import numpy as np
import scipy.integrate
import inspect

from astromodels.formula_parser import Formula
from astromodels.parameter import Parameter
from astromodels.my_yaml import my_yaml
from astromodels.named_object import NamedObject

import parser
import collections
import warnings
import sys
import functools
import pkg_resources


class WarningNoTests(Warning):
    pass


class WarningUnitsAreSlow(Warning):
    pass


class FunctionDefinitionError(Exception):
    pass


class DesignViolation(Exception):
    pass


class TestSpecificationError(Exception):
    pass


class TestFailed(Exception):
    pass


class DocstringIsNotRaw(ValueError):

    pass


class UnknownFunction(ValueError):

    pass

def input_as_array(method):
    """
    Decorator which allows the decorated functions to be coded as having always a np.array as input. However, the
    decorated function will work with any input among a list, a single float or a np.array. For example::

      > def myfunc(x):
          idx = x < 1.0
          out = np.zeros_like(x)
          out[idx] = 0.0
          out[~idx] = 1.0
          return out

    This defines a function which returns 0 for all elements below 1 and 1 above that. As it is written it only works
    with a np.array as input. Indeed::

      > print( myfunc( np.array([0,2]) ) )
      [0,1]
      > print( myfunc(1.0) )
      IndexError: too many indices for array
      > print( myfunc([0,1,2,3]) )
      array([0, 0, 0, 1])

    Calling the function with a single float as input fails, while calling it with a list does not fail but returns the
    wrong result. If we decorate the function instead:

      > dec_func = input_always_array(myfunc)
      > print( dec_func(1.0) )
      1.0
      > print( dec_func([0,1,2,3]) )
      [0 1 1 1]
      > dec_func(np.array([0,2]))
      [0,1]

    :param method: method to decorate
    :return: same as the original method, but with the dimensions squeezed to the minimum. For example, an array with
    only one element will become a single number.
    """

    # Lookup in local scope is always much faster than in module scope. Hence, to decrease the performance hit of the
    # decorator, we store the functions we are going to use in the wrapper as local variables. Note that this part
    # of the code will be executed only once during the loading of the module using the decorator, while the wrapper
    # will be executed every time the method is called

    np_array = np.array
    np_squeeze = np.squeeze

    def wrapper(input_value, *args, **kwargs):

        # Transform the input in a numpy array, if needed.
        # If the input was a single float, this will become an array with shape
        # (1,), otherwise it will keep the shape of the input.

        if isinstance(input_value, np.ndarray):

            return method(input_value, *args, **kwargs)

        else:

            new_input = np_array(input_value, ndmin=1, copy=False)

            result = method(new_input, *args, **kwargs)

            # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.

            return np_squeeze(result)

    return wrapper


def clip(**myargs):
    """
    Force the output of the method to be always between min_value and max_value. For example, to constrain a function
    to be positive (i.e. all negative results will be substituted with 0):

      @clip(min_value=0)
      def myfun(in_value):
        return -1

    To constrain a function to be between 1e-20 and 1e20:
      @clip( min_value=1e-20, max_value=1e20 )
      def function(in):
        return -20 #will return 1e-20

    :param **myargs
    :keyword min_value: minimum allowed value (default: no boundary)
    :keyword max_value: maximum allowed value (default: no boundary)
    :return: same as the input method, but the array is clipped to min_value and max_value
    """

    min_value = None if 'min_value' not in myargs.keys() else myargs['min_value']
    max_value = None if 'max_value' not in myargs.keys() else myargs['max_value']

    def clip_generator(method):
        def wrapper(*args, **kwargs):
            result = method(*args, **kwargs)

            # This will clip the result array in-place

            result.clip(min_value, max_value, out=result)

            return result

        return wrapper

    return clip_generator


def output_always_finite(method):
    """
    This calls numpy.nan_to_num on the output of method, i.e., substitute nan with 0 and inf with a large number.

    :param method: method to decorate
    :return: same as the original method, but with nan transformed to 0 and inf to a large number.
    """

    np_nan_to_num = np.nan_to_num

    def wrapper(*args, **kwargs):
        result = method(*args, **kwargs)

        return np_nan_to_num(result)

    return wrapper


# This dictionary will contain the known function by name, so that the model_parser can instance
# them by looking into this dictionary. It will be filled by the FunctionMeta meta-class.

_known_functions = {}

# The following is a metaclass for all the functions

import scipy.integrate
import inspect
import yaml as my_yaml
from yaml.reader import ReaderError
import collections
from astromodels.parameter import Parameter
import sys
import numpy as np
import functools


# The following is a metaclass for all the functions
class FunctionMeta(type):
    """
    A metaclass for the models, which takes care of setting up the parameters and the other attributes
    according to the definition given in the documentation of the function class.

    The rationale for using this solution, instead of a more usual class-inheritance mechanism, is to avoid
    re-reading and re-testing every time a new instance is created. Using a meta-class ensure that this is
    executed only once during import.
    """

    def __init__(cls, name, bases, dct):

        # This method is called at the very beginning during the import of the module,
        # after the class cls has been built.
        # We use it to perform simple checks on the function class to verify that it is
        # implemented well, and add properties and methods according to the definition
        # given in the docstring of cls

        # Enforce the presence of the evaluate method

        if 'evaluate' not in dct:

            raise AttributeError("You have to implement the 'evaluate' method in %s" % name)

        else:

            # Set the __call__ attribute to a wrapper which allows to call the function as
            # function(x), instead of function(x, par1, par2...)

            cls.__call__ = FunctionMeta.instance_call_wrapper

        # If the 'integral' method is not in there, add the default integration
        # method to the function class

        if 'integral' not in dct:

            # Use predefined numerical integral

            cls.integral = FunctionMeta.numerical_integrator

        # Minimal check of the 'evaluate' function

        variables, parameters_in_calling_sequence = FunctionMeta.check_calling_sequence(name, 'evaluate',
                                                                                         cls.evaluate,['x','y','z'])

        # Figure out the dimensionality of this function

        n_dim = len(variables)

        # Add a property to the *type* indicating the number of variables, which will be
        # propagated to the class instance as a property

        cls.n_dim = property(lambda x: n_dim, doc="Return the number of dimensions for this function.")

        # Now parse the documentation of the function which contains the parameter specification

        # The doc is a YAML document containing among other things the definition of the parameters

        # First substitute all '\' characters (which might be used in the latex formula)
        # with '\\', as otherwise the yaml load will fail

        escaped_docstring = cls.__doc__.replace(chr(92),r'\\')

        # Parse it

        try:

            function_definition = my_yaml.safe_load(escaped_docstring)

        except ReaderError:

            raise DocstringIsNotRaw("Docstring parsing has failed. " \
                                    "Did you remember to specify the docstring of %s as raw? " \
                                    "To do that, you have to put a r before the docstring, " \
                                    '''like in \n\nr"""\n(docstring)\n"""\n\ninstead of just\n\n''' \
                                    '''"""\ndocstring\n"""'''% name)

        # Enforce the presence of a description and of a parameters dictionary

        assert "description" in function_definition.keys(),"You have to provide a 'description' token in the " \
                                                           "documentation of class %s" % name

        assert "parameters" in function_definition.keys(),"You have to provide a 'parameters' token in the " \
                                                          "documentation of class %s" % name

        # If there is a latex formula, store it in the type

        if 'latex' in function_definition:

            # First remove the escaping we did to overcome the limitation of the YAML parser

            latex_formula = function_definition['latex'].replace(r"\\",chr(92))

        else:

            latex_formula = '(not provided)'

        # Add a property with the latex formula

        cls.latex = property(lambda x:latex_formula, doc="Get the formula in LaTEX (if provided).")

        # Parse the parameters' dictionary

        assert isinstance(function_definition['parameters'],dict),"Wrong syntax in 'parameters' token. It must be a " \
                                                                  "dictionary. Refers to the documentation."

        # Add the parameters as attribute of the *type* (not the instance of course, since we are working
        # on the type). During the __call__ method below this dictionary will be used to create a copy
        # of each parameter which will be made available as attribute of the instance.

        cls._parameters = collections.OrderedDict()

        for parameter_name, parameter_definition in function_definition['parameters'].iteritems():

            this_parameter = FunctionMeta.parse_parameter_definition(name, parameter_name, parameter_definition)

            cls._parameters[this_parameter.name] = this_parameter

        # Now check that all the parameters used in 'evaluate' are part of the documentation,
        # and that there are no unused parameters

        set1 = set(cls._parameters.keys())
        set2 = set(parameters_in_calling_sequence)

        if set1 != set2:

            # The parameters are different. Figure out who is missing and raise an exception accordingly

            if set1 > set2:

                missing = set1 - set2

                msg = "Parameters %s have init values but are not used in 'evaluate' in %s" % (",".join(missing), name)

            else:

                missing = set2 - set1

                msg = "Parameters %s are used in 'evaluate' but do not have init values in %s" % (",".join(missing), name)

            raise FunctionDefinitionError(msg)

        # Add the name of the function as a property

        cls.name = property(lambda x: cls.__name__)

        # Add a property returning the parameters dictionary
        cls.parameters = property(lambda instance: instance._parameters)

        # Finally add the to_dict method to serialize the class

        cls.to_dict = FunctionMeta.to_dict

        # We now proceed with the testing

        if 'tests' not in function_definition:

            warnings.warn("The function class %s contains no tests." % name, WarningNoTests)

        else:

            # Let's instance and test the class

            test_instance = cls()

            # Gather the test specifications and execute them

            for test in function_definition['tests']:

                FunctionMeta.test_simple_function(name, test, test_instance)

        # All went well, add this as a known function

        _known_functions[name] = cls

    def __call__(cls, *args, **kwargs):

        # Note that this is actually called when the class
        # cls is instanced, as a sort of decorator for the
        # constructor

        # Create the instance

        instance = type.__call__(cls, *args)

        # Create the dictionary for the parameters

        instance._parameters = collections.OrderedDict()

        # Loop over the parameters in the type, create a copy and put it in the parameters
        # dictionary. Moreover, if the value for some or all the
        # parameters have been specified in the constructor, use that so that a function can be instanced
        # as:
        # my_powerlaw = powerlaw(logK=1.2, index=-3)

        for key, value in cls._parameters.iteritems():

            # Create a copy and add to the parameters dictionary
            # Then we will add a property with the name of the parameter.
            # This is done so that the user cannot overwrite by mistake
            # the parameter

            instance._parameters[key] = value.duplicate()

            # Now add a property with the name of the parameter.

            this_setter = functools.partial(FunctionMeta.set_parameter, parameter_name=key)
            this_getter = functools.partial(FunctionMeta.get_parameter, parameter_name=key)

            setattr(cls, key, property(this_getter, this_setter, doc="Get or set %s" % key))

            if key in kwargs:

                instance._parameters[key].value = kwargs[key]

        return instance

    @staticmethod
    def to_dict(instance):

        data = collections.OrderedDict()

        for par_name, parameter in instance._parameters.iteritems():

            data[par_name] = parameter.to_dict()

        return {instance.name: data}

    @staticmethod
    def instance_call_wrapper(instance, x):

        # Gather the current parameters' values

        values = {parameter_name: parameter.value for parameter_name, parameter in instance._parameters.iteritems()}

        return instance.evaluate(x, **values)

    @staticmethod
    def check_calling_sequence(name, function_name, function, possible_variables):
        """
        Check the calling sequence for the function looking for the variables specified.
        One or more of the variables can be in the calling sequence. Note that the
        order of the variables will be enforced.
        It will also enforce that the first parameter in the calling sequence is called 'self'.

        :param function: the function to check
        :param possible_variables: a list of variables to check, The order is important, and will be enforced
        :return: a tuple containing the list of found variables, and the name of the other parameters in the calling
        sequence
        """

        # Get calling sequence

        calling_sequence = inspect.getargspec(function).args

        assert calling_sequence[0] == 'self',"Wrong syntax for 'evaluate' in %s. The first argument " \
                                             "should be called 'self'." % name

        # Figure out how many variables are used

        variables = filter(lambda var: var in possible_variables, calling_sequence)

        # Check that they actually make sense. They must be used in the same order
        # as specified in possible_variables

        assert len(variables) > 0, "The name of the variables for 'evaluate' in %s must be one or more "\
                                   "among %s" % (name,','.join(possible_variables))

        if variables != possible_variables[:len(variables)]:

            raise AssertionError("The variables %s are out of order in '%s' of %s. Should be %s."
                                 % (",".join(variables), function_name, name, possible_variables[:len(variables)]))

        other_parameters = filter(lambda var: var not in variables and var != 'self', calling_sequence)

        return variables, other_parameters

    @staticmethod
    def numerical_integrator(instance, e1, e2):

        return scipy.integrate.quad(instance.evaluate, e1, e2)[0]

    @staticmethod
    def get_parameter(instance, parameter_name):
        """
        Generic getter for a parameter
        """
        return instance._parameters[parameter_name]

    @staticmethod
    def set_parameter(instance, value, parameter_name):
        """
        Set a parameter to a new value, performing conversions if necessary
        """

        try:

            # This works if value is a astropy.Quantity

            new_value = value.to(instance._parameters[parameter_name].unit).value

        except AttributeError:

            # We get here if instead value is a simple number

            instance._parameters[parameter_name].value = value

        else:

            # Even if the to() method works, we need to warn the user that this is
            # very slow, and should only be used in interactive sessions for convenience

            warnings.warn("Using units is convenient but slow. Do not use them during computing-intensive work.",
                          WarningUnitsAreSlow)

    @staticmethod
    def parse_parameter_definition(func_name, par_name, definition):

        # Parse definition of parameter

        # Enforce the presence of attributes 'value' and 'desc'

        if 'initial value' not in definition:
            raise FunctionDefinitionError("Error for parameter %s of function %s: value for parameter must be"
                                          " specified" % (par_name, func_name))

        if 'desc' not in definition:
            raise FunctionDefinitionError("Error for parameter %s of function %s: desc for parameter must be"
                                          " specified" % (par_name, func_name))

        # Fetch attributes

        value = definition['initial value']
        desc = definition['desc']

        # Optional attributes are either None if not specified, or the value specified

        min_value = (None if 'min' not in definition else definition['min'])
        max_value = (None if 'max' not in definition else definition['max'])
        delta = (None if 'delta' not in definition else definition['delta'])
        unit = ('' if 'unit' not in definition else definition['unit'])

        # A parameter can be fixed by using fix=yes, otherwise it is free by default

        free = (True if 'fix' not in definition else not bool(definition['fix']))

        return Parameter(par_name, value, min_value=min_value, max_value=max_value,
                         delta=delta, desc=desc, free=free, unit=unit)

    @staticmethod
    def test_simple_function(name, test_specification, new_class_instance):

        # Check that all required variables are in the test

        var_names = ['x','y','z']

        for var_name in var_names[:new_class_instance.n_dim]:

            if var_name not in test_specification:
                raise TestSpecificationError("Variable %s not specified in one of the tests for %s" % (var_name, name))

        # Check that we have the minimum amount of specifications

        if 'function value' not in test_specification or \
            'tolerance' not in test_specification:

            raise TestSpecificationError("Test specification for %s lacks 'function value' or 'tolerance' attribute" %
                                         name)

        # Run the test
        # Build the point dictionary with the right number of variables

        point = {}

        for var_name in var_names[:new_class_instance.n_dim]:

            point[var_name] = float(test_specification[var_name])

        # Make a string representing the point, to be used in warnings or exceptions

        point_repr = ', '.join("{!s}={!r}".format(k, v) for (k, v) in point.iteritems())

        # Test that the value of the function is what is expected within the tolerance

        try:

            value = new_class_instance(**point)

        except TypeError:

            raise TestFailed("Cannot call function %s at point %s. Did you remember to use .value "\
                             "to access the value of a parameter in the 'evaluate' method?" % (name, point_repr))

        except:

            exc_type, value, traceback = sys.exc_info()

            raise TestFailed("Error in 'evaluate' for %s at point %s. Exception %s: '%s'" % (name, point_repr,
                                                                                             exc_type,
                                                                                             value))

        distance = value - test_specification['function value']

        if abs(distance) > float(test_specification['tolerance']):

            raise TestFailed("The function %s has value of %s instead of %s in point %s, and the difference "
                             "is larger than the tolerance %g" % (name, value,
                                                                  test_specification['function value'],
                                                                  point_repr,
                                                                  float(test_specification['tolerance'])))

        else:

            # Do nothing, test is ok
            pass


class powerlaw(object):
    r"""
    description :

        A simple power-law with normalization expressed as
        a logarithm

    latex : \frac{dN}{dx} = 10^{logK}~\frac{x}{piv}^{index}

    parameters :

        logK :

            desc : Logarithm of normalization
            initial value : 0
            min : -40
            max : 40
            unit : "1 / (keV cm2 s)"

        piv :

            desc : Pivot energy
            initial value : 1
            fix : yes
            unit: keV

        index :

            desc : Photon index
            initial value : -2
            min : -10
            max : 10

    tests :
        - { x : 10, function value: 0.01, tolerance: 1e-20}
        - { x : 100, function value: 0.0001, tolerance: 1e-20}

    """

    __metaclass__ = FunctionMeta

    def evaluate(self, x, logK, piv, index):

        return 10**logK * np.power(x / piv, index)

    def integral(self, e1, e2):

        integral =  ( 10**self.logK * np.power(e2, self.index+1) -
                      10**self.logK * np.power(e1, self.index+1) )

        return integral


def get_function(function_name):
    """
    Returns the function class "name", which must be among the known functions.

    :param function_name: the name of the function
    :return: the class (note: this is not the instance!). You have to build it yourself, like::

      my_powerlaw = get_function('powerlaw')()

    """

    if function_name in _known_functions:

        return _known_functions[function_name]

    else:

        raise UnknownFunction("Function %s is not known. Known functions are: %s" %
                              (function_name,",".join(_known_functions.keys())))