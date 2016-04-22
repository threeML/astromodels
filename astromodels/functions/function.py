from astromodels.parameter import Parameter
from astromodels.my_yaml import my_yaml
from astromodels.utils.pretty_list import dict_to_list
from astromodels.tree import Node
from astromodels.utils.table import dict_to_table
from astromodels.units import get_units
import astropy.units as u


import collections
import warnings
import sys
import os
import copy
import uuid
import re
import ast
from yaml.reader import ReaderError
import numpy as np
import scipy.integrate
import inspect

__author__ = 'giacomov'


try:

    from IPython.display import display, HTML

except:

    has_ipython = False

else:

    has_ipython = True


class WarningNoTests(ImportWarning):
    pass


class FunctionDefinitionError(Exception):
    pass


class DesignViolation(Exception):
    pass


class WrongDimensionality(Exception):
    pass


class TestSpecificationError(Exception):
    pass


class TestFailed(Exception):
    pass


class DocstringIsNotRaw(ValueError):
    pass


class UnknownFunction(ValueError):
    pass


class UnknownParameter(ValueError):
        pass


# Value to indicate that no latex formula has been given
NO_LATEX_FORMULA = '(no latex formula available)'

# Codes to indicate to Composite Function the operation between two functions
_operations = {'+': np.add,
               '-': np.subtract,
               '*': np.multiply,
               '/': np.divide,
               '**': np.power,
               'abs': np.abs,
               'of': 'compose'}


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

        # Minimal check of the 'evaluate' function

        variables, parameters_in_calling_sequence = FunctionMeta.check_calling_sequence(name, 'evaluate',
                                                                                        cls.evaluate, ['x', 'y', 'z'])

        #print("Function %s: \n variables: %s\n parameters: %s\n" % (cls.__name__, ",".join(variables),
        #                                                            ",".join(parameters_in_calling_sequence)))

        # Figure out the dimensionality of this function

        n_dim = len(variables)

        # Store the dimensionality in the *type*

        cls._n_dim = n_dim

        # Now parse the documentation of the function which contains the parameter specification

        # The doc is a YAML document containing among other things the definition of the parameters

        # Parse it

        try:

            function_definition = my_yaml.load(cls.__doc__)

        except ReaderError:

            raise DocstringIsNotRaw("Docstring parsing has failed. "
                                    "Did you remember to specify the docstring of %s as raw? "
                                    "To do that, you have to put a r before the docstring, "
                                    '''like in \n\nr"""\n(docstring)\n"""\n\ninstead of just\n\n'''
                                    '''"""\ndocstring\n"""''' % name)

        else:

            # Store the function definition in the type

            cls._function_definition = function_definition

        # Enforce the presence of a description and of a parameters dictionary

        assert "description" in function_definition.keys(), "You have to provide a 'description' token in the " \
                                                            "documentation of class %s" % name

        assert "parameters" in function_definition.keys(), "You have to provide a 'parameters' token in the " \
                                                           "documentation of class %s" % name

        # If there is a latex formula, store it in the type

        if 'latex' in function_definition:

            # First remove the escaping we did to overcome the limitation of the YAML parser

            latex_formula = function_definition['latex'].replace(r"\\", chr(92))

        else:

            latex_formula = NO_LATEX_FORMULA

        # Store latex formula in the type
        cls._latex = latex_formula

        # Store the name in the type
        cls._name = cls.__name__

        # Parse the parameters' dictionary
        assert isinstance(function_definition['parameters'], dict), "Wrong syntax in 'parameters' token. It must be " \
                                                                    "a dictionary. Refers to the documentation."

        # Add the parameters as attribute of the *type*. During the __call__ method below this dictionary will be used
        # to create a copy of each parameter which will be made available as attribute of the instance.

        cls.__parameters = collections.OrderedDict()

        for parameter_name, parameter_definition in function_definition['parameters'].iteritems():
            this_parameter = FunctionMeta.parse_parameter_definition(name, parameter_name, parameter_definition)

            cls.__parameters[this_parameter.name] = this_parameter

        # Now check that all the parameters used in 'evaluate' are part of the documentation,
        # and that there are no unused parameters

        set1 = set(cls.__parameters.keys())
        set2 = set(parameters_in_calling_sequence)

        if set1 != set2:

            # The parameters are different. Figure out who is missing and raise an exception accordingly

            if set1 > set2:

                missing = set1 - set2

                msg = "Parameters %s have init values but are not used in 'evaluate' in %s" % (",".join(missing), name)

            else:

                missing = set2 - set1

                msg = "Parameters %s are used in 'evaluate' but do not have init values in %s" % \
                      (",".join(missing), name)

            raise FunctionDefinitionError(msg)

        # Now add the constructor to the class
        cls.__init__ = FunctionMeta.class_init

        # Finally, add the info() method to the type so that it can be called even without instancing the class

        def info():

            repr_dict = collections.OrderedDict()

            repr_dict['description'] = function_definition['description']

            if 'latex' in function_definition:
                repr_dict['formula'] = function_definition['latex']

            # Add the description of each parameter and their current value
            repr_dict['default parameters'] = collections.OrderedDict()

            for parameter_name in cls.__parameters.keys():

                repr_dict['default parameters'][parameter_name] = cls.__parameters[parameter_name].to_dict()

            if has_ipython:

                display(HTML(dict_to_list(repr_dict, html=True)))

            else:

                print(dict_to_list(repr_dict, html=False))

        cls.info = staticmethod(info)

        # We now proceed with the testing

        if 'tests' not in function_definition:

            warnings.warn("The function class %s contains no tests." % name, WarningNoTests)

        else:

            # Let's instance and test the class

            test_instance = cls()

            # First add the default parameters which will be used during the test

            #test_instance._parameters = cls.__parameters

            # Gather the test specifications and execute them

            for test in function_definition['tests']:

                FunctionMeta.test_simple_function(name, test, test_instance)

        # All went well, add this as a known function

        _known_functions[name] = cls

        # Finally call the type init

        super(FunctionMeta, cls).__init__(name, bases, dct)

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

        assert calling_sequence[0] == 'self', "Wrong syntax for 'evaluate' in %s. The first argument " \
                                              "should be called 'self'." % name

        # Figure out how many variables are used

        variables = filter(lambda var: var in possible_variables, calling_sequence)

        # Check that they actually make sense. They must be used in the same order
        # as specified in possible_variables

        assert len(variables) > 0, "The name of the variables for 'evaluate' in %s must be one or more " \
                                   "among %s" % (name, ','.join(possible_variables))

        if variables != possible_variables[:len(variables)]:
            raise AssertionError("The variables %s are out of order in '%s' of %s. Should be %s."
                                 % (",".join(variables), function_name, name, possible_variables[:len(variables)]))

        other_parameters = filter(lambda var: var not in variables and var != 'self', calling_sequence)

        return variables, other_parameters

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

        var_names = ['x', 'y', 'z']

        for var_name in var_names[:new_class_instance.n_dim]:

            if var_name not in test_specification:
                raise TestSpecificationError("Variable %s not specified in one of the tests for %s" % (var_name, name))

        # Check that we have the minimum amount of specifications

        if 'function value' not in test_specification or 'tolerance' not in test_specification:

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

            raise TestFailed("Cannot call function %s at point %s. Did you remember to use .value "
                             "to access the value of a parameter in the 'evaluate' method?" % (name, point_repr))

        except:

            exc_type, value, traceback = sys.exc_info()

            raise TestFailed("Error in 'evaluate' for %s at point %s. Exception %s: '%s'" % (name, point_repr,
                                                                                             exc_type,
                                                                                             value))

        # The eval(str()) bit is needed so that a function value can be specified as np.inf or np.pi

        distance = value - eval(str(test_specification['function value']))

        if abs(distance) > float(test_specification['tolerance']):

            raise TestFailed("The function %s has value of %s instead of %s in point %s, and the difference "
                             "%s is larger than the tolerance %g" % (name, value,
                                                                     test_specification['function value'],
                                                                     point_repr,
                                                                     abs(distance),
                                                                     float(test_specification['tolerance'])))

        else:

            # Do nothing, test is ok
            pass

    @staticmethod
    def class_init(instance, **kwargs):

        # Create a copy of the parameters dictionary which is in the type,
        # otherwise every instance would share the same dictionary

        copy_of_parameters = collections.OrderedDict()

        # Fill it by duplicating the parameters contained in the dictionary in the type

        for key, value in type(instance).__parameters.iteritems():

            copy_of_parameters[key] = value.duplicate()

            # If the user has specified a value in the constructor, update the
            # corresponding parameters value. This allow to use a constructor as:
            # my_powerlaw = powerlaw(logK=1.0, index=-2)

            if key in kwargs:

                copy_of_parameters[key].value = kwargs[key]

        # Now check that all the parameters specified in the kwargs are actually parameters of this function
        for key in kwargs.keys():

            try:

                copy_of_parameters[key]

            except KeyError:

                raise UnknownParameter("You specified an init value for %s, which is not a "
                                       "parameter of function %s" % (key, type(instance)._name))

        # Now call the parent class

        Function.__init__(instance,
                          type(instance)._name,
                          type(instance)._function_definition,
                          copy_of_parameters,
                          type(instance)._n_dim)



class Function(Node):

    def __init__(self, name, function_definition, parameters, n_dim):

        # (this is called by the constructor defined in the metaclass)

        # Store name, number of dimensions and the latex formula

        # Note; in a normal situation these are stored in the type already. Thus, this is a small waste of memory.
        # However, doing this will allow to subclass the Function class without using the FunctionMeta meta-class.

        self._n_dim = n_dim

        # Store also the function definition

        assert 'description' in function_definition,"Function definition must contain a description"

        assert 'latex' in function_definition, "Function definition must contain a latex formula"

        self._function_definition = function_definition

        # Set up the node

        Node.__init__(self, name)

        # Add the parameters as children. Since the name of the key in the dictionary might
        # be different than the actual name of the parameter, use the .add_child method instead
        # of the add_children method

        for child_name, child in parameters.iteritems():

            self.add_child(child, child_name)

        # Now generate a unique identifier (UUID) in a thread safe, multi-processing safe
        # way. This is used for example in the CompositeFunction class to keep track of the different
        # instances of the same function
        self._uuid = "{" + str(self._generate_uuid()) + "}"

        # This will contain the units for the independent variables x(,y,z)
        self._independent_variables_unit = [None] * self.n_dim

    def set_independent_variables_unit(self, *units):
        """
        Set the independent variables for this function

        :param units: as many strings (like 'keV') or astropy.units.Unit instances as needed depending on the number
        of dimensions of this function.
        :return: none
        """

        assert len(units) == self.n_dim, "Wrong number of variables for function with %s dimension(s)" % self.n_dim

        self._independent_variables_unit = map(lambda x:u.Unit(x), units)

    @property
    def independent_variables_unit(self):
        """
        Returns a list with the units of the independent variables for this function

        :return: a list
        """
        return self._independent_variables_unit

    @property
    def free_parameters(self):
        """
        Returns a dictionary of free parameters for this function

        :return: dictionary of free parameters
        """

        free_parameters = collections.OrderedDict([(k,v) for k, v in self.parameters.iteritems() if v.free])

        return free_parameters

    def get_wrapper(self):
        """
        Returns a python function which can be used to call this function with the parameters in the calling sequence.
        In other words, if you can call this function with f(x), with the returned wrapper you can call it with
        f(x, parameter1, parameter2...), where parameter1, parameter2... are the *free* parameters.

        :return: a python function
        """

        # Build a list of free parameters
        free_parameters = [k for k,v in self.parameters.iteritems() if v.free]

        # Prepare the method to set the parameters to their current value

        def set_parameters(*args):

            [self.get_child(free_parameters[i]).set_value(args[i]) for i in range(len(args))]

        # Prepare the variable description

        variables = None

        if self.n_dim == 1:

            variables = 'x'

        elif self.n_dim == 2:

            variables = 'x,y'

        elif self.n_dim == 3:

            variables = 'x,y,z'

        # Prepare the free parameters string

        free_parameters_string = ",".join(free_parameters)

        # Build some code to generate a wrapper which takes care of updating the value of the parameters
        # and return the value of the function

        wrapper_code = '''

        def wrapper(%s, %s):

            set_parameters(%s)

            return self(%s)

        ''' % (variables, free_parameters_string, free_parameters_string, variables)

        #print(wrapper_code)

        exec(wrapper_code.replace("        ","")) in locals()

        return wrapper


    @staticmethod
    def _generate_uuid():
        """
        Generate a unique identifier for this function.

        :return: the UUID
        """
        return uuid.UUID(bytes=os.urandom(16), version=4)

    def numerical_integrator(self, e1, e2):

        return scipy.integrate.quad(self.__call__, e1, e2)[0]

    @property
    def description(self):
        """
        Returns a description for this function
        """

        return self._function_definition['description']

    # Add a property returning the parameters dictionary
    @property
    def parameters(self):
        """
        Returns a dictionary of parameters
        """
        return self.children

    @property
    def n_dim(self):
        """
        Returns the number of dimensions for this function (1, 2 or 3)
        """
        return self._n_dim

    @property
    def latex(self):
        """
        Returns the LaTEX formula for this function
        """
        return self._function_definition['latex']

    def evaluate(self, *args, **kwargs):

        raise NotImplementedError("You have to re-implement this")

    def get(self, *args, **kwargs):
        """
        Evaluate the function with units

        :return:
        """

        # Get the active internal units
        internal_units = get_units()

        converted = []

        for i in range(len(args)):

            if not isinstance(args[i],u.Quantity):

                raise TypeError("If you use .get() you have to provide astropy quantities (with units)")

            # This is either time or energy, normally
            try:

                physical_type = args[i].physical_type

            except AttributeError:

                raise TypeError("You cannot use get() with units without physical type")

            try:

                internal_unit = internal_units.__getattribute__(physical_type)

            except AttributeError:

                raise TypeError("The physical unit %s is not an elementary unit for astromodels" % physical_type)

            this_converted = args[i].to(internal_unit)

            converted.append(this_converted)

        # Gather the current parameters' values

        kwargs.update({parameter_name: parameter.value * parameter.unit for parameter_name, parameter
                       in self.children.iteritems()})

        return self.evaluate(*args, **kwargs)


    def __call__(self, *args, **kwargs):

        # Gather the current parameters' values

        kwargs.update({parameter_name: parameter.value for parameter_name, parameter in self.children.iteritems()})

        return self.evaluate(*args, **kwargs)

    # Define now all the operators which allow to combine functions. Each operator will return a new
    # instance of a CompositeFunction, which can then be used as a function on its own

    @staticmethod
    def __get_second_uuid(other_instance):
        """
        Return a name for the object. If the object is a function instance, return its name. Otherwise, return
        the object itself.

        :param other_instance:
        :return:
        """
        if hasattr(other_instance,'uuid'):

            # Another function

            second_uuid = other_instance.uuid

        else:

            # A number

            second_uuid = '%s' % other_instance

        return second_uuid

    def of(self, another_function):
        """
        Compose this function with another as in this_function(another_function(x))

        :param another_function: another function to compose with the current one
        :return: a composite function instance
        """
        return CompositeFunction('of', self, another_function)

    def __neg__(self):

        return CompositeFunction('-', self)

    def __abs__(self):

        return CompositeFunction('abs', self)

    def __pow__(self, other_instance):

        second_uuid = self.__get_second_uuid(other_instance)

        return CompositeFunction('**', self, other_instance)

    def __rpow__(self, other_instance):

        return CompositeFunction('**', other_instance, self)

    def __add__(self, other_instance):
        """
        Return a composite function where the current instance is summed with the given instance

        :param other_instance: the other instance. This can be either a number, or a Function instance.
        """

        return CompositeFunction('+', self, other_instance)

    __radd__ = __add__

    def __sub__(self, other_instance):

        return CompositeFunction('-', self, other_instance)

    __rsub__ = __sub__

    def __mul__(self, other_instance):

        return CompositeFunction('*', self, other_instance)

    __rmul__ = __mul__

    def __div__(self, other_instance):

        return CompositeFunction('/', self, other_instance)

    def __rdiv__(self, other_instance):

        return CompositeFunction('/', other_instance, self)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def _repr__base(self, rich_output):

        repr_dict = collections.OrderedDict()

        repr_dict['description'] = self._function_definition['description']

        if 'latex' in self._function_definition:

            repr_dict['formula'] = self._function_definition['latex']

        # Add the description of each parameter and their current value
        repr_dict['parameters'] = collections.OrderedDict()

        for parameter_name in self.children.keys():

            repr_dict['parameters'][parameter_name] = self.children[parameter_name].to_dict()

        return dict_to_list(repr_dict, rich_output)

    @property
    def uuid(self):
        """
        Returns the ID of the current function. The ID is used by the CompositeFunction class to keep track of the
        unique instances of each function. It should not be used by the user for any specific purpose.

        :return: (none)
        """
        return self._uuid

    def duplicate(self):
        """
        Create a copy of the current function with all the parameters equal to the current value

        :return: a new copy of the function
        """

        # Create a copy

        function_copy = copy.deepcopy(self)

        return function_copy


class CompositeFunction(Function):

    def __init__(self, operation, function_or_scalar_1, function_or_scalar_2=None):

        assert operation in _operations,"Do not know operation %s" % operation

        # Set the new evaluate

        if function_or_scalar_2 is None:

            # Unary operation

            self.set_evaluate(self._composite_function_factory_unary(function_or_scalar_1, _operations[operation]))

        else:

            # Binary operation

            # Check if this is a function composition with the .of method
            if _operations[operation] == 'compose':

                # Can only compose functions of 1 variable

                assert hasattr(function_or_scalar_2, 'evaluate'), "Second member of .of cannot be a scalar"

                assert function_or_scalar_1.n_dim == 1 and function_or_scalar_2.n_dim == 1, "Can only compose " \
                                                                                            "with .of functions of " \
                                                                                            "1 variable"

                def new_evaluate(*args):

                    value = function_or_scalar_2(*args)

                    return function_or_scalar_1(value, *(args[1:]))

                self.set_evaluate(new_evaluate)

            else:

                self.set_evaluate(self._composite_function_factory_binary(function_or_scalar_1, _operations[operation],
                                                                          function_or_scalar_2))

        # Save a description, but using the unique IDs of the functions involved, to keep track
        # of where they appear in the expression

        self._uuid_expression = self._get_uuid_expression(operation, function_or_scalar_1, function_or_scalar_2)

        # Makes the list of unique functions which compose this composite function.

        self._functions = []

        for function in [function_or_scalar_1, function_or_scalar_2]:

            # Check whether this is already a composite function. If it is, add the functions contained
            # in it

            if isinstance(function, CompositeFunction):

                for sub_function in function.functions:

                    if sub_function not in self._functions:

                        self._functions.append(sub_function)

            elif isinstance(function, Function):

                # This is a simple function. Add it only if it is not there already (avoid duplicate)

                if function not in self._functions:

                    self._functions.append(function)

            else:

                # This is a scalar, no need to add it among the functions

                pass

        # Check that the functions have all the same dimensionality

        n_dim = self._functions[0].n_dim

        for function in self._functions[1:]:

            if function.n_dim != n_dim:

                raise WrongDimensionality("Dimensionality mismatch when composing functions.")

        # Now assign a unique name to all the functions, to make clear which is which in the definition
        # and give an easy way for the user to understand which parameter belongs to which function

        self._id_to_uid = {}

        expression = self._uuid_expression

        for i,function in enumerate(self._functions):

            self._id_to_uid[i+1] = function.uuid

            expression = expression.replace(function.uuid, "%s{%s}" % (function.name, i+1))

        # Save the expression
        self._expression = expression

        # Build the parameters dictionary assigning a new name to each parameter to account for possible
        # duplicates.

        parameters = collections.OrderedDict()

        for i, function in enumerate(self._functions):

            for parameter_name, parameter in function.parameters.iteritems():

                # New name to avoid possible duplicates

                new_name = "%s_%i" % (parameter.name, i+1)

                # Store the parameter under the new name (obviously this is a reference to the
                # parameter, not a copy, as always in python)

                parameters[new_name] = parameter
                parameter.change_name(new_name)

        # Now build a meaningful description

        _function_definition = {'description': self.expression, 'latex': NO_LATEX_FORMULA}

        Function.__init__(self, 'composite', _function_definition, parameters, n_dim)

        self._uuid = self._uuid_expression

    def _set_units(self, x_unit, y_unit):

        # Just rely on the single functions to adjust themselves.

        for function in self.functions:

            function._set_units(x_unit, y_unit)

    @property
    def expression(self):
        return self._expression

    @staticmethod
    def _get_uuid_expression(operation, name_1, name_2=None):

        if name_2 is None:

            return '(%s %s)' % (operation, name_1.uuid)

        if hasattr(name_1, 'uuid'):

            name_1_uuid = name_1.uuid

        else:

            name_1_uuid = '%s' % name_1

        if hasattr(name_2,'uuid'):

            name_2_uuid = name_2.uuid

        else:

            name_2_uuid = '%s' % name_2

        return '(%s %s %s)' % (name_1_uuid, operation, name_2_uuid)

    @staticmethod
    def _composite_function_factory_unary(instance, numpy_operator):

        def new_evaluate(*args, **kwargs):

            return numpy_operator(instance.__call__(*args, **kwargs))

        return new_evaluate

    @staticmethod
    def _composite_function_factory_binary(first_instance, numpy_operator, second_instance):

        # Check whether the second member is a function, or a number

        if hasattr(first_instance, 'evaluate'):

            if hasattr(second_instance, 'evaluate'):

                # Get number of arguments for first function


                def new_evaluate(*args):

                    return numpy_operator(first_instance.__call__(*args),
                                          second_instance.__call__(*args))

                return new_evaluate

            else:

                def new_evaluate(*args):

                    return numpy_operator(first_instance.__call__(*args),
                                          second_instance)

                return new_evaluate

        else:

            if hasattr(second_instance, 'evaluate'):

                def new_evaluate(*args):

                    return numpy_operator(first_instance,
                                          second_instance.__call__(*args))

                return new_evaluate

            else:

                # Should never get here!

                raise RuntimeError("Should never get here")

    @property
    def functions(self):
        "A list containing the function used to build this composite function"
        return self._functions

    def evaluate(self):

        raise NotImplementedError("You cannot instance and use a composite function by itself. Use the factories.")

    def set_evaluate(self, new_evaluate_method):
        """
        This is called by the factory which create the composite function, and set the evaluate method
        to the new method which will collect the results from all the functions
        """

        self.evaluate = new_evaluate_method

    # Override the __call__ method of the Function class because the single functions in _functions
    # will handle their own collection of parameters

    def __call__(self, *args, **kwargs):

        return self.evaluate(*args, **kwargs)

    # Override the to_dict method of the Node class to add the expression to re-build this
    # composite function
    def to_dict(self, minimal=False):

        data = super(CompositeFunction, self).to_dict(minimal)

        if not minimal:

            data['expression'] = self._expression

        return data


def get_function(function_name, composite_function_expression=None):
    """
    Returns the function class "name", which must be among the known functions.

    :param function_name: the name of the function (use 'composite' if the function is a composite function)
    :param composite_function_expression: composite function specification such as
    ((((powerlaw{1} + (sin{2} * 3)) + (sin{2} * 25)) - (powerlaw{1} * 16)) + (sin{2} ** 3.0))
    :return: the class (note: this is not the instance!). You have to build it yourself, like::

      my_powerlaw = get_function('powerlaw')()

    """

    # Check whether this is a composite function or a simple function
    if composite_function_expression is not None:

        # Composite function

        return _parse_function_expression(composite_function_expression)

    else:

        if function_name in _known_functions:

            return _known_functions[function_name]()

        else:

            raise UnknownFunction("Function %s is not known. Known functions are: %s" %
                                  (function_name, ",".join(_known_functions.keys())))


def list_functions():

    # Gather all defined functions and their descriptions

    functions_and_descriptions = {key:{'Description': value._function_definition['description']}
                                  for key,value in _known_functions.iteritems()}

    # Order by key (i.e., by function name)

    ordered = collections.OrderedDict(sorted(functions_and_descriptions.items()))

    # Format in a table

    table = dict_to_table(ordered)

    return table


def _parse_function_expression(function_specification):
    """
    Parse a complex function expression like:

    ((((powerlaw{1} + (sin{2} * 3)) + (sin{2} * 25)) - (powerlaw{1} * 16)) + (sin{2} ** 3.0))

    and return a composite function instance

    :param function_specification:
    :return: a composite function instance
    """

    # NOTE FOR SECURITY
    # This function has some security concerns. Security issues could arise if the user tries to read a model
    # file which has been maliciously formatted to contain harmful code. In this function we close all the doors
    # to a similar attack, except for those attacks which assume that the user has full access to a python environment.
    # Indeed, if that is the case, then the user can already do harm to the system, and so there is no point in
    # safeguard that from here. For example, the user could format a subclass of the Function class which perform
    # malicious operations in the constructor, add that to the dictionary of known functions, and then interpret
    # it with this code. However, if the user can instance malicious classes, then why would he use astromodels to
    # carry out the attack? Instead, what we explicitly check is the content of the function_specification string,
    # so that it cannot by itself do any harm (by for example containing instructions such as os.remove).

    # This can be a arbitrarily complex specification, like
    # ((((powerlaw{1} + (sin{2} * 3)) + (sin{2} * 25)) - (powerlaw{1} * 16)) + (sin{2} ** 3.0))

    # Use regular expressions to extract the set of functions like function_name{number},
    # then build the set of unique functions by using the constructor set()

    unique_functions = set(re.findall(r'\b([a-zA-Z0-9]+)\{([0-9]?)\}',function_specification))

    # NB: unique functions is a set like:
    # {('powerlaw', '1'), ('sin', '2')}

    # Create instances of the unique functions

    instances = {}

    # Loop over the unique functions and create instances

    for (unique_function, number) in unique_functions:

        complete_function_specification = "%s{%s}" % (unique_function, number)

        # As first safety measure, check that the unique function is in the dictionary of _known_functions.
        # This could still be easily hacked, so it won't be the only check

        if unique_function in _known_functions:

            # Get the function class and check that it is indeed a proper Function class

            function_class = _known_functions[unique_function]

            if issubclass(function_class, Function):

                # Ok, let's create the instance

                instance = function_class()

                # Append the instance to the list

                instances[complete_function_specification] = instance

            else:

                raise FunctionDefinitionError("The function specification %s does not contain a proper function"
                                              % unique_function )

        else:

            raise UnknownFunction("Function %s in expression %s is unknown" % (unique_function, function_specification))

    # Check that we have found at least one instance.

    if len(instances)==0:

        raise DesignViolation("No known function in function specification")

    # The following presents a slight security problem if the model file that has been parsed comes from an untrusted
    # source. Indeed, the use of eval could make possible to execute things like os.remove.
    # In order to avoid this, first we substitute the function instances with numbers and remove the operators like
    # +,-,/ and so on. Then we try to execute the string with ast.literal_eval, which according to its documentation:

    # Safely evaluate an expression node or a Unicode or Latin-1 encoded string containing a Python literal or
    # container display. The string or node provided may only consist of the following Python literal structures:
    # strings, numbers, tuples, lists, dicts, booleans, and None.This can be used for safely evaluating strings
    # containing Python values from untrusted sources without the need to parse the values oneself.
    # It is not capable of evaluating arbitrarily complex expressions, for example involving operators or indexing.

    # If literal_eval cannot parse the string, it means that it contains unsafe input.

    # Create a copy of the function_specification

    string_for_literal_eval = function_specification

    # Remove from the function_specification all the known operators and function_expressions, and substitute them
    # with a 0 and a space

    # Let's start from the function expression

    for function_expression in instances.keys():

        string_for_literal_eval = string_for_literal_eval.replace(function_expression, '0 ')

    # Now remove all the known operators

    for operator in _operations.keys():

        string_for_literal_eval = string_for_literal_eval.replace(operator,'0 ')

    # The string at this point should contains only numbers and parenthesis separated by one or more spaces

    if re.match('''([a-zA-Z]+)''', string_for_literal_eval):

        raise DesignViolation("Extraneous input in function specification")

    # By using split() we separate all the numbers and parenthesis in a list, then we join them
    # with a comma, to end up with a comma-separated list of parenthesis and numbers like:
    # ((((0,0,(0,0,3)),0,(0,0,25)),0,(0,0,16)),0,(0,0,0,3.0))
    # This string can be parsed by literal_eval as a tuple containing other tuples, which is fine.
    # If the user has inserted some malicious content, like os.remove or more weird stuff like code objects,
    # the parsing will fail

    string_for_literal_eval = ",".join(string_for_literal_eval.split())

    print(string_for_literal_eval)

    # At this point the string should be just a comma separated list of numbers

    # Now try to execute the string
    try:

        ast.literal_eval(string_for_literal_eval)

    except (ValueError, SyntaxError):

        raise DesignViolation("The given expression is not a valid function expression")

    else:

        # The expression is safe, let's eval it

        # First substitute the reference to the functions (like 'powerlaw{1}') with a string
        # corresponding to the instance dictionary

        sanitized_function_specification = function_specification

        for function_expression in instances.keys():

            sanitized_function_specification = sanitized_function_specification.replace(function_expression,
                                                                                        'instances["%s"]' %
                                                                                        function_expression)

        # Now eval it. For safety measure, I remove all globals, and the only local is the 'instances' dictionary

        composite_function = eval(sanitized_function_specification, {}, {'instances': instances})

        return composite_function