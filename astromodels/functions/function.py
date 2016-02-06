import numpy as np
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


class FunctionDefinitionError(Exception):
    pass


class DesignViolation(Exception):
    pass


class TestSpecificationError(Exception):
    pass


class TestFailed(Exception):
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


class FunctionContainer(object):

    def add_function(self, name, function_class):

        self.__setattr__(name, function_class)

    def get_function(self, name):

        return self.__getattribute__(name)


class DefinitionParser(object):
    def __init__(self, yaml_file):
        """

        :param yaml_file: resource defining the functions
        :return:
        """

        self._yaml_file = yaml_file

        # Parse the YAML file
        # (we don't catch any exception on purpose because if the file cannot be read there is nothing we can do)

        definitions = my_yaml.load(pkg_resources.resource_string('astromodels', yaml_file))

        # Now read the functions defined there

        self.functions = {}

        for func_name, func_definition in definitions.iteritems():

            # Fill the dictionary of parameters

            these_parameters = collections.OrderedDict()

            for par_name, definition in func_definition['parameters'].iteritems():
                these_parameters[par_name] = self.parse_parameter_definition(func_name, par_name, definition)

            # Now instance the function

            if 'formula' in func_definition:

                # This is a simple function

                if 'expression' in func_definition['formula']:

                    expression = func_definition['formula']['expression']

                else:

                    raise FunctionDefinitionError("No 'expression' attribute in function %s" % func_name)

                if 'latex' in func_definition['formula']:

                    latex = func_definition['formula']['latex']

                else:

                    latex = None

                if 'tests' in func_definition:

                    # Note that the proper specification of the tests will be
                    # made in the SimpleFunction constructor

                    tests = func_definition['tests']

                else:

                    tests = None

                new_class = self.simple_function_generator(func_name,
                                                           func_definition['description'],
                                                           expression,
                                                           these_parameters,
                                                           latex=latex,
                                                           tests=tests)

                self.functions[func_name] = new_class

            else:

                raise NotImplemented("Only simple functions with formulas are implemented at the moment")

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

        # A parameter can be fixed by using fix=yes, otherwise it is free by default

        free = (True if 'fix' not in definition else not bool(definition['fix']))

        return Parameter(par_name, value, min_value=min_value, max_value=max_value, delta=delta, desc=desc, free=free)

    def simple_function_generator(self, function_name, description, formula, parameters, tests=None, latex=None):
        """

        :param function_name:
        :param description:
        :param formula:
        :param parameters:
        :param tests: a dictionary containing the tests
        :param latex:
        :return:
        """

        # Parse the formula

        parsed_formula = Formula(formula)

        # Check that the parameters match between the definition received and the one needed by the formula

        set1 = set(parameters.keys())
        set2 = set(parsed_formula.parameters)

        if set1 != set2:

            # The parameters are different. Figure out who is missing and raise an exception accordingly

            if set1 > set2:

                missing = set1 - set2

                msg = "Parameters %s have init values but are not used in function" % ",".join(missing)

            else:

                missing = set2 - set1

                msg = "Parameters %s are used in function but do not have init values" % ",".join(missing)

            raise FunctionDefinitionError(msg)

        # Now create the dictionary which will be used as the locals dictionary by the __call__ method

        class_locals = {}

        # Transfer to the locals all the functions used by the formula

        for func in parsed_formula.functions:
            # Note that the Formula class already checked that all functions are members in numpy and that they
            # are ufuncs

            class_locals[func] = getattr(np, func)

        # Finally, compile the formula and test it
        # NOTE: we don't catch exceptions because we already sanitized the input, and if we fail here we want
        # to crash

        parsed = parser.expr(parsed_formula.formula)

        compiled_formula = parsed.compile()

        # Generate the new class type

        # Members

        members = collections.OrderedDict()

        # Create a duplicate of the parameters and store them in a private member. This private member will be
        # used during the SimpleFunction construction

        members['_original_parameters'] = collections.OrderedDict()

        for k, v in parameters.iteritems():
            # Duplicate the parameter to remove any ties with this function (avoiding creating a closure)

            members['_original_parameters'][k] = v.duplicate()

        # Copy the parsed formula in a private member

        members['_parsed_formula'] = parsed_formula

        # Copy the compiled formula in a private member

        members['_compiled_formula'] = compiled_formula

        # Store the name of the function

        members['_function_name'] = function_name

        # If the LATEX expression has been given, store it in the formula,
        # otherwise store None

        if latex is not None:

            members['_latex'] = latex

        else:

            members['_latex'] = None

        # Add locals

        members['_class_locals'] = class_locals

        # Now for all the parameters for this function, add a property to the class so that the user can
        # set a parameter by doing powerlaw.logK = 5.0 and get a parameter by doing powerlaw.logK

        def my_setter(name, cls, value):

            # Update the parameter in the dictionary

            cls.parameters[name].value = value

            # Update the value of the parameter in the locals

            cls._instance_locals[name] = value

        def my_getter(name, cls):
            return cls.parameters[name]

        for par_name, par in members['_original_parameters'].iteritems():
            this_setter = functools.partial(my_setter, par_name)
            this_getter = functools.partial(my_getter, par_name)

            members[par_name] = property(this_getter,
                                         this_setter,
                                         doc='Set the value of parameter %s or get its instance' % par_name)

        # This create a new class named the content of name, derived from SimpleFunction, with the members in the
        # members dictionary

        n_var = len(parsed_formula.variables)

        members['ndim'] = n_var

        if n_var == 1:

            new_class = type(function_name, (SimpleFunction1D,), members)

        elif n_var == 2:

            new_class = type(function_name, (SimpleFunction2D,), members)

        elif n_var == 3:

            new_class = type(function_name, (SimpleFunction3D,), members)

        else:

            raise FunctionDefinitionError("Only up to 3 variables are permitted")


        # Let's test it

        # First create an instance

        new_class_instance = new_class()

        # Now perform the tests, if defined, otherwise issue a warning

        if tests is None:

            warnings.warn("The function %s contains no tests." % function_name, WarningNoTests)

        else:

            # Loop over the tests and execute them

            for test in tests:
                # this function will raise if the test is failed

                self.test_simple_function(function_name, test, parsed_formula.variables, new_class_instance)

        return new_class

    @staticmethod
    def test_simple_function(name, test_specification, variables, new_class_instance):

        # Check that all required variables are in the test

        for var in variables:

            if var not in test_specification:
                raise TestSpecificationError("Variable %s not specified in one of the tests for %s" % (var, name))

        # Check that we have the minimum amount of specifications

        if 'function value' not in test_specification or \
                        'tolerance' not in test_specification:
            raise TestSpecificationError("Test specification for %s lacks 'function value' or 'tolerance' attribute" %
                                         name)

        # Run the test
        # Build the point dictionary with the right number of variables

        point = {}

        for var in variables:

            point[var] = test_specification[var]

        # Make a string representing the point, to be used in warnings or exceptions

        point_repr = ', '.join("{!s}={!r}".format(k, v) for (k, v) in point.iteritems())

        # Test that the value of the function is what is expected within the tolerance

        try:

            value = new_class_instance(**point)

        except:

            exc_type, value, traceback = sys.exc_info()

            raise TestFailed("Cannot call function at point %s because of exception %s: '%s'" % (point_repr,
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


class Function(NamedObject):

    def __init__(self, function_name):

        super(Function, self).__init__(function_name, allow_spaces=False)

    @property
    def parameters(self):
        return self._parameters

    def __getitem__(self, item):

        return self._parameters[item]

    def to_dict(self):

        data = collections.OrderedDict()

        for par_name, parameter in self.parameters.iteritems():

            data[par_name] = parameter.to_dict()

        return {self.name: data}


class SimpleFunction(Function):

    def __init__(self):

        # This class must only be created by the simple_function_generator function. Such function for example
        # add a private member _original_parameters which contains the parameters specified in the yaml file
        # describing the function. In this constructor we need to duplicate the parameters so that we have copies
        # that we can change

        # Enforce proper use

        if not hasattr(self, '_original_parameters'):
            raise DesignViolation("You cannot instance SimpleFunction directly, but only through "
                                  "the simple_function_generator function.")

        # Duplicate the parameters so we remove the link between the different instances of this function,
        # which share the _original_parameters

        self._parameters = collections.OrderedDict()

        for k, v in self._original_parameters.iteritems():
            self._parameters[k] = v.duplicate()

        # Now create a new local dictionary to keep the locals for this instance

        self._instance_locals = {}

        # Add the members of the __class_locals for a faster access during the __call__ method

        for k,v in self._class_locals.iteritems():
            self._instance_locals[k] = v

        # Now set the parameters in the locals to the default values

        for par_name, par in self._parameters.iteritems():

            self._instance_locals[par_name] = par.value

        super(SimpleFunction, self).__init__(self._function_name)

    def __repr__(self):

        representation = "Function %s:\n" % self.name
        representation += "    -parameters: %s\n" % ",".join(self.parameters.keys())

        return representation


class SimpleFunction1D(SimpleFunction):
    def __call__(self, x):
        # Transfer variables to the local dictionary used by eval

        self._instance_locals['x'] = x

        # Note that when parameters change values, they are updated in the _class_locals dictionary by
        # the setter defined in the simple_function_generator

        # Using an empty globals dictionary speed things up, because looking into the locals dictionary is much
        # faster than looking in the globals one

        return eval(self._compiled_formula, {}, self._instance_locals)


class SimpleFunction2D(SimpleFunction):
    def __call__(self, x, y, z):

        # Transfer variables to the local dictionary used by eval

        self._instance_locals['x'] = x
        self._instance_locals['y'] = y

        # Note that when parameters change values, they are updated in the _class_locals dictionary by
        # the setter defined in the simple_function_generator

        # Using an empty globals dictionary speed things up, because looking into the locals dictionary is much
        # faster than looking in the globals one

        return eval(self._compiled_formula, {}, self._instance_locals)


class SimpleFunction3D(SimpleFunction):
    def __call__(self, x, y, z):
        # Transfer variables to the local dictionary used by eval

        self._instance_locals['x'] = x
        self._instance_locals['y'] = y
        self._instance_locals['z'] = z

        # Note that when parameters change values, they are updated in the _class_locals dictionary by
        # the setter defined in the simple_function_generator

        # Using an empty globals dictionary speed things up, because looking into the locals dictionary is much
        # faster than looking in the globals one

        return eval(self._compiled_formula, {}, self._instance_locals)


# Init this here to be a global variable, so it can be imported with "from astromodels.functions.function import f1d"

f1d = None
f2d = None
f3d = None


def build_function_containers():
    global f1d, f2d, f3d

    # Build the function containers

    f1d = FunctionContainer()

    yaml_files = pkg_resources.resource_listdir('astromodels', 'data/functions/')

    for yaml_file in yaml_files:

        df = DefinitionParser('data/functions/' + yaml_file)

        for func_name, func_class in df.functions.iteritems():

            if func_class.ndim==1:

                f1d.add_function(func_name, func_class)

            elif func_class.ndim==2:

                f2d.add_function(func_name, func_class)

            elif func_class.ndim==3:

                f3d.add_function(func_name, func_class)

            else:

                raise NotImplementedError("Cannot handle a function of type %s" % type(func_class))



# Run this on import
build_function_containers()

# Now add all the functions to the dictionary of this module, so pickle will find them

# __dict__['powerlaw'] = f1d.powerlaw