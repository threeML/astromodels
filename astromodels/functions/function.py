from astromodels.parameter import Parameter
from astromodels.my_yaml import my_yaml
from astromodels.utils.pretty_list import dict_to_list
from astromodels.tree import Node
from astromodels.utils.table import dict_to_table
from astromodels.units import get_units
import astropy.units as u
import functools


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


class ModelAssertionViolation(Exception):
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


# This dictionary will contain the known function by name, so that the model_parser can instance
# them by looking into this dictionary. It will be filled by the FunctionMeta meta-class.

_known_functions = {}


def memoize(method):
    """
    A decorator for the 2d functions which memoize the results of 1 call (useful when the minimizer is taking partial
    derivatives and calls the 2d function several times with the same arguments)

    :param method: method to be memoized
    :return: the decorated method
    """

    # Need more testing for this
    # TODO: test memoization

    # return method

    cache = method.cache = collections.OrderedDict()

    @functools.wraps(method)
    def memoizer(instance, x, *args, **kwargs):

        # Create a tuple because a tuple is hashable

        unique_id = tuple(x.value for x in instance.parameters.values()) + (x.shape[0],x.min(),x.max())

        # Create a unique identifier for this combination of inputs

        key = hash(unique_id)

        # Let's do it this way so we only look into the dictionary once

        result = cache.get(key)

        if result is not None:

            return result

        else:

            result = method(instance, x, *args, **kwargs)

            cache[key] = result

            if len(cache) > 10:

                # Remove the 10 elements that were put in first

                [cache.popitem(False) for i in range(5)]

            return result

    # Add the function as a "attribute" so we can access it
    memoizer.input_object = method

    return memoizer


# The following is a metaclass for all the functions
class FunctionMeta(type):
    """
    A metaclass for the models, which takes care of setting up the parameters and the other attributes
    according to the definition given in the documentation of the function class.
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

        # We also need the method _set_units
        if '_set_units' not in dct:

            raise AttributeError("You have to implement the '_set_units' method in %s" % name)

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

        # If the function has been memoized, it will have a "input_object" member

        try:

            calling_sequence = inspect.getargspec(function.input_object).args

        except AttributeError:

            # This might happen if the function is with memoization

            calling_sequence = inspect.getargspec(function).args

        assert calling_sequence[0] == 'self', "Wrong syntax for 'evaluate' in %s. The first argument " \
                                              "should be called 'self'." % name

        # Figure out how many variables are used

        variables = filter(lambda var: var in possible_variables, calling_sequence)

        # Check that they actually make sense. They must be used in the same order
        # as specified in possible_variables

        assert len(variables) > 0, "The name of the variables for 'evaluate' in %s must be one or more " \
                                   "among %s, instead of %s" % (name, ','.join(possible_variables), ",".join(variables))

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

        # Use unitless parameters when building the function, if no unit is specified, otherwise
        # use that unit
        if 'unit' not in definition or definition['unit'] is None or definition['unit'] == '':

            du = u.dimensionless_unscaled

        else:

            du = u.Unit(definition['unit'])

        def _parse_value(val):

            if val is None:

                return None

            if isinstance(val, str):

                return eval(val)

            else:

                return float(val)

        value = _parse_value(definition['initial value'])
        desc = definition['desc']

        # Optional attributes are either None if not specified, or the value specified

        min_value = (None if 'min' not in definition else _parse_value(definition['min']))
        max_value = (None if 'max' not in definition else _parse_value(definition['max']))
        delta = (None if 'delta' not in definition else _parse_value(definition['delta']))
        unit = du

        # A parameter can be fixed by using fix=yes, otherwise it is free by default

        free = (True if 'fix' not in definition else not bool(definition['fix']))

        new_parameter = Parameter(par_name, value, min_value=min_value, max_value=max_value,
                                  delta=delta, desc=desc, free=free, unit=unit)

        return new_parameter

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
        n_dim = type(instance)._n_dim

        if n_dim == 1:

            Function1D.__init__(instance,
                              type(instance)._name,
                              type(instance)._function_definition,
                              copy_of_parameters)

        elif n_dim == 2:

            Function2D.__init__(instance,
                                type(instance)._name,
                                type(instance)._function_definition,
                                copy_of_parameters)

        elif n_dim == 3:

            Function3D.__init__(instance,
                                type(instance)._name,
                                type(instance)._function_definition,
                                copy_of_parameters)

        # Last, if the class provides a setup method, call it
        if hasattr(instance, "_setup"):

            instance._setup()


class Function(Node):

    def __init__(self, name, function_definition, parameters):

        # (this is called by the constructor defined in the metaclass)

        # Store name, number of dimensions and the latex formula

        # Store also the function definition

        assert 'description' in function_definition,"Function definition must contain a description"

        if 'latex' not in function_definition:

            function_definition['latex'] = '$n.a.$'

        self._function_definition = function_definition

        # Set up the node

        Node.__init__(self, name)

        # Add the parameters as children. Since the name of the key in the dictionary might
        # be different than the actual name of the parameter, use the .add_child method instead
        # of the add_children method

        for child_name, child in parameters.iteritems():

            self._add_child(child, child_name)

        # Now generate a unique identifier (UUID) in a thread safe, multi-processing safe
        # way. This is used for example in the CompositeFunction class to keep track of the different
        # instances of the same function
        self._uuid = "{" + str(self._generate_uuid()) + "}"

        # Normal functions are able to handle units
        self._handle_units = True

    @property
    def n_dim(self):

        return type(self)._n_dim

    @property
    def free_parameters(self):
        """
        Returns a dictionary of free parameters for this function

        :return: dictionary of free parameters
        """

        free_parameters = collections.OrderedDict([(k,v) for k, v in self.parameters.iteritems() if v.free])

        return free_parameters

    @staticmethod
    def _generate_uuid():
        """
        Generate a unique identifier for this function.

        :return: the UUID
        """
        return uuid.UUID(bytes=os.urandom(16), version=4)

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
        return self._children

    @property
    def latex(self):
        """
        Returns the LaTEX formula for this function
        """
        return self._function_definition['latex']

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

        for parameter_name in self._children.keys():

            repr_dict['parameters'][parameter_name] = self._children[parameter_name].to_dict()

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

    def get_boundaries(self):
        """
        Returns the boundaries of this function. By default there is no boundary, but subclasses can
        override this.

        :return: a tuple of tuples containing the boundaries for each coordinate, or None if there are no boundaries
        """

        return None


class Function1D(Function):

    def __init__(self, name, function_definition, parameters):

        Function.__init__(self, name, function_definition, parameters)

        self._x_unit = None
        self._y_unit = None

    def get_wrapper(self):
        """
        Returns a python function which can be used to call this function with the parameters in the calling sequence.
        In other words, if you can call this function with f(x), with the returned wrapper you can call it with
        f(x, parameter1, parameter2...), where parameter1, parameter2... are the *free* parameters.

        :return: a python function
        """

        # Build a list of free parameters
        free_parameters = [k for k, v in self.parameters.iteritems() if v.free]

        # Prepare the method to set the parameters to their current value

        def set_parameters(*args):

            for i in range(len(args)):
                self._get_child(free_parameters[i]).value = args[i]

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

        # print(wrapper_code)

        exec (wrapper_code.replace("        ", "")) in locals()

        return wrapper

    def evaluate(self, x, *args, **kwargs):

        raise NotImplementedError("You have to re-implement this")

    def set_units(self, in_x_unit, in_y_unit):

        try:

            in_x_unit = u.Unit(in_x_unit)
            in_y_unit = u.Unit(in_y_unit)

        except:

            raise TypeError("Could not get a Unit instance from provided units when setting units "
                            "for function %s" % self.name)

        self._x_unit = in_x_unit
        self._y_unit = in_y_unit

        # Now call the underlying method to set units, which is defined by each function
        self._set_units(self._x_unit, self._y_unit)

    def _set_units(self, x_unit, y_unit):

        # This will be overridden by derived classes

        raise NotImplementedError("You have to implement the method _set_units for function %s" % self.name)

    @property
    def x_unit(self):
        return self._x_unit

    @property
    def y_unit(self):
        return self._y_unit

    def __call__(self, x, *args, **kwargs):

        #TODO args and kwargs are most probably useless here (check that!)

        # This method's code violates explicitly duck typing. The reason is that astropy.units introduce a very
        # significant overload on any computation. For this reason we treat differently the case with units from
        # the case without units, so that the latter case remains fast. Also, transforming an input
        # which is not an array into an array introduce a significant overload (10 microseconds or so), so we perform
        # this transformation only when strictly required

        if isinstance(x, np.ndarray):

            # We have an array as input

            if not isinstance(x, u.Quantity):

                # This is a normal array, let's use the fast call (without units)

                return self._call_without_units(x, *args, **kwargs)

            else:

                # This is an array with units, let's use the slow call which preserves units

                assert self.y_unit is not None, "In order to use units you need to use the function as a spectrum or " \
                                                "as something else," \
                                                "or you need to explicitly set the units."

                if self._handle_units:

                    results = self._call_with_units(x, *args, **kwargs)

                    # Now convert to the expected y unit
                    return results.to(self.y_unit)

                else:

                    # No support for units, add it artificially

                    new_input = np.array(x.value, dtype=float, ndmin=1, copy=False)

                    results = self._call_without_units(new_input, *args, **kwargs)

                    return np.squeeze(results) * self.y_unit

        else:

            # This is either a single number or a list
            if not isinstance(x, u.Quantity):

                # Transform the input to an array of floats

                new_input = np.array(x, dtype=float, ndmin=1, copy=False)

                # Compute the function

                result = self._call_without_units(new_input, *args, **kwargs)

                # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.

                return np.squeeze(result)

            else:

                if self._handle_units:

                    # This is a single number with units, let's transform it to an array with units

                    new_input = np.array(x, dtype=float, ndmin=1, copy=False) * x.unit

                    # Compute the function with units

                    result = self._call_with_units(new_input, *args, **kwargs)

                    # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.
                    # Let's also convert the result to the expected units

                    return np.squeeze(result).to(self.y_unit)

                else:

                    new_input = np.array(x.value, dtype=float, ndmin=1, copy=False)

                    # Compute the function without units

                    result = self._call_without_units(new_input, *args, **kwargs)

                    # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.
                    # Let's also convert the result to the expected units

                    return np.squeeze(result) * self.y_unit

    def _call_with_units(self, x, *args, **kwargs):

        # Gather the current parameters' values with units

        for parameter_name, parameter in self._children.iteritems():

            kwargs[parameter_name] = parameter.value * parameter.unit

        try:

            results = self.evaluate(x, *args, **kwargs)

        except u.UnitsError:

            raise

            raise u.UnitsError("Looks like you didn't provide all the units, or you provided the wrong ones, when "
                               "calling function %s" % self.name)

        else:

            return results

    def _call_without_units(self, x, *args, **kwargs):

        # Gather the current parameters' values without units, which means that the whole computation
        # will be without units, with a big speed gain (~10x)

        for parameter_name, parameter in self._children.iteritems():
            kwargs[parameter_name] = parameter.value

        return self.evaluate(x, *args, **kwargs)


class Function2D(Function):

    def __init__(self, name, function_definition, parameters):

        Function.__init__(self, name, function_definition, parameters)

        self._x_unit = None
        self._y_unit = None
        self._z_unit = None

    def evaluate(self, x, y, *args, **kwargs):

        raise NotImplementedError("You have to re-implement this")

    def set_units(self, in_x_unit, in_y_unit, in_z_unit):

        try:

            in_x_unit = u.Unit(in_x_unit)
            in_y_unit = u.Unit(in_y_unit)
            in_z_unit = u.Unit(in_z_unit)

        except:

            raise TypeError("Could not get a Unit instance from provided units when setting units "
                            "for function %s" % self.name)

        self._x_unit = in_x_unit
        self._y_unit = in_y_unit
        self._z_unit = in_z_unit

        # Now call the underlying method to set units, which is defined by each function
        self._set_units(self._x_unit, self._y_unit, self._z_unit)

    def _set_units(self, x_unit, y_unit, z_unit):

        # This will be overridden by derived classes

        raise NotImplementedError("You have to implement the method _set_units for function %s" % self.name)

    @property
    def x_unit(self):
        return self._x_unit

    @property
    def y_unit(self):
        return self._y_unit

    @property
    def z_unit(self):
        return self._z_unit

    def __call__(self, x, y, *args, **kwargs):

        # This method's code violates explicitly duck typing. The reason is that astropy.units introduce a very
        # significant overload on any computation. For this reason we treat differently the case with units from
        # the case without units, so that the latter case remains fast. Also, transforming an input
        # which is not an array into an array introduce a significant overload (10 microseconds or so), so we perform
        # this transformation only when strictly required

        assert type(x) == type(y), "You have to use the same type for x and y"

        if isinstance(x, np.ndarray):

            # We have an array as input

            if not isinstance(x, u.Quantity):

                # This is a normal array, let's use the fast call (without units)

                return self._call_without_units(x, y, *args, **kwargs)

            else:

                # This is an array with units, let's use the slow call which preserves units

                results = self._call_with_units(x, y, *args, **kwargs)

                # Now convert to the expected y unit
                return results.to(self.z_unit)

        else:

            # This is either a single number or a list
            if not isinstance(x, u.Quantity):

                # Transform the input to an array of floats

                new_x = np.array(x, dtype=float, ndmin=1, copy=False)
                new_y = np.array(y, dtype=float, ndmin=1, copy=False)

                # Compute the function

                result = self._call_without_units(new_x, new_y, *args, **kwargs)

                # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.

                return np.squeeze(result)

            else:

                # This is a single number with units, let's transform it to an array with units

                new_x = np.array(x, dtype=float, ndmin=1, copy=False) * x.unit
                new_y = np.array(y, dtype=float, ndmin=1, copy=False) * y.unit

                # Compute the function with units

                result = self._call_with_units(new_x, new_y, *args, **kwargs)

                # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.
                # Let's also convert the result to the expected units

                return np.squeeze(result).to(self.z_unit)

    def _call_with_units(self, x, y, *args, **kwargs):

        # Gather the current parameters' values with units

        for parameter_name, parameter in self._children.iteritems():
            kwargs[parameter_name] = parameter.value * parameter.unit

        try:

            results = self.evaluate(x, y, *args, **kwargs)

        except u.UnitsError:

            raise u.UnitsError("Looks like you didn't provide all the units, or you provided the wrong ones, when "
                               "calling function %s" % self.name)

        else:

            return results

    @memoize
    def _call_without_units(self, x, y, *args, **kwargs):

        # Gather the current parameters' values without units, which means that the whole computation
        # will be without units, with a big speed gain (~10x)

        for parameter_name, parameter in self._children.iteritems():
            kwargs[parameter_name] = parameter.value

        return self.evaluate(x, y, *args, **kwargs)


class Function3D(Function):

    def __init__(self, name, function_definition, parameters):

        Function.__init__(self, name, function_definition, parameters)

        self._x_unit = None
        self._y_unit = None
        self._z_unit = None
        self._w_unit = None

    def evaluate(self, x, y, z, *args, **kwargs):

        raise NotImplementedError("You have to re-implement this")

    def set_units(self, in_x_unit, in_y_unit, in_z_unit, in_w_unit):

        try:

            in_x_unit = u.Unit(in_x_unit)
            in_y_unit = u.Unit(in_y_unit)
            in_z_unit = u.Unit(in_z_unit)
            in_w_unit = u.Unit(in_w_unit)

        except:

            raise TypeError("Could not get a Unit instance from provided units when setting units "
                            "for function %s" % self.name)

        self._x_unit = in_x_unit
        self._y_unit = in_y_unit
        self._z_unit = in_z_unit
        self._w_unit = in_w_unit

        # Now call the underlying method to set units, which is defined by each function
        self._set_units(self._x_unit, self._y_unit, self._z_unit, self._w_unit)

    def _set_units(self, x_unit, y_unit, z_unit, w_unit):

        # This will be overridden by derived classes

        raise NotImplementedError("You have to implement the method _set_units for function %s" % self.name)

    @property
    def x_unit(self):
        return self._x_unit

    @property
    def y_unit(self):
        return self._y_unit

    @property
    def z_unit(self):
        return self._z_unit

    @property
    def w_unit(self):
        return self._w_unit

    def __call__(self, x, y, z, *args, **kwargs):

        # This method's code violates explicitly duck typing. The reason is that astropy.units introduce a very
        # significant overload on any computation. For this reason we treat differently the case with units from
        # the case without units, so that the latter case remains fast. Also, transforming an input
        # which is not an array into an array introduce a significant overload (10 microseconds or so), so we perform
        # this transformation only when strictly required

        assert type(x) == type(y) and type(y) == type(z), "You have to use the same type for x, y and z"

        if isinstance(x, np.ndarray) and x.shape != ():

            # We have an array as input

            if not isinstance(x, u.Quantity):

                # This is a normal array, let's use the fast call (without units)

                return self._call_without_units(x, y, z, *args, **kwargs)

            else:

                # This is an array with units, let's use the slow call which preserves units

                results = self._call_with_units(x, y, z, *args, **kwargs)

                # Now convert to the expected y unit
                return results.to(self.w_unit)

        else:

            # This is either a single number or a list
            if not isinstance(x, u.Quantity):

                # Transform the input to an array of floats

                new_x = np.array(x, dtype=float, ndmin=1, copy=False)
                new_y = np.array(y, dtype=float, ndmin=1, copy=False)
                new_z = np.array(z, dtype=float, ndmin=1, copy=False)

                # Compute the function

                result = self._call_without_units(new_x, new_y, new_z, *args, **kwargs)

                # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.

                return np.squeeze(result)

            else:

                # This is a single number with units, let's transform it to an array with units

                new_x = np.array(x, dtype=float, ndmin=1, copy=False) * x.unit
                new_y = np.array(y, dtype=float, ndmin=1, copy=False) * y.unit
                new_z = np.array(z, dtype=float, ndmin=1, copy=False) * z.unit

                # Compute the function with units

                result = self._call_with_units(new_x, new_y, new_z, *args, **kwargs)

                # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.
                # Let's also convert the result to the expected units

                return np.squeeze(result).to(self.w_unit)

    def _call_with_units(self, x, y, z, *args, **kwargs):

        # Gather the current parameters' values with units

        for parameter_name, parameter in self._children.iteritems():

            kwargs[parameter_name] = parameter.value * parameter.unit

        try:

            results = self.evaluate(x, y, z, *args, **kwargs)

        except u.UnitsError:

            raise u.UnitsError("Looks like you didn't provide all the units, or you provided the wrong ones, when "
                               "calling function %s" % self.name)

        else:

            return results

    def _call_without_units(self, x, y, z, *args, **kwargs):

        # Gather the current parameters' values without units, which means that the whole computation
        # will be without units, with a big speed gain (~10x)

        for parameter_name, parameter in self._children.iteritems():

            kwargs[parameter_name] = parameter.value

        return self.evaluate(x, y, z, *args, **kwargs)


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
                parameter._change_name(new_name)

        # Now build a meaningful description

        _function_definition = {'description': self.expression, 'latex': NO_LATEX_FORMULA}

        Function.__init__(self, 'composite', _function_definition, parameters)

        self._uuid = self._uuid_expression

    def set_units(self, x_unit, y_unit):

        # Just rely on the single functions to adjust themselves.

        for function in self.functions:

            function.set_units(x_unit, y_unit)

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
    Returns the function "name", which must be among the known functions or a composite function.

    :param function_name: the name of the function (use 'composite' if the function is a composite function)
    :param composite_function_expression: composite function specification such as
    ((((powerlaw{1} + (sin{2} * 3)) + (sin{2} * 25)) - (powerlaw{1} * 16)) + (sin{2} ** 3.0))
    :return: the an instance of the requested class

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

def get_function_class(function_name):
    """
    Return the type for the requested function

    :param function_name: the function to return
    :return: the type for that function (i.e., this is a class, not an instance)
    """

    if function_name in _known_functions:

        return _known_functions[function_name]

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

    unique_functions = set(re.findall(r'\b([a-zA-Z0-9_]+)\{([0-9]?)\}',function_specification))

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

    #print(string_for_literal_eval)

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