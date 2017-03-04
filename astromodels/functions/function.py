import ast
import collections
import copy
import inspect
import uuid
from operator import attrgetter

import astropy.units as u
import numpy as np
import os
import re
from yaml.reader import ReaderError

from astromodels.core.my_yaml import my_yaml
from astromodels.core.parameter import Parameter
from astromodels.core.tree import Node
from astromodels.utils.pretty_list import dict_to_list
from astromodels.utils.table import dict_to_table
from astromodels.core.memoization import memoize

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



# This dictionary will contain the known function by name, so that the model_parser can instance
# them by looking into this dictionary. It will be filled by the FunctionMeta meta-class.

_known_functions = {}


# The following is a metaclass for all the functions
class FunctionMeta(type):
    """
    A metaclass for the models, which takes care of setting up the parameters and the other attributes
    according to the definition given in the documentation of the function class.
    """

    def __new__(mcs, name, bases, dct):

        # We do the parsing of the parameters in the __new__ instead of the __init__ so this is the first thing
        # that runs when importing astromodels

        # Enforce the presence of the evaluate method

        if 'evaluate' not in dct:

            raise AttributeError("You have to implement the 'evaluate' method in %s" % name)

        # We also need the method _set_units

        if '_set_units' not in dct:
            raise AttributeError("You have to implement the '_set_units' method in %s" % name)

        # Now parse the documentation of the function which contains the parameter specification

        # The doc is a YAML document containing among other things the definition of the parameters

        # Parse it

        try:

            function_definition = my_yaml.load(dct['__doc__'])

        except ReaderError:  # pragma: no cover

            raise DocstringIsNotRaw("Docstring parsing has failed. "
                                    "Did you remember to specify the docstring of %s as raw? "
                                    "To do that, you have to put a r before the docstring, "
                                    '''like in \n\nr"""\n(docstring)\n"""\n\ninstead of just\n\n'''
                                    '''"""\ndocstring\n"""''' % name)

        else:

            # Store the function definition in the type

            dct['_function_definition'] = function_definition

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
        dct['_latex'] = latex_formula

        # Parse the parameters' dictionary
        assert isinstance(function_definition['parameters'], dict), "Wrong syntax in 'parameters' token. It must be " \
                                                                    "a dictionary. Refer to the documentation."

        # Add the parameters as attribute of the *type*. During the __call__ method below this dictionary will be used
        # to create a copy of each parameter which will be made available as child of the *instance*.

        dct['_parameters'] = collections.OrderedDict()

        for parameter_name, parameter_definition in function_definition['parameters'].iteritems():

            this_parameter = FunctionMeta.parse_parameter_definition(name, parameter_name, parameter_definition)

            dct['_parameters'][this_parameter.name] = this_parameter

        # Now perform a minimal check of the 'evaluate' function

        variables, parameters_in_calling_sequence = FunctionMeta.check_calling_sequence(name, 'evaluate',
                                                                                        dct['evaluate'],
                                                                                        ['x', 'y', 'z'])

        # Now check that all the parameters used in 'evaluate' are part of the documentation,
        # and that there are no unused parameters

        set1 = set(dct['_parameters'].keys())
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

        # Figure out the dimensionality of this function

        n_dim = len(variables)

        # Store the dimensionality in the *type*

        dct['_n_dim'] = n_dim

        # Now add the constructor to the class, if it does not provide one
        # You shouldn't usually provide a constructor, that's only for advance uses
        # like the TemplateModel

        if '_custom_init_' in dct:

            dct['__init__'] = dct['_custom_init_']

        else:

            dct['__init__'] = FunctionMeta.class_init

        # Finally, add the info() method to the type so that it can be called even without instancing the class

        def info():

            repr_dict = collections.OrderedDict()

            repr_dict['description'] = function_definition['description']

            if 'latex' in function_definition:
                repr_dict['formula'] = function_definition['latex']

            # Add the description of each parameter and their current value
            repr_dict['default parameters'] = collections.OrderedDict()

            for parameter_name in dct['_parameters'].keys():

                repr_dict['default parameters'][parameter_name] = dct['_parameters'][parameter_name].to_dict()

            if has_ipython:

                display(HTML(dict_to_list(repr_dict, html=True)))

            else:

                print(dict_to_list(repr_dict, html=False))

        dct['info'] = staticmethod(info)

        # Now call the __new__ of the "type" class (which then will call the __init__ of this metaclass)

        return super(FunctionMeta, mcs).__new__(mcs, name, bases, dct)

    def __init__(cls, name, bases, dct):

        # This is the MetaClass init, which is called after the __new__ is done

        # Store the name of the function in the type
        cls._name = cls.__name__

        # Add this as a known function

        _known_functions[name] = cls

        # Finally call the init of the type class

        super(FunctionMeta, cls).__init__(name, bases, dct)

    @staticmethod
    def class_init(instance, **kwargs):

        # This is what is going to be called as the __init__ of the class, every time a new instance
        # is created

        # Create a copy of the parameters dictionary which is in the type,
        # otherwise every instance would share the same dictionary

        copy_of_parameters = collections.OrderedDict()

        # Fill it by duplicating the parameters contained in the dictionary in the type

        for key, value in type(instance)._parameters.iteritems():

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

        # Now call the init of the corresponding class
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

            if isinstance(val, str):

                return eval(val)

            elif val is None:

                return None

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


class Function(Node):
    """
    Generic Function class. Will be subclassed in Function1D, Function2D and Function3D.

    """
    def __init__(self, name=None, function_definition=None, parameters=None):

        # I use default values only to avoid warnings from pycharm and other software about the
        # calling sequence of this contructor. We actually need to enforce its proper use,
        # with this assert

        assert name is not None and function_definition is not None and parameters is not None

        # Set up the node

        super(Function, self).__init__(name)

        # Store name, number of dimensions and the latex formula

        # Store also the function definition

        assert 'description' in function_definition,"Function definition must contain a description"

        if 'latex' not in function_definition:

            function_definition['latex'] = '$n.a.$'

        self._function_definition = function_definition

        # Add the parameters as children. Since the name of the key in the dictionary might
        # be different than the actual name of the parameter, use the .add_child method instead
        # of the add_children method

        self._parameters = collections.OrderedDict()

        for child_name, child in parameters.iteritems():

            self._parameters[child_name] = child

            # Add the parameter as a child of the function

            self._add_child(child)

        # Now generate a unique identifier (UUID) in a thread safe, multi-processing safe
        # way. This is used for example in the CompositeFunction class to keep track of the different
        # instances of the same function
        self._uuid = "{" + str(self._generate_uuid()) + "}"

        # Normal functions are able to change units, while some specific ones (such as the one from XSpec) are not.
        # In this second case, this variable will contain a tuple (x_unit, y_unit), but by default it should be
        # None

        self._fixed_units = None

    @property
    def n_dim(self):
        """
        :return: number of dimensions for this function (1, 2 or 3)
        """
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

    def has_fixed_units(self):
        """
        Returns True if this function cannot change units, which is the case only for very specific functions (like
        models from foreign libraries like Xspec)

        :return: True or False
        """

        return not (self._fixed_units is None)

    @property
    def fixed_units(self):
        """
            Returns the fixed units if has_fixed_units is True (see has_fixed_units)

            :return: None, or a tuple (x_unit, y_unit)
        """

        return self._fixed_units

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
        return self._parameters

    @property
    def latex(self):
        """
        Returns the LaTEX formula for this function
        """
        return self._function_definition['latex']

    # Define now all the operators which allow to combine functions. Each operator will return a new
    # instance of a CompositeFunction, which can then be used as a function on its own

    def of(self, another_function):
        """
        Compose this function with another as in this_function(another_function(x))
        :param another_function: another function to compose with the current one
        :return: a composite function instance
        """
        return CompositeFunction('of', self, another_function)

    def __neg__(self):

        return CompositeFunction('*-', self)

    def __abs__(self):

        return CompositeFunction('abs', self)

    def __pow__(self, other_instance):

        return CompositeFunction('**', self, other_instance)

    def __rpow__(self, other_instance):

        return CompositeFunction('**', other_instance, self)

    def __add__(self, other_instance):

        return CompositeFunction('+', self, other_instance)

    __radd__ = __add__

    def __sub__(self, other_instance):

        return CompositeFunction('-', self, other_instance)

    __rsub__ = __sub__

    def __mul__(self, other_instance):

        c = CompositeFunction('*', self, other_instance)

        # If the other instance is a function (and not a number), flag it so its units will be made dimensionless
        # in the set_units method of the composite function

        if isinstance(other_instance, Function):

            if self.has_fixed_units():

                # This is likely a XSpec model. The multiplication of two models with units will give the wrong
                # units to the results. So, depending on the type of the first and second model, we need to adjust
                # their units so that the result will keep the right units.

                if not self.fixed_units[1] == u.dimensionless_unscaled:

                    # This function has fixed unit and is not dimensionless (likely an additive XSpec model).
                    # We need to make the other function dimensionless so that the multiplication of them will keep
                    # the right units

                    other_instance._make_dimensionless = True

                else:

                    # This function has fixed unit, but it is dimensionless (likely a multiplicative XSpec model)
                    # The other function should keep its units, so we flag self instead

                    self._make_dimensionless = True

                    # However, if also the other function has fixed dimension and it is dimensionless
                    # (likely another XSpec multiplicative model) we need to flag that as well
                    if other_instance.has_fixed_units() and (other_instance.fixed_units[1] == u.dimensionless_unscaled):
                        other_instance._make_dimensionless = True

                        # Now we need to flag the composite function as fixed units and dimensionless
                        c._fixed_units = other_instance.fixed_units

            else:

                # We need to make the other instance dimensionless so that this function (which is not dimensionless)
                # multiplied by the other function (which we will make dimensionless) will give the right units as
                # the results

                other_instance._make_dimensionless = True

        return c

    __rmul__ = __mul__

    def __div__(self, other_instance):

        c = CompositeFunction('/', self, other_instance)

        # If the other instance is a function (and not a number), flag it so its units will be made dimensionless
        # in the set_units method of the composite function

        if isinstance(other_instance, Function):

            if self.has_fixed_units():

                # This is likely a XSpec model. Division is not supported.

                raise NotImplementedError("Division for XSpec models is not supported")

            else:

                # We need to make the other instance dimensionless

                other_instance._make_dimensionless = True

        return c

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

        for parameter in self._get_children():

            repr_dict['parameters'][parameter.name] = parameter.to_dict()

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

    def get_boundaries(self):  # pragma: no cover
        """
        Returns the boundaries of this function. By default there is no boundary, but subclasses can
        override this.

        :return: a tuple of tuples containing the boundaries for each coordinate (ra_min, ra_max), (dec_min, dec_max)
        """

        raise NotImplementedError("You have to implement this")

    def __call__(self, *args):  # pragma: no cover

        raise NotImplementedError("You have to implement this")

    def fast_call(self, *args):  # pragma: no cover

        raise NotImplementedError("You have to implement this")

    def evaluate_at(self, *args, **parameter_specification):  # pragma: no cover
        """
        Evaluate the function at the given x(,y,z) for the provided parameters, explicitly provided as part of the
        parameter_specification keywords.

        :param *args:
        :param **parameter_specification:
        :return:
        """

        # Set the parameters to the provided values
        for parameter in parameter_specification:

            self._get_child(parameter).value = parameter_specification[parameter]

        return self(*args)


class Function1D(Function):

    def __init__(self, name=None, function_definition=None, parameters=None):

        Function.__init__(self, name, function_definition, parameters)

        self._x_unit = None
        self._y_unit = None

    def evaluate(self, x, *args, **kwargs):  # pragma: no cover

        raise NotImplementedError("You have to re-implement this")

    def set_units(self, in_x_unit, in_y_unit):

        # Transform None in input to '', so that u.Unit() will generate a dimensionless unit

        in_x_unit = in_x_unit if in_x_unit is not None else ''
        in_y_unit = in_y_unit if in_y_unit is not None else ''

        # Get a Unit instance from the inputs

        try:

            in_x_unit = u.Unit(in_x_unit)
            in_y_unit = u.Unit(in_y_unit)

        except:

            raise TypeError("Could not get a Unit instance from provided units %s when setting units "
                            "for function %s" % ((in_x_unit, in_y_unit), self.name))

        # Now call the underlying method to set units, which is defined by each function
        new_units = self._set_units(in_x_unit, in_y_unit)

        # Store the units.
        # NOTE: the previous call to _set_units might return new units in special cases
        # (for example Xspec functions). So it is critical that we store them in the class' attributes
        # after the call to _set_units

        if new_units is not None:

            new_x_unit, new_y_unit = new_units

            self._x_unit = new_x_unit
            self._y_unit = new_y_unit

        else:

            self._x_unit = in_x_unit
            self._y_unit = in_y_unit

    def _set_units(self, x_unit, y_unit):  # pragma: no cover

        # This will be overridden by derived classes

        raise NotImplementedError("You have to implement the method _set_units for function %s" % self.name)

    @property
    def x_unit(self):
        """
        The unit of the independent variable
        :return: a astropy.Unit instance
        """
        return self._x_unit

    @property
    def y_unit(self):
        """
        The unit of the dependent variable
        :return: a astropy.Unit instance
        """
        return self._y_unit

    def __call__(self, x):

        # This method's code violates explicitly duck typing. The reason is that astropy.units introduce a very
        # significant overload on any computation. For this reason we treat differently the case with units from
        # the case without units, so that the latter case remains fast. Also, transforming an input
        # which is not an array into an array introduce a significant overload (10 microseconds or so), so we perform
        # this transformation only when strictly required

        # NOTE: for a single quantity such as q = (1.0 * u.keV), isinstance(q, np.ndarray) returns True

        if isinstance(x, np.ndarray):

            # We have an array as input (or a single quantity)

            if not isinstance(x, u.Quantity):

                # This is a normal array, let's use the fast call (without units)

                return self.fast_call(x)

            else:

                # This is an array with units or a single quantity, let's use the slow call which preserves units

                assert self.y_unit is not None, "In order to use units you need to use the function as a spectrum or " \
                                                "as something else, or you need to explicitly set the units."

                results = self._call_with_units(x)

                # Now convert to the expected y unit and return a astropy.Quantity by multiplying by the right unit
                return np.squeeze(results.to(self.y_unit).value) * self.y_unit

        else:

            # This is either a single number or a list

            # Transform the input to an array of floats. If x is a single number, this will be an array of size 1

            new_input = np.array(x, dtype=float, ndmin=1, copy=False)

            # Compute the function

            result = self.fast_call(new_input)

            # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number,
            # so that if the input was a single number, also the output will be a single number

            return np.squeeze(result)

    def _call_with_units(self, x):

        # Gather the current parameters' values with units
        values = map(attrgetter("as_quantity"), self._get_children())

        try:

            results = self.evaluate(x, *values)

        except u.UnitsError:  # pragma: no cover

            raise u.UnitsError("Looks like you didn't provide all the units, or you provided the wrong ones, when "
                               "calling function %s" % self.name)

        else:

            return results

    @memoize
    def fast_call(self, x):

        # Gather the current parameters' values without units, which means that the whole computation
        # will be without units, with a big speed gain (~10x)

        # NOTE: it is important to use value, and not _value, to support linking

        values = map(attrgetter("value"), self._get_children())

        return self.evaluate(x, *values)

    def get_boundaries(self):
        """
        Returns the boundaries of this function. By default there is no boundary, but subclasses can
        override this.

        :return: a tuple of tuples containing the boundaries for each coordinate (ra_min, ra_max), (dec_min, dec_max)
        """

        raise DesignViolation("Cannot call get_boundaries() on a 1d function")


class Function2D(Function):

    def __init__(self, name=None, function_definition=None, parameters=None):

        Function.__init__(self, name, function_definition, parameters)

        self._x_unit = None
        self._y_unit = None
        self._z_unit = None

    def evaluate(self, x, y, *args):  # pragma: no cover

        raise NotImplementedError("You have to re-implement this")

    def set_units(self, in_x_unit, in_y_unit, in_z_unit):

        # Change None to '' for the inputs so that the following u.Unit construction will generate a dimensionless
        # unit (it would fail with None)

        in_x_unit = in_x_unit if in_x_unit is not None else ''
        in_y_unit = in_y_unit if in_y_unit is not None else ''
        in_z_unit = in_z_unit if in_z_unit is not None else ''

        try:

            in_x_unit = u.Unit(in_x_unit)
            in_y_unit = u.Unit(in_y_unit)
            in_z_unit = u.Unit(in_z_unit)

        except:

            raise TypeError("Could not get a Unit instance from provided units when setting units "
                            "for function %s" % self.name)

        # Store the Unit instances

        self._x_unit = in_x_unit
        self._y_unit = in_y_unit
        self._z_unit = in_z_unit

        # Now call the underlying method to set units, which is defined by each function
        self._set_units(self._x_unit, self._y_unit, self._z_unit)

    def _set_units(self, x_unit, y_unit, z_unit):  # pragma: no cover

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

            # We have an array or a quantity as input

            if not isinstance(x, u.Quantity):

                # This is a normal array, let's use the fast call (without units)

                return self._call_without_units(x, y)

            else:

                # This is an array with units or a single quantity, let's use the slow call which preserves units

                results = self._call_with_units(x, y)

                # Now convert to the expected z unit and remove useless dimensions
                return np.squeeze(results.to(self.z_unit).value) * self.z_unit

        else:

            # This is either a single number or a list

            # Transform the input to an array of floats

            new_x = np.array(x, dtype=float, ndmin=1, copy=False)
            new_y = np.array(y, dtype=float, ndmin=1, copy=False)

            # Compute the function

            result = self._call_without_units(new_x, new_y)

            # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.

            return np.squeeze(result)

    def _call_with_units(self, x, y):

        # Gather the current parameters' values with units

        values = map(attrgetter("as_quantity"), self._get_children())

        try:

            results = self.evaluate(x, y, *values)

        except u.UnitsError:  # pragma: no cover

            raise u.UnitsError("Looks like you didn't provide all the units, or you provided the wrong ones, when "
                               "calling function %s" % self.name)

        else:

            return results

    @memoize
    def _call_without_units(self, x, y):

        # Gather the current parameters' values without units, which means that the whole computation
        # will be without units, with a big speed gain (~10x)

        values = map(attrgetter("value"), self._get_children())

        return self.evaluate(x, y, *values)


class Function3D(Function):

    def __init__(self, name=None, function_definition=None, parameters=None):

        Function.__init__(self, name, function_definition, parameters)

        self._x_unit = None
        self._y_unit = None
        self._z_unit = None
        self._w_unit = None

    def evaluate(self, x, y, z, *args, **kwargs):  # pragma: no cover

        raise NotImplementedError("You have to re-implement this")

    def set_units(self, in_x_unit, in_y_unit, in_z_unit, in_w_unit):

        # Change None to '' for the inputs so that the following u.Unit construction will generate a dimensionless
        # unit (it would fail with None)

        in_x_unit = in_x_unit if in_x_unit is not None else ''
        in_y_unit = in_y_unit if in_y_unit is not None else ''
        in_z_unit = in_z_unit if in_z_unit is not None else ''
        in_w_unit = in_w_unit if in_w_unit is not None else ''

        # Get instances of Unit

        try:

            in_x_unit = u.Unit(in_x_unit)
            in_y_unit = u.Unit(in_y_unit)
            in_z_unit = u.Unit(in_z_unit)
            in_w_unit = u.Unit(in_w_unit)

        except:

            raise TypeError("Could not get a Unit instance from provided units when setting units "
                            "for function %s" % self.name)

        # Store the Unit instances

        self._x_unit = in_x_unit
        self._y_unit = in_y_unit
        self._z_unit = in_z_unit
        self._w_unit = in_w_unit

        # Now call the underlying method to set units, which is defined by each function

        self._set_units(self._x_unit, self._y_unit, self._z_unit, self._w_unit)

    def _set_units(self, x_unit, y_unit, z_unit, w_unit):  # pragma: no cover

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

    def __call__(self, x, y, z):

        # This method's code violates explicitly duck typing. The reason is that astropy.units introduce a very
        # significant overload on any computation. For this reason we treat differently the case with units from
        # the case without units, so that the latter case remains fast. Also, transforming an input
        # which is not an array into an array introduce a significant overload (10 microseconds or so), so we perform
        # this transformation only when strictly required

        assert type(x) == type(y) and type(y) == type(z), "You have to use the same type for x, y and z"

        if isinstance(x, np.ndarray):

            # We have an array as input

            if not isinstance(x, u.Quantity):

                # This is a normal array, let's use the fast call (without units)

                return self._call_without_units(x, y, z)

            else:

                # This is an array with units or a single quantity, let's use the slow call which preserves units

                results = self._call_with_units(x, y, z)

                # Now convert to the expected w unit and remove useless dimensions

                return np.squeeze(results.to(self.w_unit).value) * self.w_unit

        else:

            # This is either a single number or a list
            # Transform the input to an array of floats

            new_x = np.array(x, dtype=float, ndmin=1, copy=False)
            new_y = np.array(y, dtype=float, ndmin=1, copy=False)
            new_z = np.array(z, dtype=float, ndmin=1, copy=False)

            # Compute the function

            result = self._call_without_units(new_x, new_y, new_z)

            # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.

            return np.squeeze(result)

    def _call_with_units(self, x, y, z):

        # Gather the current parameters' values with units

        values = map(attrgetter("as_quantity"), self._get_children())

        try:

            results = self.evaluate(x, y, z, *values)

        except u.UnitsError:  # pragma: no cover

            raise u.UnitsError("Looks like you didn't provide all the units, or you provided the wrong ones, when "
                               "calling function %s" % self.name)

        else:

            return results

    @memoize
    def _call_without_units(self, x, y, z):

        # Gather the current parameters' values without units, which means that the whole computation
        # will be without units, with a big speed gain (~10x)

        values = map(attrgetter("value"), self._get_children())

        return self.evaluate(x, y, z, *values)


##########################
# Composite function stuff
##########################

# Codes to indicate to Composite Function the operation between two functions
_operations = {'+': np.add,
               '-': np.subtract,
               '*-': np.negative,
               '*': np.multiply,
               '/': np.divide,
               '**': np.power,
               'abs': np.abs,
               'of': 'compose'}


# These methods need to be here to overcome the limitation of pickle with methods of classes

def _cf_evaluate_func_func(np_operator, f1, f2, *args):
    # Evaluate for when both elements are functions

    return np_operator(f1(*args), f2(*args))


def _cf_evaluate_func_number(np_operator, f1, f2, *args):
    # Evaluate for when element 1 is a function and element 2 is a number

    return np_operator(f1(*args), f2)


def _cf_evaluate_number_func(np_operator, f1, f2, *args):
    # Evaluate for when element 2 is a function and element 1 is a number

    return np_operator(f1, f2(*args))


def _cf_evaluate_func_of_func(np_operator, f1, f2, *args):
    value = f2(*args)

    return f1(value, *(args[1:]))


class CompositeFunction(Function):

    def __init__(self, operation, function_or_scalar_1, function_or_scalar_2=None):

        assert operation in _operations, "Do not know operation %s" % operation

        # Save this to make the class pickeable (see the __setstate__ and __getstate__ methods)

        self._calling_sequence = (operation, function_or_scalar_1, function_or_scalar_2)

        self._requested_x_unit = None
        self._requested_y_unit = None

        self._operation = operation

        # Set the new __call__ according to the type of the elements in the expression
        self._decide_evaluate_type()

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

                match = re.match("(.+)_[0-9]+$", parameter_name)

                if match is not None:

                    original_name = match.groups()[0]

                else:

                    original_name = parameter_name

                new_name = "%s_%i" % (original_name, i+1)

                # Store the parameter under the new name (obviously this is a reference to the
                # parameter, not a copy, as always in python)

                parameters[new_name] = parameter
                parameter._change_name(new_name)

        # Now build a meaningful description

        _function_definition = {'description': self.expression, 'latex': NO_LATEX_FORMULA}

        Function.__init__(self, 'composite', _function_definition, parameters)

        self._uuid = self._uuid_expression

    def set_units(self, x_unit, y_unit, relaxed=False):

        if relaxed and (x_unit is None) and (y_unit is None):

            # This can happen when rebuilding a composite function during unpickling, when
            # there are more than two functions composed together. We do not need to to anything in that case
            pass

        else:

            self._requested_x_unit = x_unit
            self._requested_y_unit = y_unit

            # Just rely on the single functions to adjust themselves.

            for function in self.functions:

                if hasattr(function, '_make_dimensionless'):

                    function.set_units(x_unit, u.dimensionless_unscaled)

                else:

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

    def _decide_evaluate_type(self):

        # Assign to __call__ the right evaluate according to the type of the elements in the expression

        operation, self._f1, self._f2 = self._calling_sequence

        np_operator = _operations[operation]

        self._np_operator = np_operator

        if np_operator == "compose":

            assert hasattr(self._f2, 'evaluate'), "Second member of .of cannot be a scalar"

            assert self._f1.n_dim == 1 and self._f2.n_dim == 1, "Can only compose with .of functions of 1 variable"

            self.evaluate = _cf_evaluate_func_of_func

        else:

            # Check whether the second member is a function, or a number

            if hasattr(self._f1, 'evaluate'):

                if hasattr(self._f2, 'evaluate'):

                    self.evaluate = _cf_evaluate_func_func

                else:

                    self.evaluate = _cf_evaluate_func_number

            else:

                if hasattr(self._f2, 'evaluate'):

                    self.evaluate = _cf_evaluate_number_func

                else:  # pragma: no cover

                    # Should never get here!

                    raise RuntimeError("Should never get here")

    @property
    def functions(self):
        "A list containing the function used to build this composite function"
        return self._functions

    def evaluate(self):  # pragma: no cover

        raise NotImplementedError("You cannot instance and use a composite function by itself. Use the factories.")

    # This dumb function must be here because it is not possible to override at runtime __call__ (nor any other
    # special method)
    def __call__(self, x):

        return self.evaluate(self._np_operator, self._f1, self._f2, x)

    # For composite function, fast_call is the same as __call__ (because the call will be forwarded to the
    # inner functions)

    fast_call = __call__

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

            # Maybe this is a template

            # NOTE: import here to avoid circular import

            from astromodels.functions.template_model import TemplateModel, MissingDataFile

            try:

                instance = TemplateModel(function_name)

            except MissingDataFile:

                raise UnknownFunction("Function %s is not known. Known functions are: %s" %
                                      (function_name, ",".join(_known_functions.keys())))

            else:

                return instance


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

            # It might be a template

            # This import is here to avoid circular dependency between this module and TemplateModel.py
            import astromodels.functions.template_model

            try:

                instance = astromodels.functions.template_model.TemplateModel(unique_function)

            except astromodels.functions.template_model.MissingDataFile:

                # It's not a template

                raise UnknownFunction("Function %s in expression %s is unknown. If this is a template model, you are "
                                      "probably missing the data file" % (unique_function, function_specification))

            else:

                # It's a template

                instances[complete_function_specification] = instance

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