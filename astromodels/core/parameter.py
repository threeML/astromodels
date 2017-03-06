__author__ = 'giacomov'

__doc__ = """"""

import collections
import copy
import exceptions

import astropy.units as u
import numpy as np
import scipy.stats
import warnings

from astromodels.core.tree import Node


def _behaves_like_a_number(obj):
    """

    :param obj:
    :return: True if obj supports addition, subtraction, multiplication and division. False otherwise.
    """

    try:

        obj + 1
        obj * 2
        obj / 2
        obj - 1

    except TypeError:

        return False

    else:

        return True


# Exception for when a parameter is out of its bounds
class SettingOutOfBounds(RuntimeError):
    pass


# Exception for when an object should be callable but it isn't
class NotCallableOrErrorInCall(RuntimeError):
    pass


class WarningUnitsAreSlow(Warning):
    pass


class IndependentVariableCannotBeLinked(RuntimeError):
    pass


class CannotUnderstandUnit(RuntimeError):
    pass


class CannotConvertValueToNewUnits(RuntimeError):
    pass


class ParameterMustHaveBounds(RuntimeError):
    pass


def accept_quantity(input_type=float, allow_none=False):
    """
        A class-method decorator which allow a given method (typically the set_value method) to receive both a
        astropy.Quantity or a simple float, but to be coded like it's always receiving a pure float in the right units.
        This is to give a way to avoid the huge bottleneck that are astropy.units

        :param input_type: the expected type for the input (float, int)
        :param allow_none : whether to allow or not the passage of None as argument (default: False)
        :return: a decorator for the particular type
        """

    def accept_quantity_wrapper(method):

        def handle_quantity(instance, value, *args, **kwargs):

            # For speed reasons, first run the case where the input is not a quantity, and fall back to the handling
            # of quantities if that fails. The parts that fails if the input is a Quantity is the conversion
            # input_type(value). This could have been handled more elegantly with a "finally" clause, but that would
            # have a 40 percent speed impact...

            try:

                new_value = input_type(value)

                return method(instance, new_value, *args, **kwargs)

            except TypeError:

                # Slow for slow, check that we actually have a quantity or None (if allowed)

                if isinstance(value, u.Quantity):

                    new_value = value.to(instance.unit).value

                    return method(instance, new_value, *args, **kwargs)

                elif value is None:

                    if allow_none:

                        return method(instance, None, *args, **kwargs)

                    else: # pragma: no cover

                        raise TypeError("You cannot pass None as argument for "
                                        "method %s of %s" % (method.__name__, instance.name))

                else: # pragma: no cover

                    raise TypeError("You need to pass either a %s or a astropy.Quantity "
                                    "to method %s of %s" % (input_type.__name__, method.__name__, instance.name))

        return handle_quantity

    return accept_quantity_wrapper


class ParameterBase(Node):

    def __init__(self, name, value, min_value=None, max_value=None, desc=None, unit=u.dimensionless_unscaled):

        # Make this a node

        Node.__init__(self, name)

        # Make a static name which will never change (not even after a _change_name call)
        self._static_name = str(name)

        # Callbacks are executed any time the value for the parameter changes (i.e., its value changes)

        # We start from a empty list of callbacks.
        self._callbacks = []

        # Assign to members

        # Store the units as an astropy.units.Unit instance

        self._unit = u.Unit(unit)

        # Let's store the init value

        # If the value is a Quantity, deal with that

        if isinstance(value, u.Quantity):

            # If the user did not specify an ad-hoc unit, use the unit
            # of the Quantity

            if self._unit == u.dimensionless_unscaled:

                self._unit = value.unit

            # Convert the value to the provided unit (if necessary)

            self._value = value.to(self._unit).value

        else:

            self._value = value

        # Set minimum if provided, otherwise use default
        # (use the property so the checks that are there are performed also on construction)

        self._min_value = None # this will be overwritten immediately in the next line
        self.min_value = min_value

        # Set maximum if provided, otherwise use default

        self._max_value = None  # this will be overwritten immediately in the next line
        self.max_value = max_value

        # Store description

        self._desc = desc

        # Make the description the documentation as well

        self.__doc__ = desc

        # Now perform a very lazy check that we can perform math operations on value, minimum and maximum
        # (i.e., they are numbers)

        if not _behaves_like_a_number(self._value):

            raise TypeError("The provided initial value is not a number")

    def _repr__base(self, rich_output): # pragma: no cover

        raise NotImplementedError("You need to implement this for the actual Parameter class")

    # Define the property 'description' and make it read-only

    @property
    def static_name(self):
        """
        Returns a name which will never change, even if the name of the parameter does (for example in composite
        functions)

        :return : the static name
        :type : str
        """

        return self._static_name

    @property
    def description(self):
        """
        Return a description of this parameter

        :return: a string cointaining a description of the meaning of this parameter
        """
        return self._desc

    # Define the property 'unit'
    def _set_unit(self, input_unit):

        # This will fail if the input is not valid

        new_unit = u.Unit(input_unit)

        # Now transform the current _value in the new unit, unless the current unit is dimensionless, in which
        # case there is no transformation to make

        # (self._unit is the OLD unit here)

        if self._unit != u.dimensionless_unscaled:

            # This will fail if the new unit is not compatible with the old one

            try:

                self._value = (self._value * self._unit).to(new_unit).value

            except u.UnitConversionError:

                if new_unit == u.dimensionless_unscaled:

                    new_unit_name = '(dimensionless)'

                else:

                    new_unit_name = new_unit

                raise CannotConvertValueToNewUnits("Cannot convert the value %s from %s to the "
                                                   "new units %s" % (self._value, self._unit, new_unit_name))

        else:

            # No need to transform the value
            pass

        # Finally store the new unit

        self._unit = new_unit

    def _get_unit(self):

        return self._unit

    unit = property(_get_unit, _set_unit,
                    doc="""Gets or sets the unit for the parameter""")

    # Define the property "value" with a control that the parameter cannot be set
    # outside of its bounds

    @property
    def value(self):
        """Return current parameter value"""

        return self._value

    def has_auxiliary_variable(self):

        # ParameterBase cannot have auxiliary variable (only Parameter has)

        return False

    # I use the decorator here (instead of the usual style for getters and setters) because
    # children classes need to override the getter, and using decorators is the only way
    # to achieve that
    @value.setter
    @accept_quantity(float, allow_none=False)
    def value(self, value):
        """Sets the current value of the parameter, ensuring that it is within the allowed range."""

        if self._min_value is not None and value < self._min_value:

            raise SettingOutOfBounds(
                "Trying to set parameter {0} = {1}, which is less than the minimum allowed {2}".format(
                    self.name, value, self.min_value))

        if self._max_value is not None and value > self._max_value:
            raise SettingOutOfBounds(
                "Trying to set parameter {0} = {1}, which is more than the maximum allowed {2}".format(
                    self.name, value, self.max_value))

        # Issue a warning if there is an auxiliary variable, as the setting does not have any effect
        if self.has_auxiliary_variable():

            with warnings.catch_warnings():

                warnings.simplefilter("always", RuntimeWarning)

                warnings.warn("You are trying to assign to a parameter which is either linked or "
                              "has auxiliary variables. The assignment has no effect.", RuntimeWarning)

        else:

            # Save the value as a pure floating point to avoid the overhead of the astropy.units machinery when
            # not needed

            self._value = value

        # Call the callbacks (if any)

        for callback in self._callbacks:

            try:

                callback(self)

            except:

                raise NotCallableOrErrorInCall(
                    "Could not use callback for parameter %s with value %s" % (self.name, value))

    @property
    def as_quantity(self):
        """
        Return the current value with its units (as an astropy.Quantity instance)

        :return: an instance of astropy.Quantity)
        """
        return self._value * self._unit

    def in_unit_of(self, unit, as_quantity=False):
        """
        Return the current value transformed to the new units

        :param unit: either an astropy.Unit instance, or a string which can be converted to an astropy.Unit
        instance, like "1 / (erg cm**2 s)"
        :param as_quantity: if True, the method return an astropy.Quantity, if False just a floating point number.
        Default is False
        :return: either a floating point or a astropy.Quantity depending on the value of "as_quantity"
        """

        new_unit = u.Unit(unit)

        new_quantity = self.as_quantity.to(new_unit)

        if as_quantity:

            return new_quantity

        else:

            return new_quantity.value

    # Define the property "min_value"

    def _get_min_value(self):
        """Return current minimum allowed value"""

        return self._min_value

    @accept_quantity(float, allow_none=True)
    def _set_min_value(self, min_value):

        """Sets current minimum allowed value"""

        # Store the minimum as a pure float

        self._min_value = min_value

        # Check that the current value of the parameter is still within the boundaries. If not, issue a warning

        if self._min_value is not None and self.value < self._min_value:

            warnings.warn("The current value of the parameter %s (%s) "
                          "was below the new minimum %s." % (self.name, self.value, self._min_value),
                          exceptions.RuntimeWarning)

            self._value = self._min_value

    min_value = property(_get_min_value, _set_min_value,
                         doc='Gets or sets the minimum allowed value for the parameter')

    def remove_minimum(self):
        """
        Remove the minimum from this parameter (i.e., it becomes boundless in the negative direction)
        """
        self._min_value = None

    # Define the property "max_value"

    def _get_max_value(self):

        """Return current minimum allowed value"""

        return self._max_value

    @accept_quantity(float, allow_none=True)
    def _set_max_value(self, max_value):
        """Sets current minimum allowed value"""

        self._max_value = max_value

        # Check that the current value of the parameter is still within the boundaries. If not, issue a warning

        if self._max_value is not None and self.value > self._max_value:

            warnings.warn("The current value of the parameter %s (%s) "
                          "was above the new maximum %s." % (self.name, self.value, self._max_value),
                          exceptions.RuntimeWarning)
            self._value = self._max_value

    max_value = property(_get_max_value, _set_max_value,
                         doc='Gets or sets the maximum allowed value for the parameter')

    def remove_maximum(self):
        """
        Remove the maximum from this parameter (i.e., it becomes boundless in the positive direction)
        """
        self._max_value = None

    def _set_bounds(self, bounds):
        """Sets the boundaries for this parameter to min_value and max_value"""

        # Use the properties so that the checks and the handling of units are made automatically

        min_value, max_value = bounds

        self.min_value = min_value

        self.max_value = max_value

    def _get_bounds(self):
        """Returns the current boundaries for the parameter"""

        return self.min_value, self.max_value

    bounds = property(_get_bounds, _set_bounds, doc="Gets or sets the boundaries (minimum and maximum) for this "
                                                    "parameter")

    def add_callback(self, callback):
        """Add a callback to the list of functions which are called immediately after the value of the parameter
        is changed. The callback must be a function accepting the current parameter as input. The return value of the
        callback is ignored. More than one callback can be specified. In that case, the callbacks will be called in the
        same order they have been entered."""

        self._callbacks.append(callback)

    def get_callbacks(self):
        """
        Returns the list of callbacks currently defined

        :return:

        """
        return self._callbacks

    def empty_callbacks(self):
        """Remove all callbacks for this parameter"""
        self._callbacks = []

    def duplicate(self):
        """
        Returns an exact copy of the current parameter
        """

        # Deep copy everything to make sure that there are no ties between the new instance and the old one

        new_parameter = copy.deepcopy(self)

        return new_parameter

    def to_dict(self, minimal=False):

        """Returns the representation for serialization"""

        data = collections.OrderedDict()

        if minimal:

            # In the minimal representation we just output the value

            data['value'] = self._to_python_type(self.value)

        else:

            # In the complete representation we output everything is needed to re-build the object

            data['value'] = self._to_python_type(self._value)
            data['desc'] = str(self.description)
            data['min_value'] = self._to_python_type(self._min_value)
            data['max_value'] = self._to_python_type(self._max_value)
            data['unit'] = str(self.unit.to_string())

        return data

    @staticmethod
    def _to_python_type(variable):
        """
        Returns the value in the variable handling also np.array of one element

        :param variable: input variable
        :return: the value of the variable having a python type (int, float, ...)
        """

        # Assume variable is a np.array, fall back to the case where variable is already a primitive type

        try:

            return variable.item()

        except AttributeError:

            return variable


class Parameter(ParameterBase):
    """

    Implements a numerical parameter. Optionally the parameter can vary according to an auxiliary law (see below).

    :param name: Name for the parameter
    :param value: Initial value
    :param min_value: minimum allowed value for the parameter (default: None)
    :param max_value: maximum allowed value for the parameter (default: None)
    :param delta: initial step used by some fitting engines (default: 0.1 * value)
    :param desc: description of parameter (default: '')
    :param free: whether the parameter is free or not (default: True)
    :param unit: the parameter units (default: dimensionless)
    :param prior: the parameter's prior (default: None)
    """

    def __init__(self, name=None, value=None, min_value=None, max_value=None, delta=None, desc=None, free=True, unit='',
                 prior=None):
        # NOTE: we need to set up _aux_variable immediately because we are overriding the value getter which
        # needs this

        # by default we have no auxiliary variable

        self._aux_variable = {}

        # This extends ParameterBase by adding the possibility for free/fix, and a delta for fitting purposes, as
        # well as a prior

        super(Parameter, self).__init__(name, value, min_value=min_value, max_value=max_value, desc=desc, unit=unit)

        self._free = bool(free)

        # Set delta if provided, otherwise use default

        if delta is not None:

            self._delta = delta

        else:

            # Default is 10% of the value, unless the value is zero, in which case the delta is 0.1

            if self._value == 0:

                self._delta = 0.1

            else:

                self._delta = abs(0.1 * self._value)

        # pre-defined prior is no prior
        self._prior = None

        # override the default if a prior has been given

        if prior is not None:

            # Use the property on purpose, so all checks and setups are applied

            self.prior = prior

        # Now perform a very lazy check that we can perform math operations on the delta

        if not _behaves_like_a_number(self._delta):
            raise TypeError("The provided delta is not a number")

        # Create a backup copy of the status of the parameter (useful when an auxiliary variable
        # is removed)
        self._old_free = self._free



    # Define the new get_value which accounts for the possibility of auxiliary variables
    @ParameterBase.value.getter
    def value(self):
        """Return current parameter value"""

        # If this parameter has a law of variation, use it to compute the current value

        # (remember: an empty dictionary test as False, a non-empty as true)

        if self._aux_variable:

            self._value = self._aux_variable['law'](self._aux_variable['variable'].value)

        return self._value

    # Define the property "delta"

    def _get_delta(self):
        return self._delta

    @accept_quantity(float, allow_none=False)
    def _set_delta(self, delta):
        """Sets the current delta for the parameter. The delta is used as initial step by some fitting engines"""

        self._delta = delta

    delta = property(_get_delta, _set_delta,
                     doc='''Gets or sets the delta for the parameter''')

    # Define the property "prior"

    def _get_prior(self):

        return self._prior

    def _set_prior(self, prior):
        """Set prior for this parameter. The prior must be a function accepting the current value of the parameter
        as input and giving the probability density as output."""

        # Try and call the prior with the current value of the parameter
        try:

            _ = prior(self._value)

        except:

            raise NotCallableOrErrorInCall("Could not call the provided prior. " +
                                           "Is it a function accepting the current value of the parameter?")

        try:

            prior.set_units(self.unit, u.dimensionless_unscaled)

        except AttributeError:

            raise NotCallableOrErrorInCall("It looks like the provided prior is not a astromodels function.")

        self._prior = prior

    prior = property(_get_prior, _set_prior,
                     doc='Gets or sets the current prior for this parameter. The prior must be a callable function '
                         "accepting the current value of the parameter as input and returning the probability "
                         "density as output")

    def has_prior(self):
        """
        Check whether the current parameter has a defined prior

        :return: True or False
        """

        return self._prior is not None

    def set_uninformative_prior(self, prior_class):
        """
        Sets the prior for the parameter to a uniform prior between the current minimum and maximum, or a
        log-uniform prior between the current minimum and maximum.

        NOTE: if the current minimum and maximum are not defined, the default bounds for the prior class will be used.

        :param prior_class : the class to be used as prior (either Log_uniform_prior or Uniform_prior, or a class which
        provide a lower_bound and an upper_bound properties)
        :return: (none)
        """

        prior_instance = prior_class()

        if self.min_value is None:

            raise ParameterMustHaveBounds("Parameter %s does not have a defined minimum. Set one first, then re-run "
                                          "set_uninformative_prior" % self.path)

        else:

            try:

                prior_instance.lower_bound = self.min_value

            except SettingOutOfBounds:

                raise SettingOutOfBounds("Cannot use minimum of %s for prior %s" % (self.min_value,
                                                                                    prior_instance.name))

        if self.max_value is None:

            raise ParameterMustHaveBounds("Parameter %s does not have a defined maximum. Set one first, then re-run "
                                          "set_uninformative_prior" % self.path)

        else: # pragma: no cover

            try:

                prior_instance.upper_bound = self.max_value

            except SettingOutOfBounds:

                raise SettingOutOfBounds("Cannot use maximum of %s for prior %s" % (self.max_value,
                                                                                    prior_instance.name))

        assert np.isfinite(prior_instance.upper_bound.value),"The parameter %s must have a finite maximum" % self.name
        assert np.isfinite(prior_instance.lower_bound.value),"The parameter %s must have a finite minimum" % self.name

        self._set_prior(prior_instance)

    # Define property "free"

    def _set_free(self, value=True):

        self._free = value

    def _get_free(self):

        return self._free

    free = property(_get_free, _set_free,
                    doc="Gets or sets whether the parameter is free or not. Use booleans, like: 'p.free = True' "
                        " or 'p.free = False'. ")

    # Define property "fix"

    def _set_fix(self, value=True):

        self._free = (not value)

    def _get_fix(self):

        return not self._free

    fix = property(_get_fix, _set_fix,
                   doc="Gets or sets whether the parameter is fixed or not. Use booleans, like: 'p.fix = True' "
                       " or 'p.fix = False'. ")

    def add_auxiliary_variable(self, variable, law):

        # Assign units to the law
        law.set_units(variable.unit, self.unit)

        # Test the law
        try:

            _ = law(variable.value)

        except: # pragma: no cover

            raise NotCallableOrErrorInCall("The provided law for the auxiliary variable failed on call")

        self._aux_variable['law'] = law
        self._aux_variable['variable'] = variable

        # Now add the law as an attribute (through the mother class DualAccessClass),
        # so the user will be able to access its parameters as this.name.parameter_name

        #self.add_attribute(law.name, law)

        # Now add the nodes

        if self._has_child(law.name):

            self._remove_child(law.name)

        self._add_child(law)

        # This parameter is not free anymore

        # First make a backup of the old status, so that it will be restored if the
        # law is removed
        self._old_free = self._free

        self.free = False

    def remove_auxiliary_variable(self):
        """
        Remove an existing auxiliary variable

        :return:
        """

        if not self.has_auxiliary_variable():

            # do nothing, but print a warning

            warnings.warn("Cannot remove a non-existing auxiliary variable", RuntimeWarning)

        else:

            # Remove the law from the children

            self._remove_child(self._aux_variable['law'].name)

            # Clean up the dictionary

            self._aux_variable = {}

            # Set the parameter to the status it has before the auxiliary variable was created

            self.free = self._old_free

    def has_auxiliary_variable(self):
        """
        Returns whether the parameter is linked to an auxiliary variable
        """
        if self._aux_variable:

            return True

        else:

            return False

    @property
    def auxiliary_variable(self):
        """
        Returns a tuple with the auxiliary variable and the law

        :return: tuple (variable, law)
        """
        return self._aux_variable['variable'], self._aux_variable['law']

    def _repr__base(self, rich_output=False):

        if not self.has_auxiliary_variable():

            representation = "Parameter %s = %s [%s]\n" \
                             "(min_value = %s, max_value = %s, delta = %s, free = %s)" % (self.name,
                                                                                          self.value,
                                                                                          self.unit,
                                                                                          self.min_value,
                                                                                          self.max_value,
                                                                                          self.delta,
                                                                                          self.free)

            if self._prior is not None:

                representation += " [prior: %s]" % self.prior.name

        else:

            representation = "Parameter %s = %s [%s]\n" \
                             "(linked to auxiliary variable '%s' with law '%s')" % (self.name, self.value,
                                                                                    self.unit,
                                                                                    self._aux_variable['variable'].name,
                                                                                    self._aux_variable['law'].name)

        return representation

    def to_dict(self, minimal=False):

        """Returns the representation for serialization"""

        data = super(Parameter, self).to_dict()

        if minimal:

            # No need to add anything
            pass

        else:

            # In the complete representation we output everything is needed to re-build the object

            if self.has_auxiliary_variable():

                # Store the function and the auxiliary variable

                data['value'] = 'f(%s)' % self._aux_variable['variable']._get_path()
                #data['function_of'] = {self._aux_variable['variable'].name: self._aux_variable['variable'].to_dict()}

                aux_variable_law_data = collections.OrderedDict()
                aux_variable_law_data[ self._aux_variable['law'].name ] = self._aux_variable['law'].to_dict()

                data['law'] = aux_variable_law_data

            # delta and free are attributes of Parameter, but not of ParameterBase

            data['delta'] = self._to_python_type(self._delta)
            data['free'] = self.free

            if self.has_prior():

                data['prior'] = {self.prior.name: self.prior.to_dict()}

        return data

    def get_randomized_value(self, variance=0.1):

        # Get a value close to the current value, but not identical
        # (used for the inizialization of Bayesian samplers)

        if (self._min_value is not None) or (self._max_value is not None):

            # If _value is zero, then std will be zero, which doesn't make sense
            assert self._value != 0, "You cannot randomize parameter %s because its value is exactly zero" % self.path

            # Bounded parameter. Use a truncated normal so we are guaranteed
            # to have a random value within the boundaries

            std = abs(variance * self._value)

            if self._min_value is not None:

                a = ( self._min_value - self._value ) / std

            else:

                a = - np.inf

            if self._max_value is not None:

                b = (self._max_value - self._value) / std

            else:

                b = np.inf

            sample = scipy.stats.truncnorm.rvs( a, b, loc = self._value, scale = std, size = 1)

            if (self._min_value is not None and sample < self._min_value) or \
               (self._max_value is not None and sample > self._max_value): # pragma: no cover

                # This should never happen

                raise AssertionError("Got a sample outside of the boundaries of the truncated normal distribution")

            return sample[0]

        else:

            #The parameter has no boundaries

            return np.random.normal( self._value, abs(variance * self._value) )


class IndependentVariable(ParameterBase):
    """
    An independent variable like time or energy.
    """

    # Override the constructor to make the unit specification mandatory

    def __init__(self, name, value, unit, min_value=None, max_value=None, desc=None):

        super(IndependentVariable, self).__init__(name, value, unit=unit,
                                                  min_value=min_value,
                                                  max_value=max_value,
                                                  desc=desc)

    def _repr__base(self, rich_output=False):

        return "IndependentVariable %s = %s\n" \
               "(min_value = %s, max_value = %s)" % (self.name,
                                                     self.value,
                                                     self.min_value,
                                                     self.max_value)
