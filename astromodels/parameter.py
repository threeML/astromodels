__author__ = 'giacomov'

__doc__ = """
=============
Description
=============

Create a parameter with::

  > p = Parameter("param1",1.0)

This will create a parameter called "param1" with initial value 1.0 and no boundaries. It will also have a
default delta of 10% the value, in this case 0.1 * 1 = 0.1.

You can use optional keywords to define the boundaries of the parameter as well as its delta::

  > p = Parameter("param1",1.0, min = -10.0, max = 35, delta = 3)

This will create a parameter bounded between -10.0 and 35, with initial delta = 3.

You can set the current value of the parameter as::

  > p.value = 2.5

or get its current value as::

  > current_value = p.value

You can set its boundaries with::

  > p.set_bounds(-10, 35)

A value of None for either the minimum or the maximum will remove the boundary in that direction.

You can get the current boundaries as::

  > min_value, max_value = p.get_bounds()

You can set or get the value for delta as::

  > p.delta = 0.3
  > current_delta = p.delta


==========
Callbacks
==========

The Parameter class provides also a functionality to set callbacks which are called whenever the value of the
parameter changes. This can be used for tracing purposes or other things. This code for example sets a simple
callback which just prints the value of the parameter any time it is changed::

  > def callback( value ):
        print("The parameter changed value. It is now %s" % value )
  > p.add_callback( callback )
  > p.value = 10.0
  The parameter changed value. It is now 10.0

This is instead a more elaborated example where a class is used to keep track of all the values the parameter has
assumed with time::

  > class MyRecorder( object ):
      def __init__(self):
          self.values = []
      def __call__(self, value):
          self.values.append(value)
      def print_values(self):
          print("This is the history of all the values this parameter has been set to:")
          print(self.values)
  > recorder = MyRecorder()
  > p.add_callback( recorder )
  > p.value = 5.0
  > p.value = -2.0
  > p.value = 18.0
  > recorder.print_values()
  This is the history of all the values this parameter has been set to:
  [5.0, -2.0, 18.0]

More than one callback can be registered. The callbacks will be executed in the order they are entered. To clear all
callbacks use the method empty_callbacks.
.. _parameter_auxvar:
===================
Auxiliary variable
===================

Sometimes it is useful to specify the parameter as a function of another variable ("auxiliary variable"). For
example, we want the parameter to vary with time according to a linear law. This can be done as follows.

First define the auxiliary variable "time" starting at 0, and a linear law 3.2 t + 5.6::

  > t = IndependentVariable("time",0.0)
  > law = lambda x: 3.2 * x + 5.6

For example, in t=0 our law evaluate to::

  > print( law(t.value) )
  5.6

Now we link the value of the parameter p to the law we just created::

  > p = Parameter("test",1.0)
  > p.add_auxiliary_variable(t,law)

Since t is still at 0, we can immediately see that the value of the parameter is now law(0) = 3.2 * 0 + 5.6 = 5.6::

  > p.value
  5.6

Note hence that the init values we gave during the constructor looses its meaning in this case, by design.

If we now change the value of t, the value of the parameter will change automatically according to the law::

  > t.value = 10.0
  > p.value
  37.6

This becomes even more powerful when the law contains parameters itself::

  > a = Parameter("a",3.2)
  > b = Parameter("b",5.6)
  > law = lambda x: a.value * x + b.value
  > p.add_auxiliary_variable(t,law)

If we now change the value of one of the parameters, we see the parameter p changing accordingly::

  > t.value = 10.0
  > print(p.value)
  5.6
  > a.value = 1.23
  > print(p.value)
  17.9
  > b.value = -2.12
  > print(p.value)
  10.18

"""

import collections
import copy
import exceptions
import warnings

import astropy.units as u
import scipy.stats

# The following import is necessary so that min_value and max_value can be specified as things like
# '2 * np.pi'
import numpy as np

from math import log10

from astromodels.tree import Node


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
class SettingOutOfBounds(exceptions.Exception):
    pass


# Exception for when an object should be callable but it isn't
class NotCallableOrErrorInCall(exceptions.Exception):
    pass


class WarningUnitsAreSlow(Warning):
    pass


class IndependentVariableCannotBeLinked(exceptions.Exception):
    pass


class CannotUseLogScale(exceptions.Exception):
    pass


class ParameterBase(Node):

    def __init__(self, name, value, min_value=None, max_value=None, desc=None, unit=u.dimensionless_unscaled):

        # Make this a node

        Node.__init__(self, name)

        # NOTE: we avoid enforcing a particular type or even that initial_values, max, min and delta are python numbers.
        # Indeed, as long as they behave as numbers we are going to be fine. We want to keep the possibility to use
        # numbers coming from C, C++ or other sources, for example. Note however that string
        # parameters are not supported

        # Callbacks are executed any time the value for the parameter changes (i.e., its value changes)

        # We start from a empty list of callbacks.
        self._callbacks = []

        # Assign to members

        # NOTE: the internal value stored in _value is always unscaled.
        # In the Parameter children, if the user sets this a logarithmic
        # parameter, _value will still be in linear scale. The same is true for the _min_value and _max_value.
        # (see note in the constructor of the Parameter class)

        # NOTE2: we use the eval of the str so that the value can be things like 1/2.0, np.pi and similar

        self._value = float(eval(str(value)))

        # Store the units as an astropy.units.Unit instance

        self._unit = u.Unit(unit)

        # Set minimum if provided, otherwise use default
        # (use the property so the checks that are there are performed also on construction)

        if min_value is not None:

            self.min_value = min_value

        else:

            self.min_value = None

        # Set maximum if provided, otherwise use default

        if max_value is not None:

            self.max_value = max_value

        else:

            self.max_value = None

        # Store description

        self._desc = desc

        # Make the description the documentation as well

        self.__doc__ = desc

        # Now perform a very lazy check that we can perform math operations on all members

        if not _behaves_like_a_number(self._value):
            raise TypeError("The provided initial value is not a number")

        if self._min_value is not None:

            if not _behaves_like_a_number(self._min_value):
                raise TypeError("The provided minimum value is not a number")

        if self._max_value is not None:

            if not _behaves_like_a_number(self._max_value):
                raise TypeError("The provided maximum value is not a number")

    # Define the property 'description' and make it read-only

    @property
    def description(self):
        return self._desc

    # Define the property 'unit'
    def set_unit(self, new_unit):

        self._unit = u.Unit(new_unit)

    def get_unit(self):

        return self._unit

    unit = property(get_unit, set_unit,
                    doc="""Gets or sets the unit for the parameter""")

    # Define the property "value" with a control that the parameter cannot be set
    # outside of its bounds

    @property
    def value(self):
        """Return current parameter value"""

        return self._value

    # I use the decorator here (instead of the usual style for getters and setters) because
    # children classes need to override the getter, and using decorators is the only way
    # to achieve that

    @value.setter
    def value(self, value):
        """Sets the current value of the parameter, ensuring that it is within the allowed range"""

        try:
            if self._min_value is not None and value < self._min_value:
                raise SettingOutOfBounds(
                    "Trying to set parameter {0} = {1}, which is less than the minimum allowed {2}".format(
                        self.name, value, self._min_value))

            if self._max_value is not None and value > self._max_value:
                raise SettingOutOfBounds(
                    "Trying to set parameter {0} = {1}, which is more than the maximum allowed {2}".format(
                        self.name, value, self._max_value))
        except u.UnitsError:

            raise ValueError("If you want to use astropy units, you need to use the set() method for the parameter.")

        # Call the callbacks (if any)

        for callback in self._callbacks:

            try:

                callback(value)

            except:

                raise NotCallableOrErrorInCall(
                    "Could not get callback for parameter %s with value %s" % (self.name, value))

        self._value = value

    def set(self, quantity):
        """
        Set the value of the parameter to the given value, given as a astropy.Quantity (with units).
        Note that this is quite slow due to the astropy machinery, so do not use this functionality in computing
        intensive situations. In such situations you are likely to use always the same unit, so you'd be better
        off by performing your conversion directly on the values and use my_parameter.value = your_value.

        :parameter quantity: an astropy.Quantity instance
        :return: (none)
        """

        try:

            # This works if value is a astropy.Quantity

            new_value = quantity.to(self.unit).value

        except AttributeError:

            # We get here if instead value is a simple number

            raise ValueError("You need to use a astropy.quantity object for the set() method.")

        else:

            # Even if the to() method works, we need to warn the user that this is
            # very slow, and should only be used in interactive sessions for convenience

            warnings.warn("Using units is convenient but slow. Do not use them during computing-intensive work.",
                          WarningUnitsAreSlow)

            # Set the value using the property, so the usual check about boundaries and so on is still
            # performed

            self.value = new_value

    def get(self, new_unit):

        try:

            value = (self.value * self.unit).to(new_unit)

        except AttributeError:

            # We get here if instead value is a simple number

            raise ValueError("You need to use a astropy.quantity object for the get() method.")

        else:

            # Even if the to() method works, we need to warn the user that this is
            # very slow, and should only be used in interactive sessions for convenience

            warnings.warn("Using units is convenient but slow. Do not use them during computing-intensive work.",
                          WarningUnitsAreSlow)

            return value

    # Define the property "min_value"

    def get_min_value(self):
        """Return current minimum allowed value"""
        return self._min_value

    def set_min_value(self, min_value):
        """Sets current minimum allowed value"""

        if isinstance(min_value, str):

            # This allows to specify things like np.pi as min_value

            min_value = eval(min_value)

        if min_value is None:

            self._min_value = None

        else:

            if not _behaves_like_a_number(min_value):
                raise ValueError("Provided minimum value cannot be interpreted as a number nor None")

            self._min_value = float(min_value)

            # Check that the current value of the parameter is still within the boundaries. If not, issue a warning

            if self.value < self.min_value:
                warnings.warn("The current value of the parameter %s is below the new minimum." % self.name,
                              exceptions.RuntimeWarning)

    min_value = property(get_min_value, set_min_value,
                         doc='Gets or sets the minimum allowed value for the parameter')

    # Define the property "max_value"

    def get_max_value(self):

        """Return current minimum allowed value"""

        return self._max_value

    def set_max_value(self, max_value):
        """Sets current minimum allowed value"""

        if isinstance(max_value, str):

            # This allows to specify things like 'np.pi / 3' as max_value

            max_value = eval(max_value)

        if max_value is None:

            self._max_value = None

        else:

            if not _behaves_like_a_number(max_value):
                raise ValueError("Provided maximum value is not a number nor None")

            self._max_value = float(max_value)

            # Check that the current value of the parameter is still within the boundaries. If not, issue a warning

            if self.value > self.max_value:
                warnings.warn("The current value of the parameter %s is above the new maximum." % self.name,
                              exceptions.RuntimeWarning)

    max_value = property(get_max_value, set_max_value,
                         doc='Gets or sets the maximum allowed value for the parameter')

    def set_bounds(self, min_value, max_value):
        """Sets the boundaries for this parameter to min_value and max_value"""

        # Use the properties so that the checks are made automatically

        self.min_value = min_value

        self.max_value = max_value

    def get_bounds(self):
        """Returns the current boundaries for the parameter"""

        return self.min_value, self.max_value

    def change_name(self, new_name):
        """
        Change the name of this parameter. Use with caution!

        :return: none
        """
        self._name = new_name

    def add_callback(self, callback):
        """Add a callback to the list of functions which are called whenever the value of the parameter is changed.
        The callback must be a function accepting the current value as input. The return value of the callback is
        ignored. More than one callback can be specified. In that case, the callbacks will be called in the same order
        they have been entered."""

        self._callbacks.append(callback)

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

            data['value'] = self._to_python_type(self.value) * self.unit

        else:

            # In the complete representation we output everything is needed to re-build the object

            data['value'] = self._to_python_type(self.value)

            data['min_value'] = self._to_python_type(self.min_value)
            data['max_value'] = self._to_python_type(self.max_value)

            # Try to use the unit as an astropy.Unit instance
            # If that doesn't work, treat the unit as a simple string
            try:

                data['unit'] = str(self.unit.to_string())

            except AttributeError:

                data['unit'] = str(self.unit)

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
    """

    def __init__(self, name, value, min_value=None, max_value=None, delta=None, desc=None, free=True, unit=''):

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

        # Now perform a very lazy check that we can perform math operations on the delta

        if not _behaves_like_a_number(self._delta):
            raise TypeError("The provided delta is not a number")

        # Create a backup copy of the status of the parameter (useful when an auxiliary variable
        # is removed)
        self._old_free = self._free

        # Now create the pieces needed for the mechanism which allow the user to transform a parameter
        # in the log scale only for the fitting or bayesian engine. This is how it works:
        # - if log is False, .value and .scaled_value are the same
        # - if log is True, .scaled_value is the log10 of .value. Scaled_value should be used by fitting engines
        #   and bayesian samplers, which will then work on the parameter in log scale. The user however will not
        #   notice this and will continue to see and work with the parameter in linear scale.

        # By default the log scale is False

        self._log = False

    # Define the new get_value which accounts for the possibility of auxiliary variables
    @ParameterBase.value.getter
    def value(self):
        """Return current parameter value"""

        # If this parameter has a law of variation, use it to compute the current value

        # (remember: an empty dictionary test as False, a non-empty as true)

        if self._aux_variable:

            return self._aux_variable['law'](self._aux_variable['variable'].value)

        else:

            return self._value

    def _get_log(self):

        return self._log

    def _set_log(self, new_value):

        # Check that the minimum and the maximum are larger than zero. It is not possible to use the log scale
        # on a parameter which can take zero or negative values, or which have no boundaries

        if (self.min_value is None or self.min_value <= 0) or (self.max_value is None or self.max_value <= 0):

            raise CannotUseLogScale("Cannot use log scale for parameter %s. It needs to have both a defined minimum "
                                    "and maximum, and both must be strictly larger than 0" % self.name)

        self._log = bool(new_value)

    log = property(_get_log, _set_log, doc="Sets or gets the status of the log scale switch. If log is True, "
                                           "the value, the minimum and the maximum for this parameter returned by "
                                           "scaled_value, scaled_min_value and scaled_max_value will be equal to the "
                                           "log10 of respectively the value, min_value and max_value. This is useful "
                                           "for certain fitting engine such as MINUIT if the parameter has a wide "
                                           "dynamical range.")

    # Now define the machinery needed for the log scale mechanism (see the note in the constructor)
    def _get_scaled_value(self):

        if self.log:

            return log10(self.value)

        else:

            return self.value

    def _set_scaled_value(self, new_value):

        if self.log:

            # We assume that the new value is actually the logarithm of the desired value

            self.value = pow(10, new_value)

        else:

            self.value = new_value

    scaled_value = property(_get_scaled_value, _set_scaled_value, 
                            doc="Get the current value of the parameter. You should use this in place of .value "
                                "for fitting engines and Bayesian samplers")

    def _get_scaled_min_value(self):

        if self.log:

            return log10(self.min_value)

        else:

            return self.min_value

    def _set_scaled_min_value(self, new_minimum):

        if self.log:

            # We assume that the new minimum is actually the logarithm of the desired minimum

            self.min_value = pow(10, new_minimum)

        else:

            self.min_value = new_minimum

    scaled_min_value = property(_get_scaled_min_value, _set_scaled_min_value, 
                                doc="Get the current minimum for the parameter. You should use this in place of "
                                    ".min_value for fitting engines and Bayesian samplers")

    def _get_scaled_max_value(self):

        if self.log:

            return log10(self.max_value)

        else:

            return self.max_value

    def _set_scaled_max_value(self, new_maximum):

        if self.log:

            # We assume that the new maximum is actually the logarithm of the desired maximum

            self.max_value = pow(10, new_maximum)

        else:

            self.max_value = new_maximum

    scaled_max_value = property(_get_scaled_max_value, _set_scaled_max_value,
                                doc="Get the current maximum for the parameter. You should use this in place of "
                                    ".max_value for fitting engines and Bayesian samplers")

    # Define the property "delta"

    def get_delta(self):
        """Gets the current value for the delta of the parameter"""
        return self._delta

    def set_delta(self, delta):
        """Sets the current delta for the parameter. The delta is used as initial step by some fitting engines"""

        if not _behaves_like_a_number(delta):
            raise ValueError("Provided delta is not a number" % delta)

        self._delta = delta

    delta = property(get_delta, set_delta,
                     doc='''Gets or sets the delta for the parameter''')

    # Define the property scaled_delta to be used by fitting engines and Bayesian samplers
    def _get_scaled_delta(self):

        if self.log:

            # Keep the delta appropriate on the log scale

            return log10(self.value + self.delta) - log10(self.value)

        else:

            return self.delta

    def _set_scaled_delta(self, new_delta):

        if self.log:

            # The new delta is in the log scale. This equation comes from inverting the equation in the
            # getter

            self.delta = self.value * (pow(10, new_delta) - 1)

        else:

            self.delta = new_delta

    scaled_delta = property(_get_scaled_delta, _set_scaled_delta,
                            doc="Get the current maximum for the parameter. You should use this in place of "
                                ".max_value for fitting engines and Bayesian samplers")

    # Define the property "prior"

    def get_prior(self):

        if self._prior is None:
            raise RuntimeError("There is no defined prior for parameter %s" % self.name)

        return self._prior

    def set_prior(self, prior):
        """Set prior for this parameter. The prior must be a function accepting the current value of the parameter
        as input and giving the probability density as output."""

        # Try and call the prior with the current value of the parameter
        try:

            _ = prior(self._value)

        except:

            raise NotCallableOrErrorInCall("Could not call the provided prior. " +
                                           "Is it a function accepting the current value of the parameter?")

        self._prior = prior

    prior = property(get_prior, set_prior,
                     doc='Gets or sets the current prior for this parameter. The prior must be a callable function '
                         "accepting the current value of the parameter as input and returning the probability "
                         "density as output")

    def has_prior(self):
        """
        Check whether the current parameter has a defined prior

        :return: True or False
        """

        return self._prior is not None

    # Define property "free"

    def set_free(self, value=True):

        self._free = value

    def get_free(self):

        return self._free

    free = property(get_free, set_free,
                    doc="Gets or sets whether the parameter is free or not. Use booleans, like: 'p.free = True' "
                        " or 'p.free = False'. ")

    # Define property "fix"

    def set_fix(self, value=True):

        self._free = (not value)

    def get_fix(self):

        return not self._free

    fix = property(get_fix, set_fix,
                   doc="Gets or sets whether the parameter is fixed or not. Use booleans, like: 'p.fix = True' "
                       " or 'p.fix = False'. ")

    def add_auxiliary_variable(self, variable, law):

        # Test the law
        try:

            _ = law(variable.value)

        except AttributeError:

            raise NotCallableOrErrorInCall("Cannot access the .value attribute of the aux. variable. "
                                           "Is it of the proper class?")

        except:

            raise NotCallableOrErrorInCall("The provided law for the auxiliary variable failed on call")

        self._aux_variable['law'] = law
        self._aux_variable['variable'] = variable

        # Now add the law as an attribute (through the mother class DualAccessClass),
        # so the user will be able to access its parameters as this.name.parameter_name

        #self.add_attribute(law.name, law)

        # Now add the nodes

        if law.name in self.get_children():

            self.remove_child(law.name)

        self.add_child(law)

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

            warnings.warn("Cannot remove a non-existing auxiliary variable")

        else:

            # Remove the law from the children

            self.remove_child(self._aux_variable['law'].name)

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

        representation = ""

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
                             "(linked to auxiliary variable '%s' with law '%s')" % (self.name, self.value, unit,
                                                                                    self._aux_variable['variable'].name,
                                                                                    self._aux_variable['law'].name)

        return representation

    def to_dict(self, minimal=False):

        """Returns the representation for serialization"""

        if minimal:

            data = collections.OrderedDict()

            # In the minimal representation we just output the value

            data['value'] = "%s [%s]" % (self._to_python_type(self.value), self.unit)

        else:

            data = super(Parameter, self).to_dict()

            # In the complete representation we output everything is needed to re-build the object

            if self.has_auxiliary_variable():

                # Store the function and the auxiliary variable

                data['value'] = 'f(%s)' % ".".join(self._aux_variable['variable'].get_path())
                #data['function_of'] = {self._aux_variable['variable'].name: self._aux_variable['variable'].to_dict()}

                aux_variable_law_data = collections.OrderedDict()
                aux_variable_law_data[ self._aux_variable['law'].name ] = self._aux_variable['law'].to_dict()

                data['law'] = aux_variable_law_data

            # delta and free are attributes of Parameter, but not of ParameterBase

            data['delta'] = self._to_python_type(self.delta)
            data['free'] = self.free

            if self.has_prior():

                data['prior'] = {self.prior.name: self.prior.to_dict()}

        return data

    def get_randomized_value(self, variance=0.1):

        # Get a value close to the current value, but not identical
        # (used for the inizialization of Bayesian samplers)

        if (self.min_value is not None) or (self.max_value is not None):

            # Bounded parameter. Use a truncated normal so we are guaranteed
            # to have a random value within the boundaries

            std = abs(variance * self.value)

            if self.min_value is not None:

                a = ( self.min_value - self.value ) / std

            else:

                a = - np.inf

            if self.max_value is not None:

                b = (self.max_value - self.value) / std

            else:

                b = np.inf

            sample = scipy.stats.truncnorm.rvs( a, b, loc = self.value, scale = std, size = 1)

            if (self.min_value is not None and sample < self.min_value) or \
               (self.max_value is not None and sample > self.max_value):

                raise RuntimeError("This is a bug!!")

            return sample[0]

        else:

            #The parameter has no boundaries

            return np.random.normal( self.value, variance * self.value )


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
