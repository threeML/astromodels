__author__ = "giacomov"

__doc__ = """"""

import collections
import contextlib
import copy
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import astropy.units as u
import numpy as np
import scipy.stats

from astromodels.core.parameter_transformation import ParameterTransformation
from astromodels.utils.configuration import astromodels_config
from astromodels.utils.logging import setup_logger

from .tree import Node

from .thread_safe_unit_format import ThreadSafe

log = setup_logger(__name__)



@contextlib.contextmanager
def turn_off_parameter_transforms() -> None:
    """
    deactivate parameter transforms temporarily
    :return: None
    """

    # store the old status
    
    old_status = bool(astromodels_config.modeling.use_parameter_transforms)

    # turn off the configuration value
    # which will cause the transforms to deactivate in all parameters
    
    astromodels_config.modeling.use_parameter_transforms = False

    # do your thing
    yield

    # now reset the status to the old value

    astromodels_config.modeling.use_parameter_transforms = old_status
    


def _behaves_like_a_number(obj) -> bool:
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

                    else:  # pragma: no cover

                        log.exception("You cannot pass None as argument for "
                                      "method %s of %s" %
                                      (method.__name__, instance.name))

                        raise TypeError()

                else:  # pragma: no cover

                    log.exception(
                        "You need to pass either a %s or a astropy.Quantity "
                        "to method %s of %s" %
                        (input_type.__name__, method.__name__, instance.name))

                    raise TypeError()

        return handle_quantity

    return accept_quantity_wrapper


class ParameterBase(Node):
    """

    :param name: name for parameter
    :param value: initial value
    :param min_value: minimum
    :param max_value: maximum
    :param desc: description
    :param unit: units (string or astropy.Unit)
    :param transformation: a class which implements a .forward and a .backward method to transform forth and back from
        face value (the value exposed to the user) to the internal value (the value exposed to the fitting engine)
    """
    def __init__(
        self,
        name: str,
        value: float,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None,
        desc: Optional[str] = None,
        unit: u.Unit = u.dimensionless_unscaled,
        transformation: Optional[ParameterTransformation] = None,
    ):

        # Make this a node

        Node.__init__(self, name)

        # Make a static name which will never change (not even after a _change_name call)
        self._static_name: str = str(name)

        # Callbacks are executed any time the value for the parameter changes (i.e., its value changes)

        # We start from a empty list of callbacks.
        self._callbacks = []

        # Assign to members

        # Store the units as an astropy.units.Unit instance

        self._unit = self._safe_assign_unit(unit)

        # A ParameterBase instance cannot have auxiliary variables
        self._aux_variable: Optional[Dict[str, Any]] = None

        # Set min and max to None first so that the .value setter will work,
        # we will override them later if needed
        self._external_min_value: Optional[float] = None
        self._external_max_value: Optional[float] = None

        # Store the transformation. This allows to disentangle the value of the parameter which the user interact
        # width with the value of the parameter the fitting engine (or the Bayesian sampler) interact with
        if transformation is not None:
            if not isinstance(transformation, ParameterTransformation):

                log.error("transformation is not of ParameterTransform ")

                raise AssertionError()

        self._transformation: Optional[
            ParameterTransformation] = transformation

        # save the transformation so that it can be restored
        
        self._original_transformation: Optional[
            ParameterTransformation] = transformation

        
        # Let's store the init value

        # NOTE: this will be updated immediately by the _set_value method of the "value" property
        self._internal_value: Optional[float] = None
        # If the value is a Quantity, deal with that

        if isinstance(value, u.Quantity):

            # If the user did not specify an ad-hoc unit, use the unit
            # of the Quantity

            if self._unit == u.dimensionless_unscaled:

                self._unit = value.unit

            # Convert the value to the provided unit (if necessary)

            self.value: float = value.to(self._unit).value

        else:

            self.value = value

        # Set minimum if provided, otherwise use default
        # (use the property so the checks that are there are performed also on construction)

        self.min_value: Optional[float] = min_value

        # Set maximum if provided, otherwise use default

        # this will be overwritten immediately in the next line
        self.max_value: Optional[float] = max_value

        # Store description

        self._desc: Optional[float] = desc

        # Make the description the documentation as well

        self.__doc__ = desc

        # Now perform a very lazy check that we can perform math operations on value, minimum and maximum
        # (i.e., they are numbers)

        if not _behaves_like_a_number(self.value):

            log.exception("The provided initial value is not a number")

            raise TypeError()

    def _repr__base(self, rich_output):  # pragma: no cover

        raise NotImplementedError(
            "You need to implement this for the actual Parameter class")

    # Define the property 'description' and make it read-only

    @property
    def static_name(self) -> str:
        """
        Returns a name which will never change, even if the name of the parameter does (for example in composite
        functions)

        :return : the static name
        :type : str
        """

        return self._static_name

    @property
    def description(self) -> Optional[str]:
        """
        Return a description of this parameter

        :return: a string cointaining a description of the meaning of this parameter
        """
        return self._desc

    @staticmethod
    def _safe_assign_unit(input_unit) -> u.Unit:

        # We first try to use our own, thread-safe format, if we fail then we try the astropy one

        try:

            new_unit = u.Unit(input_unit, format="threadsafe")

        except ValueError:

            # Try with the default format of astropy
            new_unit = u.Unit(input_unit)

        return new_unit

    # Define the property 'unit'
    def _set_unit(self, input_unit) -> None:

        # This will fail if the input is not valid

        new_unit = self._safe_assign_unit(input_unit)

        # Now transform the current _value in the new unit, unless the current unit is dimensionless, in which
        # case there is no transformation to make

        # (self._unit is the OLD unit here)

        if self._unit != u.dimensionless_unscaled:

            # This will fail if the new unit is not compatible with the old one

            try:

                self.value = self.as_quantity.to(new_unit).value

            except u.UnitConversionError:

                if new_unit == u.dimensionless_unscaled:

                    new_unit_name = "(dimensionless)"

                else:

                    new_unit_name = new_unit

                log.exception("Cannot convert the value %s from %s to the "
                              "new units %s" %
                              (self.value, self._unit, new_unit_name))

                raise CannotConvertValueToNewUnits()

        else:

            # This is possibly the first time the unit is set
            pass

        # Finally store the new unit

        self._unit = new_unit

    def _get_unit(self):

        return self._unit

    unit = property(_get_unit,
                    _set_unit,
                    doc="""Gets or sets the unit for the parameter""")

    @property
    def as_quantity(self) -> u.Quantity:
        """
        Return the current value with its units (as an astropy.Quantity instance)

        :return: an instance of astropy.Quantity)
        """
        return self.value * self._unit

    def in_unit_of(self, unit, as_quantity=False) -> u.Quantity:
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

    @property
    def has_auxiliary_variable(self) -> bool:

        return self._aux_variable is not None


    #@property
    def has_transformation(self) -> bool:

        if (self._transformation is None) or (not astromodels_config.modeling.use_parameter_transforms):

            return False

        else:

            return True

    @property
    def transformation(self) -> Optional[ParameterTransformation]:

        if astromodels_config.modeling.use_parameter_transforms:
        
            return self._transformation

        else:

            return None
        
    def remove_transformation(self) -> None:
        """
        
        remove any transform on the parameter

        useful in bayesian fits where 
        we do not care about 
        transformations but want speed

        :returns: 

        """
        # first get the real value
        old_value = self.value

        # now erase the transform
        
        self._transformation = None

        # restore the value which
        # will now be done without
        # the transformation
        
        self.value = old_value

    def restore_transformation(self) -> None:
        """
        
        restore the original transformation
        if it had been removed

        :returns: 

        """
        # first get the real value
        old_value = self.value

        # now reset
        
        self._transformation = self._original_transformation

        # restore the value which
        # will now be done without
        # the transformation
        
        self.value = old_value

        
    
    def internal_to_external_delta(self, internal_value, internal_delta):
        """
        Transform an interval from the internal to the external reference (through the transformation). It is useful
        if you have for example a confidence interval in internal reference and you want to transform it to the
        external reference

        :param interval_value: value in internal reference
        :param internal_delta: delta in internal reference
        :return: value and delta in external reference
        """

        external_value = self.transformation.backward(internal_value)
        bound_internal = internal_value + internal_delta
        bound_external = self.transformation.backward(bound_internal)
        external_delta = bound_external - external_value

        return external_value, external_delta

    # Define the property "value" with a control that the parameter cannot be set
    # outside of its bounds

    def _get_value(self) -> float:
        """Return current parameter value"""

        # This is going to be true (possibly) only for derived classes. It is here to make the code cleaner
        # and also to avoid infinite recursion

        if self._aux_variable:

            return self._aux_variable["law"](
                self._aux_variable["variable"].value)

        if self._transformation is None:

            return self._internal_value

        else:

            # A transformation is set. Transform back from internal value to true value
            #
            # print("Interval value is %s" % self._internal_value)
            # print("Returning %s" % self._transformation.backward(self._internal_value))

            return self._transformation.backward(self._internal_value)

    # NOTE: this function should only be used by the user. Fitting engines should only deal with
    # _get_internal_value and _set_internal_value
    @accept_quantity(
        float, allow_none=False
    )  # This means that the method will always receive a float
    def _set_value(self, new_value) -> None:
        """Sets the current value of the parameter, ensuring that it is within the allowed range."""

        if (self.min_value is not None) and (new_value < self.min_value) and not astromodels_config.modeling.ignore_parameter_bounds:

            log.error(
                "Trying to set parameter {0} = {1}, which is less than the minimum allowed {2}"
                .format(self.name, new_value, self.min_value))

            raise SettingOutOfBounds()

        if (self.max_value is not None) and (new_value > self.max_value) and not astromodels_config.modeling.ignore_parameter_bounds:

            log.error(
                "Trying to set parameter {0} = {1}, which is more than the maximum allowed {2}"
                .format(self.name, new_value, self.max_value))

            raise SettingOutOfBounds()

        # Issue a warning if there is an auxiliary variable, as the setting does not have any effect
        if self.has_auxiliary_variable:

            with warnings.catch_warnings():

                warnings.simplefilter("always", RuntimeWarning)

                log.warning(
                    "You are trying to assign to a parameter which is either linked or "
                    "has auxiliary variables. The assignment has no effect.")

        # Save the value as a pure floating point to avoid the overhead of the astropy.units machinery when
        # not needed

        if self._transformation is None:

            new_internal_value = new_value

        else:

            new_internal_value = self._transformation.forward(new_value)

        # If the parameter has changed, update its value and call the callbacks if needed

        if new_internal_value != self._internal_value:

            # Update
            self._internal_value = new_internal_value

            # Call the callbacks (if any)
            for callback in self._callbacks:

                try:

                    callback(self)

                except:

                    log.exception("Could not call callback for parameter %s" %
                                  self.name)

                    raise NotCallableOrErrorInCall()

    value = property(
        _get_value,
        _set_value,
        doc=
        "Get and sets the current value for the parameter, with or without units",
    )

    def _get_internal_value(self) -> float:
        """
        This is supposed to be only used by fitting engines at the beginning to get the starting value for free
        parameters. From then on, only the _set_internal_value should be used

        :return: the internal value
        """

        # NOTE: we don't need here to deal with auxiliary variables because if one is defined, the parameter is not
        # free thus it will not be touched by the fitting engine
        if not self._aux_variable is None:

            log.error("You cannot get the internal value of a parameter which has an auxiliary "
            "variable")

            raise AssertionError()

        return self._internal_value

    def _set_internal_value(self, new_internal_value) -> None:
        """
        This is supposed to be only used by fitting engines

        :param new_internal_value: new value in internal representation
        :return: none
        """

        if new_internal_value != self._internal_value:

            self._internal_value = new_internal_value

            # Call callbacks if any

            for callback in self._callbacks:

                callback(self)

    # Define the property "min_value"
    # NOTE: the min value is always in external representation

    def _get_min_value(self):
        """Return current minimum allowed value"""

        return self._external_min_value

    @accept_quantity(float, allow_none=True)
    def _set_min_value(self, min_value) -> None:
        """Sets current minimum allowed value"""

        # Check that the min value can be transformed if a transformation is present
        if self.has_transformation():

            if min_value is not None:

                if self._transformation.is_positive:

                    assert min_value > 0.0, (
                        "The transformation %s is postive definite and the min_value was set to a negative number for %s "
                        % (type(self._transformation), self.path))

                try:

                    _ = self._transformation.forward(min_value)

                except FloatingPointError:

                    log.exception(
                        "The provided minimum %s cannot be transformed with the transformation %s which "
                        "is defined for the parameter %s" %
                        (min_value, type(self._transformation), self.path))

                    raise ValueError()
            else:

                if self._transformation.is_positive:

                    # set it by default to for the user
                    min_value = 1e-99

                    log.warning(
                        "We have set the min_value of %s to 1e-99 because there was a postive transform"
                        % self.path)

        # Store the minimum as a pure float

        self._external_min_value = min_value

        # Check that the current value of the parameter is still within the boundaries. If not, issue a warning

        if (self._external_min_value is not None
                and self.value < self._external_min_value):

            log.warning("The current value of the parameter %s (%s) "
                        "was below the new minimum %s." %
                        (self.name, self.value, self._external_min_value))

            self.value = self._external_min_value

    min_value = property(
        _get_min_value,
        _set_min_value,
        doc="Gets or sets the minimum allowed value for the parameter",
    )

    def remove_minimum(self) -> None:
        """
        Remove the minimum from this parameter (i.e., it becomes boundless in the negative direction)
        """
        self._external_min_value = None

    def _set_internal_min_value(self):

        log.exception(
            "You should never attempt to change the internal representation of the minimum"
        )

        raise NotCallableOrErrorInCall()

    def _get_internal_min_value(self) -> float:
        """
        This is supposed to be only used by fitting engines to get the minimum value in internal representation.
        It is supposed to be called only once before doing the minimization/sampling, to set the range of the parameter

        :return: minimum value in internal representation (or None if there is no minimum)
        """

        if self.min_value is None:

            # No minimum set

            return None

        else:

            # There is a minimum. If there is a transformation, use it, otherwise just return the minimum

            if self._transformation is None:

                return self._external_min_value

            else:

                return self._transformation.forward(self._external_min_value)

    # Define the property "max_value"

    def _get_max_value(self) -> float:
        """Return current maximum allowed value"""

        return self._external_max_value

    @accept_quantity(float, allow_none=True)
    def _set_max_value(self, max_value) -> None:
        """Sets current maximum allowed value"""

        self._external_max_value = max_value

        # Check that the current value of the parameter is still within the boundaries. If not, issue a warning

        if (self._external_max_value is not None
                and self.value > self._external_max_value):

            log.warning("The current value of the parameter %s (%s) "
                        "was above the new maximum %s." %
                        (self.name, self.value, self._external_max_value))
            self.value = self._external_max_value

    max_value = property(
        _get_max_value,
        _set_max_value,
        doc="Gets or sets the maximum allowed value for the parameter",
    )

    def remove_maximum(self) -> None:
        """
        Remove the maximum from this parameter (i.e., it becomes boundless in the positive direction)
        """
        self._external_max_value = None

    def _set_internal_max_value(self) -> None:

        log.exception(
            "You should never attempt to change the internal representation of the minimum"
        )

        raise NotCallableOrErrorInCall()

    def _get_internal_max_value(self) -> float:
        """
        This is supposed to be only used by fitting engines to get the maximum value in internal representation.
        It is supposed to be called only once before doing the minimization/sampling, to set the range of the parameter

        :return: maximum value in internal representation (or None if there is no minimum)
        """

        if self.max_value is None:

            # No minimum set

            return None

        else:

            # There is a minimum. If there is a transformation, use it, otherwise just return the minimum

            if self._transformation is None:

                return self._external_max_value

            else:

                return self._transformation.forward(self._external_max_value)

    def _set_bounds(self, bounds) -> None:
        """Sets the boundaries for this parameter to min_value and max_value"""

        # Use the properties so that the checks and the handling of units are made automatically

        min_value, max_value = bounds

        # Remove old boundaries to avoid problems with the new one, if the current value was within the old boundaries
        # but is not within the new ones (it will then be adjusted automatically later)
        self.min_value = None
        self.max_value = None

        self.min_value = min_value

        self.max_value = max_value

    def _get_bounds(self) -> Tuple[float]:
        """Returns the current boundaries for the parameter"""

        return self.min_value, self.max_value

    bounds = property(
        _get_bounds,
        _set_bounds,
        doc="Gets or sets the boundaries (minimum and maximum) for this "
        "parameter",
    )

    def add_callback(self, callback) -> None:
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

    def empty_callbacks(self) -> None:
        """Remove all callbacks for this parameter"""
        self._callbacks = []

    def duplicate(self) -> "Parameter":
        """
        Returns an exact copy of the current parameter
        """

        # Deep copy everything to make sure that there are no ties between the new instance and the old one

        new_parameter = copy.deepcopy(self)

        return new_parameter

    def to_dict(self, minimal=False) -> Dict[str, Any]:
        """Returns the representation for serialization"""

        data = collections.OrderedDict()

        if minimal:

            # In the minimal representation we just output the value

            data["value"] = self._to_python_type(self.value)

        else:

            # In the complete representation we output everything is needed to re-build the object

            data["value"] = self._to_python_type(self.value)
            data["desc"] = str(self.description)
            data["min_value"] = self._to_python_type(self.min_value)
            data["max_value"] = self._to_python_type(self.max_value)
            # We use our own thread-safe format for the unit
            data["unit"] = self.unit.to_string(format="threadsafe")

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
    :param is_normalization: True or False, wether the parameter is a normalization or not (default: False)
    :param transformation: a transformation to be used between external value (the value the user interacts with) and
        the value the fitting/sampling engine interacts with (internal value). It is an instance of a class implementing a
        forward(external_value) and a backward(internal_value) method returning respectively the transformation of the
        external value in the internal value and viceversa. This is useful because for example the logarithm of a parameter
        with a large range of possible values (say from 1e-12 to 1e20) is handled much better by fitting engines than the
        raw value. The log transformation indeed makes the gradient much easier to compute.
    """

    def __init__(
        self,
        name: Optional[str]=None,
        value: Optional[float] = None,
        min_value: Optional[float]=None,
        max_value: Optional[float]=None,
        delta: Optional[float]=None,
        desc: Optional[str]=None,
        free: bool=True,
        unit="",
        prior=None,
        is_normalization: bool=False,
        transformation: Optional[ParameterTransformation]=None,
    ):

        # This extends ParameterBase by adding the possibility for free/fix, and a delta for fitting purposes, as
        # well as a prior

        super(Parameter, self).__init__(
            name,
            value,
            min_value=min_value,
            max_value=max_value,
            desc=desc,
            unit=unit,
            transformation=transformation,
        )

        self._free: bool = bool(free)

        # Set delta if provided, otherwise use default

        if delta is not None:

            self._delta: float = delta

        else:

            # Default is 10% of the value, unless the value is zero, in which case the delta is 0.1

            if self.value == 0:

                self._delta = 0.1

            else:

                self._delta = abs(0.1 * self.value)

        # pre-defined prior is no prior
        self._prior = None

        # override the default if a prior has been given

        if prior is not None:

            # Use the property on purpose, so all checks and setups are applied

            self.prior = prior

        # Now perform a very lazy check that we can perform math operations on the delta

        if not _behaves_like_a_number(self._delta):
            log.exception("The provided delta is not a number")
            raise TypeError()

        # Create a backup copy of the status of the parameter (useful when an auxiliary variable
        # is removed)
        self._old_free = self._free

        self._is_normalization = bool(is_normalization)

    @property
    def is_normalization(self) -> bool:

        return self._is_normalization

    # Define the property "delta"

    def _get_delta(self):
        return self._delta

    @accept_quantity(float, allow_none=False)
    def _set_delta(self, delta):
        """Sets the current delta for the parameter. The delta is used as initial step by some fitting engines"""

        self._delta = delta

    delta = property(
        _get_delta, _set_delta, doc="""Gets or sets the delta for the parameter"""
    )

    def _get_internal_delta(self):
        """
        This is only supposed to be used by fitting/sampling engine, to get the initial step in internal representation

        :return: initial delta in internal representation
        """

        if self._transformation is None:

            return self._delta

        else:

            delta_int = None

            for i in range(2):

                # Try using the low bound

                low_bound_ext = self.value - self.delta

                # Make sure we are within the margins

                if low_bound_ext > self.min_value:

                    # Ok, let's use that for the delta
                    low_bound_int = self._transformation.forward(low_bound_ext)

                    delta_int = abs(low_bound_int - self._get_internal_value())

                    break

                else:

                    # Nope, try with the hi bound

                    hi_bound_ext = self.value + self._delta

                    if hi_bound_ext < self.max_value:

                        # Ok, let's use it
                        hi_bound_int = self._transformation.forward(hi_bound_ext)
                        delta_int = abs(hi_bound_int - self._get_internal_value())

                        break

                    else:

                        # Fix delta
                        self.delta = abs(self.value - self.min_value) / 4.0

                        if self.delta == 0:

                            # Parameter at the minimum
                            self.delta = abs(self.value - self.max_value) / 4.0

                        # Try again
                        continue

            assert delta_int is not None, "Bug"

            return delta_int

    # Define the property "prior"

    def _get_prior(self):

        return self._prior

    def _set_prior(self, prior):
        """Set prior for this parameter. The prior must be a function accepting the current value of the parameter
        as input and giving the probability density as output."""

        if prior is None:

            # Removing prior

            self._prior = None

        else:

            # Try and call the prior with the current value of the parameter
            try:

                _ = prior(self.value)

            except:

                log.exception(
                    "Could not call the provided prior. "
                    + "Is it a function accepting the current value of the parameter?"
                )
                raise NotCallableOrErrorInCall()

            try:

                prior.set_units(self.unit, u.dimensionless_unscaled)

            except AttributeError:

                log.exception(
                    "It looks like the provided prior is not a astromodels function."
                )

                raise NotCallableOrErrorInCall()

            self._prior = prior

    prior = property(
        _get_prior,
        _set_prior,
        doc="Gets or sets the current prior for this parameter. The prior must be a callable function "
        "accepting the current value of the parameter as input and returning the probability "
        "density as output. Set to None to remove prior.",
    )

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

            log.error("Parameter %s does not have a defined minimum. Set one first, then re-run "
                "set_uninformative_prior" % self.path)
            
            raise ParameterMustHaveBounds( )

        else:

            try:

                prior_instance.lower_bound = self.min_value

            except SettingOutOfBounds:

                log.error( "Cannot use minimum of %s for prior %s"
                    % (self.min_value, prior_instance.name)
                )
                
                raise SettingOutOfBounds( )

        if self.max_value is None:

            log.error( "Parameter %s does not have a defined maximum. Set one first, then re-run "
                "set_uninformative_prior" % self.path)
            
            raise ParameterMustHaveBounds( )

        else:  # pragma: no cover

            try:

                prior_instance.upper_bound = self.max_value

            except SettingOutOfBounds:

                log.error("Cannot use maximum of %s for prior %s"
                    % (self.max_value, prior_instance.name))
                
                raise SettingOutOfBounds( )

        if not np.isfinite(prior_instance.upper_bound.value):

            log.error(
            "The parameter %s must have a finite maximum" % self.name
        )

            raise AssertionError()
            
        if not np.isfinite(prior_instance.lower_bound.value):

            log.error(
            "The parameter %s must have a finite minimum" % self.name
        )

            raise AssertionError()
            
        self._set_prior(prior_instance)

    # Define property "free"

    def _set_free(self, value=True):

        self._free = value

    def _get_free(self):

        return self._free

    free = property(
        _get_free,
        _set_free,
        doc="Gets or sets whether the parameter is free or not. Use booleans, like: 'p.free = True' "
        " or 'p.free = False'. ",
    )

    # Define property "fix"

    def _set_fix(self, value=True):

        self._free = not value

    def _get_fix(self):

        return not self._free

    fix = property(
        _get_fix,
        _set_fix,
        doc="Gets or sets whether the parameter is fixed or not. Use booleans, like: 'p.fix = True' "
        " or 'p.fix = False'. ",
    )

    def add_auxiliary_variable(self, variable, law) -> None:
        """TODO describe function

        :param variable: 
        :type variable: 
        :param law: 
        :type law: 
        :returns: 

        """
        
        # Assign units to the law
        law.set_units(variable.unit, self.unit)

        # Test the law
        try:

            _ = law(variable.value)

        except:  # pragma: no cover

            log.exception("The provided law for the auxiliary variable failed on call")

            raise NotCallableOrErrorInCall()

        self._aux_variable = dict(law=law, variable=variable)
        
        # Now add the law as an attribute
        # so the user will be able to access its parameters as this.name.parameter_name

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

        if not self.has_auxiliary_variable:

            # do nothing, but print a warning

            log.warning("Cannot remove a non-existing auxiliary variable")

        else:

            # Remove the law from the children

            self._remove_child(self._aux_variable["law"].name)

            # Clean up the dictionary

            self._aux_variable = None

            # Set the parameter to the status it has before the auxiliary variable was created

            self.free = self._old_free

    @property
    def has_auxiliary_variable(self) -> bool:
        """
        Returns whether the parameter is linked to an auxiliary variable
        """
        return  self._aux_variable is not None

    @property
    def auxiliary_variable(self) -> Tuple:
        """
        Returns a tuple with the auxiliary variable and the law

        :return: tuple (variable, law)
        """
        return self._aux_variable["variable"], self._aux_variable["law"]

    def _repr__base(self, rich_output=False):

        if not self.has_auxiliary_variable:

            representation = (
                "Parameter %s = %s [%s]\n"
                "(min_value = %s, max_value = %s, delta = %s, free = %s)"
                % (
                    self.name,
                    self.value,
                    self.unit,
                    self.min_value,
                    self.max_value,
                    self.delta,
                    self.free,
                )
            )

            if self._prior is not None:

                representation += " [prior: %s]" % self.prior.name

        else:

            representation = (
                "Parameter %s = %s [%s]\n"
                "(linked to auxiliary variable '%s' with law '%s')"
                % (
                    self.name,
                    self.value,
                    self.unit,
                    self._aux_variable["variable"].name,
                    self._aux_variable["law"].name,
                )
            )

        return representation

    def to_dict(self, minimal=False):

        """Returns the representation for serialization"""

        data = super(Parameter, self).to_dict()

        # Add wether is a normalization or not
        data["is_normalization"] = self._is_normalization

        if minimal:

            # No need to add anything
            pass

        else:

            # In the complete representation we output everything is needed to re-build the object

            if self.has_auxiliary_variable:

                # Store the function and the auxiliary variable

                data["value"] = "f(%s)" % self._aux_variable["variable"]._get_path()

                aux_variable_law_data = collections.OrderedDict()
                aux_variable_law_data[
                    self._aux_variable["law"].name
                ] = self._aux_variable["law"].to_dict()

                data["law"] = aux_variable_law_data

            # delta and free are attributes of Parameter, but not of ParameterBase

            data["delta"] = self._to_python_type(self._delta)
            data["free"] = self.free

            if self.has_prior():

                data["prior"] = {self.prior.name: self.prior.to_dict()}

        return data

    def get_randomized_value(self, variance=0.1):

        # Get a value close to the current value, but not identical
        # (used for the inizialization of Bayesian samplers)

        min_value = self.min_value
        max_value = self.max_value
        value = self.value

        if (min_value is not None) or (max_value is not None):

            # If _value is zero, then std will be zero, which doesn't make sense
            assert value != 0, (
                "You cannot randomize parameter %s because its value is exactly zero"
                % self.path
            )

            # Bounded parameter. Use a truncated normal so we are guaranteed
            # to have a random value within the boundaries

            std = abs(variance * value)

            if min_value is not None:

                a = (min_value - value) / std

            else:

                a = -np.inf

            if max_value is not None:

                b = (max_value - value) / std

            else:

                b = np.inf

            sample = scipy.stats.truncnorm.rvs(a, b, loc=value, scale=std, size=1)

            if (min_value is not None and sample < min_value) or (
                max_value is not None and sample > max_value
            ):  # pragma: no cover

                # This should never happen

                raise AssertionError(
                    "Got a sample outside of the boundaries of the truncated normal distribution"
                )

            return sample[0]

        else:

            # The parameter has no boundaries

            return np.random.normal(value, abs(variance * value))


class IndependentVariable(ParameterBase):
    """
    An independent variable like time or energy.
    """

    # Override the constructor to make the unit specification mandatory

    def __init__(self, name, value, unit, min_value=None, max_value=None, desc=None):

        super(IndependentVariable, self).__init__(
            name, value, unit=unit, min_value=min_value, max_value=max_value, desc=desc
        )

    def _repr__base(self, rich_output=False):

        return "IndependentVariable %s = %s\n" "(min_value = %s, max_value = %s)" % (
            self.name,
            self.value,
            self.min_value,
            self.max_value,
        )
