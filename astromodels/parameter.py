__author__ = 'giacomov'

"""
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
callbacks use the method empty_callbacks().

.. _parameter_auxvar:
===================
Auxiliary variable
===================

Sometimes it is useful to specify the parameter as a function of another variable ("auxiliary variable"). For
example, we want the parameter to vary with time according to a linear law. This can be done as follows.

First define the auxiliary variable "time" starting at 0, and a linear law 3.2 t + 5.6::

  > t = AuxiliaryVariable("time",0.0)
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

If we now change the value of one of the parameters, we see the parameter p changing accordingly:
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


import warnings
import exceptions


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

class Parameter(object):
    """

    Implements a numerical parameter. Optionally the parameter can vary according to an auxiliary law (see below).

    :param name: Name for the parameter
    :param initial_value: Initial value
    :keyword min: minimum allowed value for the parameter (default: None)
    :keyword max: maximum allowed value for the parameter (default: None)
    :keyword delta: initial step used by some fitting engines (default: 0.1 * value)

    """

    def __init__(self, name, initial_value, **kwargs):

        # NOTE: we avoid enforcing a particular type or even that initial_values, max, min and delta are python numbers.
        # Indeed, as long as they behave as numbers we are going to be fine. We want to keep the possibility to use
        # numbers coming from C, C++ or other sources, for example. Note however that string
        # parameters are not supported

        # Callbacks are executed any time the value for the parameter changes (i.e., its value changes)

        # We start from a empty list of callbacks.
        self._callbacks = []

        # Assign to members

        self.__name = name

        self._value = initial_value

        # Set minimum if provided, otherwise use default

        if 'min' in kwargs.keys():

            self._min_value = kwargs['min']

        else:

            self._min_value = None

        # Set maximum if provided, otherwise use default

        if 'max' in kwargs.keys():

            self._max_value = kwargs['max']

        else:

            self._max_value = None

        # Set delta if provided, otherwise use default

        if 'delta' in kwargs.keys():

            self._delta = kwargs['delta']

        else:

            self._delta = 0.1 * self._value

        # pre-defined prior is no prior
        self._prior = None

        # by default we have no auxiliary variable

        self._auxVariable = {}

        # Now perform a very lazy check that we can perform math operations on all members

        if not _behaves_like_a_number(self._value):
            raise TypeError("The provided initial value is not a number")

        if not _behaves_like_a_number(self._delta):
            raise TypeError("The provided delta is not a number")

        if self._min_value is not None:

            if not _behaves_like_a_number(self._min_value):
                raise TypeError("The provided minimum value is not a number")

        if self._max_value is not None:

            if not _behaves_like_a_number(self._max_value):
                raise TypeError("The provided maximum value is not a number")

    # Define the property "name" but make it read-only

    @property
    def name(self):
        return self.__name

    # Define the property "value" with a control that the parameter cannot be set
    # outside of its bounds

    def get_value(self):
        """Return current parameter value"""

        # If this parameter has a law of variation, use it to compute the current value

        # (remember: an empty dictionary test as False, a non-empty as true)

        if self._auxVariable:

            self.value = self._auxVariable['law'](self._auxVariable['variable'].value)

        return self._value

    def set_value(self, value):
        """Sets the current value of the parameter, ensuring that it is within the allowed range"""

        if self._min_value is not None and value < self._min_value:
            raise SettingOutOfBounds(
                "Trying to set parameter {0} = {1}, which is less than the minimum allowed {2}".format(
                    self.name, value, self._min_value))

        if self._max_value is not None and value > self._max_value:
            raise SettingOutOfBounds(
                "Trying to set parameter {0} = {1}, which is more than the maximum allowed {2}".format(
                    self.name, value, self._max_value))

        # Call the callbacks (if any)

        for callback in self._callbacks:

            try:

                callback(value)

            except:

                raise NotCallableOrErrorInCall("Could not get callback for parameter %s with value %s" % (self.name, value))

        self._value = value

    value = property(get_value, set_value,
                     doc="""Gets or sets the value for the parameter""")

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

    # Define the property "min_value"

    def get_min_value(self):
        """Return current minimum allowed value"""
        return self._min_value

    def set_min_value(self, min_value):
        """Sets current minimum allowed value"""

        if min_value is not None and not _behaves_like_a_number(min_value):
            raise ValueError("Provided minimum value is not a number nor None")

        self._min_value = min_value

        # Check that the current value of the parameter is still within the boundaries. If not, issue a warning

        if self.value < self.min_value:
            warnings.warn("The current value of the parameter is below the new minimum.", exceptions.RuntimeWarning)

    min_value = property(get_min_value, set_min_value,
                         doc='Gets or sets the minimum allowed value for the parameter')

    # Define the property "max_value"

    def get_max_value(self):

        """Return current minimum allowed value"""

        return self._max_value

    def set_max_value(self, max_value):

        """Sets current minimum allowed value"""

        if max_value is not None and not _behaves_like_a_number(max_value):
            raise ValueError("Provided maximum value is not a number nor None")

        self._max_value = max_value

        # Check that the current value of the parameter is still within the boundaries. If not, issue a warning

        if self.value > self.max_value:
            warnings.warn("The current value of the parameter is above the new maximum.", exceptions.RuntimeWarning)

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

    def add_callback(self, callback):
        """Add a callback to the list of functions which are called whenever the value of the parameter is changed.
        The callback must be a function accepting the current value as input. The return value of the callback is
        ignored. More than one callback can be specified. In that case, the callbacks will be called in the same order
        they have been entered."""

        self._callbacks.append(callback)

    def empty_callbacks(self):
        """Remove all callbacks for this parameter"""
        self._callbacks = []

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
                     doc="Gets or sets the current prior for this parameter. The prior must be a callable function " +
                         "accepting the current value of the parameter as input and returning the probability " +
                         "density as output")

    def add_auxiliary_variable(self, variable, law):

        assert isinstance(variable, AuxiliaryVariable)

        #Test the law
        try:

            this_value = law(variable.value)

        except:

            raise NotCallableOrErrorInCall("The provided law for the auxiliary variable failed on call")


        self._auxVariable['law'] = law
        self._auxVariable['variable'] = variable

        # Set the value of the parameter to the current value of the law

        self.value = this_value

class AuxiliaryVariable(Parameter):
    """
    A variable which can be used to express a parameter as a function of it. This is at the moment essentially another
    name of the :Parameter: class, used for clarity to differentiate between a parameter and a variable. See the
    documentation of the :Parameter: class (section :ref:`parameter_auxvar`) for an example of use.
    """

    pass