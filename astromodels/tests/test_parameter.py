import astropy.units as u
import pytest
from astromodels.functions.functions import Uniform_prior, Log_uniform_prior, Powerlaw

__author__ = 'giacomov'

from astromodels.core.parameter import Parameter, SettingOutOfBounds, \
    CannotConvertValueToNewUnits, NotCallableOrErrorInCall, IndependentVariable, ParameterMustHaveBounds
from astromodels.functions.functions import Line


def test_default_constructor():

    p = Parameter('test_parameter', 1.0, desc='Description')

    assert p.min_value is None
    assert p.max_value is None
    assert p.value == 1.0
    assert isinstance(p.delta, float)
    assert p.name == 'test_parameter'
    assert p.description == 'Description'
    assert p.fix == False
    assert p.free == True
    assert p.has_prior() == False
    assert p.prior is None

    assert p.unit == u.dimensionless_unscaled

    # Test that we cannot call a parameter with a name with spaces in it
    with pytest.raises(AssertionError):

        _ = Parameter('test parameter 2', 1.0)

    # Test some failures cases
    with pytest.raises(TypeError):

        _ = Parameter('test', '1.0')

    with pytest.raises(ValueError):
        _ = Parameter('test', 1.0, min_value='a')

    with pytest.raises(ValueError):

        _ = Parameter('test', 1.0, max_value='b')

    with pytest.raises(TypeError):

        _ = Parameter('test', 1.0, delta='b')

    p.display()


def test_default_constructor_units():

    p = Parameter('test_parameter', 1.0 * u.keV, desc='Description')

    assert p.min_value is None
    assert p.max_value is None
    assert p.value == 1.0
    assert isinstance(p.delta, float)
    assert p.name == 'test_parameter'
    assert p.description == 'Description'
    assert p.fix == False
    assert p.free == True
    assert p.has_prior() == False
    assert p.prior is None

    assert p.unit == u.keV

    p.display()


def test_constructor_complete():

    p = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test',
                  free=False, unit=u.MeV, prior=Uniform_prior())

    assert p.min_value == -5.0
    assert p.max_value == 5.0
    assert p.value == 1.0
    assert p.delta == 0.2
    assert p.name == 'test_parameter'
    assert p.description == 'test'
    assert p.fix == True
    assert p.free == False
    assert p.has_prior() == True

    assert p.unit == u.MeV

    p.display()


def test_conflicting_units_in_initial_value_and_unit_keyword():

    p = Parameter('test_parameter', 1.0 * u.keV, desc='Description', unit=u.MeV)

    assert p.min_value is None
    assert p.max_value is None
    assert p.value == 1.0e-3
    assert isinstance(p.delta, float)
    assert p.name == 'test_parameter'
    assert p.description == 'Description'
    assert p.fix == False
    assert p.free == True
    assert p.has_prior() == False
    assert p.prior is None

    assert p.unit == u.MeV

    p.display()


def test_constructor_with_boundaries():

    p = Parameter('test_parameter', 1.0, min_value=-5, max_value=5)

    assert p.min_value == -5
    assert p.max_value == 5

    p.display()


def test_constructor_with_delta():

    p = Parameter('test_parameter', 1.0, delta=0.3)

    assert p.delta == 0.3

    p.display()


def test_constructor_with_units():

    p = Parameter('test_parameter', 1.0, unit=u.keV)

    assert p.unit == u.keV

    p.display()


def test_set_no_units():

    p = Parameter('test_parameter',1.0)

    p.value = 25.4

    assert p.value == 25.4

    p.display()


def test_set_within_bounds_no_units():

    p = Parameter('test_parameter',1.0, min_value = -2.0, max_value = 2.0)

    p.value = 1.5

    assert p.value == 1.5

    p.display()


def test_set_outside_bounds_no_units():

    p = Parameter('test_parameter',1.0, min_value = -2.0, max_value = 2.0)

    with pytest.raises(SettingOutOfBounds):

        p.value = -10.0

    with pytest.raises(SettingOutOfBounds):

        p.value = 10.0

    p.display()


def test_set_units():

    p = Parameter('test_parameter',1.0, unit=u.keV)

    p.value = 3.0 * u.MeV

    assert p.value == 3000.0

    with pytest.raises(u.UnitConversionError):

        p.value = 3.0 * u.cm

    with pytest.raises(CannotConvertValueToNewUnits):

        p.unit = u.cm

    with pytest.raises(CannotConvertValueToNewUnits):

        p.unit = u.dimensionless_unscaled

    p.unit = u.MeV

    assert p.unit == u.MeV

    p.display()


def test_set_within_bounds_units():

    p = Parameter('test_parameter',1.0 * u.keV, min_value = -2.0 * u.MeV, max_value = 2.0 * u.MeV, unit=u.keV)

    p.value = 1.2 * u.MeV

    assert p.value == 1200.0

    p.display()


def test_set_outside_bounds_units():

    p = Parameter('test_parameter', 1.0 * u.keV, min_value = -2.0 * u.MeV, max_value = 2.0 * u.MeV, unit=u.keV)

    with pytest.raises(SettingOutOfBounds):

        p.value = -10.0 * u.MeV

    with pytest.raises(SettingOutOfBounds):

        p.value = 10.0 * u.MeV

    p.display()


def test_set_bounds_nounits():

    p = Parameter('test_parameter', 1.0)

    p.bounds = (-2.0 ,2.0)

    assert p.min_value == -2.0
    assert p.max_value == 2.0

    p.display()

    with pytest.warns(RuntimeWarning):

        p.value = 1.0
        p.min_value = 2.0


def test_set_bounds_units():

    p = Parameter('test_parameter', 1.0 * u.keV)

    p.bounds = (-2.0 * u.MeV, 2.0 * u.eV)

    assert p.min_value == -2000
    assert p.max_value == 2e-3

    p.display()


def test_set_delta_nounits():

    p = Parameter('test_parameter', 1.0)

    p.delta = 0.5

    assert p.delta == 0.5

    p.display()


def test_set_delta_units():

    p = Parameter('test_parameter', 1.0, unit='GeV')

    p.delta = 500 * u.MeV

    assert p.delta == 0.5

    p.display()


def test_duplicate():

    p1 = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test', free=False, unit='MeV')

    p2 = p1.duplicate()

    assert p1.to_dict() == p2.to_dict()

    p1.display()
    p2.display()


def test_get_randomized_value():

    # Test randomization no boundaries (normal distribution)
    p1 = Parameter('test_parameter', 1.0)

    val2 = p1.get_randomized_value(0.1)

    assert isinstance(val2, float)

    # Test the randomized value with truncated normal, i.e., with boundaries

    p2 = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test', free=False, unit='MeV')

    val1 = p2.get_randomized_value(0.1)

    assert p2.min_value <= val1 <= p2.max_value

    # Test the same but with a large variance

    val1 = p2.get_randomized_value(10.0)

    assert p2.min_value <= val1 <= p2.max_value

    p2.min_value = None

    val1 = p2.get_randomized_value(10.0)
    assert val1 <= p2.max_value

    p2.min_value = -5.0
    p2.max_value = None

    val1 = p2.get_randomized_value(10.0)
    assert val1 >= p2.min_value


def test_set_remove_minimum():
    p1 = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test', free=False, unit='MeV')

    p1.remove_minimum()

    assert p1.min_value == None

    p1.value = -1000.0

    assert p1.value == -1000.0


def test_set_remove_maximum():

    p1 = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test', free=False, unit='MeV')

    p1.remove_maximum()

    assert p1.max_value == None

    p1.value = 1000.0

    assert p1.value == 1000.0


def test_set_auxiliary_variable():

    p1 = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test', free=False, unit='MeV')

    x = Parameter('aux_variable', 1.0)

    # ax + b

    law = Line()
    law.a = 1.0
    law.b = 2.0

    p1.add_auxiliary_variable(x, law)

    assert p1.has_auxiliary_variable() == True

    assert p1.value == 3.0

    assert p1.free == False

    x.value = 4.0

    assert p1.value == 6.0

    # Check that assigning to the parameter doesn't produce any effect
    p1.value = -1.0

    assert p1.value == 6.0

    # Now check errors reporting
    with pytest.raises(AttributeError):

        p1.add_auxiliary_variable(1.0, law)

    # Now add it twice to verify that it overwrites it
    p1.add_auxiliary_variable(x, law)
    p1.add_auxiliary_variable(x, law)

    p1.display()



def test_remove_auxiliary_variable():

    p1 = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test', free=False, unit='MeV')

    x = Parameter('aux_variable', 1.0)

    # ax + b

    law = Line()
    law.a = 1.0
    law.b = 2.0

    p1.add_auxiliary_variable(x, law)

    assert p1.value == 3.0

    x.value = 4.0

    assert p1.value == 6.0

    p1.remove_auxiliary_variable()

    assert p1.has_auxiliary_variable() == False

    p1.value = -1.0

    assert p1.value == -1.0

    with pytest.warns(RuntimeWarning):

        p1.remove_auxiliary_variable()


def test_callback():

    p1 = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test', free=False, unit='MeV')

    class Callback(object):

        def __init__(self):

            self._control_value = None

        def __call__(self, p):

            assert p == p1

            self._control_value = p1.value

    working_callback = Callback()

    p1.add_callback(working_callback)

    # Test the callback
    p1.value = 2.0
    assert working_callback._control_value == p1.value

    def not_working_callback():

        # Wrong calling sequence
        pass

    p1.add_callback(not_working_callback)

    with pytest.raises(NotCallableOrErrorInCall):

        p1.value = 2.0

    p1.empty_callbacks()

    assert len(p1._callbacks) == 0


def test_to_dict():

    p1 = IndependentVariable('time', 1.0, min_value=-5.0, max_value=5.0, desc='test', unit='MeV')

    repr = p1.to_dict(minimal=False)

    assert len(repr.keys()) == 5

    repr2 = p1.to_dict(minimal=True)

    assert len(repr2.keys()) == 1
    assert 'value' in repr2

    assert repr2['value'] == p1.value

    p = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0,
                   delta=0.2, desc='test', free=False, unit='MeV')

    p.to_dict()
    p.to_dict(minimal=True)

    p = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0,
                  delta=0.2, desc='test', free=False, unit='MeV', prior=Log_uniform_prior())

    p.to_dict()
    p.to_dict(minimal=True)


def test_independent_variable_representation():

    p1 = IndependentVariable('time', 1.0, min_value=-5.0, max_value=5.0, desc='test', unit='MeV')

    print(p1._repr__base(False))
    print(p1._repr__base(True))


def test_prior():

    p1 = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0,
                   delta=0.2, desc='test', free=False, unit='MeV')

    my_prior = Uniform_prior()

    p1.prior = my_prior

    assert my_prior == p1.prior

    custom_prior = lambda x:x**2

    with pytest.raises(NotCallableOrErrorInCall):

        p1.prior = custom_prior

    invalid_prior = lambda x,y: x*y

    with pytest.raises(NotCallableOrErrorInCall):

        p1.prior = invalid_prior

    # Test the set_uninformative_prior method
    p1.min_value = None
    p1.max_value = 100.0

    with pytest.raises(ParameterMustHaveBounds):

        p1.set_uninformative_prior(Uniform_prior)

    p1.min_value = 0.0
    p1.max_value = None

    with pytest.raises(ParameterMustHaveBounds):

        p1.set_uninformative_prior(Uniform_prior)

    p1.min_value = 0.0
    p1.max_value = 100.0

    p1.set_uninformative_prior(Uniform_prior)

    # Log-uniform cannot be used if minimum is 0.0

    with pytest.raises(SettingOutOfBounds):

        p1.set_uninformative_prior(Log_uniform_prior)

    p1.min_value = 1.0
    p1.set_uninformative_prior(Log_uniform_prior)


def test_as_quantity():

    p = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0,
                  delta=0.2, desc='test', free=False, unit='MeV')

    assert isinstance(p.as_quantity, u.Quantity)
    assert p.as_quantity.to("keV").value == 1000.0


def test_in_unit_of():

    p = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0,
                  delta=0.2, desc='test', free=False, unit='MeV')

    assert p.in_unit_of(u.keV) == 1000.0
    assert p.in_unit_of(u.keV, as_quantity=True).to("MeV").value == 1.0


class Callback(object):

    def __init__(self):

        self._control_value = None

    def __call__(self, p):

        self._control_value = p.value


def test_pickle():

    import cPickle

    p_orig = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test',
                       free=False, unit=u.MeV, prior=Uniform_prior())

    # Add a callback

    working_callback = Callback()

    p_orig.add_callback(working_callback)

    # Now pickle and unpickle

    d = cPickle.dumps(p_orig)

    p = cPickle.loads(d)

    # Check that everything is fine

    assert p.min_value == -5.0
    assert p.max_value == 5.0
    assert p.value == 1.0
    assert p.delta == 0.2
    assert p.name == 'test_parameter'
    assert p.description == 'test'
    assert p.fix == True
    assert p.free == False
    assert p.has_prior() == True

    assert p.unit == u.MeV

    # Test the callback
    p.value = 2.0

    callback = p.get_callbacks()[0]

    assert callback._control_value == p.value

def test_links_and_pickle():

    import cPickle

    p_orig = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test',
                       free=False, unit=u.MeV, prior=Uniform_prior())

    # Test the linkinking and pickle

    # Add a link
    x = Parameter('aux_variable', 1.0)

    # ax + b

    law = Line()
    law.a = 1.0
    law.b = 2.0

    p_orig.add_auxiliary_variable(x, law)

    # Now pickle and unpickle

    d = cPickle.dumps(p_orig)

    p = cPickle.loads(d)

    assert p.has_auxiliary_variable() == True

    assert p.value == 3.0

    assert p.free == False

    p.auxiliary_variable[0].value = 4.0

    assert p.value == 6.0

    # Check that assigning to the parameter doesn't produce any effect
    p.value = -1.0

    assert p.value == 6.0
