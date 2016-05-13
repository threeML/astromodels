import pytest

import astropy.units as u
import numpy as np

__author__ = 'giacomov'

from astromodels.parameter import Parameter, SettingOutOfBounds, IndependentVariable
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

    with pytest.raises(RuntimeError):

        _ = p.prior

    assert p.unit == u.dimensionless_unscaled

    # Test that we cannot call a parameter with a name with spaces in it
    with pytest.raises(AssertionError):

        _ = Parameter('test parameter 2', 1.0)


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

    with pytest.raises(RuntimeError):

        _ = p.prior

    assert p.unit == u.keV


def test_constructor_complete():

    p = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test', free=False, unit=u.MeV)

    assert p.min_value == -5.0
    assert p.max_value == 5.0
    assert p.value == 1.0
    assert p.delta == 0.2
    assert p.name == 'test_parameter'
    assert p.description == 'test'
    assert p.fix == True
    assert p.free == False
    assert p.has_prior() == False

    with pytest.raises(RuntimeError):
        _ = p.prior

    assert p.unit == u.MeV


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

    with pytest.raises(RuntimeError):
        _ = p.prior

    assert p.unit == u.MeV


def test_constructor_with_boundaries():

    p = Parameter('test_parameter', 1.0, min_value=-5, max_value=5)

    assert p.min_value == -5
    assert p.max_value == 5


def test_constructor_with_delta():

    p = Parameter('test_parameter', 1.0, delta=0.3)

    assert p.delta == 0.3


def test_constructor_with_units():

    p = Parameter('test_parameter', 1.0, unit=u.keV)

    assert p.unit == u.keV


def test_set_no_units():

    p = Parameter('test_parameter',1.0)

    p.value = 25.4

    assert p.value == 25.4


def test_set_within_bounds_no_units():

    p = Parameter('test_parameter',1.0, min_value = -2.0, max_value = 2.0)

    p.value = 1.5

    assert p.value == 1.5


def test_set_outside_bounds_no_units():

    p = Parameter('test_parameter',1.0, min_value = -2.0, max_value = 2.0)

    with pytest.raises(SettingOutOfBounds):

        p.value = -10.0


def test_set_units():

    p = Parameter('test_parameter',1.0, unit=u.keV)

    p.value = 3.0 * u.MeV

    assert p.value == 3000.0


def test_set_within_bounds_units():

    p = Parameter('test_parameter',1.0 * u.keV, min_value = -2.0 * u.MeV, max_value = 2.0 * u.MeV, unit=u.keV)

    p.value = 1.2 * u.MeV

    assert p.value == 1200.0


def test_set_outside_bounds_units():

    p = Parameter('test_parameter', 1.0 * u.keV, min_value = -2.0 * u.MeV, max_value = 2.0 * u.MeV, unit=u.keV)

    with pytest.raises(SettingOutOfBounds):

        p.value = -10.0 * u.MeV


def test_set_bounds_nounits():

    p = Parameter('test_parameter', 1.0)

    p.bounds = (-2.0 ,2.0)

    assert p.min_value == -2.0
    assert p.max_value == 2.0


def test_set_bounds_units():

    p = Parameter('test_parameter', 1.0 * u.keV)

    p.bounds = (-2.0 * u.MeV, 2.0 * u.eV)

    assert p.min_value == -2000
    assert p.max_value == 2e-3


def test_set_delta_nounits():

    p = Parameter('test_parameter', 1.0)

    p.delta = 0.5

    assert p.delta == 0.5


def test_set_delta_units():

    p = Parameter('test_parameter', 1.0, unit='GeV')

    p.delta = 500 * u.MeV

    assert p.delta == 0.5


def test_duplicate():

    p1 = Parameter('test_parameter', 1.0, min_value=-5.0, max_value=5.0, delta=0.2, desc='test', free=False, unit='MeV')

    p2 = p1.duplicate()

    assert p1.to_dict() == p2.to_dict()


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

    x.value = 4.0

    assert p1.value == 6.0

    # Check that assigning to the parameter doesn't produce any effect
    p1.value = -1.0

    assert p1.value == 6.0


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