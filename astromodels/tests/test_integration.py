import astropy.units as u
import numpy as np
import pytest

from astromodels.functions import Constant, Line, Powerlaw, get_polynomial


def test_analytical_integral():
    pl = Powerlaw()
    pl.index.value = -2
    assert np.isclose(pl.integrate(0.1, 1), 9)
    assert pl.integral_numerical_error is None
    c = Constant()
    c.k.value = 1
    assert np.isclose(c.integrate(0, 1), 1)

    for order in range(2, 5, 1):
        pol = get_polynomial(order)
        assert np.isclose(
            pol.integrate(0, 1), np.sum([1 / (i + 1) for i in range(order + 1)])
        )


def test_numerical_integral():
    line = Line()
    line.a.value = 0
    line.b.value = 1
    assert np.isclose(line.integrate(0, 1, use_scipy=True), 0.5)
    assert line.integral_numerical_error is not None
    assert np.isclose(line.integrate(1, 0, use_scipy=True), -0.5)


def test_quantities():
    line = Line()
    line.set_units(u.keV, u.Unit("keV-1 cm-2 s-1"))
    line.a.value = 0
    line.b.value = 1
    int_val = line.integrate(0 * u.keV, 1 * u.keV)
    assert np.isclose(int_val.value, 0.5)
    assert int_val.unit == u.Unit("cm-2 s-1")

    with pytest.raises(ValueError):
        _ = line.integrate(1 * u.keV, 1 * u.m)
    with pytest.raises(TypeError):
        _ = line.integrate(1 * u.keV, 1)
