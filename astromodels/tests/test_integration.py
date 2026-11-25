import astropy.units as u
import warnings
import numpy as np
import pytest
from astromodels.functions.function import _known_functions, Function1D


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


def test_analytical_integral_vectorized():
    pl = Powerlaw()
    pl.index.value = -2
    assert np.all(
        np.isclose(pl.integrate(np.ones(10) * 0.1, np.ones(10)), np.ones(10) * 9)
    )
    assert pl.integral_numerical_error is None

    pl = Powerlaw()
    pl.index.value = -1
    assert np.all(
        np.isclose(
            pl.integrate(np.ones(10) * 0.1, np.ones(10)), np.ones(10) * 2.302585093
        )
    )
    assert pl.integral_numerical_error is None

    c = Constant()
    c.k.value = 1
    assert np.all(np.isclose(c.integrate(np.zeros(11), np.ones(11)), np.ones(11)))

    for order in range(2, 5, 1):
        pol = get_polynomial(order)
        assert np.isclose(
            pol.integrate(0, 1), np.sum([1 / (i + 1) for i in range(order + 1)])
        )


def test_integrate_function1D():
    for k, v in _known_functions.items():
        # exclude the functions that need external input for now
        if issubclass(v, Function1D):
            if k in [
                "GenericFunction",
                "_ComplexTestFunction",
                "TemplateModel",
            ] or (
                k.startswith("XS_") and k not in ["XS_agauss"]
            ):  # do not test the XSPEC functions
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # test that the function exists
                _ = v().integrate(1, 2)


def test_numerical_integral():
    line = Line()
    line.a.value = 0
    line.b.value = 1
    assert np.isclose(line.integrate(0, 1, use_scipy=True), 0.5)
    assert line.integral_numerical_error is not None
    assert np.isclose(line.integrate(1, 0, use_scipy=True), -0.5)

    lower = np.array([0.0, 1.0, 2.0])
    upper = np.array([1.0, 2.0, 0.0])
    assert np.allclose(
        line.integrate(lower, upper, use_scipy=True), np.array([0.5, 1.5, -2.0])
    )
    assert np.shape(line.integral_numerical_error) == lower.shape


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
