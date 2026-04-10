from astromodels.functions import Powerlaw, Line, Constant, get_polynomial
import numpy as np


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
