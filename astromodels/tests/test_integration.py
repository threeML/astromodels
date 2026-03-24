from astromodels.functions import Powerlaw, Line
import numpy as np


def test_analytical_integral():
    pl = Powerlaw()
    pl.index.value = -2
    assert np.isclose(pl.integrate(0.1, 1), 9)


def test_numerical_integral():
    line = Line()
    line.a.value = 0
    line.b.value = 1
    assert np.isclose(line.integrate(0, 1), 0.5)
