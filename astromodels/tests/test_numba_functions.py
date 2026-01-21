from astromodels.functions.numba_functions import (
    _expm1,
    _exp,
    _sqrt,
    _pow,
    _log,
    _log10,
    plaw_eval,
    plaw_flux_norm,
    cplaw_eval,
)
import numpy as np


def test_vectorized_functions():
    assert (_expm1(np.array([0, 0, 0])) == np.array([0, 0, 0])).all
    assert (_exp(np.array([0, 0, 0])) == np.array([1, 1, 1])).all
    assert (_sqrt(np.array([4, 4, 4])) == np.array([2, 2, 2])).all
    assert (_pow(np.array([2, 2, 2]), 2) == np.array([4, 4, 4])).all
    assert (_log(np.array([1, 1, 1])) == np.array([0, 0, 0])).all
    assert (_log10(np.array([1, 1, 1])) == np.array([0, 0, 0])).all


def test_jitted_functions():
    x = np.array([0.1, 1, 10, 100])
    assert (plaw_eval.py_func(x, 1, -2, 1) == np.array([100, 1, 0.01, 1e-4])).all
    assert np.isclose(plaw_flux_norm.py_func(-2.0, 1, 10), 2.302585)
    assert np.isclose(plaw_flux_norm.py_func(-4.0, 1, 10), 0.495)
    assert np.isclose(
        cplaw_eval.py_func(x, 1, 100, -2, 1),
        np.array([99.90005, 0.99005, 9.095 * 1e-3, 3.6788 * 1e-5]),
    ).all
