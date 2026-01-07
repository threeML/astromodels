import pytest
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
    cplaw_inverse_eval,
    super_cplaw_eval,
    band_eval,
    bplaw_eval,
    sbplaw_eval,
    bb_eval,
    mbb_eval,
    ggrb_int_pl,
    non_diss_photoshere_generic,
    dbl_sbpl,
)
import numpy as np


def test_vectorized_functions():
    assert (_expm1(np.array([0, 0, 0])) == np.array([0, 0, 0])).all
    assert (_exp(np.array([0, 0, 0])) == np.array([1, 1, 1])).all
    assert (_sqrt(np.array([4, 4, 4])) == np.array([2, 2, 2])).all
    assert (_pow(np.array([2, 2, 2]), 2) == np.array([4, 4, 4])).all
    assert (_log(np.array([1, 1, 1])) == np.array([0, 0, 0])).all
    assert (_log10(np.array([1, 1, 1])) == np.array([0, 0, 0])).all
