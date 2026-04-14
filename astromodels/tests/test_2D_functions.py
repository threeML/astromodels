from astromodels.functions import Latitude_galactic_diffuse
import numpy as np
import astropy.units as u
import pytest


def test_latitude_galactic_diffuse():

    lgd = Latitude_galactic_diffuse()
    lgd.set_frame("galactic")
    assert np.isclose(lgd(1, 0), 1)
    with pytest.raises(NotImplementedError):
        lgd.evaluate(0 * u.deg, 1, 1, 1, 1, 1)
