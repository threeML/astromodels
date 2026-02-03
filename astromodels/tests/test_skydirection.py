import pytest
import numpy as np
from astromodels.core.sky_direction import SkyDirection
from astromodels.core.units import get_units
from astropy.coordinates import SkyCoord
import astropy.units as u


def test_from_skycoord():

    skycoord = SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
    skydirection = SkyDirection(position=skycoord)

    assert np.isclose(skydirection.get_ra() * get_units().angle, 10 * u.deg)
    assert np.isclose(skydirection.get_dec() * get_units().angle, 20 * u.deg)

    skycoord = SkyCoord(ra=1 * u.rad, dec=0.2 * u.rad)
    skydirection = SkyDirection(position=skycoord)

    assert np.isclose(skydirection.get_ra() * get_units().angle, 1 * u.rad)
    assert np.isclose(skydirection.get_dec() * get_units().angle, 0.2 * u.rad)

    with pytest.raises(AssertionError):
        SkyDirection(ra=0, dec=0, l=0, b=0)

    sd = SkyDirection(ra=10, dec=20, unit="deg")
    assert np.isclose(sd.get_ra() * get_units().angle, 10 * u.deg)
    assert np.isclose(sd.get_dec() * get_units().angle, 20 * u.deg)

    skydirection = SkyDirection(l=0, b=0, unit="rad")

    assert np.isclose(skydirection.get_l() * get_units().angle, 0 * u.rad)
    assert np.isclose(skydirection.get_b() * get_units().angle, 0 * u.rad)
