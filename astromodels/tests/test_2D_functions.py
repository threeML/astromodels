from astromodels.functions import Latitude_galactic_diffuse
from astropy.coordinates import SkyCoord
import numpy as np
import astropy.units as u
import pytest


def test_latitude_galactic_diffuse():

    lgd = Latitude_galactic_diffuse()
    lgd.set_frame("galactic")
    assert np.isclose(lgd(1, 0), 1)
    with pytest.raises(NotImplementedError):
        lgd.evaluate(0 * u.deg, 1, 1, 1, 1, 1)

    assert np.isclose(
        lgd(SkyCoord(ra=266.40498829, dec=-28.93617776, unit="deg", frame="icrs")), 1
    )
    lgd.set_frame("icrs")
    lgd.set_units(u.deg, u.deg, u.Unit("deg-2"))

    print(lgd)
    assert np.isclose(
        lgd(266.40498829 * u.deg, -28.93617776 * u.deg), 1 * u.Unit("deg-2")
    )

    assert np.all(
        np.isclose(
            lgd(
                SkyCoord(
                    ra=[266.40498829] * 10,
                    dec=[-28.93617776] * 10,
                    unit="deg",
                    frame="icrs",
                )
            ),
            1,
        )
    )
