from __future__ import print_function
from __future__ import division
import pytest
import astropy.units as u

try:

    from astromodels.xspec import *

except:

    has_XSPEC = False

else:

    has_XSPEC = True


# This defines a decorator which can be applied to single tests to
# skip them if the condition is not met
skip_if_xspec_is_not_available = pytest.mark.skipif(not has_XSPEC,
                                                    reason="XSPEC not available")


@skip_if_xspec_is_not_available
def test_xspec_load():

    # no need to do anything really
    s = XS_phabs() * XS_powerlaw() + XS_bbody()
    print(s(1.0))
    s.set_units(u.keV, 1 / (u.keV * u.cm**2 * u.s))
    print(s(1.0 * u.keV)) 
