import pytest
import numpy as np

try:

    from astromodels.xspec import XS_bbody
    from astromodels.xspec.xspec_settings import xspec_abund, xspec_cosmo, xspec_xsect

except (ImportError, ModuleNotFoundError):

    has_XSPEC = False

else:
    XS_bbody()
    has_XSPEC = True


# This defines a decorator which can be applied to single tests to
# skip them if the condition is not met
skip_if_xspec_is_not_available = pytest.mark.skipif(
    not has_XSPEC, reason="XSPEC not available"
)


@skip_if_xspec_is_not_available
def test_xspec_abund():
    xspec_abund("lodd")
    current_abund = xspec_abund()
    assert current_abund == "lodd", "Setting abundances failed"


@skip_if_xspec_is_not_available
def test_xspec_cosmo():
    set1 = [80, 0.1, 0.73]
    set2 = [80, None, 0.73]
    xspec_cosmo(*set1)
    assert np.isclose(xspec_cosmo(), set1).all, "Setting all failed"
    xspec_cosmo(*set2)
    assert np.isclose(xspec_cosmo(), [80, 0.1, 0.73]).all, "Setting parts failed"


@skip_if_xspec_is_not_available
def test_xspec_xsect():
    xspec_xsect("vern")
    assert xspec_xsect() == "vern"
