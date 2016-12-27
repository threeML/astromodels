import pytest

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
    pass
