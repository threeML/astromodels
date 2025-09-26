from pathlib import Path

import astropy.units as u
import pytest

from astromodels import Model, PointSource, clone_model, load_model

try:

    from astromodels.xspec import XS_bbody, XS_phabs, XS_powerlaw

except (ImportError, ModuleNotFoundError):

    has_XSPEC = False

else:

    has_XSPEC = True


# This defines a decorator which can be applied to single tests to
# skip them if the condition is not met
skip_if_xspec_is_not_available = pytest.mark.skipif(
    not has_XSPEC, reason="XSPEC not available"
)


@skip_if_xspec_is_not_available
def test_xspec_load():

    # no need to do anything really
    s = XS_phabs() * XS_powerlaw() + XS_bbody()
    print(s(1.0))
    s.set_units(u.keV, 1 / (u.keV * u.cm**2 * u.s))
    print(s(1.0 * u.keV))


@skip_if_xspec_is_not_available
def test_xspec_saving():

    s = XS_powerlaw() + XS_bbody()

    ps = PointSource("test", 0, 0, spectral_shape=s)

    model = Model(ps)

    _ = clone_model(model)

    filename = "_test_xspec_model.yml"

    model.save(filename)

    _ = load_model(filename)

    p = Path(filename)

    p.unlink()
