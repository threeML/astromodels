from pathlib import Path
import os
import collections

import astropy.units as u
import pytest

from astromodels import Model, PointSource, clone_model, load_model

try:

    from astromodels.xspec import XS_bbody, XS_phabs, XS_powerlaw
    from astromodels.xspec.factory import (
        find_model_dat,
        get_models,
    )

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


@skip_if_xspec_is_not_available
def test_find_model_dat():
    model_dat_path = find_model_dat()
    assert Path(model_dat_path).is_file(), "model.dat file path is not a file"

    headas_env = os.environ.get("HEADAS")
    os.environ.update({"HEADAS": "/"})  # set HEADAS to root dir
    with pytest.raises(FileNotFoundError):
        find_model_dat()
    os.environ.update({"HEADAS": headas_env})


@skip_if_xspec_is_not_available
def test_get_models():
    model_definitions = get_models(find_model_dat())
    assert isinstance(model_definitions, collections.OrderedDict)

    # the gaussian line should be available in all versions supported
    assert (
        "agauss",
        "C_agauss",
        "add",
    ) in model_definitions.keys(), "agauss not in model.dat"
