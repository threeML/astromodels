from pathlib import Path
import os
import collections
import logging

import astropy.units as u
import pytest

from astromodels import Model, PointSource, clone_model, load_model

try:

    from astromodels.xspec import XS_bbody, XS_phabs, XS_powerlaw
    from astromodels.xspec.factory import (
        find_model_dat,
        get_models,
        generate_xs_model_file,
        xspec_model_factory,
    )

except (ImportError, ModuleNotFoundError):

    has_XSPEC = False

else:

    has_XSPEC = True

log = logging.getLogger(__name__)
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


@skip_if_xspec_is_not_available
def test_xspec_model_factory(monkeypatch):
    # first test the cutoff:
    with open("XS_agauss_factory.py", "w+") as f:
        f.write("# dummy_file")

    model_definitions = get_models(find_model_dat())

    monkeypatch.setattr(os.path, "getctime", lambda _: 500)
    monkeypatch.setattr("astromodels.xspec.factory.get_user_data_path", lambda: "")

    xspec_model_factory(
        "agauss_factory",
        "C_agauss",
        "add",
        model_definitions[("agauss", "C_agauss", "add")],
    )
    with open("XS_agauss_factory.py", "r") as f:
        cont = f.read()

    if cont == "# dummy_file":
        raise AssertionError("File was not overridden\n" + cont)
    os.remove("XS_agauss_factory.py")

    monkeypatch.undo()
    # now lets test the non-overriding version
    with open("XS_agauss_factory.py", "w+") as f:
        f.write("# dummy_file")

    # still need that
    monkeypatch.setattr("astromodels.xspec.factory.get_user_data_path", lambda: "")

    xspec_model_factory(
        "agauss_factory",
        "C_agauss",
        "add",
        model_definitions[("agauss", "C_agauss", "add")],
    )
    with open("XS_agauss_factory.py", "r") as f:
        cont = f.read()

    if cont != "# dummy_file":
        raise AssertionError("File was overridden\n" + cont)

    os.remove("XS_agauss_factory.py")

    # and the create from scratch
    xspec_model_factory(
        "agauss_factory",
        "C_agauss",
        "add",
        model_definitions[("agauss", "C_agauss", "add")],
    )
    os.remove("XS_agauss_factory.py")


@skip_if_xspec_is_not_available
def test_generate_xs_model_file():
    with pytest.raises(AssertionError):
        generate_xs_model_file("empty", "test_model", "test_model", "con", {})
    # get the actual model definitions
    model_definitions = get_models(find_model_dat())

    # create a temporary agauss file
    generate_xs_model_file(
        "test_agauss.py",
        "agauss",
        "C_agauss",
        "add",
        model_definitions[("agauss", "C_agauss", "add")],
    )
    os.remove("test_agauss.py")
