import pickle
import shutil

import numpy as np
import numpy.testing as npt
import pytest

from astromodels import Model, PointSource, clone_model, load_model
from astromodels.functions import (
    Band,
    Powerlaw,
    TemplateModel,
    TemplateModelFactory,
    XSPECTableModel,
)
from astromodels.functions.template_model import (
    convert_old_table_model,
    UnivariateSpline,
)
from astromodels.utils import _get_data_file_path
from astromodels.utils.file_utils import get_user_data_path

__author__ = "giacomov"


def get_comparison_function():

    mo = Band()
    mo.K = 1

    return mo


@pytest.mark.slow
def test_template_factory_1D():

    mo = get_comparison_function()

    energies = np.logspace(1, 3, 50)

    t = TemplateModelFactory("__test1D", "A test template", energies, ["alpha"])

    alpha_grid = np.linspace(-1.5, 1, 15)
    # beta_grid = np.linspace(-3.5, -1.6, 15)
    # xp_grid = np.logspace(1, 3, 20)

    t.define_parameter_grid("alpha", alpha_grid)

    for a in alpha_grid:
        mo.alpha = a
        mo.beta = -2.5
        mo.xp = 300.0

        t.add_interpolation_data(mo(energies), alpha=a)

    print("Data has been prepared")

    t.save_data(overwrite=True)


@pytest.mark.slow
def test_template_factory(tmp_path):

    mo = get_comparison_function()

    energies = np.logspace(1, 3, 50)

    t = TemplateModelFactory(
        "__test", "A test template", energies, ["alpha", "xp", "beta"]
    )

    alpha_grid = np.linspace(-1.5, 1, 15)
    beta_grid = np.linspace(-3.5, -1.6, 15)
    xp_grid = np.logspace(1, 3, 20)

    t.define_parameter_grid("alpha", alpha_grid)
    t.define_parameter_grid("beta", beta_grid)
    t.define_parameter_grid("xp", xp_grid)

    for a in alpha_grid:

        for b in beta_grid:

            for xp in xp_grid:
                mo.alpha = a
                mo.beta = b
                mo.xp = xp

                t.add_interpolation_data(mo(energies), alpha=a, xp=xp, beta=b)

    t.save_data(overwrite=True)

    tm = TemplateModel("__test")

    tm(energies)

    tm.clean()

    tm1 = TemplateModel("__test")

    mo1 = get_comparison_function()

    new_alpha_grid = np.linspace(-1.5, 1, 15)
    new_beta_grid = np.linspace(-3.5, -1.6, 15)
    new_xp_grid = np.logspace(1, 3, 15)

    new_energies = np.logspace(1, 3, 40)

    tm1.K = 1

    mo1.K = 1

    for a in new_alpha_grid:

        for b in new_beta_grid:

            for xp in new_xp_grid:

                mo1.alpha = a
                mo1.beta = b
                mo1.xp = xp

                tm1.alpha = a
                tm1.beta = b
                tm1.xp = xp

                res1 = mo1(new_energies)
                res2 = tm1(new_energies)

                deltas = np.abs((res2 - res1) / res1)

                idx = deltas > 0.1

                if np.any(idx):

                    raise AssertionError(
                        "Interpolation precision @ %s is %s, "
                        "worse than 10 percent, "
                        "with parameters %s!"
                        % (new_energies[idx], deltas[idx], [a, b, xp])
                    )

    tm = TemplateModel("__test")
    tm.alpha = -0.95
    tm.beta = -2.23

    fake_source = PointSource("test", ra=0.0, dec=0.0, spectral_shape=tm)

    fake_model = Model(fake_source)

    clone = clone_model(fake_model)

    assert clone.get_number_of_point_sources() == 1
    assert tm.data_file == clone.test.spectrum.main.shape.data_file

    assert clone.test.spectrum.main.shape.alpha.value == tm.alpha.value
    assert clone.test.spectrum.main.shape.beta.value == tm.beta.value

    xx = np.linspace(1, 10, 100)

    assert np.allclose(
        clone.test.spectrum.main.shape(xx),
        fake_model.test.spectrum.main.shape(xx),
    )

    # Test pickling
    dump = pickle.dumps(clone)

    clone2 = pickle.loads(dump)

    assert clone2.get_number_of_point_sources() == 1
    assert tm.data_file == clone2.test.spectrum.main.shape.data_file
    assert np.allclose(
        clone2.test.spectrum.main.shape(xx),
        fake_model.test.spectrum.main.shape(xx),
    )

    new_shape = tm + Powerlaw()

    new_shape.index_2 = -2.256

    dump2 = pickle.dumps(new_shape)

    clone3 = pickle.loads(dump2)

    assert clone3.index_2.value == new_shape.index_2.value

    # Now save to disk and reload
    fake_source2 = PointSource("test", ra=0.0, dec=0.0, spectral_shape=new_shape)

    fake_model2 = Model(fake_source2)

    fake_model2.save(tmp_path / "__test.yml", overwrite=True)

    # Now try to reload
    reloaded_model = load_model(tmp_path / "__test.yml")

    assert reloaded_model.get_number_of_point_sources() == 1
    assert np.allclose(
        fake_model2.test.spectrum.main.shape(xx),
        reloaded_model.test.spectrum.main.shape(xx),
    )

    new_shape2 = Powerlaw() + tm

    new_shape2.index_1 = -2.256

    dump2 = pickle.dumps(new_shape2)

    clone3 = pickle.loads(dump2)

    assert clone3.index_1.value == new_shape2.index_1.value

    # Now save to disk and reload
    fake_source2 = PointSource("test", ra=0.0, dec=0.0, spectral_shape=new_shape2)

    fake_model2 = Model(fake_source2)

    fake_model2.save(tmp_path / "__test.yml", overwrite=True)

    # Now try to reload
    reloaded_model = load_model(tmp_path / "__test.yml")

    assert reloaded_model.get_number_of_point_sources() == 1
    assert np.allclose(
        fake_model2.test.spectrum.main.shape(xx),
        reloaded_model.test.spectrum.main.shape(xx),
    )

    old_table_file = _get_data_file_path("tests/old_table.h5")

    p = get_user_data_path() / "old_table.h5"

    shutil.copy(old_table_file, p)

    # convert the table

    convert_old_table_model("old_table")

    # now load the old table

    old_table = TemplateModel("old_table")

    # should be the same as in test

    test = TemplateModel("__test")

    xx = np.logspace(1, 3, 50)

    npt.assert_almost_equal(test(xx), old_table(xx))

    p.unlink()


def test_xspec_table_model():

    test_table = _get_data_file_path("tests/test_xspec_table_model.fits")

    xtm = XSPECTableModel(test_table)

    xtm.to_table_model("xspectm_test", "xspec model", overwrite=True)


def test_univariate_spline():
    univ = UnivariateSpline([0, 1, 2], [3, 4, 5])
    res = univ(1.5)
    assert np.isclose(res, 4.5)

    univ_np = UnivariateSpline(np.array([0, 1, 2]), np.array([3, 4, 5]))
    res_np = univ_np(1.5)
    assert np.isclose(res_np, 4.5)
