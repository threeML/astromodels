import pytest
import os
import numpy as np

from astromodels.functions.template_model import TemplateModel, TemplateModelFactory, MissingDataFile
from astromodels.functions.functions import Band, Powerlaw
from astromodels import Model, PointSource, clone_model, load_model
import pickle

__author__ = 'giacomov'


def get_comparison_function():

    mo = Band()
    mo.K = 1

    return mo

@pytest.mark.slow
def test_template_factory():

    mo = get_comparison_function()

    energies = np.logspace(1, 3, 50)

    t = TemplateModelFactory('__test', 'A test template', energies, ['alpha', 'xp', 'beta'])

    alpha_grid = np.linspace(-1.5, 1, 15)
    beta_grid = np.linspace(-3.5, -1.6, 15)
    xp_grid = np.logspace(1, 3, 20)

    t.define_parameter_grid('alpha', alpha_grid)
    t.define_parameter_grid('beta', beta_grid)
    t.define_parameter_grid('xp', xp_grid)

    for a in alpha_grid:

        for b in beta_grid:

            for xp in xp_grid:
                mo.alpha = a
                mo.beta = b
                mo.xp = xp

                t.add_interpolation_data(mo(energies), alpha=a, xp=xp, beta=b)

    print("Data has been prepared")

    t.save_data(overwrite=True)


# This will be run second, so the template will exist
@pytest.mark.slow
def test_template_function():

    tm = TemplateModel('__test')

    mo = get_comparison_function()

    new_alpha_grid = np.linspace(-1.5, 1, 20)
    new_beta_grid = np.linspace(-3.5, -1.6, 20)
    new_xp_grid = np.logspace(1, 3, 30)

    new_energies = np.logspace(1, 3, 40)

    tm.K = 1

    mo.K = 1

    for a in new_alpha_grid:

        for b in new_beta_grid:

            for xp in new_xp_grid:

                mo.alpha = a
                mo.beta = b
                mo.xp = xp

                tm.alpha = a
                tm.beta = b
                tm.xp = xp

                res1 = mo(new_energies)
                res2 = tm(new_energies)

                deltas = np.abs((res2 - res1) / res1)

                idx = deltas > 0.1

                if np.any(idx):

                    raise AssertionError("Interpolation precision @ %s is %s, "
                                         "worse than 10 percent, "
                                         "with parameters %s!" % (new_energies[idx], deltas[idx], [a,b,xp]))


def test_input_output():

    tm = TemplateModel('__test')
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

    assert np.allclose(clone.test.spectrum.main.shape(xx), fake_model.test.spectrum.main.shape(xx))

    # Test pickling
    dump = pickle.dumps(clone)

    clone2 = pickle.loads(dump)

    assert clone2.get_number_of_point_sources() == 1
    assert tm.data_file == clone2.test.spectrum.main.shape.data_file
    assert np.allclose(clone2.test.spectrum.main.shape(xx), fake_model.test.spectrum.main.shape(xx))

    # Test pickling with other functions
    new_shape = tm * Powerlaw()

    new_shape.index_2 = -2.256

    dump2 = pickle.dumps(new_shape)

    clone3 = pickle.loads(dump2)

    assert clone3.index_2.value == new_shape.index_2.value

    # Now save to disk and reload
    fake_source2 = PointSource("test", ra=0.0, dec=0.0, spectral_shape=new_shape)

    fake_model2 = Model(fake_source2)

    fake_model2.save("__test.yml", overwrite=True)

    # Now try to reload
    reloaded_model = load_model("__test.yml")

    assert reloaded_model.get_number_of_point_sources() == 1
    assert np.allclose(fake_model2.test.spectrum.main.shape(xx), reloaded_model.test.spectrum.main.shape(xx))

    os.remove("__test.yml")

