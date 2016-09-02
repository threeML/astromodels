import pytest

import numpy as np

from astromodels.functions.template_model import TemplateModel, TemplateModelFactory, MissingDataFile
from astromodels.functions.functions import Band

__author__ = 'giacomov'


def get_comparison_function():

    mo = Band()
    mo.K = 1

    return mo


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


