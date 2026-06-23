import h5py
import numpy as np
import numpy.testing as npt

from astromodels.utils import _get_data_file_path
from astromodels.utils.list_functions import (
    get_function_class,
    list_function_names,
)

_multiplicative_models = [
    "PhAbs",
    "TbAbs",
    "WAbs",
    "APEC",
    "VAPEC",
    "EBLattenuation",
]


def test_function_values_have_not_changed():

    with h5py.File(_get_data_file_path("tests/past_1D_values.h5"), "r") as f:

        eval_x = f["eval_values"][()]

    for key in list_function_names():

        this_function = get_function_class(key)

        # Test only the power law of XSpec, which is the only one we know we can test
        # at 1 keV

        if key.find("XS") == 0 or (key in _multiplicative_models):

            # An XSpec model OR EBLattenuation function. Test it only if it's a power
            # law (the others might need other parameters during initialization)

            continue

        if key.find("TemplateModel") == 0:

            # The TemplateModel function has its own test

            continue

        if key.find("Synchrotron") == 0:

            #    Naima Synchtron function should have its own test

            continue

        if key.find("_ComplexTestFunction") == 0:

            #    Naima Synchtron function should have its own test

            continue

        if this_function._n_dim == 1:

            print("testing %s ..." % key)

            func = this_function()
            if key == "GenericFunction":

                def _f(x):
                    return x**2

                func.set_function(_f)

            new_values = np.atleast_1d(func(eval_x))

            with h5py.File(_get_data_file_path("tests/past_1D_values.h5"), "r") as f:
                if key not in f.keys():
                    msg = f"the function {key} does not exist in the past data. You"
                    msg += " must run a script to add it"
                    print(msg)

                else:

                    old_values = f[key][()]

                    npt.assert_almost_equal(new_values, old_values)

            if key == "Cutoff_powerlaw_Ep":
                func = this_function()
                func.index.value = -3
                func.xp.value = 10
                assert np.isclose(func(1), 1.10517)


def test_priors():
    np.random.seed = 123
    for key in list_function_names():
        this_function = get_function_class(key)

        if (
            hasattr(this_function, "from_unit_cube")
            and this_function().name
            != "Cauchy"  # Cauchy is too heavy tailed apparaently
        ):
            func = this_function()
            low = func.from_unit_cube(1e-9)
            high = func.from_unit_cube(1 - 1e-9)
            assert low > -np.inf
            assert high < np.inf
