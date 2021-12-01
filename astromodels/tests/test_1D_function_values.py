
import h5py
import numpy as np
import numpy.testing as npt

from astromodels.functions.function import _known_functions
from astromodels.functions.priors import *
from astromodels.utils import _get_data_file_path

_multiplicative_models = ["PhAbs", "TbAbs", "WAbs", "APEC", "VAPEC", "EBLattenuation" ]

def test_function_values_have_not_changed():

    with h5py.File(_get_data_file_path("past_1D_values.h5"), "r") as f:

        eval_x = f["eval_values"][()]

    
    for key in _known_functions:

        this_function = _known_functions[key]

        # Test only the power law of XSpec, which is the only one we know we can test at 1 keV

        if key.find("XS")==0 or (key in _multiplicative_models):

            # An XSpec model OR EBLattenuation function. Test it only if it's a power law (the others might need other parameters during
            # initialization)

            continue

        if key.find("TemplateModel")==0:

            # The TemplateModel function has its own test

            continue

        if key.find("Synchrotron")==0:

        #    Naima Synchtron function should have its own test

            continue

        if key.find("_ComplexTestFunction")==0:

        #    Naima Synchtron function should have its own test

            continue

        
        if this_function._n_dim == 1:

            print("testing %s ..." % key)


            
            
            func = this_function()
            
            new_values = np.atleast_1d(func(eval_x))

            with h5py.File(_get_data_file_path("past_1D_values.h5"), "r") as f:
                if key not in f.keys():

                    
                    print("the function %s does not exist in the past data. You must run a script to add it" %key)

                else:
                
                    old_values = f[key][()]

            
            
                    npt.assert_almost_equal(new_values, old_values)
            
