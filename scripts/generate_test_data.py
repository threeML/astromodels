# this generates test data to make sure changes in the
# code do not destroy the values. If a new function is added
# this needs to be re run
#
#

from astromodels.functions.priors import *
from astromodels.functions.function import _known_functions
import h5py

eval_x = np.logspace(-1,3, 10)
_multiplicative_models = ["PhAbs", "TbAbs", "WAbs", "APEC", "VAPEC"]


with h5py.File("past_1D_values.h5", "w") as f:
    
    f.create_dataset("eval_values", data=eval_x, compression="lzf")
    
    for key in _known_functions:

        this_function = _known_functions[key]

        # Test only the power law of XSpec, which is the only one we know we can test at 1 keV

        if key.find("XS")==0 and key != "XS_powerlaw" or (key in _multiplicative_models):

            # An XSpec model. Test it only if it's a power law (the others might need other parameters during
            # initialization)

            continue

        if key.find("TemplateModel")==0:

            # The TemplateModel function has its own test

            continue

    #        if key.find("Synchrotron")==0:

            # Naima Synchtron function should have its own test

    #            continue

        if this_function._n_dim == 1:

            print("testing %s ..." % key)

            func = this_function()

            data=func(eval_x)
            
            print(data)
            
            f.create_dataset(key, data=np.atleast_1d(data), compression="lzf")


