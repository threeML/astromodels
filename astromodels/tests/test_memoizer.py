from astromodels.functions.functions import Powerlaw
import numpy as np


def test_memoizer():

    po = Powerlaw()

    a = np.random.uniform(-3,-1, 2000)
    b = np.random.uniform(0.1,10, 2000)

    for aa, bb in zip(a, b):

        po.index = aa
        po.K = bb

        po(1.0)


