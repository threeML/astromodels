from astromodels.functions.priors import Uniform_prior, Log_uniform_prior, Gaussian
import numpy as np


def test_bounded_prior():
    logp = Log_uniform_prior(lower_bound=0.1, upper_bound=1)
    assert logp(0.2) == logp.scipy_dist.pdf(0.2)

    bounds = logp.scipy_dist.support()
    assert bounds[0] == logp.lower_bound.value
    assert bounds[1] == logp.upper_bound.value

    uni = Uniform_prior(lower_bound=-1, upper_bound=1)
    rvs = uni.scipy_dist.rvs(size=100000)
    assert np.isclose(rvs.mean(), 0.0, atol=rvs.std() * 2)

    bounds = uni.scipy_dist.support()
    assert bounds[0] == uni.lower_bound.value
    assert bounds[1] == uni.upper_bound.value
