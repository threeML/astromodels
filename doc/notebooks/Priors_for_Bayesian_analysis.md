---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.8.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Priors for Bayesian analysis

Astromodels supports the definition of priors for all parameters in
your model. You can use as prior any function (although of course not
all functions should be used this way, but the choice is up to you).

First let’s define a simple model containing one point source (see the
“Model tutorial” for more info):

```python
%%capture
from astromodels import *

# Create a point source named "pts1"
pts1 = PointSource('pts1',ra=125.23, dec=17.98, spectral_shape=Powerlaw())

# Create the model
my_model = Model(pts1)
```


Now let’s assign uniform priors to the parameters of the powerlaw
function. The function uniform_prior is defined like this:


```python
Uniform_prior.info()
```

We can use it as such:

```python
# Set 'lower_bound' to 0, 'upper bound' to 10, and leave the 'value' parameter
# to the default value
pts1.spectrum.main.Powerlaw.K.prior = Uniform_prior(lower_bound = 0, upper_bound=10)

# Display it
pts1.spectrum.main.Powerlaw.K.display()

```

Now, lets's set a Gaussian prior on the spectral index

```python

pts1.spectrum.main.Powerlaw.index.prior = Gaussian(mu=-2, sigma=1)

pts1.spectrum.main.Powerlaw.index.display()
```

```python
# Let's get 500 points uniformly distributed between -20 and 20
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)





random_points = np.random.uniform(-10,20,100)

fig, ax = plt.subplots()

ax.plot(random_points,pts1.spectrum.main.Powerlaw.K.prior(random_points), '.' )

ax.set_ylim([-0.1,1.2])
ax.set_xlabel("value of K")
ax.set_ylabel("Prior")
```

```python
random_points = np.random.uniform(-4,0,100)

fig, ax = plt.subplots()

ax.plot(random_points,pts1.spectrum.main.Powerlaw.index.prior(random_points), 'r.' )

ax.set_ylim([-0.1,0.6])
ax.set_xlabel("value of K")
ax.set_ylabel("Prior")
```

```python

```
