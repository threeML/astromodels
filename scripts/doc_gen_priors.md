---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# func_title

```python nbsphinx="hidden" tags=[]
%%capture

import numpy as np

import matplotlib.pyplot as plt

import warnings
warnings.simplefilter("ignore")

from astromodels.functions.function import _known_functions


from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)
%matplotlib inline
```

```python nbsphinx="hidden" tags=["parameters"]
func_name = "TbAbs"

positive_prior = False

```

```python nbsphinx="hidden" tags=[]
func = _known_functions[func_name]()

if not positive_prior:

	energy_grid = np.linspace(-5,5,1000)

else:
    
    energy_grid = np.linspace(0,1,1000)
    
    
    
    
blue = "#4152E3"
red = "#E3414B"
green = "#41E39E"
```
## Description
```python
func.display()
```

## Shape 

The shape of the function. 

*If this is not a photon model but a prior or linear function then ignore the units as these docs are auto-generated*

```python tags=["nbsphinx-thumbnail"]
fig, ax = plt.subplots()


ax.plot(energy_grid, func(energy_grid), color=blue, lw=3)

ax.set_xlabel("x")
ax.set_ylabel("probability")

```

## Random Number Generation

This is how we can generate random numbers from the prior.


```python


u = np.random.uniform(0,1, size=5000)

draws = [func.from_unit_cube(x) for x in u]


fig, ax = plt.subplots()


ax.hist(draws, color=green, bins=50)

ax.set_xlabel("value")
ax.set_ylabel("N")


```
