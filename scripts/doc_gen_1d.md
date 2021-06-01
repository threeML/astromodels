---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Python 3
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

x_scale="log"
y_scale="log"

linear_range = False

wide_energy_range = False
```

```python nbsphinx="hidden" tags=[]
func = _known_functions[func_name]()

if wide_energy_range:

    energy_grid = np.geomspace(1e2,1e4,500)
    
else:
    
    energy_grid = np.geomspace(2e-1,1e1,1000)

if linear_range:

	energy_grid = np.linspace(-5,5,1000)

    
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


ax.plot(energy_grid, func(energy_grid), color=blue)

ax.set_xlabel("energy (keV)")
ax.set_ylabel("photon flux")
ax.set_xscale(x_scale)
ax.set_yscale(y_scale)

```

## F$_{\nu}$

The F$_{\nu}$ shape of the photon model
*if this is not a photon model, please ignore this auto-generated plot*
```python
fig, ax = plt.subplots()

ax.plot(energy_grid, energy_grid * func(energy_grid), red)


ax.set_xlabel("energy (keV)")
ax.set_ylabel(r"energy flux (F$_{\nu}$)")
ax.set_xscale(x_scale)
ax.set_yscale(y_scale)


```

## $\nu$F$_{\nu}$

The $\nu$F$_{\nu}$ shape of the photon model
*if this is not a photon model, please ignore this auto-generated plot*

```python
fig, ax = plt.subplots()

ax.plot(energy_grid, energy_grid**2 * func(energy_grid), color=green)


ax.set_xlabel("energy (keV)")
ax.set_ylabel(r"$\nu$F$_{\nu}$")
ax.set_xscale(x_scale)
ax.set_yscale(y_scale)

```
