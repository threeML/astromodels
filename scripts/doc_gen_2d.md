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

import healpy as hp


from astropy.coordinates import ICRS, Galactic


from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=False, grid=False)
%matplotlib inline




```

```python nbsphinx="hidden" tags=["parameters"]
func_name = "Latitude_galactic_diffuse"

```

```python nbsphinx="hidden" tags=[]
func = _known_functions[func_name]()

NSIDE = 2**7
NPIX = hp.nside2npix(NSIDE)

ra, dec  = hp.pix2ang(nside=NSIDE, ipix=np.arange(NPIX), lonlat=True)

```
## Description
```python
func.display()
```

## Shape

The shape of the function on the sky.
```python tags=["nbsphinx-thumbnail"]


m=func(ra, dec)
hp.mollview(m, title=func_name, cmap="magma")
hp.graticule(color="grey", lw=2)



```
