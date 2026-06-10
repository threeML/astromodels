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

# Units in `astromodels`

In `astromodels` there are conceptually two different kinds of units:

- external units: provided to and by you, the user
- internal units: used for calculating spectra, functions, etc.

As the intention behind `astromodels` is to be used with some separate package (e.g. [https://threeml.readthedocs.io](threeML))
the functions are evaluated only with bare values to remove any overhead from using `astropy.units`. This significantly reduces the time spent in `astromodels` functionalities during a fit.
The internal units are therefore set and fixed beforehand in the `astromodels.configuration` and not explicitly passed during the evaluation.

For this reason a couple of gotchas right before diving deeper:

## Gotchas

When using `astromodels` for simply plotting a `Powerlaw` or just wanting to use it as a calculator, use the `__call__` function of the Function, not the `evaluate` one.

For the `__call__` function there are two versions implemented in the backend:
- one without units (the fast one)
- one with units


Let's look at an example. First we create a simple Powerlaw()
```python
from astromodels.functions import Powerlaw

pl = Powerlaw()
pl.info()
```
