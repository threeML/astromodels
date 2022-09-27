# astromodels

![CI](https://github.com/threeML/astromodels/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/threeML/astromodels/branch/master/graph/badge.svg)](https://codecov.io/gh/threeML/astromodels)
[![Documentation Status](https://readthedocs.org/projects/astromodels/badge/?version=latest)](http://astromodels.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub contributors](https://img.shields.io/github/contributors/threeML/astromodels)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5646925.svg)](https://doi.org/10.5281/zenodo.5646925)

![GitHub pull requests](https://img.shields.io/github/issues-pr/threeML/astromodels)
![GitHub issues](https://img.shields.io/github/issues/threeML/astromodels)
## PyPi

[![PyPI version fury.io](https://badge.fury.io/py/astromodels.svg)](https://pypi.python.org/pypi/astromodels/)
![PyPI - Downloads](https://img.shields.io/pypi/dm/astromodels)

## Conda
![Conda](https://img.shields.io/conda/pn/threeml/astromodels)
![Conda](https://img.shields.io/conda/dn/threeml/astromodels)


![alt text](https://raw.githubusercontent.com/threeml/astromodels/master/docs/media/large_logo.png)

Astromodels is a very flexible framework to define models for likelihood or Bayesian analysis of astrophysical data.

Even though it has been designed having in mind analysis in the spectral domain, it can be used also as a toolbox containing functions of any variable.

Astromodels is not a modeling package, it only gives you the tools to build a model as complex as you need. You then need a separate package (such as [3ML](github.com/threeML/threeML)) to fit that model to the data.

Some of the features which distinguish astromodels from other similar packages are: * a model can contain an arbitrary number of sources at different positions in the sky * parameters can be linked through any function (not only identity) * parameters can vary with auxiliary variables such as time. For example, you can build a model where some parameters vary with time, and you can fit the parameters of the function which describe this variability. Similary you can build models where parameters vary with the phase of a pulsar, and so on. * models can be saved in and loaded from YAML file (a human-readable format) * physical units are fully supported in input, but they are handled so that they don’t slow down the actualy computation of the models.

Astromodels has been designed with performance as priority, and is considerably faster than other python-based solution for the same problem, such as astropy.modeling and the modeling part of sherpa.
Documentation: http://astromodels.readthedocs.org/en/latest/
