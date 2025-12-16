# Installation

If you plan to install `threeML` as well please also refer to its [install instructions](inv:threeML#installation) which are more detailed.


## TL;DR
**careful with dependencies!**

with `conda` without `XSPEC`
```bash
conda install -c threeml -c conda-forge astromodels 
```

with `conda` with `XSPEC`
```bash
conda install -c https://heasarc.gsfc.nasa.gov/FTP/software/conda/ -c conda-forge xspec python=3.11
conda install -c threeml -c conda-forge astromodels 
```

with `pip` (please take care of dependencies beforehand!)
```bash
pip install astromodels 
```

in case you want the development versions:
```bash
conda install -c threeml/label/dev -c conda-forge astromodels
```
or
```bash
pip install --upgrade --pre astromodels 
```


## Without XSPEC

Installing `astromodels` without `XSPEC` is straigt forward.

### conda
We here assume that you have a basic knowledge of `conda`. If this is not the case and 
you have no `conda` installation

Simply create 

### pip
Directly install from PyPI:
```bash
pip install astromodels
```

and in case you want the latest development version use

```bash
pip install --upgrade --pre astromodels
```


## With XSPEC
The easiest way to install `astromodels` with `XSPEC` support is via `conda`.
Note that you first need to install `XSPEC` before installing `astromodels`.


First let's create a new environment
```bash
conda create -n threeml
conda activate threeml
conda install -c https://heasarc.gsfc.nasa.gov/FTP/sofware/conda -c conda-forge xspec
conda install -c https://heasarc.gsfc.nasa.gov/FTP/sofware/conda -c conda-forge xspec-data
conda install -c threeml astromodels
```
