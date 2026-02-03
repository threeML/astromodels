# Installation

If you plan to install `threeML` as well please also refer to its [installation instructions](inv:threeML#installation) which are more detailed.


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
you have no `conda` installation please first of all refer to 
[installing conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
Usually Miniconda or Miniforge are a good start. Follow their installation instructions.

After installing `conda` you are good to go:
1. create an environment and activate it
```bash
conda create -n <name_of_your_environment> python=3.11
conda activate <name_of_your_environment>
```
2. install `astromodels`
```bash
conda install -c threeml -c conda-forge astromodels 
```

or in case you want the latest development version
```bash
conda install -c threeml/label/dev -c conda-forge astromodels
```

Please make sure that the `threeml` channel has a higher priority than the `conda-forge`
one.


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
conda create -n <name_of_your_environment>
conda activate <name_of_your_environment>
```

and then install `XSPEC`. You also might want the `xspec-data` package, please check 
with the official manual.

```bash
conda install -c https://heasarc.gsfc.nasa.gov/FTP/sofware/conda -c conda-forge xspec
conda install -c https://heasarc.gsfc.nasa.gov/FTP/sofware/conda -c conda-forge xspec-data # optional
```

Install `astromodels` using either `conda` or `pip`
```bash
conda install -c threeml -c conda-forge astromodels
```
and finally test it by importing all the models
```bash
python -c "from astromodels.xspec import *"
```

The last step will build all the models and save it to your disk, usually on 
`$HOME/.astromodels/data`. You can also get this path by running 
```bash
python -c "from astromodels.utils.file_utils import get_user_data_path; print(get_user_data_path())"
```
