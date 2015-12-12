#!/usr/bin/env python

from setuptools import setup
import glob

setup(
    name="astromodels",

    packages=['astromodels',
              'astromodels/functions',
              'astromodels/tests'
              ],

    data_files=[('astromodels/data/functions', glob.glob('astromodels/data/functions/*.yaml'))],

    version='v0.1',

    description="Astromodels contains models to be used in likelihood or Bayesian analysis in astronomy",

    author='Giacomo Vianello',

    author_email='giacomo.vianello@gmail.com',

    url='https://github.com/giacomov/astromodels',

    download_url='https://github.com/giacomov/astromodels/archive/v0.1',

    keywords=['Likelihood', 'Models', 'fit'],

    classifiers=[],

    install_requires=[
        'numpy >= 1.6',
        'pyyaml']

)
