.. Astromodels documentation master file, created by
   sphinx-quickstart on Mon Nov 23 22:10:27 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Astromodels's documentation!
=======================================

About astromodels
-----------------

Astromodels is a very flexible framework to define models for likelihood or Bayesian analysis of astrophysical data. 

Even though it has been designed having in mind analysis in the spectral domain, it can be used also as a 
toolbox containing functions of any variable.

Astromodels is *not* a modeling package, it only gives you the tools to build a model as complex as you need. 
You then need a separate package (such as 3ML, github.com/giacomov/3ML) to fit that model to the data.

Some of the features which distinguish astromodels from other similar packages are:
* a model can contain an arbitrary number of sources at different positions in the sky
* parameters can be linked through any function (not only identity)
* parameters can vary with *auxiliary variables* such as time. For example, you can build a model where 
some parameters vary with time, and you can fit the parameters of the function which describe this variability. 
Similary you can build models where parameters vary with the phase of a pulsar, and so on.
* models can be saved in and loaded from YAML file (a human-readable format)
* physical units are fully supported in input, but they are handled so that they don't slow down the actualy 
computation of the models.

Astromodels has been designed with performance as priority, and is considerably faster than other python-based 
solution for the same problem, such as astropy.modeling and the modeling part of sherpa.

Contents:
=========

.. toctree::
   :maxdepth: 3

   Quick_start
   Functions_tutorials
   Point_source_tutorial
   Extended_sources_tutorial
   Model_tutorial
   Priors_for_Bayesian_analysis
   Additional_features_for_scripts_and_applications
   API
