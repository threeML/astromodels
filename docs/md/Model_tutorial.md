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

# Model tutorial

In this tutorial we show how to build a simple model with two point sources, how to save it for later use, and re-load it back. We will also plot the spectra of the two point sources, with their components.

See the documents about creating and getting information about functions, point sources and extended sources for details about these operations. Here we only focus on the global model.

```python
%%capture
from astromodels import *

# We also import astropy units to show the unit-conversion
# feature
import astropy.units as u
```

## Define sources

Now let’s define a point source (see “Creating point sources” for details and alternative way to accomplish this):

```python
# Let's start by defining a simple point source with a power law spectrum
powerlaw = Powerlaw()

pts1 = PointSource('source_1', ra=125.6, dec=-75.3,
                   spectral_shape=powerlaw)

# Get some info about what we just created
pts1.display()

# Have a quicker look
pts1.plot_tree()

```

Now let’s define another source, this time at Galactic Coordinates l = 11.25, b = -22.5, and with two spectral components:

```python
# Another point source with two spectral components

spectrum1 = Powerlaw(K=0.2, index=-0.75)
component1 = SpectralComponent('synchrotron',spectrum1)

spectrum2 = Powerlaw(K=0.8, index=-1.0)
component2 = SpectralComponent('IC',spectrum2)

point_source2 = PointSource('source_2', l=11.25, b=-22.5, components=[component1,component2])

# Have a look at what we just created

point_source2.display()
```

## Create a model

Now let’s create our model, which comprises our two sources:

```python
# Build a model with the two point sources

my_model = Model(pts1, point_source2)
```

Of course you can use as many sources as needed, like my_model = Model(pts1, pts2, pts3…)

## Getting information about a model

Using the ```.display()``` method we can see all free parameters currently in the model:

```python
my_model.display()
```

The model tree can be shown as:

```python
my_model.plot_tree()
```

A dictionary of free parameters can be obtained like this:

```python
free_parameters = my_model.free_parameters
```

We can use such dictionary to loop over all free parameters:

```python
for parameter_name, parameter in free_parameters.items():

    print("Parameter %s is free" % parameter_name)
```

More information on a particular source can be obtained like:

```python
my_model.source_1.display()
```

More information about a particular instance of a function can be obtained like:

```python
my_model.source_1.spectrum.main.Powerlaw.display()
```

## Accessing and modifying sources and parameters from the model instance

### Fully-qualified paths

Each source and each parameter has a precise path within the model. These paths are displayed by the .display() method of the model instance (see above), and can be used like ```my_model.[path]```. For example:

```python
my_model.display()
```

```python
# Access the logK parameters of the powerlaw spectrum of the main component for source 1:

my_model.source_1.spectrum.main.Powerlaw.logK = -0.5

# Access the logK parameters of the spectrum of the IC component of source 2:

my_model.source_2.spectrum.IC.Powerlaw.logK = -0.32
```

The structure of these paths is easy to understand. The model is a tree-like structure. The root of the tree is always the model instance itself. The second level is constituted by the various sources. The structure within a source can be understood by calling the ```.display``` method:

```python
my_model.source_1.display()
```

Each indentation represents one level, so to access the “ra” element we can follow the levels shown by the ```.display()``` method:

```python
ra_parameter = my_model.source_1.position.ra

ra_parameter.display()
```

<div class="alert alert-info">

**Note:** this is a Parameter instance. To get the position of the source as a floating point number, use:
    ```my_model.source_1.position.get_ra()```
    which will work for any source

</div>


while to access the index parameter of the power law function we can do

```python
K_parameter = my_model.source_1.spectrum.main.Powerlaw.K

K_parameter.display()
```

You can find much more information in the document “Additional features for scripts and applications”.

These fully-qualified paths are unique to each element, are very descriptive and easy to understand. They can always be used and are encouraged in general, but especially in scripts, when the effort spent writing them is paid off in terms of clarity. However, there is an alternative way which might be more convenient in certain situation, especially when models are simple and the chances of getting confused are low. This alternative method is described below.

### Using shortcuts

Exploiting the feature of the python language, we can create names (“shortcuts”) for objects:

```python
# Create a "shortcut" for the spectrum of a source

powerlaw_1 = my_model.source_1.spectrum.main.Powerlaw

# Now we can change the values of that power law as:
powerlaw_1.K = 1.2

# GOTCHA: while it is possible to create shortcuts for parameters, it is not encouraged
# Indeed, this will not work:
# logK_1 = my_model.source_1.spectrum.main.powerlaw.logK
# logK_1 = -1.2 # WILL NOT WORK
# In order to use a shortcut for a parameter to change its value, you have to explicitly
# set its property 'value':
# logK_1.value = -1.2 # This will work
```

<div class="alert alert-info">

**GOTACH:** while it is possible to create shortcuts for parameters, it is not encouraged. Indeed, this will not work:
   
</div>

K_1 = my_model.source_1.spectrum.main.powerlaw.K

K_1 = 1.2  *WILL NOT WORK*
    
In order to use a shortcut for a parameter to change its value, you have to explicitly
set its property 'value':
K_1.value = 1.2  *This will work*



Shortcut can point at any point of the tree:


```python
# Create a shortcut of a source
source_1 = my_model.source_1

# Now we can do:
source_1.spectrum.main.Powerlaw.index = -2.3

# Create a shortcut for a component

main_component = my_model.source_1.spectrum.main

# Now we can do:
main_component.Powerlaw.index = -1.3
```

If you are ever in doubt of what a particular shortcut stands for, you can always retrieve the full path of the element the shortcut is pointing to like this:

```python
print(main_component.path)
```

## Saving a model to file

An existing model can be saved to a file with:

```python
# Save the model to a file, overwriting it if already existing

my_model.save('my_model.yml', overwrite=True)
```

The content of the file is YAML code, which is human-readable and very easy to understand. Let’s have a look:

```python
with open('my_model.yml') as yaml_file:

    print("".join(yaml_file.readlines()))
```

## Load a model from a file

Now suppose that you want to load back a file you created in a previous session. You can do it with:

```python
my_model = load_model('my_model.yml')

# Explore the model we just loaded back

my_model.display()
```

```python
# Now evaluate and plot our models. You need matplotlib for this

import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)



fig, ax = plt.subplots()
# Energies where we want to evaluate the model

e = np.geomspace(1,1000,100)

# Loop over the sources

for src_name, src in my_model.point_sources.items():

    # Loop over the components of each source

    for comp_name, component in src.components.items():

        # Get the differential flux (in ph/cm2/s)

        flux = component.shape(e)

        # this can also be accomplished with:
        # flux = component.powerlaw(e)
        # but this requires to know the name of the
        # spectral shape which was used

        # Plot this component for this source

        ax.loglog(e, flux,label="%s of %s" % (component.name, src.name))

ax.legend(loc=0)

ax.set_xlabel("Energy")
ax.set_ylabel(r"Flux (ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$")
```

## Linking parameters

Sometimes you want to link two parameters of a model so that they have the same value. This can be easily accomplished in astromodels:

```python
# Link the photon index of the first source with the
# photon index of the IC component of the second source

my_model.link(my_model.source_2.spectrum.IC.Powerlaw.index,
              my_model.source_1.spectrum.main.Powerlaw.index)

my_model.display()
```

### Advanced use of linking: arbitrary functions

Astromodels takes this a step further. Parameters can be linked to each other through any function. The parameters of the linking function become parameters of the model like any other, and can be left free to vary or fixed. For example, let’s consider the case where we want the photon index of the IC component of the second source (p2) to be equal to the photon index of the first source (p1) plus a constant. We can link the two parameters with the ‘bias’ function $f(x) = x + k$, so that $p2(p1) = p1 + k$:

```python
# Link the photon indexes through the 'linear' function, i.e.,
# the photon index of the IC component of the second source is fixed to be the
# photon index of the first source plus a constant k

link_function = Line(a=2, b=1)
link_function.b.fix=True

my_model.link(my_model.source_2.spectrum.IC.Powerlaw.index,
              my_model.source_1.spectrum.main.Powerlaw.index,
              link_function)

# The parameters of the linking function become parameters
# of the model, and are put in the model tree under the parameter they are
# linking.
# In this case the only parameter of the 'linear' function ('k') becomes then
# my_model.source_2.spectrum.IC.powerlaw.logK.bias.k

my_model.display()
```

If we want to fix say $p2 = p1 - 1.2$, we can fix k to that:

```python
my_model.source_2.spectrum.IC.Powerlaw.index.Line.a = -1.2
my_model.source_2.spectrum.IC.Powerlaw.index.Line.a.fix = True

my_model.display()
```

As another example, we might link the two parameters using a power law function:

```python
my_model.link(my_model.source_2.spectrum.IC.Powerlaw.index,
              my_model.source_1.spectrum.main.Powerlaw.index,
              Powerlaw())

my_model.display()
```

We can use arbitrarily complex functions as link function, if needed (see “Creating and modifying functions” for more info on how to create composite functions):

```python
# A random composite function (see "Creating and modifying functions" for more info)
crazy_link = Powerlaw() + Gaussian()

my_model.link(my_model.source_2.spectrum.IC.Powerlaw.index,
              my_model.source_1.spectrum.main.Powerlaw.index,
              crazy_link)

my_model.display()
```

## Time-varying models and other independent variables

In astromodels parameters can become functions of independent variables such as time. This is accomplished in a way which is similar to the procedure to link parameters described above. First, let’s create an independent variable. An IndependentVariable instance is created and added to the model like this:

```python
# Add the time as an independent variable of the model

time = IndependentVariable("time",0.0, unit='s')

my_model.add_independent_variable(time)
```

The IndependentVariable instance is inserted at the root of the model tree. In this case, can be accessed as:

```python
my_model.time.display()
```

We can now link any parameter to be a function of time, like this:

```python
# First define the function. In this case, a linear function ax+b
law = Line(a=-2,b=-0.02)
```

Now link the index of the sync. component of source_2 to be law(t)
i.e., $index = law(t) = a*t + b$

```python
my_model.link(my_model.source_2.spectrum.synchrotron.Powerlaw.index,
             time,
             law)
```

```python
my_model.save("time_dependent.yml",overwrite=True)

my_model = load_model("time_dependent.yml")
```

This would show the link:
```my_model.display()```

Now changing the value of time will change the value of the parameter
according to the law. For example, let's loop over 10 s and print
the value of the parameter


```python
# Reset time
my_model.time = 0.0

for i in range(10):

    my_model.time = my_model.time.value + 1.0

    print("At time %s s the value of the parameter is %s" % (my_model.time.value,
                    my_model.source_2.spectrum.synchrotron.Powerlaw.index.value))
```

Now plot the synch. spectrum of the source at different times
(you will need matplotlib for this)

```python
fig, ax = plt.subplots()
# Prepare 100 logarithmically distributed energies between 1 and 100 keV
energies = np.geomspace(1,100,100)

# Compute and plot the sync. spectrum every 10 s between 0 and 50 seconds

times = np.linspace(0,50,6)

my_model.time = 0.0

for tt in times:

    my_model.time = tt

    ax.loglog(energies, my_model.source_2.spectrum.synchrotron(energies),label='t = %s' % my_model.time.value)

ax.legend(loc=1,ncol=2)
ax.set_xlabel("Energy (keV)")
ax.set_ylabel(r"Differential flux (ph. cm$^{-2}$ s$^{-1}$ keV$^{-1}$)")
```

```python

```
