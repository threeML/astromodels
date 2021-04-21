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

# Additional features for scripts and applications

In this document we describe some features of the astromodels package which are useful for non-interactive environment such as scripts or applications

First let’s import astromodels and let’s load a model from a file, which we will use as example:

```python
%%capture
from astromodels import *
my_model = load_model("my_model.yml")

```

## Get dictionaries of point and extended sources

If you don’t know the details (such as names) of the sources contained in the model, you can obtain dictionaries of point sources and extended sources like:

```python
point_sources = my_model.point_sources
extended_sources = my_model.extended_sources

# Print the names of the point sources
print(point_sources.keys())

# Print the names of the extended sources
print(extended_sources.keys())
```

You can use these dictionaries as usual. For example, you can loop over all point sources and print their position:

```python
for source_name, point_source in point_sources.items():

    print("The model contain point source %s at %s" % (source_name, point_source.position))
```

## Accessing components and spectral shapes with no previous information

Similarly you can access components and their spectral shapes (i.e., functions) without knowing the names in advance. A dictionary containing the components of a given source can be obtained with:

```python
components = my_model.source_2.components
print(components.keys())
```

So now we can loop over all the sources and print their components:

```python
for source_name, point_source in my_model.point_sources.items():

    print("Point source %s has components %s" % (source_name, point_source.components.keys()))
```

With a fully-qualified path, you would need to know the name of the function to access its parameters. Instead, you can use the generic name “shape”. For example these two statements point to the same function instance:

```python
my_model.source_1.spectrum.main.Powerlaw == my_model.source_1.spectrum.main.shape

```

Once you have a function instance, you can obtain a dictionary of its parameters as:

```python
parameters = my_model.source_1.spectrum.main.Powerlaw.parameters
print(parameters.keys())
```

Putting it all together, let’s loop over all sources in our model, then over each component in each source, then over each parameter in each component:

```python
for source_name, point_source in my_model.point_sources.items():

    print("Found source %s" % source_name)

    print("  Position of point source: %s" % point_source.position)

    for component_name, component in point_source.components.items():

        print("    Found component %s" % component_name)

        for parameter_name, parameter in component.shape.parameters.items():

            print("      Found parameter %s" % parameter_name)
```

Let’s now plot the differential flux between 1 and 100 keV of all components from all sources:



```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)

fig, ax = plt.subplots()


# Prepare 100 energies logarithmicall spaced between 1 and 100 keV
energies = np.geomspace(1,100,100)

# Now loop over all point sources and plot them
for source_name, point_source in my_model.point_sources.items():

    # Plot the sum of all components for this source

    ax.loglog(energies, point_source(energies),label=source_name)

    # If there is more than one component, plot them also separately

    if len(point_source.components) > 1:

        for component_name, component in point_source.components.items():

            ax.loglog(energies,component.shape(energies),'--',label="%s of %s" %(component_name, source_name))

# Add a legend
ax.legend(loc=0,frameon=False)

ax.set_xlabel("Energy (keV)")
ax.set_ylabel(r"Flux (ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$")
```

## Getting the path of an element and using it programmatically

Whenever you have an element from the model, you can get its fully-qualified path by using the .path property. This for example will print the path of all the parameters in the model:

```python
for source_name, point_source in my_model.point_sources.items():

    for component_name, component in point_source.components.items():

        for parameter_name, parameter in component.shape.parameters.items():

            print(parameter.path)
```

If you have a path of an element in a string, you can use it to access the element by using the [] operator of the Model class like this:

```python
my_path = 'source_2.spectrum.IC.Powerlaw.K'

K = my_model[my_path]

print(K)
```

## Alternative way of accessing the information in the model

We present here an alternative way to get information from the model without using dictionaries, and using instead source IDs. A source ID is just an ordinal number, separate for point sources and extended sources. Hence, the first point source has ID 0, the second point source has ID 1, and so on. Similarly, the first extended source has ID 0, the second has ID 1 and so on:

```python
# Get the number of point sources and of extended sources

n_pts = my_model.get_number_of_point_sources()
n_ext = my_model.get_number_of_extended_sources()

# Get the name of the first point source

print("The first point source is called %s" % my_model.get_point_source_name(0))
print("The second point source is called %s" % my_model.get_point_source_name(1))

# Of course you can achieve the same in a loop

for id in range(n_pts):

    print("Point source ID %s has name %s" % (id, my_model.get_point_source_name(id)))
```

Once you have the ID of a source, you can obtain information about it with these methods of the Model class:

```python
src_id = 1

src_name = my_model.get_point_source_name(src_id)

ra, dec = my_model.get_point_source_position(src_id) # This will always return ra,dec

# Prepare 100 energies logarithmically spaced between 1 and 100 keV
energies = np.logspace(0,2,100)

differential_flux = my_model.get_point_source_fluxes(src_id, energies)

# Similar methods exist for extended sources (to be completed)
```

Depending on your application you might find these methods more convenient that interrogating the sources directly. Note however that some features are not available through this interface. For example, it is not possible to get information about different components through these methods.

```python

```
