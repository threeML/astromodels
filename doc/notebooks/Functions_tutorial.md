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

# Functions tutorial

In astromodels functions can be used as spectral shapes for sources, or to describe time-dependence, phase-dependence, or links among parameters.

To get the list of available functions just do:

```python
%%capture
from astromodels import *
```

```python
list_functions()
```

If you need more info about a function, you can obtain it by using:

```python
Gaussian.info()
```

Note that you don’t need to create an instance in order to call the info() method.


## Creating functions

Functions can be created in two different ways. We can create an instance with the default values for the parameters like this:

```python
powerlaw_instance = Powerlaw()
```

or we can specify on construction specific values for the parameters:

```python
powerlaw_instance = Powerlaw(K=0.01, index=-2.2)
```

If you don’t remember the names of the parameters just call the .info() method as in powerlaw.info() as demonstrated above.


## Getting information about an instance

Using the ```.display()``` method we get a representation of the instance which exploits the features of the environment we are using. If we are running inside a IPython notebook, a rich representation with the formula of the function will be displayed (if available). Otherwise, in a normal terminal, the latex formula will not be rendered:

```python
powerlaw_instance.display()
```

It is also possible to get the text-only representation by simply printing the object like this:

```python
print(powerlaw_instance)
```

<div class="alert alert-info">

**Note:** the ```.display()``` method of an instance displays the current values of the parameters, while the .info() method demonstrated above (for which you don’t need an instance) displays the default values of the parameters.

</div>




## Modifying parameters

Modifying a parameter of a function is easy:

```python
# Modify current value

powerlaw_instance.K = 1.2

# Modify minimum
powerlaw_instance.K.min_value = 0.5

# Modify maximum
powerlaw_instance.K.max_value = 15

# We can also modify minimum and maximum at the same time
powerlaw_instance.K.bounds = (0.5, 15)

# Modifying the delta for the parameter
# (which can be used by downstream software for fitting, for example)
powerlaw_instance.K.delta = 0.25

# Fix the parameter
powerlaw_instance.K.fix = True

# or equivalently
powerlaw_instance.K.free = False

# Free it again
powerlaw_instance.K.fix = False

# or equivalently
powerlaw_instance.K.free = True

# We can verify what we just did by printing again the whole function as shown above,
# or simply printing the parameter:
powerlaw_instance.K.display()
```

## Using physical units

Astromodels uses the facility defined in astropy.units to make easier to convert between units during interactive analysis, when assigning to parameters. In order for functions to be aware of their units, they must be part of a ```Source``. Let's create one:

```python
powerlaw_instance = Powerlaw()

point_source = PointSource("my_point_source",ra=0,dec=0, spectral_shape=powerlaw_instance)
```

Now we can see the units

```python
powerlaw_instance.display()
```

```python
powerlaw_instance.x_unit
```

```python
powerlaw_instance.y_unit
```

```python
import astropy.units as u

# Express the differential flux at the pivot energy in 1 / (MeV cm2 s)

powerlaw_instance.K = (122.3 / (u.MeV * u.cm * u.cm * u.s))

print(powerlaw_instance.K)

# Express the differential flux at the pivot energy in 1 / (GeV m2 s)
powerlaw_instance.K = (122.3 / (u.GeV * u.m * u.m * u.s))


print(powerlaw_instance.K)
```

We see that astromodels does the unit conversion for us in the background!

However, astropy units are **very slow** and we would not want to deal with them or set parameters with units during a fit. Thus, if you do not specify units when setting a parameter, the value is assumed to have the units specified in the construction of the function. 

```python
%timeit powerlaw_instance.K = 1
```

```python
%timeit powerlaw_instance.K = (122.3 / (u.MeV * u.cm * u.cm * u.s))
```

As you can see using **astropy.units requires about 10x more than using a plain assignment**. In an interactive analysis you are unlikely to notice the difference, but if you use units in a loop or during a fit this slow-down will add up an become very noticeable. **Note that this is a feature of astropy.units, not of astromodels.**


## Composing functions

We can create arbitrary complex functions by combining “primitive” functions using the normal math operators:

```python
composite = Gaussian() + Powerlaw()

# Instead of the usual .display(), which would print all the many parameters,
# let's print just the description of the new composite functions:
print(composite.description)
```

These expressions can be as complex as needed. For example:

```python
crazy_function = 3 * Sin() + Powerlaw()**2 * (5+Gaussian()) / 3.0

print(crazy_function.description)
```

The numbers between ```{}``` enumerate the unique functions which constitute a composite function. This is useful because composite functions can be created starting from pre-existing instances of functions, in which case the same instance can be used more than once. For example:

```python
a_powerlaw = Powerlaw()
a_sin = Sin()

another_composite = 2 * a_powerlaw + (3 + a_powerlaw) * a_sin

print(another_composite.description)
```

In this case the same instance of a power law has been used twice. Changing the value of the parameters for “a_powerlaw” will affect also the second part of the expression. Instead, by doing this:

```python
another_composite2 = 2 * Powerlaw() + (3 + Powerlaw()) * Sin()

print(another_composite2.description)
```

we will end up with two independent sets of parameters for the two power laws. The difference can be seen immediately from the number of parameters of the two composite functions:

```python
print(len(another_composite.parameters)) # 6 parameters
print(len(another_composite2.parameters)) # 9 parameters
```

## Creating custom functions

```python

```
