---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.7.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- #region -->
# Point sources #

In astromodels, a point source is described by its position in the sky
and its spectral features.


## Creating a point source

A simple source with a power law spectrum can be created like this, using J2000 R.A. and Dec (ICRS), which is the default coordinate system:
<!-- #endregion -->

```python
from astromodels import *
simple_source_icrs = PointSource('simple_source', ra=123.2, dec=-13.2, spectral_shape=Powerlaw())
```

We can also use Galactic coordinates:

```python
simple_source_gal = PointSource('simple_source', l=234.320573, b=11.365142, spectral_shape=Powerlaw())
```

As spectral shape we can use any function or any composite function (see
"Creating and modifying functions")


## Getting info about a point source

Info about a point source can be obtained with the
`.display()` method (which will use the richest representation available),
or by printing it which will display a text-only representation:


```python
simple_source_icrs.display()
```

```python
print(simple_source_icrs)
```

As you can see we have created a point source with one component automatically named "main", with a power law spectrum, at the
specified position.


## Converting between coordinate systems

By default the coordinates of the point source are displayed in the same
system used during creation. However, you can always obtain R.A, Dec or
L,B like this:

```python
simple_source_icrs.position.display()
l = simple_source_icrs.position.get_l()
b = simple_source_icrs.position.get_b()
print(l,b)

simple_source_gal.position.display()
ra = simple_source_gal.position.get_ra()
dec = simple_source_gal.position.get_dec()
print(ra,dec)
```

For more control on the output and many more options, such as transform
to local frames or other equinoxes, you can obtain an instance of
[astropy.coordinates.SkyCoord](http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html) by using the `sky_coord` property of the
position object:

```python
sky_coord_instance = simple_source_icrs.position.sky_coord
ra = sky_coord_instance.transform_to('icrs').ra
dec = sky_coord_instance.transform_to('icrs').dec
print(ra.deg)
```

## Gotcha while accessing coordinates

Please note that using `get_ra()` and `.ra` (or the equivalent methods for
the other coordinates) is not the same. While `get_ra()` will always
return a single float value corresponding to the R.A. of the source, the
`.ra` property will exist only if the source has been created using R.A,
Dec as input coordinates and will return a Parameter instance:

```python
parameter_ra = simple_source_icrs.position.ra
parameter_dec = simple_source_icrs.position.dec

print( type(parameter_ra) )
parameter_ra.display()
parameter_dec.display()
```

The following would instead throw `AttributeError`, since `simple_source_icrs` was instanced using R.A. and Dec. and hence does not have the `l`, `b` parameters:

```python
try:
    print( simple_source_icrs.position.l)
except Exception as e:
    print(e)
```

In all cases, independently on how the source was instanced, you can obtain the *values* of coordinates in degrees
as  floating point numbers using `get_ra()`, `get_dec()`, `get_l()`, `get_b()`. However, you can only directly *assign* coordinates in the same system that the source direction was originally created, e.g.:


```python
simple_source_icrs.position.display()
simple_source_icrs.position.ra =   simple_source_icrs.position.ra.value + 1.0
simple_source_icrs.position.dec = simple_source_icrs.position.dec.value - 1.0
simple_source_icrs.position.display()
```

## Fitting the source position

Source coordinates, like any parameters, can be set to be free or fixed during the fit. By default, coordinates are set to be fixed. If you would like to fit them as free parameters during the likelihood fit, they can be freed as any other parameter. Note that `param.free = True` and `param.fix = False` are equivalent. 

```python
print("Free parameters (before freeing position):", simple_source_icrs.free_parameters.keys())
simple_source_icrs.position.ra.free = True
simple_source_icrs.position.dec.fix = False
print("Free parameters (after freeing position):", simple_source_icrs.free_parameters.keys())

```

For a source created in Galactic coordinates, instead use the following:

```python
print("Free parameters (before freeing position):", simple_source_gal.free_parameters.keys())
simple_source_gal.position.l.free = True
simple_source_gal.position.b.fix = False
print("Free parameters (after freeing position):", simple_source_gal.free_parameters.keys())

```

By default, the allowed range for the Right Ascension is from 0˚ to 360˚ and allowed declination values range from -90˚ to 90˚. **If fitting the source position, it is strongly recommended to restrict the coordinates to be inside the region of interest (ROI) at all times.** The source moving far enough from the ROI during the fit can lead to issues such as the minimizer getting "stuck" due to the likelihood surface being flat. For example:

```python
simple_source_icrs.position.ra.bounds = ( simple_source_icrs.position.ra.value - 5.0, simple_source_icrs.position.ra.value + 5.0 )
simple_source_icrs.position.dec.bounds = ( simple_source_icrs.position.dec.value - 5.0, simple_source_icrs.position.dec.value + 5.0 )
simple_source_icrs.position.ra.display()
simple_source_icrs.position.dec.display()

```

## "Calling" a point source

Both the point source object itself as well as the compontents are callable functions who take as argument(s) an array of energies and return the differential flux dN/dE at those energies. Energies can be provided with or without units. **If no units are provided, energies are assumed to be in keV and fluxes are returned in 1/(cm2 keV s).**

```python
from astropy import units as u
E = [1, 10, 100]*u.keV

print("Energy in keV:")
print( "With units:", simple_source_icrs(E) )
print( "Without units:", simple_source_icrs(E.value) )

print( "With units:", simple_source_icrs.spectrum.main.shape(E) )
print( "Without units:", simple_source_icrs.spectrum.main.shape(E.value) )

print("")
print("Energy in TeV:")
E_TeV = E.to(u.TeV)
print( "With units:", simple_source_icrs(E_TeV) )
print( "With units:", simple_source_icrs(E_TeV).to(1/u.cm**2/u.TeV/u.s) )
print( "Without units:", simple_source_icrs(E_TeV.value) )



```

```python

```
