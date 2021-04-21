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

# Extended source tutorial


## Introduction

Extended sources with constant (not energy-dependent) morphology can be instanced by providing a spectral shape (1D function) and a spatial shape (2D function):

```python
from astromodels import *
disk_source = ExtendedSource('simple_source', spectral_shape=Powerlaw(), spatial_shape = Disk_on_sphere() )
disk_source.display()
```

Instead of the `position` property, extended sources have a `spatial_shape` property, which contains the function that was specified as the `shape` argument earlier. Information can be printed either with the regular `print` function or using `.display()`:

```python
print("Name", disk_source.spatial_shape.name)
print("Class", type(disk_source.spatial_shape) )
disk_source.spatial_shape.display()
```

Spatial shape parameters can be listed, accessed, changed as any function parameters, for example:

```python
import astropy.units as u
print( disk_source.spatial_shape.parameters.keys() )
print( disk_source.parameters.keys() )
disk_source.spatial_shape.radius.display()
disk_source.spatial_shape.radius = 0.3 * u.deg
disk_source.spatial_shape.radius.bounds = (0.1, 1.0)*u.deg
disk_source.spatial_shape.radius.display()

```

## Available shapes


Here's a list of available 2D functions:

```python
print([c.__name__, for c in functions.function.Function2D.__subclasses__()])
```

Check the API or use the `display()` function to find out more about the shapes and their free parameters. 


### "Moveable" Functions

`Gaussian_on_sphere`, `Asymm_Gaussian_on_sphere`, `Disk_on_sphere`, `Ellipse_on_sphere`, and `Power_law_on_sphere` each have parameters called `lat0` and `lon0`, describing the longitude (R.A.) and latitude (Declination) of the function's "center" (center of mass/center of symmetry). These parameters can be set free to allow the source position to vary during the fit, or fixed to a known position. **If fitting the source position, it is strongly recommended to restrict these parameters so that source overlaps the region of interest (ROI) at all times.** The source moving far enough from the ROI during the fit can lead to issues such as the minimizer getting "stuck" due to the likelihood surface being flat.

Each of these functions also has other parameter(s) related to the source's extend and (if applicable) its orientation. 

All of these functions are normalized so that their integral is 1. The 1D `spectrum` function of an extended source can thus be interpreted as the spectrum that would be observed by a non-imaging instrument with a field of view much larger than the source in question.



### Fixed-position Functions

`Latitude_galactic_diffuse` and `SpatialTemplate_2D` do not have parameters named `lat0` and `lon0`. Their positions are always fixed throughout the fit. 

`Latitude_galactic_diffuse` is a function that is constant in Galactic longitude between `l_min` and `l_max` and follows a Gaussian shape in Galactic latitude, with a "width" parameter (`sigma_b`). 

```python
diff = Latitude_galactic_diffuse()
diff.display()
```

`SpatialTemplate_2D` is designed to read in a user-provided fits file with an image (in WCS coordinates) using the `load_file()` function. The function value will be 0 outside the WCS and equal to the value of the pixel containing the given coordinates inside the WCS. The prodided template should be in units of 1/deg2 and normalized so that its integral is 1. If that is not the case, the normalization parameter `K` may be set accordingly so that the overall function is normalized as expected.

```python
temp = SpatialTemplate_2D()
temp.display()
```

## "Calling" extended sources

Extended sources can be called as functions. They take 3 arguments: A list or `np.array` of right ascensions, a list of declinations, and a list of energies. The first two lists must have the same dimensions. The result contains the value(s) of the double differential flux dN/dE/dOmega at all combinations of coordinates and energies. See the example below. There, we picked two positions (one at the center of the disk source, one outside of it) and three energy values. **If no units are provided, energies are assumed to be in keV, coordinates in degrees, and fluxes are returned in 1/(cm2 keV s deg2).**

```python
import numpy as np
ra_center, dec_center = disk_source.spatial_shape.lon0.value, disk_source.spatial_shape.lat0.value

ra = [ra_center, ra_center+2]*u.deg
dec = [dec_center, dec_center+2]*u.deg

E = [1, 10, 100]*u.keV

print( "Whole source:" )
print( disk_source(ra, dec, E) )

print ("")
print ("Spatial only:" )
spatial = disk_source.spatial_shape(ra, dec)
print( spatial )

print ("")
print ("Spectral only:")
spectral = disk_source.spatial_shape.get_total_spatial_integral( E ) * disk_source.spectrum.main.shape( E )
print ( spectral ) 

print ("")
print ( "Spatial transposed times spectral")
print( np.atleast_2d(spatial).T  * np.atleast_2d(spectral) )
```

## Caveats and gotchas

### Large extended sources

The Gaussian shapes provided here (`Gaussian_on_sphere`, `Asymm_Gaussian_on_sphere`) are only valid as long as their width (`sigma`/`a` parameter) is much smaller than 360˚. For very large extended sources, we would need to implement the [Kent distribution](https://en.wikipedia.org/wiki/Kent_distribution), the equivalent of a Gaussian on the sphere. 

### Small extended sources

Most threeML plugins supporting imaging instruments perform a convolution of the `astromodels` source(s) with the instrument's point spread function (PSF). These convolutions are typically performed on a grid, and features in the source morphology smaller than the grid size are lost. One should take care that the extent of the source(s) in the model is constrained to be no smaller than the grid spacing used by the plugin.

### Source support

Each 2D function class has a member function called `get_boundaries`, which returns the (minimal) range in RA and Dec over which the convolution with the PSF has to be computed (see example below). The exact calculation depends on the function but the range typically depends on the centroid and the **maximum value of the width/radius parameter**. That means that even if the radius or gaussian width parameter is fixed during the fit, it is helpful to set a maximum that is close to the fixed value. This reduces computation time during the likelihood fit. 

```python
disk_source.spatial_shape.lon0 = 80*u.deg
disk_source.spatial_shape.lat0 = 20*u.deg
disk_source.spatial_shape.radius.bounds = (0.1, 2)*u.deg
bounds =  (disk_source.spatial_shape.get_boundaries() )
print ("RA range:", bounds[0])
print ("Dec range:", bounds[1])
```

## Energy-Dependent Morphology

`astromodels` functions support energy-dependent morphology in two different ways. 

1. Extended sources can be instanced with a 3D function spatial shape and a separate energy spectrum. In that case, the 3D function is expected to be a function of RA, Dec, Energy (in that order) and have units of 1/deg2. The double-differential flux dN/dE/dOmega is again given by the product of the spatial and spectral parts. See for example [`Continuous_injection_diffusion`](https:link).
2. Extended sources can be instanced with just a 3D function and no energy spectrum. In that case, the 3D function is interpreted as the double-differential flux dN/dE/dOmega.

```python

```
