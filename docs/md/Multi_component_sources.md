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

# Multi-component sources


A multi-component source is a point source (or extended source) which has different spectral
components. For example, in a Gamma-Ray Bursts you can have a
Synchrotron component and an Inverse Compton component, which come from
different zones and are described by different spectra. Depending on the
needs of your analysis, you might model this situation using a single
component constituted by the sum of the two spectra, or you might want
to model them independently. Each component has its own
polarization, which can be useful when studying polarized sources (to be
implemented). The latter choice allows you to measure for instance the
fluxes from the two components independently. See below for an example of how to create a source with two named spectral components (of course, the spectral shapes can be arbitrary).

**Sources with multiple *spatial* components are not currently supported. As a workaround, you can create two functions with the same spectral shape and link any relevant parameters. There is no restriction against overlapping extended sources.**

```python
from astromodels import *
component1 = SpectralComponent('synchrotron',shape=Powerlaw())
component2 = SpectralComponent('IC',shape=Powerlaw())
multicomp_source = PointSource('multicomp_source', ra=123.2, dec=-13.2, components=[component1,component2])
    
multicomp_source.display()

```

## Modifying features of the source and modify parameters of its spectrum

Starting from the source instance you can modify any of its components,
or its position, in a straightforward way.

Changing position:

```python
multicomp_source.position.ra = 124.5
multicomp_source.position.dec = -11.5
```

Change values for the parameters:

```python
multicomp_source.spectrum.synchrotron.shape.K = 1e-2
multicomp_source.spectrum.IC.shape.index = -1.0
multicomp_source.display()
```

Spectral components can be assigned to python variables to modify several parameters without too much repetition:

```python
po = multicomp_source.spectrum.synchrotron.shape
po.K = 0.1/u.keV/u.cm**2/u.s
po.K.min_value = 1e-10/u.keV/u.cm**2/u.s
multicomp_source.display()
```

Be careful when creating a shortcut directly to a parameter.

```python
p1 = multicomp_source.spectrum.synchrotron.shape.K
print (type(p1), p1)
p1 = 0.3
print (type(p1), p1)

```

This did **not** change the value of K, but instead assigned the float 0.3 to p1 (i.e., destroy the shortcut).   
 
However you can change the value of p1 like this:

```python
p1 = multicomp_source.spectrum.synchrotron.shape.K
p1.value = 0.5
multicomp_source.display()
```

```python

```

```python

```
