
Point sources
=============

In astromodels a point source is described by its position in the sky
and its spectral features.

Creating a point source
-----------------------

A simple source with a power law spectrum can be created like this:

.. code:: python

    from astromodels import *

.. code:: python

    # Using J2000 R.A. and Dec (ICRS), which is the default coordinate system:
    
    simple_source_icrs = PointSource('simple_source', ra=123.2, dec=-13.2, spectral_shape=powerlaw())

We can also use Galactic coordinates:

.. code:: python

    simple_source_gal = PointSource('simple_source', l=234.320573, b=11.365142, spectral_shape=powerlaw())

As spectral shape we can use any function or any composite function (see
"Creating and modifying functions")

Getting info about a point source
---------------------------------

Info about a point source can easily be obtained with the usual
.display() method (which will use the richest representation available),
or by printing it which will display a text-only representation:

.. code:: python

    simple_source_icrs.display()
    
    # or print(simple_source_icrs) for a text-only representation



.. raw:: html

    <ul>
    
    <li>simple_source (point source): 
    <ul>
    
    <li>position: 
    <ul>
    
    <li>ra: 123.2</li>
    
    <li>dec: -13.2</li>
    
    <li>equinox: J2000</li>
    
    </ul>
    
    </li>
    
    <li>spectrum: 
    <ul>
    
    <li>main: 
    <ul>
    
    <li>powerlaw: 
    <ul>
    
    <li>logK: 0.0</li>
    
    <li>piv: 1.0</li>
    
    <li>index: -2.0</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>



As you can see we have created a point source with one component (see
below) automatically named "main", with a power law spectrum, at the
specified position.

Converting between coordinates systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default the coordinates of the point source are displayed in the same
system used during creation. However, you can always obtain R.A, Dec or
L,B like this:

.. code:: python

    l = simple_source_icrs.position.get_l()
    b = simple_source_icrs.position.get_b()
    ra = simple_source_gal.position.get_ra()
    dec = simple_source_gal.position.get_dec()

For more control on the output and many more options, such as transform
to local frames or other equinoxes, you can obtain an instance of
astropy.coordinates.SkyCoord by using the sky\_coord property of the
position object:

.. code:: python

    # Refer to the transform_to() method of the astropy.coordinates.SkyCoord class:
    # http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html
    
    # For example, get the ICRS position for the source defined in Galactic coordinates:
    
    sky_coord_instance = simple_source_icrs.position.sky_coord
    
    ra = sky_coord_instance.transform_to('icrs').ra
    dec = sky_coord_instance.transform_to('icrs').dec
    
    print ra.deg


.. parsed-literal::

    123.2


Gotcha while accessing coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please note that using get\_ra() and .ra (or the equivalent methods for
the other coordinates) is not the same. While get\_ra() will always
return a single float value corresponding to the R.A. of the source, the
.ra property will exist only if the source has been created using R.A,
Dec as input coordinates and will return a Parameter instance:

.. code:: python

    # These will return two Parameter instances corresponding to the parameters ra and dec
    # NOT the corresponding floating point numbers:
    parameter_ra = simple_source_icrs.position.ra
    parameter_dec = simple_source_icrs.position.dec
    
    # This would instead throw AttributeError, since simple_source_icrs was instanced using
    # R.A. and Dec. and hence does not have the l,b parameters:
    # error = simple_source_icrs.position.l
    # error = simple_source_icrs.position.b
    
    # Similarly this will throw AttributeError, because simple_source_gal was instanced using
    # Galactic coordinates:
    # error = simple_source_gal.position.ra
    # error = simple_source_gal.position.dec
    
    # In all cases, independently on how the source was instanced, you can obtain the coordinates
    # as normal floating point numbers using:
    ra1 = simple_source_icrs.position.get_ra()
    dec1 = simple_source_icrs.position.get_dec()
    l1 = simple_source_icrs.position.get_l()
    b1 = simple_source_icrs.position.get_b()
    
    ra2 = simple_source_gal.position.get_ra()
    dec2 = simple_source_gal.position.get_dec()
    l2 = simple_source_gal.position.get_l()
    b2 = simple_source_gal.position.get_b()

Multi-component sources
-----------------------

A multi-component source is a point source which has different spectral
components. For example, in a Gamma-Ray Bursts you can have a
Synchrotron component and a Inverse Compton component, which come from
different zones and are described by different spectra. Depending on the
needs of your analysis, you might model this situation using a single
component constituted by the sum of the two spectra, or you might want
to model them independently. Also, each components has its own
polarization, which can be useful when studying polarized sources (to be
implemented). The latter choice allows you to measure for instance the
fluxes from the two components independently. Representing a source with
more than one component is easy in astromodels:

.. code:: python

    # Create the two different components 
    #(of course the shape can be any function, or any composite function)
    
    component1 = SpectralComponent('synchrotron',shape=powerlaw())
    component2 = SpectralComponent('IC',shape=powerlaw())
    
    # Create a multi-component source
    multicomp_source = PointSource('multicomp_source', ra=123.2, dec=-13.2, components=[component1,component2])
    
    multicomp_source.display()



.. raw:: html

    <ul>
    
    <li>multicomp_source (point source): 
    <ul>
    
    <li>position: 
    <ul>
    
    <li>ra: 123.2</li>
    
    <li>dec: -13.2</li>
    
    <li>equinox: J2000</li>
    
    </ul>
    
    </li>
    
    <li>spectrum: 
    <ul>
    
    <li>synchrotron: 
    <ul>
    
    <li>powerlaw: 
    <ul>
    
    <li>logK: 0.0</li>
    
    <li>piv: 1.0</li>
    
    <li>index: -2.0</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    <li>IC: 
    <ul>
    
    <li>powerlaw: 
    <ul>
    
    <li>logK: 0.0</li>
    
    <li>piv: 1.0</li>
    
    <li>index: -2.0</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>



Modifying features of the source and modify parameters of its spectrum
----------------------------------------------------------------------

Starting from the source instance you can modify any of its components,
or its position, in a straightforward way:

.. code:: python

    # Change position
    
    multicomp_source.position.ra = 124.5
    multicomp_source.position.dec = -11.5
    
    # Change values for the parameters
    multicomp_source.spectrum.synchrotron.powerlaw.logK = -1.2
    
    multicomp_source.spectrum.IC.powerlaw.index = -1.0
    
    # To avoid having to write that much, you can create a "shortcut" for a function
    po = multicomp_source.spectrum.synchrotron.powerlaw
    
    # Now you can modify its parameters more easily 
    # (see "Creating and modifying functions" for more info on what you can to with a parameter)
    po.logK = -1.3
    po.logK.min_value = -3.0
    
    # GOTCHA
    # Creating a shortcut directly to the parameter will not work:
    
    # p1 = multicomp_source.spectrum.synchrotron.powerlaw.logK
    # p1 = -1.3 # this does NOT change the value of logK, but instead assign -1.3 to p1 (i.e., destroy the shortcut)
    
    # However you can change the value of p1 like this:
    # p1.value = -1.3 # This will work
    
    multicomp_source.display()



.. raw:: html

    <ul>
    
    <li>multicomp_source (point source): 
    <ul>
    
    <li>position: 
    <ul>
    
    <li>ra: 124.5</li>
    
    <li>dec: -11.5</li>
    
    <li>equinox: J2000</li>
    
    </ul>
    
    </li>
    
    <li>spectrum: 
    <ul>
    
    <li>synchrotron: 
    <ul>
    
    <li>powerlaw: 
    <ul>
    
    <li>logK: -1.3</li>
    
    <li>piv: 1.0</li>
    
    <li>index: -2.0</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    <li>IC: 
    <ul>
    
    <li>powerlaw: 
    <ul>
    
    <li>logK: 0.0</li>
    
    <li>piv: 1.0</li>
    
    <li>index: -1.0</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>



