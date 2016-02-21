Quick tutorial
==============

.. raw:: html

   <h1>

Quick tutorial

.. raw:: html

   </h1>

In this quick tutorial we show how to build a simple model with two
point sources, how to save it for later use, and re-load it back. We
will also plot the spectra of the two point sources, with their
components

First let's import the modules we need.

.. code:: python

    from astromodels import *
    
    # We also import astropy units to show the unit-conversion
    # feature
    import astropy.units as u

Now let's define the first point source. In astromodels a point source
is characterized by a position (which can be given in Equatorial or
Galactic coordiantes) and one or more ``components``. For example, a
source featuring both a synchrotron and a Inverse Compton emission would
have two components.

Let's start from a point source located at (R.A., Dec) = (125.6, -75.3)
deg (ICRS, J2000), and having just one component with a power law
spectrum.

The spectrum can be defined as follows:

.. code:: python

    # Let's start by defining the spectrum
    
    pts1_shape = powerlaw()
    
    # Have a look at it
    
    pts1_shape.display()



.. raw:: html

    <ul>
    
    <li>description: A simple power-law with normalization expressed as a logarithm</li>
    
    <li>formula: $ \frac{dN}{dx} = 10^{logK}~\frac{x}{piv}^{index} $</li>
    
    <li>parameters: 
    <ul>
    
    <li>logK: 
    <ul>
    
    <li>value: 0</li>
    
    <li>min_value: -40</li>
    
    <li>max_value: 40</li>
    
    <li>delta: 0.0</li>
    
    <li>free: True</li>
    
    <li>unit: dex(1 / (cm2 keV s))</li>
    
    </ul>
    
    </li>
    
    <li>piv: 
    <ul>
    
    <li>value: 1</li>
    
    <li>min_value: None</li>
    
    <li>max_value: None</li>
    
    <li>delta: 0.1</li>
    
    <li>free: False</li>
    
    <li>unit: keV</li>
    
    </ul>
    
    </li>
    
    <li>index: 
    <ul>
    
    <li>value: -2</li>
    
    <li>min_value: -10</li>
    
    <li>max_value: 10</li>
    
    <li>delta: -0.2</li>
    
    <li>free: True</li>
    
    <li>unit: </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>



We can change the value of the parameters as:

.. code:: python

    pts1_shape.logK = 1.2
    pts1_shape.index = -2.1

We can also provide the parameters in different units, and astromodels
will convert automatically to the one used internally and reported
above:

.. code:: python

    # Express the differential flux at the pivot energy in 1 / (MeV cm2 s)
    
    pts1_shape.logK = 122.3 / (u.MeV * u.cm * u.cm * u.s)
    
    # Verify how things have changed
    
    pts1_shape.display()


.. parsed-literal::

    /home/giacomov/software/canopy-env/lib/python2.7/site-packages/astromodels-0.1-py2.7.egg/astromodels/functions/function.py:559: WarningUnitsAreSlow: Using units is convenient but slow. Do not use them during computing-intensive work.
      WarningUnitsAreSlow)



.. raw:: html

    <ul>
    
    <li>description: A simple power-law with normalization expressed as a logarithm</li>
    
    <li>formula: $ \frac{dN}{dx} = 10^{logK}~\frac{x}{piv}^{index} $</li>
    
    <li>parameters: 
    <ul>
    
    <li>logK: 
    <ul>
    
    <li>value: -0.912573542964</li>
    
    <li>min_value: -40</li>
    
    <li>max_value: 40</li>
    
    <li>delta: 0.0</li>
    
    <li>free: True</li>
    
    <li>unit: dex(1 / (cm2 keV s))</li>
    
    </ul>
    
    </li>
    
    <li>piv: 
    <ul>
    
    <li>value: 1</li>
    
    <li>min_value: None</li>
    
    <li>max_value: None</li>
    
    <li>delta: 0.1</li>
    
    <li>free: False</li>
    
    <li>unit: keV</li>
    
    </ul>
    
    </li>
    
    <li>index: 
    <ul>
    
    <li>value: -2.1</li>
    
    <li>min_value: -10</li>
    
    <li>max_value: 10</li>
    
    <li>delta: -0.2</li>
    
    <li>free: True</li>
    
    <li>unit: </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>



The warning you see the first time you execute this command warns you
that using conversions like these is slow, so it is fine in an
interactive environment but you need to be careful if you are using this
feature in a computing-intensive situation (like during fitting).

Now let's define the point source:

.. code:: python

    pts1 = PointSource('source_1', ra=125.6, dec=-75.3, 
                       spectral_shape=pts1_shape)
    
    # Get some info about what we just created
    pts1.display()



.. raw:: html

    <ul>
    
    <li>source_1 (point source): 
    <ul>
    
    <li>position: 
    <ul>
    
    <li>ra: 125.6</li>
    
    <li>dec: -75.3</li>
    
    <li>equinox: J2000</li>
    
    </ul>
    
    </li>
    
    <li>components: 
    <ul>
    
    <li>main: 
    <ul>
    
    <li>shape: 
    <ul>
    
    <li>powerlaw: 
    <ul>
    
    <li>logK: -0.912573542964</li>
    
    <li>piv: 1</li>
    
    <li>index: -2.1</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>



Now let's define another source, this time at Galactic Coordinates l =
11.25, b = -22.5, and with two spectral components:

.. code:: python

    # Another point source with two spectral components,
    
    # Define component 1, also showing another way to
    # instance the spectrum, which allows to specify
    # directly the parameters' values
    
    spectrum1 = powerlaw(logK=0.2, index=-0.75)
    component1 = SpectralComponent('synchrotron',spectrum1)
    
    spectrum2 = powerlaw(logK=-1, index=-1.7)
    component2 = SpectralComponent('Inverse_Compton',spectrum2)
    
    point_source2 = PointSource('source_2', l=11.25, b=-22.5, components=[component1,component2])
    
    # Have a look at what we just created
    
    point_source2.display()



.. raw:: html

    <ul>
    
    <li>source_2 (point source): 
    <ul>
    
    <li>position: 
    <ul>
    
    <li>l: 11.25</li>
    
    <li>b: -22.5</li>
    
    <li>equinox: J2000</li>
    
    </ul>
    
    </li>
    
    <li>components: 
    <ul>
    
    <li>synchrotron: 
    <ul>
    
    <li>shape: 
    <ul>
    
    <li>powerlaw: 
    <ul>
    
    <li>logK: 0.2</li>
    
    <li>piv: 1</li>
    
    <li>index: -0.75</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    <li>Inverse_Compton: 
    <ul>
    
    <li>shape: 
    <ul>
    
    <li>powerlaw: 
    <ul>
    
    <li>logK: -1</li>
    
    <li>piv: 1</li>
    
    <li>index: -1.7</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>



Now let's build the final model, which comprises our two sources:

.. code:: python

    # Build a model with the two point sources
    
    my_model = Model(pts1, point_source2)
    
    # Let's have a look at what we just created
    
    my_model.display()



.. raw:: html

    Point sources: source_1,source_2<br><br>Extended sources: (none)<br><br>Parameters:<br><br>source_1<br><table id="table71521104">
    <thead><tr><th>name</th><th>value</th><th>min_value</th><th>max_value</th><th>delta</th><th>free</th><th>unit</th></tr></thead>
    <tr><td>position.ra</td><td>125.6</td><td>0.0</td><td>360.0</td><td>12.56</td><td>False</td><td></td></tr>
    <tr><td>position.dec</td><td>-75.3</td><td>-90.0</td><td>90.0</td><td>-7.53</td><td>False</td><td></td></tr>
    <tr><td>main.powerlaw.logK</td><td>-0.912573542964</td><td>-40</td><td>40</td><td>0.0</td><td>True</td><td>dex(1 / (cm2 keV s))</td></tr>
    <tr><td>main.powerlaw.piv</td><td>1.0</td><td>None</td><td>None</td><td>0.1</td><td>False</td><td>keV</td></tr>
    <tr><td>main.powerlaw.index</td><td>-2.1</td><td>-10</td><td>10</td><td>-0.2</td><td>True</td><td></td></tr>
    </table><br><br>source_2<br><table id="table71524304">
    <thead><tr><th>name</th><th>value</th><th>min_value</th><th>max_value</th><th>delta</th><th>free</th><th>unit</th></tr></thead>
    <tr><td>position.l</td><td>11.25</td><td>0.0</td><td>360.0</td><td>1.125</td><td>False</td><td></td></tr>
    <tr><td>position.b</td><td>-22.5</td><td>-90.0</td><td>90.0</td><td>-2.25</td><td>False</td><td></td></tr>
    <tr><td>synchrotron.powerlaw.logK</td><td>0.2</td><td>-40</td><td>40</td><td>0.0</td><td>True</td><td>dex(1 / (cm2 keV s))</td></tr>
    <tr><td>synchrotron.powerlaw.piv</td><td>1.0</td><td>None</td><td>None</td><td>0.1</td><td>False</td><td>keV</td></tr>
    <tr><td>synchrotron.powerlaw.index</td><td>-0.75</td><td>-10</td><td>10</td><td>-0.2</td><td>True</td><td></td></tr>
    <tr><td>Inverse_Compton.powerlaw.logK</td><td>-1.0</td><td>-40</td><td>40</td><td>0.0</td><td>True</td><td>dex(1 / (cm2 keV s))</td></tr>
    <tr><td>Inverse_Compton.powerlaw.piv</td><td>1.0</td><td>None</td><td>None</td><td>0.1</td><td>False</td><td>keV</td></tr>
    <tr><td>Inverse_Compton.powerlaw.index</td><td>-1.7</td><td>-10</td><td>10</td><td>-0.2</td><td>True</td><td></td></tr>
    </table><br>


Now let's have a look at different ways to change parameters values or
characteristics:

.. code:: python

    # Method 1: universal method. When you have available an instance of a model
    # you can refer to a parameter with the 'extended' syntax:
    
    my_model.source_2.synchrotron.powerlaw.logK = 0.3
    
    # Method 2: universal method 2. You can obtain the same result without
    # explicitly using the name of the spectral function, by using the
    # generic name 'shape':
    
    my_model.source_2.synchrotron.shape.logK = 0.3
    
    # Method 3: using instances created previously.
    
    point_source2.synchrotron.shape.logK = -0.3
    component2.shape.logK = -0.3
    spectrum2.logK = -0.3
    
    my_model.display()



.. raw:: html

    Point sources: source_1,source_2<br><br>Extended sources: (none)<br><br>Parameters:<br><br>source_1<br><table id="table71524240">
    <thead><tr><th>name</th><th>value</th><th>min_value</th><th>max_value</th><th>delta</th><th>free</th><th>unit</th></tr></thead>
    <tr><td>position.ra</td><td>125.6</td><td>0.0</td><td>360.0</td><td>12.56</td><td>False</td><td></td></tr>
    <tr><td>position.dec</td><td>-75.3</td><td>-90.0</td><td>90.0</td><td>-7.53</td><td>False</td><td></td></tr>
    <tr><td>main.powerlaw.logK</td><td>-0.912573542964</td><td>-40</td><td>40</td><td>0.0</td><td>True</td><td>dex(1 / (cm2 keV s))</td></tr>
    <tr><td>main.powerlaw.piv</td><td>1.0</td><td>None</td><td>None</td><td>0.1</td><td>False</td><td>keV</td></tr>
    <tr><td>main.powerlaw.index</td><td>-2.1</td><td>-10</td><td>10</td><td>-0.2</td><td>True</td><td></td></tr>
    </table><br><br>source_2<br><table id="table71558672">
    <thead><tr><th>name</th><th>value</th><th>min_value</th><th>max_value</th><th>delta</th><th>free</th><th>unit</th></tr></thead>
    <tr><td>position.l</td><td>11.25</td><td>0.0</td><td>360.0</td><td>1.125</td><td>False</td><td></td></tr>
    <tr><td>position.b</td><td>-22.5</td><td>-90.0</td><td>90.0</td><td>-2.25</td><td>False</td><td></td></tr>
    <tr><td>synchrotron.powerlaw.logK</td><td>-0.3</td><td>-40</td><td>40</td><td>0.0</td><td>True</td><td>dex(1 / (cm2 keV s))</td></tr>
    <tr><td>synchrotron.powerlaw.piv</td><td>1.0</td><td>None</td><td>None</td><td>0.1</td><td>False</td><td>keV</td></tr>
    <tr><td>synchrotron.powerlaw.index</td><td>-0.75</td><td>-10</td><td>10</td><td>-0.2</td><td>True</td><td></td></tr>
    <tr><td>Inverse_Compton.powerlaw.logK</td><td>-0.3</td><td>-40</td><td>40</td><td>0.0</td><td>True</td><td>dex(1 / (cm2 keV s))</td></tr>
    <tr><td>Inverse_Compton.powerlaw.piv</td><td>1.0</td><td>None</td><td>None</td><td>0.1</td><td>False</td><td>keV</td></tr>
    <tr><td>Inverse_Compton.powerlaw.index</td><td>-1.7</td><td>-10</td><td>10</td><td>-0.2</td><td>True</td><td></td></tr>
    </table><br>


Now we can save the model for later use:

.. code:: python

    # Save the model to a file
    
    my_model.save('my_model.yml', overwrite=True)

The content of the file is YAML code, which is human-readable and very
easy to understand. Let's have a look:

.. code:: python

    with open('my_model.yml') as yaml_file:
        
        print("".join(yaml_file.readlines()))


.. parsed-literal::

    source_1:
    
      point source:
    
        position:
    
          ra: {value: 125.6, min_value: 0.0, max_value: 360.0, delta: 12.56, free: false,
    
            unit: ''}
    
          dec: {value: -75.3, min_value: -90.0, max_value: 90.0, delta: -7.53, free: false,
    
            unit: ''}
    
          equinox: J2000
    
        spectrum:
    
          main:
    
            shape:
    
              powerlaw:
    
                logK: {value: -0.9125735429637145, min_value: -40, max_value: 40, delta: 0.0,
    
                  free: true, unit: dex(1 / (cm2 keV s))}
    
                piv: {value: 1, min_value: null, max_value: null, delta: 0.1, free: false,
    
                  unit: keV}
    
                index: {value: -2.1, min_value: -10, max_value: 10, delta: -0.2, free: true,
    
                  unit: ''}
    
            polarization: {}
    
    source_2:
    
      point source:
    
        position:
    
          l: {value: 11.25, min_value: 0.0, max_value: 360.0, delta: 1.125, free: false,
    
            unit: ''}
    
          b: {value: -22.5, min_value: -90.0, max_value: 90.0, delta: -2.25, free: false,
    
            unit: ''}
    
          equinox: J2000
    
        spectrum:
    
          synchrotron:
    
            shape:
    
              powerlaw:
    
                logK: {value: -0.3, min_value: -40, max_value: 40, delta: 0.0, free: true,
    
                  unit: dex(1 / (cm2 keV s))}
    
                piv: {value: 1, min_value: null, max_value: null, delta: 0.1, free: false,
    
                  unit: keV}
    
                index: {value: -0.75, min_value: -10, max_value: 10, delta: -0.2, free: true,
    
                  unit: ''}
    
            polarization: {}
    
          Inverse_Compton:
    
            shape:
    
              powerlaw:
    
                logK: {value: -0.3, min_value: -40, max_value: 40, delta: 0.0, free: true,
    
                  unit: dex(1 / (cm2 keV s))}
    
                piv: {value: 1, min_value: null, max_value: null, delta: 0.1, free: false,
    
                  unit: keV}
    
                index: {value: -1.7, min_value: -10, max_value: 10, delta: -0.2, free: true,
    
                  unit: ''}
    
            polarization: {}
    
    


Now suppose that you want to load back a file you created in a previous
session. You can do it with:

.. code:: python

    # Re-load the model from the file, just to show how it is done
    
    my_model_2 = load_model('my_model.yml')

.. code:: python

    # Explore the model we just loaded back
    
    my_model_2.display()



.. raw:: html

    Point sources: source_1,source_2<br><br>Extended sources: (none)<br><br>Parameters:<br><br>source_1<br><table id="table71521552">
    <thead><tr><th>name</th><th>value</th><th>min_value</th><th>max_value</th><th>delta</th><th>free</th><th>unit</th></tr></thead>
    <tr><td>position.ra</td><td>125.6</td><td>0</td><td>360</td><td>12.56</td><td>False</td><td></td></tr>
    <tr><td>position.dec</td><td>-75.3</td><td>-90</td><td>90</td><td>-7.53</td><td>False</td><td></td></tr>
    <tr><td>main.powerlaw.logK</td><td>-0.912573542964</td><td>-40</td><td>40</td><td>0.0</td><td>True</td><td>dex(1 / (cm2 keV s))</td></tr>
    <tr><td>main.powerlaw.piv</td><td>1.0</td><td>None</td><td>None</td><td>0.1</td><td>False</td><td>keV</td></tr>
    <tr><td>main.powerlaw.index</td><td>-2.1</td><td>-10</td><td>10</td><td>-0.2</td><td>True</td><td></td></tr>
    </table><br><br>source_2<br><table id="table71560592">
    <thead><tr><th>name</th><th>value</th><th>min_value</th><th>max_value</th><th>delta</th><th>free</th><th>unit</th></tr></thead>
    <tr><td>position.l</td><td>11.25</td><td>0</td><td>360</td><td>1.125</td><td>False</td><td></td></tr>
    <tr><td>position.b</td><td>-22.5</td><td>-90</td><td>90</td><td>-2.25</td><td>False</td><td></td></tr>
    <tr><td>synchrotron.powerlaw.logK</td><td>-0.3</td><td>-40</td><td>40</td><td>0.0</td><td>True</td><td>dex(1 / (cm2 keV s))</td></tr>
    <tr><td>synchrotron.powerlaw.piv</td><td>1.0</td><td>None</td><td>None</td><td>0.1</td><td>False</td><td>keV</td></tr>
    <tr><td>synchrotron.powerlaw.index</td><td>-0.75</td><td>-10</td><td>10</td><td>-0.2</td><td>True</td><td></td></tr>
    <tr><td>Inverse_Compton.powerlaw.logK</td><td>-0.3</td><td>-40</td><td>40</td><td>0.0</td><td>True</td><td>dex(1 / (cm2 keV s))</td></tr>
    <tr><td>Inverse_Compton.powerlaw.piv</td><td>1.0</td><td>None</td><td>None</td><td>0.1</td><td>False</td><td>keV</td></tr>
    <tr><td>Inverse_Compton.powerlaw.index</td><td>-1.7</td><td>-10</td><td>10</td><td>-0.2</td><td>True</td><td></td></tr>
    </table><br>


.. code:: python

    # Now evaluate and plot our models. You need matplotlib for this
    
    import matplotlib.pyplot as plt
    
    %matplotlib inline
    
    # Energies where we want to evaluate the model
    
    e = np.logspace(0,3,100)
    
    # Loop over the sources
    
    for src_name, src in my_model.point_sources.iteritems():
        
        # Loop over the components of each source
        
        for comp_name, component in src.components.iteritems():
            
            # Get the differential flux (in ph/cm2/s)
            
            flux = component.shape(e)
            
            # this can also be accomplished with:
            # flux = component.powerlaw(e)
            # but this requires to know the name of the
            # spectral shape which was used
            
            # Plot this component for this source
            
            plt.plot(e,e * e * flux,label="%s of %s" % (component.name, src.name))
    
    plt.legend(loc=0)
    plt.loglog()
    plt.xlabel("Energy")
    plt.ylabel(r"$\nu$F$_{\nu}$")




.. parsed-literal::

    <matplotlib.text.Text at 0x3d9cd90>




.. image:: create_model_files/create_model_23_1.png


