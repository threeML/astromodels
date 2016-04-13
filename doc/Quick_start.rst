
Quick start
===========

In this quick tutorial we will learn how to create a simple model with
one point source, with a power law spectrum. You can of course use any
function instead of the power law. Use "list\_models()" to obtain a list
of available models.

Let's define the model:

.. code:: python

    from astromodels import *
    
    test_source = PointSource('test_source',ra=123.22, dec=-13.56, spectral_shape=powerlaw_flux())
    
    my_model = Model(test_source)


.. parsed-literal::

    /home/giacomov/software/canopy-env/lib/python2.7/site-packages/astropy/utils/introspection.py:154: UserWarning: Module errno was already imported from None, but /home/giacomov/software/canopy-env/lib/python2.7/site-packages is being added to sys.path
      from pkg_resources import parse_version


Now let's use it:

.. code:: python

    # Get and print the differential flux at 1 keV:
    
    differential_flux_at_1_keV = my_model.test_source(1.0)
    
    print("Differential flux @ 1 keV : %.3f photons / ( cm2 s keV)" % differential_flux_at_1_keV)


.. parsed-literal::

    Differential flux @ 1 keV : 1.010 photons / ( cm2 s keV)


.. code:: python

    # Evaluate the model on an array of 100 energies logarithmically distributed between 1 and 100 keV
    
    # Set up the energies
    
    energies = np.logspace(0,2,100)
    
    # Get the differential flux
    
    differential_flux = my_model.test_source(energies)

.. code:: python

    # Plot it with the help of matplotlib 
    # Matplotlib is not required by astromodels, although you need it
    # for this part of the example
    
    import matplotlib.pyplot as plt
    
    # Do not execute this if you are not using the IPython notebook
    
    ########
    
    %matplotlib inline
    
    ########
    
    plt.loglog(energies, differential_flux)
    
    plt.xlabel("Energy (keV)")
    plt.ylabel("Differential flux (ph./cm2/s/keV)")




.. parsed-literal::

    <matplotlib.text.Text at 0x5862910>




.. image:: Quick_start_files/Quick_start_7_1.png


