
Functions tutorial
==================

In astromodels functions can be used as spectral shapes for sources, or
to describe time-dependence, phase-dependence, or links among
parameters.

To get the list of available functions just do:

.. code:: python

    from astromodels import *
    
    list_functions()




.. raw:: html

    <table id="table48641168">
    <thead><tr><th>name</th><th>Description</th></tr></thead>
    <tr><td>band</td><td>The Band model from Band et al. 1993, implemented however in a way which reduces the covariances between the parameters (Calderone et al., MNRAS, 448, 403C, 2015)</td></tr>
    <tr><td>bias</td><td>Return x plus a bias</td></tr>
    <tr><td>gaussian</td><td>A Gaussian function</td></tr>
    <tr><td>identity</td><td>Return x</td></tr>
    <tr><td>line</td><td>A linear function</td></tr>
    <tr><td>log_parabola</td><td>A log-parabolic function</td></tr>
    <tr><td>log_uniform_prior</td><td>A function which is 1/x on the interval lower_bound - upper_bound and 0 outside the interval. The extremes of the interval are NOT counted as part of the interval. Lower_bound must be strictly positive.</td></tr>
    <tr><td>powerlaw</td><td>A simple power-law with normalization expressed as a logarithm</td></tr>
    <tr><td>sin</td><td>A sinusodial function</td></tr>
    <tr><td>synchrotron</td><td>Synchrotron spectrum from an input particle distribution, using Naima (naima.readthedocs.org)</td></tr>
    <tr><td>uniform_prior</td><td>A function which is constant on the interval lower_bound - upper_bound and 0 outside the interval. The extremes of the interval are counted as part of the interval.</td></tr>
    </table>



If you need more info about a function, you can obtain it by using:

.. code:: python

    gaussian.info()



.. raw:: html

    <ul>
    
    <li>description: A Gaussian function</li>
    
    <li>formula: $ 10^{logK} \frac{1}{\sigma \sqrt{2 \pi}}\exp{\frac{(x-\mu)^2}{2~\sigma^2}} $</li>
    
    <li>default parameters: 
    <ul>
    
    <li>logK: 
    <ul>
    
    <li>value: 0.0</li>
    
    <li>min_value: -40</li>
    
    <li>max_value: 40</li>
    
    <li>unit: dex(1 / (cm2 keV s))</li>
    
    <li>delta: 0.1</li>
    
    <li>free: True</li>
    
    </ul>
    
    </li>
    
    <li>mu: 
    <ul>
    
    <li>value: 0.0</li>
    
    <li>min_value: None</li>
    
    <li>max_value: None</li>
    
    <li>unit: </li>
    
    <li>delta: 0.1</li>
    
    <li>free: True</li>
    
    </ul>
    
    </li>
    
    <li>sigma: 
    <ul>
    
    <li>value: 1.0</li>
    
    <li>min_value: None</li>
    
    <li>max_value: None</li>
    
    <li>unit: </li>
    
    <li>delta: 0.1</li>
    
    <li>free: True</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>



Note that you don't need to create an instance in order to call the
info() method.

Creating functions
------------------

Functions can be created in two different ways. We can create an
instance with the default values for the parameters like this:

.. code:: python

    powerlaw_instance = powerlaw()

or we can specify on construction specific values for the parameters:

.. code:: python

    powerlaw_instance = powerlaw(logK=-2.0, index=-2.2)

If you don't remember the names of the parameters just call the .info()
method as in powerlaw.info() as demonstrated above.

Getting information about an instance
-------------------------------------

Using the .display() method we get a representation of the instance
which exploits the features of the environment we are using. If we are
running inside a IPython notebook, a rich representation with the
formula of the function will be displayed (if available). Otherwise, in
a normal terminal, the latex formula will not be rendered:

.. code:: python

    powerlaw_instance.display()



.. raw:: html

    <ul>
    
    <li>description: A simple power-law with normalization expressed as a logarithm</li>
    
    <li>formula: $ \frac{dN}{dx} = 10^{logK}~\frac{x}{piv}^{index} $</li>
    
    <li>parameters: 
    <ul>
    
    <li>logK: 
    <ul>
    
    <li>value: -2.0</li>
    
    <li>min_value: -40</li>
    
    <li>max_value: 40</li>
    
    <li>unit: dex(1 / (cm2 keV s))</li>
    
    <li>delta: 0.1</li>
    
    <li>free: True</li>
    
    </ul>
    
    </li>
    
    <li>piv: 
    <ul>
    
    <li>value: 1.0</li>
    
    <li>min_value: None</li>
    
    <li>max_value: None</li>
    
    <li>unit: keV</li>
    
    <li>delta: 0.1</li>
    
    <li>free: False</li>
    
    </ul>
    
    </li>
    
    <li>index: 
    <ul>
    
    <li>value: -2.2</li>
    
    <li>min_value: -10</li>
    
    <li>max_value: 10</li>
    
    <li>unit: </li>
    
    <li>delta: 0.2</li>
    
    <li>free: True</li>
    
    </ul>
    
    </li>
    
    </ul>
    
    </li>
    
    </ul>



It is also possible to get the text-only representation by simply
printing the object like this:

.. code:: python

    print(powerlaw_instance)


.. parsed-literal::

      * description: A simple power-law with normalization expressed as a logarithm
      * formula: $ \frac{dN}{dx} = 10^{logK}~\frac{x}{piv}^{index} $
      * parameters: 
        * logK: 
          * value: -2.0
          * min_value: -40
          * max_value: 40
          * unit: dex(1 / (cm2 keV s))
          * delta: 0.1
          * free: True
        * piv: 
          * value: 1.0
          * min_value: None
          * max_value: None
          * unit: keV
          * delta: 0.1
          * free: False
        * index: 
          * value: -2.2
          * min_value: -10
          * max_value: 10
          * unit: 
          * delta: 0.2
          * free: True
    
    


**NOTE**: the .display() method of an instance displays the *current*
values of the parameters, while the .info() method demonstrated above
(for which you don't need an instance) displays the *default* values of
the parameters.

Modifying parameters
--------------------

Modifying a parameter of a function is easy:

.. code:: python

    # Modify current value
    
    powerlaw_instance.logK = 1.2
    
    # Modify minimum 
    powerlaw_instance.logK.min_value = -10
    
    # Modify maximum
    powerlaw_instance.logK.max_value = 15
    
    # We can also modify minimum and maximum at the same time
    powerlaw_instance.logK.set_bounds(-10, 15)
    
    # Modifying the delta for the parameter 
    # (which can be used by downstream software for fitting, for example)
    powerlaw_instance.logK.delta = 0.25
    
    # Fix the parameter
    powerlaw_instance.logK.fix = True
    
    # or equivalently
    powerlaw_instance.logK.free = False
    
    # Free it again
    powerlaw_instance.logK.fix = False
    
    # or equivalently
    powerlaw_instance.logK.free = True
    
    # We can verify what we just did by printing again the whole function as shown above, 
    # or simply printing the parameter:
    powerlaw_instance.logK.display()



.. raw:: html

    Parameter logK = 1.2
    (min_value = -10, max_value = 15, delta = 0.25, free = True)


Using physical units
~~~~~~~~~~~~~~~~~~~~

Astromodels uses the facility defined in astropy.units to make easier to
convert between units during interactive analysis, when assigning to
parameters. To use this feature you have to use the .set() method of the
parameter class:

.. code:: python

    import astropy.units as u
    
    # Express the differential flux at the pivot energy in 1 / (MeV cm2 s)
    
    powerlaw_instance.logK.set(122.3 / (u.MeV * u.cm * u.cm * u.s))
    
    # Express the differential flux at the pivot energy in 1 / (GeV m2 s)
    powerlaw_instance.logK.set(122.3 / (u.GeV * u.m * u.m * u.s))


.. parsed-literal::

    /home/giacomov/software/canopy-env/lib/python2.7/site-packages/astromodels-0.1-py2.7.egg/astromodels/parameter.py:363: WarningUnitsAreSlow: Using units is convenient but slow. Do not use them during computing-intensive work.
      WarningUnitsAreSlow)


The first time you use this feature in a interactive session or in a
script you will see a warning telling you that using this facility is
convenient but slow. Indeed, it is discouraged to use the .set() method
in a computing-intensive situation. If you are running IPython, you can
verify how much .set() is slower than just using a plain number with
these instructions:

.. code:: python

    # NOTE: These requires IPython
    
    %timeit powerlaw_instance.logK = -1.3
    %timeit powerlaw_instance.logK.set(122.3 / (u.GeV * u.m * u.m * u.s))


.. parsed-literal::

    The slowest run took 4.07 times longer than the fastest. This could mean that an intermediate result is being cached 
    100000 loops, best of 3: 5.38 µs per loop
    1000 loops, best of 3: 786 µs per loop


As you can see using astropy.units requires almost 100x more than using
a plain assignment. In an interactive analysis you are unlikely to
notice the difference, but if you use .set() in a loop or during a fit
this slow-down will add up an become very noticeable. Note that this is
a feature of astropy.units, not of astromodels.

Composing functions
-------------------

We can create arbitrary complex functions by combining "primitive"
functions using the normal math operators:

.. code:: python

    composite = gaussian() + powerlaw()
    
    # Instead of the usual .display(), which would print all the many parameters,
    # let's print just the description of the new composite functions:
    print(composite.description)


.. parsed-literal::

    (gaussian{1} + powerlaw{2})


These expressions can be as complicated as needed. For example:

.. code:: python

    crazy_function = 3 * sin() + powerlaw()**2 * (5+gaussian()) / 3.0
    
    print(crazy_function.description)


.. parsed-literal::

    ((sin{1} * 3) + (((powerlaw{2} ** 2) * (gaussian{3} + 5)) / 3.0))


The numbers between {} enumerate the unique functions which constitute a
composite function. This is useful because composite functions can be
created starting from pre-existing instances of functions, in which case
the same instance can be used more than once. For example:

.. code:: python

    a_powerlaw = powerlaw()
    a_sin = sin()
    
    another_composite = 2 * a_powerlaw + (3 + a_powerlaw) * a_sin
    
    print(another_composite.description)


.. parsed-literal::

    ((powerlaw{1} * 2) + ((powerlaw{1} + 3) * sin{2}))


In this case the same instance of a power law has been used twice.
Changing the value of the parameters for "a\_powerlaw" will affect also
the second part of the expression. Instead, by doing this:

.. code:: python

    another_composite2 = 2 * powerlaw() + (3 + powerlaw()) * sin()
    
    print(another_composite2.description)


.. parsed-literal::

    ((powerlaw{1} * 2) + ((powerlaw{2} + 3) * sin{3}))


we will end up with two independent sets of parameters for the two power
laws. The difference can be seen immediately from the number of
parameters of the two composite functions:

.. code:: python

    print(len(another_composite.parameters)) # 6 parameters
    print(len(another_composite2.parameters)) # 9 parameters


.. parsed-literal::

    6
    9


