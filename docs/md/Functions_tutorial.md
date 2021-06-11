---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
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

One of the most powerful aspects of astromodels is the ability to quickly build custom functions on the fly. The source code for a function can be pure python, FORTRAN linked via f2py, C++ linked via cython, etc. Anything that provides a python function can be used to fit data. 

To build a custom spectral 1D function in astromodels, we need to import a few things that will allow astromodels to recognize your model.

```python
from astromodels.functions.function import Function1D, FunctionMeta, ModelAssertionViolation
```

<!-- #region -->
```Function1D``` is the base class for 1D spectral models and ```FunctionMeta``` is a python meta type class that ensures all the needed parts of a model are in the class as well as making the class function as it should.


There are three basic parts to declaring a model:

* the docstring
* the units setter
* the evaluate function

Let's look at the simple case of the power law already define in astromodels.

<!-- #endregion -->

```python
class Powerlaw(Function1D, metaclass=FunctionMeta):
        r"""
        description :
            A  power-law
        latex : $ K~\frac{x}{piv}^{index} $
        parameters :
            K :
                desc : Normalization (differential flux at the pivot value)
                initial value : 1.0
                is_normalization : True
                transformation : log10
                min : 1e-30
                max : 1e3
                delta : 0.1
            piv :
                desc : Pivot value
                initial value : 1
                fix : yes
            index :
                desc : Photon index
                initial value : -2
                min : -10
                max : 10	
        """


        def _set_units(self, x_unit, y_unit):
            # The index is always dimensionless
            self.index.unit = astropy_units.dimensionless_unscaled

            # The pivot energy has always the same dimension as the x variable
            self.piv.unit = x_unit

            # The normalization has the same units as the y

            self.K.unit = y_unit


        def evaluate(self, x, K, piv, index):

            xx = np.divide(x, piv)

            return K * np.power(xx, index)


```

<!-- #region -->
### The docstring

We have used the docstring interface to provide a YAML description of the function. This sets up the important information used in the fitting process and record keeping. The docstring has three parts:

- description
    - The description is a text string that provides readable info about the model. Nothing fancy, but good descriptions help to inform the user.
- latex
    - If the model is analytic, a latex formula can be included
- parameters
    - For each parameter, a description and initial value must be included. Transformations for fitting, min/max values and fixing the parameter can also be described here.

Optionally, there can be an additional ```properties``` category that we will cover later.

    
Keep in mind that this is in YAML format.

### Set units

astromodels keeps track of units for you. However, a model must be set up to properly describe the units with astropy's unit system. Keep in mind that models are fit with a differential photon flux, 

$$\frac{d N_p}{dA dt dE}$$

so your units should reflect this convention. Therefore, proper normalizations should be taken into account.


### Evaluate
This is where the function is evaluated. The first argument **must be called x** and the parameter names and ordering must reflect what is in the docstring. Any number of operations can take place inside the evaluate call, but remember that the return must be in the form of a differential photon flux. For 2D and 3D functions, the functions have **y** and **z** for their first arguments as well. ** x, y, and z  are reserved names in functions**. 


A functions is defined in a python session. **If you save the results of a fit to an AnalysisResults file and try to load this file without loading this model, you will get a error** Thus, remember to import any local models you used for an analysis before trying to reload that analysis. 


## Custom functions in other langauges

What if your model is built from a C++ function and you want to fit that directly to the data? using Cython, pybind11, f2py, etc, you can wrap these models and call them easily.
<!-- #endregion -->

```python

def cpp_function_wrapper(a):
    # we could wrap a c++ function here
    # with cython, pybind11, etc
    
    return a

```

```python
cpp_function_wrapper(2.)
```

Now we will define a astromodels function that will handle both the unit and non-unit call.

```python
import astropy.units as astropy_units

class CppModel(Function1D,metaclass=FunctionMeta):
        r"""
        description :
            A spectral model wrapping a cython function
        latex : $$
        parameters :
            a :
                desc : Normalization (differential flux)
                initial value : 1.0
                is_normalization : True
                min : 1e-30
                max : 1e3
                delta : 0.1
        """

        def _set_units(self, x_unit, y_unit):

            # The normalization has the same units as the y

            self.a.unit = y_unit

        
        def evaluate(self, x, a):
            
            # check is the function is being called with units
            
            if isinstance(a, astropy_units.Quantity):
                
                # get the values
                a_ = a.value
                
                # save the unit
                unit_ = self.y_unit
                
            else:
                
                # we do not need to do anything here
                a_ = a
                
                # this will basically be ignored
                unit_ = 1.

            # call the cython function
            flux = cpp_function_wrapper(a_)

            # add back the unit if needed
            return flux * unit_
```

We can check the unit and non-unit call by making a point source and evaluating it

```python
cpp_spectrum = CppModel()

from astromodels import PointSource

point_source = PointSource('ps',0,0,spectral_shape=cpp_spectrum)

print(point_source(10.))
point_source(10. * astropy_units.keV)
```

## Advanced functions

We are not limited to functions that can only take numerical parameters as arguments. We can also link other astromodels function into a function to expand its abilities. 
### Properties

Let's create a function that uses text based switches to alter its functionality


```python
class SwitchFunction(Function1D,metaclass=FunctionMeta):
        r"""
        description :
            A demo function that can alter its state
        latex : $$
        parameters :
            a :
                desc : Normalization (differential flux)
                initial value : 1.0
                is_normalization : True
                min : 1e-30
                max : 1e3
                delta : 0.1
	    propertes:
		    switch:
			    initial value: powerlaw
				allowed values:
				    - powerlaw
					- cosine
				function: _say_hello
        """
        def _say_hello(self):
		    # called when we set the value of switch
		
 		    print(self.switch.value)


        def _set_units(self, x_unit, y_unit):

            # The normalization has the same units as the y

            self.a.unit = y_unit

        
        def evaluate(self, x, a):
            
			if self.switch.value == "powerlaw":
				
				return a * np.powerlaw(x,-2)
				
			elif self.switch.value == "cosine":
				
				return a * np.cos(x)
		
```

We have added a text parameter called switch and specified the allowed values that it can take on. We can of course allow it to take on any value. Additionally, we have specified a function to call whenever we change the value. This allows use to do things like read in a table from a dictionary or load a file, etc. Properties can be set in the constructor of a function. They behave just like parameters except that they do not participate in the function call. Thus, their state is saved whenever you serialize the model to disk. 

```python


f = SwitchFunction()

f.switch = 'cosine'

```


In the docstring, one can also specify ```defer: True``` which allows you to not set a value until instancing an object. This is useful if you have a model that reads in file at runtime, but the file name is not known until then. Check out the source code of astromodels to see how properties can be used to expand the functionality of your custom models. For example, the absorption models such as ```TbAbs``` take advantage of this to set their abundance tables.

### Linking functions

What if you need to call another astromodels function from inside your custom function? This is achieved by linking functions. For example, let's create a convolutional model that redshifts the energies of and model linked to it:

```python
class Redshifter(Function1D, metaclass=FunctionMeta):
    r"""
    description :
		a function that can redshift the energies of any 1D function

    latex: not available

    parameters :
        redshift :
            desc : the redshift
            initial value : 0
            min : 0

    """

    
    def _set_units(self, x_unit, y_unit):

        self.redshift.unit = astropy_units.dimensionless_unscaled 


    def set_linked_function(self, function):
	     # this is an optional helper to
		 # ease in the setting of the function

	

        if "func" in self._external_functions:
   
            # since we only want to link one function
			# we unlink if we have linked before 
			
            self.unlink_external_function("func")
 
        # this allows use to link in a function 
		# with an internal name 'func'
		
        self.link_external_function(function, "func")
        
        self._linked_function = function

    def get_linked_function(self):

	    # linked functions are stored in a dictionary
		# but you 

        return self._external_functions["func"]
        

    linked_function = property(
        get_linked_function,
        set_linked_function,
        doc="""Get/set linked function""",
    )


    def evaluate(self, x, redshift):


        # we can call the function here
        return self._linked_function(x * (1 + redshift))

```

With this function, whenever we set the linked_function property to another astromodels function, its call will return that function redshifted.

```python
p = Powerlaw()
rs = Redshifter()

rs.linked_function = p

print(p(10.))

print(rs(10.))


```


We have added a lot of syntax sugar to make it easier for users to handle the function, but every function in astromodels has the members ```f.link_external_function(func, 'internal_name')```  and ```f.unlink_external_function('internal_name')```. You can link as many functions as needed and they are accessed via an internal dictionary ```self._extranal_functions```. As long as all functions used are part of the model, all the linking is saved when a model is saved to disk allowing you to restore all the complexity you built. 
