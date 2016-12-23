import pytest

import astropy.units as u
import numpy as np

from astromodels.functions.function import FunctionMeta, Function1D, Function2D, FunctionDefinitionError, \
    UnknownParameter, Function3D
from astromodels.functions.functions import Powerlaw
from astromodels.functions.functions_3D import Continuous_injection_diffusion
from astromodels.functions import function as function_module

__author__ = 'giacomov'


def get_a_function_class():

    # Try to create a function inheriting from Function with meta FunctionMeta
    class Test_function(Function1D):
        r"""
        description :

            A test function

        latex : $ a * x + b $

        parameters :

            a :

                desc : linear coefficient
                initial value : 1

            b :

                desc : intercept
                initial value : 1

        """

        __metaclass__ = FunctionMeta

        def _set_units(self, x_unit, y_unit):

            # a has units of y_unit / x_unit, so that a*x has units of y_unit
            self.a.unit = y_unit / x_unit

            # b has units of y
            self.b.unit = y_unit

        def evaluate(self, x, a, b):

            return a * x + b

    return Test_function


def test_function_meta():

    with pytest.raises(AttributeError):

        # .evaluate is lacking, ._set_units is lacking, docstring is lacking

        class Wrong_test_function1():
            __metaclass__ = FunctionMeta

    with pytest.raises(AttributeError):

        # .evaluate is lacking, ._set_units is lacking

        class Wrong_test_function2(Function1D):
            r"""
            description :

                A test function

            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                b :

                    desc : intercept
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

    with pytest.raises(AttributeError):
        # _set_units is lacking

        class Wrong_test_function3(Function1D):
            r"""
            description :

                A test function

            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                b :

                    desc : intercept
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def evaluate(self, x, a, b):

                return a * x + b

    with pytest.raises(AssertionError):
        # Signature of evaluate is wrong

        class Wrong_test_function4(Function1D):
            r"""
            description :

                A test function

            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                b :

                    desc : intercept
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit):

                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self):

                return self.a * self.x + self.b

    with pytest.raises(FunctionDefinitionError):
        # Signature of evaluate is wrong

        class Wrong_test_function5(Function1D):
            r"""
            description :

                A test function

            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                b :

                    desc : intercept
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x):
                return self.a * x + self.b

    with pytest.raises(FunctionDefinitionError):
        # Signature of evaluate is wrong

        class Wrong_test_function6(Function1D):
            r"""
            description :

                A test function

            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                b :

                    desc : intercept
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a):
                return a * x + self.b

    with pytest.raises(FunctionDefinitionError):
        # Signature of evaluate does not match docstring

        class Wrong_test_function7(Function1D):
            r"""
            description :

                A test function

            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                c :

                    desc : intercept
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(FunctionDefinitionError):
        # Definition of parameter b is not legal

        class Wrong_test_function8(Function1D):
            r"""
            description :

                A test function

            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                b :

                    desc : intercept

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(FunctionDefinitionError):
        # Parameter c declared but not used

        class Wrong_test_function9(Function1D):
            r"""
            description :

                A test function

            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                b :

                    desc : intercept
                    initial value : 1

                c :

                    desc : dumb
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(FunctionDefinitionError):
        # Parameter c used but not declared

        class Wrong_test_function10(Function1D):
            r"""
            description :

                A test function

            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                b :

                    desc : intercept
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b, c):
                return a * x + b + c


    with pytest.raises(AssertionError):
        # Docstring lacking description

        class Wrong_test_function11(Function1D):
            r"""
            latex : $ a * x + b $

            parameters :

                a :

                    desc : linear coefficient
                    initial value : 1

                b :

                    desc : intercept
                    initial value : 1

                c :

                    desc : dumb
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(FunctionDefinitionError):
        # Parameter lacking description

        class Wrong_test_function12(Function1D):
            r"""

            description: useless

            latex : $ a * x + b $

            parameters :

                a :

                    initial value : 1

                b :

                    desc : intercept
                    initial value : 1

                c :

                    desc : dumb
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(AssertionError):
        # Parameters out of order in evaluate

        class Wrong_test_function13(Function2D):
            r"""

            description: useless

            latex : $ a * x + b $

            parameters :

                a :

                    desc : blah
                    initial value : 1


                b :

                    desc : intercept
                    initial value : 1

                c :

                    desc : dumb
                    initial value : 1

            """

            __metaclass__ = FunctionMeta

            def _set_units(self, x_unit, y_unit, z_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, y, x, a, b, c):
                return a * x + b

    # A function with no latex formula (which is optional)

    class NoLatex_test_function11(Function1D):
        r"""

        description:

            A function without latex (should be fine)

        parameters :

            a :

                desc : linear coefficient
                initial value : 1

            b :

                desc : intercept
                initial value : 1

        """

        __metaclass__ = FunctionMeta

        def _set_units(self, x_unit, y_unit):
            self.a.unit = y_unit / x_unit
            self.b.unit = y_unit

        def evaluate(self, x, a, b):
            return a * x + b


def test_function_constructor():

    Test_function = get_a_function_class()

    # Instance with default parameters
    my_function = Test_function()

    assert my_function.a.value == 1.0
    assert my_function.b.value == 1.0

    # Instance with explicit parameters' values
    my_function = Test_function(b=3.2, a=-2.5)

    assert my_function.a.value == -2.5
    assert my_function.b.value == 3.2

    Test_function.info()

    function_module.has_ipython = False

    Test_function.info()

    print(my_function.free_parameters)

    with pytest.raises(UnknownParameter):

        f = Test_function(d=3.5)

    f = Test_function()

    print(f.description)

    print(f.latex)

    assert f.fixed_units is None

    assert f.has_fixed_units() == False


def test_function3D():

    c = Continuous_injection_diffusion()

    _ = c(1, 1, 1)

def test_function_values():

    # Get a function class, and test the various exceptions

    Test_function = get_a_function_class()

    # Try to instance it
    my_function = Test_function()

    # Test basic functionalities

    assert my_function(1.0)==2

    my_function.a = 2.5

    assert my_function(10.0) == 26.0

    my_function.b = -1.0

    assert my_function(10.0) == 24.0

    # Now test with list and np.array

    my_function.a.value = 2.0
    my_function.b.value = 1.0

    assert np.all(my_function([1,2,3]) == np.array([3.0, 5.0, 7.0]))
    assert np.all(my_function(np.array([3, 4, 5])) == np.array([7.0, 9.0, 11.0]))


def test_function_values_units():

    # Test units functionality
    # Get a function class, and test the various exceptions

    Test_function = get_a_function_class()

    # Try to instance it
    my_function = Test_function()

    diff_flux = 1.0 / (u.keV * u.cm**2 * u.s)

    my_function.set_units(u.keV, diff_flux)

    # Test basic functionalities

    assert my_function(1.0 * u.keV) == 2 * diff_flux

    my_function.a = 2.5 * diff_flux / u.keV

    assert my_function(10.0 * u.keV) == 26.0 * diff_flux

    my_function.b = -1.0 * diff_flux

    assert my_function(10.0 * u.keV) == 24.0 * diff_flux

    # Now test with list and np.array

    my_function.a.value = 2.0 * diff_flux / u.keV
    my_function.b.value = 1.0 * diff_flux

    assert np.all(my_function([1, 2, 3] * u.keV) == np.array([3.0, 5.0, 7.0]) * diff_flux)

    # Using one unit for each element will fail

    with pytest.raises(ValueError):

        _ = my_function([1 * u.keV, 2 * u.keV, 3 * u.keV])

    assert np.all(my_function(np.array([3, 4, 5]) * u.keV) == np.array([7.0, 9.0, 11.0]) * diff_flux)


def test_function_composition():

    Test_function = get_a_function_class()

    line = Test_function()
    powerlaw = Powerlaw()

    composite = powerlaw + line

    composite.set_units(u.keV, 1.0 / (u.keV * u.cm**2 * u.s))

    for x in ([1,2,3,4],[1,2,3,4] * u.keV, 1.0, np.array([1.0, 2.0, 3.0, 4.0])):

        assert np.all(composite(x) == line(x) + powerlaw(x))