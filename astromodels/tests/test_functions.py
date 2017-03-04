import pytest

import astropy.units as u
import numpy as np
import pickle

from astromodels.functions.function import FunctionMeta, Function1D, Function2D, FunctionDefinitionError, \
    UnknownParameter, DesignViolation, get_function, get_function_class, UnknownFunction, list_functions
from astromodels.functions.functions import Powerlaw, Line
from astromodels.functions.functions_2D import Gaussian_on_sphere
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

    with pytest.raises(DesignViolation):

        _ = f.get_boundaries()


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

    # Now test that an error is raised if units are not intelligible
    with pytest.raises(TypeError):

        _ = my_function.set_units("non_existent","non_existent")


def test_function_composition():

    Test_function = get_a_function_class()

    line = Test_function()
    powerlaw = Powerlaw()

    composite = powerlaw + line

    composite.set_units(u.keV, 1.0 / (u.keV * u.cm**2 * u.s))

    for x in ([1,2,3,4],[1,2,3,4] * u.keV, 1.0, np.array([1.0, 2.0, 3.0, 4.0])):

        assert np.all(composite(x) == line(x) + powerlaw(x))

    # Test -
    po = Powerlaw()
    li = Line()
    composite = po - li

    assert composite(1.0) == (po(1.0) - li(1.0))

    # test *
    composite = po * li

    assert composite(2.25) == po(2.25) * li(2.25)

    # test /
    composite = po / li

    assert composite(2.25) == po(2.25) / li(2.25)

    # test .of
    composite = po.of(li)

    assert composite(2.25) == po(li(2.25))

    # test power
    composite = po**li

    assert composite(2.25)  == po(2.25)**li(2.25)

    # test negation
    neg_po = -po

    assert neg_po(2.25) == -po(2.25)

    # test abs
    new_li = Line()
    new_li.b = -10.0

    abs_new_li = abs(new_li)

    assert new_li(1.0) < 0
    assert abs_new_li(1.0) == abs(new_li(1.0))

    # test rpower
    composite = 2.0**new_li

    assert composite(2.25) == 2.0**(new_li(2.25))

    # test multiplication by a number
    composite = 2.0 * po

    assert composite(2.25) == 2.0 * po(2.25)

    # Number divided by
    composite = 1.0 / li

    assert composite(2.25) == 1.0 / li(2.25)

    # Composite of composite
    composite = po*li + po - li + 2*po / li

    assert composite(2.25) == po(2.25) * li(2.25) + po(2.25) - li(2.25) + 2*po(2.25) / li(2.25)

    print(composite)


def test_duplicate():

    instance = Powerlaw()
    instance.index = -2.25
    instance.K = 0.5

    # Duplicate it

    duplicate = instance.duplicate()

    # Check that we have the same results

    assert duplicate(2.25) == instance(2.25)

    # Check that the parameters are not linked anymore
    instance.index = -1.12

    assert instance.index.value != duplicate.index.value

    print(instance)
    print(duplicate)


def test_pickling_unpickling():

    # 1d function
    po = Powerlaw()

    po.K = 5.35

    new_po = pickle.loads(pickle.dumps(po))

    assert new_po.K.value == po.K.value

    # 2d function
    gs = Gaussian_on_sphere()

    _ = pickle.loads(pickle.dumps(gs))

    # 3d function
    c = Continuous_injection_diffusion()

    _ = pickle.loads(pickle.dumps(c))

    # composite function
    po2 = Powerlaw()
    li = Line()
    composite = po2*li + po2 - li + 2*po2 / li  # type: Function1D

    # Change some parameter
    composite.K_1 = 3.2
    composite.a_2 = 1.56

    dump = pickle.dumps(composite)

    new_composite = pickle.loads(dump)

    assert new_composite.K_1.value == composite.K_1.value
    assert new_composite.a_2.value == composite.a_2.value


def test_get_function():

    po = get_function("Powerlaw")

    _ = po(1.0)

    with pytest.raises(UnknownFunction):

        _ = get_function("not_existant")


def test_get_function_class():

    po_class = get_function_class("Powerlaw")

    assert po_class == Powerlaw

    with pytest.raises(UnknownFunction):

        _ = get_function_class("not_existant")


def test_list_functions():

    print list_functions()


def test_function2D():

    c = Gaussian_on_sphere()

    _ = c(1, 1)

    a = np.array([1.0, 2.0])

    _ = c(a, a)

    c.set_units(u.deg, u.deg, 1.0 / u.deg**2)

    _ = c(1 * u.deg, 1.0 * u.deg)

    _ = c(a * u.deg, a * u.deg)

    print c.x_unit
    print c.y_unit
    print c.z_unit

    with pytest.raises(TypeError):

        c.set_units("not existent", u.deg, u.keV)


def test_function3D():

    c = Continuous_injection_diffusion()

    _ = c(1, 1, 1)

    a = np.array([1.0, 2.0])

    _ = c(a, a, a)

    c.set_units(u.deg, u.deg, u.keV, 1.0 / u.deg**2)

    _ = c(1 * u.deg, 1.0 * u.deg, 1.0 * u.keV)

    _ = c(a * u.deg, a * u.deg, a * u.keV)

    print c.x_unit
    print c.y_unit
    print c.z_unit
    print c.w_unit

    with pytest.raises(TypeError):

        c.set_units("not existent", u.deg, u.keV, 1.0 / (u.keV * u.s * u.deg**2 * u.cm**2))