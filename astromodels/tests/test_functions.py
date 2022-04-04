from __future__ import division, print_function

import os
import pickle
from builtins import object

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits
from future.utils import with_metaclass

import astromodels
from astromodels import update_logging_level
from astromodels.core.property import SettingUnknownValue
from astromodels.functions import (Continuous_injection_diffusion,
                                   Gaussian_on_sphere, Line, Powerlaw,
                                   SpatialTemplate_2D)
from astromodels.functions import function as function_module
from astromodels.functions.function import (DesignViolation, Function1D,
                                            Function2D,
                                            FunctionDefinitionError,
                                            FunctionInstanceError,
                                            FunctionMeta, UnknownFunction,
                                            UnknownParameter, get_function,
                                            get_function_class, list_functions)
from astromodels.functions.functions_1D.absorption import phabs, tbabs, wabs
from astromodels.functions.functions_1D.functions import _ComplexTestFunction
from astromodels.utils.configuration import astromodels_config

update_logging_level("DEBUG")

__author__ = 'giacomov'


def get_a_function_class():

    # Try to create a function inheriting from Function with meta FunctionMeta
    class Test_function(with_metaclass(FunctionMeta, Function1D)):
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

        class Wrong_test_function1(with_metaclass(FunctionMeta, object)):
            pass

    with pytest.raises(AttributeError):

        # .evaluate is lacking, ._set_units is lacking

        class Wrong_test_function2(with_metaclass(FunctionMeta, Function1D)):
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

    with pytest.raises(AttributeError):
        # _set_units is lacking

        class Wrong_test_function3(with_metaclass(FunctionMeta, Function1D)):
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

            def evaluate(self, x, a, b):

                return a * x + b

    with pytest.raises(FunctionDefinitionError):
        # Signature of evaluate is wrong

        class Wrong_test_function4(with_metaclass(FunctionMeta, Function1D)):
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

            def _set_units(self, x_unit, y_unit):

                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self):

                return self.a * self.x + self.b

    with pytest.raises(FunctionDefinitionError):
        # Signature of evaluate is wrong

        class Wrong_test_function5(with_metaclass(FunctionMeta, Function1D)):
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

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x):
                return self.a * x + self.b

    with pytest.raises(FunctionDefinitionError):
        # Signature of evaluate is wrong

        class Wrong_test_function6(with_metaclass(FunctionMeta, Function1D)):
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

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a):
                return a * x + self.b

    with pytest.raises(FunctionDefinitionError):
        # Signature of evaluate does not match docstring

        class Wrong_test_function7(with_metaclass(FunctionMeta, Function1D)):
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

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(FunctionDefinitionError):
        # Definition of parameter b is not legal

        class Wrong_test_function8(with_metaclass(FunctionMeta, Function1D)):
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

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(FunctionDefinitionError):
        # Parameter c declared but not used

        class Wrong_test_function9(with_metaclass(FunctionMeta, Function1D)):
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

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(FunctionDefinitionError):
        # Parameter c used but not declared

        class Wrong_test_function10(with_metaclass(FunctionMeta, Function1D)):
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

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b, c):
                return a * x + b + c

    with pytest.raises(AssertionError):
        # Docstring lacking description

        class Wrong_test_function11(with_metaclass(FunctionMeta, Function1D)):
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

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(FunctionDefinitionError):
        # Parameter lacking description

        class Wrong_test_function12(with_metaclass(FunctionMeta, Function1D)):
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

            def _set_units(self, x_unit, y_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, x, a, b):
                return a * x + b

    with pytest.raises(AssertionError):
        # Parameters out of order in evaluate

        class Wrong_test_function13(with_metaclass(FunctionMeta, Function2D)):
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

            def _set_units(self, x_unit, y_unit, z_unit):
                self.a.unit = y_unit / x_unit
                self.b.unit = y_unit

            def evaluate(self, y, x, a, b, c):
                return a * x + b

    # A function with no latex formula (which is optional)

    class NoLatex_test_function11(with_metaclass(FunctionMeta, Function1D)):
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

    assert my_function(1.0) == 2

    my_function.a = 2.5

    assert my_function(10.0) == 26.0

    my_function.b = -1.0

    assert my_function(10.0) == 24.0

    # Now test with list and np.array

    my_function.a.value = 2.0
    my_function.b.value = 1.0

    assert np.all(my_function([1, 2, 3]) == np.array([3.0, 5.0, 7.0]))
    assert np.all(my_function(
        np.array([3, 4, 5])) == np.array([7.0, 9.0, 11.0]))


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

    assert np.all(my_function([1, 2, 3] * u.keV) ==
                  np.array([3.0, 5.0, 7.0]) * diff_flux)

    # Using one unit for each element will fail

    # (depending on the version of astropy, it might raise ValueError or TypeError)
    with pytest.raises((ValueError, TypeError)):

        _ = my_function([1 * u.keV, 2 * u.keV, 3 * u.keV])

    assert np.all(my_function(
        np.array([3, 4, 5]) * u.keV) == np.array([7.0, 9.0, 11.0]) * diff_flux)

    # Now test that an error is raised if units are not intelligible
    with pytest.raises(TypeError):

        _ = my_function.set_units("non_existent", "non_existent")


def test_function_composition():

    Test_function = get_a_function_class()

    line = Test_function()
    powerlaw = Powerlaw()

    composite = powerlaw + line

    composite.set_units(u.keV, 1.0 / (u.keV * u.cm**2 * u.s))

    for x in ([1, 2, 3, 4], [1, 2, 3, 4] * u.keV, 1.0, np.array([1.0, 2.0, 3.0, 4.0])):

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

    assert composite(2.25) == po(2.25)**li(2.25)

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

    assert composite(2.25) == po(2.25) * li(2.25) + \
        po(2.25) - li(2.25) + 2*po(2.25) / li(2.25)

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


    # now try with previously used function1d

    # composite function
    po3 = Powerlaw()
    
    composite2 = po3*li

    # Change some parameter
    composite2.K_1 = 3.2
    composite2.a_2 = 1.56

    dump2 = pickle.dumps(composite2)

    new_composite2 = pickle.loads(dump2)


    
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

    print(list_functions())


def test_function2D():

    c = Gaussian_on_sphere()

    f1 = c(1, 1)
    assert np.isclose(f1, 5.17276409, rtol=1e-10)

    a = np.array([1.0, 2.0])

    fa = c(a, a)
    assert np.isclose(fa, [5.17276409, 5.01992404], rtol=1e-10).all()

    c.set_units(u.deg, u.deg, 1.0 / u.deg**2)

    f1d = c(1 * u.deg, 1.0 * u.deg)
    assert np.isclose(f1d.value, 5.17276409 , rtol=1e-10)
    assert f1d.unit == u.deg**-2

    assert c.x_unit == u.deg
    assert c.y_unit == u.deg
    assert c.z_unit == u.deg**-2

    assert c.get_total_spatial_integral(1) == 1
    assert np.isclose(c.get_total_spatial_integral(
        [1, 1]),  [1, 1], rtol=1e-10).all()

    with pytest.raises(TypeError):

        c.set_units("not existent", u.deg, u.keV)


def test_function3D():

    c = Continuous_injection_diffusion()

    f1 = c(1, 1, 1)
    assert np.isclose(f1, 134.95394313247866, rtol=1e-10)

    a = np.array([1.0, 2.0])

    fa = c(a, a, a)
    assert np.isclose(fa,  [[134.95394313, 132.19796573], [
                      25.40751507, 27.321443]], rtol=1e-10).all()

    c.set_units(u.deg, u.deg, u.keV, 1.0 / u.deg**2)

    f1d = c(1 * u.deg, 1.0 * u.deg, 1.0 * u.keV)
    assert np.isclose(f1d.value, 134.95394313247866, rtol=1e-10)
    assert f1d.unit == u.deg**-2

    assert c.x_unit == u.deg
    assert c.y_unit == u.deg
    assert c.z_unit == u.keV
    assert c.w_unit == u.deg**-2

    assert c.get_total_spatial_integral(1) == 1
    assert np.isclose(c.get_total_spatial_integral(
        [1, 1]),  [1, 1], rtol=1e-10).all()

    with pytest.raises(TypeError):

        c.set_units("not existent", u.deg, u.keV, 1.0 /
                    (u.keV * u.s * u.deg**2 * u.cm**2))


def test_spatial_template_2D():

    # make the fits files with templates to test.
    cards = {
        "SIMPLE": "T",
        "BITPIX": -32,
        "NAXIS": 2,
        "NAXIS1": 360,
        "NAXIS2": 360,
        "DATE": '2018-06-15',
        "CUNIT1": 'deg',
        "CRVAL1":  83,
        "CRPIX1": 0,
        "CDELT1": -0.0166667,
        "CUNIT2": 'deg',
        "CRVAL2": -2.0,
        "CRPIX2": 0,
        "CDELT2": 0.0166667,
        "CTYPE1": 'GLON-CAR',
        "CTYPE2": 'GLAT-CAR'}

    data = np.zeros([400, 400])
    data[0:100, 0:100] = 1
    hdu = fits.PrimaryHDU(data=data, header=fits.Header(cards))
    hdu.writeto("test1.fits", overwrite=True)

    data[:, :] = 0
    data[200:300, 200:300] = 1
    hdu = fits.PrimaryHDU(data=data, header=fits.Header(cards))
    hdu.writeto("test2.fits", overwrite=True)

    # Now load template files and test their evaluation
    shape1 = SpatialTemplate_2D(fits_file="test1.fits")
    
    shape1.K = 1

    shape2 = SpatialTemplate_2D(fits_file="test2.fits")
    
    shape2.K = 1

    assert shape1.hash != shape2.hash

    assert np.all(shape1.evaluate(
        [312, 306], [41, 41], [1, 1], [40, 2], 0) == [1., 0.])
    assert np.all(shape2.evaluate(
        [312, 306], [41, 41], [1, 1], [40, 2], 0) == [0., 1.])
    assert np.all(shape1.evaluate(
        [312, 306], [41, 41], [1, 10], [40, 2], 0) == [1., 0.])
    assert np.all(shape2.evaluate(
        [312, 306], [41, 41], [1, 10], [40, 2], 0) == [0., 10.])

    shape1.K = 1
    shape2.K = 1
    assert np.all(shape1([312, 306], [41, 41], 0) == [1., 0.])
    assert np.all(shape2([312, 306], [41, 41], 0) == [0., 1.])

    shape1.K = 1
    shape2.K = 10
    assert np.all(shape1([312, 306], [41, 41], 0) == [1., 0.])
    assert np.all(shape2([312, 306], [41, 41], 0) == [0., 10.])
    
    
    os.remove("test1.fits")
    os.remove("test2.fits")

def test_linking_external_functions():

    p = Powerlaw()
    p2 = Powerlaw()


    # nothing there yet
    assert not p.external_functions
    
    p.link_external_function(p2, "p2")

    # should be there
    
    assert "p2" in p.external_functions

    with pytest.raises(RuntimeError):

        p.link_external_function(p2, "p2")


    p.unlink_external_function("p2")
    
    # nothing there now
    assert not p.external_functions

    with pytest.raises(RuntimeError):

        p.link_external_function("dummy", "p2")


    with pytest.raises(RuntimeError):

        p.unlink_external_function("p2")

    p.link_external_function(p2, "p2")

    p.unlink_all_external_functions()

    assert not p.external_functions

    p.link_external_function(p2, "p2")
    

    assert p2 == p.external_functions["p2"]

    data = p.to_dict(minimal=False)

    assert "external_functions" in data


    p3 = Powerlaw()

    p4 = p3 + p

    
    data = p4.to_dict(minimal=False)

    assert "external_functions" in data


    assert data["external_functions"][1]["p2"] == p2.path


def test_function_properties():
    
    with pytest.raises(FunctionInstanceError):

        c = _ComplexTestFunction()

    c = _ComplexTestFunction(file_name="lost.txt", dummy="test")
    

    with pytest.raises(SettingUnknownValue):

        c = _ComplexTestFunction(file_name="f.txt", dummy="wrong")


def test_abs_model():

    for i, m in enumerate([astromodels.TbAbs, astromodels.WAbs, astromodels.PhAbs]):

        instance = m()

        instance.info()

        instance.abundance_table_info
        
        if i != 1:
        
            instance.abundance_table = "AG89"

        if i ==0:

            assert tbabs._current_table == "AG89"

        if i == 2:

            assert phabs._current_table == "AG89"



            
            
    
