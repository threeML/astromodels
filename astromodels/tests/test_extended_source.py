import astropy.units as u
import numpy as np
import pytest

from astromodels.core.spectral_component import SpectralComponent
from astromodels.functions.functions import Powerlaw
from astromodels.functions.functions_2d import Gaussian_on_sphere
from astromodels.sources.extended_source import ExtendedSource


from astromodels.functions.priors import *
from astromodels.functions.function import _known_functions

__author__ = 'henrikef'


def test_constructor():

    # RA, Dec and L,B of the same point in the sky

    ra, dec = (125.6, -75.3)
    l, b = (288.44190139183564, -20.717313145391525)

    # This should throw as we are using Powerlaw instead of Powerlaw()
    with pytest.raises(TypeError):

        _ = ExtendedSource("my_source", Gaussian_on_sphere, Powerlaw)

    # Init with RA, Dec

    shape = Gaussian_on_sphere()
    source1 = ExtendedSource('my_source', shape, Powerlaw())
    shape.lon0=ra*u.degree
    shape.lat0=dec*u.degree

    assert source1.shape.lon0 == ra
    assert source1.shape.lat0 == dec

    # Verify that the position is fixed by default
    assert source1.shape.lon0.fix
    assert source1.shape.lon0.fix

test_constructor()


def test_call():

    # Multi-component

    po1 = Powerlaw()
    po2 = Powerlaw()

    c1 = SpectralComponent("component1", po1)
    c2 = SpectralComponent("component2", po2)

    point_source = PointSource("test_source", 125.4, -22.3, components=[c1, c2])

    assert np.all(point_source.spectrum.component1([1, 2, 3]) == po1([1, 2, 3]))
    assert np.all(point_source.spectrum.component2([1, 2, 3]) == po2([1, 2, 3]))

    one = point_source.spectrum.component1([1, 2, 3])
    two = point_source.spectrum.component2([1, 2, 3])

    assert np.all( np.abs(one + two - point_source([1,2,3])) == 0 )


def test_call_with_units():

    po = Powerlaw()

    result = po(1.0)

    assert result.ndim == 0

    with pytest.raises(AssertionError):

        # This raises because the units of the function have not been set up

        _ = po(1.0 * u.keV)

    # Use the function as a spectrum
    ps = PointSource("test",0,0,po)

    result = po(1.0 * u.keV)

    assert isinstance(result, u.Quantity)

    result = po(np.array([1,2,3])* u.keV)

    assert isinstance(result, u.Quantity)

    # Now test all the functions
    def test_one(class_type):

        instance = class_type()

        if not instance.is_prior:

            # if we have fixed x_units then we will use those
            # in the test

            if instance.has_fixed_units():

                x_unit_to_use = instance.fixed_units[0]

            else:

                x_unit_to_use = u.keV



            # Use the function as a spectrum
            ps = PointSource("test", 0, 0, instance)


            result = ps(1.0 * x_unit_to_use)


            assert isinstance(result, u.Quantity)

            result = ps(np.array([1, 2, 3]) * x_unit_to_use)

            assert isinstance(result, u.Quantity)

            result = ps(1.0)

            assert isinstance(result, float)

        else:

            print ('Skipping prior function')



    for key in _known_functions:

        this_function = _known_functions[key]

        # Test only the power law of XSpec, which is the only one we know we can test at 1 keV

        if key.find("XS")==0 and key != "XS_powerlaw":

            # An XSpec model. Test it only if it's a power law (the others might need other parameters during
            # initialization)

            continue

        if key.find("TemplateModel")==0:

            # The TemplateModel function has its own test

            continue

        if this_function._n_dim == 1:

            print("testing %s ..." % key)

            test_one(_known_functions[key])


def test_call_with_composite_function_with_units():

    def one_test(spectrum):

        print("Testing %s" % spectrum.expression)

        # # if we have fixed x_units then we will use those
        # # in the test
        #
        # if spectrum.expression.has_fixed_units():
        #
        #     x_unit_to_use, y_unit_to_use = spectrum.expression.fixed_units[0]
        #
        # else:

        x_unit_to_use = u.keV

        pts = PointSource("test", ra=0, dec=0, spectral_shape=spectrum)

        res = pts([100, 200] * x_unit_to_use)

        # This will fail if the units are wrong
        res.to(1 / (u.keV * u.cm**2 * u.s))

    # Test a simple composition

    spectrum = Powerlaw() * Exponential_cutoff()

    one_test(spectrum)

    spectrum = Band() + Blackbody()

    one_test(spectrum)

    # Test a more complicate composition

    spectrum = Powerlaw() * Exponential_cutoff() + Blackbody()

    one_test(spectrum)

    spectrum = Powerlaw() * Exponential_cutoff() * Exponential_cutoff() + Blackbody()

    one_test(spectrum)

    if has_xspec:

        spectrum = XS_phabs() * Powerlaw()

        one_test(spectrum)

        spectrum = XS_phabs() * XS_powerlaw()

        one_test(spectrum)

        spectrum = XS_phabs() * XS_powerlaw() * XS_phabs()

        one_test(spectrum)

        spectrum = XS_phabs() * XS_powerlaw() * XS_phabs() + Blackbody()

        one_test(spectrum)

        spectrum = XS_phabs() * XS_powerlaw() * XS_phabs() + XS_powerlaw()

        one_test(spectrum)


def test_free_param():

    spectrum = Log_parabola()
    source = PointSource("test_source", ra=123.4, dec=56.7, spectral_shape=spectrum)

    parameters = [spectrum.alpha, spectrum.beta, spectrum.piv, spectrum.K, source.position.ra, source.position.dec]

    for param in parameters:
        param.free = False

    assert len(source.free_parameters) == 0

    for i, param in enumerate(parameters):
        param.free = True
        assert len(source.free_parameters) == i+1

