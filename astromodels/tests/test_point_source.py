import astropy.units as u
import numpy as np
import pytest

from astromodels.core.spectral_component import SpectralComponent
from astromodels.functions.functions import Powerlaw, Exponential_cutoff, Blackbody, Band
from astromodels.sources.point_source import PointSource

try:

    from astromodels.xspec import XS_phabs, XS_powerlaw

except:

    has_xspec = False

else:

    has_xspec = True


from astromodels.functions.function import _known_functions

__author__ = 'giacomov'


def test_constructor():

    # RA, Dec and L,B of the same point in the sky

    ra, dec = (125.6, -75.3)
    l, b = (288.44190139183564, -20.717313145391525)

    # This should throw as we are using Powerlaw instead of Powerlaw()
    with pytest.raises(TypeError):

        _ = PointSource("my_source", ra, dec, Powerlaw)

    # Init with RA, Dec

    point_source1 = PointSource('my_source',ra, dec, Powerlaw())

    assert point_source1.position.get_ra() == ra
    assert point_source1.position.get_dec() == dec

    assert abs(point_source1.position.get_l() - l) < 1e-7
    assert abs(point_source1.position.get_b() - b) < 1e-7

    assert point_source1.position.ra.value == ra
    assert point_source1.position.dec.value == dec

    # Verify that the position is fixed by default
    assert point_source1.position.ra.fix
    assert point_source1.position.dec.fix

    # Init with l,b

    point_source2 = PointSource('my_source', l=l, b=b, spectral_shape=Powerlaw())

    assert point_source2.position.get_l() == l
    assert point_source2.position.get_b() == b

    assert abs(point_source2.position.get_ra() - ra) < 1e-7
    assert abs(point_source2.position.get_dec() - dec) < 1e-7

    assert point_source2.position.l.value == l
    assert point_source2.position.b.value == b

    # Verify that the position is fixed by default
    assert point_source2.position.l.fix
    assert point_source2.position.b.fix

    # Multi-component

    po1 = Powerlaw()
    po2 = Powerlaw()

    c1 = SpectralComponent("component1", po1)
    c2 = SpectralComponent("component2", po2)

    point_source3 = PointSource("test_source", ra, dec, components=[c1, c2])

    assert np.all(point_source3.spectrum.component1([1,2,3]) == po1([1,2,3]))
    assert np.all(point_source3.spectrum.component2([1,2,3]) == po2([1,2,3]))

    with pytest.raises(AssertionError):

        # Illegal RA

        _ = PointSource("test",720.0, -15.0, components=[c1,c2])

    with pytest.raises(AssertionError):
        # Illegal Dec

        _ = PointSource("test", 120.0, 180.0, components=[c1, c2])

    with pytest.raises(AssertionError):
        # Illegal l

        _ = PointSource("test", l=-195, b=-15.0, components=[c1, c2])

    with pytest.raises(AssertionError):
        # Illegal b

        _ = PointSource("test", l=120.0, b=-180.0, components=[c1, c2])


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

        # Use the function as a spectrum
        ps = PointSource("test", 0, 0, instance)

        result = ps(1.0 * u.keV)

        assert isinstance(result, u.Quantity)

        result = ps(np.array([1, 2, 3]) * u.keV)

        assert isinstance(result, u.Quantity)

        result = ps(1.0)

        assert isinstance(result, float)

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

        pts = PointSource("test", ra=0, dec=0, spectral_shape=spectrum)

        res = pts([100, 200] * u.keV)

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
