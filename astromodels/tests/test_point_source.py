from __future__ import division, print_function

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pytest

from astromodels.core.spectral_component import SpectralComponent
from astromodels.functions import (Band, Blackbody, Exponential_cutoff,
                                   Log_parabola, Powerlaw)
from astromodels.functions.functions_1D.functions import _ComplexTestFunction

try:
    from astromodels.functions import PhAbs, TbAbs, WAbs

    has_abs_models = True

except:

    has_abs_models = False

    
from astromodels.core.model import Model
from astromodels.core.model_parser import clone_model, load_model
from astromodels.sources.particle_source import ParticleSource
from astromodels.sources.point_source import PointSource

try:

    from astromodels.xspec import XS_phabs, XS_powerlaw

except:

    has_xspec = False

else:

    has_xspec = True

try:

    from astromodels.functions import EBLattenuation

except:

    has_ebl = False

else:

    has_ebl = True

from astromodels.functions.function import _known_functions
from astromodels.functions.priors import *

__author__ = 'giacomov'



_multiplicative_models = ["PhAbs", "TbAbs", "WAbs", "EBLattenuation", "ZDust"]

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

        if class_type == _ComplexTestFunction:
        
            instance = class_type(file_name="test.txt")
            
        else:

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

            if instance.name in [ "Synchrotron", "_ComplexTestFunction" ]:

                # we should not do it this way

                p = Powerlaw()

                particleSource = ParticleSource("particles", p)

                

                instance.set_particle_distribution(p)


            # elif instance.name in ["PhAbs", "TbAbs"]:

            #     instance
                

            result = ps(1.0)

            assert isinstance(result, float)

            result = ps(1.0 * x_unit_to_use)

            assert isinstance(result, u.Quantity)

            result = ps(np.array([1, 2, 3]) * x_unit_to_use)

            assert isinstance(result, u.Quantity)
            
            if instance.name in [ "Synchrotron", "_ComplexTestFunction" ]:
              model = Model( particleSource, ps)
                  
                  
            else:
              model = Model( ps )
            
            new_model = clone_model( model )
            
            new_result =  new_model["test"](np.array([1, 2, 3]) * x_unit_to_use)
            
            assert np.all(new_result==result)

            model.save("__test.yml", overwrite=True)
            
            new_model = load_model("__test.yml")

            new_result =  new_model["test"](np.array([1, 2, 3]) * x_unit_to_use)
            
            assert np.all(new_result==result)

        else:

            print ('Skipping prior function')



    for key in _known_functions:

        this_function = _known_functions[key]

        # Test only the power law of XSpec, which is the only one we know we can test at 1 keV

        if key.find("XS")==0 and key != "XS_powerlaw" or (key in _multiplicative_models):

            # An XSpec model. Test it only if it's a power law (the others might need other parameters during
            # initialization)

            continue

        if key.find("TemplateModel")==0:

            # The TemplateModel function has its own test

            continue

        if key.find("Synchrotron")==0:

            # Naima Synchtron function should have its own test

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

        res = pts(np.array([100., 200.]) * x_unit_to_use)

        # This will fail if the units are wrong
        res.to(old_div(1, (u.keV * u.cm**2 * u.s)))

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

    # test the absorption models

    if has_abs_models:
    
        spectrum = PhAbs() * Powerlaw()


        one_test(spectrum)

        spectrum = TbAbs() * Powerlaw()

        one_test(spectrum)


        spectrum = WAbs() * Powerlaw()

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
        
    if has_ebl:
    
        spectrum = Powerlaw() * EBLattenuation()
        
        one_test(spectrum)


def test_free_param():

    spectrum = Log_parabola()
    source = PointSource("test_source", ra=123.4, dec=56.7, spectral_shape=spectrum)

    parameters = [spectrum.alpha, spectrum.beta, spectrum.piv, spectrum.K, source.position.ra, source.position.dec]

    for param in parameters:
        param.free = False

    assert not source.has_free_parameters
        
    assert len(source.free_parameters) == 0

    assert len(source.parameters) > len(source.free_parameters)
    
    for i, param in enumerate(parameters):
        param.free = True
        assert len(source.free_parameters) == i+1


    assert source.has_free_parameters


def test_local_deriv():

    p = Powerlaw(index=-2.)


    npt.assert_allclose(-2., p.local_spectral_index(np.logspace(1,3,10)))
        

        
