import astropy.units as u
import numpy as np
import pytest

from astromodels.core.spectral_component import SpectralComponent
from astromodels.functions.functions import Powerlaw, Log_parabola
from astromodels.functions.functions_2D import Gaussian_on_sphere, Disk_on_sphere
from astromodels.sources.extended_source import ExtendedSource


from astromodels.functions.priors import *
from astromodels.functions.function import _known_functions

__author__ = 'henrikef'


def test_constructor():

    # RA, Dec and L,B of the same point in the sky

    ra, dec = (125.6, -75.3)
    l, b = (288.44190139183564, -20.717313145391525)

    # This should throw an error as we are using Powerlaw instead of Powerlaw()
    with pytest.raises(RuntimeError):

        _ = ExtendedSource("my_source", Gaussian_on_sphere, Powerlaw)

    # This should throw an error because we should use a 2D function for the spatial shape
    with pytest.raises(RuntimeError):

        _ = ExtendedSource("my_source", Powerlaw(), Powerlaw())

    # Init with RA, Dec

    shape = Gaussian_on_sphere()
    source1 = ExtendedSource('my_source', shape, Powerlaw())
    shape.lon0=ra*u.degree
    shape.lat0=dec*u.degree

    assert source1.spatial_shape.lon0.value == ra
    assert source1.spatial_shape.lat0.value == dec

    # Verify that the position is free by default
    assert source1.spatial_shape.lon0.free
    assert source1.spatial_shape.lon0.free



def test_call():

    # Multi-component

    po1 = Powerlaw()
    po2 = Powerlaw()

    c1 = SpectralComponent("component1", po1)
    c2 = SpectralComponent("component2", po2)

    ra, dec = (125.6, -75.3)

    shape = Gaussian_on_sphere()
    source = ExtendedSource('test_source', shape, components=[c1, c2])
    shape.lon0=ra*u.degree
    shape.lat0=dec*u.degree

    assert np.all(source.spectrum.component1([1, 2, 3]) == po1([1, 2, 3]))
    assert np.all(source.spectrum.component2([1, 2, 3]) == po2([1, 2, 3]))

    one = source.spectrum.component1([1, 2, 3])
    two = source.spectrum.component2([1, 2, 3])

    #check spectral components
    assert np.all( np.abs(one + two - source.get_spatially_integrated_flux([1,2,3])) == 0 )
    
    #check spectral and spatial components
    total = source( [ra, ra, ra], [dec, dec, dec], [1,2,3])
    spectrum = one + two
    spatial = source.spatial_shape( [ra, ra, ra], [dec, dec, dec] )
    assert np.all( np.abs( total - spectrum*spatial ) == 0 )
    
    total = source( [ra*1.01, ra*1.01, ra*1.01], [dec*1.01, dec*1.01, dec*1.01], [1,2,3])
    spectrum = one + two
    spatial = source.spatial_shape( [ra*1.01, ra*1.01, ra*1.01], [dec*1.01, dec*1.01, dec*1.01] )
    assert np.all( np.abs( total - spectrum*spatial ) == 0 )
    


def test_call_with_units():

    po = Powerlaw()
    ga = Gaussian_on_sphere()

    result = ga(1.0, 1.0)
    
    assert result.ndim == 0

    with pytest.raises(u.UnitsError):

        # This raises because the units of the function have not been set up

        _ = ga(1.0 * u.deg, 1.0*u.deg)

    # Use the function as a spectrum
    source = ExtendedSource("test",ga,po)

    result = ga(1.0 * u.deg, 1.0*u.deg)

    assert isinstance(result, u.Quantity)

    result = ga(np.array([1,2,3]) * u.deg, np.array([1,2,3]) * u.deg) 

    assert isinstance(result, u.Quantity)

    # Now test all the functions
    def test_one(class_type):

        instance = class_type()

        if not instance.is_prior:

            # if we have fixed x_units then we will use those
            # in the test

            if instance.has_fixed_units():

                x_unit_to_use = instance.fixed_units[0]
                y_unit_to_use = instance.fixed_units[1]
                z_unit_to_use = instance.fixed_units[2]

            else:

                x_unit_to_use = u.deg
                y_unit_to_use = u.deg
                z_unit_to_use = u.keV



            # Use the function as a spectrum
            source = ExtendedSource("test",ga,po)


            #result = source(np.atleast_1d([1.0]) * x_unit_to_use,  np.atleast_1d([1.0]) * y_unit_to_use, np.atleast_1d([1.0])*z_unit_to_use)

            #assert isinstance(result, u.Quantity)

            result = source(np.array([1, 2, 3]) * x_unit_to_use, np.array([1, 2, 3]) * y_unit_to_use, np.array([1, 2, 3]) * z_unit_to_use)

            assert isinstance(result, u.Quantity)

            result = source(1.0, 1.0, 1.0)
            
            assert result.dtype == np.dtype('float64' )

        else:

            print ('Skipping prior function')


    for key in _known_functions:

        this_function = _known_functions[key]

        if this_function._n_dim in [2,3]:

            print("testing %s ..." % key)

            test_one(_known_functions[key])


def test_free_param():

    spectrum = Log_parabola()
    source = ExtendedSource("test_source", spatial_shape=Gaussian_on_sphere(), spectral_shape=spectrum)

    parameters = [spectrum.alpha, spectrum.beta, spectrum.piv, spectrum.K, source.spatial_shape.lat0, source.spatial_shape.lon0, source.spatial_shape.sigma]

    for param in parameters:
        param.free = False

    assert len(source.free_parameters) == 0

    for i, param in enumerate(parameters):
        param.free = True
        assert len(source.free_parameters) == i+1


