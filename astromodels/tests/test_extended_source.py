import astropy.units as u
import numpy as np
import pytest

from astromodels.core.spectral_component import SpectralComponent
from astromodels.functions.functions import Powerlaw, Log_parabola
from astromodels.functions.functions_2D import *
from astromodels.functions.functions_3D import *
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

    def test_one(class_type, name ):
        
      print("testing %s ..." % name)

      shape = class_type()
      source = ExtendedSource('test_source_%s' % name, shape, components=[c1, c2])
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
    
      total = source( [ra*1.01]*3, [dec*1.01]*3, [1,2,3])
      spectrum = one + two
      spatial = source.spatial_shape( [ra*1.01]*3, [dec*1.01]*3 )
      assert np.all( np.abs( total - spectrum*spatial ) == 0 )
    
    for key in _known_functions:

        if key in [ "SpatialTemplate_2D", "Latitude_galactic_diffuse"]:
        #not testing spatial template  as we don't have a file to read in right now.
            continue
        
        this_function = _known_functions[key]

        if this_function._n_dim == 2  and not this_function().is_prior:

            test_one(this_function, key)
            
    with pytest.raises(AssertionError):
        #this will fail because the Latitude_galactic_diffuse function isn't normalized.
        test_one(_known_functions["Latitude_galactic_diffuse"], "Latitude_galactic_diffuse")
   
@pytest.mark.xfail #fails for some functions.
def test_call_with_units():

    # Multi-component

    po1 = Powerlaw()
    po2 = Powerlaw()

    c1 = SpectralComponent("component1", po1)
    c2 = SpectralComponent("component2", po2)

    ra, dec = (125.6, -75.3)

    def test_one(class_type, name ):
        
      print("testing %s ..." % name)


      shape = class_type()
      source = ExtendedSource('test_source_%s' % name, spatial_shape = shape, components=[c1, c2])

      shape.lon0=ra*u.degree
      shape.lat0=dec*u.degree

      assert np.all(source.spectrum.component1([1, 2, 3]*u.keV) == po1([1, 2, 3]*u.keV))
      assert np.all(source.spectrum.component2([1, 2, 3]*u.keV) == po2([1, 2, 3]*u.keV))

      one = source.spectrum.component1([1, 2, 3]*u.keV)
      two = source.spectrum.component2([1, 2, 3]*u.keV)

      #check spectral components
      assert np.all( np.abs(one + two - source.get_spatially_integrated_flux([1,2,3]*u.keV)) == 0 )
    
      #check spectral and spatial components
      spatial = source.spatial_shape( ra*u.deg,dec*u.deg )
      spatial = source.spatial_shape( [ra, ra, ra]*u.deg, [dec, dec, dec]*u.deg )

      total = source( [ra, ra, ra]*u.deg, [dec, dec, dec]*u.deg, [1,2,3]*u.keV)
      spectrum = one + two
      assert np.all( np.abs( total - spectrum*spatial ) == 0 )
    
      total = source( [ra*1.01]*3*u.deg, [dec*1.01]*3*u.deg, [1,2,3]*u.keV)
      spectrum = one + two
      spatial = source.spatial_shape( [ra*1.01]*3*u.deg, [dec*1.01]*3*u.deg )
      assert np.all( np.abs( total - spectrum*spatial ) == 0 )
    
    for key in _known_functions:

        if key in [ "SpatialTemplate_2D", "Latitude_galactic_diffuse"]:
        #not testing spatial template  as we don't have a file to read in right now.
            continue
        
        this_function = _known_functions[key]

        if this_function._n_dim == 2  and not this_function().is_prior:

            test_one(this_function, key)
            
    with pytest.raises(AssertionError):
        #this will fail because the Latitude_galactic_diffuse function isn't normalized.
        test_one(_known_functions["Latitude_galactic_diffuse"], "Latitude_galactic_diffuse")
   

  
 
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


