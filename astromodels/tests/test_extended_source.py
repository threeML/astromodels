from __future__ import print_function

# this prevent a crash in macos. If does not import threeML first the code crashes
# with a segmantiation violation (Need to investigate more)s
try:
    from threeML import *
except:
    pass
import astropy.io.fits as fits
import astropy.units as u
from astropy import wcs
import numpy as np
import pytest

from astromodels.core.model import Model
from astromodels.core.model_parser import clone_model
from astromodels.core.spectral_component import SpectralComponent
from astromodels.functions import *
from astromodels.functions import Log_parabola, Powerlaw
from astromodels.functions.function import _known_functions
from astromodels.sources.extended_source import ExtendedSource

__author__ = 'henrikef'


def make_test_template(ra, dec, fitsfile):

    # Test template function: 40 pixel (0.8 deg) wide square centered approximately around a given ra, dec.
    test_wcs=False
    if (test_wcs):
        # this is an alternative way to build the header from WCS:
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [100, 100]
        w.wcs.cdelt = np.array([-0.02, 0.02])
        w.wcs.crval = [ra, dec]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        dOmega = (abs(w.wcs.cdelt[0] * w.wcs.cdelt[1]) *
                 u.degree * u.degree).to(u.steradian).value
        header=w.to_header()
    else:
        cards = {
            "SIMPLE": "T",
            "BITPIX": -32,
            "NAXIS": 2,
            "NAXIS1": 200,
            "NAXIS2": 200,
            "DATE": '2018-11-13',
            "CUNIT1": 'deg',
            "CRVAL1":  ra,
            "CRPIX1": 100,
            "CDELT1": -0.02,
            "CUNIT2": 'deg',
            "CRVAL2": dec,
            "CRPIX2": 100,
            "CDELT2": 0.02,
            "CTYPE1": 'RA---TAN',
            "CTYPE2": 'DEC--TAN'}

        dOmega = (abs(cards["CDELT1"] * cards["CDELT2"]) *
                  u.degree*u.degree).to(u.steradian).value
        header=fits.Header(cards)

    data = np.zeros([200, 200])
    data[80:120, 80:120] = 1

    total = np.sum(data)

    data = data / total / dOmega

    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(fitsfile, overwrite=True)


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
    shape.lon0 = ra*u.degree
    shape.lat0 = dec*u.degree

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

    def test_one(class_type, name):

        print("testing %s ..." % name)


        if name != "SpatialTemplate_2D":

            shape = class_type()
            source = ExtendedSource('test_source_%s' %
                                name, shape, components=[c1, c2])

            
            shape.lon0 = ra*u.degree
            shape.lat0 = dec*u.degree

        else:
            make_test_template(ra, dec, "__test.fits")
            shape = class_type(fits_file="__test.fits")
            source = ExtendedSource('test_source_%s' %
                                name, shape, components=[c1, c2])

            
            shape.K = 1.0

        assert np.all(source.spectrum.component1([1, 2, 3]) == po1([1, 2, 3]))
        assert np.all(source.spectrum.component2([1, 2, 3]) == po2([1, 2, 3]))

        one = source.spectrum.component1([1, 2, 3])
        two = source.spectrum.component2([1, 2, 3])

        # check spectral components
        assert np.all(
            np.abs(one + two - source.get_spatially_integrated_flux([1, 2, 3])) == 0)

        # check spectral and spatial components
        total = source([ra, ra, ra], [dec, dec, dec], [1, 2, 3])
        spectrum = one + two
        spatial = source.spatial_shape([ra, ra, ra], [dec, dec, dec])
        assert np.all(np.abs(total - spectrum*spatial) == 0)

        total = source([ra*1.01]*3, [dec*1.01]*3, [1, 2, 3])
        spectrum = one + two
        spatial = source.spatial_shape([ra*1.01]*3, [dec*1.01]*3)
        assert np.all(np.abs(total - spectrum*spatial) == 0)

    for key in _known_functions:

        if key in ["Latitude_galactic_diffuse"]:
            # not testing latitude galactic diffuse for now.
            continue

        this_function = _known_functions[key]

        if key in ["SpatialTemplate_2D"]:

            test_one(this_function, key)
        
        elif this_function._n_dim == 2 and not this_function().is_prior:

            test_one(this_function, key)

    with pytest.raises(AssertionError):
        # this will fail because the Latitude_galactic_diffuse function isn't normalized.
        test_one(
            _known_functions["Latitude_galactic_diffuse"], "Latitude_galactic_diffuse")


def test_call_with_units():

    # Multi-component

    po1 = Powerlaw()
    po2 = Powerlaw()

    c1 = SpectralComponent("component1", po1)
    c2 = SpectralComponent("component2", po2)

    ra, dec = (125.6, -75.3)

    def test_one(class_type, name):

        print("testing %s ..." % name)

        if name != "SpatialTemplate_2D":

            shape = class_type()
            source = ExtendedSource('test_source_%s' %
                                    name, spatial_shape=shape, components=[c1, c2])


            
            shape.lon0 = ra*u.degree
            shape.lat0 = dec*u.degree

        else:
            make_test_template(ra, dec, "__test.fits")

            shape = class_type(fits_file="__test.fits")
            source = ExtendedSource('test_source_%s' %
                                name, spatial_shape=shape, components=[c1, c2])

    
            shape.K = 1.0

        assert np.all(source.spectrum.component1(
            [1, 2, 3]*u.keV) == po1([1, 2, 3]*u.keV))
        assert np.all(source.spectrum.component2(
            [1, 2, 3]*u.keV) == po2([1, 2, 3]*u.keV))

        one = source.spectrum.component1([1, 2, 3]*u.keV)
        two = source.spectrum.component2([1, 2, 3]*u.keV)

        # check spectral components
        assert np.all(
            np.abs(one + two - source.get_spatially_integrated_flux([1, 2, 3]*u.keV)) == 0)

        # check spectral and spatial components
        #spatial = source.spatial_shape( ra*u.deg,dec*u.deg )
        spatial = source.spatial_shape(
            [ra, ra, ra]*u.deg, [dec, dec, dec]*u.deg)

        total = source([ra, ra, ra]*u.deg, [dec, dec, dec]
                       * u.deg, [1, 2, 3]*u.keV)
        spectrum = one + two
        assert np.all(np.abs(total - spectrum*spatial) == 0)

        total = source([ra*1.01]*3*u.deg, [dec*1.01]*3*u.deg, [1, 2, 3]*u.keV)
        spectrum = one + two
        spatial = source.spatial_shape([ra*1.01]*3*u.deg, [dec*1.01]*3*u.deg)
        assert np.all(np.abs(total - spectrum*spatial) == 0)

        model = Model(source)
        new_model = clone_model(model)

        new_total = new_model['test_source_%s' % name](
            [ra*1.01]*3*u.deg, [dec*1.01]*3*u.deg, [1, 2, 3]*u.keV)
        assert np.all(np.abs(total - new_total) == 0)

    for key in _known_functions:

        if key in ["Latitude_galactic_diffuse"]:
            # not testing latitude galactic diffuse for now.
            continue

        this_function = _known_functions[key]

        if key in ["SpatialTemplate_2D"]:

            test_one(this_function, key)
        
        elif this_function._n_dim == 2 and not this_function().is_prior:

            test_one(this_function, key)

    with pytest.raises(AssertionError):
        # this will fail because the Latitude_galactic_diffuse function isn't normalized.
        test_one(
            _known_functions["Latitude_galactic_diffuse"], "Latitude_galactic_diffuse")


def test_free_param():

    spectrum = Log_parabola()
    source = ExtendedSource(
        "test_source", spatial_shape=Gaussian_on_sphere(), spectral_shape=spectrum)

    parameters = [spectrum.alpha, spectrum.beta, spectrum.piv, spectrum.K,
                  source.spatial_shape.lat0, source.spatial_shape.lon0, source.spatial_shape.sigma]

    for param in parameters:
        param.free = False

    assert len(source.free_parameters) == 0

    for i, param in enumerate(parameters):
        param.free = True
        assert len(source.free_parameters) == i+1
