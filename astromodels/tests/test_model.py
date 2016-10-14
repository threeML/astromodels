import pytest

import astropy.units as u

__author__ = 'giacomov'

from astromodels.model import Model, InvalidInput
from astromodels.sources.point_source import PointSource
from astromodels.sources.extended_source import ExtendedSource
from astromodels.sources.particle_source import ParticleSource
from astromodels.functions.functions import Powerlaw
from astromodels.functions.functions_2D import Gaussian_on_sphere
from astromodels.parameter import Parameter


def _get_point_source(name="test"):

    pts = PointSource(name, ra=0, dec=0, spectral_shape=Powerlaw())

    return pts


def _get_extended_source(name="test_ext"):

    ext = ExtendedSource(name, Gaussian_on_sphere(), Powerlaw())

    return ext


def _get_particle_source(name="test_part"):

    part = ParticleSource(name, Powerlaw())

    return part



class ModelGetter(object):

    def __init__(self):

        # 2 point sources and 3 ext sources
        pts1 = _get_point_source("one") # 5 params, 2 free
        pts2 = _get_point_source("two")
        ext1 = _get_extended_source("ext_one") # 6 params, 5 free
        ext2 = _get_extended_source("ext_two")
        ext3 = _get_extended_source("ext_three")

        # Two particle sources
        part1 = _get_particle_source("part_one")
        part2 = _get_particle_source("part_two")

        self._m = Model(pts1, pts2, ext1, ext2, ext3, part1, part2)

    @property
    def model(self):

        return self._m

    @property
    def n_sources(self):

        return 4

    @property
    def n_point_sources(self):

        return 2

    @property
    def n_extended_sources(self):

        return 3

    @property
    def n_particle_sources(self):

        return 2

    @property
    def total_number_of_parameters(self):

        return 5 * self.n_point_sources + 6 * self.n_extended_sources + 3 * self.n_particle_sources

    @property
    def number_of_free_parameters(self):

        return 2 * self.n_point_sources + 5 * self.n_extended_sources + 2 * self.n_particle_sources


def test_default_constructor():

    # Test that we cannot build a model with no sources

    with pytest.raises(AssertionError):

        _ = Model()


def test_constructor_1source():

    # Test with one point source
    pts = _get_point_source()

    m = Model(pts)

    # Test with a point source with an invalid name
    pts = PointSource("name", 0, 0, Powerlaw())

    with pytest.raises(InvalidInput):

        _ = Model(pts)

    # Test with two identical point sources, which should raise, as sources must have unique names
    many_sources = [pts] * 2

    with pytest.raises(InvalidInput):

        _ = Model(*many_sources)


def test_constructor_with_many_point_sources():

    # Test with 200 point sources

    many_p_sources = map(lambda x:_get_point_source("pts_source%i" %x), range(200))

    m = Model(*many_p_sources)

    assert m.get_number_of_point_sources()==200


def test_constructor_with_many_extended_sources():

    # Test with one extended source
    ext = _get_extended_source()

    _ = Model(ext)

    # Test with 200 extended sources
    many_e_sources = map(lambda x: _get_extended_source("ext_source%i" %x), range(200))

    m = Model(*many_e_sources)

    assert m.get_number_of_extended_sources() == 200


def test_constructor_with_mix():

    many_p_sources = map(lambda x: _get_point_source("pts_source%i" % x), range(200))

    many_e_sources = map(lambda x: _get_extended_source("ext_source%i" % x), range(200))

    many_part_sources = map(lambda x: _get_particle_source("part_source%i" % x), range(200))

    all_sources = []
    all_sources.extend(many_p_sources)
    all_sources.extend(many_e_sources)
    all_sources.extend(many_part_sources)

    m = Model(*all_sources)

    assert m.get_number_of_point_sources() == 200
    assert m.get_number_of_extended_sources() == 200
    assert m.get_number_of_particle_sources() == 200


def test_parameters_property():

    mg = ModelGetter()

    m = mg.model

    assert isinstance(m.parameters, dict)
    assert isinstance(m.free_parameters, dict)

    assert len(m.parameters) == mg.total_number_of_parameters
    assert len(m.free_parameters) == mg.number_of_free_parameters


def test_accessors():

    mg = ModelGetter()

    m = mg.model

    # Test __getitem__

    _ = m['ext_three.Gaussian_on_sphere.lat0']

    for p in m.parameters:

        assert isinstance(m[p], Parameter)

    # Test __iter__
    for p in m:

        assert isinstance(p, Parameter)

    # Test __contains__

    for p in m:

        assert p.path in m


def test_links():

    mg = ModelGetter()

    m = mg.model

    n_free_before_link = len(m.free_parameters)

    # Link as equal (default)
    m.link(m.one.spectrum.main.Powerlaw.K, m.two.spectrum.main.Powerlaw.K)

    assert len(m.free_parameters) == n_free_before_link -1

    # Try to display it just to make sure it works

    m.display()

    # Now test the link

    # This should print a warning, as trying to change the value of a linked parameters does not have any effect
    with pytest.warns(RuntimeWarning):

        m.one.spectrum.main.Powerlaw.K = 1.23456

    # This instead should work
    new_value = 1.23456
    m.two.spectrum.main.Powerlaw.K.value = new_value

    assert m.one.spectrum.main.Powerlaw.K.value == new_value

    # Now try to remove the link

    # First we remove it from the wrong parameters, which should issue a warning
    with pytest.warns(RuntimeWarning):

        m.unlink(m.two.spectrum.main.Powerlaw.K)

    # Remove it from the right parameter

    m.unlink(m.one.spectrum.main.Powerlaw.K)

    assert len(m.free_parameters) == n_free_before_link

    # Redo the same, but with a powerlaw law
    link_law = Powerlaw()

    link_law.K.value = 1.0
    link_law.index.value = -1.0

    n_free_before_link = len(m.free_parameters)

    m.link(m.one.spectrum.main.Powerlaw.K, m.two.spectrum.main.Powerlaw.K, link_law)

    # The power law adds two parameters, but the link removes one, so
    assert len(m.free_parameters) - 1 == n_free_before_link

    # Check that the link works
    new_value = 1.23456
    m.two.spectrum.main.Powerlaw.K.value = new_value

    predicted_value = link_law(new_value)

    assert m.one.spectrum.main.Powerlaw.K.value == predicted_value

    # Remove the link
    m.unlink(m.one.spectrum.main.Powerlaw.K)


def test_3ML_interface():

    pass

