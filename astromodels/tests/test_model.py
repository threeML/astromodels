import os
import pytest

__author__ = 'giacomov'

from astromodels.core.model import Model, DuplicatedNode, ModelFileExists, CannotWriteModel
from astromodels.sources.point_source import PointSource
from astromodels.sources.extended_source import ExtendedSource
from astromodels.sources.particle_source import ParticleSource
from astromodels.functions.functions import Powerlaw, _ComplexTestFunction, Uniform_prior
from astromodels.functions.functions_2D import Gaussian_on_sphere
from astromodels.core.parameter import Parameter, IndependentVariable
from astromodels.core.model_parser import *
from astromodels import u
import numpy as np


def _get_point_source(name="test"):

    pts = PointSource(name, ra=0, dec=0, spectral_shape=Powerlaw())

    return pts

def _get_point_source_gal(name):

    pts = PointSource(name, l=0, b=0, spectral_shape=Powerlaw())

    return pts

def _get_point_source_composite(name):

    spectral_shape = Powerlaw() + Powerlaw()

    pts = PointSource(name, l=0, b=0, spectral_shape=spectral_shape)

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
        pts2 = _get_point_source_gal("two")
        pts3 = _get_point_source_composite('three')
        ext1 = _get_extended_source("ext_one") # 6 params, 5 free
        ext2 = _get_extended_source("ext_two")
        ext3 = _get_extended_source("ext_three")

        # Two particle sources
        part1 = _get_particle_source("part_one")
        part2 = _get_particle_source("part_two")

        self._m = Model(pts1, pts2, pts3, ext1, ext2, ext3, part1, part2)

    @property
    def model(self):

        return self._m

    @property
    def n_sources(self):

        return 4

    @property
    def n_point_sources(self):

        return 3

    @property
    def n_extended_sources(self):

        return 3

    @property
    def n_particle_sources(self):

        return 2

    @property
    def total_number_of_parameters(self):

        return 5 * self.n_point_sources + 6 * self.n_extended_sources + 3 * self.n_particle_sources + 3

    @property
    def number_of_free_parameters(self):

        return 2 * self.n_point_sources + 5 * self.n_extended_sources + 2 * self.n_particle_sources + 2

    @property
    def center_of_extended_source(self):

        return (Gaussian_on_sphere.lon0.value, Gaussian_on_sphere.lat0.value)


def test_default_constructor():

    # Test that we cannot build a model with no sources

    with pytest.raises(AssertionError):

        _ = Model()


def test_constructor_1source():

    # Test with one point source
    pts = _get_point_source()

    m = Model(pts)

    # Test with a point source with an invalid name
    with pytest.raises(AssertionError):

        _ = PointSource("name", 0, 0, Powerlaw())

    assert len(m.sources) == 1


def test_constructor_duplicated_sources():

    # Test with one point source
    pts = _get_point_source()

    with pytest.raises(DuplicatedNode):

        m = Model(pts, pts)


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


def test_accessors_failures():

    mg = ModelGetter()

    m = mg.model

    for p in m:

        assert not((p.path + ".not") in m)


def test_display():

    mg = ModelGetter()
    m = mg.model

    m.display()

    m.display(complete=True)

    # Now display a model without free parameters
    m = Model(PointSource("test",0.0, 0.0, Powerlaw()))

    for parameter in m.parameters.values():

        parameter.fix = True

    m.display()

    # Now display a model without fixed parameters (very unlikely)
    for parameter in m.parameters.values():

        parameter.free = True

    m.display()


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


def test_external_parameters():

    mg = ModelGetter()

    m = mg.model

    n_free_before_external = len(m.free_parameters)

    # Create parameter
    fake_parameter = Parameter("external_parameter",1.0, min_value=-1.0, max_value=1.0, free=True)

    # Link as equal (default)
    m.add_external_parameter(fake_parameter)

    assert len(m.free_parameters) - 1 == n_free_before_external

    # Try to display it just to make sure it works

    m.display()

    # Now test that we can remove it
    m.remove_external_parameter("external_parameter")

    assert len(m.free_parameters) == n_free_before_external

    # Now test that adding it twice will not fail, but will issue a warning

    m.add_external_parameter(fake_parameter)

    with pytest.warns(RuntimeWarning):

        m.add_external_parameter(fake_parameter)


def test_input_output_basic():

    mg = ModelGetter()
    m = mg.model

    temp_file = "__test.yml"

    m.save(temp_file, overwrite=True)

    # Now reload it again
    m_reloaded = load_model(temp_file)

    # Check that all sources have been recovered
    assert m_reloaded.sources.keys() == m.sources.keys()

    os.remove(temp_file)

    # Now check that saving twice with the same name issues an exception if overwrite is not True
    m.save(temp_file, overwrite=True)

    with pytest.raises(ModelFileExists):

        m.save(temp_file, overwrite=False)

    # Try to write in an impossible location and verify that we throw an exception
    with pytest.raises(CannotWriteModel):

        m.save("/dev/null/ciaps", overwrite=True)

    # Add a prior to one of the parameters
    m.free_parameters.values()[0].prior = Uniform_prior()

    m.save(temp_file, overwrite=True)

    new_m = load_model(temp_file)

    assert m.free_parameters.values()[0].prior.to_dict() == new_m.free_parameters.values()[0].prior.to_dict()

    os.remove(temp_file)

    # Free and fix parameters
    m.one.position.ra.free = True
    m.two.spectrum.main.shape.K.fix = True

    m.save(temp_file, overwrite=True)

    new_m = load_model(temp_file)

    assert new_m.one.position.ra.free == True
    assert new_m.two.spectrum.main.shape.K.fix == True

    os.remove(temp_file)


def test_input_output_with_links():

    mg = ModelGetter()
    m = mg.model

    # Make a link
    m.link(m.one.spectrum.main.Powerlaw.K, m.two.spectrum.main.Powerlaw.K)

    temp_file = "__test.yml"

    m.save(temp_file, overwrite=True)

    # Now reload it again
    m_reloaded = load_model(temp_file)

    os.remove(temp_file)

    # Check that all sources have been recovered
    assert m_reloaded.sources.keys() == m.sources.keys()

    # Check that the link have been recovered
    new_value = 0.987
    m_reloaded.two.spectrum.main.Powerlaw.K.value = new_value

    assert m_reloaded.one.spectrum.main.Powerlaw.K.value == new_value


def test_input_output_with_external_parameters():

    mg = ModelGetter()
    m = mg.model

    # Create an external parameter
    fake_parameter = Parameter("external_parameter", 1.0, min_value=-1.0, max_value=1.0, free=True)

    # Link as equal (default)
    m.add_external_parameter(fake_parameter)

    # Save model

    temp_file = "__test.yml"

    m.save(temp_file, overwrite=True)

    # Now reload it again
    m_reloaded = load_model(temp_file)

    os.remove(temp_file)

    # Check that all sources have been recovered
    assert m_reloaded.sources.keys() == m.sources.keys()

    # Check that the external parameter have been recovered
    assert 'external_parameter' in m_reloaded

    # Remove external parameter
    m.remove_external_parameter("external_parameter")

    # Add a prior to the external parameter
    fake_parameter.prior = Uniform_prior()

    fake_parameter.value = -0.1

    m.add_external_parameter(fake_parameter)

    # Save model

    temp_file = "__test.yml"

    m.save(temp_file, overwrite=True)

    # Now reload it again
    m_reloaded = load_model(temp_file)

    os.remove(temp_file)

    # Check that all sources have been recovered
    assert m_reloaded.sources.keys() == m.sources.keys()

    # Check that the external parameter have been recovered
    assert 'external_parameter' in m_reloaded

    assert m.external_parameter.value == m_reloaded.external_parameter.value


def test_input_output_with_complex_functions():

    my_particle_distribution = Powerlaw()

    my_particle_distribution.index = -1.52

    electrons = ParticleSource('electrons', distribution_shape=my_particle_distribution)

    # Now set up the synch. spectrum for our source and the source itself

    synch_spectrum = _ComplexTestFunction()

    # Use the particle distribution we created as source for the electrons
    # producing synch. emission

    synch_spectrum.particle_distribution = my_particle_distribution

    synch_source = PointSource('synch_source', ra=12.6, dec=-13.5, spectral_shape=synch_spectrum)

    my_model = Model(electrons, synch_source)

    my_model.display()

    my_model.save("__test.yml")

    new_model = load_model("__test.yml")

    assert len(new_model.sources) == len(my_model.sources)

    assert my_particle_distribution.index.value == new_model.electrons.spectrum.main.shape.index.value


def test_add_remove_sources():

    mg = ModelGetter()
    m = mg.model

    # Point source

    new_src = _get_point_source("new")
    m.add_source(new_src)

    assert m.get_number_of_point_sources() == mg.n_point_sources + 1

    m.remove_source(new_src.name)

    assert m.get_number_of_point_sources() == mg.n_point_sources

    # Extended

    new_src = _get_extended_source("new")
    m.add_source(new_src)

    assert m.get_number_of_extended_sources() == mg.n_extended_sources + 1

    m.remove_source(new_src.name)

    assert m.get_number_of_extended_sources() == mg.n_extended_sources

    # Particle
    new_src = _get_particle_source("new")
    m.add_source(new_src)

    assert m.get_number_of_particle_sources() == mg.n_particle_sources + 1

    m.remove_source(new_src.name)

    assert m.get_number_of_particle_sources() == mg.n_particle_sources


def test_add_and_remove_independent_variable():

    mg = ModelGetter()
    m = mg.model

    # Create an independent variable
    independent_variable = IndependentVariable("time", 1.0, u.s)

    # Try to add it
    m.add_independent_variable(independent_variable)

    # Try to remove it
    m.remove_independent_variable("time")

    with pytest.raises(AssertionError):

        m.add_independent_variable(Parameter("time", 1.0))

    # Try to add it twice, which shouldn't fail
    m.add_independent_variable(independent_variable)
    m.add_independent_variable(independent_variable)

    # Try to display it just to make sure it works

    m.display()

    # Now try to use it
    link_law = Powerlaw()

    link_law.K.value = 1.0
    link_law.index.value = -1.0

    n_free_before_link = len(m.free_parameters)

    m.link(m.one.spectrum.main.Powerlaw.K, independent_variable, link_law)

    # The power law adds two parameters, but the link removes one, so
    assert len(m.free_parameters) - 1 == n_free_before_link

    # Now see if it works

    for t in np.linspace(0,10,100):

        independent_variable.value = t

        assert m.one.spectrum.main.Powerlaw.K.value == link_law(t)


def test_input_output_with_independent_variable():

    mg = ModelGetter()
    m = mg.model

    # Create an independent variable
    independent_variable = IndependentVariable("time", 1.0, u.s)

    # Try to add it
    m.add_independent_variable(independent_variable)

    # Save model

    temp_file = "__test.yml"

    m.save(temp_file, overwrite=True)

    # Now reload it again
    m_reloaded = load_model(temp_file)

    os.remove(temp_file)

    # Check that all sources have been recovered
    assert m_reloaded.sources.keys() == m.sources.keys()

    # Check that the external parameter have been recovered
    assert 'time' in m_reloaded


def test_3ML_interface():

    mg = ModelGetter()
    m = mg.model

    assert m.get_number_of_point_sources() == len(m.point_sources)
    assert m.get_number_of_particle_sources() == len(m.particle_sources)
    assert m.get_number_of_extended_sources() == len(m.extended_sources)

    # Test point source interface

    ra, dec = m.get_point_source_position(0)

    assert ra == m.point_sources.values()[0].position.ra.value
    assert dec == m.point_sources.values()[0].position.dec.value

    energies = np.logspace(1,2,10)
    fluxes = m.get_point_source_fluxes(0, energies)

    assert np.all(fluxes == m.point_sources.values()[0](energies))

    assert m.get_point_source_name(0) == m.point_sources.values()[0].name

    # Test extended source interface
    ra = np.random.uniform(0,1.0, 100)
    dec = np.random.uniform(0.0, 1.0, 100)
    energies = np.logspace(3,4,100)
    fluxes = m.get_extended_source_fluxes(0, ra, dec, energies)

    assert np.all(fluxes == m.extended_sources.values()[0](ra, dec, energies))

    assert m.get_extended_source_name(0) == m.extended_sources.values()[0].name

    res = m.get_extended_source_boundaries(0)

    res1 = m.extended_sources.values()[0].get_boundaries()

    assert np.all(np.array(res).flatten() == np.array(res1).flatten())

    assert m.is_inside_any_extended_source(0.0, 0.0) == True
    assert m.is_inside_any_extended_source(0.0, -90.0) == False

    # Test particle source interface

    energies = np.logspace(1, 2, 10)
    fluxes = m.get_particle_source_fluxes(0, energies)

    assert np.all(fluxes == m.particle_sources.values()[0](energies))

    assert m.get_particle_source_name(0) == m.particle_sources.values()[0].name


def test_clone_model():

    mg = ModelGetter()
    m1 = mg.model

    m2 = clone_model(m1)

    for source in m2.sources:

        assert source in m1.sources

    # Test that changing the parameter in one model does not changes the other

    m2.free_parameters.values()[0].value = m2.free_parameters.values()[0].value / 2.0

    assert m2.free_parameters.values()[0].value != m1.free_parameters.values()[0].value


def test_model_parser():

    mg = ModelGetter()
    m1 = mg.model

    m1.save("__test.yml")

    mp = ModelParser("__test.yml")

    with pytest.raises(ModelIOError):

        _ = ModelParser("__not_existing.yml")

    # Corrupt the yaml file
    with open("__test.yml", "a") as f:

        f.write("this is made to break the yaml parser")

    with pytest.raises(ModelYAMLError):

        _ = ModelParser("__test.yml")

    os.remove("__test.yml")
