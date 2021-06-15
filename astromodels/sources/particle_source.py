from __future__ import division
from past.utils import old_div
__author__ = 'giacomov'

import collections

import numpy

from astromodels.core.spectral_component import SpectralComponent
from astromodels.core.tree import Node
from astromodels.core.units import get_units
from astromodels.sources.source import Source, SourceType
from astromodels.utils.pretty_list import dict_to_list
from astromodels.utils.logging import setup_logger


log = setup_logger(__name__)

class ParticleSource(Source, Node):
    """
    A source of particles with a certain distribution. This source does not produce any electromagnetic signal,
    but it is useful if used in conjuction with spectral shapes which need a particle distribution as input, like
    for example a Synchrotron kernel.

    :param name: name for the source
    :param distribution_shape: a function describing the energy distribution of the particles
    :param components: a list of SpectralComponents instances
    :return:

    """

    def __init__(self, name, distribution_shape=None, components=None):

        Node.__init__(self, name)

        if components is None:

            if distribution_shape is None:

                log.error("You have to either provied a list of components, or a " \
                                                   "distribution shape")

                raise AssertionError()

            components = [SpectralComponent("main", distribution_shape)]

        Source.__init__(self, components, SourceType.PARTICLE_SOURCE)

        # Add a node called 'spectrum'

        spectrum_node = Node('spectrum')
        spectrum_node._add_children(list(self._components.values()))

        self._add_child(spectrum_node)

        type(self).__call__ = type(self).get_flux

        # Set the units
        # Now sets the units of the parameters for the energy domain

        current_units = get_units()

        # energy as x and particle flux as y
        x_unit = current_units.energy
        y_unit = old_div(1, current_units.energy)

        # Now set the units of the components
        for component in list(self._components.values()):
            component.shape.set_units(x_unit, y_unit)

    def get_flux(self, energies):

        """Get the total flux of this particle source at the given energies (summed over the components)"""

        results = [component.shape(energies) for component in list(self.components.values())]

        return numpy.sum(results, 0)

    def _repr__base(self, rich_output=False):
        """
        Representation of the object

        :param rich_output: if True, generates HTML, otherwise text
        :return: the representation
        """

        # Make a dictionary which will then be transformed in a list

        repr_dict = collections.OrderedDict()

        key = '%s (particle source)' % self.name

        repr_dict[key] = collections.OrderedDict()
        repr_dict[key]['spectrum'] = collections.OrderedDict()

        for component_name, component in list(self.components.items()):

            repr_dict[key]['spectrum'][component_name] = component.to_dict(minimal=True)

        return dict_to_list(repr_dict, rich_output)

    @property
    def free_parameters(self):
        """
        Returns a dictionary of free parameters for this source.
        We use the parameter path as the key because it's 
        guaranteed to be unique, unlike the parameter name.

        :return:
        """
        free_parameters = collections.OrderedDict()

        for component in self._components.values():

            for par in component.shape.parameters.values():

                if par.free:

                    free_parameters[par.path] = par

        return free_parameters

    @property
    def parameters(self):
        """
        Returns a dictionary of all parameters for this source.
        We use the parameter path as the key because it's 
        guaranteed to be unique, unlike the parameter name.

        :return:
        """
        all_parameters = collections.OrderedDict()

        for component in self._components.values():

            for par in component.shape.parameters.values():

                all_parameters[par.path] = par

        return all_parameters


