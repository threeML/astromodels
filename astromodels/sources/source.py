__author__ = 'giacomov'

import collections
import exceptions

from astromodels.units import get_units

PARTICLE_SOURCE = 'particle source'
POINT_SOURCE = 'point source'
EXTENDED_SOURCE = 'extended source'


class UnknownSourceType(exceptions.Exception):
    pass

class Source(object):

    def __init__(self, list_of_components, src_type):

        # Make the dictionary of components
        self._components = collections.OrderedDict()

        for component in list_of_components:

            self._components[component.name] = component

        # Store the type string
        self._src_type = str(src_type)

        # Now sets the units of the parameters for the energy domain

        current_units = get_units()

        if self._src_type == POINT_SOURCE:

            # Components in this case have energy as x and differential flux as y

            x_unit = current_units.energy
            y_unit = (current_units.energy * current_units.area * current_units.time) ** (-1)

        elif self._src_type == PARTICLE_SOURCE:

            # energy as x and particle flux as y
            x_unit = current_units.energy
            y_unit = 1 / current_units.energy

        elif self._src_type == EXTENDED_SOURCE:

            pass

        else:

            raise UnknownSourceType("Source of type %s is unknown" % self._src_type)

        # Now set the units of the components
        for component in self._components.values():

            component.shape.set_units(x_unit, y_unit)

    @property
    def components(self):
        """
        Return the dictionary of components

        :return: dictionary of components
        """

        return self._components

    @property
    def source_type(self):
        """
        Return the type of the source ('point source' or 'extended source')

        :return: type of the source
        """

        return self._src_type
