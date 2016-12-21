__author__ = 'giacomov'

import collections
import exceptions

PARTICLE_SOURCE = 'particle source'
POINT_SOURCE = 'point source'
EXTENDED_SOURCE = 'extended source'


class UnknownSourceType(exceptions.Exception):
    pass


class Source(object):

    def __init__(self, list_of_components, src_type, spatial_shape=None):

        # Make the dictionary of components
        self._components = collections.OrderedDict()

        for component in list_of_components:

            self._components[component.name] = component

        if src_type not in (PARTICLE_SOURCE, POINT_SOURCE, EXTENDED_SOURCE):

            raise UnknownSourceType("Source of type %s is unknown" % src_type)

        else:

            # Store the type string
            self._src_type = str(src_type)

    def has_free_parameters(self):

        raise NotImplementedError("You need to override this")

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
