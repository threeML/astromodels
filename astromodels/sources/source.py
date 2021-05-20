
__author__ = 'giacomov'

from enum import Enum, unique
from astromodels.utils.logging import setup_logger
import collections

log = setup_logger(__name__)

@unique
class SourceType(Enum):
    PARTICLE_SOURCE = 'particle source'
    POINT_SOURCE = 'point source'
    EXTENDED_SOURCE = 'extended source'

    def __str__(self):
        return f"{self.value}"


class UnknownSourceType(Exception):
    pass


class Source(object):

    def __init__(self, list_of_components, src_type, spatial_shape=None):

        # Make the dictionary of components
        self._components = collections.OrderedDict()

        for component in list_of_components:

            self._components[component.name] = component

        if src_type not in SourceType:

            log.error(f"Source of type {src_type} is unknown")
            
            raise UnknownSourceType()

        else:

            # Store the type string
            self._src_type = src_type

    def has_free_parameters(self):

        raise NotImplementedError("You need to override this")

    @property
    def free_parameters(self):
        """
        Returns a dictionary of free parameters for this source

        :return:
        """

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
