from builtins import str
from builtins import object
__author__ = 'giacomov'

from typing import Dict, Optional, List, Any

from astromodels.core.parameter import Parameter

import collections

PARTICLE_SOURCE = 'particle source'
POINT_SOURCE = 'point source'
EXTENDED_SOURCE = 'extended source'


class UnknownSourceType(Exception):
    pass


class Source(object):

    def __init__(self, list_of_components: List[Any], src_type: str, spatial_shape=None):

        # Make the dictionary of components
        self._components: Dict[str, Any] = collections.OrderedDict()

        for component in list_of_components:

            self._components[component.name] = component

        if src_type not in (PARTICLE_SOURCE, POINT_SOURCE, EXTENDED_SOURCE):

            raise UnknownSourceType("Source of type %s is unknown" % src_type)

        else:

            # Store the type string
            self._src_type: str = str(src_type)

    def has_free_parameters(self):

        raise NotImplementedError("You need to override this")

    @property
    def free_parameters(self) -> Dict[str, Parameter]:
        """
        Returns a dictionary of free parameters for this source

        :return:
        """

        raise NotImplementedError("You need to override this")

    @property
    def components(self) -> Dict[str, Any]:
        """
        Return the dictionary of components

        :return: dictionary of components
        """

        return self._components

    @property
    def source_type(self) -> str:
        """
        Return the type of the source ('point source' or 'extended source')

        :return: type of the source
        """

        return self._src_type
