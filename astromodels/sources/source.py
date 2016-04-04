__author__ = 'giacomov'

import collections

from astromodels.named_object import NamedObject
from astromodels.dual_access_class import DualAccessClass


class Source(NamedObject):

    def __init__(self, source_name, list_of_components, src_type):

        NamedObject.__init__(self, source_name, allow_spaces=False)

        # Make the dictionary of components
        self._components = collections.OrderedDict()

        for component in list_of_components:
            self._components[component.name] = component

        # Store the type string
        self._src_type = str(src_type)

        # This will allow to access the components as instance.component, instead of instance.components['component']

        #DualAccessClass.__init__(self, "component", self._components)

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
