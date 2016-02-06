__author__ = 'giacomov'

import collections

from astromodels.named_object import NamedObject
from astromodels.dual_access_class import DualAccessClass


class Source(NamedObject, DualAccessClass):

    def __init__(self, source_name, list_of_components):

        NamedObject.__init__(self, source_name, allow_spaces=False)

        # Make the dictionary of components
        components = collections.OrderedDict()

        for component in list_of_components:

            components[component.name] = component

        self._components = components

        # This will allow to access the components as instance.component, instead of instance.components['component']

        DualAccessClass.__init__(self, "component", self._components)

    @property
    def components(self):

        return self._components
