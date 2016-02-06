__author__ = 'giacomov'

import collections

from astromodels.named_object import NamedObject
from astromodels.dual_access_class import ProtectedAttribute


class SpectralComponent(NamedObject):

    # This is needed to avoid problems when constructing the class, due to the fact that we are overriding
    # the __setattr__ and __getattr__ attributes

    _spectral_shape = None

    def __init__(self, name, shape, polarization=None):

        NamedObject.__init__(self, name, allow_spaces=False)

        # Check that we can call the shape (i.e., it is a function)

        assert hasattr(shape, '__call__'), "The shape must be callable (i.e., behave like a function)"

        self._spectral_shape = shape

        self._polarization = polarization

    @property
    def shape(self):
        return self._spectral_shape

    @property
    def polarization(self):
        return self._polarization

    def __repr__(self):

        representation = "Spectral component %s\n" % self.name
        representation += "    -shape: %s\n" % self.shape.name

        # Print the polarization only if it's not None

        if self._polarization is not None:

            representation += "    -polarization: %s\n" % self.polarization.name

        return representation

    # Allow access to the spectral shape also by name. For example, if the shape is a powerlaw,
    # the user can access it either through component.shape or through component.powerlaw
    # Remember that __getattr__ is called only if __getattribute__ fails to find the requested attribute,
    # thus the performance impact of this is minimal

    def __getattr__(self, item):

        # The first check is needed to avoid infinite recursion during the constructor

        if self._spectral_shape is not None and item==self._spectral_shape.name:

            return self._spectral_shape

        else:

            raise AttributeError("Attribute %s does not exist" % item)

    # Implement this setattr to forbid the user from overwriting the attribute with the name of the shape.
    # For example, if the shape is a powerlaw, this will fail: component.powerlaw = 'something'.

    def __setattr__(self, key, value):

        # The first check is needed to avoid infinite recursion during the constructor

        if self._spectral_shape is not None and key==self.shape.name:

            raise ProtectedAttribute("You cannot overwrite the shape %s" % key)

        else:

            super(SpectralComponent, self).__setattr__(key, value)

    def to_dict(self):

        data = collections.OrderedDict()

        data['shape'] = self.shape.to_dict()

        if self.polarization is not None:

            data['polarization'] = self.polarization.to_dict()

        else:

            data['polarization'] = {}

        return data