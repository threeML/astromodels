__author__ = 'giacomov'

from astromodels.core.polarization import Polarization
from astromodels.core.tree import Node


class SpectralComponent(Node):

    # This is needed to avoid problems when constructing the class, due to the fact that we are overriding
    # the __setattr__ and __getattr__ attributes

    _spectral_shape = None

    def __init__(self, name, shape, polarization=None):

        # Check that we can call the shape (i.e., it is a function)

        assert hasattr(shape, '__call__'), "The shape must be callable (i.e., behave like a function)"

        self._spectral_shape = shape

        # Store the polarization

        if polarization is None:

            self._polarization = Polarization()

        else:

            self._polarization = polarization

        # Add shape and polarization as children

        Node.__init__(self, name)

        try:

            self._add_children([self._spectral_shape, self._polarization])

        except TypeError:

            raise TypeError("Couldn't instance the spectral component. Please verify that you are using an "
                            "*instance* of a function, and not a class. For example, you need to use "
                            "'%s()', and not '%s'." % (shape.__name__, shape.__name__))

    def _repr__base(self, rich_output):

        representation = "Spectral component %s\n" % self.name
        representation += "    -shape: %s\n" % self._spectral_shape.name

        # Print the polarization only if it's not None

        if self._polarization is not None:
            representation += "    -polarization: %s\n" % self._polarization.name

        return representation

    @property
    def shape(self):
        """
        A generic name for the spectral shape.

        :return: the spectral shape instance
        """
        return self._spectral_shape

    def __call__(self, energies):

        return self._spectral_shape(energies)
