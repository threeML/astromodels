from __future__ import division

import collections
from typing import Any, Dict, List, Optional, Union

import astropy.units as u
import numba as nb
import numpy
import scipy.integrate
from past.utils import old_div

from astromodels.core.memoization import use_astromodels_memoization
from astromodels.core.parameter import Parameter
from astromodels.core.sky_direction import SkyDirection
from astromodels.core.spectral_component import SpectralComponent
from astromodels.core.tree import Node
from astromodels.core.units import get_units
from astromodels.functions.function import Function1D
from astromodels.sources.source import Source, SourceType
from astromodels.utils.logging import setup_logger
from astromodels.utils.pretty_list import dict_to_list

__author__ = 'giacomov'

__all__ = ["PointSource"]

log = setup_logger(__name__)

class PointSource(Source, Node):
    """
    A point source. You can instance this class in many ways.

    - with Equatorial position and a function as spectrum (the component will be automatically called 'main')::

        >>> from astromodels import *
        >>> point_source = PointSource('my_source', 125.6, -75.3, Powerlaw())

    - with Galactic position and a function as spectrum (the component will be automatically called 'main')::

        >>> point_source = PointSource('my_source', l=15.67, b=80.75, spectral_shape=Powerlaw())

    - with Equatorial position or Galactic position and a list of spectral components::

        >>> c1 = SpectralComponent("component1", Powerlaw())
        >>> c2 = SpectralComponent("component2", Powerlaw())
        >>> point_source = PointSource("test_source",125.6, -75.3,components=[c1,c2])

        Or with Galactic position:

        >>> point_source = PointSource("test_source",l=15.67, b=80.75,components=[c1,c2])
    
    NOTE: by default the position of the source is fixed (i.e., its positional parameters are fixed)
    
    :param source_name: name for the source
    :param ra: Equatorial J2000 Right Ascension (ICRS)
    :param dec: Equatorial J2000 Declination (ICRS)
    :param spectral_shape: a 1d function representing the spectral shape of the source
    :param l: Galactic latitude
    :param b: Galactic longitude
    :param components: list of spectral components (instances of SpectralComponent)
    :param sky_position: an instance of SkyDirection
    :return:
    """
    def __init__(self,
                 source_name: str,
                 ra: Optional[float] = None,
                 dec: Optional[float] = None,
                 spectral_shape: Optional[Function1D] = None,
                 l: Optional[float] = None,
                 b: Optional[float] = None,
                 components=None,
                 sky_position: Optional[SkyDirection]=None):

        # Check that we have all the required information

        # (the '^' operator acts as XOR on booleans)

        # Check that we have one and only one specification of the position

        if not ((ra is not None and dec is not None) ^
                (l is not None and b is not None) ^
                (sky_position is not None)):

            log.error(
                "You have to provide one and only one specification for the position"
            )

            raise AssertionError()

        # Gather the position

        if not isinstance(sky_position, SkyDirection):

            if (ra is not None) and (dec is not None):

                # Check that ra and dec are actually numbers

                try:

                    ra = float(ra)
                    dec = float(dec)

                except (TypeError, ValueError):

                    log.error("RA and Dec must be numbers. If you are confused by this message, you "
                                         "are likely using the constructor in the wrong way. Check the documentation.")
                    
                    raise AssertionError()

                sky_position = SkyDirection(ra=ra, dec=dec)

            else:

                sky_position = SkyDirection(l=l, b=b)

        self._sky_position: SkyDirection = sky_position

        # Now gather the component(s)

        # We need either a single component, or a list of components, but not both
        # (that's the ^ symbol)

        if not (spectral_shape is not None) ^ (components is not None):

            log.error("You have to provide either a single component, or a list of components (but not both).")

            raise AssertionError()

        # If the user specified only one component, make a list of one element with a default name ("main")

        if spectral_shape is not None:

            components = [SpectralComponent("main", spectral_shape)]

        Source.__init__(self, components, src_type=SourceType.POINT_SOURCE)

        # A source is also a Node in the tree

        Node.__init__(self, source_name)

        # Add the position as a child node, with an explicit name

        self._add_child(self._sky_position)

        # Add a node called 'spectrum'

        spectrum_node = Node('spectrum')
        spectrum_node._add_children(list(self._components.values()))

        self._add_child(spectrum_node)

        # Now set the units
        # Now sets the units of the parameters for the energy domain

        current_units = get_units()

        # Components in this case have energy as x and differential flux as y

        x_unit = current_units.energy
        y_unit = (current_units.energy * current_units.area * current_units.time) ** (-1)

        # Now set the units of the components
        for component in list(self._components.values()):

            component.shape.set_units(x_unit, y_unit)

    def __call__(self, x, tag=None):

        if tag is None:

            # No integration nor time-varying or whatever-varying

            if isinstance(x, u.Quantity):

                # Slow version with units

                results = [component.shape(x) for component in list(self.components.values())]

                # We need to sum like this (slower) because using np.sum will not preserve the units
                # (thanks astropy.units)

                return sum(results)

            else:

                # Fast version without units, where x is supposed to be in the same units as currently defined in
                # units.get_units()

                results = numpy.array([component.shape(x) for component in list(self.components.values())])

                return _sum(results)

        else:

            # Time-varying or energy-varying or whatever-varying

            integration_variable, a, b = tag

            if b is None:

                # Evaluate in a, do not integrate

                # Suspend memoization because the memoization gets confused when integrating
                with use_astromodels_memoization(False):

                    integration_variable.value = a

                    res = self.__call__(x, tag=None)

                return res

            else:

                # Integrate between a and b

                integrals = numpy.zeros(len(x))

                # TODO: implement an integration scheme avoiding the for loop

                # Suspend memoization because the memoization gets confused when integrating
                with use_astromodels_memoization(False):

                    reentrant_call = self.__call__

                    for i, e in enumerate(x):

                        def integral(y):

                            integration_variable.value = y

                            return reentrant_call(e, tag=None)

                        # Now integrate
                        integrals[i] = scipy.integrate.quad(integral, a, b, epsrel=1e-5)[0]

                return old_div(integrals, (b - a))

    @property
    def has_free_parameters(self) -> bool:
        """
        Returns True or False whether there is any parameter in this source

        :return:
        """

        for component in list(self._components.values()):

            for par in list(component.shape.parameters.values()):

                if par.free:

                    return True

        for par in list(self.position.parameters.values()):

            if par.free:

                return True

        return False

    @property
    def free_parameters(self) -> Dict[str, Parameter]:
        """
        Returns a dictionary of free parameters for this source.
        We use the parameter path as the key because it's 
        guaranteed to be unique, unlike the parameter name.

        :return:
        """
        free_parameters = collections.OrderedDict()

        for component in list(self._components.values()):

            for par in list(component.shape.parameters.values()):

                if par.free:

                    free_parameters[par.path] = par

        for par in list(self.position.parameters.values()):

            if par.free:

                free_parameters[par.path] = par

        return free_parameters

    @property
    def parameters(self) -> Dict[str, Parameter]:
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

        for par in self.position.parameters.values():

            all_parameters[par.path] = par

        return all_parameters

    def _repr__base(self, rich_output=False):
        """
        Representation of the object

        :param rich_output: if True, generates HTML, otherwise text
        :return: the representation
        """

        # Make a dictionary which will then be transformed in a list

        repr_dict = collections.OrderedDict()

        key = '%s (point source)' % self.name

        repr_dict[key] = collections.OrderedDict()
        repr_dict[key]['position'] = self._sky_position.to_dict(minimal=True)
        repr_dict[key]['spectrum'] = collections.OrderedDict()

        for component_name, component in list(self.components.items()):

            repr_dict[key]['spectrum'][component_name] = component.to_dict(minimal=True)

        return dict_to_list(repr_dict, rich_output)



@nb.njit(fastmath=True)
def _sum(x):
    return numpy.sum(x, axis=0)
    
