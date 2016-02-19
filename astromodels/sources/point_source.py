__author__ = 'giacomov'

import collections
import numpy

from astromodels.sources.source import Source
from astromodels.sky_direction import SkyDirection
from astromodels.spectral_component import SpectralComponent


class PointSource(Source):
    """
    A point source. You can instance this class in many ways.

    - with Equatorial position and a function as spectrum (the component will be automatically called 'main')::

        point_source = PointSource('my_source', 125.6, -75.3, f1d.powerlaw())

    - with Galactic position and a function as spectrum (the component will be automatically called 'main')::

        point_source = PointSource('my_source', l=15.67, b=80.75, spectral_shape=f1d.powerlaw())

    - with Equatorial position or Galactic position and a list of spectral components::

        c1 = SpectralComponent("component1", f1d.powerlaw())
        c2 = SpectralComponent("component2", f1d.powerlaw())
        point_source = PointSource("test_source",125.6, -75.3,components=[c1,c2])
        # Or with Galactic position:
        point_source = PointSource("test_source",l=15.67, b=80.75,components=[c1,c2])

    - with a SkyDirection instance and either a single or a list of components::

        c1 = SpectralComponent("component1", f1d.powerlaw())
        c2 = SpectralComponent("component2", f1d.powerlaw())
        sky_dir = SkyDirection(RA=125.6, Dec=-75.3)
        point_source = PointSource("test_source",sky_direction=sky_dir,components=[c1,c2])
    
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

    def __init__(self, source_name, ra=None, dec=None, spectral_shape=None,
                 l=None, b=None, components=None, sky_position=None):

        # Check that we have all the required information

        # (the '^' operator acts as XOR on booleans)

        # Check that we have one and only one specification of the position

        assert ((ra is not None and dec is not None) ^
                (l is not None and b is not None) ^
                (sky_position is not None)), "You have to provide one and only one specification for the position"

        # Gather the position

        if not isinstance(sky_position, SkyDirection):

            if (ra is not None) and (dec is not None):

                # Check that ra and dec are actually numbers

                try:

                    ra = float(ra)
                    dec = float(dec)

                except (TypeError, ValueError):

                    raise AssertionError("RA and Dec must be numbers. If you are confused by this message, you "
                                         "are likely using the constructor in the wrong way. Check the documentation.")

                sky_position = SkyDirection(ra=ra, dec=dec)

            else:

                sky_position = SkyDirection(l=l, b=b)

        self._sky_position = sky_position

        # Fix the position by default

        self._sky_position.fix()

        # Now gather the component(s)

        # We need either a single component, or a list of components, but not both

        assert (spectral_shape is not None) ^ (components is not None), "You have to provide either a single " \
                                                                        "component, or a list of components " \
                                                                        "(but not both)."

        # If the user specified only one component, make a list of one element

        if spectral_shape is not None:
            components = [SpectralComponent("main", spectral_shape)]

        super(PointSource, self).__init__(source_name, components)

    @property
    def position(self):
        return self._sky_position

    def get_flux(self, energies):

        """Get the total flux of this point source at the given energies (summed over the components)"""

        results = [component.shape(energies) for component in self.components.values()]

        return numpy.sum(results, 0)

    def __repr__(self):

        representation = ''
        representation += 'Point source %s\n' % self.name
        representation += '    -position: (R.A., Dec) = (%s, %s)\n' % (self.position.ra, self.position.dec)
        representation += '    -components: %s\n' % ",".join(self.components.keys())

        return representation

    def to_dict(self):

        key = 'point source'

        data = collections.OrderedDict()

        # Serialize position

        position_key = 'position'

        position_dict = self.position.to_dict()

        data[position_key] = position_dict

        # Serialize components

        components_key = 'spectrum'

        components_dict = collections.OrderedDict()

        for component_name, component in self.components.iteritems():
            components_dict[component_name] = component.to_dict()

        data[components_key] = components_dict

        return {key: data}
