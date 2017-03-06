__author__ = 'giacomov'

import collections

from astropy import coordinates

from astromodels.core.parameter import Parameter
from astromodels.core.tree import Node


class WrongCoordinatePair(ValueError):
    pass


class IllegalCoordinateValue(ValueError):
    pass


class WrongCoordinateSystem(ValueError):
    pass


class SkyDirection(Node):
    """
    This is essentially a wrapper around astropy.coordinates.SkyCoord with a possibility for
    being serialized and deserialized with YAML.
    """

    def __init__(self, ra=None, dec=None, l=None, b=None, equinox='J2000'):
        """

        :param ra: Right Ascension in degrees
        :param dec: Declination in degrees
        :param l: Galactic latitude in degrees
        :param b: Galactic longitude in degrees
        :param equinox: string
        :return:
        """

        self._equinox = equinox

        # Create the node

        Node.__init__(self, 'position')

        # Check that we have the right pairs of coordinates

        if ra is not None and dec is not None:

            # This goes against duck typing, but it is needed to provide a means of initiating this class
            # with either Parameter instances or just floats

            # Try to transform it to float, if it works than we transform it to a parameter

            ra = self._get_parameter_from_input(ra, 0, 360, 'ra','Right Ascension')

            dec = self._get_parameter_from_input(dec, -90, 90, 'dec','Declination')

            self._coord_type = 'equatorial'

            self._add_child(ra)
            self._add_child(dec)

        elif l is not None and b is not None:

            # This goes against duck typing, but it is needed to provide a means of initiating this class
            # with either Parameter instances or just floats

            # Try to transform it to float, if it works than we transform it to a parameter

            l = self._get_parameter_from_input(l, 0, 360, 'l','Galactic longitude')

            b = self._get_parameter_from_input(b, -90, 90, 'b','Galactic latitude')

            self._coord_type = 'galactic'
            self._add_child(l)
            self._add_child(b)

        else:

            raise WrongCoordinatePair("You have to specify either (ra, dec) or (l, b).")

    @staticmethod
    def _get_parameter_from_input(number_or_parameter, minimum, maximum, what, desc):

        # Try to transform it to float, if it works than we transform it to a parameter

        try:

            number_or_parameter = float(number_or_parameter)

        except TypeError:

            assert isinstance(number_or_parameter, Parameter), "%s must be either a number or a " \
                                                               "parameter instance" % what

            # So this is a Parameter instance already. Enforce that it has the right maximum and minimum

            parameter = number_or_parameter

            assert parameter.min_value == minimum, "%s must have a minimum of %s" % (what, minimum)
            assert parameter.max_value == maximum, "%s must have a maximum of %s" % (what, maximum)

        else:

            # This was a float. Enforce that it has a legal value

            assert minimum <= number_or_parameter <= maximum, "%s cannot have a value of %s, " \
                                                              "it must be %s <= %s <= %s" % (what, number_or_parameter,
                                                                                             minimum, what, maximum)

            parameter = Parameter(what, number_or_parameter,
                                  desc=desc, min_value=minimum, max_value=maximum, unit='deg', free=False)

        return parameter

    def get_ra(self):
        """
        Get R.A. corresponding to the current position (ICRS, J2000)

        :return: Right Ascension
        """

        try:

            return self.ra.value

        except AttributeError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to('icrs').ra.value

    def get_dec(self):
        """
        Get Dec. corresponding to the current position (ICRS, J2000)

        :return: Declination
        """

        try:

            return self.dec.value

        except AttributeError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to('icrs').dec.value

    def get_l(self):
        """
        Get Galactic Longitude (l) corresponding to the current position

        :return: Galactic Longitude
        """

        try:

            return self.l.value

        except AttributeError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to('galactic').l.value

    def get_b(self):
        """
        Get Galactic latitude (b) corresponding to the current position

        :return: Latitude
        """

        try:

            return self.b.value

        except AttributeError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to('galactic').b.value

    def _get_sky_coord(self):

        if self._coord_type == 'galactic':

            l = self.l.value
            b = self.b.value

            return coordinates.SkyCoord(l=l, b=b,
                                        frame='galactic', equinox=self._equinox,
                                        unit="deg")

        else:

            ra = self.ra.value
            dec = self.dec.value

            return coordinates.SkyCoord(ra=ra, dec=dec,
                                        frame='icrs', equinox=self._equinox,
                                        unit="deg")

    @property
    def sky_coord(self):
        """
        Return an instance of astropy.coordinates.SkyCoord which can be used to make all transformations supported
        by it

        :return: astropy.coordinates.SkyCoord
        """
        return self._get_sky_coord()

    @property
    def parameters(self):
        """
        Get the dictionary of parameters (either ra,dec or l,b)

        :return: dictionary of parameters
        """

        if self._coord_type == 'galactic':

            return collections.OrderedDict(('l', self.l), ('b', self.b))

        else:

            return collections.OrderedDict(('ra', self.ra), ('dec', self.dec))

    @property
    def equinox(self):
        """
        Returns the equinox for the coordinates.

        :return:
        """
        return self._equinox

    def to_dict(self, minimal=False):

        data = collections.OrderedDict()

        if self._coord_type == 'equatorial':

            data['ra'] = self.ra.to_dict(minimal)
            data['dec'] = self.dec.to_dict(minimal)
            data['equinox'] = self._equinox

        else:

            data['l'] = self.l.to_dict(minimal)
            data['b'] = self.b.to_dict(minimal)
            data['equinox'] = self._equinox

        return data

    def fix(self):
        """
        Fix the parameters with the coordinates (either ra,dec or l,b depending on how the class
        has been instanced)
        
        """

        if self._coord_type == 'equatorial':

            self.ra.fix = True
            self.dec.fix = True

        else:

            self.l.fix = True
            self.b.fix = True

    def free(self):
        """
        Free the parameters with the coordinates (either ra,dec or l,b depending on how the class
        has been instanced)
        
        """

        if self._coord_type == 'equatorial':

            self.ra.fix = False
            self.dec.fix = False

        else:

            self.l.fix = False
            self.b.fix = False

    @classmethod
    def from_dict(cls, data):

        return cls(**data)

    def _repr__base(self, rich_output):

        if self._coord_type == 'equatorial':

            representation = 'Sky direction (R.A., Dec.) = (%.5f, %.5f) (%s)' % (self.ra.value,
                                                                                 self.dec.value,
                                                                                 self.equinox)

        else:

            representation = 'Sky direction (l, b) = (%.5f, %.5f) (%s)' % (self.l.value,
                                                                           self.b.value,
                                                                           self.equinox)

        return representation