__author__ = 'giacomov'

from astropy import coordinates
from astropy import units as u
import collections

from astromodels.parameter import Parameter
from astromodels.tree import Node


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

            if isinstance(ra, float):
                ra = Parameter('ra', ra, min_value=0.0, max_value=360.0, unit='deg')

            if isinstance(dec, float):
                dec = Parameter('dec', dec, min_value=-90.0, max_value=90.0, unit='deg')

            assert 0 <= ra.value <= 360, "R.A. cannot have a value of %s, it must be 0 <= ra <= 360" % ra
            assert -90 <= dec.value <= 90, "dec cannot have a value of %s, it must be -90 <= dec <= 90" % dec

            self._coord_type = 'equatorial'

            self.add_child(ra)
            self.add_child(dec)

        elif l is not None and b is not None:

            # This goes against duck typing, but it is needed to provide a means of initiating this class
            # with either Parameter instances or just floats

            if isinstance(l, float):
                l = Parameter('l', l, min_value=0.0, max_value=360.0, unit='deg')

            if isinstance(b, float):
                b = Parameter('b', b, min_value=-90.0, max_value=90.0, unit='deg')

            assert 0 <= l.value <= 360, "L cannot have a value of %s, it must be 0 <= L <= 360" % l
            assert -90 <= b.value <= 90, "B cannot have a value of %s, it must be -90 <= B <= 90" % b

            self._coord_type = 'galactic'
            self.add_child(l)
            self.add_child(b)

        else:

            raise WrongCoordinatePair("You have to specify either (ra, dec) or (l, b).")

    def get_ra(self):
        """
        Get R.A. corresponding to the current position (ICRS, J2000)

        :return: Right Ascension
        """

        try:

            return self.children['ra'].value

        except KeyError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to('icrs').ra.deg

    def get_dec(self):
        """
        Get Dec. corresponding to the current position (ICRS, J2000)

        :return: Declination
        """

        try:

            return self.children['dec'].value

        except KeyError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to('icrs').dec.deg

    def get_l(self):
        """
        Get Galactic Longitude (l) corresponding to the current position

        :return: Galactic Longitude
        """

        try:

            return self.children['l'].value

        except KeyError:

            # Transform from L,B to R.A., Dec

            return float(self.sky_coord.transform_to('galactic').l.deg)

    def get_b(self):
        """
        Get Galactic latitude (b) corresponding to the current position

        :return: Latitude
        """

        try:

            return self.children['b'].value

        except KeyError:

            # Transform from L,B to R.A., Dec

            return float(self.sky_coord.transform_to('galactic').b.deg)

    def _get_sky_coord(self):

        if self._coord_type == 'galactic':

            l = self.children['l'].value
            b = self.children['b'].value

            return coordinates.SkyCoord(l=l * u.deg, b=b * u.deg,
                                        frame='galactic', equinox=self._equinox,
                                        unit="deg")

        else:

            ra = self.children['ra'].value
            dec = self.children['dec'].value

            return coordinates.SkyCoord(ra=ra * u.deg, dec=dec * u.deg,
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

        return self.children

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

            data['ra'] = self.children['ra'].to_dict(minimal)
            data['dec'] = self.children['dec'].to_dict(minimal)
            data['equinox'] = self._equinox

        else:

            data['l'] = self.children['l'].to_dict(minimal)
            data['b'] = self.children['b'].to_dict(minimal)
            data['equinox'] = self._equinox

        return data

    def fix(self):
        """
        Fix the parameters with the coordinates (either ra,dec or l,b depending on how the class
        has been instanced)
        
        """

        if self._coord_type == 'equatorial':

            self.children['ra'].fix = True
            self.children['dec'].fix = True

        else:

            self.children['l'].fix = True
            self.children['b'].fix = True

    def free(self):
        """
        Free the parameters with the coordinates (either ra,dec or l,b depending on how the class
        has been instanced)
        
        """

        if self._coord_type == 'equatorial':

            self.children['ra'].fix = False
            self.children['dec'].fix = False

        else:

            self.children['l'].fix = False
            self.children['b'].fix = False

    @classmethod
    def from_dict(cls, data):

        return cls(**data)

    def _repr__base(self, rich_output):

        if self._coord_type == 'equatorial':

            representation = 'Sky direction (R.A., Dec.) = (%.5f, %.5f) (%s)' % (self.children['ra'].value,
                                                                                 self.children['dec'].value,
                                                                                 self.equinox)

        else:

            representation = 'Sky direction (l, b) = (%.5f, %.5f) (%s)' % (self.children['l'].value,
                                                                           self.children['b'].value,
                                                                           self.equinox)

        return representation