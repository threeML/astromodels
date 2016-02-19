__author__ = 'giacomov'

from astropy import coordinates
from astropy import units as u
import yaml
import collections
from astromodels.parameter import Parameter


class WrongCoordinatePair(ValueError):
    pass


class IllegalCoordinateValue(ValueError):
    pass


class WrongCoordinateSystem(ValueError):
    pass


class SkyDirection(object):
    """
    This is essentially a wrapper around astropy.coordinates.SkyCoord with a possibility for
    being serialized and deserialized with YAML.
    """

    def __init__(self, ra=None, dec=None, l=None, b=None, equinox='J2000'):
        """

        :param ra: float or Parameter
        :param dec: float or Parameter
        :param l: float or Parameter
        :param b: float or Parameter
        :param equinox: string
        :return:
        """

        self._equinox = equinox

        self._parameters = collections.OrderedDict()

        # Check that we have the right pairs of coordinates

        if ra is not None and dec is not None:

            # This goes against duck typing, but it is needed to provide a means of initiating this class
            # with either Parameter instances or just floats

            if isinstance(ra, float):

                ra = Parameter('ra',ra,min_value=0.0,max_value=360.0)

            if isinstance(dec, float):

                dec = Parameter('dec',dec,min_value=-90.0, max_value=90.0)

            assert 0 <= ra.value <= 360, "R.A. cannot have a value of %s, it must be 0 <= ra <= 360" % ra
            assert -90 <= dec.value <= 90, "dec cannot have a value of %s, it must be -90 <= dec <= 90" % dec

            self._coord_type = 'equatorial'
            self._parameters['ra'] = ra
            self._parameters['dec'] = dec

            # self._sky_pos = coordinates.SkyCoord(ra=ra.value * u.deg, dec=dec.value * u.deg,
            #                                    frame='icrs', equinox=equinox)

            # Pre-compute galactic coordinates

            # self._parameters['l'] = Parameter('l',self._sky_pos.galactic.l.value,min_value=0.0,max_value=360.0)
            # self._parameters['b'] = Parameter('b',self._sky_pos.galactic.b.value,min_value=-90,max_value=+90.0)

        elif l is not None and b is not None:

            # This goes against duck typing, but it is needed to provide a means of initiating this class
            # with either Parameter instances or just floats

            if isinstance(l, float):

                l = Parameter('l',l,min_value=0.0,max_value=360.0)

            if isinstance(b, float):

                b = Parameter('b',b,min_value=-90.0, max_value=90.0)

            assert 0 <= l.value <= 360, "L cannot have a value of %s, it must be 0 <= L <= 360" % l
            assert -90 <= b.value <= 90, "B cannot have a value of %s, it must be -90 <= B <= 90" % b

            self._coord_type = 'galactic'
            self._parameters['l'] = l
            self._parameters['b'] = b

            # self._sky_pos = coordinates.SkyCoord(l=l.value * u.deg, b=b.value * u.deg,
            #                                     frame='galactic', equinox=equinox)

            # self._parameters['ra'] = Parameter('ra',self._sky_pos.icrs.ra.value,min_value=0.0,max_value=360.0)
            # self._parameters['dec'] = Parameter('dec',self._sky_pos.icrs.dec.value,min_value=-90,max_value=+90.0)

        else:

            raise WrongCoordinatePair("You have to specify either (ra, dec) or (l, b).")

    def _get_sky_coord(self):

        if self._coord_type == 'galactic':

            l = self._parameters['l'].value
            b = self._parameters['b'].value

            return coordinates.SkyCoord(l=l* u.deg, b=b * u.deg,
                                        frame='galactic', equinox=self._equinox)

        else:

            ra = self._parameters['ra'].value
            dec = self._parameters['dec'].value

            return coordinates.SkyCoord(ra=ra * u.deg, dec=dec * u.deg,
                                        frame='icrs', equinox=self._equinox)

    @property
    def parameters(self):
        """
        Get the dictionary of parameters (either ra,dec or l,b)

        :return: dictionary of parameters
        """

        return self._parameters

    def set_ra(self, value):
        """
        Set the new R.A. for this sky direction

        :return: none
        """

        if self._coord_type == 'equatorial':

            self._parameters['ra'].value = value

        else:

            raise WrongCoordinateSystem('You cannot set R.A. since you have instanced the coordinates with l,b')

    def get_ra(self):
        """
        Get R.A. for this sky direction

        :return: Right Ascension
        """

        if self._coord_type == 'equatorial':

            return self._parameters['ra']

        else:

            sky_pos = self._get_sky_coord()

            return sky_pos.icrs.ra.value

    ra = property(get_ra,set_ra,doc="Get/set the new Right Ascension. Note that you can set R.A. only if you have "
                                    "instanced the object with the pair (R.A., Dec)")

    def set_dec(self, value):
        """
        Set the new Dec. for this sky direction

        :return: none
        """

        if self._coord_type == 'equatorial':

            self._parameters['dec'].value = value

        else:

            raise WrongCoordinateSystem('You cannot set Dec. since you have instanced the coordinates with l,b')

    def get_dec(self):
        """
        Get Declination

        :return: Declination
        """

        if self._coord_type == 'equatorial':

            return self._parameters['dec']

        else:

            sky_pos = self._get_sky_coord()

            return sky_pos.icrs.dec.value

    dec = property(get_dec,set_dec,doc="Get/set the new Declination. Note that you can set Dec. only if you have "
                                       "instanced the object with the pair (R.A., Dec)")

    def set_l(self, value):
        """
        Set the new L for this sky direction

        :return: none
        """

        if self._coord_type == 'galactic':

            self._parameters['l'].value = value

        else:

            raise WrongCoordinateSystem('You cannot set L since you have instanced the coordinates with R.A., Dec.')

    def get_l(self):
        """
        Get Galactic latitude

        :return: L
        """

        if self._coord_type == 'galactic':

            return self._parameters['l']

        else:

            sky_pos = self._get_sky_coord()

            return sky_pos.galactic.l.value

    l = property(get_l,set_l,doc="Get/set the new L. Note that you can set L only if you have "
                                 "instanced the object with the pair (L,B)")

    def set_b(self, value):
        """
        Set the new L for this sky direction

        :return: none
        """

        if self._coord_type == 'galactic':

            self._parameters['b'].value = value

        else:

            raise WrongCoordinateSystem('You cannot set B since you have instanced the coordinates with R.A., Dec.')

    def get_b(self):
        """
        Get Galactic longitude

        :return: b
        """

        if self._coord_type == 'galactic':

            return self._parameters['b']

        else:

            sky_pos = self._get_sky_coord()

            return sky_pos.galactic.b.value

    b = property(get_b,set_b,doc="Get/set the new B. Note that you can set B only if you have "
                                 "instanced the object with the pair (L,B)")

    def to_dict(self):

        data = collections.OrderedDict()

        if self._coord_type == 'equatorial':

            data['ra'] = self._parameters['ra'].to_dict()
            data['dec'] = self._parameters['dec'].to_dict()
            data['equinox'] = self._equinox

        else:

            data['l'] = self._parameters['l'].to_dict()
            data['b'] = self._parameters['b'].to_dict()
            data['equinox'] = self._equinox

        return data
    
    def fix(self):
        """
        Fix the parameters with the coordinates (either ra,dec or l,b depending on how the class
        has been instanced)
        
        """
        
        if self._coord_type == 'equatorial':
            
            self._parameters['ra'].fix = True
            self._parameters['dec'].fix = True
        
        else:
            
            self._parameters['l'].fix = True
            self._parameters['b'].fix = True
    
    def free(self):
        """
        Free the parameters with the coordinates (either ra,dec or l,b depending on how the class
        has been instanced)
        
        """
        
        if self._coord_type == 'equatorial':
            
            self._parameters['ra'].fix = False
            self._parameters['dec'].fix = False
        
        else:
            
            self._parameters['l'].fix = False
            self._parameters['b'].fix = False
        
    
    @classmethod
    def from_dict(cls, data):

        return cls(**data)
