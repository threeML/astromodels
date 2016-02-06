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

        self.parameters = collections.OrderedDict()

        # Check that we have the right pairs of coordinates

        if ra is not None and dec is not None:

            # This goes against duck typing, but it is needed to provide a means of initiating this class
            # with either Parameter instances or just floats

            if isinstance(ra, float):

                ra = Parameter('RA',ra,min_value=0.0,max_value=360.0)

            if isinstance(dec, float):

                dec = Parameter('Dec',dec,min_value=-90.0, max_value=90.0)

            assert 0 <= ra.value <= 360, "R.A. cannot have a value of %s, it must be 0 <= RA <= 360" % ra
            assert -90 <= dec.value <= 90, "Dec cannot have a value of %s, it must be -90 <= Dec <= 90" % dec

            self._coord_type = 'equatorial'
            self.parameters['RA'] = ra
            self.parameters['Dec'] = dec

            self._sky_pos = coordinates.SkyCoord(ra=ra.value * u.deg, dec=dec.value * u.deg,
                                                 frame='icrs', equinox=equinox)

            # Pre-compute galactic coordinates

            self.parameters['l'] = Parameter('l',self._sky_pos.galactic.l.value,min_value=0.0,max_value=360.0)
            self.parameters['b'] = Parameter('b',self._sky_pos.galactic.b.value,min_value=-90,max_value=+90.0)

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
            self.parameters['l'] = l
            self.parameters['b'] = b

            self._sky_pos = coordinates.SkyCoord(l=l.value * u.deg, b=b.value * u.deg,
                                                 frame='galactic', equinox=equinox)

            self.parameters['RA'] = Parameter('RA',self._sky_pos.icrs.ra.value,min_value=0.0,max_value=360.0)
            self.parameters['Dec'] = Parameter('Dec',self._sky_pos.icrs.dec.value,min_value=-90,max_value=+90.0)

        else:

            raise WrongCoordinatePair("You have to specify either (RA, Dec) or (l, b).")

    @property
    def ra(self):
        return self.parameters['RA'].value

    @property
    def dec(self):
        return self.parameters['Dec'].value

    @property
    def l(self):
        return self.parameters['l'].value

    @property
    def b(self):
        return self.parameters['b'].value

    def to_dict(self):

        data = collections.OrderedDict()

        if self._coord_type == 'equatorial':

            data['ra'] = self.parameters['RA'].to_dict()
            data['dec'] = self.parameters['Dec'].to_dict()
            data['equinox'] = self._equinox

        else:

            data['l'] = self.parameters['l'].to_dict()
            data['b'] = self.parameters['b'].to_dict()
            data['equinox'] = self._equinox

        return data

    @classmethod
    def from_dict(cls, data):

        return cls(**data)
