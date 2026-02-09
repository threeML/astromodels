__author__ = "giacomov"
from typing import Optional, Union

import collections

import astropy
from astropy import coordinates
import astropy.units as u

from astromodels.core.parameter import Parameter
from astromodels.core.tree import Node
from astromodels.core.units import get_units
from astromodels.utils.logging import setup_logger

log = setup_logger(__name__)


class WrongCoordinatePair(ValueError):
    pass


class IllegalCoordinateValue(ValueError):
    pass


class WrongCoordinateSystem(ValueError):
    pass


class SkyDirection(Node):
    """This is essentially a wrapper around astropy.coordinates.SkyCoord with a
    possibility for being serialized and deserialized with YAML."""

    def __init__(
        self,
        position: Optional[astropy.coordinates.SkyCoord] = None,
        ra: Optional[float] = None,
        dec: Optional[float] = None,
        l: Optional[float] = None,
        b: Optional[float] = None,
        unit: Union[astropy.units.Unit, str, None] = None,
        equinox: str = "J2000",
    ) -> None:
        """
        :param position: the position on the sky
        :type position: astropy.coordinates.SkyCoord
        :param ra: Right Ascension
        :type ra: float
        :param dec: Declination
        :type dec: float
        :param l: Galactic latitude
        :type l: float
        :param b: Galactic longitude
        :type b: float
        :param unit: astropy.units.Unit or string of the unit for the coordinates
        :type unit: astropy.units.Unit or str
        :param equinox: the equinox for the coordinates
        :type equinox: str

        :return:
        """

        self._equinox = equinox

        # Create the node

        Node.__init__(self, "position")

        # Check that we have the right pairs of coordinates

        if ra is not None or dec is not None:
            assert l is None and b is None, WrongCoordinatePair(
                "You have to specify either (ra, dec) or (l, b)."
            )

        frame = get_units().frame  # get the current coordinate_frame
        if isinstance(position, coordinates.SkyCoord):
            position = position.transform_to(frame)
            if frame == "icrs":
                ra = self._get_parameter_from_input(
                    position.ra.to(get_units().angle).value,
                    get_units.lon_bounds.lower_bound,
                    get_units.lon_bounds.upper_bound,
                    "ra",
                    "Right Ascension",
                )
                dec = self._get_parameter_from_input(
                    position.dec.to(get_units().angle).value,
                    get_units.lat_bounds.lower_bound,
                    get_units.lat_bounds.upper_bound,
                    "dec",
                    "Declination",
                )
                self._coord_type = "equatorial"
                self._add_child(ra)
                self._add_child(dec)
            elif frame == "galactic":

                lon = self._get_parameter_from_input(
                    position.l.to(get_units().angle).value,
                    get_units.lon_bounds.lower_bound,
                    get_units.lon_bounds.upper_bound,
                    "l",
                    "Galactic longitude",
                )
                lat = self._get_parameter_from_input(
                    position.b.to(get_units().angle).value,
                    get_units.lat_bounds.lower_bound,
                    get_units.lat_bounds.upper_bound,
                    "b",
                    "Galactic latitude",
                )

                self._coord_type = "galactic"
                self._add_child(lon)
                self._add_child(lat)
            else:
                raise NotImplementedError()

        elif ra is not None and dec is not None:

            # This goes against duck typing, but it is needed to provide a means of
            # initiating this class with either Parameter instances or just floats

            # Try to transform it to float, if it works than we transform it to a
            # parameter
            if unit is not None:
                if not isinstance(unit, u.Unit):
                    unit = u.Unit(unit)
                if not isinstance(ra, Parameter):
                    ra = (ra * unit).to(get_units().angle).value
                else:
                    ra = (ra.value * ra.unit).to(get_units().angle).value
                if not isinstance(dec, Parameter):
                    dec = (dec * unit).to(get_units().angle).value
                else:
                    dec = (dec.value * dec.unit).to(get_units().angle).value

            ra = self._get_parameter_from_input(
                ra,
                get_units.lon_bounds.lower_bound,
                get_units.lon_bounds.upper_bound,
                "ra",
                "Right Ascension",
            )

            dec = self._get_parameter_from_input(
                dec,
                get_units.lat_bounds.lower_bound,
                get_units.lat_bounds.upper_bound,
                "dec",
                "Declination",
            )

            self._coord_type = "equatorial"

            self._add_child(ra)
            self._add_child(dec)

        elif l is not None and b is not None:

            # This goes against duck typing, but it is needed to provide a means of
            # initiating this class with either Parameter instances or just floats

            # Try to transform it to float, if it works than we transform it to a
            # parameter
            if unit is not None:
                if not isinstance(unit, u.Unit):
                    unit = u.Unit(unit)
                if not isinstance(l, Parameter):
                    l = (l * unit).to(get_units().angle).value
                else:
                    l = (l.value * l.unit).to(get_units().angle).value
                if not isinstance(b, Parameter):
                    b = (b * unit).to(get_units().angle).value
                else:
                    b = (b.value * b.unit).to(get_units().angle).value

            l = self._get_parameter_from_input(
                l,
                get_units.lon_bounds.lower_bound,
                get_units.lon_bounds.upper_bound,
                "l",
                "Galactic longitude",
            )

            b = self._get_parameter_from_input(
                b,
                get_units.lat_bounds.lower_bound,
                get_units.lat_bounds.upper_bound,
                "b",
                "Galactic latitude",
            )

            self._coord_type = "galactic"
            self._add_child(l)
            self._add_child(b)

    @staticmethod
    def _get_parameter_from_input(number_or_parameter, minimum, maximum, name, desc):

        # Try to transform it to float, if it works than we transform it to a parameter

        try:

            number_or_parameter = float(number_or_parameter)

        except TypeError:

            assert isinstance(number_or_parameter, Parameter), (
                "%s must be either a number or a " "parameter instance" % name
            )

            # So this is a Parameter instance already. Enforce that it has the right
            # maximum and minimum

            parameter = number_or_parameter

            assert (
                parameter.min_value >= minimum
            ), "%s must have a minimum greater than or equal to %s" % (
                name,
                minimum,
            )
            assert (
                parameter.max_value <= maximum
            ), "%s must have a maximum less than or equal to %s" % (
                name,
                maximum,
            )

        else:

            # This was a float. Enforce that it has a legal value

            assert (
                minimum <= number_or_parameter <= maximum
            ), "%s cannot have a value of %s, " "it must be %s <= %s <= %s" % (
                name,
                number_or_parameter,
                minimum,
                name,
                maximum,
            )

            parameter = Parameter(
                name,
                number_or_parameter,
                desc=desc,
                min_value=minimum,
                max_value=maximum,
                unit=get_units().angle,
                free=False,
            )

        return parameter

    def get_ra(self):
        """Get R.A. corresponding to the current position (ICRS, J2000)

        :return: Right Ascension
        """

        try:

            return self.ra.value

        except AttributeError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to("icrs").ra.value

    def get_dec(self):
        """Get Dec. corresponding to the current position (ICRS, J2000)

        :return: Declination
        """

        try:

            return self.dec.value

        except AttributeError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to("icrs").dec.value

    def get_l(self):
        """Get Galactic Longitude (l) corresponding to the current position.

        :return: Galactic Longitude
        """

        try:

            return self.l.value

        except AttributeError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to("galactic").l.value

    def get_b(self):
        """Get Galactic latitude (b) corresponding to the current position.

        :return: Latitude
        """

        try:

            return self.b.value

        except AttributeError:

            # Transform from L,B to R.A., Dec

            return self.sky_coord.transform_to("galactic").b.value

    def _get_sky_coord(self):

        if self._coord_type == "galactic":

            l = self.l.value
            b = self.b.value

            return coordinates.SkyCoord(
                l=l, b=b, frame="galactic", equinox=self._equinox, unit="deg"
            )

        else:

            ra = self.ra.value
            dec = self.dec.value

            return coordinates.SkyCoord(
                ra=ra, dec=dec, frame="icrs", equinox=self._equinox, unit="deg"
            )

    @property
    def sky_coord(self):
        """Return an instance of astropy.coordinates.SkyCoord which can be used
        to make all transformations supported by it.

        :return: astropy.coordinates.SkyCoord
        """
        return self._get_sky_coord()

    @property
    def parameters(self):
        """Get the dictionary of parameters (either ra,dec or l,b)

        :return: dictionary of parameters
        """

        if self._coord_type == "galactic":

            return collections.OrderedDict((("l", self.l), ("b", self.b)))

        else:

            return collections.OrderedDict((("ra", self.ra), ("dec", self.dec)))

    @property
    def equinox(self):
        """Returns the equinox for the coordinates.

        :return:
        """
        return self._equinox

    def to_dict(self, minimal=False):

        data = collections.OrderedDict()

        if self._coord_type == "equatorial":

            data["ra"] = self.ra.to_dict(minimal)
            data["dec"] = self.dec.to_dict(minimal)
            data["equinox"] = self._equinox

        else:

            data["l"] = self.l.to_dict(minimal)
            data["b"] = self.b.to_dict(minimal)
            data["equinox"] = self._equinox

        return data

    def fix(self):
        """Fix the parameters with the coordinates (either ra,dec or l,b
        depending on how the class has been instanced)"""

        if self._coord_type == "equatorial":

            self.ra.fix = True
            self.dec.fix = True

        else:

            self.l.fix = True
            self.b.fix = True

    def free(self):
        """Free the parameters with the coordinates (either ra,dec or l,b
        depending on how the class has been instanced)"""

        if self._coord_type == "equatorial":

            self.ra.fix = False
            self.dec.fix = False

        else:

            self.l.fix = False
            self.b.fix = False

    @classmethod
    def from_dict(cls, data):

        return cls(**data)

    def _repr__base(self, rich_output):

        if self._coord_type == "equatorial":

            representation = "Sky direction (R.A., Dec.) = (%.5f, %.5f) (%s)" % (
                self.ra.value,
                self.dec.value,
                self.equinox,
            )

        else:

            representation = "Sky direction (l, b) = (%.5f, %.5f) (%s)" % (
                self.l.value,
                self.b.value,
                self.equinox,
            )

        return representation
