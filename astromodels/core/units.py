from builtins import object
__author__ = 'giacomov'

import collections

import astropy.units as u

from astromodels.core.tree import Node
from astromodels.utils.pretty_list import dict_to_list

# This module keeps the configuration of the units used in astromodels

# Pre-defined values

_ENERGY = u.keV
_TIME = u.s
_ANGLE = u.deg
_AREA = u.cm**2


class UnknownUnit(Exception):
    pass


class UnitMismatch(Exception):
    pass


def _check_unit(new_unit, old_unit):
    """
    Check that the new unit is compatible with the old unit for the quantity described by variable_name

    :param new_unit: instance of astropy.units.Unit
    :param old_unit: instance of astropy.units.Unit
    :return: nothin
    """

    try:

        new_unit.physical_type

    except AttributeError:

        raise UnitMismatch("The provided unit (%s) has no physical type. Was expecting a unit for %s"
                           % (new_unit, old_unit.physical_type))

    if new_unit.physical_type != old_unit.physical_type:

        raise UnitMismatch("Physical type mismatch: you provided a unit for %s instead of a unit for %s"
                           % (new_unit.physical_type, old_unit.physical_type))


class _AstromodelsUnits(object):
    """
    Store the fundamental units of time, energy, angle and area to be used in astromodels.
    """
    def __init__(self, energy_unit=None, time_unit=None, angle_unit=None, area_unit=None):

        if energy_unit is None: energy_unit = _ENERGY
        if time_unit is None: time_unit = _TIME
        if angle_unit is None: angle_unit = _ANGLE
        if area_unit is None: area_unit = _AREA

        self._units = collections.OrderedDict()

        self._units['energy'] = energy_unit
        self._units['time'] = time_unit
        self._units['angle'] = angle_unit
        self._units['area'] = area_unit

    # This __new__ method add the properties to the class. We could have achieved the same with a metaclass,
    # but this method is more clearer, for a tiny performance penalty. Consider also that under normal circumstances
    # this class will be only created once per session

    def __new__(cls, *args, **kwargs):

        cls.energy = property(*(cls._create_property('energy')))
        cls.time = property(*(cls._create_property('time')))
        cls.angle = property(*(cls._create_property('angle')))
        cls.area = property(*(cls._create_property('area')))

        obj = super(_AstromodelsUnits, cls).__new__(cls)

        return obj

    def get(self, what):

        return self._get_unit(what)

    def _set_unit(self, what, new_unit):

        try:

            old_unit = self._units[what]

        except KeyError:

            raise UnknownUnit("You can only assign units for energy, time, angle and area. Don't know "
                              "anything about %s" % what)

        # This allows to use strings in place of Unit instances as new_unit

        new_unit = u.Unit(new_unit)

        # Check that old and new unit are for the appropriate quantity

        _check_unit(new_unit, old_unit)

        # set the new unit
        self._units[what] = new_unit

    def _get_unit(self, what):

        try:

            return self._units[what]

        except KeyError:

            raise UnknownUnit("%s is not a fundamental unit" % what)

    # This is a function which generates the elements needed to make a property.photon
    # Just a trick to avoid having to duplicate the same code for each unit.
    # It is called in __new__

    @staticmethod
    def _create_property(what):

        return (lambda self: self._get_unit(what), lambda self, new_unit: self._set_unit(what, new_unit),
                "Sets or gets the unit for %s" % what)

    # Add the == and != operators
    def __eq__(self, other):

        return other.to_dict() == self.to_dict()

    def __ne__(self, other):

        return other.to_dict() != self.to_dict()

    def _repr__base(self, rich_output):

        return dict_to_list(self._units, html=rich_output)

    def to_dict(self, minimal=False):

        return self._units


# This is a factory which will always return the same instance of the _AstromodelsUnits class

class _AstromodelsUnitsFactory(object):

        _instance = None

        def __call__(self, *args, **kwds):

            if self._instance is None:

                # Create and return a new instance

                self._instance = _AstromodelsUnits(*args, **kwds)

                return self._instance

            else:

                # Use the instance already created

                return self._instance


# Create the factory to be used in the program

get_units = _AstromodelsUnitsFactory()