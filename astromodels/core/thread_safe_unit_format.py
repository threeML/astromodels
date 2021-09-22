# NOTE: the astropy.units.format pre-defined formatters use ply.yacc
# for parsing, which unfortunately is not thread safe. This causes
# all sort of weird problems when trying to use ipyparallel features.
# Here we implement a format which is very simple, does not use
# ply.yacc and is thread safe

from builtins import str
from builtins import map
from builtins import zip
import re
from astropy.units.format.base import Base
import astropy.units as u
from functools import reduce


# NOTE: the metaclass in Base will take care of registering
# this format, which will be available in the u.Unit
# constructor and in the .to_string method of the Unit and
# Quantity classes as format='threadsafe' (case insensitive)

def _format_one(xxx_todo_changeme):
    # This is for example 'cm-2' if base=cm and power=-2,
    # but only 'cm' if base=cm and power=1

    (base, power) = xxx_todo_changeme
    return "%s%s" % (base.to_string(), power if power != 1 else '')


class ThreadSafe(Base):

    @classmethod
    def parse(cls, s):

        # Since often we get these strings from a YAML file,
        # we want to make sure there is no new line within
        # the string.
        assert "\n" not in s

        # We assume that the string is something like
        # [unit1]{power} [unit2]{power} ..., like for example:
        # "m s-1" (a speed) or "keV-1 cm-2 s-1"
        # (a differential flux in 1 / (keV cm**2 s))
        # This is of course the format that is the output of
        # our to_string method(). See there for details

        tokens = re.findall('([a-zA-z]+)(-?\+?[0-9]+)?', s)

        # tokens is a list of tuples of the type [(unit name, power), ...]
        # Here we build a list like [u.m, u.s**(-1), ...]
        r = []

        for (unit, power) in tokens:

            # Get the primitive unit from the units module
            thisr = getattr(u, unit)

            # If needed, raise it to the power
            if power != '':

                thisr = thisr**power

            r.append(thisr)

        # Now we multiply the units
        # so that we get something like u.m * u.s**(-1)

        return reduce(lambda x, y: x*y, r)

    @classmethod
    def to_string(cls, unit):

        # The Unit class has two lists: bases (which are primitive units)
        # and powers (the correspondent powers). Thus we just need to use
        # those and encode them in a string like [unit1]{power} [unit2]{power} ...

        tokens = list(map(_format_one, list(zip(unit.bases, unit.powers))))

        return str(" ".join(tokens))
