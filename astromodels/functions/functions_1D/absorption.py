import os
import sys
from functools import lru_cache, wraps

import astropy.units as astropy_units
import numpy as np
import numba as nb
import six
from astropy.io import fits

from astromodels.functions.function import Function1D, FunctionMeta
from astromodels.utils import configuration
from astromodels.utils.data_files import _get_data_file_path

from interpolation import interp


def cache_array_method(*args, **kwargs):
    """
    LRU cache implementation for methods whose FIRST parameter is a numpy array
    modified from: https://gist.github.com/Susensio/61f4fee01150caaac1e10fc5f005eb75
    """
    def decorator(function):
        @wraps(function)
        def wrapper(s, np_array, *args, **kwargs):
            hashable_array = tuple(np_array)
            return cached_wrapper(s, hashable_array, *args, **kwargs)

        @lru_cache(*args, **kwargs)
        def cached_wrapper(s, hashable_array, *args, **kwargs):
            array = np.array(hashable_array)
            return function(s, array, *args, **kwargs)

        # copy lru_cache attributes over too
        wrapper.cache_info = cached_wrapper.cache_info
        wrapper.cache_clear = cached_wrapper.cache_clear
        return wrapper

    return decorator


_abs_tables = {
    "phabs": {
        "AG89": "angr",
        "ASPL": "aspl"
    },
    "tbabs": {
        "AG89": "angr",
        "ASPL": "aspl",
        "WILM": "wilm"
    },
    "wabs": {
        "AG89": "angr"
    },
}
_abund_info = {}
_abund_info[
    "WILM"] = "wilms\nfrom Wilms, Allen & McCray (2000), ApJ 542, 914 \n except for elements not listed which are given zero abundance)\n https://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html "
_abund_info[
    "AG89"] = "angr\nfrom Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197)\n https://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html"
_abund_info[
    "ASPL"] = "aspl\nfrom Asplund M., Grevesse N., Sauval A.J. & Scott P. (2009, ARAA, 47, 481)\nhttps://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html"


def _get_xsect_table(model, abund_table):
    """
    contructs the abundance table from the values given
    """

    assert model in _abs_tables, "the model %s does not exist" % model
    assert abund_table in _abs_tables[model], ("the table %s does not exist" %
                                               abund_table)

    path_to_xsect = _get_data_file_path(
        os.path.join(
            "xsect",
            "xsect_%s_%s.fits" % (model, _abs_tables[model][abund_table])))

    fxs = fits.open(path_to_xsect)
    dxs = fxs[1].data
    xsect_ene = dxs["ENERGY"]
    xsect_val = dxs["SIGMA"]

    return np.array(xsect_ene,dtype=np.float64), np.array(xsect_val, dtype=np.float64)


# PhAbs class


class PhAbs(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        Photometric absorption (phabs implementation), f(E) = exp(- NH * sigma(E))
        contributed by Dominique Eckert
    parameters :
        NH :
            desc : absorbing column density in units of 1e22 particles per cm^2
            initial value : 1.0
            is_normalization : False
            transformation : log10
            min : 1e-4
            max : 1e4
            delta : 0.1

        redshift :
            desc : the redshift of the source
            initial value : 0.
            is_normalization : False
            min : 0
            max : 15
            delta : 0.1
            fix: True


    """
    def _setup(self):
        self._fixed_units = (astropy_units.keV,
                             astropy_units.dimensionless_unscaled)
        self.init_xsect()

    def _set_units(self, x_unit, y_unit):
        self.NH.unit = astropy_units.cm**(-2)
        self.redshift.unit = astropy_units.dimensionless_unscaled

    def init_xsect(self, abund_table="AG89"):
        """
        Set the abundance table

        :param abund_table: "ASPL", "AG89" 
        :returns: 
        :rtype: 

        """

        # load cross section data

        try:
            self.xsect_ene, self.xsect_val = _get_xsect_table(
                "phabs", abund_table)
            self._abund_table = abund_table

        except:

            print("defaulting to AG89")
            self.xsect_ene, self.xsect_val = _get_xsect_table(
                "phabs", abund_table)

            self._abund_table = "AG89"

    @cache_array_method()
    def _cached_interp(self, x):

        return np.interp(x, self.xsect_ene, self.xsect_val)

    def evaluate(self, x, NH, redshift):

        if isinstance(x, astropy_units.Quantity):

            _unit = astropy_units.cm**2
            _y_unit = astropy_units.dimensionless_unscaled
            _x = x.value
            _redshift = redshift.value
        else:

            _unit = 1.0
            _y_unit = 1.0
            _redshift = redshift
            _x = x

        xsect_interp = interp( self.xsect_ene, self.xsect_val,_x * (1 + _redshift))

        spec = np.exp(-NH * xsect_interp * _unit) * _y_unit

        return spec


# TbAbs class


class TbAbs(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        Photometric absorption (Tbabs implementation), f(E) = exp(- NH * sigma(E))
        contributed by Dominique Eckert
    parameters :
        NH :
            desc : absorbing column density in units of 1e22 particles per cm^2
            initial value : 1.0
            is_normalization : True
            transformation : log10
            min : 1e-4
            max : 1e4
            delta : 0.1

        redshift :
            desc : the redshift of the source
            initial value : 0.
            is_normalization : False
            min : 0
            max : 15
            delta : 0.1
            fix: True


    """
    def _setup(self):

        self.init_xsect()

        self._fixed_units = (astropy_units.keV,
                             astropy_units.dimensionless_unscaled)

    def _set_units(self, x_unit, y_unit):
        self.NH.unit = astropy_units.cm**(-2)
        self.redshift.unit = astropy_units.dimensionless_unscaled

    def init_xsect(self, abund_table="WILM"):
        """
        Set the abundance table

        :param abund_table: "WILM", "ASPL", "AG89" 
        :returns: 
        :rtype: 

        """

        try:
            self.xsect_ene, self.xsect_val = _get_xsect_table(
                "tbabs", abund_table)
            self._abund_table = abund_table

        except:

            print("defaulting to WILM")
            self.xsect_ene, self.xsect_val = _get_xsect_table(
                "tbabs", abund_table)

            self._abund_table = "WILM"

    @property
    def abundance_table(self):
        print(_abund_info[self._abund_table])

    @cache_array_method()
    def _cached_interp(self, x):

        return np.interp(x, self.xsect_ene, self.xsect_val)

    def evaluate(self, x, NH, redshift):

        if isinstance(x, astropy_units.Quantity):

            _unit = astropy_units.cm**2
            _y_unit = astropy_units.dimensionless_unscaled
            _x = x.value
            _redshift = redshift.value
            
        else:

            _unit = 1.0
            _y_unit = 1.0
            _redshift = redshift
            _x = x

        xsect_interp = interp( self.xsect_ene, self.xsect_val,_x * (1 + _redshift))

        spec = _numba_eval(NH, xsect_interp) * _y_unit

        return spec


# WAbs class


class WAbs(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        Photometric absorption (Wabs implementation), f(E) = exp(- NH * sigma(E))
        contributed by Dominique Eckert
    parameters :
        NH :
            desc : absorbing column density in units of 1e22 particles per cm^2
            initial value : 1.0
            is_normalization : True
            transformation : log10
            min : 1e-4
            max : 1e4
            delta : 0.1
        redshift :
            desc : the redshift of the source
            initial value : 0.
            is_normalization : False
            min : 0
            max : 15
            delta : 0.1
            fix: True


    """
    def _setup(self):
        self._fixed_units = (astropy_units.keV,
                             astropy_units.dimensionless_unscaled)
        self.init_xsect()

    def _set_units(self, x_unit, y_unit):
        self.NH.unit = astropy_units.cm**(-2)
        self.redshift.unit = astropy_units.dimensionless_unscaled

    def init_xsect(self):
        """
        Set the abundance table

        :returns:
        :rtype:

        """

        self.xsect_ene, self.xsect_val = _get_xsect_table("wabs", "AG89")

        self._abund_table = "AG89"

    @property
    def abundance_table(self):
        print(_abund_info[self._abund_table])

    @cache_array_method()
    def _cached_interp(self, x):

        return np.interp(x, self.xsect_ene, self.xsect_val)

    def evaluate(self, x, NH, redshift):

        if isinstance(x, astropy_units.Quantity):

            _unit = astropy_units.cm ** 2
            _y_unit = astropy_units.dimensionless_unscaled
            _x = x.value
            _redshift = redshift.value
        else:

            _unit = 1.0
            _y_unit = 1.0
            _redshift = redshift
            _x = x


        xsect_interp = interp( self.xsect_ene, self.xsect_val, _x * (1+ _redshift))

        spec = _numba_eval(NH,xsect_interp) * _y_unit

        return spec

@nb.njit(fastmath=True)
def _numba_eval(nh, xsect_interp):

    return np.exp(-nh * xsect_interp )
