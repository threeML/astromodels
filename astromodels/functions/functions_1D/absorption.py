import math
import os
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

import astropy.units as astropy_units
import numba as nb
import numpy as np
from astropy.io import fits
from interpolation import interp

from astromodels.functions.function import Function1D, FunctionMeta
from astromodels.utils import _get_data_file_path
from astromodels.utils.configuration import astromodels_config
from astromodels.utils.logging import setup_logger

log = setup_logger(__name__)

try:

    # ebltable is a Python packages to read in and interpolate tables for the photon density of
    # the Extragalactic Background Light (EBL) and the resulting opacity for high energy gamma
    # rays.

    import ebltable.tau_from_model as ebltau

    has_ebltable = True

except ImportError:

    if astromodels_config.logging.startup_warnings:

        log.warning(
            "The ebltable package is not available. Models that depend on it will not be available"
        )

    has_ebltable = False


class EBLTableNotAvailable(ImportWarning):
    pass


class InvalidUsageForFunction(Exception):
    pass


@dataclass(frozen=False)
class AbundanceTable:
    name: str
    tables: List[str]
    _current_table: str

    def set_table(self, table: str):

        """
        set the current table from AG89, WILM or ASPL

        :param table:
        :type table: str
        :returns:

        """

        old_table = self._current_table

        self._current_table = table

        if self.current_table not in self.tables:

            log.error(
                f"{self.name} does not contain {table} choose {','.join(self.table)}"
            )

            self._current_table = old_table

            raise AssertionError()

    @property
    def current_table(self) -> str:

        convert = {"AG89": "angr", "ASPL": "aspl", "WILM": "wilm"}

        return convert[self._current_table]

    @property
    def info(self) -> str:

        _abund_info = {}
        _abund_info[
            "WILM"
        ] = "wilms\nfrom Wilms, Allen & McCray (2000), ApJ 542, 914 \n except for elements not listed which are given zero abundance)\n https://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html "
        _abund_info[
            "AG89"
        ] = "angr\nfrom Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197)\n https://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html"
        _abund_info[
            "ASPL"
        ] = "aspl\nfrom Asplund M., Grevesse N., Sauval A.J. & Scott P. (2009, ARAA, 47, 481)\nhttps://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html"

        return _abund_info[self._current_table]

    @property
    def xsect_table(self) -> np.ndarray:

        """
        returns the XSECT table for the current model

        :returns:

        """
        _path: Path = (
            Path("xsect") / f"xsect_{self.name}_{self.current_table}.fits"
        )

        path_to_xsect: Path = _get_data_file_path(_path)

        with fits.open(path_to_xsect) as fxs:

            dxs = fxs[1].data
            xsect_ene = dxs["ENERGY"]
            xsect_val = dxs["SIGMA"]

        return np.array(xsect_ene, dtype=np.float64), np.array(
            xsect_val, dtype=np.float64
        )


phabs = AbundanceTable("phabs", ["angr", "aspl"], "AG89")
tbabs = AbundanceTable("tbabs", ["angr", "aspl", "wilm"], "WILM")
wabs = AbundanceTable("wabs", ["angr"], "AG89")


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

    properties:
         abundance_table:
            desc: the abundance table for the model
            initial value: AG89
            allowed values:
              - AG89
              - ASPL
            function: _init_xsect

    """

    def _setup(self):
        self._fixed_units = (
            astropy_units.keV,
            astropy_units.dimensionless_unscaled,
        )

        # # astromodels_config.absorption_models.phabs_table.value
        # self.init_xsect(self.abundance_table.value)

    def _set_units(self, x_unit, y_unit):
        self.NH.unit = astropy_units.cm ** (-2)
        self.redshift.unit = astropy_units.dimensionless_unscaled

    def _init_xsect(self):
        """
        Set the abundance table

        :param abund_table: "ASPL", "AG89"
        :returns:
        :rtype:

        """

        # load cross section data

        phabs.set_table(self.abundance_table.value)

        self.xsect_ene, self.xsect_val = phabs.xsect_table

    @property
    def abundance_table_info(self):
        print(phabs.info)

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

        xsect_interp = interp(
            self.xsect_ene, self.xsect_val, _x * (1 + _redshift)
        )

        # evaluate the exponential with numba

        spec = _numba_eval(NH, xsect_interp) * _y_unit

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

    properties:
         abundance_table:
            desc: the abundance table for the model
            initial value: WILM
            allowed values:
             - WILM
             - AG89
             - ASPL
            function: _init_xsect



    """

    def _setup(self):

        # self.init_xsect(self.abundance_table)

        self._fixed_units = (
            astropy_units.keV,
            astropy_units.dimensionless_unscaled,
        )

    def _set_units(self, x_unit, y_unit):

        self.NH.unit = astropy_units.cm ** (-2)
        self.redshift.unit = astropy_units.dimensionless_unscaled

    def _init_xsect(self):
        """
        Set the abundance table

        :param abund_table: "WILM", "ASPL", "AG89"
        :returns:
        :rtype:

        """

        tbabs.set_table(self.abundance_table.value)

        self.xsect_ene, self.xsect_val = tbabs.xsect_table

    @property
    def abundance_table_info(self):
        print(tbabs.info)

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

        xsect_interp = interp(
            self.xsect_ene, self.xsect_val, _x * (1 + _redshift)
        )

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
        self._fixed_units = (
            astropy_units.keV,
            astropy_units.dimensionless_unscaled,
        )
        self._init_xsect()

    def _set_units(self, x_unit, y_unit):
        self.NH.unit = astropy_units.cm ** (-2)
        self.redshift.unit = astropy_units.dimensionless_unscaled

    def _init_xsect(self):
        """
        Set the abundance table

        :returns:
        :rtype:

        """

        self.xsect_ene, self.xsect_val = wabs.xsect_table

    @property
    def abundance_table_info(self):
        print(wabs.info)

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

        xsect_interp = interp(
            self.xsect_ene, self.xsect_val, _x * (1 + _redshift)
        )

        spec = _numba_eval(NH, xsect_interp) * _y_unit

        return spec


if has_ebltable:

    class EBLattenuation(Function1D, metaclass=FunctionMeta):
        r"""
        description :
            Attenuation factor for absorption in the extragalactic background light (EBL) ,
            to be used for extragalactic source spectra. Based on package "ebltable" by
            Manuel Meyer, https://github.com/me-manu/ebltable .

        latex: not available

        parameters :

          redshift :
                desc : redshift of the source
                initial value : 1.0
                fix : yes

          attenuation :
                desc : scaling factor for the strength of attenuation
                initial value : 1.0
                min : 0.0
                max : 10.0
                fix : yes

        properties:
           ebl_model:
              desc: set the EBL model
              initial value: dominguez
              allowed values:
                 - dominguez
                 - franceschini
                 - kneiske
                 - inuoe
                 - gilmore
              function: _set_ebl_model


        """

        # def _setup(self):

        #     # define EBL model, use dominguez as default
        #     self._tau = ebltau.OptDepth.readmodel(model=astromodels_config.absorption_models.ebl_table.value)

        def _set_ebl_model(self):

            # passing modelname to ebltable, which will check if defined
            self._tau = ebltau.OptDepth.readmodel(model=self.ebl_model.value)

        def _set_units(self, x_unit, y_unit):

            if (
                not hasattr(x_unit, "physical_type")
                or x_unit.physical_type != "energy"
            ):

                # x should be energy
                raise InvalidUsageForFunction(
                    "Unit for x is not an energy. The function "
                    "EBLOptDepth calculates energy-dependent "
                    "absorption."
                )

            # y should be dimensionless
            if (
                not hasattr(y_unit, "physical_type")
                or y_unit.physical_type != "dimensionless"
            ):
                raise InvalidUsageForFunction(
                    "Unit for y is not dimensionless."
                )

            self.redshift.unit = astropy_units.dimensionless_unscaled
            self.attenuation.unit = astropy_units.dimensionless_unscaled

        def evaluate(self, x, redshift, attenuation):

            if isinstance(x, astropy_units.Quantity):

                # ebltable expects TeV
                eTeV = x.to(astropy_units.TeV).value
                _unit = astropy_units.dimensionless_unscaled
                _redshift = redshift.value
                _attenuation = attenuation.value

            else:

                # otherwise it's in keV
                eTeV = x / 1.0e9

                _unit = 1.0
                _redshift = redshift
                _attenuation = attenuation

            return (
                _numba_eval(self._tau.opt_depth(_redshift, eTeV), _attenuation)
                * _unit
            )


@nb.vectorize
def _exp(x):
    return math.exp(x)


@nb.njit(fastmath=True)
def _numba_eval(nh, xsect_interp):

    return _exp(-nh * xsect_interp)
