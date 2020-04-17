import astropy.units as astropy_units
import numpy as np

try:

    import pyatomdb

    has_atomdb = True

except:

    has_atomdb = False


from astromodels.functions.function import Function1D, FunctionMeta
from astromodels.utils import configuration
from astromodels.utils.data_files import _get_data_file_path

from astropy.io import fits
import os
import six


if has_atomdb:
    # APEC class
    @six.add_metaclass(FunctionMeta)
    class APEC(Function1D):
        r"""
        description :
            The Astrophysical Plasma Emission Code (APEC, Smith et al. 2001)
            contributed by Dominique Eckert
        parameters :
            K :
                desc : Normalization in units of 1e-14/(4*pi*(1+z)^2*dA*2)*EM
                initial value : 1.0
                is_normalization : True
                transformation : log10
                min : 1e-30
                max : 1e3
                delta : 0.1
            kT :
                desc : Plasma temperature
                initial value : 1.0
                min : 0.08
                max : 64
                delta : 0.1
            abund :
                desc : Metal abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            redshift :
                desc : Source redshift
                initial value : 0.1
                min : 0.0
                max : 10.0
                delta : 1e-3
                fix : yes

        """

        def _set_units(self, x_unit, y_unit):
            self.kT.unit = astropy_units.keV

            self.abund.unit = astropy_units.dimensionless_unscaled

            self.redshift.unit = astropy_units.dimensionless_unscaled

            self.K.unit = y_unit

        def init_session(self, abund_table="AG89"):
            # initialize PyAtomDB session
            self.session = pyatomdb.spectrum.CIESession(abundset=abund_table)

        def evaluate(self, x, K, kT, abund, redshift):
            sess = self.session

            nval = len(x)

            xz = x * (1.0 + redshift)

            ebplus = (np.roll(xz, -1) + xz)[: nval - 1] / 2.0

            ebounds = np.empty(nval + 1)

            ebounds[1:nval] = ebplus

            ebounds[0] = xz[0] - (ebplus[0] - xz[0])

            ebounds[nval] = xz[nval - 1] + (xz[nval - 1] - ebplus[nval - 2])

            binsize = (np.roll(ebounds, -1) - ebounds)[:nval]

            sess.set_response(ebounds, raw=True)

            sess.set_abund(
                [
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                ],
                abund,
            )

            spec = sess.return_spectrum(kT) / binsize / 1e-14

            return K * spec

    # APEC class
    @six.add_metaclass(FunctionMeta)
    class VAPEC(Function1D):
        r"""
        description :
            The Astrophysical Plasma Emission Code (APEC, Smith et al. 2001), variable
            contributed by Dominique Eckert
        parameters :
            K :
                desc : Normalization in units of 1e-14/(4*pi*(1+z)^2*dA*2)*EM
                initial value : 1.0
                is_normalization : True
                transformation : log10
                min : 1e-30
                max : 1e3
                delta : 0.1
            kT :
                desc : Plasma temperature
                initial value : 1.0
                min : 0.08
                max : 64
                delta : 0.1
            abund :
                desc : Metal abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            redshift :
                desc : Source redshift
                initial value : 0.1
                min : 0.0
                max : 10.0
                delta : 1e-3
                fix : yes

        """

        def _set_units(self, x_unit, y_unit):
            self.kT.unit = astropy_units.keV

            self.abund.unit = astropy_units.dimensionless_unscaled

            self.redshift.unit = astropy_units.dimensionless_unscaled

            self.K.unit = astropy_units.ph / astropy_units.s / astropy_units.keV

        def init_session(self, abund_table="AG89"):
            # initialize PyAtomDB session
            self.session = pyatomdb.spectrum.CIESession(abundset=abund_table)

        def evaluate(self, x, K, kT, abund, redshift):
            sess = self.session

            nval = len(x)

            xz = x * (1.0 + redshift)

            ebplus = (np.roll(xz, -1) + xz)[: nval - 1] / 2.0

            ebounds = np.empty(nval + 1)

            ebounds[1:nval] = ebplus

            ebounds[0] = xz[0] - (ebplus[0] - xz[0])

            ebounds[nval] = xz[nval - 1] + (xz[nval - 1] - ebplus[nval - 2])

            binsize = (np.roll(ebounds, -1) - ebounds)[:nval]

            sess.set_response(ebounds, raw=True)

            sess.set_abund(
                [
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                ],
                abund,
            )

            spec = sess.return_spectrum(kT) / binsize / 1e-14

            return K * spec


# PhAbs class
@six.add_metaclass(FunctionMeta)
class PhAbs(Function1D):
    r"""
    description :
        Photometric absorption (phabs implementation), f(E) = exp(- NH * sigma(E))
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
    """

    def _setup(self):
        self._fixed_units = (astropy_units.keV, astropy_units.dimensionless_unscaled)

    
    def _set_units(self, x_unit, y_unit):
        self.NH.unit = astropy_units.cm ** (-2)

    def init_xsect(self, abund_table="AG89"):
        """
        Set the abundance table

        :param abund_table: "ASPL", "AG89" 
        :returns: 
        :rtype: 

        """
        
        # load cross section data

        if abund_table == "AG89":
            path_to_xsect = _get_data_file_path(
                os.path.join("xsect", "xsect_phabs_angr.fits")
            )

        elif abund_table == "ASPL":
            path_to_xsect = _get_data_file_path(
                os.path.join("xsect", "xsect_phabs_aspl.fits")
            )

        else:
            print("Unknown abundace table %s, reverting to AG89" % (abund_table))
            path_to_xsect = _get_data_file_path(
                os.path.join("xsect", "xsect_phabs_angr.fits")
            )

        fxs = fits.open(path_to_xsect)
        dxs = fxs[1].data
        self.xsect_ene = dxs["ENERGY"]
        self.xsect_val = dxs["SIGMA"]

    def evaluate(self, x, NH):
        assert self.xsect_ene is not None and self.xsect_val is not None, "please run init_xsect(abund)"

        if isinstance(NH, astropy_units.Quantity):

            _unit = astropy_units.cm**2

        else:

            _unit = 1.
            
        xsect_interp = np.interp(x, self.xsect_ene, self.xsect_val)

        spec = np.exp(-NH * xsect_interp * _unit)

        return spec


# TbAbs class
@six.add_metaclass(FunctionMeta)
class TbAbs(Function1D):
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

    """

    def _setup(self):
        self._fixed_units = (astropy_units.keV, astropy_units.dimensionless_unscaled)
    
    def _set_units(self, x_unit, y_unit):
        self.NH.unit = astropy_units.cm ** (-2)

    def init_xsect(self, abund_table="WILM"):
        """
        Set the abundance table

        :param abund_table: "WILM", "ASPL", "AG89" 
        :returns: 
        :rtype: 

        """
        

        if abund_table == "AG89":
            path_to_xsect = _get_data_file_path(
                os.path.join("xsect", "xsect_tbabs_angr.fits")
            )

        elif abund_table == "ASPL":
            path_to_xsect = _get_data_file_path(
                os.path.join("xsect", "xsect_tbab_aspl.fits")
            )
        elif abund_table == "WILM":
            path_to_xsect = _get_data_file_path(
                os.path.join("xsect", "xsect_tbabs_wilm.fits")
            )

        else:
            print("Unknown abundace table %s, reverting to WILM" % (abund_table))
            path_to_xsect = _get_data_file_path(
                os.path.join("xsect", "xsect_tbabs_wilm.fits")
            )

        fxs = fits.open(path_to_xsect)
        dxs = fxs[1].data
        self.xsect_ene = dxs["ENERGY"]
        self.xsect_val = dxs["SIGMA"]

    def evaluate(self, x, NH):
        assert self.xsect_ene is not None and self.xsect_val is not None, "please run init_xsect(abund)"

        if isinstance(NH, astropy_units.Quantity):

            _unit = astropy_units.cm**2

        else:

            _unit = 1.


        xsect_interp = np.interp(x, self.xsect_ene, self.xsect_val)

        spec = np.exp(-NH * xsect_interp * _unit )

        return spec
