from functools import wraps


try:
    from functools import lru_cache
except ImportError:

    def lru_cache(user_function, **kwargs):
        cache = {}

        @wraps(user_function)
        def wrapper(*args):
            key = tuple(args)
            if key not in cache:
                cache[key] = user_function(*args)
            return cache[key]

        return wrapper


import astropy.units as astropy_units
import numpy as np

import sys





is_py3 = True

if sys.version_info[0] < 3:

    is_py3 = False

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
            assert self.session is not None, "please run init_session(abund)"

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

    # VAPEC class
    @six.add_metaclass(FunctionMeta)
    class VAPEC(Function1D):
        r"""
        description :
            The Astrophysical Plasma Emission Code (APEC, Smith et al. 2001), variable abundances for individual elements
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
            Fe :
                desc : Fe abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            C :
                desc : C abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            N :
                desc : N abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            O :
                desc : O abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            Ne :
                desc : Ne abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            Mg :
                desc : Mg abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            Al :
                desc : Al abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            Si :
                desc : Si abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            S :
                desc : S abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            Ar :
                desc : Ar abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            Ca :
                desc : Ca abundance
                initial value : 1
                min : 0.0
                max : 5.0
                delta : 0.01
                fix : yes
            Ni :
                desc : Ni abundance
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

            self.Fe.unit = astropy_units.dimensionless_unscaled

            self.C.unit = astropy_units.dimensionless_unscaled

            self.N.unit = astropy_units.dimensionless_unscaled

            self.O.unit = astropy_units.dimensionless_unscaled

            self.Ne.unit = astropy_units.dimensionless_unscaled

            self.Mg.unit = astropy_units.dimensionless_unscaled

            self.Al.unit = astropy_units.dimensionless_unscaled

            self.Si.unit = astropy_units.dimensionless_unscaled

            self.Ar.unit = astropy_units.dimensionless_unscaled

            self.Ca.unit = astropy_units.dimensionless_unscaled

            self.Ni.unit = astropy_units.dimensionless_unscaled

            self.redshift.unit = astropy_units.dimensionless_unscaled

            self.K.unit = y_unit

        def init_session(self, abund_table="AG89"):
            # initialize PyAtomDB session
            self.session = pyatomdb.spectrum.CIESession(abundset=abund_table)

        def evaluate(
            self, x, K, kT, Fe, C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Ni, redshift
        ):
            assert self.session is not None, "please run init_session(abund)"

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
                [6,], C,
            )

            sess.set_abund(
                [7,], N,
            )

            sess.set_abund(
                [8,], O,
            )

            sess.set_abund(
                [10,], Ne,
            )

            sess.set_abund(
                [12,], Mg,
            )

            sess.set_abund(
                [13,], Al,
            )

            sess.set_abund(
                [14,], Si,
            )

            sess.set_abund(
                [16,], S,
            )

            sess.set_abund(
                [18,], Ar,
            )

            sess.set_abund(
                [20,], Ca,
            )

            sess.set_abund(
                [26,], Fe,
            )

            sess.set_abund(
                [28,], Ni,
            )

            sess.set_abund(
                [9, 11, 15, 17, 19, 21, 22, 23, 24, 25, 27, 29, 30], Fe
            )  # Remaining elements are set to Fe

            spec = sess.return_spectrum(kT) / binsize / 1e-14

            return K * spec



if is_py3:


    _abs_tables = {
        "phabs": {"AG89": "angr", "ASPL": "aspl"},
        "tbabs": {"AG89": "angr", "ASPL": "aspl", "WILM": "wilm"},
        "wabs": {"AG89": "angr"},
    }
    _abund_info = {}
    _abund_info[
        "WILMS"
    ] = "wilms\nfrom Wilms, Allen & McCray (2000), ApJ 542, 914 \n except for elements not listed which are given zero abundance)\n https://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html "
    _abund_info[
        "AG89"
    ] = "angr\nfrom Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197)\n https://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html"
    _abund_info[
        "ASPL"
    ] = "aspl\nfrom Asplund M., Grevesse N., Sauval A.J. & Scott P. (2009, ARAA, 47, 481)\nhttps://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html"


    def _get_xsect_table(model, abund_table):
        """
        contructs the abundance table from the values given
        """

        assert model in _abs_tables, "the model %s does not exist" % model
        assert abund_table in _abs_tables[model], (
            "the table %s does not exist" % abund_table
        )

        path_to_xsect = _get_data_file_path(
            os.path.join(
                "xsect", "xsect_%s_%s.fits" % (model, _abs_tables[model][abund_table])
            )
        )

        fxs = fits.open(path_to_xsect)
        dxs = fxs[1].data
        xsect_ene = dxs["ENERGY"]
        xsect_val = dxs["SIGMA"]

        return np.array(xsect_ene), np.array(xsect_val)




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
                fixed: True


        """

        def _setup(self):
            self._fixed_units = (astropy_units.keV, astropy_units.dimensionless_unscaled)
            self.init_xsect()

        def _set_units(self, x_unit, y_unit):
            self.NH.unit = astropy_units.cm ** (-2)
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
                self.xsect_ene, self.xsect_val = _get_xsect_table("phabs", abund_table)
                self._abund_table = abund_table

            except:

                print("defaulting to AG89")
                self.xsect_ene, self.xsect_val = _get_xsect_table("phabs", abund_table)

                self._abund_table = "AG89"

        @cache_array_method()
        def _cached_interp(self, x):

            return np.interp(x, self.xsect_ene, self.xsect_val)

        def evaluate(self, x, NH, redshift):

            if isinstance(x, astropy_units.Quantity):

                _unit = astropy_units.cm ** 2
                _y_unit = astropy_units.dimensionless_unscaled
                _x = x.value

            else:

                _unit = 1.0
                _y_unit = 1.0

                _x = x

            xsect_interp = self._cached_interp(_x * (1 + redshift))

            spec = np.exp(-NH * xsect_interp * _unit) * _y_unit

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

            redshift :
                desc : the redshift of the source
                initial value : 0.
                is_normalization : False
                min : 0
                max : 15
                delta : 0.1
                fixed: True


        """

        def _setup(self):

            self.init_xsect()

            self._fixed_units = (astropy_units.keV, astropy_units.dimensionless_unscaled)

        def _set_units(self, x_unit, y_unit):
            self.NH.unit = astropy_units.cm ** (-2)
            self.redshift.unit = astropy_units.dimensionless_unscaled

        def init_xsect(self, abund_table="WILM"):
            """
            Set the abundance table

            :param abund_table: "WILM", "ASPL", "AG89" 
            :returns: 
            :rtype: 

            """

            try:
                self.xsect_ene, self.xsect_val = _get_xsect_table("tbabs", abund_table)
                self._abund_table = abund_table

            except:

                print("defaulting to WILM")
                self.xsect_ene, self.xsect_val = _get_xsect_table("tbabs", abund_table)

                self._abund_table = "WILM"

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

            else:

                _unit = 1.0
                _y_unit = 1.0

                _x = x

            xsect_interp = self._cached_interp(_x * (1 + redshift))

            spec = np.exp(-NH * xsect_interp * _unit) * _y_unit

            return spec


    # WAbs class
    @six.add_metaclass(FunctionMeta)
    class WAbs(Function1D):
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
                fixed: True


        """

        def _setup(self):
            self._fixed_units = (astropy_units.keV, astropy_units.dimensionless_unscaled)
            self.init_xsect()

        def _set_units(self, x_unit, y_unit):
            self.NH.unit = astropy_units.cm ** (-2)
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

            else:

                _unit = 1.0
                _y_unit = 1.0

                _x = x

            xsect_interp = self._cached_interp(_x * (1 + redshift))

            spec = np.exp(-NH * xsect_interp * _unit) * _y_unit

            return spec
