import os
import sys
from functools import lru_cache, wraps

import astropy.units as astropy_units
import numpy as np
import six
from astropy.io import fits

from astromodels.functions.function import Function1D, FunctionMeta
from astromodels.utils import configuration
from astromodels.utils import _get_data_file_path
import gc

try:

    import pyatomdb

    has_atomdb = True

except:

    has_atomdb = False

if has_atomdb:
    # APEC class
    
    class APEC(Function1D, metaclass=FunctionMeta):
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

        def clean(self):
            """
            Clean the current APEC session to avoid having too many open files
            :returns: 
            """
            
            self.session = None
            del self.session
            gc.collect()

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
    
    class VAPEC(Function1D, metaclass=FunctionMeta):
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

        def clean(self):
            """
            Clean the current APEC session to avoid having too many open files
            :returns: 
            """
            
            self.session = None
            del self.session
            gc.collect()

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
                [6, ], C,
            )

            sess.set_abund(
                [7, ], N,
            )

            sess.set_abund(
                [8, ], O,
            )

            sess.set_abund(
                [10, ], Ne,
            )

            sess.set_abund(
                [12, ], Mg,
            )

            sess.set_abund(
                [13, ], Al,
            )

            sess.set_abund(
                [14, ], Si,
            )

            sess.set_abund(
                [16, ], S,
            )

            sess.set_abund(
                [18, ], Ar,
            )

            sess.set_abund(
                [20, ], Ca,
            )

            sess.set_abund(
                [26, ], Fe,
            )

            sess.set_abund(
                [28, ], Ni,
            )

            sess.set_abund(
                [9, 11, 15, 17, 19, 21, 22, 23, 24, 25, 27, 29, 30], Fe
            )  # Remaining elements are set to Fe

            spec = sess.return_spectrum(kT) / binsize / 1e-14

            return K * spec


