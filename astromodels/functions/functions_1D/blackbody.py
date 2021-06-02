from __future__ import division

import astropy.units as astropy_units
import numpy as np
from past.utils import old_div

import astromodels.functions.numba_functions as nb_func
from astromodels.core.units import get_units
from astromodels.functions.function import (Function1D, FunctionMeta,
                                            ModelAssertionViolation)


class Blackbody(Function1D, metaclass=FunctionMeta):
    r"""

    description :
        A blackbody function

    latex : $f(x) = K \frac{x^2}{\exp(\frac{x}{kT}) -1}  $

    parameters :
        K :
            desc :
            initial value : 1e-4
            min : 0.
            is_normalization : True

        kT :
            desc : temperature of the blackbody
            initial value : 30.0
            min: 0.

    """

    def _set_units(self, x_unit, y_unit):
        # The normalization has the same units as y
        self.K.unit = old_div(y_unit, (x_unit ** 2))

        # The break point has always the same dimension as the x variable
        self.kT.unit = x_unit

    def evaluate(self, x, K, kT):

        if isinstance(x, astropy_units.Quantity):

            K_ = K.value
            kT_ = kT.value

            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, kT_, x_, = (
                K,
                kT,
                x,
            )

        result = nb_func.bb_eval(x_, K_, kT_)

        return result * unit_


class ModifiedBlackbody(Function1D, metaclass=FunctionMeta):
    r"""

    description :
        A blackbody function

    latex : $f(x) = K \frac{x^2}{\exp(\frac{x}{kT}) -1}  $

    parameters :
        K :
            desc :
            initial value : 1e-4
            min : 0.
            is_normalization : True

        kT :
            desc : temperature of the blackbody
            initial value : 30.0
            min: 0.

    """

    def _set_units(self, x_unit, y_unit):
        # The normalization has the same units as y
        self.K.unit = y_unit / x_unit

        # The break point has always the same dimension as the x variable
        self.kT.unit = x_unit

    def evaluate(self, x, K, kT):

        if isinstance(x, astropy_units.Quantity):

            K_ = K.value
            kT_ = kT.value

            x_ = x.value

            unit_ = self.y_unit

        else:

            unit_ = 1.0
            K_, kT_, x_, = (
                K,
                kT,
                x,
            )

        result = nb_func.mbb_eval(x_, K_, kT_)

        return result * unit_


class NonDissipativePhotosphere(Function1D, metaclass=FunctionMeta):
    r"""

    description :
        Non-dissipative photosphere of a GRB occuring above the saturation radius
        Acuner, Z., Ryde, F. & Yu, H.-F. Mon Not R Astron Soc 487, 5508–5519 (2019).
    latex : $N_{\mathrm{E}}=K\left(\frac{E}{E_{\mathrm{pivot}}}\right)^{0.4} e^{-\left(\frac{E}{E_{c}}\right)^{0.65}}$

    parameters :
        K :
            desc :
            initial value : 1e-4
            min : 0.
            is_normalization : True

        ec :
            desc : peak energy
            initial value : 200.0
            min: 0.
       
        piv :
            desc : the pivot energy
            initial value: 100. 
            fix: True


    """

    def _set_units(self, x_unit, y_unit):
        # The normalization has the same units as y
        self.K.unit = y_unit

        # The break point has always the same dimension as the x variable
        self.ec.unit = x_unit

        self.piv.unit = x_unit
        
    def evaluate(self, x, K, ec, piv):

        if isinstance(x, astropy_units.Quantity):

            K_ = K.value
            ec_ = ec.value
            piv_ = piv.value

            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, ec_, x_, piv_ = (
                K,
                ec,
                x,
                piv
            )

        result = nb_func.non_diss_photoshere_generic(x_, K_, ec_, piv_, 0.4, 0.65)

        return result * unit_

    
class NonDissipativePhotosphere_Deep(Function1D, metaclass=FunctionMeta):
    r"""

    description :
        Non-dissipative photosphere of a GRB occuring BELOW the saturation radius.
        Acuner, Z., Ryde, F. & Yu, H.-F. Mon Not R Astron Soc 487, 5508–5519 (2019).
   
    
    latex : $N_{\mathrm{E}}=K\left(\frac{E}{E_{\mathrm{pivot}}}\right)^{0.66} e^{-\left(\frac{E}{E_{c}}\right)}$

    parameters :
        K :
            desc :
            initial value : 1e-4
            min : 0.
            is_normalization : True

        ec :
            desc : peak energy
            initial value : 200.0
            min: 0.
       
        piv :
            desc : the pivot energy
            initial value: 100. 
            fix: True


    """

    def _set_units(self, x_unit, y_unit):
        # The normalization has the same units as y
        self.K.unit = y_unit

        # The break point has always the same dimension as the x variable
        self.ec.unit = x_unit

        self.piv.unit = x_unit
        
    def evaluate(self, x, K, ec, piv):

        if isinstance(x, astropy_units.Quantity):

            K_ = K.value
            ec_ = ec.value
            piv_ = piv.value

            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, ec_, x_, piv_ = (
                K,
                ec,
                x,
                piv
            )

        result = nb_func.non_diss_photoshere_generic(x_, K_, ec_, piv_, 0.66, 1.)

        return result * unit_

    
