import math
import warnings
from typing import Iterable

import astropy.units as astropy_units
import numpy as np
import six
from past.utils import old_div
from scipy.special import erfcinv, gamma, gammaincc

import astromodels.functions.numba_functions as nb_func
from astromodels.core.units import get_units
from astromodels.functions.function import (Function1D, FunctionMeta,
                                            ModelAssertionViolation)

try:
    from threeML.config.config import threeML_config

    _has_threeml = True

except ImportError:

    _has_threeml = False


from astromodels.utils.logging import setup_logger

log = setup_logger(__name__)

__author__ = 'giacomov'
# DMFitFunction and DMSpectra add by Andrea Albert (aalbert@slac.stanford.edu) Oct 26, 2016

erg2keV = 6.24151e8


class Powerlaw(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A simple power-law

    latex : $ K~\frac{x}{piv}^{index} $

    parameters :

        K :

            desc : Normalization (differential flux at the pivot value)
            initial value : 1.0
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e3
            delta : 0.1

        piv :

            desc : Pivot value
            initial value : 1
            fix : yes

        index :

            desc : Photon index
            initial value : -2.01
            min : -10
            max : 10

    tests :
        - { x : 10, function value: 0.01, tolerance: 1e-20}
        - { x : 100, function value: 0.0001, tolerance: 1e-20}

    """

    def _set_units(self, x_unit, y_unit):
        # The index is always dimensionless
        self.index.unit = astropy_units.dimensionless_unscaled

        # The pivot energy has always the same dimension as the x variable
        self.piv.unit = x_unit

        # The normalization has the same units as the y

        self.K.unit = y_unit

    # noinspection PyPep8Naming
    def evaluate(self, x, K, piv, index):

        if isinstance(x, astropy_units.Quantity):
            index_ = index.value
            K_ = K.value
            piv_ = piv.value
            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, piv_, x_, index_ = K, piv, x, index

        result = nb_func.plaw_eval(x_, K_, index_, piv_)

        return result * unit_


# noinspection PyPep8Naming


class Powerlaw_flux(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A simple power-law with the photon flux in a band used as normalization. This will reduce the correlation
        between the index and the normalization.

    latex : $ \frac{F(\gamma+1)} {b^{\gamma+1} - a^{\gamma+1}} (x)^{\gamma}$

    parameters :

        F :

            desc : Integral between a and b
            initial value : 1
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e3
            delta : 0.1

        index :

            desc : Photon index
            initial value : -2
            min : -10
            max : 10

        a :

            desc : lower bound for the band in which computing the integral F
            initial value : 1.0
            fix : yes

        b :

            desc : upper bound for the band in which computing the integral F
            initial value : 100.0
            fix : yes


    """

    def _set_units(self, x_unit, y_unit):
        # The flux is the integral over x, so:
        self.F.unit = y_unit * x_unit

        # The index is always dimensionless
        self.index.unit = astropy_units.dimensionless_unscaled

        # a and b have the same units as x

        self.a.unit = x_unit
        self.b.unit = x_unit

    # noinspection PyPep8Naming
    def evaluate(self, x, F, index, a, b):

        if isinstance(x, astropy_units.Quantity):
            index_ = index.value
            F_ = F.value
            a_ = a.value
            b_ = b.value
            x_ = x.value

            yunit_ = self.y_unit
            xunit_ = self.x_unit

        else:
            yunit_ = 1.0
            xunit_ = 1.0
            F_,  x_, index_, a_, b_ = F, x, index, a, b

        gp1 = index_ + 1

        norm = F_ * gp1 / (((b_)**gp1 - (a_)**gp1))

        return nb_func.plaw_eval(x_, norm, index_, 1.) * yunit_


class Powerlaw_Eflux(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        A  power-law where the normalization is the energy flux defined between a and b
    latex : $ F~\frac{x}{piv}^{index} $
    parameters :
        F :
            desc : Normalization (energy flux at the between a and b) erg /cm2 s
            initial value : 1.e-5
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e3
            delta : 0.1
        piv :
            desc : Pivot value
            initial value : 1
            fix : yes
        index :
            desc : Photon index
            initial value : -2
            min : -10
            max : 10


        a :
            desc : lower energy integral bound (keV)
            initial value : 1
            min : 0
            fix: yes

        b :
            desc : upper energy integral bound (keV)
            initial value : 100
            min : 0
            fix: yes

    """

    def _set_units(self, x_unit, y_unit):
        # The index is always dimensionless
        self.index.unit = astropy_units.dimensionless_unscaled

        # The pivot energy has always the same dimension as the x variable
        self.piv.unit = x_unit

        # a and b have the same units of x
        self.a.unit = x_unit
        self.b.unit = x_unit

        # The normalization has the same units as the y

        self.F.unit = y_unit * x_unit

    # noinspection PyPep8Naming
    def evaluate(self, x, F, piv, index, a, b):

        if isinstance(x, astropy_units.Quantity):
            index_ = index.value
            F_ = F.value
            piv_ = piv.value
            x_ = x.value
            a_ = a.value
            b_ = b.value

            yunit_ = self.y_unit
            xunit_ = self.x_unit

        else:
            yunit_ = 1.0
            xunit_ = 1.0
            F_, piv_, x_, index_, a_, b_ = F, piv, x, index, a, b

        intflux = nb_func.plaw_flux_norm(index_, a_, b_)

        norm = (F_ / (intflux)) * erg2keV

        out = nb_func.plaw_eval(x_, 1., index_, piv_)

        return norm * out * yunit_


class Cutoff_powerlaw(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A power law multiplied by an exponential cutoff

    latex : $ K~\frac{x}{piv}^{index}~\exp{-x/xc} $

    parameters :

        K :

            desc : Normalization (differential flux at the pivot value)
            initial value : 1.0
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e3
            delta : 0.1

        piv :

            desc : Pivot value
            initial value : 1
            fix : yes

        index :

            desc : Photon index
            initial value : -2
            min : -10
            max : 10

        xc :

            desc : Cutoff energy
            initial value : 10.0
            transformation : log10
            min: 1.0

    """

    def _set_units(self, x_unit, y_unit):
        # The index is always dimensionless
        self.index.unit = astropy_units.dimensionless_unscaled

        # The pivot energy has always the same dimension as the x variable
        self.piv.unit = x_unit

        self.xc.unit = x_unit

        # The normalization has the same units as the y

        self.K.unit = y_unit

    # noinspectionq PyPep8Naming
    def evaluate(self, x, K, piv, index, xc):

        if isinstance(x, astropy_units.Quantity):
            index_ = index.value
            K_ = K.value
            piv_ = piv.value
            xc_ = xc.value
            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, piv_, x_, index_, xc_ = K, piv, x, index, xc

        result = nb_func.cplaw_eval(x_, K_, xc_, index_, piv_)

        return result * unit_


class Cutoff_powerlaw_Ep(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A power law multiplied by an exponential cutoff parametrized with Ep

    latex : $ K~\frac{x}{piv}^{index}~\exp{-x(2+index)/xp} $

    parameters :

        K :

            desc : Normalization (differential flux at the pivot value)
            initial value : 1.0
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e3
            delta : 0.1

        piv :

            desc : Pivot value
            initial value : 1
            fix : yes

        index :

            desc : Photon index
            initial value : -2
            min : -10
            max : 10

        xp :

            desc : peak in the x * x * N (nuFnu if x is a energy)
            initial value : 500
            min : 10
            max : 1e4
            transformation : log10

    """

    def _set_units(self, x_unit, y_unit):
        # The index is always dimensionless
        self.index.unit = astropy_units.dimensionless_unscaled

        # The pivot energy has always the same dimension as the x variable
        self.piv.unit = x_unit

        self.xp.unit = x_unit

        # The normalization has the same units as the y

        self.K.unit = y_unit

    # noinspectionq PyPep8Naming
    def evaluate(self, x, K, piv, index, xp):

        if isinstance(x, astropy_units.Quantity):
            index_ = index.value
            K_ = K.value
            piv_ = piv.value
            xp_ = xp.value
            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, piv_, x_, index_, xp_ = K, piv, x, index, xp

        xc = xp / (2+index)

        result = nb_func.cplaw_eval(x_, K_, xc, index_, piv_)

        return result * unit_


class Inverse_cutoff_powerlaw(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        A power law multiplied by an exponential cutoff [Note: instead of cutoff energy energy parameter xc, b = 1/xc is used]
    latex : $ K~\frac{x}{piv}^{index}~\exp{-x~\b} $
    parameters :
        K :
            desc : Normalization (differential flux at the pivot value)
            initial value : 1.0
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e3
            delta : 0.1
        piv :
            desc : Pivot value
            initial value : 1
            fix : yes
        index :
            desc : Photon index
            initial value : -2
            min : -10
            max : 10
        b :
            desc : inverse cutoff energy i.e 1/xc
            initial value : 1
    """

    def _set_units(self, x_unit, y_unit):
        # The index is always dimensionless
        self.index.unit = astropy_units.dimensionless_unscaled

        # The pivot energy has always the same dimension as the x variable
        self.piv.unit = x_unit

        self.b.unit = 1 / x_unit

        # The normalization has the same units as the y

        self.K.unit = y_unit

    # noinspectionq PyPep8Naming
    def evaluate(self, x, K, piv, index, b):

        if isinstance(x, astropy_units.Quantity):
            index_ = index.value
            K_ = K.value
            piv_ = piv.value
            b_ = b.value
            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, piv_, x_, index_, b_ = K, piv, x, index, b

        result = nb_func.cplaw_inverse_eval(x_, K_, b_, index_, piv_)

        return result * unit_


class Super_cutoff_powerlaw(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A power law with a super-exponential cutoff

    latex : $ K~\frac{x}{piv}^{index}~\exp{(-x/xc)^{\gamma}} $

    parameters :

        K :

            desc : Normalization (differential flux at the pivot value)
            initial value : 1.0
            is_normalization : True

        piv :

            desc : Pivot value
            initial value : 1
            fix : yes

        index :

            desc : Photon index
            initial value : -2
            min : -10
            max : 10

        xc :

            desc : Cutoff energy
            initial value : 10.0
            min : 1.0

        gamma :

            desc : Index of the super-exponential cutoff
            initial value : 1.0
            min : 0.1
            max : 10.0

    """

    def _set_units(self, x_unit, y_unit):
        # The index is always dimensionless
        self.index.unit = astropy_units.dimensionless_unscaled
        self.gamma.unit = astropy_units.dimensionless_unscaled

        # The pivot energy has always the same dimension as the x variable
        self.piv.unit = x_unit

        # The cutoff has the same dimensions as x
        self.xc.unit = x_unit

        # The normalization has the same units as the y

        self.K.unit = y_unit

    # noinspection PyPep8Naming
    def evaluate(self, x, K, piv, index, xc, gamma):

        if isinstance(x, astropy_units.Quantity):
            index_ = index.value
            K_ = K.value
            piv_ = piv.value
            xc_ = xc.value
            gamma_ = gamma.value
            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, piv_, x_, index_, xc_, gamma_ = K, piv, x, index, xc, gamma

        result = nb_func.super_cplaw_eval(x_, K_, piv_, index_, xc_, gamma_)

        return result * unit_


class SmoothlyBrokenPowerLaw(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A Smoothly Broken Power Law

    Latex : $  $

    parameters :

        K :

            desc : normalization
            initial value : 1
            min : 0
            is_normalization : True


        alpha :

            desc : power law index below the break
            initial value : -1
            min : -1.5
            max : 2

        break_energy:

            desc: location of the peak
            initial value : 300
            fix : no
            min : 10

        break_scale :

            desc: smoothness of the break
            initial value : 0.5
            min : 0.
            max : 10.
            fix : yes

        beta:

            desc : power law index above the break
            initial value : -2.
            min : -5.0
            max : -1.6

        pivot:

            desc: where the spectrum is normalized
            initial value : 100.
            fix: yes


    """

    def _set_units(self, x_unit, y_unit):

        # norm has same unit as energy
        self.K.unit = y_unit

        self.break_energy.unit = x_unit

        self.pivot.unit = x_unit

        self.alpha.unit = astropy_units.dimensionless_unscaled
        self.beta.unit = astropy_units.dimensionless_unscaled
        self.break_scale.unit = astropy_units.dimensionless_unscaled

    def evaluate(self, x, K, alpha, break_energy, break_scale, beta, pivot):

        if isinstance(x, astropy_units.Quantity):
            alpha_ = alpha.value
            beta_ = beta.value
            K_ = K.value
            pivot_ = pivot.value
            break_energy_ = break_energy.value
            break_scale_ = break_scale.value
            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, pivot_, x_, alpha_, beta_, break_scale_, break_energy_ = (
                K,
                pivot,
                x,
                alpha,
                beta,
                break_scale,
                break_energy,
            )

        result = nb_func.sbplaw_eval(
            x_, K_, alpha_, break_energy, break_scale_, beta_, pivot_
        )

        return result * unit_


class Broken_powerlaw(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A broken power law function

    latex : $ f(x)= K~\begin{cases}\left( \frac{x}{x_{b}} \right)^{\alpha} & x < x_{b} \\ \left( \frac{x}{x_{b}} \right)^{\beta} & x \ge x_{b} \end{cases} $

    parameters :

        K :

            desc : Normalization (differential flux at x_b)
            initial value : 1.0
            is_normalization : True

        xb :

            desc : Break point
            initial value : 10
            min : 1.0

        alpha :

            desc : Index before the break xb
            initial value : -1.5
            min : -10
            max : 10

        beta :

            desc : Index after the break xb
            initial value : -2.5
            min : -10
            max : 10

        piv :

            desc : Pivot energy
            initial value : 1.0
            fix : yes

    """

    def _set_units(self, x_unit, y_unit):
        # The normalization has the same units as y
        self.K.unit = y_unit

        # The break point has always the same dimension as the x variable
        self.xb.unit = x_unit

        # alpha and beta are dimensionless
        self.alpha.unit = astropy_units.dimensionless_unscaled
        self.beta.unit = astropy_units.dimensionless_unscaled

        self.piv.unit = x_unit

    # noinspection PyPep8Naming
    def evaluate(self, x, K, xb, alpha, beta, piv):
        # The K * 0 is to keep the units right. If the input has unit, this will make a result
        # array with the same units as K. If the input has no units, this will have no
        # effect whatsoever

        if isinstance(x, astropy_units.Quantity):
            alpha_ = alpha.value
            beta_ = alpha.value
            K_ = K.value
            xb_ = xb.value
            piv_ = piv.value
            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            alpha_, beta_, K_, piv_, x_, xb_ = alpha, beta, K, piv, x, xb

        result = nb_func.bplaw_eval(x_, K_, xb_, alpha_, beta_, piv_)

        return result * unit_


class Band(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        Band model from Band et al., 1993, parametrized with the peak energy

    latex : $ K \begin{cases} \left(\frac{x}{piv}\right)^{\alpha} \exp \left(-\frac{(2+\alpha) x}{x_{p}}\right) & x \leq (\alpha-\beta) \frac{x_{p}}{(\alpha+2)} \\ \left(\frac{x}{piv}\right)^{\beta} \exp (\beta-\alpha)\left[\frac{(\alpha-\beta) x_{p}}{piv(2+\alpha)}\right]^{\alpha-\beta} &x>(\alpha-\beta) \frac{x_{p}}{(\alpha+2)} \end{cases} $

    parameters :

        K :

            desc : Differential flux at the pivot energy
            initial value : 1e-4
            is_normalization : True

        alpha :

            desc : low-energy photon index
            initial value : -1.0
            min : -1.5
            max : 3

        xp :

            desc : peak in the x * x * N (nuFnu if x is a energy)
            initial value : 500
            min : 10

        beta :

            desc : high-energy photon index
            initial value : -2.0
            min : -5.0
            max : -1.6

        piv :

            desc : pivot energy
            initial value : 100.0
            fix : yes
    """

    def _set_units(self, x_unit, y_unit):
        # The normalization has the same units as y
        self.K.unit = y_unit

        # The break point has always the same dimension as the x variable
        self.xp.unit = x_unit

        self.piv.unit = x_unit

        # alpha and beta are dimensionless
        self.alpha.unit = astropy_units.dimensionless_unscaled
        self.beta.unit = astropy_units.dimensionless_unscaled

    def evaluate(self, x, K, alpha, xp, beta, piv):
        E0 = old_div(xp, (2 + alpha))

        if alpha < beta:
            raise ModelAssertionViolation("Alpha cannot be less than beta")

        if isinstance(x, astropy_units.Quantity):
            alpha_ = alpha.value
            beta_ = alpha.value
            K_ = K.value
            E0_ = E0.value
            piv_ = piv.value
            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            alpha_, beta_, K_, piv_, x_, E0_ = alpha, beta, K, piv, x, E0

        return nb_func.band_eval(x_, K_, alpha_, beta_, E0_, piv_) * unit_


class Band_grbm(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        Band model from Band et al., 1993, parametrized with the cutoff energy

    latex : $  $

    parameters :

        K :

            desc : Differential flux at the pivot energy
            initial value : 1e-4
            is_normalization : True

        alpha :

            desc : low-energy photon index
            initial value : -1.0
            min : -1.5
            max : 3

        xc :

            desc : cutoff of exp
            initial value : 500
            min : 10

        beta :

            desc : high-energy photon index
            initial value : -2.0
            min : -5.0
            max : -1.6

        piv :

            desc : pivot energy
            initial value : 100.0
            fix : yes
    """

    def _set_units(self, x_unit, y_unit):
        # The normalization has the same units as y
        self.K.unit = y_unit

        # The break point has always the same dimension as the x variable
        self.xc.unit = x_unit

        self.piv.unit = x_unit

        # alpha and beta are dimensionless
        self.alpha.unit = astropy_units.dimensionless_unscaled
        self.beta.unit = astropy_units.dimensionless_unscaled

    def evaluate(self, x, K, alpha, xc, beta, piv):

        if alpha < beta:
            raise ModelAssertionViolation("Alpha cannot be less than beta")

        idx = x < (alpha - beta) * xc

        # The K * 0 part is a trick so that out will have the right units (if the input
        # has units)

        out = np.zeros(x.shape) * K * 0

        out[idx] = (
            K * np.power(old_div(x[idx], piv), alpha) *
            np.exp(old_div(-x[idx], xc))
        )
        out[~idx] = (
            K
            * np.power((alpha - beta) * xc / piv, alpha - beta)
            * np.exp(beta - alpha)
            * np.power(old_div(x[~idx], piv), beta)
        )

        return out


class Band_Calderone(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        The Band model from Band et al. 1993, implemented however in a way which reduces the covariances between
        the parameters (Calderone et al., MNRAS, 448, 403C, 2015)

    latex : $ \text{(Calderone et al., MNRAS, 448, 403C, 2015)} $

    parameters :

        alpha :
            desc : The index for x smaller than the x peak
            initial value : -1
            min : -10
            max : 10

        beta :

            desc : index for x greater than the x peak (only if opt=1, i.e., for the
                   Band model)
            initial value : -2.2
            min : -7
            max : -1

        xp :

            desc : position of the peak in the x*x*f(x) space (if x is energy, this is the nuFnu or SED space)
            initial value : 200.0
            min : 0

        F :

            desc : integral in the band defined by a and b
            initial value : 1e-6
            is_normalization : True

        a:

            desc : lower limit of the band in which the integral will be computed
            initial value : 1.0
            min : 0
            fix : yes

        b:

            desc : upper limit of the band in which the integral will be computed
            initial value : 10000.0
            min : 0
            fix : yes

        opt :

            desc : option to select the spectral model (0 corresponds to a cutoff power law, 1 to the Band model)
            initial value : 1
            min : 0
            max : 1
            fix : yes

    """

    def _set_units(self, x_unit, y_unit):

        # alpha and beta are always unitless

        self.alpha.unit = astropy_units.dimensionless_unscaled
        self.beta.unit = astropy_units.dimensionless_unscaled

        # xp has the same dimension as x
        self.xp.unit = x_unit

        # F is the integral over x, so it has dimensions y_unit * x_unit
        self.F.unit = y_unit * x_unit

        # a and b have the same units of x
        self.a.unit = x_unit
        self.b.unit = x_unit

        # opt is just a flag, and has no units
        self.opt.unit = astropy_units.dimensionless_unscaled

    @staticmethod
    def ggrb_int_cpl(a, Ec, Emin, Emax):

        # Gammaincc does not support quantities
        i1 = gammaincc(2 + a, old_div(Emin, Ec)) * gamma(2 + a)
        i2 = gammaincc(2 + a, old_div(Emax, Ec)) * gamma(2 + a)

        return -Ec * Ec * (i2 - i1)

    def evaluate(self, x, alpha, beta, xp, F, a, b, opt):

        assert opt == 0 or opt == 1, "Opt must be either 0 or 1"

        if alpha < beta:
            raise ModelAssertionViolation("Alpha cannot be smaller than beta")

        if alpha < -2:
            raise ModelAssertionViolation("Alpha cannot be smaller than -2")

        # Cutoff energy

        if alpha == -2:

            Ec = old_div(xp, 0.0001)  # TRICK: avoid a=-2

        else:

            Ec = old_div(xp, (2 + alpha))

        # Split energy

        Esplit = (alpha - beta) * Ec

        # Evaluate model integrated flux and normalization

        if isinstance(alpha, astropy_units.Quantity):

            # The following functions do not allow the use of units
            alpha_ = alpha.value
            Ec_ = Ec.value
            a_ = a.value
            b_ = b.value
            Esplit_ = Esplit.value
            beta_ = beta.value
            x_ = x.value

            unit_ = self.x_unit

        else:

            alpha_, Ec_, a_, b_, Esplit_, beta_, x_ = alpha, Ec, a, b, Esplit, beta, x
            unit_ = 1.0

        if opt == 0:

            # Cutoff power law

            intflux = self.ggrb_int_cpl(alpha_, Ec_, a_, b_)

        else:

            # Band model

            if a <= Esplit and Esplit <= b:

                intflux = self.ggrb_int_cpl(
                    alpha_, Ec_, a_, Esplit_
                ) + nb_func.ggrb_int_pl(alpha_, beta_, Ec_, Esplit_, b_)

            else:

                if Esplit < a:

                    intflux = nb_func.ggrb_int_pl(alpha_, beta_, Ec_, a_, b_)

                else:

                    intflux = nb_func.ggrb_int_cpl(alpha_, Ec_, a_, b_)

        norm = F * erg2keV / (intflux * unit_)

        if opt == 0:

            # Cutoff power law

            flux = nb_func.cplaw_eval(x_, 1.0, Ec_, alpha_, Ec_)

            # flux = norm * np.power(old_div(x, Ec), alpha) * \
            #     np.exp(old_div(- x, Ec))

        else:

            # The norm * 0 is to keep the units right

            flux = nb_func.band_eval(x_, 1.0, alpha_, beta_, Ec_, Ec_)

        return norm * flux



class DoubleSmoothlyBrokenPowerlaw(Function1D, metaclass=FunctionMeta):
    r"""
    description : A smoothly broken power law with two breaks as parameterized in Ravasio, M. E. et al. Astron Astrophys 613, A16 (2018).

    latex : $\begin{array}{l}\begin{aligned}f(x)=& A x_{\mathrm{b}}^{\alpha_{1}}\left[\left[\left(\frac{x}{x_{\mathrm{b}}}\right)^{-\alpha_{1} n_{1}}+\left(\frac{x}{x_{\mathrm{b}}}\right)^{-\alpha_{2} n_{1}}\right]^{\frac{n_{2}}{n_{1}}}\right.\\&\left.+\left(\frac{x}{x_{\mathrm{j}}}\right)^{-\beta n_{2}} \cdot\left[\left(\frac{x_{\mathrm{j}}}{x_{\mathrm{b}}}\right)^{-\alpha_{1} n_{1}}+\left(\frac{x_{\mathrm{j}}}{x_{\mathrm{b}}}\right)^{-\alpha_{2} n_{1}}\right]^{\frac{n_{2}}{n_{1}}}\right]^{-\frac{1}{n_{2}}}\end{aligned}\\\text { where }\\x_{\mathrm{j}}=x_{\mathrm{p}} \cdot\left(-\frac{\alpha_{2}+2}{\beta+2}\right)^{\frac{1}{\left.\beta-\alpha_{2}\right) n_{2}}}\end{array}$

    parameters :

        K :

            desc : Differential flux at the pivot energy
            initial value : 1e-4
            is_normalization : True

        alpha1 :

            desc : photon index below xb
            initial value : -0.66

        xb :

            desc : break energy below xp
            initial value : 100
            min : 0

        n1 :

            desc : curvature of the first break
            initial value : 2.0
            min : 0
            fix: True

        alpha2 :

            desc : photon index between xb and xp
            initial value : -1.5

        xp :

            desc : nuFnu peak
            initial value : 300
            min : 0

        n2 :

            desc : curvature of the break at xp
            initial value : 2.0
            min : 0
            fix: True

        beta :

            desc : photon index above xp
            initial value : -2.5
            max : 2


        piv :

            desc : pivot energy
            initial value : 1.
            fix : yes
    """

    def _set_units(self, x_unit, y_unit):
        # The normalization has the same units as y
        self.K.unit = y_unit

        # The break point has always the same dimension as the x variable
        self.xp.unit = x_unit
        self.xb.unit = x_unit
        self.piv.unit = x_unit

        # alpha and beta are dimensionless
        self.alpha1.unit = astropy_units.dimensionless_unscaled
        self.alpha2.unit = astropy_units.dimensionless_unscaled                
        self.beta.unit = astropy_units.dimensionless_unscaled

        self.n1.unit = astropy_units.dimensionless_unscaled
        self.n2.unit = astropy_units.dimensionless_unscaled                
   

    def _fix_units(self, x, K, alpha1, xb, n1, alpha2, xp, n2, beta, piv):

            if isinstance(x, astropy_units.Quantity):

                return ( x.value,
                         K.value,
                         alpha1.value,
                         xb.value,
                         n1.value,
                         alpha2.value,
                         xp.value,
                         n2.value,
                         beta.value,
                         piv.value,
                         self.y_unit
                )
            
            else:

                return ( x,
                         K,
                         alpha1,
                         xb,
                         n1,
                         alpha2,
                         xp,
                         n2,
                         beta,
                         piv,
                         1.
                )


        
    def free_curvature(self) -> None:

        """
        free the two curvature parameters n1, n2

        :returns: 

        """
        self.n1.free = True
        self.n2.free = True


    def fix_curvature(self) -> None:

        """
        fix the two curvature parameters n1, n2

        :returns: 

        """
        self.n1.fix = True
        self.n2.fix = True

        
    def evaluate(self, x, K, alpha1, xb, n1, alpha2, xp, n2, beta, piv):

        x_, K_, alpha1_, xb_, n1_, alpha2_, xp_, n2_, beta_, piv_, y_unit = self._fix_units(x, K, alpha1, xb, n1, alpha2, xp, n2, beta, piv)

        return nb_func.dbl_sbpl(x_,
                        K_,
                        alpha1_,
                        alpha2_,
                        beta_,
                        xp_,
                        xb_,
                        n1_,
                        n2_,
                        piv_) * y_unit
        


    
