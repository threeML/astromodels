from __future__ import division

import astropy.units as astropy_units
import numpy as np
from past.utils import old_div

import astromodels.functions.numba_functions as nb_func
from astromodels.core.units import get_units
from astromodels.functions.function import Function1D, FunctionMeta
from astromodels.utils.configuration import astromodels_config
from astromodels.utils.logging import setup_logger

log = setup_logger(__name__)

__author__ = "giacomov"
# DMFitFunction and DMSpectra add by Andrea Albert (aalbert@slac.stanford.edu) Oct 26, 2016

erg2keV = 6.24151e8


class GSLNotAvailable(ImportWarning):
    pass


class NaimaNotAvailable(ImportWarning):
    pass



class InvalidUsageForFunction(Exception):
    pass


# Now let's try and import optional dependencies

import astropy.units as u

try:

    # Naima is for numerical computation of Synch. and Inverse compton spectra in randomly oriented
    # magnetic fields


    import naima

except ImportError:

    if astromodels_config.logging.startup_warnings:

        log.warning(
            "The naima package is not available. Models that depend on it will not be available"
        )

    has_naima = False

else:

    has_naima = True

try:

    # GSL is the GNU Scientific Library. Pygsl is the python wrapper for it. It is used by some
    # functions for faster computation

    from pygsl.testing.sf import gamma_inc

except ImportError:

    if astromodels_config.logging.startup_warnings:

        log.warning(
            "The GSL library or the pygsl wrapper cannot be loaded. Models that depend on it will not be "
            "available."
        )

    has_gsl = False

else:

    has_gsl = True



class StepFunction(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A function which is constant on the interval lower_bound - upper_bound and 0 outside the interval. The
        extremes of the interval are counted as part of the interval.

    latex : $ f(x)=\begin{cases}0 & x < \text{lower_bound} \\\text{value} & \text{lower_bound} \le x \le \text{upper_bound} \\ 0 & x > \text{upper_bound} \end{cases}$

    parameters :

        lower_bound :

            desc : Lower bound for the interval
            initial value : 0

        upper_bound :

            desc : Upper bound for the interval
            initial value : 1

        value :

            desc : Value in the interval
            initial value : 1.0

    tests :
        - { x : 0.5, function value: 1.0, tolerance: 1e-20}
        - { x : -0.5, function value: 0, tolerance: 1e-20}

    """

    def _set_units(self, x_unit, y_unit):
        # Lower and upper bound has the same unit as x
        self.lower_bound.unit = x_unit
        self.upper_bound.unit = x_unit

        # value has the same unit as y
        self.value.unit = y_unit

    def evaluate(self, x, lower_bound, upper_bound, value):
        # The value * 0 is to keep the units right

        result = np.zeros(x.shape) * value * 0

        idx = (x >= lower_bound) & (x <= upper_bound)
        result[idx] = value

        return result


class StepFunctionUpper(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A function which is constant on the interval lower_bound - upper_bound and 0 outside the interval. The
        upper interval is open.

    latex : $ f(x)=\begin{cases}0 & x < \text{lower_bound} \\\text{value} & \text{lower_bound} \le x \le \text{upper_bound} \\ 0 & x > \text{upper_bound} \end{cases}$

    parameters :

        lower_bound :

            desc : Lower bound for the interval
            initial value : 0
            fix : yes

        upper_bound :

            desc : Upper bound for the interval
            initial value : 1
            fix : yes

        value :

            desc : Value in the interval
            initial value : 1.0

    tests :
        - { x : 0.5, function value: 1.0, tolerance: 1e-20}
        - { x : -0.5, function value: 0, tolerance: 1e-20}

    """

    def _set_units(self, x_unit, y_unit):
        # Lower and upper bound has the same unit as x
        self.lower_bound.unit = x_unit
        self.upper_bound.unit = x_unit

        # value has the same unit as y
        self.value.unit = y_unit

    def evaluate(self, x, lower_bound, upper_bound, value):
        # The value * 0 is to keep the units right

        result = np.zeros(x.shape) * value * 0

        idx = (x >= lower_bound) & (x < upper_bound)
        result[idx] = value

        return result


# noinspection PyPep8Naming



    

# noinspection PyPep8Naming


class Sin(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A sinusodial function

    latex : $ K~\sin{(2\pi f x + \phi)} $

    parameters :

        K :

            desc : Normalization
            initial value : 1
            is_normalization : True

        f :

            desc : frequency
            initial value : 1.0 / (2 * np.pi)
            min : 0

        phi :

            desc : phase
            initial value : 0
            min : -np.pi
            max : +np.pi
            unit: rad

    tests :
        - { x : 0.0, function value: 0.0, tolerance: 1e-10}
        - { x : 1.5707963267948966, function value: 1.0, tolerance: 1e-10}

    """

    def _set_units(self, x_unit, y_unit):
        # The normalization has the same unit of y
        self.K.unit = y_unit

        # The unit of f is 1 / [x] because fx must be a pure number. However,
        # np.pi of course doesn't have units, so we add a rad
        self.f.unit = x_unit ** (-1) * astropy_units.rad

        # The unit of phi is always the same (radians)

        self.phi.unit = astropy_units.rad

    # noinspection PyPep8Naming
    def evaluate(self, x, K, f, phi):
        return K * np.sin(2 * np.pi * f * x + phi)




class DiracDelta(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        return  at zero_point

    latex : $ value $

    parameters :

        value :

            desc : Constant value
            initial value : 0

        zero_point:

             desc: value at which function is non-zero
             initial value : 0
             fix : yes


    """

    def _set_units(self, x_unit, y_unit):

        self.value.unit = y_unit
        self.zero_point.unit = x_unit

    def evaluate(self, x, value, zero_point):

        out = np.zeros(x.shape) * value * 0

        out[x == zero_point] = value

        return out


if has_naima:

    class Synchrotron(Function1D, metaclass=FunctionMeta):
        r"""
        description :
            Synchrotron spectrum from an input particle distribution, using Naima (naima.readthedocs.org)
        latex: not available
        parameters :
            B :
                desc : magnetic field
                initial value : 3.24e-6
                unit: Gauss
            distance :
                desc : distance of the source
                initial value : 1.0
                unit : kpc
            emin :
                desc : minimum energy for the particle distribution
                initial value : 1
                fix : yes
                unit: GeV
            emax :
                desc : maximum energy for the particle distribution
                initial value : 510e3
                fix : yes
                unit: GeV
            need:
                desc: number of points per decade in which to evaluate the function
                initial value : 10
                min : 2
                max : 100
                fix : yes
        """

        def _setup(self):

            # this is needed to load from
            # yaml

            for k,v in self._children.items():

                if k not in self._parameters:

                    # ok, this is probably a
                    #  particle distribution
                    self.set_particle_distribution(v)

                    break


        
        def _set_units(self, x_unit, y_unit):

            # This function can only be used as a spectrum,
            # so let's check that x_unit is a energy and y_unit is
            # differential flux

            if hasattr(x_unit, "physical_type") and x_unit.physical_type == "energy":

                # Now check that y is a differential flux
                current_units = get_units()
                should_be_unitless = y_unit * (
                    current_units.energy * current_units.time * current_units.area
                )

                if (
                    not hasattr(should_be_unitless, "physical_type")
                    or should_be_unitless.decompose().physical_type != "dimensionless"
                ):
                    # y is not a differential flux
                    raise InvalidUsageForFunction(
                        "Unit for y is not differential flux. The function synchrotron "
                        "can only be used as a spectrum."
                    )
            else:

                raise InvalidUsageForFunction(
                    "Unit for x is not an energy. The function synchrotron can only be used "
                    "as a spectrum"
                )

                # we actually don't need to do anything as the units are already set up

        def set_particle_distribution(self, function):

            if "particle_distribution" in self._external_functions:

                self.unlink_external_function("particle_distribution")

            self.link_external_function(function, "particle_distribution")

            self._particle_distribution_name = function.name

            self._particle_distribution = function


            # Now set the units for the function

            current_units = get_units()

            self._particle_distribution.set_units(
                current_units.energy, current_units.energy ** (-1)
            )

            # Naima wants a function which accepts a quantity as x (in units of eV) and returns an astropy quantity,
            # so we need to create a wrapper which will remove the unit from x and add the unit to the return
            # value

            self._particle_distribution_wrapper = lambda x: old_div(
                function(x.value), current_units.energy
            )


        def get_particle_distribution(self):

            return self._external_functions["particle_distribution"]
 

        particle_distribution = property(
            get_particle_distribution,
            set_particle_distribution,
            doc="""Get/set particle distribution for electrons""",
        )

        def fix_units(self, x, B, distance, emin, emax):

            if isinstance(x, u.Quantity):
                return (
                    True,
                    x.to(get_units().energy),
                    B.to(u.Gauss),
                    distance.to(u.kpc),
                    emin.to(u.GeV),
                    emax.to(u.GeV),
                )
            else:
                return (
                    False,
                    x * (get_units().energy),
                    B * (u.Gauss),
                    distance * (u.kpc),
                    emin * (u.GeV),
                    emax * (u.GeV),
                )

        # noinspection PyPep8Naming
        def evaluate(self, x, B, distance, emin, emax, need):

            has_units, x, B, distance, emin, emax = self.fix_units(
                x, B, distance, emin, emax
            )

            _synch = naima.models.Synchrotron(
                self._particle_distribution_wrapper,
                B,
                Eemin=emin,
                Eemax=emax,
                nEed=need,
            )

            if has_units:
                return _synch.flux(x, distance=distance)
            else:
                return _synch.flux(x, distance=distance).value

        def to_dict(self, minimal=False):

            data = super(Function1D, self).to_dict(minimal)

            # if not minimal:
            #     data["extra_setup"] = {
            #         "particle_distribution": self.particle_distribution.path
            #     }

            return data


class _ComplexTestFunction(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        A useless function to be used during automatic tests

    latex: not available

    parameters :
        A :
            desc : none
            initial value : 3.24e-6
            min : 1e-6
            max : 1e-5

        B :
            desc : none
            initial value : -10
            min : -100
            max : 100
            delta : 0.1
    
    properties:
        file_name: 
            desc: a file name
            defer: True 
        dummy:
            desc: a dummy property
            initial value: test
            allowed values:
                - test
                - love

    """

    
    def _set_units(self, x_unit, y_unit):

        self.A.unit = y_unit
        self.B.unit = old_div(y_unit, x_unit)

    def set_particle_distribution(self, function):

        if "particle_distribution" in self._external_functions:

            self.unlink_external_function("particle_distribution")

        self.link_external_function(function, "particle_distribution")

        self._particle_distribution_name = function.name
        
        self._particle_distribution = function

    def get_particle_distribution(self):

        return self._external_functions["particle_distribution"]
        

    particle_distribution = property(
        get_particle_distribution,
        set_particle_distribution,
        doc="""Get/set particle distribution for electrons""",
    )

    # noinspection PyPep8Naming
    def evaluate(self, x, A, B):

        return A + B * x

    def to_dict(self, minimal=False):

        data = super(Function1D, self).to_dict(minimal)

        # if not minimal:

        #     data["extra_setup"] = {
        #         "particle_distribution": self.particle_distribution.name
        #     }

        return data


class Log_parabola(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A log-parabolic function. NOTE that we use the high-energy convention of using the natural log in place of the
        base-10 logarithm. This means that beta is a factor 1 / log10(e) larger than what returned by those software
        using the other convention.

    latex : $ K \left( \frac{x}{piv} \right)^{\alpha -\beta \log{\left( \frac{x}{piv} \right)}} $

    parameters :

        K :

            desc : Normalization
            initial value : 1.0
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e5

        piv :
            desc : Pivot (keep this fixed)
            initial value : 1
            fix : yes

        alpha :

            desc : index
            initial value : -2.0

        beta :

            desc : curvature (positive is concave, negative is convex)
            initial value : 1.0

    """

    def _set_units(self, x_unit, y_unit):

        # K has units of y

        self.K.unit = y_unit

        # piv has the same dimension as x
        self.piv.unit = x_unit

        # alpha and beta are dimensionless
        self.alpha.unit = astropy_units.dimensionless_unscaled
        self.beta.unit = astropy_units.dimensionless_unscaled

    def evaluate(self, x, K, piv, alpha, beta):

        # print("Receiving %s" % ([K, piv, alpha, beta]))

        xx = np.divide(x, piv)

        try:

            return K * xx ** (alpha - beta * np.log(xx))

        except ValueError:

            # The current version of astropy (1.1.x) has a bug for which quantities that have become
            # dimensionless because of a division (like xx here) are not recognized as such by the power
            # operator, which throws an exception: ValueError: Quantities and Units may only be raised to a scalar power
            # This is a quick fix, waiting for astropy 1.2 which will fix this

            xx = xx.to("")

            return K * xx ** (alpha - beta * np.log(xx))

    @property
    def peak_energy(self):
        """
        Returns the peak energy in the nuFnu spectrum

        :return: peak energy in keV
        """

        # Eq. 6 in Massaro et al. 2004
        # (http://adsabs.harvard.edu/abs/2004A%26A...413..489M)

        return self.piv.value * pow(
            10, old_div(((2 + self.alpha.value) * np.log(10)),
                        (2 * self.beta.value))
        )


if has_gsl:

    class Cutoff_powerlaw_flux(Function1D, metaclass=FunctionMeta):
        r"""
        description :

            A cutoff power law having the flux as normalization, which should reduce the correlation among
            parameters.

        latex : $ \frac{F}{T(b)-T(a)} ~x^{index}~\exp{(-x/x_{c})}~\text{with}~T(x)=-x_{c}^{index+1} \Gamma(index+1, x/C)~\text{(}\Gamma\text{ is the incomplete gamma function)} $

        parameters :

            F :

                desc : Integral between a and b
                initial value : 1e-5
                is_normalization : True

            index :

                desc : photon index
                initial value : -2.0

            xc :

                desc : cutoff position
                initial value : 50.0

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
            # K has units of y * x
            self.F.unit = y_unit * x_unit

            # alpha is dimensionless
            self.index.unit = astropy_units.dimensionless_unscaled

            # xc, a and b have the same dimension as x
            self.xc.unit = x_unit
            self.a.unit = x_unit
            self.b.unit = x_unit

        @staticmethod
        def _integral(a, b, index, ec):
            ap1 = index + 1

            def integrand(x):
                return -pow(ec, ap1) * gamma_inc(ap1, old_div(x, ec))

            return integrand(b) - integrand(a)

        def evaluate(self, x, F, index, xc, a, b):
            this_integral = self._integral(a, b, index, xc)

            return (
                F / this_integral *
                np.power(x, index) * np.exp(-1 * np.divide(x, xc))
            )


class Exponential_cutoff(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        An exponential cutoff

    latex : $ K \exp{(-x/xc)} $

    parameters :

        K :

            desc : Normalization
            initial value : 1.0
            fix : no
            is_normalization : True

        xc :
            desc : cutoff
            initial value : 100
            min : 1
    """

    def _set_units(self, x_unit, y_unit):
        # K has units of y

        self.K.unit = y_unit

        # piv has the same dimension as x
        self.xc.unit = x_unit

    def evaluate(self, x, K, xc):
        return K * np.exp(np.divide(x, -xc))
