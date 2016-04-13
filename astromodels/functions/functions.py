__author__ = 'giacomov'

import math
import numpy as np
import warnings
from scipy.special import gammaincc, gamma

from astromodels.functions.function import Function
from astromodels.functions.function import FunctionMeta

# Now let's try and import optional dependencies

try:

    # Naima is for numerical computation of Synch. and Inverse compton spectra in randomly oriented
    # magnetic fields

    import naima
    import astropy.units as u

except ImportError:

    warnings.warn("The naima package is not available. Models that depend on it will not be available")

    has_naima = False

else:

    has_naima = True

try:

    # GSL is the GNU Scientific Library. Pygsl is the python wrapper for it. It is used by some
    # functions for faster computation

    from pygsl.testing.sf import gamma_inc

except ImportError:

    warnings.warn("The GSL library or the pygsl wrapper cannot be loaded. Models that depend on it will not be "
                  "available.")

    has_gsl = False

else:

    has_gsl = True

# noinspection PyPep8Naming
class powerlaw(Function):
    r"""
    description :

        A simple power-law with normalization expressed as
        a logarithm

    latex : $ \frac{dN}{dx} = 10^{logK}~\frac{x}{piv}^{index} $

    parameters :

        logK :

            desc : Logarithm of normalization
            initial value : 0
            min : -40
            max : 40
            unit : "1 / (keV cm2 s)"

        piv :

            desc : Pivot energy
            initial value : 1
            fix : yes
            unit: keV

        index :

            desc : Photon index
            initial value : -2
            min : -10
            max : 10

    tests :
        - { x : 10, function value: 0.01, tolerance: 1e-20}
        - { x : 100, function value: 0.0001, tolerance: 1e-20}

    """

    __metaclass__ = FunctionMeta

    # noinspection PyPep8Naming
    def evaluate(self, x, logK, piv, index):

        return 10 ** logK * np.power(x / piv, index)

# noinspection PyPep8Naming
class powerlaw_flux(Function):
        r"""
        description :

            A simple power-law with the photon flux in a band used as normalization. This will reduce the correlation
            between the index and the normalization.

        latex : $ \frac{10^{logF}(\gamma+1)} {b^{\gamma+1} - a^{\gamma+1}} x^{\gamma}$

        parameters :

            logF :

                desc : Logarithm of the photon flux between a and b
                initial value : 0
                min : -40
                max : 40
                unit : "1 / (cm2 s)"

            gamma :

                desc : Photon index
                initial value : -2
                min : -10
                max : 10

            a :

                desc : lower bound for the band in which computing the flux F
                initial value : 1.0
                fix : yes

            b :

                desc : upper bound for the band in which computing the flux F
                initial value : 100.0
                fix : yes

        """

        __metaclass__ = FunctionMeta

        # noinspection PyPep8Naming
        def evaluate(self, x, logF, gamma, a, b):

            gp1 = gamma + 1

            return 10 ** logF * gp1 / (b**gp1 - a**gp1) * np.power(x, gamma)

# noinspection PyPep8Naming
class gaussian(Function):
    r"""
    description :

        A Gaussian function

    latex : $ 10^{logK} \frac{1}{\sigma \sqrt{2 \pi}}\exp{\frac{(x-\mu)^2}{2~\sigma^2}} $

    parameters :

        logK :

            desc : Logarithm of normalization
            initial value : 0
            min : -40
            max : 40
            unit : "1 / (keV cm2 s)"

        mu :

            desc : Central value
            initial value : 0.0

        sigma :

            desc : Standard deviation
            initial value : 1.0

    tests :
        - { x : 0.0, function value: 0.3989422804014327, tolerance: 1e-10}
        - { x : -1.0, function value: 0.24197072451914337, tolerance: 1e-9}

    """

    __metaclass__ = FunctionMeta

    __norm = 1.0 / (math.sqrt(2 * np.pi))

    # noinspection PyPep8Naming
    def evaluate(self, x, logK, mu, sigma):

        norm = self.__norm / sigma

        return 10**logK * norm * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.)))


class uniform_prior(Function):
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

    __metaclass__ = FunctionMeta

    def evaluate(self, x, lower_bound, upper_bound, value):

        return np.where( (x >= lower_bound) & (x <= upper_bound), value, 0.0)


class log_uniform_prior(Function):
    r"""
    description :

        A function which is 1/x on the interval lower_bound - upper_bound and 0 outside the interval. The
        extremes of the interval are NOT counted as part of the interval. Lower_bound must be strictly positive.

    latex : $ f(x)=\begin{cases}0 & x \le \text{lower_bound} \\\frac{1}{x} & \text{lower_bound} < x < \text{upper_bound} \\ 0 & x \ge \text{upper_bound} \end{cases}$

    parameters :

        lower_bound :

            desc : Lower bound for the interval
            initial value : 0
            min_value : 0

        upper_bound :

            desc : Upper bound for the interval
            initial value : 100

    tests :
        - { x : 50, function value: (1.0 / 50.0), tolerance: 1e-20}
        - { x : 200, function value: 0, tolerance: 1e-20}

    """

    __metaclass__ = FunctionMeta

    def evaluate(self, x, lower_bound, upper_bound):

        return np.where((x >= lower_bound) & (x <= upper_bound), 1/x, 0.0)


# noinspection PyPep8Naming
class sin(Function):
    r"""
    description :

        A sinusodial function

    latex : $ 10^logK~\sin{2~\pi~f~x + \phi} $

    parameters :

        logK :

            desc : Logarithm of normalization
            initial value : 0
            min : -40
            max : 40
            unit : "1 / (keV cm2 s)"

        f :

            desc : frequency
            initial value : 0.15915494309189535
            min : 0
            max : None
            unit: Hz

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

    __metaclass__ = FunctionMeta

    # noinspection PyPep8Naming
    def evaluate(self, x, logK, f, phi):

        return 10**logK * np.sin(2 * np.pi * f * x + phi)


if has_naima:

    class synchrotron(Function):
        r"""
        description :

            Synchrotron spectrum from an input particle distribution, using Naima (naima.readthedocs.org)

        latex: not available

        parameters :

            B :

                desc : magnetic field
                initial value : 3.24e-6
                unit: Gauss

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

        __metaclass__ = FunctionMeta

        def set_particle_distribution(self, function):

            self._particle_distribution = function

            self._particle_distribution_wrapper = lambda x: function(x / u.eV) / u.eV

        def get_particle_distribution(self):

            return self._particle_distribution

        particle_distribution = property(get_particle_distribution, set_particle_distribution,
                                         doc="""Get/set particle distribution for electrons""")

        # noinspection PyPep8Naming
        def evaluate(self, x, B, emin, emax, need):

            _synch = naima.models.Synchrotron(self._particle_distribution_wrapper, B * u.Gauss,
                                              Eemin = emin * u.GeV, Eemax = emax * u.GeV, nEed = need)

            return _synch._spectrum(x * u.keV)

        def to_dict(self, minimal=False):

            data = super(Function, self).to_dict(minimal)

            if not minimal:

                data['extra_setup'] = {'particle_distribution': self.particle_distribution.path}

            return data


class line(Function):
    r"""
    description :

        A linear function

    latex : $ 10^logK~a * x + 10^logK~b) $

    parameters :

        logK :

            desc : Logarithm of scale factor
            initial value : 0
            min : -40
            max : 40
            unit : "1 / (keV cm2 s)"

        a :

            desc : linear coefficient
            initial value : 1
            min : None
            max : None

        b :

            desc : intercept
            initial value : 0

    """

    __metaclass__ = FunctionMeta

    def evaluate(self, x, logK, a, b):

        k = 10**logK

        return k * a * x + k * b

class identity(Function):
    r"""
    description :

        Return x

    latex : $ x $

    parameters : {}

    """

    __metaclass__ = FunctionMeta

    def evaluate(self, x):

        return x

class bias(Function):
    r"""
    description :

        Return x plus a bias

    latex : $ x + k$

    parameters :

        k :

            desc : Constant value
            initial value : 0

    """

    __metaclass__ = FunctionMeta

    def evaluate(self, x, k):

        return x + k

class band(Function):
    r"""
    description :

        The Band model from Band et al. 1993, implemented however in a way which reduces the covariances between
        the parameters (Calderone et al., MNRAS, 448, 403C, 2015)

    latex : $ \text{(Calderone et al., MNRAS, 448, 403C, 2015)} $

    parameters :

        alpha :
            desc : The photon index for energies smaller than the peak energy
            initial value : -1
            min : -10
            max : 10

        beta :

            desc : photon index for energies greater than the peak energy (only if opt=1, i.e., for the
                   Band model)
            initial value : -2.2
            min : -7
            max : -1

        log_Ep :

            desc : the logarithm of the energy of the spectrum peak in the nuFnu representation
            initial value : 2.2
            min : 0
            max : 3
            unit : keV

        log_Flux :

            desc : the logarithm of the integrated flux in the energy band defined by Emin and Emax
            initial value : -6
            min : -40
            max : 40
            unit : erg / (cm2 * s)

        Emin:

            desc : lower limit of the energy band in which the flux will be computed
            initial value : 10
            min : 0
            max : None
            unit : keV
            fix : yes

        Emax:

            desc : upper limit of the energy band in which the flux will be computed.
            initial value : 10000
            min : 0
            max : None
            unit : keV
            fix : yes

        opt :

            desc : option to select the spectral model (0 corresponds to a cutoff power law, 1 to the Band model)
            initial value : 1
            min : 0
            max : 1
            fix : yes

    """

    __metaclass__ = FunctionMeta


    @staticmethod
    def ggrb_int_cpl( a, Ec, Emin, Emax):

        i1 = gammaincc(2 + a, Emin / Ec) * gamma(2+a)
        i2 = gammaincc(2 + a, Emax / Ec) * gamma(2+a)

        return -Ec * Ec * (i2 - i1)

    @staticmethod
    def ggrb_int_pl(a, b, Ec, Emin, Emax):

        pre = pow(a-b, a-b) * math.exp(b-a) / pow(Ec, b)

        if b != -2:

            return pre / (2+b) * (pow(Emax, 2+b) - pow(Emin, 2+b))

        else:

            return pre * math.log(Emax / Emin)

    def evaluate(self, x, alpha, beta, log_Ep, log_Flux, Emin, Emax, opt):

        assert opt == 0 or opt == 1, "Opt must be either 0 or 1"

        # Cutoff energy

        if alpha == -2:

            Ec = pow(10, log_Ep) / 0.0001 #TRICK: avoid a=-2

        else:

            Ec = pow(10, log_Ep) / (2 + alpha)

        # Split energy

        Esplit = (alpha-beta) * Ec

        # Evaluate model integrated flux and normalization

        if opt==0:

            # Cutoff power law

            intflux = self.ggrb_int_cpl(alpha, Ec, Emin, Emax)

        else:

            # Band model

            if Emin <= Esplit and Esplit <= Emax:

                intflux = ( self.ggrb_int_cpl(alpha,    Ec, Emin,   Esplit) +
                            self.ggrb_int_pl (alpha, beta, Ec, Esplit, Emax) )

            else:

                if Esplit < Emin:

                    intflux = self.ggrb_int_pl(alpha, beta, Ec, Emin, Emax)

                else:

                    raise RuntimeError("Esplit > emax!")

        erg2keV = 6.24151e8

        norm = pow(10, log_Flux) * erg2keV / intflux

        if opt==0:

            # Cutoff power law

            flux = norm * np.power(x / Ec, alpha) * np.exp( - x / Ec)

        else:

            idx = (x < Esplit)

            flux = np.zeros_like( x )

            flux[idx] = ( norm * np.power(x[idx] / Ec, alpha) *
                          np.exp(-x[idx] / Ec) )

            nidx = ~idx

            flux[nidx] = ( norm * pow(alpha-beta, alpha-beta) * math.exp(beta-alpha) *
                           np.power(x[nidx] / Ec, beta) )

        return flux


class log_parabola(Function):
    r"""
    description :

        A log-parabolic function

    latex : $ 10^{logK} \left( \frac{x}{piv} \right)^{\alpha -\beta \log{\left( \frac{x}{piv} \right)}} $

    parameters :

        logK :

            desc : Logarithm of scale factor
            initial value : 0
            min : -40
            max : 40
            unit : "1 / (keV cm2 s)"

        piv :
            desc : Pivot energy (keep this fixed)
            initial value : 1
            fix : yes
            unit : keV

        alpha :

            desc : photon index
            initial value : -2.0
            min : -10
            max : 0

        beta :

            desc : curvature
            initial value : 1.0

    """

    __metaclass__ = FunctionMeta

    def evaluate(self, x, logK, piv, alpha, beta):

        xx = x/piv

        return pow(10,logK) * xx**(alpha - beta * np.log10(xx))

    @property
    def peak_energy(self):
        """
        Returns the peak energy in the nuFnu spectrum

        :return: peak energy in keV
        """

        # Eq. 6 in Massaro et al. 2004
        # (http://adsabs.harvard.edu/abs/2004A%26A...413..489M)

        return self.piv.value * pow(10, (2 + self.alpha.value) / (2 * self.beta.value) )

if has_gsl:

    class cutoff_powerlaw_flux(Function):
        r"""
            description :

                A cutoff power law having the flux as normalization, which should reduce the correlation among
                parameters.

            latex : $ \frac{10^{logK}}{T(b)-T(a)} ~x^{\alpha}~\exp{(x/10^{logC})}~\text{with}~T(x)=-x_{c}^{\alpha+1} \Gamma(\alpha+1, x/10^{logC})~\text{(}\Gamma\text{ is the incomplete gamma function)} $

            parameters :

                logK :

                    desc : Logarithm of the flux F between a and b
                    initial value : 0
                    min : -40
                    max : 40
                    unit : "1 / (cm2 s)"

                alpha :

                    desc : photon index
                    initial value : -2.0
                    min : -10
                    max : 10

                logC :

                    desc : logarithm of the cutoff value
                    initial value : 1.0
                    min : -5
                    max : 10

                a :

                    desc : lower bound for the band in which computing the flux F
                    initial value : 1.0
                    fix : yes

                b :

                    desc : upper bound for the band in which computing the flux F
                    initial value : 100.0
                    fix : yes
            """

        __metaclass__ = FunctionMeta


        @staticmethod
        def _integral(a,b, alpha, ec):

            integrand = lambda x: -ec * pow(x, alpha) * pow(x / ec, -alpha) * gamma_inc(alpha + 1, x / ec)

            return integrand(b) - integrand(a)

        def evaluate(self, x, logK, alpha, logC, a, b):

            ec = pow(10, logC)

            this_integral = self._integral(a,b, alpha, ec)

            return pow(10, logK) / this_integral * np.power(x, alpha) * np.exp(-x / ec)