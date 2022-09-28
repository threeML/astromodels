__author__ = "giacomov"


from astromodels.core.tree import Node
from astromodels.core.parameter import Parameter
import numpy as np


class Polarization(Node):
    def __init__(self, polarization_type="linear"):

        assert polarization_type in [
            "linear",
            "stokes",
        ], "polarization must be linear or stokes"

        self._polarization_type = polarization_type

        Node.__init__(self, "polarization")

    @staticmethod
    def _get_parameter_from_input(
        number_or_parameter, minimum, maximum, what, desc, unit
    ):

        # Try to transform it to float, if it works than we transform it to a parameter

        try:

            number_or_parameter = float(number_or_parameter)

        except TypeError:

            assert isinstance(number_or_parameter, Parameter), (
                "%s must be either a number or a " "parameter instance" % what
            )

            # So this is a Parameter instance already. Enforce that it has the right maximum and minimum

            parameter = number_or_parameter

            assert parameter.min_value == minimum, "%s must have a minimum of %s" % (
                what,
                minimum,
            )
            assert parameter.max_value == maximum, "%s must have a maximum of %s" % (
                what,
                maximum,
            )

        else:

            # This was a float. Enforce that it has a legal value

            assert (
                minimum <= number_or_parameter <= maximum
            ), "%s cannot have a value of %s, " "it must be %s <= %s <= %s" % (
                what,
                number_or_parameter,
                minimum,
                what,
                maximum,
            )

            parameter = Parameter(
                what,
                number_or_parameter,
                desc=desc,
                min_value=minimum,
                max_value=maximum,
                unit=unit,
                free=True,
            )

        return parameter


# TODO: add transform between polarizations

class LinearPolarization(Polarization):
    def __init__(self, degree, angle):
        """
        Linear parameterization of polarization

        :param degree: The polarization degree
        :param angle: The polarization angle
        """
        super(LinearPolarization, self).__init__(polarization_type="linear")

        if callable(degree):
            self.degree = LinearParameter('degree', degree)
        else:
            self.degree = self._get_parameter_from_input(degree, 0, 100, 'degree', 'Polarization degree', 'dimensionless_unscaled')

        if callable(angle):
            self.angle = LinearParameter('angle', angle)
        else:
            self.angle = self._get_parameter_from_input(angle, 0, 180, 'angle', 'Polarization angle', 'deg')

        self._add_child(self.degree)
        self._add_child(self.angle)

    def __call__(self, energies, stokes):

        if stokes == 'Q':
            return self.degree(energies) * np.cos(2.0 * np.radians(self.angle(energies)))
        elif stokes == 'U':
            return self.degree(energies) * np.sin(2.0 * np.radians(self.angle(energies)))
        return 1


class StokesPolarization(Polarization):
    def __init__(self, I=None, Q=None, U=None, V=None):
        """
        Stokes parameterization of polarization
        """
        super(StokesPolarization, self).__init__(polarization_type='stokes')
        self._Q = StokesParameter('Q', Q)
        self._add_child(self._Q)
        self._U = StokesParameter('U', U)
        self._add_child(self._U)

    def __call__(self, energies, stokes):
        if stokes == 'Q':
            return self._Q(energies)
        elif stokes == 'U':
            return self._U(energies)
        return 1

#     def to_linear_polarization(self):
#         # polarization angle
# #        psi = 0.5 * np.arctan2(U_bin, Q_bin)
#
#         # polarization fraction
# #        frac = np.sqrt(Q_bin ** 2 + U_bin ** 2) / I_bin
#
#         pass
#
#         #angle = 0.5 * np.arctan2(se)
#
#


class StokesParameter(Node):

    def __init__(self, name, value):

        assert name in ['I', 'Q', 'U', 'V']
        Node.__init__(self, name)
        self._add_child(value)
        self.value = value

    def __call__(self, energies):
        return self.value(energies)


class LinearParameter(Node):
    def __init__(self, name, value):
        assert name in ['degree', 'angle']
        Node.__init__(self, name)
        self._add_child(value)
        self.value = value

    def __call__(self, energies):
        return self.value(energies)
