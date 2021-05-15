import astropy.units as astropy_units
import numpy as np
from past.utils import old_div

import astromodels.functions.numba_functions as nb_func
from astromodels.core.units import get_units
from astromodels.functions.function import (Function1D, FunctionMeta,
                                            ModelAssertionViolation)


def get_polynomial(order: int) -> Function1D:
    """
    get a polynomial function of order

    :param order: the order of the polynomical
    :type order: int
    :returns: 

    """
    return [Constant(), Line(),Quadratic(), Cubic(), Quartic()][order]
    
    
    
    


class Constant(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        Return k

    latex : $ k $

    parameters :

        k :

            desc : Constant value
            initial value : 0

    """
    def _set_units(self, x_unit, y_unit):
        self.k.unit = y_unit

    def evaluate(self, x, k):

        return k * np.ones(np.shape(x))


class Line(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A linear function

    latex : $ b * x + a $

    parameters :

        a :

            desc :  intercept
            initial value : 0

        b :

            desc : coeff
            initial value : 1

    """
    def _set_units(self, x_unit, y_unit):
        # a has units of y_unit / x_unit, so that a*x has units of y_unit
        self.a.unit = y_unit

        # b has units of y
        self.b.unit = y_unit / x_unit

    def evaluate(self, x, a, b):
        return b * x + a


class Quadratic(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A Quadratic function

    latex : $ a + b \cdot x + c \cdot x^2 $

    parameters :

        a :

            desc : coefficient
            initial value : 1

        b :

            desc : coefficient
            initial value : 1

        c :

            desc : coefficient
            initial value : 1


    """
    def _set_units(self, x_unit, y_unit):
        # a has units of y_unit / x_unit, so that a*x has units of y_unit
        self.a.unit = y_unit

        # b has units of y
        self.b.unit = y_unit / x_unit

        self.c.unit = y_unit / (x_unit) ** 2

    def evaluate(self, x, a, b, c):
        return a + b * x + c * x * x


class Cubic(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A cubic function

    latex : $ a + b \cdot x + c \cdot x^2 + d \cdot x^3$

    parameters :

        a :

            desc : coefficient
            initial value : 1

        b :

            desc : coefficient
            initial value : 1

        c :

            desc : coefficient
            initial value : 1

        d :

            desc : coefficient
            initial value : 1


    """

    def _set_units(self, x_unit, y_unit):
        # a has units of y_unit / x_unit, so that a*x has units of y_unit
        self.a.unit = y_unit

        # b has units of y
        self.b.unit = y_unit / x_unit

        self.c.unit = y_unit / (x_unit) ** 2

        self.d.unit = y_unit / (x_unit) ** 3

    def evaluate(self, x, a, b, c, d):

        x2 = x * x

        x3 = x2 * x

        return a + b * x + c * x2 + d * x3


class Quartic(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A quartic function

    latex : $ a + b \cdot x + c \cdot x^2 + d \cdot x^3 + e \cdot x^4$

    parameters :

        a :

            desc : coefficient
            initial value : 1

        b :

            desc : coefficient
            initial value : 1

        c :

            desc : coefficient
            initial value : 1

        d :

            desc : coefficient
            initial value : 1

        e :

            desc : coefficient
            initial value : 1


    """

    def _set_units(self, x_unit, y_unit):
        # a has units of y_unit / x_unit, so that a*x has units of y_unit
        self.a.unit = y_unit

        # b has units of y
        self.b.unit = y_unit / x_unit

        self.c.unit = y_unit / (x_unit) ** 2

        self.d.unit = y_unit / (x_unit) ** 3

        self.e.unit = y_unit / (x_unit) ** 4

    def evaluate(self, x, a, b, c, d, e):

        x2 = x * x

        x3 = x2 * x

        x4 = x3 * x

        return a + b * x + c * x2 + d * x3 + e * x4
