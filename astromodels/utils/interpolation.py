import numpy as np
from interpolation import interp
from interpolation.splines import eval_linear
from typing import Tuple

class GridInterpolate(object):
    def __init__(self, grid: Tuple[np.ndarray], values: np.ndarray) -> None:
        """
        A numba version of the RegularGridInterpolator from scipy.
        The call is the same as the scipy version

        :param grid: 
        :type grid: np.ndarray
        :param values: 
        :type values: np.ndarray
        :returns: 

        """
        self._grid: np.ndarray = grid
        self._values: np.ndarray() = np.ascontiguousarray(values)
        
        def __call__(self, v) -> np.ndarray:

            return eval_linear(self._grid, self._values, v)

class UnivariateSpline(object):
    
    def __init__(self, x: np.ndarray, y: np.ndarray) -> None:

        """
        A numba version of the InterpolatedUnivariateSpline from scipy.
        The call is the same as the scipy version


        :param x: 
        :type x: np.ndarray
        :param y: 
        :type y: np.ndarray
        :returns: 

        """
        self._x: np.ndarray = x
        self._y: np.ndarray = y

    def __call__(self, v) -> np.ndarray:

        return interp(self._x, self._y, v)



__all__ = ["GridInterpolate", "UnivariateSpline"]
