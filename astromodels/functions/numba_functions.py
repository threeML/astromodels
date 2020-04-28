from math import pow
import numpy as np
import numba as nb


@nb.njit()
def band_eval(x, K, alpha, beta, E0, piv):

    n = x.shape[0]
    out = np.empty(n)

    factor_ab = np.exp(beta - alpha) *  pow((alpha - beta) * E0 / piv, alpha - beta)
    break_point = (alpha - beta) * E0
    
    for idx in range(n):

        if x[idx] < break_point:
            out[idx] = K * pow(x[idx]/piv, alpha) * np.exp(-x[idx]/E0)

        else:

            out[idx] = K *  factor_ab  * pow(x[idx]/piv, beta)

    return out
