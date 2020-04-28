import math

import numba as nb
import numpy as np


@nb.njit(fastmath=True, cache=True)
def plaw_eval(x,K,index,piv):

    n = x.shape[0]
    out = np.empty(n)

    for idx in range(n):

        out[idx] = K * math.pow(x[idx]/piv, idx)

    return out



@nb.njit(fastmath=True, cache=True)
def band_eval(x, K, alpha, beta, E0, piv):

    n = x.shape[0]
    out = np.empty(n)

    factor_ab = np.exp(beta - alpha) * math.pow((alpha - beta) * E0 / piv, alpha - beta)
    break_point = (alpha - beta) * E0

    for idx in range(n):

        if x[idx] < break_point:
            out[idx] = K * math.pow(x[idx] / piv, alpha) * np.exp(-x[idx] / E0)

        else:

            out[idx] = K * factor_ab * math.pow(x[idx] / piv, beta)

    return out


@nb.njit(fastmath=True, cache=True)
def bplaw_eval(x, K, xb, alpha, beta, piv):

    n = x.shape[0]
    out = np.empty(n)

    factor = math.pow(xb / piv, alpha - beta)

    for idx in range(n):

        if x[idx] < xb:

            out[idx] = K * math.pow(x[idx] / piv, alpha)
        else:

            out[idx] = K * factor * math.pow(x[idx] / piv, beta)

    return out
