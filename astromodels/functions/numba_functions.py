import math

import numba as nb
import numpy as np


@nb.njit(fastmath=True, cache=True)
def plaw_eval(x, K, index, piv):

    n = x.shape[0]
    out = np.empty(n)

    for idx in range(n):

        out[idx] = K * math.pow(x[idx] / piv, idx)

    return out


@nb.njit(fastmath=True, cache=True)
def cplaw_eval(x, K, xc, index, piv):

    n = x.shape[0]
    out = np.empty(n)

    for i in range(n):
        # Compute it in logarithm to avoid roundoff errors, then raise it
        log_v = index * np.log(x[i] / piv) - (x[i] / xc)
        out[i] = K * np.exp(log_v)

    return out


@nb.njit(fastmath=True, cache=True)
def cplaw_inverse_eval(x, K, b, index, piv):

    n = x.shape[0]
    out = np.empty(n)

    for i in range(n):
        # Compute it in logarithm to avoid roundoff errors, then raise it
        log_v = index * np.log(x[i] / piv) - x[i] * b
        out[i] = K * np.exp(log_v)

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


@nb.njit(fastmath=True, cache=True)
def sbplaw_eval(x, K, alpha, be, bs, beta, piv):

    n = x.shape[0]
    out = np.zeros(n)

    B = 0.5 * (alpha + beta)
    M = 0.5 * (beta - alpha)

    arg_piv = np.log10(piv / be) / bs

    log2 = np.log(2.0)

    Mbs = M * bs

    if arg_piv < -6.0:
        pcosh_piv = M * bs * (-arg_piv - log2)

    elif arg_piv > 4.0:

        pcosh_piv = M * bs * (arg_piv - log2)
    else:

        pcosh_piv = M * bs * np.log((np.exp(arg_piv) + np.exp(-arg_piv)) / 2.0)

    ten_pcosh_piv = pow(10., pcosh_piv)
        
    for idx in range(n):

        arg = np.log10(x[idx] / be) / bs

        if arg < -6.0:

            pcosh = Mbs * (-arg - log2)

        elif arg > 4.0:

            pcosh = Mbs * (arg - log2)

        else:

            pcosh = Mbs * np.log(0.5 * ((np.exp(arg) + np.exp(-arg))))

        out[idx] = K * pow(x[idx]/piv, B) * pow(10., pcosh)/ten_pcosh_piv


    return out
