import ctypes
import math

import numba as nb
import numpy as np

# from numba.extending import get_cython_function_address

# addr1 = get_cython_function_address(
#     "scipy.special.cython_special", "gammaincc")
# functype1 = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
# gammaincc_fn = functype1(addr1)


# @nb.vectorize('float64(float64, float64)')
# def vec_gammaincc(x, y):
#     return gammaincc_fn(x, y)


# addr2 = get_cython_function_address("scipy.special.cython_special", "gamma")
# functype2 = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
# gamma_fn = functype2(addr2)


# @nb.vectorize('float64(float64)')
# def vec_gamma(x):
#     return gamma_fn(x)


@nb.njit(fastmath=True, cache=True)
def plaw_eval(x, K, index, piv):

    n = x.shape[0]
    out = np.empty(n)

    for idx in range(n):

        out[idx] = K * math.pow(x[idx] / piv, index)

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
def super_cplaw_eval(x, K, piv, index, xc, gamma):

    n = x.shape[0]
    out = np.empty(n)

    for i in range(n):

        log_v = index * np.log(x[i] / piv) - gamma*(x[i] / xc)

        out[i] = K * np.exp(log_v)

    return out


@nb.njit(fastmath=True, cache=True)
def band_eval(x, K, alpha, beta, E0, piv):

    n = x.shape[0]
    out = np.empty(n)

    break_point = (alpha - beta) * E0

    factor_ab = np.exp(beta - alpha) * \
        math.pow(break_point / piv, alpha - beta)

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
        pcosh_piv = Mbs * (-arg_piv - log2)

    elif arg_piv > 4.0:

        pcosh_piv = Mbs * (arg_piv - log2)
    else:

        pcosh_piv = Mbs * np.log((np.exp(arg_piv) + np.exp(-arg_piv)) / 2.0)

    ten_pcosh_piv = math.pow(10., pcosh_piv)

    for idx in range(n):

        arg = np.log10(x[idx] / be) / bs

        if arg < -6.0:

            pcosh = Mbs * (-arg - log2)

        elif arg > 4.0:

            pcosh = Mbs * (arg - log2)

        else:

            pcosh = Mbs * np.log(0.5 * ((np.exp(arg) + np.exp(-arg))))

        out[idx] = K * math.pow(x[idx]/piv, B) * \
            math.pow(10., pcosh)/ten_pcosh_piv

    return out


@nb.njit(fastmath=True, cache=True)
def bb_eval(x, K, kT):

    n = x.shape[0]
    out = np.empty(n)

    for idx in range(n):

        arg = x[idx]/kT
        out[idx] = K * x[idx] * x[idx] / np.expm1(arg)

    return out

# @nb.njit(fastmath=True, cache=True)
# def bbrad_eval(x, K, kT):

#     n = x.shape[0]
#     out = np.empty(n)

#     for idx in range(n):

#         arg = x[idx]/kT
#         out[idx] = K * x[idx] * x[idx] / np.expm1(arg)

#     return out



# band calderone


@nb.njit(fastmath=True, cache=True)
def ggrb_int_pl(a, b, Ec, Emin, Emax):

    pre = math.pow(a - b, a - b) * math.exp(b - a) / math.pow(Ec, b)

    if b != -2:
        b2 = 2+b

        return pre / (b2) * (math.pow(Emax, b2) - math.pow(Emin, b2))

    else:

        return pre * math.log(Emax/Emin)


# @nb.njit(fastmath=True, cache=True)
# def ggrb_int_cpl(a, Ec, Emin, Emax):

#     # Gammaincc does not support quantities
#     i1 = vec_gammaincc(2 + a, Emin/Ec) * vec_gamma(2 + a)
#     i2 = vec_gammaincc(2 + a, Emax/Ec) * vec_gamma(2 + a)

#     return -Ec * Ec * (i2 - i1)
