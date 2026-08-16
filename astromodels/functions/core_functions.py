import numpy as np

def plaw_eval(x, K, index, piv):

    return K * np.pow(x / piv, index)


def plaw_flux_norm(index, a, b):
    """Energy flux power law."""

    if index != -2:
        dp2 = 2 + index
        intflux = (np.pow(b, dp2) - np.pow(a, dp2)) / dp2
    else:
        intflux = -np.log(a / b)

    return intflux


def cplaw_eval(x, K, xc, index, piv):

    # Compute it in logarithm to avoid roundoff errors, then raise it
    log_v = index * np.log(x / piv) - (x / xc)
    return K * np.exp(log_v)


def cplaw_inverse_eval(x, K, b, index, piv):

    # Compute it in logarithm to avoid roundoff errors, then raise it
    log_v = index * np.log(x / piv) - x * b
    return K * np.exp(log_v)


def super_cplaw_eval(x, K, piv, index, xc, gamma):

    log_v = index * np.log(x / piv) - np.pow(x / xc, gamma)
    return K * np.exp(log_v)


def band_eval(x, K, alpha, beta, E0, piv):

    out = np.empty_like(x)

    break_point = (alpha - beta) * E0

    factor_ab = np.exp(beta - alpha) * np.pow(break_point / piv, alpha - beta)

    lo_mask = (x < break_point)

    xlo = x[lo_mask]
    out[lo_mask] = K * np.pow(xlo / piv, alpha) * np.exp(-xlo / E0)

    xhi = x[~lo_mask]
    out[~lo_mask] = K * factor_ab * np.pow(xhi / piv, beta)

    return out


def bplaw_eval(x, K, xb, alpha, beta, piv):

    out = np.empty_like(x)

    factor = np.pow(xb / piv, alpha - beta)

    lo_mask = (x < xb)

    xlo = x[lo_mask]
    out[lo_mask] =  K * np.pow(xlo / piv, alpha)

    xhi = x[~lo_mask]
    out[~lo_mask] = K * factor * np.pow(xhi / piv, beta)

    return out


LN2 = np.log(2.0)

def sbplaw_eval(x, K, alpha, be, bs, beta, piv):

    def log_cosh(x):
        """
        An overflow-resistant computation of
        log(0.5 * (exp(x) + exp(-x)))
        """
        absx = np.abs(x)
        return absx + np.log1p(np.exp(-2 * absx)) - LN2

    B = 0.5 * (alpha + beta)
    M = 0.5 * (beta - alpha)

    Mbs = M * bs

    # previously, the log1p/exp contribution to log_cosh was
    # ignored for arg < -6 or > +4; this introduced deviation of
    # up to ~1e-4 vs log_cosh.

    arg_piv = np.log10(piv / be) / bs
    pcosh_piv = log_cosh(arg_piv)

    arg = np.log10(x / be) / bs
    pcosh = log_cosh(arg)

    return K * np.pow(x / piv, B) * np.pow(10.0, Mbs * (pcosh - pcosh_piv))


def bb_eval(x, K, kT):

    return K * x * x / np.expm1(x / kT)


def mbb_eval(x, K, kT):

    arg = x / kT
    exp_arg = np.exp(-arg)

    out = np.pow(arg, 1.5) * exp_arg / np.sqrt(1 - exp_arg)
    return K * out / x


def ggrb_int_pl(a, b, Ec, Emin, Emax):

    pre = np.pow(a - b, a - b) * np.exp(b - a) / np.pow(Ec, b)

    if b != -2:
        b2 = 2 + b
        return pre / b2 * (np.pow(Emax, b2) - np.pow(Emin, b2))
    else:
        return pre * np.log(Emax / Emin)


def non_diss_photoshere_generic(x, K, ec, piv, a, b):

    log_v = a * np.log(x / piv) - np.pow(x / ec, b)
    return K * np.exp(log_v)


def dbl_sbpl(x, K, a1, a2, b1, xp, xb, n1, n2, xpiv):

    xj = xp * np.pow(-(a2 + 2) / (b1 + 2), 1.0 / ((b1 - a2) * n2))

    arg1 = xj / xb
    arg2 = x / xb
    arg3 = x / xj

    inner1 = np.pow(arg2, -a1 * n1) + np.pow(arg2, -a2 * n1)

    inner2 = np.pow(arg1, -a1 * n1) + np.pow(arg1, -a2 * n1)

    out = np.pow(xb / xpiv, a1) * np.pow(
        np.pow(inner1, n2 / n1) + np.pow(arg3, -b1 * n2) * np.pow(inner2, n2 / n1),
        -1 / n2,
    )

    return K * out


def plaw_integral(a, b, K, i, piv):
    if i != -1:
        M = K * piv / (i+1)
        return M * ((b / piv) ** (i + 1) - (a / piv) ** (i + 1))
    else:
        M = K * piv
        return M * (np.log(b / piv) - np.log(a / piv))
