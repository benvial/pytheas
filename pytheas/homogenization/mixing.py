import numpy as np


def maxwell_garnett(incl_vol_frac, eps_incl, eps_host):
    num = 2 * incl_vol_frac * (eps_incl - eps_host) + eps_incl + 2 * eps_host
    den = -incl_vol_frac * (eps_incl - eps_host) + eps_incl + 2 * eps_host
    return num / den * eps_host


def arithmetic_mean(incl_vol_frac, eps_incl, eps_host):
    return (1 - incl_vol_frac) * eps_host + incl_vol_frac * eps_incl


def harmonic_mean(incl_vol_frac, eps_incl, eps_host):
    return 1 / ((1 - incl_vol_frac) / eps_host + incl_vol_frac / eps_incl)


def cplx2realtand(z):
    return z.real, np.abs(z.imag / z.real)
