from functools import lru_cache

import numpy as np
from scipy.special import betaln
from scipy.special import gamma as gamma_fn

# -------------- Gaussian Constants --------------#


def gaussian_moment(q: int) -> float:
    """
    Returns q-th moment of standard normal distribution
    """
    return 2 ** (q / 2) * gamma_fn((q + 1) / 2) / np.sqrt(np.pi)


kappa_4_gaussian = 0
m4_gaussian = gaussian_moment(4)
m6_gaussian = gaussian_moment(6)

# -------------- t-distribution Constants --------------#


@lru_cache(maxsize=128)
def kappa_4_t(dof):
    return 6 / (dof - 4)


def t_moment(q: int, dof: int) -> float:
    if dof <= q:
        raise ValueError(f"E[Y^{q}] does not exist: need dof > q (got dof={dof}, q={q})")

    if q % 2 == 1:
        return 0.0

    log_num = betaln((q + 1) / 2, (dof - q) / 2)
    log_den = betaln(0.5, dof / 2)
    return (dof - 2) ** (q / 2) * np.exp(log_num - log_den)


@lru_cache(maxsize=128)
def m4_t(dof: int):
    return t_moment(4, dof)


@lru_cache(maxsize=128)
def m6_t(dof: int):
    return t_moment(6, dof)
