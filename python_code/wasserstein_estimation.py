"""Simulate the Wasserstein-1 distance between sqrt(n)(gamma - gamma_hat) and its normal limit."""

import numpy as np
from scipy.special import ndtri

from ar1_model import Process


def _gamma_hat(paths: np.ndarray, k: int) -> np.ndarray:
    """

    :param paths: shape (R,n) contains R paths of the process X^(n) each with length n
    :param k: lag parameter
    :return: shape (R) contains R realisations of the empirical autocovariance evaluated at lag k
    """
    n = paths.shape[1]
    return np.sum(paths[:, k:] * paths[:, : n - k], axis=1) / n


def estimate_W(
    process: Process,
    R: int,
    B: int,
    rng: np.random.Generator,
    normalise=False,
    return_covariance_distribution=False,
) -> float:
    """
    Draws ``R`` paths of length ``process.n`` from ``process``, computes
    empirical covariances at lag ``process.k``, standardizes them to
    ``sqrt(Sigma) * (gamma(k) - gamma_hat(k)). The sorted values of this act
    as quantiles at ``(2i-1)/(2R)``. Matches those quantiles agains the normal ones
    and averages the results over ``B``independent repetitions.

    :param process: the process instance (row and lag fixed) that is simulated
    :param R: number of covariance realisations used to approximate the quantile function
    :param B: number of Wasserstein distances are averaged over
    :param rng: random number generator
    :param normalise: normalise rate according to Corollary 6
    :return: single wasserstein distance approximation from row n and,
        optionally, the standardized covariance distribution pooled across
        all ``B`` repetitions (shape ``(B * R,)``)
    """
    n, k = process.n, process.k
    gamma_k = process.gamma(k)
    target = ndtri((2.0 * np.arange(1, R + 1) - 1.0) / (2.0 * R))
    if not normalise:
        target *= np.sqrt(process.Sigma(k))

    total = 0.0
    zs = []
    for _ in range(B):
        z = np.sqrt(n) * (_gamma_hat(process.sample_paths(R, rng), k) - gamma_k)
        if normalise:
            z /= np.sqrt(process.Sigma(k))
        z.sort()
        total += float(np.mean(np.abs(z - target)))
        if return_covariance_distribution:
            zs.append(z)
    if return_covariance_distribution:
        return total / B, np.sort(np.concatenate(zs))
    else:
        return total / B


def sweep(
    process_cls: type[Process],
    N,
    R: int,
    B: int,
    rng: np.random.Generator,
    normalise=False,
    **params,
) -> dict[int, float]:
    """
    Estimates Wasserstein distance in each row ``n`` in ``N`` using ``R`` path
    realisations and averaging over ``B`` repetitions.

    :param process_cls: the process that is simulated
    :param N: list of rows
    :param R: number of covariance realisations used to approximate the quantile function
    :param B: number of Wasserstein distances are averaged over
    :param rng: random number generator
    :param params: parameters used to intialize a row of process_cls
    :return: dict containing rows from ``N`` as keys and wassertein distance estimates as values
    """
    return {int(n): estimate_W(process_cls(int(n), **params), R, B, rng, normalise=normalise) for n in N}
