from collections.abc import Callable
from scipy.signal import lfilter
from abc import ABC, abstractmethod
import numpy as np


from consts import kappa_4_gaussian, kappa_4_t


def _burn_in(a: float, tol: float = 1e-10) -> int:
    """Number of initial samples to discard so the AR(1) recursion has
    forgotten its (non-stationary) start and is effectively stationary."""
    if a <= 0:
        return 0
    return int(np.ceil(np.log(tol) / np.log(a)))



class Process(ABC):
    """
    An instance of ``Process`` represents a single row `n` of the triangular
    array X^(n). The row number ``n`` and the lag ``k``is stored per instance.

    Sweeping a row index set ``N`` over this creates one instance ``Process(n,**params)
    per ``n``in ``N``.
    """

    def __init__(self, n: int, k: int = 1):
        """

        :param n: Row number
        :param k: lag
        """
        self.n = int(n)
        self.k = int(k)

    # ------------------------------------------------------------------
    # example-specific symbols
    # ------------------------------------------------------------------

    @abstractmethod
    def gamma(self, k: int) -> float:
        """
        :param k: lag
        :return: Autocovariance ``gamma^(n)(k)``.
        """

    @abstractmethod
    def Sigma(self, k: int) -> float:
        """
        :param k: lag
        :return: Long-run variance ``Sigma^(n)(k)`` of ``sqrt(n) gamma_hat(k)``.
        """

    @abstractmethod
    def sample(self, length: int, rng: np.random.Generator) -> np.ndarray:
        """
        Realise ``length`` observations of ``X^(n)``
        :param length: number of observations
        :param rng: random number generator
        :return: shape (length) ``X^(n)(1),...,X^(n)(length)``
        """

    def sample_paths(self, reps: int, rng: np.random.Generator) -> np.ndarray:
        """
        ``reps`` independent realisations from row X^(n) of length ``n``, shape ``(reps, n)``.
        Default loops over :meth:`sample`; override for a vectorised version.

        :param reps: number of realisations
        :param rng: random number generator
        :return: shape (reps, n)
        """

        return np.stack([self.sample(self.n, rng) for _ in range(reps)])




class AR1(Process):
    """``X_t = a_n X_{t-1} + eps_t = sum_{i=0}^infty a_n^i eps_{t-i}``,
    ``a_n in (0,1)``.

    The Regime 2 rate for this example is used via its hand-derived
    closed-form exponents (:func:`rate.closed_form`) rather than assembled
    from ``D``, ``m_n``, ``C0``, ``T`` via :func:`rate.regime2`, so only
    ``gamma``, ``Sigma`` and ``sample``/``sample_paths`` (needed by the
    simulation) are implemented; ``L`` and ``C`` are unused and raise
    ``NotImplementedError``.
    """

    def __init__(self, n: int, a: float | Callable[[int], float], kappa_2: float, kappa_4: float, k: int = 1):
        super().__init__(n, k)
        self.a = float(a(n)) if callable(a) else float(a)
        self.kappa_2 = kappa_2
        self.kappa_4 = kappa_4

    def gamma(self, k: int) -> float:
        return self.kappa_2 * self.a ** abs(k) / (1 - self.a**2)

    def Sigma(self, k: int) -> float:
        if self.k != k:
            raise NotImplementedError
        a, k = self.a, abs(k)
        kappa_2 = self.kappa_2  # was: kappa_2 = 1.0
        term1 = kappa_2**2 * (1 + a**2 + a ** (2 * k) * (1 + a**2 + 2 * k * (1 - a**2))) / (1 - a**2) ** 3
        term2 = self.kappa_4 * a ** (2 * k) / (1 - a**2) ** 2  # unchanged, already exact
        return float(term1 + term2)



class AR1Gaussian(AR1):
    def __init__(self, n: int, a: float | Callable[[int], float], kappa_2: float = 1.0, k: int = 1):
        super().__init__(n, a, kappa_2, kappa_4_gaussian, k)

    def sample(self, size: int, rng: np.random.Generator) -> np.ndarray:
        burn_in = _burn_in(self.a)
        eps = np.sqrt(self.kappa_2) * rng.standard_normal(size + burn_in)
        return lfilter([1.0], [1.0, -self.a], eps)[burn_in:]

    def sample_paths(self, reps: int, rng: np.random.Generator) -> np.ndarray:
        burn_in = _burn_in(self.a)
        eps = np.sqrt(self.kappa_2) * rng.standard_normal((reps, self.n + burn_in))
        return lfilter([1.0], [1.0, -self.a], eps, axis=1)[:, burn_in:]


class AR1T(AR1):
    def __init__(self, n: int, a: float | Callable[[int], float], dof: int, kappa_2: float = 1.0, k: int = 1):
        super().__init__(n, a, kappa_2, kappa_4_t(dof), k)
        self.dof = dof
        self.scale = np.sqrt((dof - 2) / dof)

    def sample(self, size: int, rng: np.random.Generator) -> np.ndarray:
        burn_in = _burn_in(self.a)
        eps = np.sqrt(self.kappa_2) * rng.standard_t(self.dof, size=size + burn_in) * self.scale
        return lfilter([1.0], [1.0, -self.a], eps)[burn_in:]

    def sample_paths(self, reps: int, rng: np.random.Generator) -> np.ndarray:
        burn_in = _burn_in(self.a)
        eps = np.sqrt(self.kappa_2) * rng.standard_t(self.dof, size=(reps, self.n + burn_in)) * self.scale
        return lfilter([1.0], [1.0, -self.a], eps, axis=1)[:, burn_in:]
