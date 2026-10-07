"""Fairness measures of the catalogue (Section 2) and exact ex-ante quantities."""
import numpy as np
from scipy.optimize import linprog


def utilitarian(y, w=None):
    y = np.asarray(y, float)
    return float(y.sum() if w is None else np.dot(w, y))

def maxmin(y):
    return float(np.min(y))

def leximin_key(y):
    return tuple(np.sort(y))

def generalized_lorenz(y):
    return np.cumsum(np.sort(np.asarray(y, float)))

def owa(y, weights):
    return float(np.dot(np.sort(np.asarray(y, float)), weights))

def gini(y):
    y = np.asarray(y, float); n = len(y); mu = y.mean()
    return float(np.abs(y[:, None] - y[None, :]).sum() / (2 * n * n * mu)) if mu > 0 else 0.0

def gini_sen(y):
    y = np.asarray(y, float); n = len(y); k = np.arange(1, n + 1)
    return float(np.dot(np.sort(y), (2 * (n - k) + 1) / n ** 2))

def generalized_gini(y, delta=2.0):
    y = np.sort(np.asarray(y, float)); n = len(y); k = np.arange(1, n + 1)
    w = ((n - k + 1) / n) ** delta - ((n - k) / n) ** delta
    return float(np.dot(w, y))

def atkinson_ede(y, eps):
    y = np.asarray(y, float)
    if eps == 1:
        return float(np.exp(np.mean(np.log(y))))
    return float(np.mean(y ** (1 - eps)) ** (1 / (1 - eps)))

def atkinson_index(y, eps):
    return 1.0 - atkinson_ede(y, eps) / float(np.mean(y))

def alpha_fair(y, alpha):
    y = np.asarray(y, float)
    return float(np.sum(np.log(y))) if alpha == 1 else float(np.sum(y ** (1 - alpha)) / (1 - alpha))

def nash_welfare(y):
    return float(np.prod(np.asarray(y, float)))

def kalai_smorodinsky(y, d, a):
    y, d, a = (np.asarray(z, float) for z in (y, d, a))
    return float(np.min((y - d) / (a - d)))

def variance(y):
    return float(np.var(y))

def coefficient_of_variation(y):
    return float(np.std(y) / np.mean(y))

def jain_index(y):
    y = np.asarray(y, float)
    return float(y.sum() ** 2 / (len(y) * (y ** 2).sum()))

def generalized_entropy(y, a):
    y = np.asarray(y, float); r = y / y.mean()
    if a == 0:
        return float(-np.mean(np.log(r)))
    if a == 1:
        return float(np.mean(r * np.log(r)))
    return float(np.mean(r ** a - 1) / (a * (a - 1)))


def utility_matrix(inst):
    """U[x, i] = u_i(x) for every independent set x, and the bit matrix B[x, v] = x_v."""
    F, N, n = inst["F"], inst["N"], inst["n"]
    bits = ((F[:, None] >> np.arange(N)[None, :]) & 1).astype(float)
    A = np.zeros((N, n)); A[np.arange(N), inst["owners"]] = inst["rates"]
    return bits @ A, bits


def ex_ante_maxmin(U):
    """v* = max_p min_i E_p u_i:  max t  s.t.  U^T p >= t 1,  1^T p = 1,  p >= 0 (HiGHS)."""
    m, n = U.shape
    c = np.zeros(m + 1); c[-1] = -1.0
    A_ub = np.hstack([-U.T, np.ones((n, 1))])
    A_eq = np.hstack([np.ones((1, m)), np.zeros((1, 1))])
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(n), A_eq=A_eq, b_eq=[1.0],
                  bounds=[(0, None)] * m + [(None, None)], method="highs")
    assert res.status == 0
    return float(-res.fun), res.x[:m]


def deterministic_maxmin(U):
    return float(np.max(np.min(U, axis=1)))


def support_function(U, theta):
    return float(np.max(U @ theta))
