"""Markov kernels on F: local, swap, uniform, quantum and mixtures; lazy Metropolis kernel;
spectral gaps; warm-start tracking cost. Proposal matrices are column-stochastic, Q[x, y] = Q(x | y)."""
import numpy as np


def local_proposal(sub):
    """Pick v uniformly; propose x = y + e_v (mod 2) if independent, else stay."""
    dim, N = sub.dim, sub.N
    Q = np.zeros((dim, dim))
    for v in range(N):
        lo, hi = sub.pairs[v]
        Q[hi, lo] += 1.0 / N
        Q[lo, hi] += 1.0 / N
    Q[np.diag_indices(dim)] += 1.0 - Q.sum(axis=0)
    return Q


def swap_proposal(sub):
    """Pick an edge (u, v) uniformly; if exactly one endpoint is occupied and moving the particle
    across the edge keeps the set independent, propose the move; else stay."""
    F, dim, edges = sub.F, sub.dim, sub.inst["edges"]
    Q = np.zeros((dim, dim))
    for (u, v) in edges:
        for a, b in ((u, v), (v, u)):
            src = np.nonzero((((F >> a) & 1) == 1) & (((F >> b) & 1) == 0))[0]
            tgt = F[src] ^ (1 << a) ^ (1 << b)
            j = np.minimum(np.searchsorted(F, tgt), dim - 1)
            ok = F[j] == tgt
            Q[j[ok], src[ok]] += 1.0 / len(edges)
    Q[np.diag_indices(dim)] += 1.0 - Q.sum(axis=0)
    return Q


def local_swap_proposal(sub):
    return 0.5 * local_proposal(sub) + (0.5 * swap_proposal(sub) if sub.inst["edges"] else 0.5 * np.eye(sub.dim))


def uniform_proposal(sub):
    return np.full((sub.dim, sub.dim), 1.0 / sub.dim)


def proposal_params(N, w, rng, M=3):
    """Random proposal parameters in the style of Layden et al. (2023):
    H = (1 - g) a H_theta + g H_PXP, a = sqrt(N) / ||w||, g ~ U[0.25, 0.6], t ~ U[2, 10], tau <= 0.5.
    Returns tuples (t_ev, K, gamma, h)."""
    a = np.sqrt(N) / np.linalg.norm(w)
    out = []
    for _ in range(M):
        g = rng.uniform(0.25, 0.6); t = rng.uniform(2.0, 10.0)
        out.append((t, int(np.ceil(t / 0.5)), (1 - g) * a, g))
    return out


def quantum_proposal(sub, E, params):
    return sum(sub.proposal_matrix(t, K, gam, h, E, sub.inst["colours"]) for (t, K, gam, h) in params) / len(params)


def metropolis(Q, E, beta, lazy=True):
    """Row-stochastic P[y, x] from symmetric proposal Q[x, y], target pi(x) ~ exp(beta E(x))."""
    A = np.minimum(1.0, np.exp(beta * (E[None, :] - E[:, None])))
    P = Q.T * A
    np.fill_diagonal(P, 0.0)
    P[np.diag_indices_from(P)] = 1.0 - P.sum(axis=1)
    return 0.5 * (np.eye(len(E)) + P) if lazy else P


def gibbs(E, beta):
    z = np.exp(beta * (E - E.max()))
    return z / z.sum()


def spectral_gap(P, pi):
    d = np.sqrt(pi)
    S = (d[:, None] * P) / d[None, :]
    ev = np.linalg.eigvalsh(0.5 * (S + S.T))
    return float(1.0 - np.sort(ev)[-2])


def check_reversible(P, pi):
    M = pi[:, None] * P
    return float(np.abs(M - M.T).max()), float(np.abs(pi @ P - pi).max())


def tracking_steps(P, pi_start, pi_target, eps=0.1, kmax=10 ** 9):
    """Least k with chi^2(pi_start P^k || pi_target) <= eps^2, by exact spectral decomposition."""
    d = np.sqrt(pi_target)
    S = (d[:, None] * P) / d[None, :]
    lam, V = np.linalg.eigh(0.5 * (S + S.T))
    c2 = (V.T @ (d * (pi_start / pi_target - 1.0))) ** 2
    keep = lam < 1 - 1e-12
    lam, c2 = np.clip(lam[keep], 0.0, 1.0), c2[keep]
    chi = lambda k: float(np.sum(c2 * lam ** (2 * k)))
    if chi(0) <= eps ** 2:
        return 0
    lo, hi = 0, 1
    while chi(hi) > eps ** 2:
        lo, hi = hi, 2 * hi
        if hi > kmax:
            return kmax
    while hi - lo > 1:
        m = (lo + hi) // 2
        lo, hi = (lo, m) if chi(m) <= eps ** 2 else (m, hi)
    return hi
