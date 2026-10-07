"""Exact evolution restricted to the feasible subspace H_F = span{|x> : x independent}.

Implements the palindromic Trotter step
    S(tau) = C(tau/2) D_1(tau/2) ... D_{chi-1}(tau/2) D_chi(tau) D_{chi-1}(tau/2) ... D_1(tau/2) C(tau/2)
with C(s) = prod_v P(s gamma w_v) and D_c(s) = prod_{v in C_c} exp(-i s h X_v Pi_v).
On H_F, exp(-i a X_v Pi_v) rotates each pair (x, x + e_v) with x_v = 0 and Pi_v|x> = |x> by R_X(2a)
and leaves every other amplitude unchanged; this is exact.
"""
import numpy as np


class Subspace:
    def __init__(self, inst):
        self.inst = inst
        self.F = inst["F"]
        self.N = inst["N"]
        self.dim = len(self.F)
        self.bits = ((self.F[:, None] >> np.arange(self.N)[None, :]) & 1).astype(np.int8)
        nbmask = np.array([sum(1 << u for u in inst["nb"][v]) for v in range(self.N)], dtype=np.int64)
        self.pairs = []
        for v in range(self.N):
            free = (self.F & nbmask[v]) == 0
            lo = np.nonzero(free & (((self.F >> v) & 1) == 0))[0]
            hi = np.searchsorted(self.F, self.F[lo] | (1 << v))
            assert np.all(self.F[hi] == (self.F[lo] | (1 << v)))
            self.pairs.append((lo, hi))

    def energies(self, w):
        """E(x) = sum_v w_v x_v = -<x|H_theta|x>."""
        return self.bits @ np.asarray(w, float)

    def apply_cost(self, psi, s, gamma, E):
        ph = np.exp(1j * s * gamma * E)
        return psi * (ph[:, None] if psi.ndim == 2 else ph)

    def apply_driver_vertex(self, psi, v, a):
        lo, hi = self.pairs[v]
        c, sn = np.cos(a), np.sin(a)
        A = psi[lo].copy(); B = psi[hi].copy()
        psi[lo] = c * A - 1j * sn * B
        psi[hi] = c * B - 1j * sn * A
        return psi

    def strang_step(self, psi, tau, gamma, h, E, colours):
        psi = self.apply_cost(psi, tau / 2, gamma, E)
        for cls in colours[:-1]:
            for v in cls:
                psi = self.apply_driver_vertex(psi, v, h * tau / 2)
        for v in colours[-1]:
            psi = self.apply_driver_vertex(psi, v, h * tau)
        for cls in reversed(colours[:-1]):
            for v in cls:
                psi = self.apply_driver_vertex(psi, v, h * tau / 2)
        return self.apply_cost(psi, tau / 2, gamma, E)

    def evolve(self, psi, t_ev, K, gamma, h, E, colours):
        psi = np.array(psi, dtype=complex)
        tau = t_ev / K
        for _ in range(K):
            psi = self.strang_step(psi, tau, gamma, h, E, colours)
        return psi

    def proposal_matrix(self, t_ev, K, gamma, h, E, colours):
        """Q[x, y] = |<x| S(tau)^K |y>|^2 for all x, y in F (column y is the proposal from y)."""
        U = self.evolve(np.eye(self.dim), t_ev, K, gamma, h, E, colours)
        return np.abs(U) ** 2

    def sample_proposal(self, y, t_ev, K, gamma, h, E, colours, rng):
        psi = np.zeros(self.dim, complex); psi[y] = 1.0
        p = np.abs(self.evolve(psi, t_ev, K, gamma, h, E, colours)) ** 2
        return int(rng.choice(self.dim, p=p / p.sum()))
