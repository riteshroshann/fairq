"""E4. Warm-start tracking cost along multiplicative-weights sequences -> results/tracking.json

For checkpoints theta_a -> theta_b (gap rounds apart) of an MWU run with exact Gibbs sampling,
the least k with chi^2(pi_a P_b^k || pi_b) <= eps^2 is computed exactly for each kernel."""
import json, sys, time
import numpy as np
from fairq.instances import make_instance
from fairq.subspace import Subspace
from fairq.fairness import utility_matrix
from fairq.chains import (local_proposal, local_swap_proposal, uniform_proposal, proposal_params,
                          quantum_proposal, metropolis, gibbs, tracking_steps)

METHODS = ("local", "swap", "uniform", "quantum", "mixed")
SIZES = {"udg": (8, 10, 12, 14, 16), "er": (8, 10, 12, 14), "reg3": (8, 10, 12, 14)}


def mwu_sequence(inst, U, beta, T, seed):
    rng = np.random.default_rng(seed)
    n = inst["n"]
    eta = np.sqrt(8 * np.log(n) / T) / float(U.max())
    theta = np.full(n, 1.0 / n); seq = [theta.copy()]
    for _ in range(T):
        x = rng.choice(len(U), p=gibbs(U @ (n * theta), beta))
        theta = theta * np.exp(-eta * U[x]); theta /= theta.sum(); seq.append(theta.copy())
    return np.array(seq)


def track(kind, Ns, seeds=(0, 1, 2), beta=3.0, T=300, gap=30, eps=0.1, n=3):
    rows = []
    for N in Ns:
        for sd in seeds:
            t0 = time.time()
            inst = make_instance(N, n, sd, kind=kind)
            sub = Subspace(inst); U, _ = utility_matrix(inst)
            seq = mwu_sequence(inst, U, beta, T, 500 + sd)
            Ql, Qs, Qu = local_proposal(sub), local_swap_proposal(sub), uniform_proposal(sub)
            base = proposal_params(N, inst["rates"], np.random.default_rng(100 + sd))
            rn = np.linalg.norm(inst["rates"])
            ks = {m: [] for m in METHODS}; drift = []
            for a in range(gap, T + 1 - gap, gap):
                ta, tb = n * seq[a], n * seq[a + gap]
                Ea, Eb = U @ ta, U @ tb
                pa, pb = gibbs(Ea, beta), gibbs(Eb, beta)
                wb = tb[inst["owners"]] * inst["rates"]
                par = [(t, K, gam * rn / np.linalg.norm(wb), h) for (t, K, gam, h) in base]
                Qq = quantum_proposal(sub, Eb, par)
                props = dict(local=Ql, swap=Qs, uniform=Qu, quantum=Qq, mixed=0.5 * Ql + 0.5 * Qq)
                for m in METHODS:
                    ks[m].append(tracking_steps(metropolis(props[m], Eb, beta), pa, pb, eps))
                drift.append(float(np.abs(Eb - Ea).max()))
            rows.append(dict(family=kind, N=N, seed=sd, F=int(sub.dim), beta=beta, gap_rounds=gap, eps=eps,
                             drift=float(np.mean(drift)), **{f"k_{m}": float(np.mean(v)) for m, v in ks.items()}))
            print(kind, N, sd, sub.dim, {m: round(float(np.mean(v)), 1) for m, v in ks.items()},
                  f"{time.time() - t0:.1f}s", flush=True)
    return rows


if __name__ == "__main__":
    kinds = sys.argv[1:] or list(SIZES)
    rows = []
    for kind in kinds:
        rows += track(kind, SIZES[kind])
    json.dump(rows, open("results/tracking.json", "w"), indent=1)
