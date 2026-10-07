"""E2. Spectral-gap scaling of five Metropolis kernels -> results/gaps_<family>.json, exponents_<family>.json"""
import json, sys, time
import numpy as np
from fairq.instances import make_instance
from fairq.subspace import Subspace
from fairq.chains import (local_proposal, local_swap_proposal, uniform_proposal, proposal_params,
                          quantum_proposal, metropolis, gibbs, spectral_gap)

METHODS = ("local", "swap", "uniform", "quantum", "mixed")
SIZES = {"udg": (8, 10, 12, 14, 16, 18, 20), "er": (8, 10, 12, 14, 16, 18), "reg3": (8, 10, 12, 14, 16, 18)}


def sweep(kind, Ns, seeds=(0, 1, 2, 3, 4), betas=(1.0, 3.0), n=3):
    rows = []
    for N in Ns:
        for sd in seeds:
            t0 = time.time()
            inst = make_instance(N, n, sd, kind=kind)
            sub = Subspace(inst)
            w = inst["rates"]; E = sub.energies(w)
            Qq = quantum_proposal(sub, E, proposal_params(N, w, np.random.default_rng(100 + sd)))
            Ql, Qs, Qu = local_proposal(sub), local_swap_proposal(sub), uniform_proposal(sub)
            props = dict(local=Ql, swap=Qs, uniform=Qu, quantum=Qq, mixed=0.5 * Ql + 0.5 * Qq)
            for beta in betas:
                pi = gibbs(E, beta)
                g = {m: spectral_gap(metropolis(props[m], E, beta), pi) for m in METHODS}
                rows.append(dict(family=kind, N=N, seed=sd, beta=beta, F=int(sub.dim),
                                 edges=len(inst["edges"]), **g))
            print(kind, N, sd, sub.dim, f"{time.time() - t0:.1f}s", flush=True)
    return rows


def fit(rows):
    out = {}
    for beta in sorted({r["beta"] for r in rows}):
        R = [r for r in rows if r["beta"] == beta]
        x = np.array([r["N"] for r in R], float)
        A = np.vstack([x, np.ones_like(x)]).T
        out[str(beta)] = {}
        for m in METHODS:
            y = np.log([r[m] for r in R])
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            res = y - A @ coef
            se = np.sqrt(res @ res / (len(y) - 2) / np.sum((x - x.mean()) ** 2))
            out[str(beta)][m] = dict(kappa=float(-coef[0]), se=float(se))
    return out


if __name__ == "__main__":
    kinds = sys.argv[1:] or list(SIZES)
    for kind in kinds:
        rows = sweep(kind, SIZES[kind])
        json.dump(rows, open(f"results/gaps_{kind}.json", "w"), indent=1)
        ex = fit(rows)
        json.dump(ex, open(f"results/exponents_{kind}.json", "w"), indent=1)
        print(json.dumps(ex, indent=1))
