"""E3. Repeated max-min allocation: classical (q = 0) versus quantum-mixed (q = 0.5) -> results/dynamic.json"""
import json
import numpy as np
from fairq.instances import make_instance
from fairq.subspace import Subspace
from fairq.fairness import utility_matrix, ex_ante_maxmin, deterministic_maxmin, gini
from fairq.chains import proposal_params
from fairq.dynamic import run_dynamic


def main(N=12, seed=1, n=3, beta=3.0, T=1500, k=8, k0=200, reps=5, kind="udg"):
    inst = make_instance(N, n, seed, kind=kind)
    sub = Subspace(inst); U, _ = utility_matrix(inst)
    vstar, _ = ex_ante_maxmin(U); vdet = deterministic_maxmin(U)
    params = proposal_params(N, inst["rates"], np.random.default_rng(7))
    out = dict(family=kind, N=N, seed=seed, F=int(sub.dim), n=n, beta=beta, T=T, k=k, k0=k0, reps=reps,
               v_star=vstar, v_deterministic=vdet,
               guarantee_slack=float(np.log(sub.dim) / (beta * n) + U.max() * np.sqrt(np.log(n) / (2 * T))),
               runs={})
    for q in (0.0, 0.5):
        res = [run_dynamic(inst, sub, U, beta, T, k, k0, q, params, seed=100 + r) for r in range(reps)]
        vals = [r["value"] for r in res]
        out["runs"][f"q={q}"] = dict(values=vals, mean=float(np.mean(vals)), sd=float(np.std(vals, ddof=1)),
                                     dual_upper=float(np.mean([r["dual_upper"] for r in res])),
                                     gini=float(np.mean([gini(r["avg"]) for r in res])))
    json.dump(out, open("results/dynamic.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
