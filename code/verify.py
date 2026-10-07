"""E1. Verification of the correctness lemmas on concrete instances -> results/verification.json"""
import json
import numpy as np
from fairq.instances import make_instance
from fairq.subspace import Subspace
from fairq.circuit import proposal_circuit, simulate, counts, closed_form_counts
from fairq.chains import local_proposal, local_swap_proposal, metropolis, gibbs, check_reversible


def verify(N, seed, kind="udg", n=3):
    inst = make_instance(N, n, seed, kind=kind)
    sub = Subspace(inst)
    rng = np.random.default_rng(0)
    w = rng.uniform(0.5, 1.5, N); E = sub.energies(w)
    tev, K, gam, h = 1.7, 3, 0.8, 1.0
    err, leak = 0.0, 0.0
    for y in rng.choice(sub.dim, min(5, sub.dim), replace=False):
        gates, nq = proposal_circuit(inst, int(sub.F[y]), w, tev, K, gam, h)
        psi = np.zeros(2 ** nq, complex); psi[0] = 1.0
        psi = simulate(gates, nq, psi)
        e = np.zeros(sub.dim, complex); e[y] = 1.0
        ref = sub.evolve(e, tev, K, gam, h, E, inst["colours"])
        err = max(err, float(np.abs(psi[sub.F] - ref).max()))
        leak = max(leak, float(1 - np.sum(np.abs(psi[sub.F]) ** 2)))
    one, _ = proposal_circuit(inst, 0, w, tev / K, 1, gam, h)
    cf = closed_form_counts(inst); cnt = counts(one)
    counts_match = all(cnt.get(k, 0) == v for k, v in cf.items())
    Q = sub.proposal_matrix(tev, K, gam, h, E, inst["colours"])
    beta = 2.0; pi = gibbs(E, beta)
    db, st = check_reversible(metropolis(0.5 * local_proposal(sub) + 0.5 * Q, E, beta), pi)
    Qs = local_swap_proposal(sub)
    return dict(family=kind, N=N, seed=seed, edges=len(inst["edges"]), F=int(sub.dim), qubits=nq,
                circuit_vs_subspace=err, leakage=leak, closed_form_counts_match=counts_match,
                proposal_asymmetry=float(np.abs(Q - Q.T).max()),
                proposal_column_sum_error=float(np.abs(Q.sum(0) - 1).max()),
                swap_asymmetry=float(np.abs(Qs - Qs.T).max()),
                detailed_balance=db, stationarity=st)


if __name__ == "__main__":
    out = [verify(10, 3), verify(12, 0), verify(14, 1), verify(12, 2, "er"), verify(12, 1, "reg3")]
    json.dump(out, open("results/verification.json", "w"), indent=1)
    for r in out:
        print(json.dumps(r))
