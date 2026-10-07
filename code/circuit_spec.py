"""E5. Exact gate-level specification for the benchmark instance -> results/circuit_spec.json, proposal_step.qasm"""
import json
from fairq.instances import make_instance
from fairq.circuit import proposal_circuit, counts, closed_form_counts, depth, to_qasm


def main(N=20, seed=0, n=3, kind="udg", t_ev=5.0, K=10, gamma=1.0, h=0.4):
    inst = make_instance(N, n, seed, kind=kind)
    w = inst["rates"]; deg = [len(x) for x in inst["nb"]]
    step, nq = proposal_circuit(inst, 0, w, t_ev / K, 1, gamma, h)
    full, _ = proposal_circuit(inst, 0, w, t_ev, K, gamma, h)
    spec = dict(family=kind, N=N, seed=seed, edges=len(inst["edges"]), F=int(len(inst["F"])),
                max_degree=max(deg), degrees=deg, colour_classes=inst["colours"], qubits=nq, ancillas=nq - N,
                t_ev=t_ev, K=K, gamma=gamma, h=h,
                step_counts=counts(step), step_counts_closed_form=closed_form_counts(inst),
                step_counts_toffoli_decomposed=counts(step, True),
                step_depth=depth(step, nq), step_two_qubit_depth=depth(step, nq, True),
                full_counts=counts(full), full_depth=depth(full, nq), full_two_qubit_depth=depth(full, nq, True),
                edges_list=inst["edges"], rates=[round(float(r), 6) for r in inst["rates"]],
                owners=[int(a) for a in inst["owners"]])
    json.dump(spec, open("results/circuit_spec.json", "w"), indent=1)
    open("results/proposal_step.qasm", "w").write(to_qasm(step, nq, N))
    print(json.dumps({k: v for k, v in spec.items() if k not in ("edges_list", "degrees", "rates", "owners")}, indent=1))


if __name__ == "__main__":
    main()
