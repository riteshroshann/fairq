"""Gate-level proposal circuit, gate counts, depth, OpenQASM 2.0 export and a full
state-vector simulator used to verify the gate list against the subspace evolution.

Qubits 0..N-1 hold x_v; qubits N..N+A-1 are ancillas with A = max(Delta - 1, 0).
Gate set: x, h, p(phi), rx(theta), rz(theta), cx, ccx.
"""
import numpy as np
from collections import Counter


def crx(gates, c, t, phi):
    """CR_X(phi) = (I x H) CR_Z(phi) (I x H), CR_Z(phi) = [CX, RZ(-phi/2), CX, RZ(phi/2)] in time order."""
    gates += [("h", (t,), None), ("cx", (c, t), None), ("rz", (t,), -phi / 2),
              ("cx", (c, t), None), ("rz", (t,), phi / 2), ("h", (t,), None)]


def driver_gate(gates, v, nbrs, phi, N):
    """exp(-i (phi/2) X_v Pi_v): R_X(phi) on v controlled on every neighbour being |0>."""
    d = len(nbrs)
    if d == 0:
        gates.append(("rx", (v,), phi))
        return
    for u in nbrs:
        gates.append(("x", (u,), None))
    if d == 1:
        crx(gates, nbrs[0], v, phi)
    else:
        anc = [N + k for k in range(d - 1)]
        chain = [("ccx", (nbrs[0], nbrs[1], anc[0]), None)]
        for k in range(2, d):
            chain.append(("ccx", (nbrs[k], anc[k - 2], anc[k - 1]), None))
        gates += chain
        crx(gates, anc[-1], v, phi)
        gates += list(reversed(chain))
    for u in nbrs:
        gates.append(("x", (u,), None))


def proposal_circuit(inst, y_mask, w, t_ev, K, gamma, h):
    """Gate list for preparing |y>, applying S(tau)^K, before measurement. Returns (gates, n_qubits)."""
    N, nb, colours = inst["N"], inst["nb"], inst["colours"]
    tau = t_ev / K
    gates = [("x", (v,), None) for v in range(N) if (y_mask >> v) & 1]

    def cost(s):
        for v in range(N):
            gates.append(("p", (v,), s * gamma * w[v]))

    def drive(cls, s):
        for v in cls:
            driver_gate(gates, v, nb[v], 2 * s * h, N)

    for _ in range(K):
        cost(tau / 2)
        for cls in colours[:-1]:
            drive(cls, tau / 2)
        drive(colours[-1], tau)
        for cls in reversed(colours[:-1]):
            drive(cls, tau / 2)
        cost(tau / 2)
    return gates, N + max(max(len(x) for x in nb) - 1, 0)


def closed_form_counts(inst):
    """Per-step gate counts from the closed-form expressions of Section 5.4."""
    nb, last = inst["nb"], set(inst["colours"][-1])
    m = {v: (1 if v in last else 2) for v in range(inst["N"])}
    d = {v: len(nb[v]) for v in range(inst["N"])}
    two = sum(2 * m[v] for v in d if d[v] >= 1)
    return dict(p=2 * inst["N"], x=sum(2 * m[v] * d[v] for v in d),
                ccx=sum(2 * m[v] * max(d[v] - 1, 0) for v in d), cx=two, h=two, rz=two,
                rx=sum(m[v] for v in d if d[v] == 0))


def counts(gates, decompose_toffoli=False):
    c = Counter(g[0] for g in gates)
    if decompose_toffoli:
        k = c.pop("ccx", 0)
        c["cx"] += 6 * k; c["h"] += 2 * k; c["t_tdg"] += 7 * k; c["s"] += k
    return dict(c)


def depth(gates, nq, two_qubit_only=False):
    lvl = [0] * nq; d = 0
    for _, qs, _ in gates:
        if two_qubit_only and len(qs) == 1:
            continue
        L = max(lvl[q] for q in qs) + 1
        for q in qs:
            lvl[q] = L
        d = max(d, L)
    return d


def to_qasm(gates, nq, n_meas):
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{nq}];", f"creg c[{n_meas}];"]
    for name, qs, par in gates:
        args = ",".join(f"q[{q}]" for q in qs)
        lines.append(f"{name}({par:.12f}) {args};" if par is not None else f"{name} {args};")
    lines += [f"measure q[{v}] -> c[{v}];" for v in range(n_meas)]
    return "\n".join(lines) + "\n"


_H = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
_X = np.array([[0, 1], [1, 0]], complex)


def _one(name, par):
    if name == "x":
        return _X
    if name == "h":
        return _H
    if name == "p":
        return np.diag([1, np.exp(1j * par)])
    if name == "rz":
        return np.diag([np.exp(-1j * par / 2), np.exp(1j * par / 2)])
    if name == "rx":
        c, s = np.cos(par / 2), np.sin(par / 2)
        return np.array([[c, -1j * s], [-1j * s, c]])
    raise ValueError(name)


def simulate(gates, nq, psi):
    """Full state-vector simulation; bit q of the basis index is qubit q."""
    psi = np.array(psi, complex).reshape([2] * nq)
    ax = lambda q: nq - 1 - q
    for name, qs, par in gates:
        if name in ("cx", "ccx"):
            *ctrl, t = qs
            sl = [slice(None)] * nq
            for c in ctrl:
                sl[ax(c)] = 1
            sub = psi[tuple(sl)]
            tax = ax(t) - sum(1 for c in ctrl if ax(c) < ax(t))
            psi[tuple(sl)] = np.flip(sub, axis=tax).copy()
        else:
            U = _one(name, par)
            psi = np.moveaxis(np.tensordot(U, psi, axes=([1], [ax(qs[0])])), 0, ax(qs[0]))
    return psi.reshape(-1)
