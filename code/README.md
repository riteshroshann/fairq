# fairq 1.0.0 — Guided Quantum Gibbs Dynamics for Fair Allocation

Reference implementation for the manuscript
"Fairness Is Advantage-Neutral: Guided Quantum Gibbs Dynamics for Fair Allocation of Indivisible Resources"
by Ritesh Roshan.

## Requirements
Python 3.10 or later with numpy, scipy, matplotlib and networkx:

    pip install -r requirements.txt

## Reproduce everything

    sh run_all.sh          # about 60 minutes on one CPU core

or step by step:

| Script | Experiment | Output in `results/` |
|---|---|---|
| `verify.py` | E1: circuit vs exact evolution, symmetry, detailed balance | `verification.json` |
| `circuit_spec.py` | E5: exact gate counts and OpenQASM for the benchmark | `circuit_spec.json`, `proposal_step.qasm` |
| `gaps.py udg er reg3` | E2: spectral-gap scaling, five kernels, three families | `gaps_<family>.json`, `exponents_<family>.json` |
| `tracking.py udg er reg3` | E4: exact warm-start tracking cost along MWU trajectories | `tracking.json` |
| `dynamic_run.py` | E3: repeated max-min allocation, classical vs quantum-mixed | `dynamic.json` |
| `figures.py` | Figures and summary statistics | `fig_gaps.pdf`, `fig_tracking.pdf`, `summary.json` |

Where a single command must stay short, `gaps_part.py <family> <N> <seeds...>` runs one chunk into
`results/parts/` and `gaps_merge.py <family>` merges the chunks and fits the exponents.

## Package layout
- `fairq/instances.py`  conflict graphs (unit-disk, Erdős–Rényi, random 3-regular), rates, owners, colouring, independent sets
- `fairq/fairness.py`   fairness measures of the catalogue, ex-ante max-min LP, support function
- `fairq/subspace.py`   exact palindromic Trotter evolution restricted to the feasible subspace
- `fairq/circuit.py`    gate-level circuit, closed-form counts, depth, OpenQASM 2.0 export, full state-vector simulator
- `fairq/chains.py`     local, swap, uniform and quantum proposals; lazy Metropolis; spectral gaps; tracking cost
- `fairq/dynamic.py`    multiplicative-weights fairness algorithm with warm-started sampling

All randomness is seeded. The results shipped in `results/` were produced by these scripts.
