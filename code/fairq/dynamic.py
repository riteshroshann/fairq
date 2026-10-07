"""Algorithm of Section 6: multiplicative-weights max-min fairness with warm-started Metropolis
sampling, with local and (optionally) quantum proposals."""
import numpy as np


def _local_move(sub, y, rng):
    v = rng.integers(sub.N)
    lo, hi = sub.pairs[v]
    j = np.searchsorted(lo, y)
    if j < len(lo) and lo[j] == y:
        return int(hi[j])
    j = np.searchsorted(hi, y)
    return int(lo[j]) if (j < len(hi) and hi[j] == y) else y


def run_dynamic(inst, sub, U, beta, T, k, k0, q, params, seed, eta=None):
    """Returns the schedule (indices into F), averaged utilities, fairness value and dual bound.
    params: proposal parameters (t_ev, K, gamma, h), with gamma scaled for unit-norm rates."""
    rng = np.random.default_rng(seed)
    n = inst["n"]
    rho = float(U.max())
    eta = np.sqrt(8 * np.log(n) / T) / rho if eta is None else eta
    theta = np.full(n, 1.0 / n)
    y = 0                                    # F is sorted, so F[0] = 0 is the empty set
    xs, thetas = [], []
    rnorm = np.linalg.norm(inst["rates"])
    for t in range(T):
        w = n * theta[inst["owners"]] * inst["rates"]          # normalised weights, mean scale 1
        E = sub.energies(w)
        for _ in range(k0 if t == 0 else k):
            if rng.random() < 0.5:
                continue                     # laziness
            if rng.random() < q:
                tev, K, gam, h = params[rng.integers(len(params))]
                x = sub.sample_proposal(y, tev, K, gam * rnorm / np.linalg.norm(w), h, E, inst["colours"], rng)
            else:
                x = _local_move(sub, y, rng)
            if rng.random() < min(1.0, np.exp(beta * (E[x] - E[y]))):
                y = x
        xs.append(y); thetas.append(theta.copy())
        theta = theta * np.exp(-eta * U[y]); theta /= theta.sum()
    xs = np.array(xs); avg = U[xs].mean(axis=0)
    upper = min(float(np.max(U @ th)) for th in thetas)
    return dict(xs=xs, avg=avg, value=float(avg.min()), dual_upper=upper, thetas=np.array(thetas))
