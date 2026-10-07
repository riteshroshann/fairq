"""Figures and summary statistics of the manuscript -> results/fig_gaps.pdf, fig_tracking.pdf, summary.json"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.serif": ["TeX Gyre Pagella", "Palatino", "DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})
FAM = [("udg", "unit-disk"), ("er", r"Erdős–Rényi $G(N,0.3)$"), ("reg3", "random 3-regular")]
STY = {"local": ("black", "o", "local"), "swap": ("0.45", "v", "local + swap"),
       "uniform": ("0.7", "s", "uniform on $F$"), "quantum": ("#1f4e9a", "^", "quantum"),
       "mixed": ("#9a1f1f", "D", "mixed")}
gm = lambda a: float(np.exp(np.mean(np.log(a))))

summary = {"exponents": {}, "ratios": {}, "tracking": {}}
fig, axes = plt.subplots(2, 3, figsize=(7.3, 4.6))
for j, (fam, title) in enumerate(FAM):
    rows = json.load(open(f"results/gaps_{fam}.json"))
    ex = json.load(open(f"results/exponents_{fam}.json"))
    summary["exponents"][fam] = ex
    for i, beta in enumerate((1.0, 3.0)):
        ax = axes[i, j]
        R = [r for r in rows if r["beta"] == beta]
        Ns = sorted({r["N"] for r in R})
        for m, (c, mk, lab) in STY.items():
            ax.scatter([r["N"] for r in R], [r[m] for r in R], color=c, marker=mk, s=6, alpha=0.3, linewidths=0)
            ax.plot(Ns, [gm([r[m] for r in R if r["N"] == N]) for N in Ns], color=c, marker=mk, ms=3, lw=1, label=lab)
        ax.set_yscale("log")
        if i == 0: ax.set_title(title, fontsize=9)
        if i == 1: ax.set_xlabel("$N$")
        if j == 0: ax.set_ylabel(rf"gap $\gamma$,  $\beta={beta:g}$")
        e = ex[str(beta)]
        best = min(("local", "swap", "uniform"), key=lambda m: e[m]["kappa"])
        a = e[best]["kappa"] / e["quantum"]["kappa"]
        se = a * np.hypot(e[best]["se"] / e[best]["kappa"], e["quantum"]["se"] / e["quantum"]["kappa"])
        summary["ratios"][f"{fam}_beta{beta:g}"] = dict(best_classical=best, ratio=float(a), se=float(se))
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, ncol=5, loc="lower center", fontsize=8)
fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig("results/fig_gaps.pdf")

tr = json.load(open("results/tracking.json"))
fig, axes = plt.subplots(1, 3, figsize=(7.3, 2.6), sharey=True)
for j, (fam, title) in enumerate(FAM):
    R = [r for r in tr if r["family"] == fam]; Ns = sorted({r["N"] for r in R})
    summary["tracking"][fam] = {m: [gm([r[f"k_{m}"] for r in R if r["N"] == N]) for N in Ns] for m in STY}
    summary["tracking"][fam]["N"] = Ns
    for m, (c, mk, lab) in STY.items():
        axes[j].plot(Ns, summary["tracking"][fam][m], color=c, marker=mk, ms=3, lw=1, label=lab)
    axes[j].set_yscale("log"); axes[j].set_title(title, fontsize=9); axes[j].set_xlabel("$N$")
axes[0].set_ylabel("steps per checkpoint")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, ncol=5, loc="lower center", fontsize=8)
fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig("results/fig_tracking.pdf")

allr = [r for r in tr]
summary["tracking_overall"] = {m: gm([r[f"k_{m}"] for r in allr]) for m in STY}
summary["tracking_best_classical_over_quantum"] = gm([min(r["k_local"], r["k_swap"], r["k_uniform"]) / r["k_quantum"] for r in allr])
summary["tracking_quantum_best_fraction"] = float(np.mean([r["k_quantum"] <= min(r["k_local"], r["k_swap"], r["k_uniform"]) for r in allr]))
json.dump(summary, open("results/summary.json", "w"), indent=1)
print(json.dumps({k: summary[k] for k in ("ratios", "tracking_overall", "tracking_best_classical_over_quantum", "tracking_quantum_best_fraction")}, indent=1))
for fam in summary["tracking"]:
    print(fam, {m: [round(v, 1) for v in summary["tracking"][fam][m]] for m in summary["tracking"][fam]})
