# Fairness Is Advantage-Neutral — project broken into phases

Status as of 30 Sep 2026. Each phase lists: what we have, the one-paragraph plain statement, what is
genuinely ours versus already known (with the paper that knows it), what has been verified, and the gaps
to close before moving to the next phase. Citation keys refer to `paper/refs.bib`.

---

## The whole algorithm in plain words

**One sentence.** Fairness picks the prices; the computer only has to shop at those prices, and a quantum
device can only help with the shopping.

**Five steps, repeated (read it like Grover's "flip, then invert about the mean, repeat").**

1. **Priorities.** Every agent starts with the same priority weight.
2. **Prices.** Each link (resource) gets a price = its owner's priority × its rate.
3. **Shop.** Draw one feasible schedule from a lottery that favours expensive schedules
   (probability ∝ e^{β·total price}). We draw it by a random walk: *propose* a new schedule, *accept* it
   if it is more valuable, otherwise accept with probability e^{β·(change in value)}. The proposal can be
   a single local flip or a short quantum evolution that physically cannot create a conflict.
4. **Pay.** Whoever received little this round gets a higher priority next round (multiply by
   e^{−η·utility}); whoever received a lot gets a lower one.
5. **Repeat.** After T rounds the average allocation is within ε of the fairest lottery that exists.

The weight update in step 4 plays the same role for fairness that "inversion about the mean" plays in
Grover's search. Agents above average are pushed down and agents below average are pushed up. Each
round moves the system a little toward the fair point.

**The SVM analogy.** A support vector machine maximises the margin, and at the optimum only the
*support vectors* carry dual weight. Max-min fairness maximises the worst-off agent's utility, and at the
optimum only the *worst-off agents* carry priority weight (complementary slackness). The dual weights θ
are the "support agents" of the allocation.

**The theory in three lines.**
- *Separation (Thm 2):* fairest lottery = min over prices θ of [best weighted welfare at θ − Ψ*(θ)]. The
  fairness criterion only chooses which prices are allowed. The combinatorics only answers "best schedule
  at these prices".
- *Neutrality (Thm 4):* hence fair allocation is exactly as hard as weighted-welfare allocation, on any
  computer. Fairness cannot create or destroy a quantum advantage.
- *Tracking (Thm 9, 11):* prices move slowly, so the target lottery moves slowly, and a chain that is
  already warm needs about (1/γ)·log(1/ε) steps per round. The spectral gap γ of the sampler is therefore
  the only thing a quantum device can improve.

---

## Phase 1 — Fairness criteria (the normative layer)

**What we have.** Table 2 of the manuscript: 12 criteria (utilitarian, max-min, leximin, Lorenz, OWA,
Gini–Sen, Nash, α-fair/Atkinson, Kalai–Smorodinsky, need floors, ex-ante envy-freeness, mean−variance),
each with its dual price set. Appendix A derives them. `code/fairq/fairness.py` implements the measures.

**Plain statement.** Every reasonable fairness rule is a concave function of average utilities, and each
one corresponds to a set of allowed prices on agents.

**Genuine or known?** The individual criteria and their conjugates are textbook
(Atkinson 1970; Sen 1973; Weymark 1981; Yager 1988; Mo–Walrand 2000; Kalai–Smorodinsky 1975;
Rockafellar 1970). What is ours is the *unified price table*, especially the Pareto diagnostic: a
criterion is monotone ⇔ its prices are non-negative, which is why mean−variance can hurt the well-off.
Present it as an organising contribution, not as new mathematics.

**Verified.** Derivations checked by hand: Gini–Sen weights (2(n−k)+1)/n², α-fair conjugate α/(α−1)·Σθ^{(α−1)/α},
Nash conjugate 1+log θ, and the mean−variance gradient are all correct.

**Gaps to close before Phase 2.**
- [ ] Leximin and ex-ante envy-freeness are not a single concave Ψ. Leximin is lexicographic (a sequence
  of max-min LPs), and EF adds n² linear constraints. Say this explicitly in the table.
- [ ] Nash and α-fair have *unbounded* price sets (θ_i = 1/y_i → ∞ as y_i → 0). Phase 2's complexity
  bound needs bounded prices, so either assume a utility floor y_i ≥ δ or restrict to Lipschitz Ψ.
- [x] Notation fixed: the appendix used Φ where the text uses Ψ.

**Exit criterion.** One table, one lemma ("Ψ concave ⇒ allowed prices = dom Ψ*; monotone ⇔ prices ≥ 0"),
and every row either marked Lipschitz or given its floor assumption.

---

## Phase 2 — From criteria to prices (separation and advantage neutrality)

**What we have.** Thm 2 (separation), Cor 3 (zero representation cost), Thm 4 (advantage neutrality),
Prop 5 (deterministic max-min is NP-hard, ex-ante is not; Shapley–Folkman rounding), Prop 6 (degree of
max-min and Gini; concavity of sum-of-variances), Lemma 7 (Gibbs variational bound), Thm 8 (max-entropy
lottery is a Gibbs state).

**Plain statement.** Randomising removes the representation problem. Fairness becomes a price vector,
and the hard part is always "find the best schedule at given prices", which is the same problem you had
without fairness.

**Genuine or known?**
- Thm 2 is Fenchel duality plus Sion's minimax theorem (Rockafellar 1970; Sion 1958). The mathematics is
  standard. The *reading* of it ("normative vs combinatorial") is the contribution.
- Thm 4 for **max-min** was already shown by **Kawase & Sumita (AAAI 2020)**, who compute max-min fair
  lotteries by repeated "virtual welfare" maximisation. It also follows the Grötschel–Lovász–Schrijver
  equivalence of optimisation and separation. Our extension to *every* concave Ψ, and the statement about
  quantum advantage, is new. This is now cited and positioned in the paper.
- Prop 5's NP-hardness is a textbook PARTITION reduction (Karp 1972; Garey–Johnson 1979); the rounding
  uses Shapley–Folkman (Starr 1969).
- Lemma 7 and Thm 8 are classical (Jaynes 1957; Wainwright–Jordan 2008).

**Gaps to close before Phase 3 (these are the ones a strict reader will attack).**
- [ ] **Thm 4 overstates "polynomially equivalent".** The dual method needs poly(D, ρ, 1/ε) calls, which is
  pseudo-polynomial in 1/ε, and it needs a bounded, Lipschitz Ψ. Restate it as
  "ε-FAIR reduces to O(L²ρ²·log n/ε²) calls of an ε-approximate weighted best-response oracle", with the
  Lipschitz constant L of Ψ written in.
- [ ] **Approximate oracles.** A quantum or heuristic solver returns approximate best responses. State how
  oracle error δ enters the final accuracy (additively, ε + δ).
- [ ] **Prop 6 uses "d" without defining it.** It should be the multilinear degree of the valuations.
- [ ] **Thm 8's title ("the fairest lottery").** Maximum entropy chooses *among* lotteries with the same
  mean. That is a tie-break rule, not fairness. Rename it to "the maximum-entropy fair lottery is a Gibbs state".

**Exit criterion.** Thm 4 restated with explicit oracle complexity and error propagation, and a one-line
remark crediting Kawase–Sumita for the max-min case.

---

## Phase 3 — Implementing it: dynamic fair allocation (classical algorithm)

**What we have.** Eq. (1): Hedge/multiplicative weights on agents. Thm 9: regret guarantee
v* − log|F|/(βn) − ρε − ρ√(ln n/2T). Lemma 10: target drift. Thm 11: warm-start tracking,
k ≥ (1/2γ)·ln[(e^{2βnδ}(1+ε²)−1)/ε²], independent of |F|. Code: `code/fairq/dynamic.py`,
`code/dynamic_run.py`.

**Plain statement.** Raise priority for whoever is behind, resample a schedule from the tilted lottery,
repeat. Because priorities change slowly, the sampler never has to start from scratch.

**Genuine or known?**
- Hedge and its regret bound are standard (Freund–Schapire 1997; Arora–Hazan–Kale 2012; Cesa-Bianchi–Lugosi 2006).
- **The classical pipeline itself is known in wireless networking.** Adaptive CSMA
  (**Jiang & Walrand, IEEE/ACM ToN 2010**) samples independent sets from a Gibbs distribution whose weights
  are dual prices of a utility-maximisation problem. The **network adiabatic theorem
  (Rajagopalan–Shah–Shin, SIGMETRICS 2009)** shows that the chain tracks slowly varying weights.
  Our Thm 11 is an explicit, non-asymptotic, |F|-independent χ² version of that idea. Both papers are now
  cited, and the novelty sentence was changed accordingly.

**Verified.** Proofs of Thm 9, Lemma 10 and Thm 11 checked line by line. η = ρ⁻¹√(8 ln n/T) gives regret
ρ√(T ln n/2), and the χ² contraction (1−γ)^{2k} and TV ≤ ½√χ² are correct. We reproduced the numbers by
running `dynamic_run.py` again: v* = 1.4604, best deterministic 1.2266, slack 0.6197 (paper: 0.620),
classical 1.179 ± 0.139, quantum-mixed 1.327 ± 0.0075, SD ratio 18.6 ("about nineteen").

**Gaps to close before Phase 4.**
- [ ] **Close the loop with one end-to-end theorem.** Combine Thm 9 and Thm 11: after
  T = O(ρ² ln n/ε²) rounds of k = O(γ⁻¹·(βnδ + log 1/ε)) steps, the fairness error is ≤ 3ε. The total cost
  is then a single formula, total steps ≈ (ρ² ln n/ε²)·(1/γ)·log(1/ε). This is the "Grover-style" headline
  result the paper is missing. It uses only what is already proved.
- [ ] **The entropy term dominates the slack.** log|F|/(βn) = 0.54 of the 0.62 in the benchmark. Either
  anneal β upward over rounds, or state that β must scale like log|F|/(nε).
- [ ] **The ± values come from 5 repetitions.** Report them as mean ± SD (n = 5) and add the 5 raw values.

**Exit criterion.** One theorem that gives total sampler steps as a function of (ε, n, ρ, γ). Everything
else in Phase 4–5 is then about making γ large.

---

## Phase 4 — The quantum proposal (device layer)

**What we have.** H = γ_c·H_θ + h·H_PXP. Lemma 12: the PXP driver never leaves the independent-set subspace.
Lemma 13: non-adjacent terms commute, so colour classes are exact. Eq. (2): palindromic Trotter step.
Lemma 14: U^T = U, so the proposal is exactly symmetric at *any* Trotter depth. Thm 15: the lazy Metropolis
chain is exact, irreducible and aperiodic. Gate-level circuit with closed-form counts; OpenQASM export.
Code: `fairq/subspace.py`, `fairq/circuit.py`, `verify.py`, `circuit_spec.py`.

**Plain statement.** Let the device shuffle the schedule for a short time. The Rydberg blockade makes a
conflicting schedule physically impossible, and running the shuffle forwards-then-backwards
(palindromic) makes "a→b" exactly as likely as "b→a". So a classical coin flip on the value difference
keeps the answer exactly right.

**Genuine or known?**
- Quantum proposals with a classical Metropolis filter: **Layden et al. (Nature 2023)**; follow-ups
  Nakano et al. (PRR 2024), Ferguson–Wallden (PRR 2025), Christmann et al. (PRA 2025),
  Marshall et al. (arXiv 2026).
- Rydberg-native independent sets: Jaksch 2000; Lukin 2001; Pichler 2018; Ebadi (Science 2022).
  Weighted instances: de Oliveira (PRX Quantum 2025). Embeddings: Nguyen (PRX Quantum 2023).
- **Closest prior art for sampling:** Wild, Sels, Pichler, Zanoci, Lukin (PRL 2021; PRA 2021) already sample
  weighted independent-set Gibbs states on Rydberg arrays, by *adiabatic* preparation. Our difference is
  that the quantum part is only a proposal, so correctness never depends on adiabaticity. Now cited.
- **Genuinely ours:** the three-line argument "blockade ⇒ feasible, palindrome + real generators ⇒
  symmetric, Metropolis ⇒ exact", its gate-level realisation with exact counts, and its use inside a
  fairness loop whose target moves.

**Verified.** The shipped `verification.json` matches the paper exactly (3.1×10⁻¹⁵, 2.0×10⁻¹⁴,
1.7×10⁻¹⁶, 6.5×10⁻¹⁹, 1.1×10⁻¹⁶). A fresh rerun of `verify.py` on this machine gives the same picture,
with floating-point differences in the last digit (circuit vs subspace ≤ 4.1×10⁻¹⁵, asymmetry ≤ 2.2×10⁻¹⁶),
and closed-form counts match on all 5 instances. Quote "below 10⁻¹⁴" in the paper, not the exact
machine-dependent maxima.
`circuit_spec.py` reproduces 40 P / 400 X / 324 CCX / 76 CNOT / 76 H / 76 R_Z, 29 qubits, depth 445,
and 2020 CNOT after decomposition.

**Gaps to close before Phase 5.**
- [ ] **Noise robustness (the obvious referee question).** Decoherence breaks U^T = U. Add a lemma: if the
  proposal asymmetry is at most δ in total variation, the stationary distribution moves by
  O(δ/γ). This is a standard Markov-chain perturbation bound (verify and cite, e.g. Mitrophanov,
  J Appl Probab 2005). It turns the limitation into a quantitative statement.
- [ ] **Scars and revivals.** PXP has quantum many-body scars (Turner 2018) with near-periodic revivals from
  Néel-like states. Evolution times near a revival give near-identity proposals. The code already
  randomises t_ev ∈ [2, 10]; say why in one sentence.
- [ ] **Non-unit-disk graphs.** Erdős–Rényi and 3-regular conflict graphs are *not* native to 2D arrays and
  need gadgets (Nguyen 2023), which change the timescales (Bombieri et al., PRX Quantum 2025). The
  "no gate compilation needed" sentence should be restricted to unit-disk graphs.

**Exit criterion.** Lemma 12–14 plus a noise-perturbation lemma, and a clear statement of which graph
families are hardware-native.

---

## Phase 5 — Quantum advantage: evidence and limits (experiments)

**What we have.** Spectral gaps of 5 kernels on unit-disk / ER / 3-regular graphs, N = 8…20, 5 seeds,
β ∈ {1, 3}; fitted rates κ and ratios a = κ_C/κ_Q; exact warm-start tracking cost on 39 instances; one
repeated max-min run.

**Plain statement.** The quantum proposal makes the chain mix 3× to 1000× faster than local moves on every
instance tested. Whether the advantage *grows* with size depends on the graph family, and only
Erdős–Rényi graphs show a scaling gain.

**Verified (all against the shipped `code/results/*.json`).**

| Claim in paper | Value in results | Status |
|---|---|---|
| quantum/local gap ratio 3.4 to 1116 | 3.43 to 1115.5 | OK |
| geometric means 25.6 / 20.3 / 48.6 | 25.6 / 20.3 / 48.6 | OK |
| Table 3 κ and a values | all 6 rows match `summary.json` | OK |
| tracking 21.4 / 32.3 / 120.6 / 1910; 1.44×; 82% of 39 | 21.35 / 32.33 / 120.57 / 1909.8; 1.437; 0.8205 | OK |
| feasible sets "17 to 2635" | **14** to 2635 (udg N=8 seed 4) | **fixed in paper** |
| max-min run 1.179±0.139 vs 1.327±0.007, v*=1.460, slack 0.620 | reproduced exactly by rerun | OK |

**Honest reading (what a strict reviewer will say).**
- The scaling fits use N = 8…18/20 and 5 seeds. A ratio like a = 2.28 ± 0.76 is about 1.7σ from 1.
- The idealised *uniform* proposal beats the quantum one on 34 of 190 (instance, β) cases.
- On 3-regular graphs at β = 1, the **mixed** kernel (κ = 0.080) scales better than the best classical
  kernel (κ = 0.095). The paper currently reports only the pure quantum kernel there, so this is a
  positive result it leaves out.
- Orfi–Sels (PRA 2024) bound what any unital proposal can do in the worst case. The gain must come from
  graph structure, so the paper should say *which* structure.

**Gaps to close (in order of value per effort).**
1. [ ] **Right temperature regime.** Sampling independent sets is classically hard only above the
   uniqueness threshold λ_c(Δ) = (Δ−1)^{Δ−1}/(Δ−2)^Δ (Weitz 2006; Sly 2010; Mossel–Weitz–Wormald 2009);
   for Δ = 3, λ_c = 4. Sweep fugacity e^{β w} across λ_c and show the quantum/classical gap ratio jumps
   there. This is simple, principled and a clear figure.
2. [ ] **Bigger N.** Matrix-free Lanczos (Lanczos 1950) on the symmetrised kernel gets to N ≈ 30–36
   (|F| ~ 10⁵–10⁶) on one machine. Use 20 seeds and bootstrap confidence intervals for a.
3. [ ] **Stronger classical baselines.** Add the quantum-inspired proposal of Christmann et al. (PRA 2025)
   and parallel tempering (Hukushima–Nemoto 1996; as used for QeMCMC by Marshall et al. 2026). A reviewer
   will ask for both.
4. [ ] **Cost in time, not steps.** Price one quantum proposal in µs of Rydberg evolution plus readout,
   against ns per classical flip. This reports steps-to-solution and wall-clock side by side.
5. [ ] **Report the mixed kernel's 3-regular result.**

**Exit criterion.** One figure of gap ratio vs fugacity across λ_c, one table with bootstrap CIs at
N ≤ 36, and both new baselines.

---

## Phase 6 — Write-up and positioning

**Done in this pass.**
- Rebuilt the LaTeX source from the PDF (`paper/main.tex`) and rendered `paper/main.pdf`.
- Moved the bibliography to `paper/refs.bib`: 123 entries, every one cited, entries added in this pass
  marked `% ADDED`, DOIs throughout.
- Fixed 8 wrong or outdated references:
  - [24] was missing its first author, **Aziz**.
  - Orfi–Sels is now published in PRA 110:052414.
  - Caragiannis et al. is now published at AAMAS 2024.
  - Udekwe et al. is now published in J Infrastruct Preserv Resil 2026, with a fifth author.
  - Pathak et al.: full title, pages and DOI.
  - White et al.: "multicenter", pages 503–506, DOI.
  - Anschuetz et al.: new title in arXiv v2.
  - Akrami et al.: arXiv v3 title.
- Added citations wherever a claim needed one: the QAOA paper (Farhi–Goldstone–Gutmann 2014, the PDF that
  was in this folder), Metropolis/Hastings, Hoeffding, Hedge, Trotter/Lloyd, PXP (Fendley; Lesanovsky),
  blockade (Jaksch; Lukin), hardness of independent-set counting (Provan–Ball; Dyer–Frieze–Jerrum; Sly;
  Mossel–Weitz–Wormald), graph models, and software (NumPy, SciPy, NetworkX).
- **Rewrote Related Work** as four threads, each ending with what it leaves open, plus a comparison table
  against the five closest precedents (Kawase–Sumita, Jiang–Walrand, Rajagopalan–Shah–Shin, Layden, Wild).
- **Corrected the novelty claim.** "Fairness criteria have not previously been placed in this sampling
  framework" is not true (Jiang–Walrand). It now reads that a quantum proposal has not been used inside a
  fairness-driven Gibbs sampler.
- Fixed Φ/Ψ notation and the feasible-set range (14, not 17).

**Still to do.**
- [ ] Apply the Phase 1–4 restatements (Thm 4 oracle complexity, the end-to-end theorem, the noise lemma).
- [ ] Six references are still preprints: Akrami, Anschuetz, Babaioff, Bach, Choi, Farhi 2020.
  Pichler 2018, QAOA 2014 and OpenQASM 2017 have no journal versions. Before submission, check for
  published versions again, or for Farhi 2020 cite the peer-reviewed Basso et al. (FOCS 2022) alongside it.
- [ ] Kaneko–Nakamura page range: RePEc gives 423–435, the Econometric Society gives 423–436. Check the scan.

---

## Research gaps worth a contribution (ranked: simple, principled, not already done)

1. **End-to-end complexity of fair allocation.** Total sampler steps = Õ(ρ² ln n/ε² · 1/γ). This is one
   formula that makes "the spectral gap is everything" a theorem, not a slogan. It needs no new
   experiments, only Phase 3.
2. **Uniqueness-threshold experiment.** Does the quantum advantage switch on exactly where classical
   sampling becomes hard (λ > λ_c)? It is one clean figure. If yes, it explains the family dependence (dense
   ER graphs reach the hard regime at lower β), and it answers Orfi–Sels's "structure is needed" with a
   named structure.
3. **Noise-robust exactness.** A bound on the stationary bias of O(asymmetry/γ). It turns "exact only on ideal
   hardware" into a hardware requirement.
4. **Quantum-CSMA.** Jiang–Walrand's adaptive CSMA with the PXP proposal replacing local backoff. This is
   the natural systems contribution, and the wireless community already accepts Gibbs-sampling scheduling.
5. **Beyond independent sets.** Matchings and assignment with XY mixers (Wang et al. 2020): the same three
   lemmas with a different driver. It shows the construction is general.

---

## Folder map

```
qca/
  PROJECT_PHASES.md          this file
  paper/  main.tex, refs.bib, figures/, main.pdf (rendered)
  code/   fairq package, experiment scripts, shipped results/ (from fairq_v1.0.0)
  old/    original manuscript PDF, files.zip, QAOA paper (1411.4028v1.pdf)
```
Build the PDF: `tectonic paper/main.tex` (or pdflatex + bibtex).
