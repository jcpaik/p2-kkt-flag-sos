# The $e_5(I-A_2)$ second weight: exact expansion, implementation, first measurements

*2026-08-18.  Agent W deliverable for the
[multi-weight program](MULTI_WEIGHT_PROGRAM.md) §4/§6: the label-algebra
implementation of the face-vanishing multiplier
$e_5(B)=\det(I-A_2)$, $B=I-A_2$, matched to the Jensen-saturation
escape of [Unprojected escape note](UNPROJECTED_ESCAPE_NOTE.md) §5.
Code: `sos_search.py` (`gap_elementary_vector`, `--gap-cut-e5`,
`--e5-weight`, `--e5-localized-harmonics`); exact checks:
`sdpa_runs/e5_label_check.py`, `sdpa_runs/e5_perturbations.py`,
regression tests in `test_sos_search.py`.*

## 1. Inventory (W1): what the label algebra already had, what was missing

**The reducer was never the obstruction.**  `reduce_graph_matrix` /
`canonical_connected_label` canonicalize connected multigraphs on *any*
number of vertices (minimum over all vertex permutations), emitting
`graph_5`, `graph_6`, … labels; ONB, pole–equator, and (new) uniform
evaluation are equally generic.  The production weighted export
`deg14_h2w_h2all.dat-s` (nine toggles + `--h2-weighted-target
--h2-localized-all`, **without** `--max-flag-arity`) already contains 18
`graph_5` and 21 `graph_6` labels — all of them from the unrooted
`spin2_flag` Gram block (two harmonic vertices with leaves; verified by
regenerating that block's label set), and **none** of them the 5-cycle
labels of $\operatorname{tr}A_2^5$.

**Five-sample flag squares exist and reach the cycles.**  With
`--max-flag-arity 5 --max-root-factor-degree 2` (the README's strongest
unweighted config), the `weighted_flag_5_*` family (3 roots + 2 leaves,
root-root degree $\le2$ per flag) generates 554 distinct `graph_5`
labels at degree 14, and **all 11** five-sample labels of $e_5$ are
among them (measured; §3).  A 5-cycle embeds as a two-leaf double star
with one root–root edge; the doubled 5-cycle of $(\operatorname{tr}P)^2$
needs total root-root degree 2 = the current budget.  So no new flag
family is needed for coverage; the $A_2$-word blocks
$E[v^{\mathsf T}(I-A_2)v]$, $v=A_2^d[\text{spin-2 flag}]$ remain the
*escalation* if the arity-5 family proves too weak in practice, not a
prerequisite.

**What was genuinely missing** (now implemented): (i) the exact rational
label expansion of $e_5$ itself; (ii) an objective/cut path that puts
$e_5$-labels into the exported problem; (iii) 1×1 coverage blocks so the
new label sector cannot open an uncontrolled dual direction.

## 2. Exact expansion (W2)

`gap_power_trace_vector(k)` computes
$\operatorname{tr}A_2^k=E[\chi_2(\rho_{x_1}\cdots\rho_{x_k})]$,
$\chi_2(R)=(\operatorname{tr}R)^2-\operatorname{tr}R-1$, by expanding
$\operatorname{tr}(\rho_{x_1}\cdots\rho_{x_k})
=\sum_{S\subseteq[k]}(-1)^{k-|S|}2^{|S|}\operatorname{cyc}(S)$ and
reducing every resulting multigraph through the standard reducer.
`a2_elementary_vector(k)` applies Newton's identities (label-vector
products = disconnected label products, `multiply_label_vectors`), and
`gap_elementary_vector(k)` returns
$e_k(B)=\sum_j(-1)^j\binom{5-j}{k-j}e_j(A_2)$.  All arithmetic is
`Fraction`-exact.

**The $e_5$ vector**: 35 labels — 1 constant, 2 pair, 4 triangle, 6
`graph_4`, **11 `graph_5`** (the genuinely five-sample cycle sector,
coefficients like $-512$, $1280$, $-1024/5$), 3 pair×pair and 8
pair×triangle products.  Every connected label has per-vertex degree
$\le4$: the obstruction was label arity only, never polynomial degree.

**Machine checks (all exact):**

| check | result |
|---|---|
| `gap_elementary_vector(2,3,4)` vs the audited `--gap-scalar-cuts` tables | identical dictionaries |
| pole–equator family, **symbolic** in $(w,\hat\nu(2),\hat\nu(4),\hat\nu(6))$ via the label expansion (`e5_label_check.py` V1) | identically $0$ |
| pole + regular equators (orders 3,4,5,6) and ONB, exact evaluators | $0$ |
| uniform measure (`uniform_label_value`, new exact Wick evaluator) | $(4/5)^5=1024/3125$ |
| 3-atom rational measures (three weight sets), expansion vs direct $5\times5$ $\det(I-A_2)$ over $\mathbb Q$ | equal (test + standalone script) |
| $e_4$ on the same family (V3) | $(1-r^2)(1-w^2)^2\,(\cdot)\neq0$: $e_5$ is the *first* vanishing invariant |

The symbolic family check V1 doubles as a validation of the `graph_5`
evaluation path, since the family evaluator recomputes every label from
its raw multigraph.

## 3. Implementation (W3)

* **Labels**: no reducer change needed (§1).  New exact evaluator
  `uniform_label_value` (S² Wick pairing, vertex-by-vertex elimination).
* **Cut** `--gap-cut-e5`: $e_5(B)\ge0$ as a 1×1 block (35 labels).
  Validity: all-measures (eigenvalues of $B$ in $[0,2]$).
* **Weighted target** `--e5-weight KAPPA` (rational $\kappa\ge0$,
  requires `--h2-weighted-target`): target
  $(h_2+\kappa e_5(B))\,E$.  The $e_5E$ part is the disconnected
  product of the 35 $e_5$ labels with the four energy terms
  ($-4/3,20p_2,-48p_4,32p_6$).  Reduction lemma for this weight:
  [Multi-weight program](MULTI_WEIGHT_PROGRAM.md) Fact 3.  The target
  is now assembled in exact rationals and floated once; a regression
  test verifies every coefficient survives `rationalize_float`
  round-trip for $\kappa\in\{1/4,1,4\}$.
* **Coverage** `--e5-localized-harmonics` (implied by `--e5-weight`,
  as is the cut) — designed by exact recession analysis after a first
  1×1-only attempt measured *dual unbounded*:
  * `gap_cut_e5_upper`: $(4/5)^5-e_5\ge0$ — AM–GM on the five
    nonnegative eigenvalues of $B$ ($\operatorname{tr}B=4$
    identically), with **equality exactly at the uniform measure**
    (test-asserted).  Its recession forces $e_5(r)=0$ on any dual ray.
  * `e5loc_hankel_even/odd`: $e_5\cdot E[v(t)v(t)^{\mathsf T}]\succeq0$,
    $t=X\!\cdot\!Y$, $v=(1,t^2)$ resp. $(t,t^3)$ — a nonnegative scalar
    invariant times a pair moment matrix (disconnected samples
    factor).  Given $e_5(r)=0$, the two Hankels force
    $e_5p_2(r)=e_5p_4(r)=0$, leaving only the $+32\kappa\,e_5p_6$
    objective direction ($\ge0$): the $\kappa e_5E$ part cannot open a
    recession ray.  Confirmed operationally: with only 1×1 combos the
    $\kappa$-target dual is unbounded at every $\kappa$; with the
    Hankel+upper-cut module it is bounded (§4).
  * `e5loc_harmonic_d`, `e5comp_harmonic_d` ($d=2,4,6$):
    $e_5\cdot E[P_d]\ge0$ and $((4/5)^5-e_5)\cdot E[P_d]\ge0$.

  Together with the cut these contain every label of $\kappa e_5E$
  (test-asserted).  Under `--h2-localized-all` every e5 block acquires
  an $h_2$-multiplied copy automatically ($h_2\cdot e_5\cdot(\cdot)\ge0$,
  same validity).
* **Exporter hardening**: the exact image-elimination in
  `export_sdpa_problem` now tracks the objective residual of every
  *dropped* (image-dependent) direction.  A dropped direction whose
  objective is not the matching combination of kept objectives means
  the true dual is unbounded and the exported file would hide it; this
  is now counted (`dropped_objective_inconsistent`) and warned on
  stderr.  All e5 exports below have count 0.
* **Validity ledger**: every new block is all-measures valid
  (operator bound + products of nonnegative invariants); nothing here
  is KKT-only.  The composition caveat of PLAN §5 is untouched.

## 4. Measurements (W4; MOSEK double, indicative only)

All-measures weighted cone, degree 14, flags of
[Weighted-E1 note](WEIGHTED_E1_NOTE.md) §7 (`--h2-weighted-target
--h2-localized-all --harmonics --three/four-point --two-root
--rank-relations --dual --scale-constraints`).  Baseline reproduced
exactly: $-1.9197930\times10^{-4}$.

| configuration | bound (double) | status |
|---|---:|---|
| baseline $h_2E$ (reproduced) | $-1.91979\times10^{-4}$ | optimal |
| + e5 cut module, target unchanged | $\mathbf{-1.83199\times10^{-4}}$ | optimal |
| `--e5-weight 1/4` (1×1 coverage only, first attempt) | unbounded | — |
| `--e5-weight 1/4` (Hankel + AM–GM module) | $-1.7649\times10^{-2}$ | optimal |
| `--e5-weight 1` (same) | $-7.0237\times10^{-2}$ | optimal |
| `--e5-weight 4` (same) | — | MOSEK numerical failure |
| control + arity-5 weighted flags ($\kappa=0$) | $-1.91752\times10^{-4}$ | optimal |
| `--e5-weight 1/4` + arity-5 weighted flags | $-1.37225\times10^{-2}$ | optimal |

Readings: (i) the **cut module alone improves the all-measures double
bound by 4.6%** — the first *scalar* invariant to move this bound at
all (the $e_2,e_3,e_4$ cuts moved only the projected problem's ray
structure); (ii) the $\kappa$-weighted bounds are $\approx-0.07\kappa$:
the coverage module closes the recession but cannot yet *control* the
$e_5E$ mass; (iii) adding the arity-5 weighted flags (who do carry all
11 `graph_5` labels) improves the $\kappa=1/4$ bound by 22% while
leaving the $h_2E$ control bound essentially unchanged
($-1.9175\times10^{-4}$) — the `graph_5` sector is doing real work on
the weighted target, but degree-14 blocks are far from taming
$\kappa e_5E$.  This is the $e_5$-sector analogue of $h_2E$ before
`--h2-localized-all` existed: the missing object is the
**$e_5$-localized module** ($e_5\times$flag blocks; label count
$\approx35\times$ per localized family — a degree-16-scale job), not
more scalar cuts.

GMP-ready exports (exact rational, NOT solved here; queued for the
orchestrator; both have `dropped_objective_inconsistent = 0`):

| file | m | labels | shift | note |
|---|---:|---:|---|---|
| `deg14_h2w_h2all_e5cut.dat-s` | 797 | 1456 | $+2/3$ | KKT-inclusive weighted problem + full e5 cut module, target $h_2E$ — the direct test against the $-2.7909\times10^{-5}$ wall; the deep-escape pairing (below) predicts the cut bites below $\varepsilon\sim5\times10^{-5}$ |
| `deg14_e5w_k1d4.dat-s` | 797 | 1456 | $-58/9$ | $(h_2+\tfrac14e_5)E$ target, same blocks — *reference only* (double bound $-1.76\times10^{-2}$; not competitive until an $e_5$-localized module exists) |

**Escape pairing (partial, $\le4$-sample shadow).**  The 24
$\le4$-sample labels of $e_5$ paired against stored selector growth
directions $D$ (the 11 `graph_5` labels pair as 0 — absent from those
problems — so these numbers are the *visible shadow* only,
`sdpa_runs/e5_invariant_expansion.json`):

| interval (problem) | $\langle e_5^{\le4},D\rangle$ | $\langle h_2,D\rangle$ |
|---|---:|---:|
| $10^{-3}\to10^{-4}$ (control, `D_e3e4`) | $-6.37\times10^{5}$ | $+11.9$ |
| $10^{-4}\to5.33\times10^{-5}$ (control, `D_e4e5`) | $+2.36\times10^{4}$ | $+12.0$ |
| toep2 residual (`D_toep2`) | $-1.18\times10^{4}$ | $-0.80$ |
| **deep toep3 $5.33\times10^{-5}\to2\times10^{-5}$** (`D_toep3_deep`, new) | $\mathbf{-1.26\times10^{5}}$ | — |

By the sign rule the deep-escape pairing is **cut-signed**: the visible
shadow of $e_5\ge0$ excludes the deep escape direction.  (The weight
role rests on the face geometry, §5, not on this pairing.)

## 5. Does $e_5$ grow along the saturation face? (W5)

`sdpa_runs/e5_perturbations.py`: exact sympy series of
$\det(I-A_2)$ for model families leaving the pole–equator stratum
(the stratum itself is $e_5$-flat in *all* directions —
$\det\equiv0$ identically in $(w,\hat\nu(2),\hat\nu(4))$, so mass
reshuffling on the stratum, including the equator mode-4 direction of
`e5_face_check.py` V4, never registers):

| family (off-stratum motion) | leading order of $e_5$ |
|---|---|
| F1 split latitudes: equator $\to$ Haar rings at $\pm u$ | $\dfrac{256}{81}u^2$ |
| F2 rings with mode-4 density $r$ | $\dfrac{256}{81}(1-r^2)u^2$; at the ONB corner $r=\pm1$: $+\dfrac{1024}{81}u^4$ |
| F3 equator circle tilted by $t$ | $\dfrac{64}{81}t^2$ |
| F4 pole mass spread to a polar ring (angle $u$) | $\dfrac{128}{81}u^2$ |
| F5 control: any on-stratum $w$ | identically $0$ |

**Verdict.**  $e_5$ is strictly positive-quadratic transverse to the
stratum in every model direction tested, including the
fiber-deterministic motions (F1 split latitudes is exactly a
deterministic-$|z|$-fiber deformation; F4 is the pole-side analogue),
and degenerates to quartic-but-positive only at the isolated mode-4
corner $r=\pm1$.  Combined with the cut-signed deep-escape pairing
(§4), $e_5$ acts on the surviving escape from both sides: as a valid
cut on its visible shadow today, and as a weight whose zero set
excludes exactly the saturation face once the `graph_5` sector is in
the problem.  The operational test of requirement (d) is the exported
$\kappa$-grid on the GMP queue.

## 6. Files

* Code: `sos_search.py` — `multiply_label_vectors`,
  `rho_word_trace_terms`, `gap_power_trace_vector`,
  `a2_elementary_vector`, `gap_elementary_vector`,
  `uniform_label_value`, `exact_harmonic_pair_vector`, flags
  `--gap-cut-e5 --e5-localized-harmonics --e5-weight`, exporter
  objective-residual tracking.
* Checks: `sdpa_runs/e5_label_check.py` (symbolic family = 0, uniform,
  first-vanishing), `sdpa_runs/e5_perturbations.py` (growth orders),
  7 new tests in `test_sos_search.py` (suite green).
* Data: `sdpa_runs/e5_invariant_expansion.json` (exact expansion, full
  and $\le4$ shadow), `sdpa_runs/fingerprint_D_toep3_deep.json` (deep
  growth direction), `sdpa_runs/mosek_e5_*.json` (double-precision
  measurements), exports listed in §4.
