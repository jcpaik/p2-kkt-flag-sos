# The $h_2$-weighted experiment: what was done and why it gained $161\times$

This note documents one working session: reproducing the SDPA-GMP
baseline, reformulating the target as a *rational certificate* problem,
and the sequence of controlled experiments that took the degree-14
bound from $-4.5\times10^{-3}$ to $-2.8\times10^{-5}$ — a $161\times$
improvement at the **same degree, arity, and label algebra**.  It also
explains the one negative control that makes the mechanism of the gain
unambiguous.  Current status of the follow-ups lives in
[PLAN](../PLAN.md) §4–§5; the theory write-up is
[Exact zero program](EXACT_ZERO_PROGRAM.md).

## 1. Baseline reproduction

SDPA-GMP (an arbitrary-precision interior-point semidefinite solver;
double precision provably stalls on these problems) was built from
source (`nakatamaho/sdpa-gmp`, GMP via Homebrew; two build fixes:
`-Wno-int-conversion` for the bundled SPOOLES C code, `-std=c++17` for
the C++ part).  The pruned nine-toggle hierarchy at degree 14,

```sh
python3 sos_search.py --export-sdpa deg14_pruned.dat-s \
  --degree 14 --no-pointwise-sos \
  --harmonics --three-point-flags --four-point-flags --two-root-flags \
  --gradient --potential --hessian --global-tangent-gaps --rank-relations
sdpa_gmp -ds deg14_pruned.dat-s -o deg14_pruned.result -p param_200bit.sdpa
```

reproduced the published bound to 10 significant digits:
$-4.48560259\times10^{-3}$ (m = 398, 46 iterations, 6 min, feasibility
errors $\sim10^{-26}$).  The binary, parameter files and all problem
files are kept in `sdpa_runs/`.

## 2. Why more degree cannot reach zero, and what can

The published measurements contain a structural diagnosis:

* the bounds decay like $\exp(-c\,d^2)$ — they approach $0$ but at no
  finite degree reach it;
* the minimal certificate trace for $E+\varepsilon\ge0$ diverges like
  $\operatorname{tr}C(\varepsilon)\approx1.07/\varepsilon$ (legacy scale; the pole order is scale-free) — a
  **simple pole**.

Together these say the sharp certificate does not exist in the
finite-degree cone: the $\varepsilon$-certificates escape to infinity
along a fixed recession direction (classical sum-of-squares
**non-attainment**, the Motzkin phenomenon).  The classical repair is
not more degree but a *multiplier*: find a nonnegative quantity $Q$
vanishing to first order along the escape so that $Q\cdot E$ admits a
certificate even though $E$ does not.  A simple pole predicts a single
factor suffices.

The escape variable was already identified in the moment diagnostics:
the optimal pseudo-moment pays the spin-2 harmonic energy

$$h_2=\frac{3p_2-1}{2}=\sum_m|\hat\mu_{2m}|^2\ \ge\ 0$$

linearly while earning contraction violations of order $\sqrt{h_2}$.
So the canonical multiplier is $Q=h_2$.

**Reduction lemma (nothing is lost).**  If $h_2(\mu)E(\mu)\ge0$ for
every antipodal probability measure, then $E\ge0$ for every such
measure: on $\{h_2>0\}$ divide; an isotropic $\mu$ is the weak-\* limit
of the anisotropic $\mu_t=(1-t)\mu+t\,\delta_{\pm e}$ and $E$ is
continuous.  So certifying the weighted target proves the full
conjecture — no separate isotropic branch remains.

## 3. Making $h_2E$ a target the machinery accepts

$h_2E$ is *polynomial* in the existing label algebra:

$$h_2E=\tfrac23-12p_2+24p_4-16p_6
+30\,p_2p_2-72\,p_2p_4+48\,p_2p_6 ,$$

where a product $p_2p_j$ is the expectation of a *disconnected*
four-sample graph — two independent pairs — which the moment reducer
already represents as a `("product", ("pair",2), ("pair",j))` label.
Implemented as `--h2-weighted-target` (a different rational objective
vector; every other part of the export pipeline unchanged).

## 4. The experiments, and the anatomy of the $161\times$ gain

All runs: degree 14, pruned nine-toggle base, SDPA-GMP 200-bit,
`epsilonStar 1e-25`, primal = dual to $\ge13$ digits.

| # | configuration | m | bound |
|---|---|---:|---:|
| 0 | unweighted $E$ (baseline) | 398 | $-4.4856\times10^{-3}$ |
| 1 | $h_2E$, same blocks | 398 | $-2.0840\times10^{-3}$ |
| 2 | $h_2E$ + `--h2-localized-flags` | 414 | $-1.2032\times10^{-3}$ |
| 3 | $h_2E$ + $p_2\times$(all families) — control | 790 | $-2.0840\times10^{-3}$ |
| 4 | $h_2E$ + $h_2\times$(all families) | 790 | $\mathbf{-2.7909\times10^{-5}}$ |

### Why run 1 only doubled the bound

The relaxation *linearizes* products: the label $y_{p_2p_j}$ is a free
variable, not the product $y_{p_2}y_{p_j}$.  With no constraints tying
the product labels to their factors, the pseudo-moment simply decouples
the weight from the target, and most of the intended cancellation never
happens.  (The one coupling already present — the `empty_type_flag`
block, which is the moment matrix over $\{1,p_2,p_4,p_6\}$ — supplies
the one-sided inequality $y_{p_2p_2}\ge y_{p_2}^2$, which is why run 1
still improves at all.)

### The negative control (run 3), and what it proves

Multiplying every PSD block $B\succeq0$ by $p_2$ is *valid* ($p_2\ge
\tfrac13>0$ for every measure, and disconnected samples factor), and it
doubles the problem to m = 790 — yet the bound is **identical to run 1
to eleven digits**.  Reason: a $p_2$-multiplied block touches *only*
fresh product labels $y_{p_2\cdot\ell}$.  Those variables appear in no
other constraint and not in the objective, so the pseudo-moment can
satisfy the new block trivially (e.g. by inflating its diagonal)
without constraining anything it cares about.  A localizing block is
worthless unless its entries **mix** the new product labels with the
old base labels inside one PSD constraint.

### The fix (run 4): localize by $h_2$, not by $p_2$

Since $h_2=\tfrac{3p_2-1}{2}$, the $h_2$-multiplied copy of a block
with coefficient matrices $\{A_\ell\}$ has entries

$$h_2\times A_\ell\;\longmapsto\;
\tfrac32\,A_\ell\ \text{on}\ y_{p_2\cdot\ell}
\;\;-\;\tfrac12\,A_\ell\ \text{on}\ y_{\ell},$$

and the $-\tfrac12$ part is the whole point: every localized block now
couples each product label to its base label with a definite sign.  A
pseudo-moment that tries to depress $y_{p_2\cdot\ell}$ below the true
product is caught by the PSD constraint through the $y_\ell$ term it
can no longer decouple from.  `--h2-localized-all` adds the
$h_2$-multiplied copy of **every** PSD family and the $p_2$-multiplied
copy of every scalar equality family (gradient, potential, rank; for
equalities the constant part is redundant, since the unmultiplied
relation is already present).

On the certificate (dual) side the same structure is what the
non-attainment analysis demands.  The escaping $\varepsilon$-certificates
needed a multiplier growing like $1/\varepsilon$ to fight the
$\sqrt{h_2}$ leak; in the weighted problem that role is played by the
*fixed* $h_2$ factor carried by the localized blocks, so certificates
of the form $\sigma_0+h_2\sigma_1$ (with $\sigma_i$ in the ordinary
flag-square cone) become expressible at finite size.  Formally, run 4
optimizes over the degree-14 truncation of the quadratic module
generated by $h_2$, where runs 0–1 only had the trivial module.

The result, $-2.7909\times10^{-5}$, is

* $161\times$ better than the unweighted bound at the same degree,
  arity, and label set (run 0 → run 4);
* $2.7\times$ better than the *unweighted degree-18* bound
  ($-7.5595\times10^{-5}$) with $30\%$ fewer variables (790 vs 1089) —
  i.e. the reformulation is worth more than two full degree steps;
* attributable, by the run-3 control, entirely to the $-\tfrac12$
  coupling term — not to problem size (identical m), not to the weight
  alone (run 1), not to a subset of families (run 2).

### The attainment measurement

The decisive question for an *exact* zero is whether the weighted
problem attains its certificate.  The first data point (quoted in the
legacy $(3/16)$-scale in which the selector problems were exported;
ratios are scale-free): minimal certificate trace for the legacy-scaled
$h_2E\ge-10^{-4}$ (via `sdpa_selector.py`) is $330.1$, where the
unweighted pole law $1.07/\varepsilon$ predicts $\sim10{,}700$ — a
$32\times$ collapse of the pole coefficient.  The follow-up $\varepsilon$-sweep (recorded in
[PLAN](../PLAN.md) §4) shows the trace still grows between
$\varepsilon=10^{-4}$ and $10^{-5}$: one factor of $h_2$ tames but does
not kill the pole, and the surviving escape is the pure $g_4$
($\hat K_4<0$) mode — which is what the operator-gap route in PLAN §4
now targets.

## 5. Takeaways

1. **Weighted targets are cheap and strong.**  Multiplying the target
   by a nonnegative moment quantity costs nothing (a new objective
   vector) and, with the *matching localized module*, moved the
   hierarchy further than two degree steps.  The reduction lemma makes
   the reformulation logically lossless.
2. **Localizers must mix product and base labels.**  Any multiplier
   with zero constant term ($p_2$, or the raw product $h_2\times$flag
   labels of the earlier inert experiment) produces vacuous blocks
   under linearization.  The affine part of the multiplier is what
   gives the relaxation teeth.  This is a general design rule for
   flag-algebra localization, worth remembering for the operator-gap
   blocks.
3. **The pole is measurable and informative.**  The
   $\varepsilon$-trace law cleanly separates "bound is small" from
   "certificate exists": run 4's bound looks excellent, but the trace
   law says degree 14 still does not attain — and *which* direction
   escapes ($g_4$) tells you the next constraint to add.

New code: `--h2-weighted-target`, `--h2-localized-all` in
`sos_search.py`.  Artifacts: `sdpa_runs/` (solver binary, parameter
files, problems, results).
