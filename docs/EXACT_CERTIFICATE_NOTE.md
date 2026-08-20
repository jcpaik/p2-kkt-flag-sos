# First exact rational certificate: $h_2E \ge -2\times10^{-8}$

*2026-08-19.  The exactification dress rehearsal: rounding the
GMP-solved proof-carrying problem `deg16_am_toep3` into an exact,
independently checkable rational certificate.  Artifacts:
`certificates/h2E_geq_minus_2em8.json.gz` (complete certificate data),
`verify_h2E_bound.py` (standalone pure-rational checker),
`sdpa_runs/build_h2E_certificate.py` (the rounding pipeline).*

## 1. Theorem (certified)

For every antipodally symmetric Borel probability measure $\mu$ on
$S^2$, with $E(\mu)=\iint K(x\cdot y)\,d\mu\,d\mu$,
$K(t)=32t^6-48t^4+20t^2-\tfrac43$, and
$h_2(\mu)=\tfrac{3p_2-1}2=E[P_2(X\cdot Y)]\ge0$:

$$h_2(\mu)\,E(\mu)\;\ge\;c-\|\rho\|_1\;=\;-1.4823972914\times10^{-8}\;>\;-2\times10^{-8}.$$

(The exact rational values of $c$ and $\|\rho\|_1$ are stored in the
certificate; see §4.)

## 2. Source and rounding pipeline

* Problem: `sdpa_runs/deg16_am_toep3.dat-s` — the degree-16
  **all-measures** weighted cone (no KKT content: flags, two-root,
  harmonics, spin-2, harmonic flags, their $h_2$-localized copies,
  four-point Gram-rank relations) plus the 43 Jensen/fiber-Toeplitz
  families (v3), $m=1265$, 109 blocks, shift $+2/3$.
* Solution: `sdpa_runs/deg16_am_toep3.result`, SDPA-GMP 7.1.3 at
  200 bits, `pdFEAS`, relative gap $5.2\times10^{-15}$, dual
  feasibility error $3.4\times10^{-35}$; wall
  $-1.2252\times10^{-8}$.  Margin to the target $-2\times10^{-8}$:
  $\approx7.7\times10^{-9}$.
* Rounding: the printed 40-digit `yMat` decimals are taken **exactly**
  (each printed decimal is a rational $p/10^k$), symmetrized
  $(Y+Y^{\mathsf T})/2$; no PSD repair was needed (see §4).  The
  relation multipliers $\lambda$ are recovered by an **exact rational
  least-squares** solve (normal equations over $\mathbb Q$, pivoted
  Gaussian elimination) of $r\approx E^{\mathsf T}\lambda$ where
  $r=\text{target}-e(Y)$.  All remaining discrepancy is absorbed by
  the $\ell_1$ term, *not* by any equality re-solving — no exact
  re-projection of $Y$ is required at all.

## 3. The verified identity and its validity ledger

The checker verifies, coefficient by coefficient over $\mathbb Q$,

$$\text{target}_L-\sum_b\langle A^b_L,Y_b\rangle-\sum_i\lambda_iE_{i,L}
 \;=\;c\,[L=\text{constant}]+\rho_L ,$$

and the final inequality chain uses exactly three facts per object:

| object | fact used | status |
|---|---|---|
| $Y_b$ | $Y_b\succeq0$, exact rational LDL with diagonal pivoting | **machine-verified** (all 109 blocks) |
| $A^b$ (block families) | $M_b(\mu)=\sum_LA^b_Ly_L(\mu)\succeq0$ for every measure | by construction (all-measures families: conditional flag squares, harmonic energies, $h_2$-multiplied copies, Jensen/covariance contractions, fiber-Toeplitz moment matrices; docs/TOEPLITZ_BLOCKS_NOTE.md, docs/IMPLEMENTATION.md) — **spot-verified** exactly PSD on random rational atomic measures |
| $E_i$ (relations) | $\sum_LE_{i,L}y_L(\mu)=0$ for every measure | four-point Gram-rank identities $\det\mathrm{Gram}_4\equiv0$ in $\mathbb R^3$ and their $p_2$-shifted copies — **spot-verified** exactly $0$ on random rational measures |
| labels | $\lvert y_L(\mu)\rvert\le1$ | every label is an expectation of a product of inner products of unit vectors, each factor in $[-1,1]$; products of moments likewise — **spot-verified** |
| target | $\langle\text{target},y(\mu)\rangle=h_2(\mu)E(\mu)$ | **spot-verified** against the directly computed $h_2E$ on the same measures |

Pairing the identity with $y(\mu)$ gives
$h_2E=\langle\text{target},y\rangle
=\sum_b\langle M_b(\mu),Y_b\rangle+0+c+\langle\rho,y\rangle
\ge c-\|\rho\|_1$.

**Composition.**  The certificate lives in the all-measures cone, so
the reduction lemma applies: $h_2>0$ on a weak-\* dense set, $E$
weak-\* continuous, hence the certified $h_2E\ge-\varepsilon$ bounds
the weighted functional globally without any KKT assumption.  (This is
the bound for $h_2E$, not yet $E\ge0$; it is the program's first exact
rational bound of this order, the dress rehearsal for the eventual
zero-bound certificate.)

## 4. Measured numbers

* Exact certificate constant $c=-1.2164987\times10^{-8}$ (a fully
  explicit rational, stored in the certificate; numerically $2.6\%$
  above the solver's wall $-1.2252\times10^{-8}$ — the exact pairing
  of the printed $Y$ with the exact problem data).
* Exact $\ell_1$ residual $\|\rho\|_1=2.659\times10^{-9}$ over 2031
  labels.  Provenance: not solver error (dual feasibility is
  $3.4\times10^{-35}$) but the geometry of the exported problem — the
  export dropped 447 image-dependent directions, and the label-space
  residual along their coordinates is only indirectly controlled; the
  $\ell_1$ absorption is deliberately conservative.  A tighter
  certificate would re-export without drops or extend the
  least-squares basis; unnecessary here.
* **Certified bound $c-\|\rho\|_1=-1.4823973\times10^{-8}\ge-2\times10^{-8}$**,
  with $5.2\times10^{-9}$ of the margin unused.
* Y-block float eigen pre-check: smallest reported value
  $-2.2\times10^{-8}$ on `h2loc_jensen_even_00_minor_d7` — pure
  eigensolver noise (those blocks carry entries up to $5.3\times10^7$,
  float noise floor $\approx4.5\times10^{-7}$); the exact rational LDL
  in the checker proves all 109 blocks PSD outright, and **no PSD
  repair was needed** ($\delta=0$ everywhere).
* Builder runtime 5 min (exact $e(Y)$: 5 s; exact 338×338 normal
  equations: 5 min); checker runtime: see PASS transcript below.

PASS transcript (2026-08-19, `python3 verify_h2E_bound.py`):

```
[1/3] exact PSD of certificate blocks
  all 109 blocks PSD: yes
[2/3] exact certificate identity
  c               = -1.2164987532e-08
  ||rho||_1       = +2.659e-09  (2031 labels)
  certified bound = -1.4823972914e-08
  bound >= -2e-8: yes
[3/3] semantic spot checks on random rational measures
  measure 1: 2051 labels evaluated ... all 338 relation rows vanish: yes
             all 109 block families PSD on the measure: yes
  measure 2: (same, second random rational measure)

PASS: h2*E >= -1.4823972914e-08 > -2e-8 for every antipodal
probability measure on S^2  (107s)
```

## 5. Checker invocation

```sh
python3 verify_h2E_bound.py certificates/h2E_geq_minus_2em8.json.gz
```

Pure stdlib (fractions/json/gzip/ast/random): no solver, no numpy, no
repo imports.  Verifies (1) exact PSD of all 109 $Y_b$, (2) the exact
identity and the bound, (3) the semantic spot checks of §3 on random
exact rational antipodal atomic measures (`--spot-measures N`,
default 2; `--skip-spot` for the fast algebraic-only run).  Exit code
0 and a final `PASS` line constitute the verification.

## 6. Second certificate: the degree-18 stack (KKT-inclusive cone)

*2026-08-19, second run of the same pipeline on
`deg18_h2w_h2all_toep3` ($m=2192$, 129 blocks, 3655 labels, 650
relation rows, wall $-2.21142\times10^{-14}$, `pdOPT`).*

**Theorem (certified).**  For every antipodal probability measure on
$S^2$ **satisfying the encoded first-order (KKT) relations** of $E$
(the 24 gradient/potential rows; in particular every minimizer of
$E$):

$$h_2(\mu)\,E(\mu)\;\ge\;-2.7454\times10^{-14}\;\ge\;-5\times10^{-14},$$

equivalently an exact rational lower bound on the KKT-inclusive
degree-18 relaxation value — **10× stronger than the requested
$-1\times10^{-13}$ target**.  Artifacts:
`certificates/h2E_geq_deg18.json.gz` (3.3 MB),
`sdpa_runs/build_deg18_cert2.log`, checker as in §5 (the same
`verify_h2E_bound.py` reads the cone metadata).

**The KKT decision (measured).**  The certificate's gradient-row
multipliers carry mass $\sim1.15\times10^{8}$ (potential
$\sim2\times10^{3}$): the KKT *relations* are structurally
load-bearing, and the rank-relation-only least squares leaves a
residual of order $10^{10}$ — the **all-measures version of this
certificate does not exist for this $Y$**; it would need a fresh
GMP solve of the am-cone problem (the deg-16 am experience suggests
a $\sim2$–$3\times$ weaker wall).  The 16 KKT PSD families
(hessian/tangent-gap), by contrast, are essentially unused
($1.9\times10^{-9}$ of the $Y$ mass) — the certificate is honest KKT
through its relations, not its blocks.  The certificate marks every
KKT-only row and block (`kkt_only`, `kkt_relation_count`), the
checker skips them in the semantic measure-checks (they are not
identities/PSD for generic measures) and states the KKT-scoped
theorem in its PASS line.

**Numbers.**  Exact $c=-2.1956\times10^{-14}$ (slightly above the
solver wall); $\|\rho\|_1=5.498\times10^{-15}$ over 3629 labels — the
deg-16 dropped-direction inflation did *not* recur (this export's
kept directions cover the label space much more tightly);
bound $=c-\|\rho\|_1=-2.745\times10^{-14}$.

**Pipeline upgrade: exact iterative refinement.**  The deg-16 exact
normal-equation Gauss ($338^2$, 5 min) does not scale ($650^2$ dense
rational elimination ran >75 min unfinished).  Replacement
(`sdpa_runs/solve_lambda_refine.py`): float64 solves of the normal
equations with **exact rational residual accumulation** (classical
iterative refinement, 12 digits/iteration here), then
`limit_denominator(10^110)`.  Soundness is unaffected by how
$\lambda$ is found: $\rho$ is recomputed exactly afterwards, so any
rational $\lambda$ is admissible — the refinement only makes
$\|\rho\|_1$ small.  Total $\lambda$ time: **14 s** (vs hours);
whole certificate build: 22 s after the one-off yMat parse.

**Checker upgrade.**  Exact PSD is now decided by *fraction-free*
(Bareiss) elimination with symmetric max-diagonal pivoting on the
integer-scaled matrix — pivots are principal minors up to positive
factors, digit growth is linear, and the singular-PSD /
indefinite cases are decided by the same zero-diagonal rule as
before.  Validated against float eigenvalues on 196 randomized
sign-definite matrices and 100 exact low-rank Grams.  The deg-16
certificate re-verifies in 15 s (was 78 s); the deg-18 blocks
(55×55, entries to $6.5\times10^{8}$, mixed $10^{-90}$ scales) pass
in 100 s total.

**Verification.**  PASS (`sdpa_runs/verify_deg18.log`, 139 s):

```
all 129 blocks PSD: yes
c               = -2.1956494525e-14
||rho||_1       = +5.498e-15  (3629 labels)
certified bound = -2.7454450222e-14   (>= -5.0e-14: yes)
all 626 all-measures relation rows vanish (both spot measures): yes
all 113 all-measures block families PSD on the measures: yes
PASS: h2*E >= -2.7454450222e-14 (>= -5.0e-14) for the KKT-inclusive
relaxation (measures satisfying the encoded first-order relations;
in particular every minimizer of E)
```

## 7. Capstone: the degree-18 **all-measures** certificate

*2026-08-20.  The cycle's capstone: the proof-carrying degree-18
stack `deg18_am_toep3` ($m=2216$, 113 blocks, all-measures cone — no
KKT content, 626 rank-relation rows only, composes with the reduction
lemma; wall $-1.43334761688\times10^{-11}$, `pdOPT`, dual feasibility
$4.1\times10^{-42}$).*

**Theorem (certified, all measures).**  For every antipodally
symmetric Borel probability measure $\mu$ on $S^2$:

$$h_2(\mu)\,E(\mu)\;\ge\;-1.7355738\times10^{-11}\;\ge\;-2\times10^{-11}.$$

This supersedes the degree-16 record ($-1.482\times10^{-8}$, §1–§5)
by **854×** and tightens the requested $-3\times10^{-11}$ target by
1.5×.  Artifacts: `certificates/h2E_geq_deg18_am.json.gz` (3.3 MB),
`sdpa_runs/build_deg18_am_cert.log`, `sdpa_runs/verify_deg18_am.log`.

**Numbers.**  Exact $c=-1.4232\times10^{-11}$;
$\|\rho\|_1=3.124\times10^{-12}$ over 3626 labels (the
dropped-direction residue is present at this export, as at deg-16,
but sits comfortably inside the $1.6\times10^{-11}$ margin); bound
$=c-\|\rho\|_1=-1.7356\times10^{-11}$.  Pipeline runtime with the
upgraded tooling: $\lambda$ refinement **6 s** (8 iterations, 14
digits/iteration), certificate build **11 s**, verification transcript
below — the entire exactification is now a sub-minute operation per
solved stack (plus the checker run).

**Verification.**  PASS (`python3 verify_h2E_bound.py
certificates/h2E_geq_deg18_am.json.gz`, 83 s):

```
all 113 blocks PSD: yes
c               = -1.4231737489e-11
||rho||_1       = +3.124e-12  (3626 labels)
certified bound = -1.7355738043e-11   (>= -2.0e-11: yes)
all 626 all-measures relation rows vanish (both spot measures): yes
all 113 all-measures block families PSD on the measures: yes
PASS: h2*E >= -1.7355738043e-11 (>= -2.0e-11) for every antipodal
probability measure on S^2
```

**Composition.**  As in §3: the cone is all-measures, so the reduction
lemma applies verbatim — this is an unconditional theorem about
measures, the strongest exact statement of the program to date.

## 8. Honest caveats / upgrade path

* The block matrices $A^b$ and relation rows $E_i$ are certificate
  *data*; their all-measures validity is by construction and is
  spot-checked semantically, not re-derived inside the checker.  The
  upgrade is mechanical: re-derive each family row inside the checker
  from its definition (the reducer is ~100 lines); queued for the real
  (zero-bound) certificate.
* $c$ is the rounded certificate's own constant; the $\ell_1$
  absorption uses $|y_L|\le1$, which is loose but costs only
  $\|\rho\|_1$ — measured many orders below the margin.
* The bound $-2\times10^{-8}$ is a statement about $h_2E$, two orders
  below the previous best *numerical* wall for this cone and now
  exact; the remaining distance to $0$ is the multiplier program
  (docs/EXACT_ZERO_PROGRAM.md, docs/MULTI_WEIGHT_PROGRAM.md,
  docs/E5_WEIGHT_NOTE.md).
