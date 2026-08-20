# Limit-certificate extraction from the degree-14/16/18 tower

*Agent L, 2026-08-20.  Data: the three GMP-solved duals of the
structured Jensen/fiber-Toeplitz family tower (KKT-inclusive cone),
`sdpa_runs/deg{14,16,18}_h2w_h2all_toep3.result`, plus the two
all-measures duals `deg{16,18}_am_toep3.result`.  Method: per-label
fingerprints (`e_L(Y)=\sum_b\langle A_L^b,Y_b\rangle`), per-block Gram
extraction at 40 digits (`sdpa_extract.py` conventions), corner-nested
cross-degree comparison, model fits, and closed-form identification.
Analysis scripts and per-degree extractions (`fp_deg*.json`,
`spike_antisym.json`, `converging_labels.json`) live in the session
scratchpad; everything below is reproducible from the committed
exports + result files.*

**Headline.**  The naive blockwise limit of the dual sequence does
**not** exist — dual mass keeps escaping to the truncation boundary —
but the sequence is far from structureless.  Three exact objects were
identified: (i) the KKT tower's dominant component is a rank-one
square of the *kernel-potential difference* $[K(x\cdot y)-K(z\cdot
y)]^2$ — proved to be an **exact null of the KKT-inclusive cone**
(so its huge multiplier is gauge, and the KKT-vs-am wall gap is the
feasible space's exact knowledge of the equal-potential identity);
(ii) the pair-sector
complement block converges to the square of the admissible deviatoric
leaf $(t^2-\tfrac13)$, with a transient exactly along the kernel top
$48t^4-32t^6$; (iii) the KKT walls decay *faster than any fixed*
$q^{d^2}$ — both prescribed two-parameter models overshoot zero —
while the **all-measures walls (three points, the third solved by
this session) decay only geometrically ($\sim900\times$/2 degrees)
with $c_\infty<0$ fits**: the proof-carrying tower needs a new
mechanism, not more degree.  The escaping tail is a
swap-antisymmetric boundary layer in the two-root *even* sectors
whose $X\!\cdot\!Z$-modulation index climbs slowly with degree
($k\approx1\to2\to3$).

---

## 1. Alignment across degrees (L1)

* **Labels nest exactly**: the fingerprint label sets satisfy
  $L_{14}\subset L_{16}\subset L_{18}$ ($1287\subset1944\subset3521$;
  verified as string sets from the map sidecars).
* **Base blocks correspond by name** (74 shared names); moment-type
  blocks grow (`two_root_even_00`: $30\to45\to55$; sizes follow the
  all-even/all-odd exponent count of `monomials(d/2)`), 14 harmonic
  $1\times1$ blocks stay fixed.
* **Corner nesting is exact.**  For every shared base block and for
  families 0–28 of the tower, the degree-16 label matrices restricted
  to the top-left degree-14 corner equal the degree-14 matrices
  entrywise (max deviation $0.0$ over all shared labels), and labels
  new at 16/18 have zero corner.  `sos_search.monomials` is
  degree-graded, so Gram corners are directly comparable across
  degrees.
* **Caveats.**  (a) The degree-14 v3 export has **39** families, not
  43 — `ftoep3_even_11_*`, `ftoep2_*_r4`-type entries only enter at
  16/18; family *roles* align by construction order for indices 0–28
  (the Jensen/h2comp towers, name-shifted `_d7/_d8/_d9`), while the
  fiber-Toeplitz block indices 29+ *migrate in modulation index*
  (`ftoep2_even_00_r3 -> _r4 -> _r5`, etc.) and are **not**
  corner-nested (different objects per degree).  (b) `phase.value`:
  deg-14/18 pdOPT, deg-16 pdFEAS (gap $\sim10^{-11}$) — deg-16 entries
  carry more solver noise.

## 2. The wall laws (L2)

40-digit walls (recomputed from `xVec`, shift $+2/3$):

| tower | $d=14$ | $d=16$ | $d=18$ |
|---|---|---|---|
| KKT-inclusive | $-3.6605120927\times10^{-6}$ | $-5.2504467051\times10^{-9}$ | $-2.2114191705\times10^{-14}$ |
| all-measures | $-1.2080927467\times10^{-5}$ | $-1.2251946394\times10^{-8}$ | $-1.4333476169\times10^{-11}$ |

(The $d=14$ all-measures point was produced by this session:
`deg14_am_toep3.dat-s`, $m=810$, 200-bit, pdFEAS — see §7.)

KKT ratios per 2 degrees: $697.2\times$, then $237{,}424\times$;
$\ln|W|$ increments $-6.5470$, $-12.3776$.

* **Pure geometric $r^d$: rejected** (local $r$ = 0.038 → 0.0021).
* **Theta law $q^{d^2}$ with fixed $q$: rejected.**  Local
  $q(14\!\to\!16)=0.89663$, $q(16\!\to\!18)=0.83358$ — the decay is
  *super*-$q^{d^2}$.  (Numerology: the second value is $5/6$ to
  $3\times10^{-4}$; recorded, not relied on.)
* **Prescribed model $c_\infty+Aq^{d^2}$ (ENRICHMENTS §"extraction
  fallback"), exact 3-point solve**: $q=0.896634$,
  $A=-7095.1$, $c_\infty=+3.13\times10^{-12}$.  The same for
  $c_\infty+Ar^d$: $r=0.037900$, $c_\infty=+7.53\times10^{-12}$.
  Both give $c_\infty>0$, i.e. the models *overshoot through zero*:
  the true walls fall faster than either family.  **This is the
  best evidence yet that the walls converge to exactly $0$ from
  below** (a genuinely negative limit would force $c_\infty<0$; a
  wall sequence *slower* than the model would force $c_\infty<0$ as
  well).
* **Single-power law** $|W|=\exp(-cd^{\,p})$: the three points are
  consistent with exactly one exponent, $p^\* = 6.116$
  ($c=5.07\times10^{-7}$).  The proximity to the kernel degree 6 is
  noted; with 3 points this is a fit, not a finding.
* Quadratic-exponent fit (exact): $\ln|W| = -129.94 + 18.591d
  -0.72882 d^2$; the enormous linear prefactor says this
  parameterization is not natural either.
* **All-measures tower (3 points, this session)**: ratios $986.0\times$
  then $854.8\times$ per 2 degrees — **nearly constant geometric,
  mildly decelerating**; local theta-$q$ *rises* ($0.89146\to0.90549$).
  The exact 3-point fits give
  $c_\infty+Aq^{d^2}$: $q=0.89145$, $c_\infty=-9.39\times10^{-12}$
  ($=0.65\,W_{18}$), and $c_\infty+Ar^d$: $r=0.03184$,
  $c_\infty=-1.91\times10^{-12}$ ($=0.13\,W_{18}$).  **Both give
  $c_\infty<0$**: on its current enrichment set the proof-carrying
  tower extrapolates to a *stall* around $-10^{-12}$, not to zero.
  The single-power consistency exponent is $p^\*\approx0.7$
  (sub-geometric) — a plain approximation ladder.  The KKT/am gap
  widens ($3.3\times$ at 14, $2.3\times$ at 16, $648\times$ at 18) —
  the KKT acceleration is a KKT-face effect (see §4a).

**Verdict on the law**: the two towers obey qualitatively different
laws.  KKT: super-geometric and super-theta, acceleration itself
accelerating — consistent with walls $\to0^-$ (both 2-parameter
models overshoot, $c_\infty>0$); a drifting-$q$ theta law
$q_d^{d^2}$ (pole sliding with degree, cf. the pole-law entries in
PLAN.md) or $\exp(-cd^{\,6})$ both remain viable.  All-measures:
ordinary geometric with mild deceleration and $c_\infty<0$ fits —
the family tower v3 does **not** yet contain the mechanism that
drives the KKT walls; without new all-measures-valid enrichments the
am ladder is headed for a plateau near $-10^{-12}$.

## 3. Gauge anatomy: what the raw fingerprints actually measure

Classification of the 1287 shared labels by $(e_L(Y_{14}),
e_L(Y_{16}), e_L(Y_{18}))$ with Aitken extrapolation
(`converging_labels.json`): 44 "converging", 24 decaying, 639
numerically zero, 580 unclassifiable at 3 points.  **But the raw
fingerprint is gauge-dominated** and the headline magnitudes
($|e_L|\sim10^{11}$ at $d=14$) are *relation multipliers*, not
certificate content:

* label pairs tied by ideal relations appear with equal-and-opposite
  huge values (e.g. `('product',p2,p4,p6)` $=+1.5407\times10^{11}$ vs
  `('product',p2,('triangle',0,4,6))` $=-1.5388\times10^{11}$;
  class sums are $\sim100\times$ smaller);
* the dominant Gram spike (§4a), expanded to label space
  ($\|e\|=3.4\times10^{11}$), pairs with **every** feasible direction
  of the export at $\le9.6\times10^{2}$ — a relative cancellation of
  $3\times10^{8}$.  It is, to that accuracy, an element of the
  relation-ideal cone: invisible to all feasible moment vectors of
  the KKT-inclusive problem;
* a family of *exact null corners* carries arbitrary shared mass
  $s_d$ ($1.346\times10^{8},\,2.907\times10^{8},\,1.388\times10^{6}$ —
  non-monotone, hence non-physical): the $(0,0)$ Jensen corners
  ($T_{00}-G_{00} = 1-1=0$ on probability measures — the
  `jensen_pair_d*` top eigenvector is exactly $(1,0,0,0,\dots)$),
  the `flag_1` direction $(5/24,-1,1)/\|\cdot\|$, and at $d=18$ also
  `ftoep2_even_00_r5`.  The interior point spreads one free mass over
  this orbit; **any cross-degree extrapolation must quotient it
  out.**

Consequence: *per-label extrapolation in the raw label basis is the
wrong gauge*; the meaningful objects are per-block Gram data modulo
these null directions, as below.

## 4. Identified structure (L3)

### 4a. The KKT tower's dominant object: $\big[K(x\cdot y)-K(z\cdot y)\big]^2$

The largest Gram block at every degree is `h2loc_two_root_even_00`
(trace $1.25\times10^{11}\to1.06\times10^{10}\to2.10\times10^{9}$),
and it is a *rank-one spike in the swap-antisymmetric sector*
($X\leftrightarrow Z$).  Its eigenvector, in the two-root monomial
basis $(X\!\cdot\!Y)^i(Z\!\cdot\!Y)^j(X\!\cdot\!Z)^k$, is supported on
$\{(2,0,0),(4,0,0),(6,0,0)\}$ minus the swap image, with coefficients
**exactly** $(5,-12,8)/\sqrt{466}$:

$$v \;\propto\; p(X\!\cdot\!Y)-p(Z\!\cdot\!Y),\qquad
p(t)=8t^6-12t^4+5t^2=\tfrac14\big(K(t)+\tfrac43\big).$$

Since the constant cancels in the difference, the spike **is** the
$h_2$-localized square of the kernel-potential difference at two
roots, $\big[K(X\!\cdot\!Y)-K(Z\!\cdot\!Y)\big]^2$.  Residuals of the
computed eigenvector from this exact direction: $3.6\times10^{-8}$
($d{=}14$), $0.108$ ($d{=}16$), $0.0112$ ($d{=}18$); the $d=16/18$
deviation is itself structured — the same polynomial modulated by
$(X\!\cdot\!Z)^2$, with amplitude falling $\sim10\times$ per step —
so $v_d\to v_{\text{exact}}$.

Multipliers: $\lambda_d = 1.2467\times10^{11},\ 1.0377\times10^{10},\
2.1003\times10^{9}$.  **Exact-null verdict**: pairing the *exact*
unit-normalized $[\Delta K]^2$ square against the deg-16 KKT export
gives $\max_j|\langle\text{spike},q_j\rangle| = 3.3\times10^{-20}$
and $\langle\text{spike},y_0\rangle = 0.0$ — it is a **certified
exact null of the KKT-inclusive cone** (the KKT rows encode constant
potential on the support, so the equal-potential defect square has
identically zero expectation on the feasible affine span).  Hence
$\lambda_d$ itself is pure gauge (like $s_d$), and what is
load-bearing is the *modulated copy*: the same polynomial times
$(X\!\cdot\!Z)^2$ pairs at $0.82$ with the directions — a genuine
constraint — and its amplitude in the dual, $\varepsilon_d = 0.108\to
0.0112$, is real certificate data decaying $\sim10\times$ per step.
Against the *all-measures* export the pure $[\Delta K]^2$ square pairs
at $4.3\times10^{-3}$: **not** null there.  The am duals accordingly
do not carry the spike (top eigenvalue $\sim3\text{–}5\times10^4$,
overlap $\approx0.2$).  So the entire KKT-vs-am wall gap ($648\times$
at $d=18$) is the *feasible space* knowing the equal-potential
identity exactly — information the proof-carrying cone lacks.

### 4b. The pair-sector complement block: the deviatoric leaf, exactly

`h2comp_gram_pair_d*` (the $(1-h_2)$-localized pair Gram complement,
basis $(X\!\cdot\!Y)^{0,2,4,6,\dots}$) is **rank one at every degree
and in both towers** ($2\times2$ determinant $\le10^{-10}\lambda^2$),
with converging direction and a closed-form transient:

| | $w_1/w_2$ | $w_3/w_2$ | $w_4/w_3$ | $\lambda$ |
|---|---|---|---|---|
| KKT 14 | $-0.3210884$ | $-0.13033$ | $-0.6666730$ | $0.8210$ |
| KKT 16 | $-0.3327976$ | $-0.0056968$ | $-0.6666589$ | $0.6165$ |
| KKT 18 | $-0.3333304$ | $-2.903\times10^{-5}$ | $-0.6669$ | $0.2968$ |
| am 14 | $-0.3243385$ | $-0.11621$ | $-0.6665953$ | $0.8652$ |
| am 16 | $-0.3328001$ | $-0.0062106$ | $-0.6666566$ | $0.5530$ |
| am 18 | $-0.3333155$ | $-2.116\times10^{-4}$ | $-0.6666677$ | $0.4877$ |

* $w_1/w_2 \to -1/3$ (Aitken $-0.3333558$): the limit vector is
  $(X\!\cdot\!Y)^2-\tfrac13 = \tfrac12(T_2+\tfrac13)$ — **the
  admissible spin-0 deviatoric leaf** of the (E1) classification
  (docs/SHARP_STRUCTURE.md II.1).
* $w_4/w_3 = -2/3$ at *every* degree and in *both* towers: the
  transient is exactly proportional to $48t^4-32t^6$ — **the
  top-degree part of the kernel** — with amplitude dying
  super-geometrically ($\times22.9$, then $\times196$ per step,
  mirroring the wall acceleration).
* $\lambda$: KKT copy decays with the wall; the am copy converges:
  $0.8652\to0.5529\to0.4877$, Aitken $\lambda_\infty=0.4704$
  ($\rho=0.209$) — **the deviatoric-leaf square survives in the am
  limit with a finite positive coefficient** $\approx0.47$ (no
  credible closed form at this precision; $8/17=0.4706$ is within the
  extrapolation error and is recorded as numerology only).

### 4c. A shared certificate core, and what is *not* shared

At $d=18$, the top-5 eigenspaces of `h2loc_two_root_even_00` in the
KKT dual (swap-symmetric part) and the am dual share a 3-dimensional
subspace with principal cosines $(0.967,\,0.944,\,0.613)$; the
remaining two directions are tower-specific.  The core vectors are
dense in monomials of shells 5–8 (leading components $(2,2,4)$,
$(3,3,3)$, $(3,1,3){+}(1,3,3)$, …) — no sparse closed form found at
this precision.  This core is the best-supported "limit object"
beyond §4a/4b, but it lives at the moving truncation boundary
(§5), so it should be read as the *shape of the escape*, not as a
convergent finite block.

### 4d. PSLQ / closed-form summary

Identified exactly (rational, PSLQ-grade agreement):
$(5,-12,8)/\sqrt{466}$ spike $= [\Delta K]^2$, an exact null of the
KKT cone (§4a; eigenvector agreement $3.6\times10^{-8}$ at $d{=}14$,
pairing zero at $10^{-20}$); $-1/3$ and $-2/3$ ratios (§4b,
$10^{-5}$–$10^{-6}$ after extrapolation); Jensen null corner
$(1,0,0,0)$ and `flag_1` null $(5/24,-1,1)$ (§3).  Suggestive only: local theta-$q=0.83358\approx
5/6$ at $16\!\to\!18$.  PSLQ against $\{\pi^2, (4/5)^5, 2/3$-powers$\}$
on the wall values and $\lambda$'s produced no credible hits (as
expected for solver-path-dependent magnitudes).

## 5. The escaping tail (L4)

Diagonal Gram mass by total monomial degree ("shell") in
`h2loc_two_root_even_00`, am tower:

    d=16:  0:1.7  2:5.1e2  4:6.9e3  5:9.9e3  6:1.8e4  7:2.6e4  8:1.6e4
    d=18:  0:2.0  2:3.0e2  4:5.2e3  5:1.0e4  6:1.5e4  7:4.7e4  8:4.1e4  9:1.3e4

The profile is a wave riding the truncation boundary: the peak sits
at shell $\approx d/2-1$ and its amplitude *grows* ($\times1.8$ at
shell 7, $\times2.5$ at shell 8 per 2 degrees; KKT symmetric-sector
top eigenvalue grows $1.02\times10^4\to2.51\times10^4\to
1.28\times10^5$, accelerating).  Characterization of the boundary
layer:

* **sector**: the *even* two-root sectors only (`even_00` and
  `even_11` traces grow; both `odd` sectors decay,
  $7.8\times10^4\to3.2\times10^4$);
* **symmetry**: top escape directions are swap-**antisymmetric**
  (e.g. $d{=}18$: $\pm0.39$ on $(2,4,2)-(4,2,2)$ and
  $(5,1,1)-(1,5,1)$);
* **modulation**: mass by $X\!\cdot\!Z$-exponent $k$ *climbs with
  degree*: the am peak sits at $k=0\text{–}1$ at $d{=}14$, $k=2$ at
  $d{=}16$, $k=2\text{–}3$ at $d{=}18$ ($k{=}3$: $1.4\times10^4\to
  3.2\times10^4$) — a slowly rising fiber-modulation index, *not*
  captured by the explicit `ftoep*` families, which stay empty in
  both towers (traces $\le10^{-9}$ up to gauge mass);
* the `h2comp_gram_even_*` and Jensen towers at $d=18$ hold only
  $O(10^3$–$10^4)$ am mass — the v3 families did their job
  (they moved the law), but the residual escape has moved past them.

**Verdict (L4).**  As a *finite object in the current algebra* the
$d\to\infty$ certificate does **not** exist: a boundary layer with
slowly growing amplitude escapes along the leaf-degree axis (and
secondarily the $X\!\cdot\!Z$-modulation axis) inside the even
two-root sectors.  The convergent part is small and now identified
(§4a, §4b + a 3-dim core §4c); everything else is transient or gauge.
The limit certificate exists only in an *enriched* algebra that
resums the boundary layer: finitely many converging blocks plus one
tail family indexed along the degree axis.  Confidence: high on the
gauge anatomy, the $[\Delta K]^2$ and deviatoric-leaf identifications,
and the KKT model-overshoot ($c_\infty>0$) statements; medium-high on
the am tower being a plain geometric ladder with a predicted stall
(3 points, consistent fits); medium on the tail-profile growth rates;
low on any specific KKT wall law.

## 6. Recommendation for the closing construction

1. **Port the equal-potential identity to the proof-carrying cone.**
   The 648× KKT advantage is exactly the feasible span's exact
   knowledge that $h_2\text{loc}[K(X\!\cdot\!Y)-K(Z\!\cdot\!Y)]^2$ has
   zero expectation (§4a) — the am cone lacks it and its escape
   directions ($(2,4,2)-(4,2,2)$, $(5,1,1)-(1,5,1)$, … §5) are visibly
   trying to assemble $\Delta K$-cross terms from boundary monomials.
   The all-measures-valid bridge is the *kernel-row covariance*:
   $\mathrm{Var}_\mu$-type blocks $T-G\succeq0$ conjugated by the
   kernel coefficient vector $(5,-12,8)$ on $(t^2,t^4,t^6)$ — i.e.
   dedicated small blocks $\mathrm{Cov}(K\text{-row})$, plus
   $(X\!\cdot\!Z)^{2m}$-modulated copies ($m=1,2$; the $m=1$ copy is
   load-bearing even in the KKT cone).  These are sections of the
   existing Jensen/covariance families, so they are exactly rational
   and cheap; pair-test against the current am escape first.
2. **Resum the even-sector boundary layer along the degree axis.**
   The tail is a $q$-profile in leaf degree at fixed low
   $X\!\cdot\!Z$-modulation ($k\le3$), swap-antisymmetric.  The
   generating-function labels $k_G$ of the one-sided truncation
   principle, built from the two admissible spin-0 kernels
   ($T_2+\tfrac13$, $T_6+\tfrac13$) and their products with
   $(X\!\cdot\!Z)^{2}$, are the natural finite dictionary: they carry
   exactly the $(u^m-w^m)$-difference structure the escape uses.
3. **Weight the deviatoric square.**  §4b says the am certificate's
   stable pair-sector content is $\lambda\,[(X\!\cdot\!Y)^2-\tfrac13]^2$
   under $(1-h_2)$; giving this leaf its own weighted copy (multiplier
   dictionary entry $w=(t^2-\tfrac13)^2$) aligns the certificate with
   its own limit and is cheap to pair-test against the current escape
   data before any solve.
4. **Spend future GMP budget on the am tower**, which is the one that
   composes with the reduction lemma; the third am point (deg-14+v3
   am) delivered by this session (§7) shows the current am ladder is
   geometric with $c_\infty<0$ fits — i.e. degree-20 am without new
   blocks is predicted to land near $-1.7\times10^{-14}$ and then
   flatten; measure item 1's blocks at degree 14/16 *first* (cheap
   solves) before buying more degree.

## 7. Session artifacts / in-flight

* Extraction pipeline: `extract_fp.py`, `analyze_fp.py`,
  `gram_conv.py` + per-degree `fp_deg*.json` (fingerprints, 40-digit
  Gram blocks, traces/eigs), `spike_antisym.json`,
  `converging_labels.json` — session scratchpad.
* Bug fixed in `toeplitz_export.py`: `--out` with a
  (degree, version, cone) triple absent from the default table no
  longer KeyErrors.
* **Third all-measures point (delivered)**: `deg14_am_toep3.dat-s`
  ($m=810$, 97 blocks, 39 families) exported and solved at 200-bit
  in 19.4 min, pdFEAS, objValPrimal $=-6.6667874759413337\times
  10^{-1}$, wall $=-1.2080927467\times10^{-5}$ (40-digit
  recomputation from xVec).  Capture pickle
  `toeplitz_capture_deg14_h2w_v3_am.pkl` and `fp_deg14am.json` in the
  scratchpad; 3-point am analysis in `am3_analysis.py`.  Structure of
  the $d{=}14$ am dual: no $[\Delta K]^2$ spike (overlap $0.0000$),
  spread spectrum ($\lambda_1=8.0\times10^3$,
  $\lambda_2=5.9\times10^3$), boundary-shell profile peaking at
  shells 6–7, $X\!\cdot\!Z$ modulation peak at $k=0\text{–}1$ —
  fully consistent with §5's traveling-wave picture.

**§7 addendum (orchestrator, 2026-08-20): the third all-measures
point.**  `deg14_am_toep3` (m = 810, 200-bit, pdFEAS):
objValPrimal = −0.66667874759413337 → wall = **−1.20809e−5**.
The all-measures tower is therefore

| $d$ | 14 | 16 | 18 |
|---|---:|---:|---:|
| wall | −1.2081e−5 | −1.2252e−8 | −1.4333e−11 |
| ratio | | 986× | 855× |

— **essentially pure geometric decay** ($r\approx0.033$ per two
degrees), in contrast to the KKT tower's acceleration
(697× → 237,424×).  Conclusion sharpened: the acceleration is a
KKT-face effect (§4a's equal-potential defect); the proof-carrying
cone decays geometrically toward 0 and will not land at any finite
degree by escalation alone.  The closing construction (§6 items 1–2:
adjoin the defect square's all-measures analogue and resum the
even-sector boundary layer via the one-sided-truncation
generating-function labels) is *required* for the am cone, not
optional.
