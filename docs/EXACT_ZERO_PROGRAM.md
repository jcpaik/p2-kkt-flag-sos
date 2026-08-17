# Reaching the exact zero bound: enriched certificate algebras

The assumption-free four-point hierarchy converges to zero
super-exponentially in the degree
([Numerical results](../RESULTS.md)): bounds decay like
$\exp(-c\,d^2)$, while the $\varepsilon$-canonical certificate trace
diverges like $\operatorname{tr}C(\varepsilon)\approx1.07/\varepsilon$.
Both are the signature of **non-attainment**: the target functional $E$
lies on the boundary of the closure of the finite-degree certificate
cone, but not in the cone itself, so raising the degree can approach
zero but never land on it.

This document works out three ways to enlarge the certificate algebra so
that a *finite* object can be exact — rational certificates
(denominators), square-root adjunction, and entire/hypergeometric kernel
generators — and shows that they are not equally promising.  The first
is implemented (`--h2-weighted-target`).

Throughout, labels and the expansion map $e(\cdot)$ are as in
[Implementation](IMPLEMENTATION.md);
$p_j=\iint(x\cdot y)^j\,d\mu\,d\mu$, and

$$h_2=\frac{3p_2-1}{2}=\sum_{m=-2}^{2}\big|\hat\mu_{2m}\big|^2\ \ge\ 0$$

is the spin-2 harmonic energy, zero exactly on second-moment-isotropic
measures.

## 1. What the pole is telling us

Write $C_\varepsilon$ for the minimal-trace certificate of
$E+\varepsilon\ge0$ (measured in the legacy $(3/16)E$ scale).  The measurements

* $\operatorname{tr}C(\varepsilon)\cdot\varepsilon\to1.07\pm0.04$
  (a **simple** pole),
* heavy Gram entries growing $10^6\to10^{11}$ when the inactive families
  are removed,
* the $L^2$ blow-up criterion of
  [Limit certificate](LIMIT_CERTIFICATE.md) §5

say that $\varepsilon C_\varepsilon\to Y_0$ with $Y_0\succeq0$,
$Y_0\neq0$, and $e(Y_0)=0$: a positive-semidefinite combination whose
label expansion is the zero functional.  $Y_0$ is a recession direction
of the certificate cone; the sharp certificate escapes to infinity along
it.  This is precisely the classical mechanism (Motzkin-type) by which a
nonnegative polynomial fails to be a sum of squares while $q\cdot f$ is
one for a suitable positive multiplier $q$ — and a *simple* pole
predicts a multiplier vanishing to *first* order along the leak.

The leak variable is known from the moment diagnostics: the optimal
pseudo-moment pays $h_2$ linearly and earns contraction violations of
order $\sqrt{h_2}$.  The multiplier that vanishes to first order along
the leak is $h_2$ itself.

## 2. Route A — rational certificates: the weighted target $h_2E$

### 2.1 Reduction lemma (no stratification is needed)

**Lemma.** If $h_2(\mu)\,E(\mu)\ge0$ for every antipodal probability
measure $\mu$, then $E(\mu)\ge0$ for every antipodal probability
measure.

*Proof.* On $\{h_2>0\}$ divide.  If $h_2(\mu)=0$, let
$\mu_t=(1-t)\mu+t\,\delta_{\pm e}$; its second moment
$\frac{1-t}{3}I+t\,ee^{\mathsf T}$ is anisotropic for every $t\in(0,1]$,
so $E(\mu_t)\ge0$, and $E$ is weak-\* continuous (polynomial kernel), so
$E(\mu)=\lim_{t\downarrow0}E(\mu_t)\ge0$. $\qquad\blacksquare$

So a certificate for the weighted target proves the full conjecture —
the isotropic stratum comes for free by density, and no separate
isotropic-branch argument is needed.

### 2.2 The weighted target is polynomial in the existing algebra

$$h_2E=\Big(\tfrac{3p_2-1}{2}\Big)
\Big(-\tfrac43+20p_2-48p_4+32p_6\Big)
=\tfrac23-12p_2+24p_4-16p_6
+30\,p_2p_2-72\,p_2p_4+48\,p_2p_6 .$$

The products $p_2p_j$ linearize as genuine four-sample moments — the
disconnected two-edge labels
$(\texttt{product},(\texttt{pair},2),(\texttt{pair},j))$ already present
in the four-point label algebra.  No new machinery is required: the
weighted problem is the same relaxation with a different rational
objective vector, exported by

```sh
python3 sos_search.py \
  --export-sdpa PROBLEM.dat-s --h2-weighted-target \
  --degree DEGREE --no-pointwise-sos \
  --harmonics --three-point-flags --four-point-flags --two-root-flags \
  --gradient --potential --hessian \
  --global-tangent-gaps --rank-relations
```

The natural quadratic module for a weighted certificate contains
$h_2\times(\text{flag square})$ terms; `--h2-localized-flags` provides
exactly these and should be ablated in.

### 2.3 What the outcomes would mean

| outcome at degree 14–18 | interpretation | next step |
|---|---|---|
| bound $=0$ to solver precision ($<10^{-20}$) | exact rational certificate $E=e(C)/h_2$ exists at finite degree | attainment check, then rationalization (§2.4) |
| bound negative, several orders better than $-4.5\times10^{-3}$ | right denominator family, wrong power or missing module term | try $h_2^2E$; add `--h2-localized-flags`; re-measure the $\varepsilon$-pole on the weighted target |
| bound in the same band | the leak is not $h_2$-graded | measure $Y_0$ with the strictly convex (Frobenius) selector and read the denominator off its label expansion |

### 2.4 The attainment diagnostic

The decisive follow-up measurement is the $\varepsilon$-trace law on the
*weighted* target: solve the minimal-trace selector of
$h_2E+\varepsilon\ge0$ for a decreasing $\varepsilon$-sequence.

* $\operatorname{tr}C_w(\varepsilon)=O(1)$ as $\varepsilon\downarrow0$:
  attainment is restored; the max-margin go/no-go of
  [PLAN](../PLAN.md) §4 applies verbatim to the weighted problem, and
  the exactification pipeline (rational rounding on the interior point,
  `verify_certificate.py`) can run unchanged.
* $\operatorname{tr}C_w(\varepsilon)\sim c/\varepsilon$ again: the pole
  survives one $h_2$; iterate the residue analysis — the new $Y_0^{(w)}$
  names the next denominator factor.

### 2.5 Choosing denominators beyond $h_2$

If $h_2$ alone is not the whole denominator, the candidates, in order:

1. $h_2^2$ — if the pole on the weighted target has order one again;
2. $h_2+\beta h_4+\gamma h_6$ with $h_\ell$ the higher harmonic
   energies (all are nonnegative labels) — if the residue $Y_0$ loads
   on spin-4/spin-6 blocks;
3. the "squared-form" weight $\iint T_3(x\cdot y)^2-\tfrac13+h_2$ —
   nonnegative by the frame-form identity of
   [Structure](STRUCTURE.md) §1 only *jointly*, so use with care.

In every case the reduction lemma survives as long as the weight is
nonnegative on all measures and strictly positive on a dense class.

## 3. Route B — square roots of scalars are dominated (negative result)

Adjoining $s=\sqrt{h_2}$ as a formal variable ($s^2=h_2$, $s\ge0$, and
certificates that are sums of squares in $(\text{labels},s)$) looks like
it should repair the $\sqrt{h_2}$ leak directly.  It cannot, for a
structural reason:

* A square $(a+sb)^2=a^2+h_2b^2+2s\,ab$ uses only $s^2=h_2$ and
  $s\ge0$.  On the moment side, the old optimal pseudo-moment $y^*$
  extends to the $s$-variables by the rank-one assignment
  $y^*_{s\cdot m}=\sqrt{y^*_{h_2}}\;y^*_m$ whenever only moments linear
  in $s$ are constrained; every $2\times2$ extension block
  $\begin{pmatrix}1&s\\ s&h_2\end{pmatrix}$ is then exactly singular —
  the extension *saturates* Cauchy–Schwarz, which is exactly the
  behaviour the leak already exploits.  First-order use of $\sqrt{h_2}$
  is therefore identically the Cauchy–Schwarz closure of the spin-2
  Gram block, already inside the hierarchy.
* The only place the extension is *not* automatically feasible is the
  deeper identity $y_{s^2\cdot m}=y_{h_2\cdot m}$ for nonconstant $m$,
  i.e. products of $h_2$ with other labels — but those are precisely
  the `--h2-localized-flags` products, measured inert at degree 14.

So the fix the $\sqrt{}$-intuition is pointing at is not the scalar
$\sqrt{h_2}$: it is the *ratio* structure $v/\sqrt{h_2}$ near the
minimizing face, and ratios are Route A's territory (multiply through by
the denominator).  Square roots do earn their place, but as *kernels* —
see Route C.

The genuinely missing constraint remains the spin-2 rank condition (all
spin-2 correlation vectors in one five-dimensional space).  Neither
$\sqrt{h_2}$ nor any scalar function of invariants sees it, because it
is a statement about the *joint* Gram of many correlation vectors, not
about any single scalar.

## 4. Route C — entire kernels: the hypergeometric dictionary

$\exp(-c\,d^2)$ decay of the hierarchy bounds is Gaussian in the degree
filtration: the degreewise coefficients of the asymptotic certificate
behave like $q^{d^2}$ — a **theta series** in the degree.  The limit
certificate is then an *entire* kernel of order-two type in the Gram
entries, and the finite object to search for is a certificate over a
finite dictionary of such kernels.

### 4.1 The dictionary

Candidates whose Legendre expansions are exactly computable:

| kernel | closed form | Legendre coefficients | decay |
|---|---|---|---|
| Gegenbauer / Poisson | $(1-2rt+r^2)^{-1/2}$ | $r^\ell$ | geometric |
| heat | $\sum_\ell e^{-\ell(\ell+1)\tau}\frac{2\ell+1}{4\pi}P_\ell(t)$ | $e^{-\ell(\ell+1)\tau}$ | **Gaussian — matches the observed decay** |
| Bessel / ${}_0F_1$ | $e^{\kappa t}$ | $(2\ell+1)\,i_\ell(\kappa)$ | $1/\ell!$ |

(Antipodal parity: use the even parts, e.g.
$\frac12\big[G_r(t)+G_r(-t)\big]$, whose coefficients live on even
$\ell$ only.)  The Gegenbauer kernel is an *algebraic* function — a
square root of a rational function of Gram monomials — which is where
square roots re-enter usefully: as generating functions resumming the
whole degree axis, not as scalars.

### 4.2 Validity with exact arithmetic

Positive-coefficient kernels stay valid under truncation in a one-sided
way: if $G=\sum_\ell a_\ell P_\ell$ with $a_\ell\ge0$, then
$G-G_{\le D}$ is positive-definite, so the label
$k_G=\iint G(x\cdot y)\,d\mu\,d\mu$ obeys the exact linear inequalities

$$k_G\ \ge\ \sum_{\ell\le D}a_\ell\,\iint P_\ell\,d\mu\,d\mu
\qquad\text{for every }D,$$

with rational data whenever $r$ (or $\kappa$) is rational.  Adjoining
finitely many $k_G$ labels with these one-sided cuts gives a finite,
exactly-rational SDP whose certificates genuinely use infinite series —
"hypergeometric series of Gram monomials as variables", made rigorous.

### 4.3 The selection principle for kernel leaves

Not every kernel is admissible as a flag leaf: complementary slackness
forces the master circle-mode equations (E1) of
[Limit certificate](LIMIT_CERTIFICATE.md) — a leaf $g$ must restrict to
every compatible great circle with Fourier modes $\{0,\pm2,\pm6\}$ only,
with the mean tied to the pole value by (ii).  These are *linear*
conditions on Legendre coefficients, so they can be solved inside the
dictionary span: the admissible kernels form an affine subspace, and the
dictionary should be projected onto it before the SDP runs.  This is the
principled version of coefficient extrapolation: solve (E1) in closed
form, then let the SDP choose the finitely many remaining coefficients.

### 4.4 Extraction (if the SDP route stalls)

`sdpa_extract.py` already reconstructs 40-digit certificates.  Per-label
extrapolation of the degree-$d$ solutions at $d=12,\dots,20$ against the
model $c_\infty+A\,q^{d^2}$, followed by integer-relation detection
(PSLQ) over the basis $\{1,\ r^d,\ q^{d^2},\ 1/d!,\ \Gamma\text{-ratios}\}$,
identifies which dictionary row each label family is converging to.

## 5. Iterating the multiplier heuristic (second round)

Round one (weight $h_2$) collapsed the pole residue $32\times$ but the
$\varepsilon$-sweep shows $13$–$17\times$ trace growth per decade: the
pole survives, and the surviving escape is the pure $g_4$ ray
([RESULTS](../RESULTS.md), `sdpa_runs/ray_projected_deg14.json`).  The
heuristic iterates by the same recipe — measure the ray, find a valid
quantity matched to it, multiply or cut, re-measure — but the matching
conditions now do real work.

**Requirements on a second weight $Q_2$.**  (a) $Q_2\ge0$ valid for all
measures; (b) $Q_2>0$ on a dense set (reduction lemma); (c) $Q_2=0$ on
the minimizer face (vanishing-order matching); (d)
$\langle Q_2,\mathrm{ray}\rangle>0$ — the weight must *grow* to first
order along the escape so one factor absorbs a simple pole.

**Ruled out by (c)/(d).**  $h_2$ again: the ray is $h_2$-orthogonal
($\langle g_2,\mathrm{ray}\rangle=0.005$), so $h_2^2$ vanishes to the
wrong order along the wrong direction.  $h_4=\sum_m|\hat\mu_{4m}|^2$:
pairs correctly ($\langle E[P_4],\mathrm{ray}\rangle=5.353$) but fails
(c) — the octahedron (ONB) is a 3-design, not a 4-design, so $h_4>0$ at
minimizers and the weighted certificate could not be sharp there.

**The spectral-gap invariants, measured.**  The natural face-vanishing
candidates are invariants of the gap operator $B=I-A_2\succeq0$, since
$A_2$ has eigenvalue $1$ across the whole zero family.  Exact label
expansions via $\operatorname{tr}A_2^k=E[\chi_2(\rho_{x_1}\cdots\rho_{x_k})]$,
$\chi_2(R)=(\operatorname{tr}R)^2-\operatorname{tr}R-1$,
$\rho_x=2xx^{\mathsf T}-I$, and Newton's identities
(`sdpa_runs/gap_invariant_expansions.json`; note
$\operatorname{tr}B=4$ *identically*, so $e_1$ is trivial):

| invariant | ONB | pole+Haar / $m$-gon | $\langle\cdot,\mathrm{ray}\rangle$ |
|---|---:|---:|---:|
| $e_2(B)=6+6p_2-8p_4$ | $16/3$ | $52/9$ | $-9.79$ |
| $e_3(B)$ | $64/27$ | $32/9$ | $-19.26$ |
| $e_4(B)$ | $0$ | $64/81$ | $-5.85$ |

(Calibration: $\langle g_2,\mathrm{ray}\rangle=0.0051$,
$\langle E[P_4],\mathrm{ray}\rangle=5.3531$, matching the documented
$0.005$ and $5.353$.)

**Verdict: no $e_k$ qualifies as a weight** — none vanishes on the whole
face ($e_4$ vanishes only at the ONB, where $B$ has a *two*-dimensional
kernel; the continuous strata have a one-dimensional kernel and only
$e_5=\det B$ vanishes there, which needs five-sample labels outside the
four-point algebra).  All pairings are negative, the wrong sign for (d).

**But the wrong sign for a weight is the right sign for a cut.**  Each
$e_k(B)\ge0$ is a valid *linear* inequality on (product) labels, and a
functional that strictly decreases along a ray kills it: no recession
direction of $\{e_k\ge0\}$ can pair negatively.  Moreover the negative
pairing *proves* these cuts are not implied by the projected problem
(the ray was feasible there).  So the $g_4$ escape dies not by a second
denominator but by the cheapest scalar shadow of the operator gap:

$$e_2(I-A_2)=6+6p_2-8p_4\ \ge\ 0
\qquad\Longleftrightarrow\qquad p_4\ \le\ \tfrac{3(1+p_2)}{4},$$

one $1\times1$ block on two pair labels ($-9.79$ per unit ray), with
$e_3,e_4$ as independent backups.  Implemented as `--gap-scalar-cuts`.
This is the general lesson of the round-two measurement: **whether the
matched object is a weight or a cut is decided by the sign of its ray
pairing**, and the sign came out "cut".

**Measured outcome (A/B `--find-ray`, degree 14, (E1)-projected).**
Without cuts: `optimal_inaccurate`, squared norm $55.7$,
$E[P_4]\cdot\mathrm{ray}=5.3498$ — the documented $g_4$ escape,
reproduced (`sdpa_runs/ray_nocuts_deg14.json`).  With
`--gap-scalar-cuts`: **`infeasible`**
(`sdpa_runs/ray_gapcuts_deg14.json`) — not merely the $g_4$ ray but
*every* improving recession direction is gone; the projected dual is
bounded again.  Note the deep-feature operator-gap blocks
(`--spin2-operator-gap`, depth 1) had *not* achieved this; the three
scalar invariants did, at a cost of three $1\times1$ blocks.

Next per the tree: bound sweep and the $\varepsilon$-trace law with the
cuts stacked on the weighted problem
(`--h2-weighted-target --h2-localized-all --gap-scalar-cuts`, GMP), in
both the KKT and the all-measures cones.  If a *new* pole appears there,
its ray feeds back into the same recipe; if the trace plateaus,
attainment is restored and the §4-of-PLAN exactification pipeline
(max-margin, rational rounding, independent verification) runs on the
weighted-and-cut problem.  If successive rays were ever to rotate
through the $g_\ell$ tower without terminating, that would be the
signature that scalar cuts cannot close it, and the theta/kernel
dictionary (§4) takes over as the systematic completion.

## 6. Work plan

1. **(running)** Solve the weighted target $h_2E$ at degree 14, pruned
   nine-toggle hierarchy, SDPA-GMP 200-bit — compare against the
   unweighted $-4.4856\times10^{-3}$.
2. Ablate `--h2-localized-flags` in; sweep degree 16/18 for the winner.
3. $\varepsilon$-trace law on the weighted target (attainment
   diagnostic, §2.4).  Bounded trace ⇒ max-margin ⇒ rationalization
   pipeline unchanged.
4. If the pole survives: Frobenius selector, measure $Y_0$, pick the
   next denominator factor (§2.5).
5. Fallback: kernel dictionary (§4), starting from the heat kernel at
   two rational temperatures projected onto the (E1) subspace.
