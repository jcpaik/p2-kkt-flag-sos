# Enrichments: weights, cuts, and block families beyond the plain hierarchy

This document collects every enrichment of the certificate algebra —
multiplier weights, scalar and matrix cuts, atom adjunctions, and new
block families — with, for each: the definition, the validity proof,
the exact label-expansion facts, and a one-or-two-sentence measured
verdict.  Detailed measurement narratives and run logs live in
[SUBCASES_AND_RECORD.md](SUBCASES_AND_RECORD.md); definitions of
labels, flags, blocks, and the export pipeline are in
[FOUNDATIONS.md](FOUNDATIONS.md); the sharp-ansatz/(E1) structure is
in [SHARP_STRUCTURE.md](SHARP_STRUCTURE.md).

Throughout: $\mu$ ranges over antipodal probability measures on $S^2$,
$p_j=\iint(x\cdot y)^j\,d\mu\,d\mu$, the target is
$E=-\tfrac43+20p_2-48p_4+32p_6$, and

$$h_2=\frac{3p_2-1}{2}=\sum_{m=-2}^{2}\big|\hat\mu_{2m}\big|^2\ \ge\ 0$$

is the spin-2 harmonic energy, zero exactly on second-moment-isotropic
measures ($p_2=\operatorname{tr}\Sigma^2\in[1/3,1]$, so also
$0\le h_2\le1$ on probability measures).  A *label* is an isomorphism
class of an edge-weighted graph on $k$ vertices standing for the
moment
$\int\!\cdots\!\int\prod_{i<j}(x_i\cdot x_j)^{w_{ij}}\,d\mu^{\otimes k}$;
the relaxation replaces $\mu$ by a pseudo-moment vector
$y=(y_\ell)$ constrained by PSD blocks $A_b(y)\succeq0$ and linear
equalities, and by SDP duality every bound is certified by PSD
multipliers $Q_b$ with
$E-\text{bound}=\sum_b\langle Q_b,A_b(\cdot)\rangle+(\text{equality
multipliers})$ as an identity in the label algebra.  A product of
moments over disjoint sample sets is itself a label (a disconnected
graph) — this "disconnected-label mechanism" is used by every weight
below.

---

## 1. The multiplier program

### 1.1 Non-attainment and what the pole says

The assumption-free four-point hierarchy converges to zero
super-exponentially: bounds decay like $\exp(-c\,d^2)$ in the degree
$d$, while the minimal-trace certificate of $E+\varepsilon\ge0$
diverges like $\operatorname{tr}C(\varepsilon)\approx1.07/\varepsilon$
(legacy $(3/16)E$ scale; the pole *order* is scale-free) — a
**simple pole**.  Together with heavy Gram entries growing
$10^6\to10^{11}$ when inactive families are removed and the $L^2$
blow-up criterion, this is the signature of **non-attainment**: the
target $E$ lies on the boundary of the closure of the finite-degree
certificate cone but not in the cone itself.  Concretely
$\varepsilon C_\varepsilon\to Y_0$ with $Y_0\succeq0$, $Y_0\neq0$,
$e(Y_0)=0$: a recession direction of the certificate cone along which
the sharp certificate escapes to infinity.  This is the classical
Motzkin mechanism by which a nonnegative polynomial fails to be a sum
of squares while $q\cdot f$ is one for a suitable positive multiplier
$q$ — and a *simple* pole predicts a multiplier vanishing to *first*
order along the leak.

The leak variable is known from the moment diagnostics: the optimal
pseudo-moment pays $h_2$ linearly and earns contraction violations of
order $\sqrt{h_2}$.  The multiplier vanishing to first order along
the leak is $h_2$ itself.

On the moment side, non-attainment of the dual is equivalent to
unboundedness of a *projected* primal, witnessed by a **recession
direction**: a vector $r$ with

$$r_{\mathrm{const}}=0,\qquad A_b(r)\succeq0\ \ \forall b,\qquad
\text{equalities}(r)=0,\qquad \langle E,r\rangle=-1 ,$$

so that $y+tr$ is feasible for all $t\ge0$ with
$\langle E,y+tr\rangle\to-\infty$.  `--find-ray` computes the
least-norm such $r$ on the (E1)-projected problem.  The documented ray
(stored in the legacy $(3/16)E$ normalization) pairs as
$\langle\iint P_4,r\rangle=5.35$ with all other harmonic energies
pairing to $0$ — the escape concentrates in the $\ell=4$ mode, the
unique negative Legendre coefficient of $K$.

### 1.2 The reduction lemma

**Lemma.** If $h_2(\mu)\,E(\mu)\ge0$ for every antipodal probability
measure $\mu$, then $E(\mu)\ge0$ for every antipodal probability
measure.

*Proof.* On $\{h_2>0\}$ divide.  If $h_2(\mu)=0$, let
$\mu_t=(1-t)\mu+t\,\delta_{\pm e}$; its second moment
$\frac{1-t}{3}I+t\,ee^{\mathsf T}$ is anisotropic for every
$t\in(0,1]$, so $E(\mu_t)\ge0$, and $E$ is weak-\* continuous
(polynomial kernel), so
$E(\mu)=\lim_{t\downarrow0}E(\mu_t)\ge0$. $\blacksquare$

A certificate for the weighted target therefore proves the full
conjecture; the isotropic stratum comes for free by density, with no
separate isotropic-branch argument.  The lemma survives for any
weight that is nonnegative on all measures and strictly positive on a
weak-\* dense class (see §3, Fact 3, for the general density
argument).

### 1.3 The weighted target $h_2E$ and its polynomial expansion

$h_2E$ is polynomial in the existing label algebra:

$$h_2E=\Big(\tfrac{3p_2-1}{2}\Big)
\Big(-\tfrac43+20p_2-48p_4+32p_6\Big)
=\tfrac23-12p_2+24p_4-16p_6
+30\,p_2p_2-72\,p_2p_4+48\,p_2p_6 ,$$

where each product $p_2p_j$ is the expectation of a *disconnected*
four-sample graph — two independent pairs — already represented as a
`("product", ("pair",2), ("pair",j))` label.  Implemented as
`--h2-weighted-target` (a different rational objective vector; the
rest of the export pipeline is unchanged):

```sh
python3 sos_search.py \
  --export-sdpa PROBLEM.dat-s --h2-weighted-target \
  --degree DEGREE --no-pointwise-sos \
  --harmonics --three-point-flags --four-point-flags --two-root-flags \
  --gradient --potential --hessian \
  --global-tangent-gaps --rank-relations
```

### 1.4 The localized module, and the mixing design rule

The relaxation *linearizes* products: the label $y_{p_2p_j}$ is a
free variable, not the product $y_{p_2}y_{p_j}$.  With no constraints
tying product labels to their factors, the pseudo-moment decouples
the weight from the target and most of the intended cancellation
never happens.  (The one coupling always present — the
`empty_type_flag` block, the moment matrix over $\{1,p_2,p_4,p_6\}$ —
supplies the one-sided inequality $y_{p_2p_2}\ge y_{p_2}^2$.)

The natural quadratic module for a weighted certificate contains
$h_2\times(\text{flag square})$ terms.  Since $h_2=\tfrac{3p_2-1}2$,
the $h_2$-multiplied copy of a block with coefficient matrices
$\{A_\ell\}$ has entries

$$h_2\times A_\ell\;\longmapsto\;
\tfrac32\,A_\ell\ \text{on}\ y_{p_2\cdot\ell}
\;\;-\;\tfrac12\,A_\ell\ \text{on}\ y_{\ell},$$

and the $-\tfrac12$ part is the whole point: every localized block
couples each product label to its base label with a definite sign, so
a pseudo-moment that depresses $y_{p_2\cdot\ell}$ below the true
product is caught by the PSD constraint through the $y_\ell$ term.
`--h2-localized-flags` provides the flag-family subset;
`--h2-localized-all` adds the $h_2$-multiplied copy of **every** PSD
family and the $p_2$-multiplied copy of every scalar equality family
(for equalities the constant part is redundant).  On the certificate
side this is exactly what non-attainment demands: certificates of the
form $\sigma_0+h_2\sigma_1$ ($\sigma_i$ in the ordinary flag-square
cone) become expressible at finite size — the degree-$d$ truncation
of the quadratic module generated by $h_2$.

**Negative control (the mechanism proof).**  Multiplying every PSD
block by $p_2$ instead is valid ($p_2\ge\tfrac13>0$ always,
disconnected samples factor) and doubles the problem size — yet the
bound is identical to the unlocalized weighted bound to eleven
digits.  A $p_2$-multiplied block touches *only* fresh product labels
$y_{p_2\cdot\ell}$, which appear in no other constraint and not in
the objective, so the pseudo-moment satisfies the new block trivially.
**Design rule**: a localizing block is worthless unless its entries
mix the new product labels with the old base labels inside one PSD
constraint; the affine part of the multiplier is what gives the
relaxation teeth.  Any multiplier with zero constant term produces
vacuous blocks under linearization.

**Verdict.**  At degree 14 (same degree, arity, and label set),
`--h2-weighted-target --h2-localized-all` moved the bound from
$-4.4856\times10^{-3}$ to $-2.7909\times10^{-5}$ — $161\times$, worth
more than two full degree steps ($2.7\times$ better than unweighted
degree 18 with 30% fewer variables) and attributable by the controls
entirely to the $-\tfrac12$ coupling.  The pole coefficient collapsed
$32\times$ but the $\varepsilon$-trace still grows ($13$–$17\times$
per decade at first measurement): one factor of $h_2$ tames but does
not kill the pole, and the surviving escape is the pure $g_4$ ray.

**The attainment diagnostic** (used for every weight iteration):
solve the minimal-trace selector of the weighted target
$+\varepsilon\ge0$ for decreasing $\varepsilon$.
$\operatorname{tr}C_w(\varepsilon)=O(1)$ means attainment is restored
and the max-margin/rational-rounding exactification pipeline applies
unchanged; $\operatorname{tr}C_w(\varepsilon)\sim c/\varepsilon$
means the pole survives and the new residue $Y_0^{(w)}$ names the
next factor.

### 1.5 Requirements (a)–(d) on further weights

A second multiplier $Q$ in an iterated rational certificate
$Q\,h_2\,E=\text{SOS}$ must satisfy:

* **(a)** $Q\ge0$ valid for every measure;
* **(b)** $Q>0$ on a weak-\* dense set (reduction lemma);
* **(c)** $Q=0$ on the whole zero-energy family — otherwise near a
  minimizer $\mu^*$ with $Q(\mu^*)>0$ the product vanishes to the
  same order as $E$ and inherits the same representability
  obstruction;
* **(d)** $\langle\mathrm{lin}(Q),\,r\rangle>0$ — the weight must
  grow to first order along the escape so one factor absorbs a simple
  pole (this is what $h_2$ did to its own leak).

**Ruled out by (c)/(d).**  $h_2$ again: the ray is $h_2$-orthogonal
($\langle g_2,\mathrm{ray}\rangle=0.005$), so $h_2^2$ vanishes to the
wrong order along the wrong direction.
$h_4=\sum_m|\hat\mu_{4m}|^2$: pairs correctly
($\langle E[P_4],\mathrm{ray}\rangle=5.353$) but fails (c) — the
octahedron (ONB) is a 3-design, not a 4-design, so $h_4>0$ at
minimizers and the weighted certificate could not be sharp there.
Denominator candidates beyond $h_2$, in order: $h_2^2$ (if the pole
on the weighted target were again order one);
$h_2+\beta h_4+\gamma h_6$ with $h_\ell$ the higher harmonic energies
(all nonnegative labels, if the residue loads on spin-4/spin-6
blocks); the "squared-form" weight
$\iint T_3(x\cdot y)^2-\tfrac13+h_2$ (nonnegative only *jointly* by
the frame-form identity, use with care).  The spectral-gap invariants
$e_k(I-A_2)$ were the measured face-vanishing candidates — see §2
(they fail as weights and become cuts) and §3 (the one that succeeds
as a weight, $e_5$).

### 1.6 Negative result: square-root adjunction is dominated

Adjoining $s=\sqrt{h_2}$ as a formal variable ($s^2=h_2$, $s\ge0$,
certificates sums of squares in $(\text{labels},s)$) looks like it
should repair the $\sqrt{h_2}$ leak directly.  It cannot, for a
structural reason:

* A square $(a+sb)^2=a^2+h_2b^2+2s\,ab$ uses only $s^2=h_2$ and
  $s\ge0$.  On the moment side, the old optimal pseudo-moment $y^*$
  extends to the $s$-variables by the rank-one assignment
  $y^*_{s\cdot m}=\sqrt{y^*_{h_2}}\;y^*_m$ whenever only moments
  linear in $s$ are constrained; every $2\times2$ extension block
  $\begin{pmatrix}1&s\\ s&h_2\end{pmatrix}$ is then exactly singular
  — the extension *saturates* Cauchy–Schwarz, exactly the behaviour
  the leak already exploits.  First-order use of $\sqrt{h_2}$ is
  identically the Cauchy–Schwarz closure of the spin-2 Gram block,
  already inside the hierarchy.
* The only place the extension is *not* automatically feasible is the
  deeper identity $y_{s^2\cdot m}=y_{h_2\cdot m}$ for nonconstant $m$
  — precisely the `--h2-localized-flags` products, measured inert at
  degree 14.

The fix the $\sqrt{}$-intuition points at is the *ratio* structure
$v/\sqrt{h_2}$ near the minimizing face, and ratios are the
multiplier program's territory.  Square roots earn their place as
*kernels* (§1.7), not scalars.  The genuinely missing constraint
remains the spin-2 rank condition (all spin-2 correlation vectors in
one five-dimensional space): a statement about the joint Gram of many
correlation vectors, invisible to any scalar function of invariants.

### 1.7 Route C: the kernel dictionary and the one-sided truncation principle

$\exp(-c\,d^2)$ decay is Gaussian in the degree filtration: the
degreewise coefficients of the asymptotic certificate behave like
$q^{d^2}$ — a theta series in the degree.  The limit certificate is
then an *entire* kernel of order-two type in the Gram entries, and
the finite object to search for is a certificate over a finite
dictionary of such kernels.

**The dictionary** (Legendre expansions exactly computable; antipodal
parity: use even parts, e.g. $\tfrac12[G_r(t)+G_r(-t)]$):

| kernel | closed form | Legendre coefficients | decay |
|---|---|---|---|
| Gegenbauer / Poisson | $(1-2rt+r^2)^{-1/2}$ | $r^\ell$ | geometric |
| heat | $\sum_\ell e^{-\ell(\ell+1)\tau}\frac{2\ell+1}{4\pi}P_\ell(t)$ | $e^{-\ell(\ell+1)\tau}$ | **Gaussian — matches the observed decay** |
| Bessel / ${}_0F_1$ | $e^{\kappa t}$ | $(2\ell+1)\,i_\ell(\kappa)$ | $1/\ell!$ |

The Gegenbauer kernel is an *algebraic* function — a square root of a
rational function of Gram monomials — which is where square roots
re-enter usefully: as generating functions resumming the whole degree
axis.

**One-sided truncation principle (exact validity).**
Positive-coefficient kernels stay valid under truncation in a
one-sided way: if $G=\sum_\ell a_\ell P_\ell$ with $a_\ell\ge0$, then
$G-G_{\le D}$ is positive-definite, so the adjoined label
$k_G=\iint G(x\cdot y)\,d\mu\,d\mu$ obeys the exact linear
inequalities

$$k_G\ \ge\ \sum_{\ell\le D}a_\ell\,\iint P_\ell\,d\mu\,d\mu
\qquad\text{for every }D,$$

with rational data whenever $r$ (or $\kappa$) is rational.  Adjoining
finitely many $k_G$ labels with these one-sided cuts gives a finite,
exactly-rational SDP whose certificates genuinely use infinite series.
(The theta atoms of §4 are the realized instance of this principle,
on the modulation axis rather than the degree axis.)

**Selection principle.**  Not every kernel is admissible as a flag
leaf: complementary slackness forces the master circle-mode equations
(E1) of [SHARP_STRUCTURE.md](SHARP_STRUCTURE.md) — a leaf must
restrict to every compatible great circle with Fourier modes
$\{0,\pm2,\pm6\}$ only, the mean tied to the pole value.  These are
*linear* conditions on Legendre coefficients, so the admissible
kernels form an affine subspace of the dictionary span; project the
dictionary onto it before the SDP runs.

**Extraction fallback.**  Per-label extrapolation of the degree-$d$
solutions at $d=12,\dots,20$ against the model $c_\infty+A\,q^{d^2}$,
followed by integer-relation detection (PSLQ) over the basis
$\{1,\ r^d,\ q^{d^2},\ 1/d!,\ \Gamma\text{-ratios}\}$, identifies
which dictionary row each label family converges to
(`sdpa_extract.py` reconstructs 40-digit certificates).

### 1.8 The multi-weight program is one convex problem

Because the target $E$ is fixed, a multiplier drawn from a finite
dictionary $q=\sum_j\lambda_j w_j$, $\lambda_j\ge0$, enters the
certificate identity *linearly*: each $w_jE$ is a fixed rational
label vector, so

$$e\Big(\sum_j\lambda_j\,w_jE\Big)-e(\sigma)=0,\qquad
\sigma\in\text{(certificate cone)},\ \lambda\ge0,\ \textstyle\sum_j\lambda_j=1$$

is jointly convex in $(\lambda,\sigma)$ — one SDP, not an
alternation.  The hand-iterated weight search is the special case of
singleton dictionaries; requirement (d) is then *selected
automatically* by the SDP instead of matched by hand.  The bottleneck
is the dictionary: the invariants satisfying (a)–(c) are $h_2$ and
$e_5(I-A_2)$ (§3).

---

## 2. The sign rule, and the $e_k(I-A_2)$ scalar gap cuts

### 2.1 The operator bound and the exact expansions

Embed $\mathbb{RP}^2\hookrightarrow SO(3)$ by
$x\mapsto\rho_x=2xx^{\mathsf T}-I$ (rotation by $\pi$ about $x$) and
let $\pi_2$ be the 5-dimensional spin-2 representation.  Then
$A_2=\int\pi_2(\rho_x)\,d\mu$ is a symmetric contraction
($\pi_2(\rho_x)$ orthogonal, so $\|A_2\|\le1$), hence
$B=I-A_2\succeq0$ is a valid constraint for every measure; and at
every zero-energy measure $A_2$ has eigenvalue exactly $1$
(eigenvector: the axial quadrupole $\mathrm{diag}(1,1,-2)$).  The
natural scalar invariants are the elementary symmetric polynomials
$e_k(B)\ge0$ of its eigenvalues.

Exact label expansions: the character identity
$\chi_2(R)=\tau^2-\tau-1$, $\tau=\operatorname{tr}R$, gives
$\operatorname{tr}(A_2^k)=E\big[\chi_2(\rho_{x_1}\cdots\rho_{x_k})\big]$,
and $\operatorname{tr}(\rho_{x_1}\cdots\rho_{x_k})$ expands into
cyclic dot-product monomials

$$\operatorname{tr}(\rho_{x_1}\cdots\rho_{x_k})
=\sum_{S\subseteq[k]} 2^{|S|}(-1)^{k-|S|}
\prod_{\text{cycle of }S}(x_i\cdot x_j),$$

so every $\operatorname{tr}(A_2^k)$ is an exact rational combination
of $k$-vertex labels, and Newton's identities give the $e_k$
(expansions in `sdpa_runs/gap_invariant_expansions.json`).
Byproduct: $\operatorname{tr}B=4$ **identically** (since
$\chi_2(\rho_x)=1$), so $e_1$ is trivial.  For example

$$e_2(I-A_2)=6+6p_2-8p_4\ \ge\ 0
\qquad\Longleftrightarrow\qquad p_4\ \le\ \tfrac34\,(1+p_2).$$

### 2.2 Face values and ray pairings

Measured exactly (face values) and against the documented projected
ray $r$ (pairings calibrated by $\langle g_2,r\rangle=0.0051$ and
$\langle\iint P_4,r\rangle=5.3531$, both reproducing the documented
$0.005$ and $5.353$):

| invariant | ONB | pole+Haar / $m$-gon strata | $\langle\cdot,r\rangle$ |
|---|---:|---:|---:|
| $e_2(B)=6+6p_2-8p_4$ | $16/3$ | $52/9$ | $-9.79$ |
| $e_3(B)$ | $64/27$ | $32/9$ | $-19.26$ |
| $e_4(B)$ | $0$ | $64/81$ | $-5.85$ |

**No $e_k$, $k\le4$, qualifies as a weight**: every candidate fails
(c) — on the continuous strata $B$ has a one-dimensional kernel, so
only $e_5=\det B$ vanishes identically on the face ($e_4$ vanishes
only at the ONB, where $B$ has a *two*-dimensional kernel), and
$e_5$ needs five-sample labels outside the four-point algebra (§3) —
and every candidate fails (d), with strictly negative pairings.

### 2.3 The sign rule

**Whether a matched valid quantity is a weight or a cut is decided by
the sign of its ray pairing.**  A weight must satisfy (d): positive
pairing, so a factor can absorb the pole.  A valid linear inequality
$\langle q,y\rangle\ge0$ appended as a $1\times1$ block adds
$\langle q,r\rangle\ge0$ to every recession direction; a functional
that strictly *decreases* along a ray therefore kills it — **negative
pairing kills the ray**.  Moreover the negative pairing is
simultaneously a *separation certificate*: the inequality cannot be
implied by the existing constraint set, or the ray could not have
been feasible.  Operationally: before any big solve, pair the
candidate against the current escape data; positive → weight,
negative → cut.

### 2.4 The gap scalar cuts and their verdicts

Each $e_k(B)\ge0$, expanded in labels, is a linear inequality in $y$
(products of moments over disjoint samples are disconnected labels).
Implemented as `--gap-scalar-cuts` (three $1\times1$ blocks:
$e_2,e_3,e_4$).

**Verdict (projected ray).**  A/B on the (E1)-projected degree-14
`--find-ray` problem: without cuts, feasible, $\|r\|^2=55.7$,
$\langle\iint P_4,r\rangle=5.3498$
(`sdpa_runs/ray_nocuts_deg14.json`); with the three cuts,
**infeasible** (`sdpa_runs/ray_gapcuts_deg14.json`) — not merely the
$g_4$ ray but *every* improving recession direction is gone, at the
cost of three $1\times1$ blocks.  The depth-1 operator-gap localizing
blocks `--spin2-operator-gap` had *not* achieved this.

**Verdict (finite-$\varepsilon$ optima).**  On the projected
all-measures cone at $\varepsilon=0$ the gap cuts restore boundedness
(dual $-3.18615$, from unbounded) but that value is catastrophically
weak against the unprojected cone ($-6.1\times10^{-3}$): recession
repair, not certificate strength.  On unprojected weighted problems
the cuts are slack at the optimum — bounds and selector traces
unchanged — the recurring slack-at-optimum pattern of every scalar
cut measured (see also §4, §3).

---

## 3. The $e_5(I-A_2)$ weight and the multi-weight dictionary

### 3.1 The three facts

With $A_2$ acting on traceless symmetric $S$ by
$A_2[S]=\int\rho_xS\rho_x\,d\mu$ and $B=I-A_2\succeq0$:

**Fact 1 (validity, (a)).** $e_5(B)=\det(I-A_2)\ge0$ for every
measure.  *(All eigenvalues of $B$ lie in $[0,2]$.)*

**Fact 2 (face-vanishing, (c) — strengthened).** $\det(I-A_2)=0$ for
**every** pole–equator measure $\mu=w\,\delta_{\pm e}+(1-w)\nu$, for
*all* $w$ and *all* circle measures $\nu$ — not only the zero family.
*Proof.*  The axial quadrupole $Q=I-3ee^{\mathsf T}$ obeys
$\rho_xQ\rho_x=Q$ for $x=\pm e$ (both diagonal in an adapted basis)
and for every equatorial $x$ (there $\rho_xe=-e$, hence
$\rho_x(I-3ee^{\mathsf T})\rho_x=I-3ee^{\mathsf T}$).  So $A_2[Q]=Q$
for any measure supported on $\{\pm e\}\cup e^\perp$, i.e.
$1\in\operatorname{spec}A_2$ and $\det B=0$. $\square$
Verified symbolically for the parametrized family
($w,\ \hat\nu(2),\ \hat\nu(4)$ free): `sdpa_runs/e5_face_check.py`
V1/V4 give $\det(I-A_2)\equiv0$ identically in all three parameters,
alongside exact reproduction of the documented spectra
(ONB: $\{1,1,-\tfrac13,-\tfrac13,-\tfrac13\}$; $e_{1..5}(B)$ at ONB
$=4,\ \tfrac{16}3,\ \tfrac{64}{27},\ 0,\ 0$; at pole+Haar
$=4,\ \tfrac{52}9,\ \tfrac{32}9,\ \tfrac{64}{81},\ 0$).

**Fact 3 (dense positivity, (b)).** $\{e_5>0\}$ is weak-\* dense.
*Proof.*  $A_2$ is affine in $\mu$ and
$A_2(\text{uniform})=\tfrac15I$ (the trace is
$\chi_2(\rho_x)\equiv1$, and irreducibility under the rotation
average forces the scalar).  For any $\mu$ set
$\mu_t=(1-t)\mu+t\,\text{uniform}$: then
$B(\mu_t)=(1-t)B(\mu)+\tfrac{4t}5I\succeq\tfrac{4t}5I\succ0$, so
$e_5(\mu_t)>0$ for all $t>0$ and $\mu_t\to\mu$. $\square$
The same argument re-proves the reduction lemma for any weight of the
form $q_\kappa=h_2+\kappa\,e_5(B)$, $\kappa\ge0$.

**Why $e_5$ and not $e_2,e_3,e_4$**: on the continuous strata of the
zero family $B$ has a one-dimensional kernel, so $e_5$ is the *first*
elementary symmetric invariant forced to vanish there — **the first
face-vanishing scalar shadow of the operator gap lives at word length
five**, which explains structurally why the depth-1 operator-gap
blocks and the scalar $e_k$ cuts could not serve as weights.

The combined weight $q_\kappa=h_2+\kappa\,e_5(B)$ has zero set
$\{h_2=0\}\cap\{e_5=0\}$ = isotropic measures with an
$A_2$-invariant quadrupole — strictly smaller than $\{h_2=0\}$, so
the weighted-(E1) forced-vanishing conditions for $q_\kappa E$ are
weaker than for $h_2E$ (certificates get more room) while the
reduction lemma still applies.

### 3.2 The arity obstruction and the exact 35-label expansion

$e_5$ expands via
$\operatorname{tr}A_2^k=E[\chi_2(\rho_{x_1}\cdots\rho_{x_k})]$,
$k\le5$; the genuinely five-sample cycle labels of
$\operatorname{tr}A_2^5$ are outside the four-sample algebra.  The
obstruction is **label arity, not polynomial degree**: every
connected label of $e_5$ has per-vertex degree $\le4$.  The reducer
(`reduce_graph_matrix` / `canonical_connected_label`) canonicalizes
connected multigraphs on any number of vertices, so no reducer change
was needed; the missing pieces were the expansion itself, an
objective/cut path, and coverage blocks.

Implementation (`sos_search.py`): `gap_power_trace_vector(k)`
computes $\operatorname{tr}A_2^k$ by the cyclic expansion
$\operatorname{tr}(\rho_{x_1}\cdots\rho_{x_k})
=\sum_{S\subseteq[k]}(-1)^{k-|S|}2^{|S|}\operatorname{cyc}(S)$;
`a2_elementary_vector(k)` applies Newton's identities
(label-vector products = disconnected products,
`multiply_label_vectors`); `gap_elementary_vector(k)` returns
$e_k(B)=\sum_j(-1)^j\binom{5-j}{k-j}e_j(A_2)$.  All arithmetic is
`Fraction`-exact.

**The $e_5$ vector**: 35 labels — 1 constant, 2 pair, 4 triangle, 6
`graph_4`, **11 `graph_5`** (the genuinely five-sample cycle sector,
coefficients like $-512$, $1280$, $-1024/5$), 3 pair×pair and 8
pair×triangle products (`sdpa_runs/e5_invariant_expansion.json`).

**Machine checks (all exact; `sdpa_runs/e5_label_check.py`,
regression tests in `test_sos_search.py`):**

| check | result |
|---|---|
| `gap_elementary_vector(2,3,4)` vs the audited `--gap-scalar-cuts` tables | identical dictionaries |
| pole–equator family, symbolic in $(w,\hat\nu(2),\hat\nu(4),\hat\nu(6))$ via the label expansion (V1) | identically $0$ |
| pole + regular equators (orders 3,4,5,6) and ONB, exact evaluators | $0$ |
| uniform measure (`uniform_label_value`, exact $S^2$ Wick evaluator) | $(4/5)^5=1024/3125$ |
| 3-atom rational measures (three weight sets), expansion vs direct $5\times5$ $\det(I-A_2)$ over $\mathbb Q$ | equal |
| $e_4$ on the same family (V3) | $(1-r^2)(1-w^2)^2\,(\cdot)\neq0$: $e_5$ is the *first* vanishing invariant |

Coverage note: five-sample flag squares already reach the cycles —
with `--max-flag-arity 5 --max-root-factor-degree 2` the
`weighted_flag_5_*` family generates 554 distinct `graph_5` labels at
degree 14, and **all 11** five-sample labels of $e_5$ are among them
(a 5-cycle embeds as a two-leaf double star with one root–root edge;
the doubled 5-cycle of $(\operatorname{tr}P)^2$ needs total root-root
degree 2).  The $A_2$-word blocks
$E[v^{\mathsf T}(I-A_2)v]$, $v=A_2^d[\text{spin-2 flag}]$, $d\le4$,
remain the escalation if the arity-5 family proves too weak
(Cayley–Hamilton ties $\det B$ to $B$-words of length $\le5$ against
the quadrupole eigendirection; by Fact 2 such blocks are
sharpness-compatible on the whole pole–equator stratum).

### 3.3 The two-sided cut and the coverage module

* **Cut** `--gap-cut-e5`: $e_5(B)\ge0$ as a $1\times1$ block
  (35 labels).  Validity: all-measures (eigenvalues of $B$ in
  $[0,2]$).
* **Upper cut** `gap_cut_e5_upper`: $(4/5)^5-e_5\ge0$.
  *Proof (AM–GM):* the five eigenvalues of $B$ are nonnegative with
  $\operatorname{tr}B=4$ identically, so
  $e_5=\prod\lambda_i\le(\tfrac{\sum\lambda_i}5)^5=(4/5)^5$, with
  **equality exactly at the uniform measure** (test-asserted).  Its
  recession form forces $e_5(r)=0$ on any dual ray.
* **Weighted target** `--e5-weight KAPPA` (rational $\kappa\ge0$,
  requires `--h2-weighted-target`): target
  $(h_2+\kappa e_5(B))\,E$.  The $e_5E$ part is the disconnected
  product of the 35 $e_5$ labels with the four energy terms
  ($-4/3,20p_2,-48p_4,32p_6$) — six samples counting the $E$
  factor's pair.  Assembled in exact rationals and floated once; a
  regression test verifies every coefficient survives
  `rationalize_float` round-trip for $\kappa\in\{1/4,1,4\}$.
* **Coverage** `--e5-localized-harmonics` (implied by `--e5-weight`,
  as is the cut), designed by exact recession analysis:
  * `e5loc_hankel_even/odd`:
    $e_5\cdot E[v(t)v(t)^{\mathsf T}]\succeq0$, $t=X\!\cdot\!Y$,
    $v=(1,t^2)$ resp. $(t,t^3)$ — a nonnegative scalar invariant
    times a pair moment matrix (disconnected samples factor).  Given
    $e_5(r)=0$ (from the upper cut), the two Hankels force
    $e_5p_2(r)=e_5p_4(r)=0$, leaving only the $+32\kappa\,e_5p_6$
    objective direction ($\ge0$): the $\kappa e_5E$ part cannot open
    a recession ray.  (With $1\times1$ combos only, the
    $\kappa$-target dual is unbounded at every $\kappa$; with the
    Hankel + AM–GM module it is bounded.)
  * `e5loc_harmonic_d`, `e5comp_harmonic_d` ($d=2,4,6$):
    $e_5\cdot E[P_d]\ge0$ and $((4/5)^5-e_5)\cdot E[P_d]\ge0$.

  Together with the cut these contain every label of
  $\kappa e_5E$ (test-asserted).  Under `--h2-localized-all` every
  e5 block acquires an $h_2$-multiplied copy automatically
  ($h_2\cdot e_5\cdot(\cdot)\ge0$, same validity).  Everything here
  is all-measures valid; nothing is KKT-only.

Exporter hardening (motivated by this module): the exact
image-elimination in `export_sdpa_problem` tracks the objective
residual of every dropped (image-dependent) direction; a mismatch
means the true dual is unbounded and the exported file would hide it
— counted as `dropped_objective_inconsistent` and warned on stderr.

### 3.4 Off-stratum growth and the escape pairing

Exact sympy series of $\det(I-A_2)$ for model families leaving the
pole–equator stratum (`sdpa_runs/e5_perturbations.py`; the stratum
itself is $e_5$-flat in *all* directions — $\det\equiv0$ identically
in $(w,\hat\nu(2),\hat\nu(4))$):

| family (off-stratum motion) | leading order of $e_5$ |
|---|---|
| F1 split latitudes: equator $\to$ Haar rings at $\pm u$ | $\dfrac{256}{81}u^2$ |
| F2 rings with mode-4 density $r$ | $\dfrac{256}{81}(1-r^2)u^2$; at the ONB corner $r=\pm1$: $+\dfrac{1024}{81}u^4$ |
| F3 equator circle tilted by $t$ | $\dfrac{64}{81}t^2$ |
| F4 pole mass spread to a polar ring (angle $u$) | $\dfrac{128}{81}u^2$ |
| F5 control: any on-stratum $w$ | identically $0$ |

$e_5$ is strictly positive-quadratic transverse to the stratum in
every model direction tested — including the fiber-deterministic
motions (F1 is a deterministic-$|z|$-fiber deformation, F4 the
pole-side analogue) — degenerating to quartic-but-positive only at
the isolated mode-4 corner $r=\pm1$.  So $e_5$'s zero set excludes
exactly the saturation face of §5.

**Escape pairing** (the $\le4$-sample visible shadow; the 11
`graph_5` labels are absent from the stored problems and pair as 0):
the deep selector growth direction pairs
$\langle e_5^{\le4},D\rangle=-1.26\times10^{5}$ — **cut-signed** by
the sign rule: the visible shadow of $e_5\ge0$ excludes the deep
escape direction.  (The weight role rests on the face geometry
above, not on this pairing.)

**Verdict.**  In double precision the e5 cut module alone improved
the all-measures weighted bound $-1.9198\times10^{-4}\to
-1.8320\times10^{-4}$ (4.6%) — the first *scalar* invariant to move
this bound at all — but at GMP precision the bare e5 cut is **inert**
on the degree-14 weighted problem (dual bound and $10^{-4}$ selector
identical to baseline).  The $\kappa$-weighted targets are bounded
but uncontrolled ($\approx-0.07\kappa$; arity-5 weighted flags
improve $\kappa=1/4$ by 22% while leaving the $h_2E$ control
unchanged): the weight route needs the missing
**$e_5$-localized module** ($e_5\times$flag blocks, label count
$\approx35\times$ per family — a degree-16-scale job), the
$e_5$-sector analogue of $h_2E$ before `--h2-localized-all` existed.

### 3.5 The $W$-KKT re-encoding

For $W=h_2E$ the first variation at $\mu$ in direction $\nu-\mu$ has
density

$$\Phi^W_\mu(x)\;=\;3\,(x^{\mathsf T}M_2(\mu)\,x)\,E(\mu)\;+\;2\,h_2(\mu)\,U_\mu(x),
\qquad U_\mu(x)=\int K(x\cdot y)\,d\mu(y),$$

since $h_2=(3p_2-1)/2$ gives $h_2'(\mu)(x)=3\,x^{\mathsf T}M_2x$.
The KKT/first-variation inequality $\Phi^W_\mu\ge\lambda$ (equality
on $\operatorname{supp}\mu$) is expressible in the existing label
algebra: $(x^{\mathsf T}M_2x)\cdot E$ terms are disconnected products
of pair labels with rooted spin-2 flags (root + 3 samples), and
$h_2\cdot U$ terms are $h_2$-localized potential flags — the same
disconnected-label mechanism as `--h2-localized-all`.  **Composition
is sound**: a certificate over the $W$-KKT cone proves $W\ge0$ at
every $W$-critical measure; $W$ attains its minimum (weak-\*
compactness, polynomial continuity), minimizers are critical, hence
$W\ge0$ everywhere, hence $E\ge0$ by the reduction lemma.  This
upgrades KKT-inclusive numbers to proof-carrying status once the four
KKT toggles (gradient, potential, Hessian, tangent gaps) are
re-derived for $W$ with the $\Phi^W$ density above.

---

## 4. Theta atoms: entire-kernel generators on the modulation axis

### 4.1 Setting

For a polynomial $A(t_1,t_2,s)$ the **two-root quadratic form** is

$$\mathcal Q[A](\mu)\;=\;\iint\Big(\int A(x_1\cdot y,\;x_2\cdot y,\;x_1\cdot x_2)\,d\mu(y)\Big)^{2}\,d\mu(x_1)\,d\mu(x_2)\;\ge\;0,$$

a $1\times1$ two-root Gram block; expanding the square over leaves
$y,y'$ gives the exact rational label expansion
$\mathcal Q[A]=\sum_\ell q^A_\ell\,y_\ell$ over four-vertex
multigraph labels of degree $2\deg A$.  The modulated generators are

$$\hat C_n=C_n+\tfrac13T_{|n|}(s),\qquad \hat S_n=S_n+\tfrac13T_{|n|}(s),
\qquad n\in\mathbb Z,$$

$\deg\hat C_n=|n|+1$, $\deg\hat S_n=|n|+5$ ($n\ne0$), with the
machine-checked closed forms and Gram-body bounds

$$|\hat C_n|\le\tfrac43+4|n|,\qquad |\hat S_n|\le\tfrac43+12|n|
\qquad\text{on }|t_1|,|t_2|,|s|\le1 .$$

Write $\hat G^{(2)}_n=\hat C_n$, $\hat G^{(6)}_n=\hat S_n$,
$c_2=4$, $c_6=12$, $a=\tfrac43$.

### 4.2 The diagonal atom, and why not the coherent square

**Definition (theta atom, diagonal/incoherent form).**  For
$f\in\{2,6\}$ and rational $q\in(0,1)$, the theta-atom moment of a
measure $\mu$ is

$$\tau_{f,q}(\mu)\;:=\;\sum_{n\in\mathbb Z}q^{n^2}\,\mathcal Q[\hat G^{(f)}_n](\mu),$$

absolutely convergent because
$0\le\mathcal Q[\hat G^{(f)}_n]\le(a+c_f|n|)^2$ (the leaf integral of
a probability measure is bounded by the sup of the leaf).  On the
moment side each atom adjoins one new scalar pseudo-moment $y_\tau$
— a new label $(\texttt{theta},f,q)$ — linked to the polynomial
labels by one-sided linear cuts (§4.3).  Dually, the cut multipliers
assemble into a certificate generator that may carry *either sign* on
its high-modulation part, the negative part paid for by an exactly
rational constant.  No closed-form value of $\tau_{f,q}$ is ever
needed: the atom exists in the SDP only through its cuts, exactly as
the kernel labels $k_G$ of §1.7.

**Coherent vs diagonal domination.**  The alternative "coherent" atom
— the single leaf $\Theta_q=\sum_nq^{n^2}\hat G_n$ inside one square
$\mathcal Q[\Theta_q]$ — expands with cross terms of both signs; its
truncation error is controlled only by Cauchy–Schwarz windows of
width $\sqrt T$, not rational and not one-sided.  No strength is
lost by going diagonal: by Cauchy–Schwarz on the leaf integrals, at
every measure

$$\mathcal Q\Big[\textstyle\sum_nq^{n^2}\hat G_n\Big]\;\le\;
\Big(\textstyle\sum_nq^{n^2}\Big)\cdot\sum_nq^{n^2}\,\mathcal Q[\hat G_n]
\;=\;\theta_3(q)\,\tau_{f,q}(\mu),$$

so every coherent-square certificate term is dominated by the
diagonal atom up to the explicit rational-majorizable constant
$\theta_3(q)=\sum_nq^{n^2}$.  (The coherent object stays available as
a *finite* leaf $\sum_{|n|\le N}q^{n^2}\hat G_n$ inside the ordinary
two-root blocks — polynomial, hence already in the hierarchy.)

### 4.3 The cuts and their validity

Fix $f,q$; abbreviate $L_n(y)=\sum_\ell q^{\hat G_n}_\ell\,y_\ell$.
The adjunction consists of, for a finite set of truncation orders
$N$:

* **(G$_n$) generator blocks** ($1\times1$ squares):
  $L_n(y)\ge0$ for each carried $|n|\le N_{\max}$;
* **(T$^-_N$) lower truncation cuts**:
  $y_\tau-\sum_{|n|\le N}q^{n^2}L_n(y)\ \ge\ 0$;
* **(T$^+_N$) upper truncation cuts**:
  $\sum_{|n|\le N}q^{n^2}L_n(y)+T^{(f)}_q(N)-y_\tau\ \ge\ 0$,
  with the rational tail majorant $T^{(f)}_q(N)$ of §4.4.

**Validity lemma.**  *For every antipodal probability measure $\mu$,
the assignment $y=y(\mu)$, $y_\tau=\tau_{f,q}(\mu)$ satisfies
(G$_n$), (T$^-_N$), (T$^+_N$) for all $n,N$.  The extended problem is
still a relaxation, its optimum a valid lower bound, and — every
constraint being valid for* all *measures — the adjunction composes
with the reduction lemma and the all-measures certificate track.*
*Proof.*  (G$_n$) is $\mathcal Q[\hat G_n](\mu)\ge0$, a square.
(T$^-_N$): the omitted tail is a sum of squares times positive
$q^{n^2}$.  (T$^+_N$): each omitted term is at most
$q^{n^2}(a+c_f|n|)^2$ and §4.4 majorizes the sum. $\square$

All data are rational: $q^{n^2}\in\mathbb Q$, the expansions
$q^{\hat G_n}_\ell$ are exact rationals (Chebyshev recurrence with
rational seeds), and $T^{(f)}_q(N)\in\mathbb Q$.

### 4.4 The exactly rational tail majorants

For rational $q\in(0,1)$ and $N\ge0$, use the geometric domination
$q^{n^2}\le q^{N^2}\,q^{(2N+1)(n-N)}$ for $n>N$ (because
$n^2-N^2=(n-N)(n+N)\ge(n-N)(2N+1)$ when $n-N\ge1$).  Put
$r=q^{2N+1}\in(0,1)$, $\alpha=a+c_fN$, $m=n-N$, and use
$\sum_{m\ge1}r^m=\frac r{1-r}$,
$\sum_{m\ge1}mr^m=\frac r{(1-r)^2}$,
$\sum_{m\ge1}m^2r^m=\frac{r(1+r)}{(1-r)^3}$:

**First-moment tail** (linear atom variants and diagnostics):

$$\sum_{n>N}q^{n^2}(a+c_fn)\;\le\;
q^{N^2}\Big[\alpha\,\frac{r}{1-r}+c_f\,\frac{r}{(1-r)^2}\Big]
\;=:\;t^{(f)}_q(N)\;\in\mathbb Q .$$

**Second-moment tail** (the one the atom cuts use):

$$\sum_{n>N}q^{n^2}(a+c_fn)^2\;\le\;
q^{N^2}\Big[\alpha^2\,\frac{r}{1-r}
+2\alpha c_f\,\frac{r}{(1-r)^2}
+c_f^2\,\frac{r(1+r)}{(1-r)^3}\Big] ,$$

and the two-sided-in-$n$ majorant used in (T$^+_N$) is

$$T^{(f)}_q(N)\;=\;2\,q^{N^2}\Big[\alpha^2\,\frac{r}{1-r}
+2\alpha c_f\,\frac{r}{(1-r)^2}
+c_f^2\,\frac{r(1+r)}{(1-r)^3}\Big]\;\in\;\mathbb Q ,$$

valid because $|\hat G^{(f)}_{-n}|\le a+c_fn$ as well (the closed
form $C_n=C_0T_n(s)+2t_1(t_2-st_1)U_{n-1}(s)$ extends to $n<0$ with
$T_{-n}=T_n$, $U_{-n}=-U_{n-2}$, so the same triangle bound applies;
machine-checked in `theta_atoms.py --self-test`).  Example values
(`theta_atoms.py --majorants`, exact fractions, checked against
400-term 60-digit partial sums): at $q=\tfrac12$,
$T^{(2)}_q(N)=238,\ 14.4,\ 0.733,\ 9.28\times10^{-3},\
2.72\times10^{-5},\ 1.87\times10^{-8}$ for $N=0..5$ and
$T^{(6)}_q(N)=1860,\ 108,\ 5.76,\ 7.52\times10^{-2},\
2.25\times10^{-4},\ 1.57\times10^{-7}$ — super-exponential collapse
in $N$, so deep-window cuts are numerically tight.

### 4.5 The window inequalities and the recession-cone consequence

$y_\tau$ occurs only in (T$^\pm$).  Eliminating it, the adjunction's
shadow on the polynomial labels is exactly: for every carried pair
$N'<N$,

$$\boxed{\;\sum_{N'<|n|\le N}q^{n^2}\,L_n(y)\;\le\;T^{(f)}_q(N')\;}
\tag{W$_{N',N}$}$$

(take (T$^+_{N'}$) and (T$^-_N$); pairs $N\le N'$ give consequences
of (G$_n$) alone, and any $y$ satisfying all (W$_{N',N}$) with the
(G$_n$) lifts back to a feasible $y_\tau$, e.g.
$y_\tau=\max_N\sum_{|n|\le N}q^{n^2}L_n(y)$).  (W$_{N',N}$) is a
genuinely new valid inequality: a uniform, $q$-graded bound on the
total high-modulation diagonal two-root mass — no finite-degree
polynomial relaxation couples the modulation axis to a constant.

**Recession-cone consequence.**  A recession direction $r$ of the
extended problem must satisfy the homogenized cuts (constants drop):

$$L_n(r)\ \ge\ 0\quad(|n|\le N_{\max}),\qquad
\sum_{N'<|n|\le N}q^{n^2}\,L_n(r)\ \le\ 0
\;\Longrightarrow\;
L_n(r)=0\quad\text{for all }N'<|n|\le N .$$

So the atom kills every escape ray carrying *any* strictly positive
diagonal $\hat C_n/\hat S_n$ mass in the window — and only those.  By
the sign rule the window functional pairs against a ray as
$-\sum_{N'<|n|\le N}q^{n^2}L_n(r)$: negative pairing (= positive
window mass) kills the ray; zero pairing for every $q,f$, window
means the theta ansatz misses the escape.

**Effective window at solver degree $d$.**  The top-degree
cancellation in the closed form works only forward:
$\deg\hat G_n=n+\mathrm{off}_f$ for $n\ge1$ but
$|n|+\mathrm{off}_f+1$ for $n\le-1$, $\mathrm{off}_2=1$,
$\mathrm{off}_6=5$ (verified exactly).  Within the degree-$d$
algebra the atom windows reach
$-(d/2-\mathrm{off}_f-1)\le n\le d/2-\mathrm{off}_f$; at $d=14$:
$n\in[-5,6]$ (family 2), $n\in[-1,2]$ (family 6).  Larger $|n|$ may
be carried (their labels are new but genuine moments) but are
constrained only by (G$_n$) and the cuts, so the load-bearing window
is the within-degree one.

**Certificate-side reading.**  A dual solution assigns multipliers
$\lambda_N\ge0$ to (T$^-_N$) and $\rho_{N'}\ge0$ to (T$^+_{N'}$) with
$\sum\lambda_N=\sum\rho_{N'}$ (the $y_\tau$ coefficient cancels):
certificates gain

$$\sum_{n}q^{n^2}\Big[\underbrace{\textstyle\sum_{N'\ge|n|}\rho_{N'}-\sum_{N\ge|n|}\lambda_N}_{=:c_n,\ \text{may be negative}}\Big]\,\mathcal Q[\hat G_n](\mu)
\;+\;\sum_{N'}\rho_{N'}T^{(f)}_q(N') ,$$

nonnegative for every $\mu$ by construction: a certificate may
*subtract* high-modulation squares at rate $q^{n^2}$, paying an exact
rational constant — the resummation the escaping
$\varepsilon$-certificates were building term by term.

### 4.6 Localized atoms, $h_2$-localized atoms, and sup-caps

**Localized atoms.**  The plain rows $\mathcal Q[\hat G_n]\ge0$ are
implied by the projected two-root blocks; the genuinely new content
is the windows.  The localized escalation uses leaves
$\hat G_n\cdot m_\alpha$ over *all* monomial multipliers
$m_\alpha(t_1,t_2,s)$ within the degree cap (per-multiplier atoms
$\tau_{f,q,m}=\sum_nq^{n^2}\mathcal Q[\hat G_nm]$, same majorants
since $|m_\alpha|\le1$ on the cube).  These squares are **not**
(E1)-admissible, so they reach the (E1)-dead directions the
projection removed.  Multipliers of odd leaf parity are dropped
(their $\mathcal Q$ vanishes identically on antipodal measures).

**$h_2$-localized atoms.**  For $f,q$ as before,
$$\tau^{\mathrm{loc}}_{f,q}(\mu)=\sum_nq^{n^2}\,h_2(\mu)\,
\mathcal Q[\hat G^{(f)}_n](\mu),$$
with label rows the exact $h_2$-shift of $L_n$
($h_2\cdot F=\tfrac32\,(p_2\times\text{labels})-\tfrac12\,\text{labels}$,
`theta_atoms.h2_localize`).  Validity: each term is $h_2\ge0$ times a
valid square, so lower cuts hold; and $h_2\in[0,1]$ on probability
measures, so $h_2\mathcal Q[\hat G_n]\le(a+c_f|n|)^2$ pointwise and
**the same rational tail majorants $T^{(f)}_q(N)$ apply verbatim** to
the upper cuts.  All-measures valid; composes with the reduction
lemma.

**Sup-cap cuts.**  The same sup bound gives, for every single $n$,
the valid rows
$$\mathcal Q[\hat G_n](y)\ \le\ (a+c_f|n|)^2,\qquad
h_2\mathcal Q[\hat G_n](y)\ \le\ (a+c_f|n|)^2$$
(exact rational constants) — the "$N'=n$ slice" of the atom package
carried individually.

### 4.7 Measured verdicts

* **Plain diagonal atoms do not kill the projected ray**: they cut
  the documented least-norm ray (window pairings $-0.010\dots-0.057$,
  every tested $q$, led by $\hat S_0$ with $\ell_0=0.057$) but the
  escape cone is fat — a nearby recession direction with every
  within-degree theta diagonal at exact zero survives at $+1.9\%$
  norm cost, same pure-$g_4$ escape.
* **Localized atoms re-bound the projected dual**: the full localized
  package (146 within-degree rows, $q=\tfrac12$) makes the
  (E1)-projected all-measures `--find-ray` problem **infeasible**,
  and by cone inclusion (the KKT cone only adds PSD blocks and
  equalities, both shrinking the recession cone) the KKT-inclusive
  projected dual is recession-free too.  The gap scalar cuts achieve
  the same with 3 rows, and both repairs leave the projected
  $\varepsilon=0$ bound at $\approx-3.19$: recession repair, not
  certificate strength.
* **Inert at finite-$\varepsilon$ optima and on the pole**: bound
  unchanged on top of the gap cuts ($\Delta\approx2\times10^{-6}$,
  solver noise), selector trace unchanged at $\varepsilon=0.3$
  (multipliers activate, zero trace reduction), and the window/cap
  cuts on the weighted problem inert at GMP (controls 765.0/8133.57
  unchanged).  Since the majorants collapse super-exponentially, the
  inertness is not a constant-size artifact.
* **What the theta ansatz misses** (constructive output): the
  escape's load-bearing content lives in the orientation-odd two-root
  sectors, the $g_4$ harmonic block, and (KKT cone) the
  Hessian/tangent-gap blocks — none touched by $\hat C_n/\hat S_n$.
  A resummed completion of the sharp ansatz would need odd-sector
  modulated families (orientation leaves
  $\det(x_1,x_2,y)B(t_1,t_2,s)$ with their own $q^{n^2}$ towers)
  and/or the kernel family $DQ-\tfrac13(1-s^2)Q(0,0,s)$.
* **Sup-caps** pair kill-signed against the weighted escape (13 of
  28, led by
  $\operatorname{cap}[\text{plain},6,n{=}2]:-5.56\times10^6$,
  $\operatorname{cap}[h_2\text{loc},6,n{=}0]:-5.53\times10^6$; the
  per-generator pairings are sign-alternating with the
  collision-boundary $s=\pm1$, $T_n(\pm1)=(\pm1)^n$ signature) yet
  inert at GMP: codimension-few *scalar* slices reroute — the matrix
  families of §5 are what moved the pole.

Code and artifacts: `theta_atoms.py` (exact generators, expansions,
majorants, pairing CLI, 8 self-tests including label expansion vs
direct integration to $5\times10^{-11}$), `theta_ab.py`,
`theta_ray_blocks.py`, `theta_trace.py`, `theta_export.py`,
`theta_fingerprint_pair.py`; exact GMP exports
`sdpa_runs/deg14_allm_gapcuts_theta.dat-s`,
`sdpa_runs/deg14_h2w_h2all_theta.dat-s` (every cut row verified
exactly at the ONB measure before writing).  Note: `--e1-project`
conjugates block matrices in floating point, so a *projected* export
cannot currently be certification-grade without an exact-conjugation
pass.

---

## 5. The Jensen (averaging-contraction) and fiber-Toeplitz blocks

These are the two constraint families that broke the weighted pole's
inertness — the first families measured *active* at selector optima
and the first to move the pole law.

### 5.1 Setting

Roots $x_1,x_2$ and leaves $y,y'$ sampled independently from $\mu$;
$t_1=x_1\!\cdot\!y$, $t_2=x_2\!\cdot\!y$, $s=x_1\!\cdot\!x_2$.  The
azimuth $\varphi$ of the leaf about the root-pair frame enters
polynomially through the complex frame variable

$$w:=(t_2-st_1)+i\det(x_1,x_2,y),\qquad
|w|^2=(1-t_1^2)(1-s^2),\qquad w=\sin\delta\,\sin\theta\,e^{i\varphi}$$

(frame $x_1=e_3$, $x_2$ at angle $\delta$ from $x_1$, leaf at polar
angle $\theta$).  $\operatorname{Re}w^k$ is a polynomial in
$(t_1,t_2,s)$; $\operatorname{Im}w^k$ carries one factor
$\det(x_1,x_2,y)$.  Same-leaf products reduce to 3-sample
(`triangle`) labels, two-leaf products to 4-sample (`graph_4`)
labels, $h_2$-multiplied copies to disconnected $p_2\times(\cdot)$
product labels.  Two-root sectors: leaf monomials $t_1^it_2^js^k$
split by the parities of $(i{+}j,\,i{+}k,\,j{+}k)$ into `even_00`,
`even_11` (leaf-even) and `odd_01`, `odd_10` (leaf-odd, orientation
leaves $\det\cdot B$); `_minor` carries the root weight $1-s^2$;
`h2loc_` copies are multiplied by $h_2$.

### 5.2 The Jensen (averaging-contraction) block

**Data.**  A vector $\varphi=(\varphi_\alpha)$ of polynomial leaf
functions of the rooted triple — monomials $t_1^it_2^js^k$ in the
leaf-even sectors, or $\det(x_1,x_2,y)\cdot t_1^it_2^js^k$ in the odd
sectors — and a nonnegative root weight $\rho\in\{1,\ 1-s^2\}$.
Define

$$T_{ab}=\mathbb E\big[\rho\,\varphi_a(y)\,\varphi_b(y)\big]
\quad(\text{same leaf; triangle labels}),\qquad
G_{ab}=\mathbb E\big[\rho\,\varphi_a(y)\,\varphi_b(y')\big]
\quad(\text{independent leaves; the existing two-root Gram}).$$

**Theorem (validity, all measures).**

$$\boxed{\;T-G\;=\;\mathbb E_{x_1,x_2}\big[\rho\cdot
\mathrm{Cov}_y(\varphi\mid x_1,x_2)\big]\;\succeq\;0\;}$$

for every probability measure $\mu$ on $S^2$; likewise the
$h_2$-localized copy $h_2\,(T-G)\succeq0$ and the complements
$(1-h_2)\,G\succeq0$, $(1-h_2)(T-G)\succeq0$, since $0\le h_2\le1$
for every measure.  No KKT content; composes with the reduction
lemma.
*Proof.*  $\mathrm{Cov}_y(\varphi)\succeq0$ pointwise in the roots;
$\rho\ge0$; take expectations.  Label expansions are exact: $T$
entries through 3-point monomials, $G$ entries through the standard
4-point two-root reduction, $h_2$-shifts through
$\tfrac32\,p_2\times L-\tfrac12L$. $\square$

Equivalently, with the mixed instant/averaged flag vector
$f=(\varphi(y),\,\mathbb E_y\varphi)$,
$U=\mathbb E[ff^{\top}]=\begin{pmatrix}T&G\\G&G\end{pmatrix}\succeq0$,
and since $G(I-GG^{+})=0$ automatically,
$U\succeq0\iff G\succeq0$ and $T-G\succeq0$: **one new matrix
inequality per parity sector**, the "unfolding"
$[\![A]\!]^2\le[\![A^2]\!]$ of flag calculus.  No family of the
nine-toggle cone imposes it: the two-root Grams constrain $G$ alone,
the three-point `flag_k` Grams constrain different triangle-sector
matrices (one root, two averaged leaves), `--no-pointwise-sos`
removes the same-leaf moment matrix, and no existing block mixes
same-leaf with averaged-leaf entries.

**Unrooted (pair-sample) copy.**  Taking the whole pair $(X,Y)$ as
the sample with flags $(X\!\cdot\!Y)^d$, $d\in\{0,2,4,6\}$:

$$T_{ab}=p_{a+b}\ (\text{pair labels}),\qquad
G_{ab}=p_a\,p_b\ (\text{product labels}),\qquad T-G\succeq0,$$

the covariance matrix of $((X\!\cdot\!Y)^d)_d$ — a Hankel-vs-product
coupling between the pair sector and the pair-product sector not
implied by the `empty_type_flag` Gram.

### 5.3 The spanning theorem, and why $T-G$ is the circle-pair mechanism

**Spanning theorem (per-pair Toeplitz Grams are already inside).**
Let $v=(1,A_1,\dots,A_r)$ be polynomial leaf functions with
$\deg A_i\le7$ whose monomials lie in a single parity sector, and
$M_{ab}=\mathbb E\big[(\mathbb E_yv_a)(\mathbb E_{y'}v_b)\big]$ the
two-root flag Gram.  Then $M=P^{\top}G_{\text{sector}}P$, where $P$
is the coefficient matrix of the $v_a$ over the sector's monomial
basis and $G_{\text{sector}}$ the full `two_root_*` block; every dual
multiplier for $M$ lifts as $P\Lambda P^{\top}$.  The modulated
generators fit: $\hat C_{2k},\hat S_{2k}\in$ `even_00`,
$\hat C_{2k+1},\hat S_{2k+1}\in$ `even_11` (machine-checked for all
within-degree $n$, `toeplitz_blocks.py --self-test`); mixed triples
such as $(1,\hat C_3,\hat C_6)$ split block-diagonally, since
cross-parity entries have odd degree at a root vertex and vanish
identically for antipodal measures.  So every conditional Gram over
modulated leaves with *independent* leaves — including the
$(1,\hat C_k,\hat C_{2k})$ Toeplitz Grams and the whole theta tower
— is a principal restriction of the existing `two_root_even_00`
$\oplus$ `two_root_even_11` blocks: the mechanism cannot be reached
by any two-leaf object at fixed degree.

**The double-angle identity is automatic.**
$(\operatorname{Re}w^k)^2=\tfrac12\big((1-t_1^2)^k(1-s^2)^k
+\operatorname{Re}w^{2k}\big)$ is a polynomial identity applied by
the moment reducer (`graph_expectation_label`) whenever a same-leaf
product is expanded.

**The circle-pair mechanism.**  The load-bearing inequality of the
circle-pair theorem is the fiber Cauchy–Schwarz
$1+\operatorname{Re}\zeta_{2k}=2\,\mathbb E\cos^2k\varphi\ge
2\,(\mathbb E\cos k\varphi)^2$, comparing a **two-leaf** square
(graph$_4$ labels) with a **same-leaf** average (triangle labels).
With the frame sampled as a root pair, the measure-level Toeplitz
matrix tested against latitude functions is exactly a sum of
same-leaf squares ($T$), the mode averages are the averaged flags
($G$), and the proof inequality is a diagonal entry of $T-G$ after
the automatic double-angle reduction.  Grams over modulated leaves
with independent leaves are spanned (above) — the mechanism is
reachable only by mixing same-leaf with averaged-leaf entries, which
is what $T-G$ does.

**Sharpness / saturation.**  $\mathrm{Cov}_y(\varphi)=0$ exactly on
fiber-deterministic configurations (leaf angle a function of the
roots); hence the blocks vanish identically on coherent one-orbit
measures $\{\pm u\}$ (ONB-face compatible; machine-checked exactly,
alongside strict positivity of the $C_1$-direction gap at generic
circle measures).  Dually, a deep escape direction can *saturate*
the blocks by concentrating fibers — the measured degree-14
re-steepening mechanism (**saturation phenomenon**: covariance
kernels = deterministic fibers); the richer degree-16 cone suppresses
it (no re-steepening measured).

### 5.4 The seesaw, the $h_2$-complement family, and the $s^2$ axis

The measured escape moves plain and $h_2$-localized copies of the
*same* generator in opposite directions: the pairing ratio
h2loc/plain lies in $[-6.3,-4.5]$ (mean $\approx-4.9$) for all 14
generators with pairing above $10^3$.  Since the h2loc row is
$\tfrac32(p_2\times L)-\tfrac12L$, a uniform ratio $\approx-4.9$
says the escape's $p_2\times L$ content is $\approx-2.9\times$ its
plain-$L$ content, coherently across the whole modulated tower —
maximally non-measure-like (a genuine measure near ONB has ratio
$+p_2\approx+\tfrac13$).  The plain and h2loc blocks are separate
PSD constraints, so nothing forbids this shuttling; the valid
couplings that do are the complements

$$G_{\text{plain}}-G_{\text{h2loc}}=(1-h_2)\,G\succeq0
\qquad(\text{and }(1-h_2)(T-G)\succeq0),$$

whose label support is exactly the union of the existing plain and
h2loc block labels — zero labels outside the problem.

**The $s^2$ root-weight axis** (motivated by the collision-boundary
signature of the escape): $s^2\cdot\mathrm{Cov}\succeq0$ is valid
since $s^2\ge0$ pointwise, and the weight algebra
$s^2=1-(1-s^2)$ makes $A_{s^2}=A_{\text{plain}}-A_{\text{minor}}$ on
a common basis — $s^2$-blocks relate to plain+minor exactly as the
$h_2$-complement relates to plain+h2loc: a *difference of PSD
constraints*, not implied by them.  (Pure $p_2$- or $(1-p_2)$-weights
are conic combinations of $\{1,h_2\}$-weights and add nothing;
$h_2^2$ would need $p_2^2\times$ labels outside the vocabulary —
both rejected on principle.)

**Cross-pair candidate, settled by data.**  Two disjoint root pairs
either factorize ($\mathbb E[F]\,\mathbb E[G]$: rank one, no
content) or produce graph$_5$/graph$_6$ labels (shared leaf, or
harmonic-index contraction), where the measured escape has $0.0\%$
mass; the sector where the escape lives (graph$_4$ +
$p_2\times$graph$_4$) is exactly the label support of $T-G$ and its
$h_2$-localization — the Jensen block absorbs the cross-pair intent.

### 5.5 The fiber-Toeplitz block

**Data.**  A radial depth $K\ge2$, root polynomials
$g_a(t_1,t_2,s)$, and the frame variable $w$.  Index by pairs
$(j,a)$, $j=0,\dots,K$: $V_{(j,a)}=|w|^{K}e^{2ij\varphi}\,g_a$.
Then $\mathbb E_y[VV^{H}]$ (conditional on the roots, then averaged)
is Hermitian PSD, and its real part has polynomial entries

$$\boxed{\;M_{(j,a),(k,b)}
=\mathbb E\big[(|w|^2)^{\,K-|j-k|}\,
\operatorname{Re}\!\big(w^{2(j-k)}\big)\,g_a\,g_b\big]\;\succeq\;0\;}$$

— the radially-weighted trigonometric moment matrix of the leaf's
azimuthal fiber distribution $(\sin\delta\sin\theta)^{2K}d\mu_y$
pushed to the circle.  Entries are `triangle` labels; $h_2$-localized
copies give $p_2\times$`triangle`.  Only the cosine
($\operatorname{Re}$) part is used — the imaginary parts carry
single $\det$ factors whose 3-point moments are outside the label
algebra — and $\operatorname{Re}$ of a Hermitian PSD matrix is PSD.

**Validity.**  Positivity is *moment-matrix* positivity of a
nonnegative fiber measure — not a sum of squares.  Machine-checked:
exact rational equality of the label expansion against direct
evaluation, and exact PSD at the ONB and 4-point-cross measures
(every exported family, in-run).

**Fejér–Riesz non-spanning theorem ($K\ge2$).**  Every band of $M$
carries the *constant* radial weight $|w|^{2K}$, while a polynomial
Gram $\mathbb E[v_j\bar v_k]$ with $v_j=w^{2j}q_j$ forces band
weights $|w|^{4\min(j,k)}q_jq_k$; matching them requires
$q_j=|w|^{K-2j}$ — non-polynomial for $j>K/2$ (equivalently: the
Fejér–Riesz factors of the truncated Toeplitz cone are the
non-polynomial $|w|^Ke^{2ij\varphi}$).  For $K\le1$ the content
reduces to polynomial squares; for $K\ge2$ the block is strictly
outside the within-degree square cone — and the problem carries no
same-leaf constraints beyond the $T-G$ ties, so the fiber-Toeplitz
blocks add triangle-sector PSD structure nothing else provides.
This is the flag-algebra incarnation of the pointwise $7\times7$
fiber-Toeplitz PSD identified as the honest frontier of the
cylindrical route.  Entry degree of FT$(K,r)$ is $4K+2r$, so degree
16 opens $K=4$.

**Pair-sector companions.**  The localized Hankel
$[\,p_{a+b}-p_{a+b+2}\,]\succeq0$ (moment matrix of $(1-t^2)d\nu$ on
pair moments; not implied by the $g_\ell\ge0$ rows or Grams) and the
weighted pair-Jensen $\mathrm{Cov}\big((1-t^2)t^a\big)\succeq0$,
whose $G$ side is the pair-product matrix
$[(p_a-p_{a+2})(p_b-p_{b+2})]$; $h_2$-localized copies land on
$p_2\times$pair and $p_2\times$pair$\times$pair labels.

### 5.6 Family roster and export soundness

The final family set (43 families at degree 16) comprises: the
Jensen blocks per parity sector with minors
(`jensen_even_00`/`jensen_even_11`/odd-sector and `_minor_d6`
copies, e.g. `jensen_even_00_d8` at deg-16 caps), their h2loc copies
(`h2loc_jensen_*`), the seesaw complements
(`h2comp_gram_even_00/11`, `h2comp_gram_odd_01/10`, minor copies,
`h2comp_cov_even_*_d6`), the $s^2$-weighted Jensen families
(`jensen_*_s2_d6`, h2loc copies), the pair-sample copies
(`jensen_pair`, `h2loc_jensen_pair`, `h2comp_gram_pair`), the
fiber-Toeplitz tower (`ftoep2_even_00_r3`, `ftoep2_even_11_r3`,
`ftoep3_even_00_r0`, and at degree 16 specs $(4,0),(3,2),(2,4)$ per
even sector) with h2loc copies (e.g. `h2loc_ftoep2_even_11_r3`),
and the pair companions (`pair_hankel_loc_d6`,
`pair_jensen_minor_d4` and h2loc copies).  Code:
`toeplitz_blocks.py` (exact families, 6/6 self-tests: exact
expansions = direct covariance at relative gap $<10^{-14}$, exact
PSD at random atomic measures, coherence sharpness, strict $C_1$
gap, $h_2\ge0$, generator spans; h2-shift convention matches
`--h2-localized-all` exactly), `toeplitz_export.py`,
`toeplitz_ab.py`.

**Export soundness.**  The base export drops directions whose images
on the base blocks are linearly dependent (337 at degree 14); a
block appended to the finished `.dat-s` would silently pin those
coordinates — an *invalid strengthening*.  `toeplitz_export.py`
therefore passes the families into `sos_search.export_sdpa_problem`
itself, so the exact elimination / image-selection machinery redoes
the bookkeeping; family entries are integers/halves with the float
hand-off verified round-trip-exact, and every family is re-verified
(exact equality + exact PSD at ONB and cross) inside the run.
Appending unprojected valid families to a projected cone is sound:
they are valid moment constraints for every measure and the
certificate stays inside the all-measures cone.

### 5.7 Verdicts

* **First family active at a selector optimum, and first to move the
  pole law**: the Jensen blocks reduce the minimal trace at every
  tested $\varepsilon$ (gap cuts and theta atoms were slack
  everywhere), and at GMP the pole-decade growth moved
  $10.63\times\to9.59\times$ (Jensen alone), then $\to3.53\times$
  (exponent $\varepsilon^{-1.03}\to\varepsilon^{-0.55}$) with the
  complements/minors/$s^2$ families.
* **The residual rotated to the same-leaf side**
  ($p_2\times$graph$_4$+graph$_4$ collapsed $77\%\to3.8\%$; residual
  $49.7\%$ $p_2\times$triangle): the certificate pumps the $T$-side
  that Jensen ties to $G$, and the fiber-Toeplitz blocks bound that
  side directly, carrying top-6 certificate trace when added.
* **Degree-16 wall**: with the full family set the selector ladder
  has exponent $0.207/0.224$ with **no re-steepening**, and the dual
  wall is $-5.2504\times10^{-9}$ — $250\times$ beyond plain degree
  16 ($-1.31\times10^{-6}$; degree 14 moved
  $-2.79\times10^{-5}\to-3.66\times10^{-6}$) — with the Jensen d7/d8
  towers + `h2loc_two_root_even_00` active at the wall and the wall
  adversary a perturbed zero-family with partial equator mode-4 mass
  ($p_2=\tfrac13+9.7\times10^{-5}$).  Degree 14 re-steepened to
  exponent $1.17$ below $5\times10^{-5}$ as the escape saturated the
  covariance kernels (§5.3); degree 16 does not.
* **Sharp-face repair**: on the weighted-(E1)-projected all-measures
  cone at degree 18 the unstacked selector fails at
  $\varepsilon=10^{-2}$ while the stacked one is optimal — the
  predicted $\varepsilon$-pathology of the sharp-face cone is
  repaired by the admissible-tower cuts.
* These families are **not** members of the sharp (E1) ansatz
  ($\mathrm{Cov}_y\ne0$ at the zero measures): they are
  $\varepsilon$-regime/pole instruments, not sharp-certificate
  generators.

---

## 6. Consolidated validity ledger

Status codes: **AM** = valid for all measures (composes with the
reduction lemma and the all-measures certificate track); **KKT** =
valid only at critical/minimizing measures of the stated functional.

| family / object | statement | validity | proof / check |
|---|---|---|---|
| $h_2E$ weighted target (`--h2-weighted-target`) | certifying $h_2E\ge0$ proves $E\ge0$ | AM (reduction lemma) | §1.2 proof; expansion §1.3 exact rational |
| $h_2$-localized module (`--h2-localized-flags`, `--h2-localized-all`) | $h_2\times(\text{PSD block})\succeq0$, $p_2\times$(equality)$=0$ | AM | $h_2\ge0$; disconnected samples factor; §1.4 |
| $\sqrt{h_2}$ adjunction | dominated by existing hierarchy + localized products | — (negative result) | §1.6 rank-one extension argument |
| kernel labels $k_G$ + truncation cuts (Route C) | $k_G\ge\sum_{\ell\le D}a_\ell\iint P_\ell$ for $a_\ell\ge0$ | AM | positive-definiteness of the kernel tail, §1.7 |
| operator bound $I-A_2\succeq0$ | $\|A_2\|\le1$ | AM | $\pi_2(\rho_x)$ orthogonal |
| gap scalar cuts $e_2,e_3,e_4(B)\ge0$ (`--gap-scalar-cuts`) | linear label inequalities, e.g. $e_2=6+6p_2-8p_4\ge0$ | AM | eigenvalues of $B$ in $[0,2]$; exact expansions `sdpa_runs/gap_invariant_expansions.json` |
| $e_5(B)\ge0$ (`--gap-cut-e5`) | $\det(I-A_2)\ge0$, 35-label expansion | AM | Fact 1; exact checks `sdpa_runs/e5_label_check.py` + 7 tests in `test_sos_search.py` |
| $e_5$ upper cut `gap_cut_e5_upper` | $e_5\le(4/5)^5$ | AM | AM–GM with $\operatorname{tr}B=4$; equality at uniform, test-asserted |
| $e_5$ coverage (`--e5-localized-harmonics`: `e5loc_hankel_even/odd`, `e5loc_harmonic_d`, `e5comp_harmonic_d`) | nonnegative invariant × pair moment matrices / $E[P_d]$ | AM | products of nonnegative quantities over disjoint samples; §3.3 |
| $(h_2+\kappa e_5)E$ target (`--e5-weight`) | certifying it proves $E\ge0$ | AM (reduction lemma, Fact 3) | §3.1; exact-rational assembly, round-trip test |
| $W$-KKT re-encoding ($W=h_2E$) | $\Phi^W_\mu\ge\lambda$ at $W$-critical measures | KKT (for $W$), composition sound | §3.5 compactness argument |
| theta generator rows (G$_n$) | $\mathcal Q[\hat G_n]\ge0$ | AM | squares |
| theta truncation cuts (T$^\pm_N$), windows (W$_{N',N}$) | one-sided rational cuts on $y_\tau$ | AM | §4.3 validity lemma; majorants §4.4, `theta_atoms.py --self-test` (8 checks) |
| $h_2$-localized theta atoms | $h_2\mathcal Q[\hat G_n]$ tower, same majorants | AM | $h_2\in[0,1]$; §4.6 |
| sup-cap cuts | $\mathcal Q[\hat G_n]\le(a+c_f|n|)^2$ (plain and h2loc) | AM | sup bound + Cauchy–Schwarz on the leaf integral |
| Jensen blocks $T-G$ (J1; `jensen_*`, minors, $\rho\in\{1,1-s^2\}$) | $T-G=\mathbb E[\rho\,\mathrm{Cov}_y(\varphi)]\succeq0$, all sectors | AM | §5.2 theorem; exact expansions vs direct covariance (rel. gap $<10^{-14}$), rational equality at ONB and 4-point cross |
| pair-sample Jensen copy (J2; `jensen_pair`) | $[p_{a+b}]-[p_ap_b]\succeq0$ | AM | same theorem, sample = pair |
| $h_2$ bounds (J3) | $0\le h_2\le1$ | AM | $p_2=\operatorname{tr}\Sigma^2\in[1/3,1]$ |
| h2loc / $(1-h_2)$-complement copies (J4; `h2loc_jensen_*`, `h2comp_gram_*`, `h2comp_cov_*`) | scalar $\ge0$ times PSD | AM | J1–J3 |
| exact PSD at sharp measures (J5) | every exported family PSD at ONB and cross | check | rational pivoted-elimination PSD check, in-run |
| $h_2$-shift convention (J6) | matches `--h2-localized-all` exactly | check | gap $0.0$ vs `blocks_deg14_h2w_h2all.json`, incl. minors |
| spanning / parity of $\hat C_n,\hat S_n$ (J7) | per-pair Toeplitz Grams inside `two_root_*` blocks | theorem + check | §5.3; `toeplitz_blocks.py --self-test` |
| coherence sharpness (J8) | even-sector $T-G\equiv0$ at one-orbit $\{\pm u\}$ | check | machine-checked exactly |
| $s^2$-weighted Jensen (`jensen_*_s2_d6`) | $s^2\cdot\mathrm{Cov}\succeq0$ | AM | $s^2\ge0$ pointwise |
| fiber-Toeplitz blocks (`ftoep*`, `h2loc_ftoep*`) | $M_{(j,a),(k,b)}=\mathbb E[(|w|^2)^{K-|j-k|}\operatorname{Re}(w^{2(j-k)})g_ag_b]\succeq0$ | AM | fiber-measure moment positivity (§5.5); exact expansion + exact PSD checks in-run |
| pair companions (`pair_hankel_loc_d6`, `pair_jensen_minor_*` + h2loc) | localized Hankel / weighted pair-Jensen $\succeq0$ | AM | moment matrix of $(1-t^2)d\nu$ / covariance |

KKT-only assumptions in this document: none except the $W$-KKT
re-encoding row, whose composition with the reduction lemma is proved
in §3.5.  Everything else composes with the all-measures certificate
track.
