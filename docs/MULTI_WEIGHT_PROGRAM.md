# The multi-weight convex program and the $e_5$ discovery

*2026-08-17.  Orchestrator note; exact checks in
`sdpa_runs/e5_face_check.py` (sympy, rational).  Cross-linked from
[PLAN](../PLAN.md) §5.*

## 1. Choosing the multiplier is a convex problem

The enriched-algebra program ([Exact zero program](EXACT_ZERO_PROGRAM.md))
seeks a weight $Q\ge0$ with $QE$ certifiable.  Because the target $E$ is
**fixed**, a multiplier drawn from a finite dictionary
$q=\sum_j\lambda_j w_j$, $\lambda_j\ge0$, enters the certificate identity
*linearly*: each $w_jE$ is a fixed rational label vector, so

$$e\Big(\sum_j\lambda_j\,w_jE\Big)-e(\sigma)=0,\qquad
\sigma\in\text{(certificate cone)},\ \lambda\ge0,\ \textstyle\sum_j\lambda_j=1$$

is jointly convex in $(\lambda,\sigma)$ — one SDP, not an alternation.
The hand-iterated weight search (round 1: $h_2$; round 2: the $e_k$
sign-rule measurements) is the special case of singleton dictionaries.
Requirements on each dictionary element (Exact zero program §5):
(a) $w_j\ge0$ valid for all measures; (b) $w_j>0$ on a dense set
(reduction lemma); (c) $w_j=0$ on the whole zero family (else the pole
persists near the face); (d) is then *selected automatically* by the SDP
instead of being matched by hand.

The bottleneck is the dictionary: until now the only known invariant
satisfying (a)–(c) was $h_2$ itself.

## 2. $e_5(I-A_2)$ satisfies (a)–(c) — with a stronger (c) than expected

Let $A_2=\int\pi_2(\rho_x)\,d\mu$, $\rho_x=2xx^{\mathsf T}-I$, acting on
traceless symmetric $S$ by $A_2[S]=\int\rho_xS\rho_x\,d\mu$, and
$B=I-A_2\succeq0$ (operator bound, [Structure](STRUCTURE.md) §4;
$\pi_2(\rho_x)$ is orthogonal, so $\|A_2\|\le1$).

**Fact 1 (validity, (a)).** $e_5(B)=\det(I-A_2)\ge0$ for every measure.
*(All eigenvalues of $B$ lie in $[0,2]$.)*

**Fact 2 (face-vanishing, (c) — strengthened).** $\det(I-A_2)=0$ for
**every** pole–equator measure $\mu=w\,\delta_{\pm e}+(1-w)\nu$, for
*all* $w$ and *all* circle measures $\nu$ — not only the zero family.
*Proof.* The axial quadrupole $Q=I-3ee^{\mathsf T}$ obeys
$\rho_xQ\rho_x=Q$ for $x=\pm e$ (both diagonal in an adapted basis) and
for every equatorial $x$ (there $\rho_xe=-e$, hence
$\rho_x(I-3ee^{\mathsf T})\rho_x=I-3ee^{\mathsf T}$).  So $A_2[Q]=Q$ for
any measure supported on $\{\pm e\}\cup e^\perp$, i.e. $1\in\operatorname{spec}A_2$
and $\det B=0$. $\square$
Verified symbolically for the parametrized family
($w,\ \hat\nu(2),\ \hat\nu(4)$ free): `e5_face_check.py` V1/V4 give
$\det(I-A_2)\equiv0$ identically in all three parameters, alongside
exact reproduction of the documented spectra
(ONB: $\{1,1,-\tfrac13,-\tfrac13,-\tfrac13\}$; $e_{1..5}(B)$ at ONB
$=4,\ \tfrac{16}3,\ \tfrac{64}{27},\ 0,\ 0$; at pole+Haar
$=4,\ \tfrac{52}9,\ \tfrac{32}9,\ \tfrac{64}{81},\ 0$).

**Fact 3 (dense positivity, (b)).** $\{e_5>0\}$ is weak-\* dense.
*Proof.* $A_2$ is affine in $\mu$ and $A_2(\text{uniform})=\tfrac15I$
(the trace is $\chi_2(\rho_x)\equiv1$, and irreducibility under the
rotation average forces the scalar).  For any $\mu$ set
$\mu_t=(1-t)\mu+t\,\text{uniform}$: then
$B(\mu_t)=(1-t)B(\mu)+\tfrac{4t}5I\succeq\tfrac{4t}5I\succ0$, so
$e_5(\mu_t)>0$ for all $t>0$ and $\mu_t\to\mu$. $\square$
(The same argument re-proves the reduction lemma for any weight of the
form $q=h_2+\kappa\,e_5(B)$, $\kappa\ge0$.)

**Why $e_5$ and not $e_2,e_3,e_4$.**  On the continuous strata of the
zero family $B$ has a one-dimensional kernel, so $e_5$ is the *first*
elementary symmetric invariant forced to vanish there; $e_2,e_3,e_4$ are
measured nonzero on the face (Exact zero program §5) and indeed carry
the *cut* sign, not the weight sign.  This explains structurally why the
measured depth-1 operator-gap blocks and the scalar $e_k$ cuts could not
serve as weights: **the first face-vanishing scalar shadow of the
operator gap lives at word length five.**

## 3. The combined weight $q_\kappa=h_2+\kappa\,e_5(B)$

Zero set $\{q_\kappa=0\}=\{h_2=0\}\cap\{e_5=0\}$ = isotropic measures
with an $A_2$-invariant quadrupole — *strictly smaller* than
$\{h_2=0\}$.  Consequences:

* The weighted-(E1) forced-vanishing conditions for $q_\kappa E$ are
  *weaker* than for $h_2E$ (certificates get more room), while the
  reduction lemma still applies (Fact 3 argument).
* Whether $e_5$ absorbs the surviving escape is exactly requirement
  (d): the pairing of $e_5$ against the unprojected-escape fingerprint
  (PLAN Next actions #1).  This cannot be evaluated on the stored ray
  directly (see §4); the operational test is the localized-module solve
  below.

## 4. The arity obstruction, and the implementable shadow

$e_5$ expands via $\operatorname{tr}A_2^k=E[\chi_2(\rho_{x_1}\cdots\rho_{x_k})]$,
$k\le5$.  The exported label algebra is four-sample (`graph_4`
multigraph labels; inspected in `deg14_h2w_h2all.dat-s.map.json`), so
the genuinely five-sample cycle labels of $\operatorname{tr}A_2^5$ are
outside it — the obstruction is **label arity**, exactly as PLAN §4
noted ("needs five-sample labels outside the four-point algebra"), not
polynomial degree: the per-vertex degree of every $e_5$ label is $\le4$,
well inside the degree budget.  A `graph_5` extension of the moment
reducer restricted to vertex degree $\le4$ is therefore a *bounded*
piece of work, and the weight $e_5E$ (6 samples counting the $E$
factor's pair — disconnected labels) needs `graph_5`×pair products,
i.e. the same disconnected-label mechanism that already handles
$p_2p_j$.  The ray pairing $\langle e_5,\text{ray}\rangle$ is undefined
until then.

The implementable counterpart: **$A_2$-word operator-gap blocks of word
length up to five** (PLAN §4 already names "A₂-word blocks" as the
escalation), i.e. localized blocks
$G_{\alpha\beta}=E[v_\alpha^{\mathsf T}(I-A_2)v_\beta]$ over features
$v_\alpha$ that are words $A_2^{d}[\text{spin-2 flag}]$, $d\le4$.  The
$5$-word diagonal contains the $e_5$-graded directions (Cayley–Hamilton
ties $\det B$ to $B$-words of length $\le5$ against the quadrupole
eigendirection).  Design constraint from Fact 2: the blocks are
sharpness-compatible on the whole pole–equator stratum because $B$
annihilates the adapted quadrupole there — the same mechanism that made
the depth-1 gap blocks ONB-audit-clean.

## 5. The $W$-KKT re-encoding (PLAN §5 route (ii)), derived

For $W=h_2E$ the first variation at $\mu$ in direction $\nu-\mu$ has
density

$$\Phi^W_\mu(x)\;=\;3\,(x^{\mathsf T}M_2(\mu)\,x)\,E(\mu)\;+\;2\,h_2(\mu)\,U_\mu(x),
\qquad U_\mu(x)=\int K(x\cdot y)\,d\mu(y),$$

since $h_2=(3p_2-1)/2$ gives $h_2'(\mu)(x)=3\,x^{\mathsf T}M_2x$.  The
KKT/first-variation inequality $\Phi^W_\mu\ge\lambda$ (equality on
$\operatorname{supp}\mu$) is expressible in the existing label algebra:
$(x^{\mathsf T}M_2x)\cdot E$ terms are disconnected products of pair
labels with rooted spin-2 flags (root + 3 samples), and $h_2\cdot U$
terms are $h_2$-localized potential flags (root + 3 samples) — the same
disconnected-label mechanism as `--h2-localized-all`.  Composition is
sound: a certificate over the $W$-KKT cone proves $W\ge0$ at every
$W$-critical measure; $W$ attains its minimum (weak-\* compactness,
polynomial continuity), minimizers are critical, hence $W\ge0$
everywhere, hence $E\ge0$ by the reduction lemma.  This upgrades the
strong KKT-inclusive numbers to proof-carrying status once the four
KKT toggles are re-derived for $W$ (gradient, potential, Hessian,
tangent gaps with the $\Phi^W$ density above).  Implementation queued
behind the fingerprint verdict.

## 6. Decision protocol (feeds PLAN Next actions #1)

1. Fingerprint the unprojected escape (selector `yMat` norms across
   $\varepsilon$) — in progress.
2. Expand the dominant escape direction in labels; evaluate its
   *geometric* content by fitting model measure families; test whether
   $e_5$ (equivalently the quadrupole eigenvalue $\lambda_Q(B)$) grows
   to first order along those families.  ($e_5$ is flat along the
   equator mode-4 direction $\hat\nu(4)$ — `e5_face_check.py` V4 — so
   the test discriminates genuinely.)
3. If yes: implement depth-$\le4$ $A_2$-word gap blocks (the $e_5$
   shadow) and/or the label-algebra extension for the true
   $q_\kappa=h_2+\kappa e_5$ weighted target, and run the multi-weight
   convex program of §1 with dictionary $\{h_2,\ e_5\text{-shadow}\}$.
4. If no: the fingerprint's own sign-rule verdict names the next
   object; the convex program of §1 still applies to whatever
   dictionary emerges.
