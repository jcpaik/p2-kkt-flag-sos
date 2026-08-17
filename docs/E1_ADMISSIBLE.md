# The (E1) equations solved: admissible certificate leaves in closed form

Every claim in this document is verified in exact rational arithmetic by

```sh
python3 solve_e1.py            # 54 checks
python3 solve_e1.py --table    # + admissible-dimension tables
```

This carries out the selection principle of
[Exact zero program](EXACT_ZERO_PROGRAM.md) §4.3 — *solve the
complementary-slackness equations in closed form first, and only then let
an SDP choose the remaining coefficients* — instead of extrapolating
solver output.  The master equations are the circle-mode conditions (E1)
of [Limit certificate](LIMIT_CERTIFICATE.md) §4.

## 0. Setup

The zero-energy family is the $SO(3)$-orbit of

$$\mu^*=\tfrac13\,\delta_{\pm e}+\tfrac23\,\nu,\qquad
\nu\ \text{antipodal on}\ C=e^\perp,\quad \hat\nu(2)=\hat\nu(6)=0 .$$

If a sharp certificate exists (at any degree, or over any enriched
algebra), complementary slackness forces every flag-square leaf
$\Phi(R;\cdot)$ to satisfy

$$\int \Phi(R;y)\,d\mu^*(y)=0
\qquad\text{for all } R\in(\operatorname{supp}\mu^*)^r
\text{ and all admissible } (e,\nu). \tag{E1}$$

Two reductions make (E1) finite and exactly computable:

* **Interior lemma.**  The truncated moment sets of admissible $\nu$
  have nonempty interior around Haar measure inside the slice
  $\{\hat\nu(2)=\hat\nu(6)=0\}$, so "$=0$ for every admissible $\nu$" is
  equivalent to "$=0$ identically as a polynomial in the free moments
  $\tau_m=\int e^{im\theta}d\nu$" ($m$ even, $|m|\notin\{0,2,6\}$).
* **Configuration lemma.**  Up to rotation about $e$ and reflection,
  every root configuration is a tuple from $\{\pm e\}\cup C$ with
  symbolic equator angles; rotating a compatible circle about an
  equator root multiplies the spin-$k$ leaf integral by a phase
  $e^{ik\beta}$ and creates no new condition (checked exactly for
  $k\le3$: the "tilted circle" rows change nothing).

Everything below is the exact nullspace of the resulting linear systems,
over $\mathbb Q$, with no floating point.

## 1. One-root leaves: a seven-dimensional space, at every degree

A spin-$k$ leaf about a root $x$ is $\Phi=f(t)\,w^k$ with $t=x\cdot y$
and $w=\sin\alpha\,e^{i\psi}$ (the solver's `flag_k` blocks pair these
through the Chebyshev kernel $R_k$).  Solving (E1):

| spin $k$ | admissible radial $f$ | dim |
|---|---|---:|
| 0 | $\operatorname{span}\{\,T_2+\tfrac13,\;T_6+\tfrac13\,\}$ | 2 |
| 1 | $\operatorname{span}\{\,t,\;U_5\,\}$ | 2 |
| 2 | $\operatorname{span}\{\,1,\;(4t^2-1)^2\,\}$ | 2 |
| 3 | $\operatorname{span}\{\,t^3\,\}$ | 1 |
| $\ge4$ | $\{0\}$ | 0 |

**The dimensions are independent of the degree cap** (stable from cap 8
through 24, and provably for all degrees: the conditions pin the
Chebyshev support).  Structural facts, all verified:

* $K=(T_2+\tfrac13)+(T_6+\tfrac13)$: the kernel itself is the sum of
  the two admissible spin-0 leaves.  Their rooted averages are
  $\int(T_2+\tfrac13)(x\cdot y)\,d\mu(y)=2x^{\mathsf T}(M-\tfrac13I)x$
  (the deviatoric quadratic) and the analogous hexapole average.
* The spin-1 space is exactly the $\theta$-derivative of the spin-0
  space ($\partial_\theta T_2=-2\sin\theta\,\cdot 2t$,
  $\partial_\theta T_6=-6\sin\theta\,U_5$): the admissible gradient
  flags are the gradients of the admissible potentials.
* Spin 3 survives on $t^3$ alone — $\sin^3\theta\cos^3\theta=
  \tfrac{1}{32}(3\sin2\theta-\sin6\theta)$, the $z_1z_2z_3$ structure of
  the ONB potential.  The surviving odd spins $\{1,3\}$ are precisely
  the regular orders $m\in\{1,3\}$ whose pole–equator faces do *not*
  annihilate the target ([Structure](STRUCTURE.md) §5).
* Spins $\ge4$ die entirely; for spin 4 the one candidate
  $(1-t^2)^2f=\tfrac1{16}(8-9T_2+T_6)$, i.e. $f\propto 1+2t^2$, is
  killed by the pole-rooted condition $f(0)=0$ — the $\hat K_4<0$
  obstruction in leaf form.

So the entire one-root layer of any sharp certificate is at most
**seven-dimensional**, at every degree: Gram blocks of sizes
$2,2,2,1$.  Raising the SDP degree adds nothing to this layer — the
degree escalation observed numerically must flow through the two-root
layer and the KKT multipliers.

## 2. Unrooted families

Evaluated on the whole family with symbolic $\nu$:

* $g_2=\iint P_2\,d\mu\,d\mu$ **vanishes identically** on the family
  ($\mu^*$ is isotropic): the $h_2$ face.  $g_4$ and $g_6$ do not:
  `harmonic_4`, `harmonic_6`, ... are dead in any sharp certificate.
* **Pair-power flags** (`empty_type_flag`): the only admissible
  combinations of $\{p_0,p_2,\dots,p_{12}\}$ are
  $$\operatorname{span}\{\,3p_2-1,\;\;E\,\}$$
  — the isotropy deficit and the target itself.  (Dim 2, exactly.)
* **Harmonic-weighted pair flags** (`harmonic_flag_l`): a weight
  combination survives iff its inner average
  $F_c(x)=\int\big(\sum_a c_a(x\cdot y)^a\big)d\mu^*(y)$ kills the
  order-$l$ pairing; the admissible spaces are exactly
  $$\operatorname{span}\{T_2+\tfrac13,\;T_6+\tfrac13\}\ (l\ge4),\qquad
  \operatorname{span}\{1,\;T_2+\tfrac13,\;T_6+\tfrac13\}\ (l=2),$$
  i.e. the (E1) spin-0 radials again — the leaves whose $\mu^*$-average
  vanishes on $\operatorname{supp}\mu^*$.  These blocks add no leaf
  freedom beyond §1.
* **Spin-2 deviatoric flags** (`spin2_flag`): 9 of 10 directions
  survive at degree cap 6 — the deviatoric pairing dies on the isotropic
  family, so this block passes facial reduction nearly whole.
  Consistent with $h_2$ being the leak variable of the
  $\varepsilon$-pole.
* The **uniform global-gap multiplier** is dead by (E3): at the ONB,
  $\int(U-E)\,d\sigma=\tfrac{32}{105}>0$, so its multiplier must vanish
  in any certificate sharp at the ONB (`--check-onb` shows the 0.3048).

## 3. Two-root leaves: where the infinite-dimensional freedom lives

Expository treatment with all generator formulas written out:
[Two-root generators](TWO_ROOT_GENERATORS.md).

A two-root leaf is $\Psi=A(t_1,t_2,s)+\det(x_1,x_2,y)\,B(t_1,t_2,s)$
with $t_i=x_i\cdot y$, $s=x_1\cdot x_2$ (the solver's `two_root_*`
sectors).  The conditions come from five configuration classes: both
roots on the equator at symbolic separation $\delta$ (b2), one root at a
pole (a2 and its mirror), and both roots at poles ($s=\pm1$; dropped for
the `_minor` sectors, whose root weight $1-s^2$ vanishes there).

### The even sector, solved

On the coplanar variety $\{D=0\}$
($D=\det\operatorname{Gram}(x_1,x_2,y)$) parametrize
$t_1=\cos\theta$, $t_2=\cos(\theta-\delta)$, $s=\cos\delta$.  Then:

1. **Profile structure.**  The restriction of an admissible $A$ to the
   circle has $\theta$-modes $\{0,\pm2,\pm6\}$ only, with free
   $\delta$-profiles on modes $2$ and $6$ and the mode-0 profile tied to
   the pole values.  The kernel of the profile map is exactly the family
   $$A\;=\;D\,Q\;-\;\tfrac13(1-s^2)\,Q(0,0,s),$$
   with $Q$ subject to the two pole-root slice conditions.  Verified as
   an exact dimension identity at degree 8:
   $\dim(\text{admissible})=52=26\ (\text{kernel family})+26\ (\text{profile rank})$.
2. **Canonical modulated generators.**  Let $C_n$ be the polynomial
   lift of $\cos(2\theta-n\delta)$ defined by $C_0=T_2(t_1)$,
   $C_1=2t_1t_2-s$, $C_{n+1}=2sC_n-C_{n-1}$, and $S_n$ the lift of
   $\cos(6\theta-n\delta)$ seeded by $S_0=T_6(t_1)$,
   $S_1=2T_4(t_1)C_1-C_{-1}$.  Then
   $$\boxed{\;\hat C_n=C_n+\tfrac13T_n(s)\quad\text{and}\quad
   \hat S_n=S_n+\tfrac13T_n(s)\ \text{ are (E1)-admissible for every }n\;}$$
   with $\hat C_0=(T_2+\tfrac13)(t_1)$ and $\hat S_0=(T_6+\tfrac13)(t_1)$:
   the arity-1 leaves are the $n=0$ members, and every integer
   $\delta$-modulation of the two surviving circle modes is realizable
   (degree grows linearly in $|n|$).
3. **Pure ideal part.**  $D\,Q$ is admissible iff $Q(0,0,s)\equiv0$ and
   both slices $Q(0,t,0)$, $Q(t,0,0)$ lie in
   $\operatorname{span}\{2t^4-t^2\}$ (from
   $(T_2-T_6)/(1-t^2)=8t^2(2t^2-1)$).  In particular
   $D\cdot t_1t_2\cdot(\text{anything of matching parity})$ is
   admissible.
4. **The obstruction, localized.**  Products
   $(T_2+\tfrac13)(t_1)\,(T_2+\tfrac13)(t_2)$ of admissible arity-1
   leaves are *not* admissible: their circle restriction has mode
   $4=2+2$.  This is the Delsarte/two-point obstruction reappearing as a
   mode-selection rule.

### The odd (orientation) sector is nearly free

The coplanar configuration imposes *nothing* on $B$ ($\det=0$ on the
circle, and the $\pm e$ pole terms cancel).  The only conditions are the
two pole-root slices $B(0,\cdot,0),\,B(\cdot,0,0)\in
\operatorname{span}\{U_1,U_5\}$, giving the exact codimension

$$\operatorname{codim}(D)=2\,\#\{m\ \text{even},\ m\notin\{2,6\},\
4\le m\le D+1\}$$

(verified for $D=2,\dots,10$: e.g. 4 of 80 at degree 8).  The heavy use
of the orientation-odd sectors by the numerical optimum is consistent
with this near-freedom.

### Dimension tables

`solve_e1.py --table` prints, per degree, for the combined even/odd
bases and for the four solver parity sectors.  At solver degree 14
(flag degree 7) the projection file records

| sector | admissible / total |
|---|---:|
| `two_root_even_00` | 19 / 30 |
| `two_root_even_11` | 20 / 30 |
| `two_root_odd_01` | 28 / 30 |
| `two_root_odd_10` | 28 / 30 |
| minors | 13/24, 13/20, 19/20, 19/20 |
| `spin2_flag` | 9 / 10 |

## 4. The projected hierarchy: `--e1-project`

```sh
python3 solve_e1.py --export-projection e1_projection_deg14.json --solver-degree 14
python3 sos_search.py ... --degree 14 --e1-project e1_projection_deg14.json
```

restricts every flag-square family to its admissible subspace: one-root
blocks to the closed forms of §1 (spin $\ge4$ dropped), pair flags to
$\{3p_2-1,\,E\}$, `harmonic_flag_l` to the spin-0 radials, and the
two-root sectors and `spin2_flag` to the exact rational bases in the
JSON file.  (The scalar `harmonic_l` inequalities are kept: dead
multipliers may sit at zero, and as valid inequalities they help keep
the $\varepsilon>0$ dual bounded.)  `--e1-project-families` restricts
the projection to a comma list of families for ablation.

`--check-onb` under the projection shows **every projected flag block
vanishes identically at the ONB** (entries at $10^{-16}$), as
complementary slackness demands.  This is complete facial reduction
against the entire continuous zero-energy family, performed symbolically
at the leaf level rather than numerically at the Gram level.

### Measurements (degree 14, MOSEK dual, `--scale-constraints`)

| projected families | dual bound for $E$ |
|---|---:|
| none (baseline, current tree) | $-4.49\times10^{-3}$ |
| one-root only | $-5.39\times10^{-3}$ |
| one-root + pair flags | $-4.97\times10^{-3}$ |
| two-root only | $-2.58\times10^{-2}$ |
| two-root + pair flags | $-2.56\times10^{-2}$ |
| one-root + two-root | $-8.69\times10^{-1}$ |
| everything | **unbounded** (also at degree 16) |

Readings, in order of importance:

1. **The one-root truncation is free.**  Cutting the one-root layer
   from its full degree-14 basis to the seven closed-form leaves moves
   the bound by $\sim2\times10^{-4}$ only: the hierarchy never used the
   dead one-root directions.  This is strong evidence that (E1) is the
   right ansatz layer by layer.
2. **The $\varepsilon>0$ escape runs through dead directions, jointly.**
   Removing the dead two-root directions costs a factor 6; removing
   dead one-root *and* two-root directions together collapses the bound
   and then unbounds the dual.  The recession ray of the fully
   projected dual is a clean, cheap proxy for the escape direction
   $Y_0$ of [Exact zero program](EXACT_ZERO_PROGRAM.md) §1 — it
   identifies exactly which moment directions only dead squares
   control.
3. On the zero family itself, the only certificate terms able to pay a
   positive constant are the global-gap multipliers
   ($\int(U_\mu-E)\,d\sigma = \tfrac{32}{105}$ for every member of the
   family), so any $E\ge\lambda$ certificate with $\lambda<0$ must
   route through them; at sharpness their (E3) slackness kills them.

## 5. What this says about the certificate

* The sharp certificate's square layers are: a 7-dimensional one-root
  layer (fixed for all degrees), the $h_2$/target pair flags, the
  nearly-free spin-2 deviatoric block, and the two-root layer whose
  freedom is **exactly two free profiles** $\gamma_2(\delta)$,
  $\gamma_6(\delta)$ (plus the $(M,Q)$-kernel family and the nearly
  free orientation sector).  All degree escalation is escalation in the
  $\delta$-modulation order $n$ and in the multiplier layers.
* The observed $\exp(-c\,d^2)$ convergence therefore predicts profile
  coefficients $\gamma_{k,n}\sim q^{n^2}$: the natural dictionary atoms
  for the entire-kernel route are **theta-modulated leaves**
  $$\Theta_q^{(2)}=\sum_{n\in\mathbb Z}q^{n^2}\hat C_n,\qquad
  \Theta_q^{(6)}=\sum_{n\in\mathbb Z}q^{n^2}\hat S_n,\qquad q\in\mathbb Q\cap(0,1),$$
  which stay inside the admissible cone term by term; each truncation
  tail is a sum of admissible PSD blocks, so the one-sided cut trick of
  [Exact zero program](EXACT_ZERO_PROGRAM.md) §4.2 applies verbatim.
* The recession direction $Y_0$ should now be re-measured **inside the
  projected cone**: its loading on the $(n$-graded$)$ two-root
  modulations identifies how much of each profile the escape uses, i.e.
  which finitely many dictionary atoms to adjoin.

## 6. The weighted target $h_2E$ (Part D of `solve_e1.py`)

The enriched program of [Exact zero program](EXACT_ZERO_PROGRAM.md)
replaces $E$ by $W=h_2E$.  Its zero set is larger:
$\{W=0\}=\{h_2=0\}\cup\{E=0\}$, and $\{h_2=0\}$ is **every isotropic
measure** (the pole–equator family is itself isotropic, so
$\{E=0\}\subset\{h_2=0\}$ among the known zeros).  Complementary
slackness against Haar measure and its isotropic perturbations then
forces, for every *pure* square of a sharp all-measures certificate:

> the leaf, as a function of the leaf variable $y$, is a **pure
> degree-2 spherical harmonic**, for every root configuration —

while any block carrying an explicit $h_2$ factor is exempt (the factor
already vanishes on the whole zero set).  Intersecting with the
pole–equator conditions of §1 (which turn out to be implied here):

| spin | weighted-admissible radial | vs. unweighted (§1) |
|---|---|---|
| 0 | $\operatorname{span}\{t^2-\tfrac13\}$ | was $\{T_2+\tfrac13,\ T_6+\tfrac13\}$ |
| 1 | $\operatorname{span}\{t\}$ | was $\{t,\ U_5\}$ |
| 2 | $\operatorname{span}\{1\}$ | was $\{1,\ (4t^2-1)^2\}$ |
| $\ge3$ | $\{0\}$ | spin 3 had $\{t^3\}$ |

The three survivors are exactly the spin components of the deviatoric
tensor flag $\int(yy^{\mathsf T}-\tfrac13I)\,d\mu$ rooted anywhere; the
pair-power layer collapses to $\operatorname{span}\{3p_2-1\}$ ($E$ drops
out, since $E$ does not vanish on general isotropic measures).  Hence
**`--h2-localized-all` is the forced structure of the sharp weighted
certificate** — $h_2\times$every family, plus the thin deviatoric
layer — not a heuristic.  All verified exactly (`weighted spin k` and
`weighted pair flags` checks).

**Composition caveat.**  The reduction lemma ($h_2E\ge0$ for all
measures $\Rightarrow E\ge0$) requires an **all-measures** certificate.
Runs that include the KKT-only families (`--gradient --potential
--hessian --global-tangent-gaps`) prove only "any counterexample
minimizer is isotropic", which is vacuous because the zero family is
isotropic.  Proof-carrying weighted runs must drop those toggles (all
remaining families are identities or valid squares), or re-encode the
KKT conditions for the functional $W$ itself.  Measured (MOSEK double,
degree 14, all-measures cone): plain weighted $-1.1\times10^{-2}$; with
`--h2-localized-all` $-3.8\times10^{-4}$ — the proof-carrying cone is
as strong as the KKT-inclusive one.
