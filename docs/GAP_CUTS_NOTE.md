# From failed multipliers to cutting planes: the $e_k(I-A_2)$ session

A mathematical account of the session that turned the spectral-gap
invariants of $A_2$ from rejected second-round multipliers into the
cutting planes that kill the projected escape ray.  Companion documents:
[Exact zero program](EXACT_ZERO_PROGRAM.md) (the program this belongs
to), [H2 weighted experiment](H2_WEIGHTED_EXPERIMENT.md) (round one),
[PLAN](../PLAN.md) §4 (status).

## Setup

$\mu$ ranges over antipodally symmetric probability measures on $S^2$.
With $p_j=\iint (x\cdot y)^j\,d\mu\,d\mu$, the conjecture is
$E(\mu)\ge 0$ where
$E=-\tfrac43+20p_2-48p_4+32p_6$.

The solver works with a **moment relaxation**.  A *label* is an
isomorphism class of an edge-weighted graph on $k$ vertices; the label
with weights $w_{ij}$ stands for the moment
$\int\!\cdots\!\int\prod_{i<j}(x_i\cdot x_j)^{w_{ij}}\,d\mu^{\otimes k}$.
The relaxation replaces $\mu$ by a vector $y=(y_\ell)_\ell$ of real
numbers indexed by labels ("pseudo-moments"), constrained by the
linear-algebraic consequences of $\mu$ being a measure that we can
afford to enforce: for each family of test functions (flag squares,
localized blocks, …) a matrix $A_b(y)$, linear in $y$, that must be
PSD; plus linear equalities (KKT identities, rank identities).  Every
real measure gives a feasible $y$, so

$$\text{bound}\;=\;\min\{\langle E,y\rangle : y\ \text{feasible},\
y_{\mathrm{const}}=1\}\;\le\;\inf_\mu E(\mu),$$

and by SDP duality the bound is certified by a **certificate**: PSD
matrices $Q_b$ with
$E-\text{bound}=\sum_b\langle Q_b,A_b(\cdot)\rangle+(\text{equality
multipliers})$ as an identity in the label algebra.

Prior state.  $\min\langle E,y\rangle$ at degree 14 is
$-4.5\times10^{-3}$; weighting the target by
$h_2=\tfrac{3p_2-1}{2}=\sum_m|\hat\mu_{2m}|^2$ and adding the
$h_2$-multiplied blocks gave
$\min\langle h_2E,y\rangle=-2.8\times10^{-5}$ ($h_2E$ is again linear
in labels because a product of moments over disjoint sample sets *is* a
label — a disconnected graph).  But the minimal-trace certificate of
$h_2E+\varepsilon\ge0$ still satisfies
$\operatorname{tr}C(\varepsilon)\to\infty$ as $\varepsilon\downarrow0$:
non-attainment persists.

## The escape ray

Non-attainment of the dual is equivalent to unboundedness of a
*projected* primal, witnessed by a **recession direction**: a vector
$r$ with

$$r_{\mathrm{const}}=0,\qquad A_b(r)\succeq0\ \ \forall b,\qquad
\text{equalities}(r)=0,\qquad \langle E,r\rangle=-1 .$$

If such $r$ exists, $y+tr$ is feasible for all $t\ge0$ and
$\langle E,y+tr\rangle\to-\infty$: the pseudo-moments escape to
infinity while the target decreases.  `--find-ray` computes the
least-norm such $r$ on the (E1-)projected problem.  The documented ray (stored in the legacy $(3/16)E$ normalization, in
which its pairings are quoted)
pairs as $\langle\iint P_4,\,r\rangle=5.35$ with all other harmonic
energies pairing to $0$ — the escape is concentrated in the $\ell=4$
mode, the unique negative Legendre coefficient of $K$.

## What was done

### 1. Candidate second multipliers from the operator bound

Embed $\mathbb{RP}^2\hookrightarrow SO(3)$ by
$x\mapsto\rho_x=2xx^{\mathsf T}-I$ (rotation by $\pi$ about $x$), and
let $\pi_2$ be the 5-dimensional spin-2 representation.  Then
$A_2=\int\pi_2(\rho_x)\,d\mu$ is a symmetric contraction, so
$I-A_2\succeq0$ is a *valid* constraint for every measure; and at every
zero-energy measure $A_2$ has eigenvalue exactly $1$ (eigenvector: the
axial quadrupole $\mathrm{diag}(1,1,-2)$).  The natural scalar
invariants are the elementary symmetric polynomials
$e_k(I-A_2)\ge0$ of its eigenvalues.

### 2. Exact label expansions

Using the character identity $\chi_2(R)=\tau^2-\tau-1$ with
$\tau=\operatorname{tr}R$, one gets
$\operatorname{tr}(A_2^k)=E\big[\chi_2(\rho_{x_1}\cdots\rho_{x_k})\big]$,
and $\operatorname{tr}(\rho_{x_1}\cdots\rho_{x_k})$ expands into cyclic
dot-product monomials
$\sum_S 2^{|S|}(-1)^{k-|S|}\prod_{\text{cycle of }S}(x_i\cdot x_j)$.
Hence every $\operatorname{tr}(A_2^k)$ is an exact rational combination
of $k$-vertex labels, and Newton's identities give the $e_k$.
(Byproduct: $\operatorname{tr}(I-A_2)=4$ identically, since
$\chi_2(\rho_x)=1$.)  For example

$$e_2(I-A_2)=6+6p_2-8p_4 .$$

Full expansions: `sdpa_runs/gap_invariant_expansions.json`.

### 3. The admissibility tests for a second multiplier

A second multiplier $Q$ in an iterated rational certificate
$Q\,h_2\,E=\text{SOS}$ must satisfy:

* (i) $Q\ge0$ valid for every measure;
* (ii) $Q>0$ on a dense set of measures (reduction lemma);
* (iii) $Q=0$ on the zero-energy family — otherwise near a minimizer
  $\mu^*$ with $Q(\mu^*)>0$ the product vanishes to the same order as
  $E$ and inherits the same representability obstruction;
* (iv) $\langle\mathrm{lin}(Q),\,r\rangle>0$ — $Q$ must grow to first
  order along the escape so one factor absorbs a simple pole (this is
  what $h_2$ did to its own leak).

Measured (exact face values; ray pairings calibrated by
$\langle g_2,r\rangle=0.005$ and $\langle\iint P_4,r\rangle=5.353$,
both reproducing the documented values):

| | ONB | pole–equator strata | $\langle\cdot,r\rangle$ |
|---|---:|---:|---:|
| $e_2$ | $16/3$ | $52/9$ | $-9.79$ |
| $e_3$ | $64/27$ | $32/9$ | $-19.26$ |
| $e_4$ | $0$ | $64/81$ | $-5.85$ |

Every candidate fails (iii) — on the continuous strata $I-A_2$ has a
one-dimensional kernel, so only $e_5=\det(I-A_2)$ vanishes identically
on the face, and that needs 5-vertex labels outside the arity-4
algebra — and every candidate fails (iv), with strictly *negative*
pairings.

### 4. The sign flip turns them into cuts

Each $e_k(I-A_2)\ge0$, expanded in labels, is a **linear** inequality
in $y$ (again: products of moments over disjoint samples are
disconnected labels).  Appending it to the relaxation as a $1\times1$
block adds the constraint $\langle e_k,r\rangle\ge0$ to any recession
direction.  Since the measured pairing is $-9.79<0$, the documented ray
violates the cut — and, simultaneously, the negative pairing is a
*separation certificate*: the inequality cannot be implied by the
existing constraint set, or the ray could not have been feasible.
Concretely, $e_2\ge0$ is

$$p_4\ \le\ \tfrac34\,(1+p_2).$$

Implemented as `--gap-scalar-cuts` (three $1\times1$ blocks:
$e_2,e_3,e_4$).

### 5. A/B verification

Re-solving the least-norm-ray problem on the projected relaxation:

* **without** the cuts — feasible, $\|r\|^2=55.7$,
  $\langle\iint P_4,r\rangle=5.3498$: the documented ray, reproduced
  (`sdpa_runs/ray_nocuts_deg14.json`);
* **with** the three cuts — **infeasible**
  (`sdpa_runs/ray_gapcuts_deg14.json`).

So not only that ray but *every* recession direction is excluded: the
projected relaxation is bounded again.  (The depth-1 operator-gap
localizing blocks `--spin2-operator-gap` had not achieved this.)

### 6. In progress

$\min\langle h_2E,y\rangle$ over the cut relaxation at degree 14 in
high precision (`sdpa_runs/deg14_h2w_h2all_cuts.dat-s`), then the
$\varepsilon$-trace law $\operatorname{tr}C(\varepsilon)$ on it:
bounded trace as $\varepsilon\downarrow0$ means the dual optimum is
attained, which is the precondition for rounding the certificate to
exact rationals and verifying it independently — the endgame for the
exact zero.

## Summary

The invariants of $I-A_2$ fail as multipliers because their ray
pairings are negative — and that same negativity is precisely what
makes them valid, non-redundant cutting planes that render the escape
infeasible.
