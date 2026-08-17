# Wrapper lemmas for the certificate routes

These are the statements that turn an exact SDP certificate into a proof
of the conjecture.  They are standard, but every route in
[PLAN §4](../PLAN.md) relies on them, so they are collected here with
proofs or proof sketches.  Notation: $\mathcal M$ is the set of Borel
probability measures on $S^2$, $\mathcal M_a\subset\mathcal M$ the
antipodally symmetric ones,
$E(\mu)=\iint K(x\cdot y)\,d\mu\,d\mu$,
$U_\mu(z)=\int K(z\cdot y)\,d\mu(y)$.

## L0. Reduction to antipodal measures

$K$ is even, so for the symmetrization
$\bar\mu=\tfrac12(\mu+\mu_-)$ (with $\mu_-$ the pushforward under
$x\mapsto-x$) all four sign-combinations in the double integral are
equal and $E(\bar\mu)=E(\mu)$.  Hence
$\inf_{\mathcal M}E=\inf_{\mathcal M_a}E$ and it suffices to prove
$E\ge0$ on $\mathcal M_a$.

## L1. Existence of a minimizer

$\mathcal M_a$ is weak-\* compact (Banach–Alaoglu; antipodality and
total mass are weak-\* closed constraints) and $E$ is weak-\*
continuous ($K$ is a polynomial, so $E$ is a finite sum of products
of moments).  Hence $E$ attains its minimum at some
$\mu^\dagger\in\mathcal M_a$.

## L2. First-order conditions at a minimizer

For $\nu\in\mathcal M_a$ and $t\in[0,1]$,
$$E\big((1-t)\mu^\dagger+t\nu\big)
 =E(\mu^\dagger)+2t\Big(\textstyle\int U_{\mu^\dagger}\,d\nu-E(\mu^\dagger)\Big)+O(t^2),$$
so minimality gives $\int U_{\mu^\dagger}\,d\nu\ge E(\mu^\dagger)$
for every $\nu$; taking $\nu=\tfrac12(\delta_z+\delta_{-z})$:

$$U_{\mu^\dagger}(z)\ \ge\ E(\mu^\dagger)\qquad\text{for every }z\in S^2 .$$

Moreover $\int(U_{\mu^\dagger}-E)\,d\mu^\dagger=0$ with a continuous
nonnegative integrand, so $U_{\mu^\dagger}=E(\mu^\dagger)$ on
$\operatorname{supp}\mu^\dagger$.

## L3. Second-order conditions at support points

By L2 every $x_0\in\operatorname{supp}\mu^\dagger$ is a global minimum
of the smooth function $U_{\mu^\dagger}$ on $S^2$, hence
$\nabla_{S^2}U_{\mu^\dagger}(x_0)=0$ and
$\operatorname{Hess}_{S^2}U_{\mu^\dagger}(x_0)\succeq0$.

## L4. Identities valid for every measure

The following hold identically and may be multiplied by *free-sign*
multipliers: $|x|^2=1$ on the sphere; vanishing of odd antipodal
moments; exchangeability of sample points; and
$\det\operatorname{Gram}(x_1,x_2,x_3,x_4)=0$ for any four points of
$\mathbb R^3$ (the `--rank-relations` family).  The identity
$\int U_\mu\,d\mu=E(\mu)$ is the constant-multiplier case of the
potential relation; nonconstant potential/gradient/Hessian multipliers
and the global-gap terms are **not** identities — they are valid only
through L2/L3 at a minimizer.

## L5. Validity of a certificate

Suppose rational data give the label-algebra identity

$$E-\lambda=\sum_\beta\langle Q_\beta,\,G_\beta(\mu)\rangle
+\sum_j c_j\,R_j(\mu)+\sum_k f_k\,I_k(\mu)$$

with $Q_\beta\succeq0$ (so each pairing $\ge0$ at every measure),
$c_j\ge0$ multiplying quantities $R_j\ge0$ **at minimizers** (the
L2/L3 families), and $f_k$ free multiplying the L4 identities
$I_k\equiv0$.  Evaluating at $\mu^\dagger$ gives
$E(\mu^\dagger)\ge\lambda$, hence $E\ge\lambda$ on all of
$\mathcal M_a$ by L1.  If no $R_j$ terms appear, the same display
evaluated at an arbitrary $\mu$ gives $E\ge\lambda$ directly, with no
appeal to L1–L3 (an *all-measures* certificate).

## L6. The weighted target and its composition requirement

With $h_2=\tfrac{3p_2-1}2\ge0$ (zero exactly on isotropic measures),
the reduction lemma of [Exact zero program](EXACT_ZERO_PROGRAM.md) §2.1
states: if $h_2E\ge0$ **for every** $\mu\in\mathcal M_a$, then
$E\ge0$ on $\mathcal M_a$ (divide on $\{h_2>0\}$; on $\{h_2=0\}$
perturb by $\delta_{\pm e}$ and use weak-\* continuity).

*Composition caveat.*  The hypothesis must hold for all measures.  A
certificate of $h_2E\ge0$ that uses the L2/L3 families proves the
inequality only at minimizers of $E$; combined with L1 it yields
"$h_2(\mu^\dagger)E(\mu^\dagger)\ge0$", and since a hypothetical
counterexample minimizer could be isotropic ($h_2=0$) — indeed the
conjectured zero family is isotropic — this is not a contradiction.
Valid conclusions require either (i) an all-measures weighted
certificate (L5 with no $R_j$ terms), or (ii) L2/L3 re-derived for the
functional $W=h_2E$ itself.  By the product rule its first-variation
potential is
$$\Phi^W_\mu(z)=h_2(\mu)\,U_\mu(z)+E(\mu)\,u_\mu(z)+c(\mu),\qquad
u_\mu(z)=\tfrac32\int(z\cdot x)^2\,d\mu(x),$$
with the $\mu$-dependent constant
$c(\mu)$ fixed by $\int\Phi^W_\mu\,d\mu=2W(\mu)$; every term is
polynomial in the label algebra.  With this potential, L1–L3 and L5
apply to $W$ verbatim and give $W\ge\lambda$ on all of $\mathcal M_a$.
(Not yet implemented in the solver.)

## L7. Theta atoms: convergence and truncation cuts

Let $\hat C_n=C_n+\tfrac13T_n(s)$, $\hat S_n=S_n+\tfrac13T_n(s)$ be the
modulated two-root generators
([Two-root generators](TWO_ROOT_GENERATORS.md)).  The recurrence
solves in closed form (checked exactly in `solve_e1.py`):

$$C_n=C_0\,T_n(s)+2t_1(t_2-st_1)\,U_{n-1}(s),\qquad
S_n=S_0\,T_n(s)+(S_1-sS_0)\,U_{n-1}(s).$$

On the Gram body ($|t_1|,|t_2|,|s|\le1$) we have $|T_n|\le1$,
$|U_{n-1}|\le n$, $|C_0|,|S_0|\le1$, $|2t_1(t_2-st_1)|\le4$,
$|S_1-sS_0|\le12$, so

$$|\hat C_n|\le\tfrac43+4n,\qquad |\hat S_n|\le\tfrac43+12n .$$

Hence for rational $q\in(0,1)$ the theta atoms
$\Theta_q=\sum_{n\in\mathbb Z}q^{n^2}\hat C_n$ (and the $\hat S$
analogue) converge absolutely and uniformly, and every Gram pairing of
$\Theta_q$-leaves is the absolutely convergent sum of the pairings of
its terms.  Since each $\hat C_n$-block pairing is PSD at every
measure, the tails are PSD and the one-sided cuts

$$\operatorname{Pair}(\Theta_q)\;\succeq\;\sum_{|n|\le N}q^{n^2}\operatorname{Pair}(\hat C_n)
\qquad(N=0,1,2,\dots)$$

are valid at every measure, with exactly rational data ($q^{n^2}$ and
the $\hat C_n$ coefficients are rational).  The adjoined atom labels
never need closed-form values: they enter the SDP only as variables
bounded by these cuts.

## Status

| lemma | status |
|---|---|
| L0–L5 | written above; routine, to be included in the final paper |
| L6(i) all-measures cone | measured strong (see RESULTS); the proof-carrying track |
| L6(ii) $W$-KKT encoding | not yet implemented |
| L7 | bounds proved above; closed form machine-checked (`solve_e1.py`) |
