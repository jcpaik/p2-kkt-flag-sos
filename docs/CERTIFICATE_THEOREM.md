# What an exact certificate would prove

This document separates the exact mathematical implication from the current
floating-point evidence.

## 1. KKT notation

For a probability measure \(\mu\) on \(\mathbb{RP}^2\), write

\[
E(\mu)=\iint K(x\cdot y)\,d\mu(x)d\mu(y),
\]

\[
q_\mu(z)=U_\mu(z)-E(\mu),
\qquad
g_\mu(x)=\nabla_{S^2}U_\mu(x),
\qquad
H_\mu(x)=\operatorname{Hess}_{S^2}U_\mu(x).
\]

At a global minimizer,

\[
q_\mu(z)\ge0\quad\text{for every }z,
\tag{1.1}
\]

and

\[
q_\mu(x)=0,\qquad
g_\mu(x)=0,\qquad
H_\mu(x)\succeq0
\tag{1.2}
\]

for \(\mu\)-almost every \(x\).

In the isotropic branch,

\[
E(\mu)=\frac{16}{3}-48p_4+32p_6.
\tag{1.3}
\]

## 2. Abstract exact certificate

An exact KKT-infused flag certificate is an identity

\[
\begin{aligned}
E(\mu)
={}&
\sum_\alpha
\left\|\Phi_{\alpha,\mu}\right\|_{L^2(\mu^{r_\alpha})}^2
\\
&+
\int
\operatorname{tr}\!\left(B_\mu(x)H_\mu(x)\right)d\mu(x)
\\
&+
\int_{\mathbb{RP}^2}
\rho_\mu(z)q_\mu(z)\,d\sigma(z)
\\
&+
\int a_\mu(x)q_\mu(x)\,d\mu(x)
\\
&+
\int v_\mu(x)\cdot g_\mu(x)\,d\mu(x)
\\
&+
\mathcal Z(\mu).
\end{aligned}
\tag{2.1}
\]

The datum in (2.1) must satisfy:

- every \(\Phi_{\alpha,\mu}\) is a rooted flag evaluation;
- \(B_\mu(x)\succeq0\) pointwise on the tangent plane;
- \(\rho_\mu(z)\ge0\);
- \(a_\mu\) and \(v_\mu\) are arbitrary real equality multipliers;
- \(\mathcal Z(\mu)=0\) follows exactly from antipodal symmetry, isotropy,
  exchangeability, sphere equations, and Gram-rank identities;
- expansion of both sides into moment labels agrees exactly.

At finite degree, the first three types of nonnegative multipliers are encoded
by rational positive semidefinite Gram matrices.

## 3. Certificate theorem for the encoded branch

**Theorem.** Suppose (2.1) is an exact identity for every isotropic probability
measure satisfying the KKT conditions (1.1)–(1.2). Then every such measure has
\(E(\mu)\ge0\).

**Proof.**

Each squared norm in (2.1) is nonnegative. Since \(B_\mu(x)\succeq0\) and
\(H_\mu(x)\succeq0\), their Frobenius contraction is nonnegative:

\[
\operatorname{tr}(B_\mu(x)H_\mu(x))\ge0.
\]

Likewise, \(\rho_\mu(z)q_\mu(z)\ge0\) by (1.1). The two equality-multiplier
integrals vanish by (1.2), and \(\mathcal Z(\mu)=0\) identically. Therefore
the right side of (2.1) is nonnegative, so \(E(\mu)\ge0\).

\(\square\)

## 4. From a KKT certificate to global copositivity

The space of probability measures on the compact space \(\mathbb{RP}^2\) is
weak-* compact. Because \(K\) is continuous, \(E(\mu)\) is weak-* continuous.
Hence \(E\) has a global minimizer.

The following is the precise contradiction argument.

**Global reduction theorem.** Assume:

1. every global minimizer with negative energy lies in the class covered by
   (2.1); and
2. (2.1) has been verified exactly on that class.

Then \(K\) is copositive.

**Proof.**

If copositivity failed, some probability measure would have negative energy.
A global minimizer \(\mu_*\) would then exist and satisfy

\[
E(\mu_*)<0.
\]

Global minimality gives the KKT conditions (1.1)–(1.2). Assumption 1 places
\(\mu_*\) in the class covered by the certificate. The certificate theorem
then gives \(E(\mu_*)\ge0\), a contradiction.

\(\square\)

## 5. Coverage without the isotropy assumption

The current program makes no isotropy assumption. Degree-two sampled
vertices are retained in the canonical moment labels, so the certificate
algebra covers every antipodally symmetric probability measure satisfying
the KKT conditions. Antipodal symmetry itself is exact: \(K\) is even, so
symmetrizing \(\mu\mapsto\tfrac12(\mu+\mu^-)\) preserves \(E\), and a
negative minimizer may be taken antipodal without loss of generality.

Consequently a rationally verified certificate in the current algebra
would prove unrestricted copositivity directly — assumption 1 of the
global reduction theorem holds for the full class.

The isotropy deficit is represented inside the hierarchy by harmonic flag
squares: \(h_2=\sum_m|\hat\mu_{2m}|^2\ge0\) and spin-2 Gram blocks
containing the deviatoric second moment. However, as recorded in
[Numerical results](../RESULTS.md), the numerical dual value of the
non-isotropic hierarchy at degree 14 is approximately
\(-4.5\times10^{-3}\), far from zero: no such certificate has been found,
and the convex relaxation of the spin-2 sector provably leaks at order
\(\sqrt{h_2}\) against the linear cost of \(h_2\). Coverage is therefore
no longer the obstacle; existence of a finite-degree certificate is.

## 6. What “rationalized and exactified” requires

A publishable finite certificate must contain all of the following.

### Exact Gram matrices

Every SOS and positive multiplier block must be a matrix over
\(\mathbb Q\), or an explicitly described algebraic number field, and must be
proved positive semidefinite. Acceptable checks include an exact
\(LDL^T\) decomposition with nonnegative diagonal entries or an exact
factorization \(Q=R^TR\).

### Exact coefficient identity

After expanding every block, the coefficient of every canonical moment label
must equal the target coefficient exactly. A small floating-point residual is
not sufficient.

### Exact equality identities

All free terms must be explicit combinations of identities that truly vanish:

- odd-degree antipodal moments;
- isotropic degree-two contractions;
- permutation/exchangeability relations;
- sphere relations;
- \(\det\operatorname{Gram}_4=0\) and its polynomial multiples;
- the KKT equalities on the support.

### Correctly signed KKT multipliers

Multipliers of \(q_\mu(z)\ge0\) must be nonnegative, and multipliers contracted
with \(H_\mu(x)\succeq0\) must be positive semidefinite tangent matrices.
Only multipliers of equality constraints may be unrestricted.

### Coverage of the minimizing class

Finally, the certificate must apply to every possible negative global
minimizer. For the current code this is exactly where the isotropy caveat
enters.

## 7. Infinite-degree alternative

It may happen that no finite exact identity exists, while for every
\varepsilon>0\) there is a certificate

\[
E(\mu)+\varepsilon=C_\varepsilon(\mu),
\]

where \(C_\varepsilon\) is a valid KKT/SOS expression. Letting
\(\varepsilon\downarrow0\) still proves \(E(\mu)\ge0\), provided the
certificate is valid on the full minimizing class.

Such a family is an asymptotic certificate. Its datum is a convergent family
of positive rooted kernels and matrix-valued Hessian multipliers—the
higher-point analogue of a magic function.
