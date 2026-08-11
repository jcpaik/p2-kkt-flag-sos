# Collective transport second variation

This note adds a stronger positional optimality condition for a hypothetical
minimizer of the P2 energy.  The existing rooted Hessian constraints move one
support point against a fixed background measure.  A true minimizer is also
stable when **all support points move simultaneously** according to a tangent
vector field.  The resulting second variation contains off-diagonal terms
coupling the motions at two sampled points.

## 1. Transport variation

Let

\[
E(\mu)=\iint K(x\cdot y)\,d\mu(x)d\mu(y),
\qquad
K(t)=32t^6-48t^4+20t^2-\frac43.
\]

Fix a tangent field \(V(x)\perp x\).  Move each point along the spherical
geodesic with initial velocity \(V(x)\):

\[
x_t=\cos(t\|V(x)\|)x+
\frac{\sin(t\|V(x)\|)}{\|V(x)\|}V(x),
\qquad
\mu_t=(x\mapsto x_t)_*\mu.
\]

At a local minimizer,

\[
E'(0)=0,\qquad E''(0)\ge0.
\]

For \(s=x\cdot y\), \(v=V(x)\), and \(w=V(y)\),

\[
\frac12E''(0)=\iint\Big[
K''(s)(v\cdot y)^2-sK'(s)\|v\|^2
+K''(s)(v\cdot y)(x\cdot w)
+K'(s)v\cdot w
\Big]d\mu(x)d\mu(y).
\]

The first two terms are the averaged pointwise support Hessian already used by
the current hierarchy.  The last two are the new **collective cross term**.
They measure whether two support points can lower the energy by moving in a
correlated way.

## 2. Polynomial vector fields

The initial implementation uses the projectively well-defined fields

\[
V_r(x)=\mathbb E_Z[(x\cdot Z)^r P_{x^\perp}Z],
\qquad r=1,3,5,\ldots,
\]

where

\[
P_{x^\perp}Z=Z-(x\cdot Z)x.
\]

Odd \(r\) is forced by antipodal symmetry: the lifted field satisfies
\(V_r(-x)=-V_r(x)\), as required for a tangent vector field on
\(\mathbb{RP}^2\).

For

\[
V=\sum_r c_rV_r,
\]

the second variation is a quadratic form in the coefficients \(c_r\).
Therefore its moment matrix must be positive semidefinite.

## 3. Four-point Gram formula

Use independent samples \(X,Y,Z,W\) and the edge order

\[
XY=a,\quad XZ=c,\quad XW=d,\quad
YZ=b,\quad YW=e,\quad ZW=f.
\]

For the local term, set

\[
v=P_{X^\perp}Z,\qquad u=P_{X^\perp}W.
\]

Then

\[
L=
K''(a)(b-ac)(e-ad)
-aK'(a)(f-cd).
\]

For the collective cross term, set

\[
v=P_{X^\perp}Z,\qquad w=P_{Y^\perp}W.
\]

The relevant contractions are

\[
v\cdot Y=b-ac,
\qquad
X\cdot w=d-ae,
\]

and

\[
v\cdot w=f-eb-cd+ace.
\]

Hence

\[
C=
K''(a)(b-ac)(d-ae)
+K'(a)(f-eb-cd+ace).
\]

For basis powers \(r,s\), the matrix entry is the expectation of

\[
(X\cdot Z)^r(X\cdot W)^s L
+
(X\cdot Z)^r(Y\cdot W)^s C.
\]

This is a polynomial in the six Gram entries, so it passes through the same
antipodal, isotropic-contraction, exchangeability, and rank-three reduction
machinery as the existing flag blocks.

## 4. Degree bookkeeping

Both base kernels have total Gram degree at most 8.  Therefore, with an SDP
degree cap \(D\), a square matrix indexed by powers \(r,s\) is admissible when

\[
r+s+8\le D.
\]

The default symmetric choice is

\[
r,s\in\{1,3,5,\ldots\},
\qquad
r,s\le\left\lfloor\frac{D-8}{2}\right\rfloor.
\]

Thus the first nontrivial block appears at degree 10, and degree 14 gives the
basis \([1,3]\).

## 5. Code

`collective_transport.py` provides

```python
collective_transport_polynomials()
collective_transport_expectation_matrix(auxiliary_degrees)
collective_transport_degrees(total_degree)
```

The expectation-matrix output has the same type as the existing
`four_point_hessian_expectation_matrix` output and is intended to be appended
to the solver's PSD `blocks` list:

```python
from collective_transport import (
    collective_transport_degrees,
    collective_transport_expectation_matrix,
)

transport_degrees = collective_transport_degrees(args.degree)
if args.collective_transport and transport_degrees:
    label_matrices = collective_transport_expectation_matrix(transport_degrees)
    variable = cp.Variable(
        (len(transport_degrees), len(transport_degrees)),
        symmetric=True,
        name="collective_transport_hessian",
    )
    blocks.append(("collective_transport_hessian", variable, label_matrices))
    constraints.append(variable >> 0)
```

The new block is a necessary condition for every local minimizer under smooth
transport variations.  It does not require isotropy analytically, although the
current moment reducer still uses isotropy when canonicalizing its entries.

## 6. Tests

`test_collective_transport.py` checks:

1. the symbolic local and cross Gram polynomials against direct tangent-vector
   contractions for explicit points on the sphere;
2. symmetry of the reduced moment matrices;
3. positive semidefiniteness on the uniform ONB equality measure;
4. degree scheduling at degrees 8, 10, and 14.

A natural next extension is to replace the one-leaf fields \(V_r\) by general
rooted vector-valued flags with additional conditioning roots.  That would be
the collective analogue of the current `--max-hessian-arity` hierarchy.
