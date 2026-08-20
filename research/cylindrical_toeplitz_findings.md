# Cylindrical / Toeplitz findings (exact audit)

This note records exact identities obtained independently of `PLAN.md`.  It is
not itself a proof of copositivity; the final missing inequality is isolated in
§6.

## 1. Global irreducible identity

Put

\[
  r(x)=xx^{\mathsf T}-I/3\in V:=\operatorname{Sym}_0(3),\qquad
  m=\mathbb E r,\quad C=\mathbb E r^{\otimes2},\quad
  T=\mathbb E r^{\otimes3}.
\]

On `V`, use the traceless Jordan product

\[
 A*B=(AB+BA)/2-\operatorname{tr}(AB)I/3.
\]

For every rank-one projector, `r*r=r/3`.  Decompose under `SO(3)` as

\[
 C=C_0+C_2+C_4,\qquad T=T_0+T_2+T_4+T_6.
\]

The odd spin-3 component of `T` is absent because `r(x)^3` is an even
polynomial in `x`.  Exact Jordan contractions give

\[
 \begin{aligned}
 \|C_0\|^2&=4/45,&\quad \|C_2\|^2&={4\over21}\|m\|^2,\\
 \|T_0\|^2&=16/945,& \|T_2\|^2&={4\over21}\|m\|^2,
 &\|T_4\|^2&={4\over11}\|C_4\|^2.
 \end{aligned}
\]

For factor auditing: the squared singular values of the fully symmetric
contraction `Sym^3(V_2) -> V_2 tensor V_2` on spins `0,2,3,4` are

\[
       {7\over12},\quad {1\over9},\quad {5\over12},\quad {11\over36},
\]

and spin 6 is the kernel.  A symbolic check is in
`research/cylindrical_sym3_coupling.py`.

Writing `s=(x.y)^2-1/3=<r(x),r(y)>`, one has

\[
 K=32s^3-16s^2-{4\over3}s+{32\over27}.
\]

Consequently the following is an **exact identity for every measure**:

\[
 \boxed{
 {E(\mu)\over16}={2\over105}+{3\over28}\|m\|^2
       -{3\over11}\|C_4\|^2+2\|T_6\|^2. }
\]

Equivalently, in normalized Legendre energies
`I_l=double integral P_l(x.y)`, the target is

\[
 I_6+{33\over64}I_2+{11\over80}-{9\over20}I_4\ge0.
\]

Thus the *entire* cylindrical difficulty is one coupled spin-4 deficit.

## 2. Exact rank-one splits of the cylindrical mode matrices

With the notation of `docs/SUBCASES_AND_RECORD.md`, and using the principal
axis constraint `z_{1,0}=0`,

\[
 Q_1=792\left|z_{1,2}-{10\over11}z_{1,1}\right|^2
             -{336\over11}|z_{1,1}|^2.
\]

For `(x,y,z)=(z_{2,0},z_{2,1},z_{2,2})`,

\[
 Q_2={66\over17}\left|x-{52\over11}y+{60\over11}z\right|^2
       +{1080\over11}|y-2z|^2
       -{1\over17}|7x-30y+15z|^2.
\]

For `(A,B)=(z_{3,0},z_{3,1})`,

\[
 Q_3=270|B-A/3|^2-2|3A-5B|^2.
\]

Finally,

\[
 Q_4=66|z_{4,1}-z_{4,0}/11|^2-{6\over11}|z_{4,0}|^2.
\]

Dropping the displayed positive squares destroys sharpness on pole--equator
measures (in particular the mode-4 square cancels its negative term when the
height is zero).

## 3. A sharp measure-valued `3 -> 6` Toeplitz block

Set `s=a^2`, `r^2=1-s`, and

\[
 Z=\int ar^3(3-5s)\,d\mu_3,quad
 C_6=\int r^6\,d\mu_6,quad
 R=\int r^6\,d\sigma,quad
 H=\int s(3-5s)^2\,d\sigma.
\]

The integrated Gram matrix of
`(r^3, a(3-5s)e^{3i phi}, r^3e^{6i phi})` is

\[
 \begin{pmatrix}
 R&\bar Z&\bar C_6\\ Z&H&\bar Z\\ C_6&Z&R
 \end{pmatrix}\succeq0.
\]

It is sharp at the ONB height `s=1/3`, where `H=2R`.  In particular it must
be used as a matrix Schur complement; separate scalar Cauchy bounds lose the
ONB equality.

## 4. Tangent spin-4 determinant and ONB potential

For a root `x` and oriented tangent frame `(u,v)`, let

\[
 A_x=\int t^2(1-t^2)^2d\mu(y),\qquad
 B_x=\int t^2((u+iv)\mathbin\cdot y)^4d\mu(y),\quad t=x\mathbin\cdot y.
\]

Then `A_x >= |B_x|`, and

\[
 D_x=A_x^2-|B_x|^2,qquad F=\int D_xd\mu(x)
 ={8\over3}\mathbb E\!\left[D\sum_{cyc}a^2b^2(c-ab)^2\right]\ge0.
\]

For the ON frame

\[
 R(x,y)=\left(x,{y-(x\cdot y)x\over\sqrt{1-(x\cdot y)^2}},
 {x\times y\over\sqrt{1-(x\cdot y)^2}}\right)
\]

and `J_R(mu)=integral 32 product_j(R_j.y)^2 dmu(y)`,

\[
 F={1\over4}\mathbb E_{x,y}\left[(x\cdot y)^2(1-(x\cdot y)^2)^2
 J_{R(x,y)}(\mu)\right].
\]

More generally, among ONBs containing a fixed root `x`,

\[
 J_{\min}(x)=4(A_x-|B_x|),\quad J_{\max}(x)=4(A_x+|B_x|),
 \quad 16D_x=J_{\min}(x)J_{\max}(x).
\]

Thus the numerically plausible universal strengthening `E >= 16F` is
exactly

\[
 E\ge\int J_{\min}(x)J_{\max}(x)d\mu(x).
\]

## 5. Compound-operator form of `F`

Let `H` be the `H_3 -> H_3` block of the raw cubic moment matrix
`M=E |x^3><x^3|`.  At a root `x`, define harmonic cubics

\[
 f_c(y)=(x\cdot y)((u\cdot y)^2-(v\cdot y)^2),\quad
 f_s(y)=2(x\cdot y)(u\cdot y)(v\cdot y).
\]

Their coefficient tensors are orthogonal and both have squared norm `2/3`.
For the normalized basis `e_{c,s}=sqrt(3/2)f_{c,s}`,

\[
 \det(H|_{\operatorname{span}(e_c,e_s)})={9\over16}D_x.
\]

Hence, on `wedge^2 H_3` (dimension 21),

\[
 \boxed{F={16\over9}\operatorname{tr}((\wedge^2H)W_\mu)},\qquad
 W_\mu=\int|e_c(x)\wedge e_s(x)\rangle\langle\cdot|d\mu(x)\succeq0,
 \quad\operatorname{tr}W_\mu=1.
\]

This is the precise full measure-valued Toeplitz coupling represented by the
tangent determinant.

## 6. Exact remaining obstruction

The following attractive reductions are false:

1. `E >= 4 integral(A-|B|) dmu`.  The exact pole--equator perturbation
   described by the tensor strategy violates it already to leading order.
2. A pointwise bound of `D_x` by a fixed multiple of the rooted potential.
   Rooted potentials can be negative while `D_x >= 0`.
3. A sharp isotropic bound `E >= (105/2)F` does not extend to arbitrary
   measures: pole--equator perturbations have `E/F -> 24`.

No factor or sign issue remains in the reduction.  A complete proof by this
route is equivalent to either of the following exact global inequalities:

\[
 {2\over105}+{3\over28}\|m\|^2+2\|T_6\|^2
       \ge {3\over11}\|C_4\|^2,
\]

or the stronger (still unproved here) compound inequality

\[
 E/16\ge F={16\over9}\operatorname{tr}((\wedge^2H)W_\mu).
\]

Any scalarization of the `1 -> 2 -> 4` and `3 -> 6` chains before summing
fails on the documented equality faces; the needed certificate must retain
the operator-valued Schur complement through the final step.

