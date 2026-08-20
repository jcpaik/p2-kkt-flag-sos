# Exact KKT/contact-set reductions

This note records exact facts for the variational/contact-sextic route.  It is
independent of `PLAN.md`.  None of the results below is an epsilon estimate.

Throughout,

\[
 K(t)=32t^6-48t^4+20t^2-\frac43,
 \qquad E(\mu)=\iint K(x\mathbin\cdot y)\,d\mu(x)d\mu(y).
\]

## 1. What a hypothetical negative minimizer gives

Weak-* compactness gives a global minimizer.  If its energy is negative, write
it as \(\mu\), put \(E=E(\mu)\), and define

\[
 q(z)=U_\mu(z)-E,\qquad U_\mu(z)=\int K(z\mathbin\cdot y)\,d\mu(y).
\]

Then \(q\geq0\) on the sphere and \(q=0\) on \(\operatorname{supp}\mu\).
The homogeneous gap is the psd ternary sextic

\[
 Q(z)=32\int(z\mathbin\cdot y)^6d\mu(y)
 -48|z|^2\int(z\mathbin\cdot y)^4d\mu(y)
 +20|z|^4\int(z\mathbin\cdot y)^2d\mu(y)
 -(\tfrac43+E)|z|^6.                                      \tag{1}
\]

There is also a second variation that is sometimes omitted.  If \(\eta\) is
any finite signed measure supported on \(\operatorname{supp}\mu\), with
\(\eta(1)=0\), then \(\mu+t\eta\) is feasible for sufficiently small positive
and negative \(t\) after first restricting to bounded densities and then
approximating.  Hence

\[
 \iint K(x\mathbin\cdot y)\,d\eta(x)d\eta(y)\geq0.          \tag{2}
\]

For an atomic measure \(\mu=\sum_iw_i\delta_{x_i}\), (2) says that the kernel
matrix \(K_{ij}=K(x_i\cdot x_j)\) is positive semidefinite on
\(\mathbf1^\perp\).  The equilibrium equations are

\[
 Kw=E\mathbf1.                                             \tag{3}
\]

In harmonic feature notation, with

\[
 K(x\cdot y)=a_0+\langle u(x),u(y)\rangle
             -\langle v(x),v(y)\rangle,
\]

where \(u=\sqrt{8/7}\,Y_2\oplus\sqrt{512/231}\,Y_6\) and
\(v=\sqrt{384/385}\,Y_4\), condition (2) has the following exact form.  There
is a contraction \(T\) such that

\[
 v(x)-\bar v=T\bigl(u(x)-\bar u\bigr)
 \quad\text{for }\mu\text{-almost every }x.                \tag{4}
\]

This is the Douglas factorization lemma applied to the two feature-integration
maps on the mean-zero subspace of \(L^2(\mu)\).  Thus the full weight second
variation is strictly stronger than contact of the sextic.

## 2. Finite versus infinite contact sets

Every real zero of a psd form is singular.  A nonzero psd ternary sextic with
a finite real projective zero set has at most ten zeros (the sharp
Choi--Lam--Reznick bound).  Consequently, if \(Z_{\mathbb R}(Q)\) is finite,
then \(\mu\) is atomic on at most ten projective lines, and (2)--(3) apply.

If the real zero set is infinite, its Zariski closure contains a real curve.
Let \(g\) be a square-free irreducible equation of one such component.  Since
all partial derivatives of \(Q\) vanish on a Zariski-dense subset of
\(V(g)\), writing \(Q=gh\) and differentiating shows \(g\mid h\).  Therefore

\[
 g^2\mid Q,\qquad \deg g\leq3.                             \tag{5}
\]

The quotient is nonnegative off \(V(g)\), hence everywhere by continuity.
This is the elementary factorization behind the infinite-contact branch.
It does **not** by itself imply pole--equator structure; Section 4 gives an
exact counterexample.

If the support is Zariski dense in \(V(g)\), (2) gives one further algebraic
condition.  Regard every harmonic as a homogeneous sextic by multiplying by
the appropriate power of \(|x|^2\), and put

\[
 V_+=H_0\oplus H_2\oplus H_6,\qquad V_-=H_4.
\]

If a signed measure on the support annihilates \(V_+\), (2) forces it to
annihilate \(V_-\).  Finite-dimensional duality and Zariski density therefore
give

\[
 H_4\subset (H_0\oplus H_2\oplus H_6)+g\,\mathbb R[x,y,z]_{6-\deg g}. \tag{6}
\]

The normed version of (6) is exactly the contraction condition (4).  Thus an
infinite-support proof has to use an operator-norm restriction on the factor
\(g\), not factorization alone.

## 3. The ten-zero Robinson contact sextic is exact KKT, but not stable

Take the six edge lines

\[
 [1:\pm1:0],\quad[1:0:\pm1],\quad[0:1:\pm1]
\]

and the four diagonal lines

\[
 [1:\pm1:\pm1].
\]

Give each edge line weight \(184/2265\) and each diagonal line weight
\(387/3020\).  The weights sum to one.  Orbit averaging gives

\[
 E_{EE}=1,\qquad E_{ED}=-\frac{16}{27},\qquad
 E_{DD}=\frac{224}{243},
\]

and the equal-potential equation gives total edge weight \(368/755\).
The energy is

\[
 E=\frac{416}{2265}>0.                                    \tag{7}
\]

Direct exact expansion of (1) gives

\[
 Q(x,y,z)=\frac{64}{151}R(x,y,z),                           \tag{8}
\]

where

\[
 R=x^6+y^6+z^6-\sum_{\rm sym}x^4y^2+3x^2y^2z^2.
\]

Writing \(a=x^2,b=y^2,c=z^2\), this is Schur's form
\(\sum a(a-b)(a-c)\), so it is nonnegative.  Its ten projective zeros are
exactly the six edge lines and four diagonal lines above.  Thus a KKT proof
cannot assume that the contact sextic is SOS, reducible, or has at most nine
zeros.

This example is eliminated exactly by the weight second variation.  In the
ordering

\[
 [0:1:1],[0:1:-1],[1:0:1],[1:0:-1],[1:1:0],[1:-1:0],
 \text{ four diagonals},
\]

the kernel matrix has characteristic polynomial

\[
 (\lambda+1)^2
 (243\lambda^2-2354\lambda+3328)
 (729\lambda^2-4612\lambda+3584)^3/94143178827.
\]

In particular, the mass-zero vectors

\[
 (-1,-1,1,1,0,0,0,0,0,0),\qquad
 (-1,-1,0,0,1,1,0,0,0,0)
\]

are eigenvectors of eigenvalue \(-1\).  Hence (2) fails.

## 4. An infinite contact curve that is not pole--equator

Fix a pole \(e\), and let \(\mu_a\) be Haar measure on the projective latitude

\[
 C_a=\{[x]:(x\cdot e)^2=a\}.
\]

If \(v=(z\cdot e)^2\), its potential is

\[
\begin{aligned}
F(a,v)=\frac23(&693a^3v^3-945a^3v^2+315a^3v-15a^3
-945a^2v^3+1260a^2v^2-405a^2v+18a^2\\
&+315av^3-405av^2+126av-6a
-15v^3+18v^2-6v+1).
\end{aligned}                                               \tag{9}
\]

The contact/stationarity equation at its own latitude is

\[
 P(a):=693a^5-1575a^4+1260a^3-420a^2+54a-2=0.              \tag{10}
\]

There is exactly one root \(a_*\) of (10) in
\([421/1000,8/19]\), since the endpoint signs are respectively negative and
positive and a Sturm count is one.  Numerically only for identification,
\(a_*=0.4210393869\ldots\).  Modulo (10),

\[
 F(a,v)-F(a,a)=(v-a)^2(C(a)v+D(a)),                         \tag{11}
\]

where

\[
 C=2(231a^3-315a^2+105a-5),
\]
\[
 D=2(462a^4-945a^3+630a^2-145a+6).
\]

On the isolating interval, Bernstein coefficients show
\(C>0\) and \(D-C>0\).  Hence the third root \(-D/C\) is less than \(-1\),
so (11) is nonnegative for every \(v\in[0,1]\).  The energy is

\[
 F(a_*,a_*)>0
\]

(positivity again follows from positive Bernstein coefficients on the same
interval).  Thus this is an exact first-order KKT/contact example with an
infinite non-pole--equator contact curve.

It is rejected by (2).  On the latitude, the Fourier coefficient of mode two
of the restricted kernel is

\[
 \widehat K_2=\frac{(1-a)^2}{2}
 (495a^4-540a^3+162a^2-12a+1).                             \tag{12}
\]

All Bernstein coefficients of the quartic in (12) are negative on
\([421/1000,8/19]\), so \(\widehat K_2<0\).  A mode-two signed density has
zero mass and negative quadratic energy, contradicting (2).

## 5. A sharp geometric endpoint: the spin-four determinant

For a root \(x\), choose an oriented tangent orthonormal frame \((u,v)\).  For
a leaf \(y\), write

\[
 t=x\cdot y,\qquad \zeta=(u+iv)\cdot y,
\]

and define

\[
 A(x)=\int t^2|\zeta|^4d\mu(y),\qquad
 B(x)=\int t^2\zeta^4d\mu(y).
\]

The matrix

\[
 \begin{pmatrix}A&B\\\bar B&A\end{pmatrix}
\]

is psd, since it is the integral of rank-one psd matrices.  Therefore

\[
 \mathcal F(\mu):=\int(A(x)^2-|B(x)|^2)d\mu(x)\geq0.        \tag{13}
\]

For iid \(X,Y,Z\), put

\[
 a=X\cdot Y,\quad b=X\cdot Z,\quad c=Y\cdot Z,\quad
 D=1+2abc-a^2-b^2-c^2.
\]

The exact leaf expansion is

\[
 \boxed{\mathcal F(\mu)=8\,\mathbb E[
 a^2b^2(c-ab)^2D]}.                                       \tag{14}
\]

Indeed, if \(R=c-ab\) and \(J=\det(X,Y,Z)\), then
\(\zeta_Y\bar\zeta_Z=R+iJ\), and

\[
 (R^2+J^2)^2-\Re(R+iJ)^4=8R^2J^2.
\]

The zero locus of (13) has a complete support classification.

**Proposition.** If \(\mathcal F(\mu)=0\), then the support of \(\mu\) is
contained in a great circle, or in the union of a great circle and its
perpendicular pole.  Consequently \(E(\mu)\geq0\).

**Proof.**  The integrand in (14) is continuous and nonnegative.  If its
integral is zero, it vanishes on all of
\((\operatorname{supp}\mu)^3\).  For a linearly independent triple, \(D>0\).
Applying (14) in each cyclic rooting shows

\[
 a^2b^2(c-ab)^2=b^2c^2(a-bc)^2=c^2a^2(b-ca)^2=0.           \tag{15}
\]

If exactly two of \(a,b,c\) were nonzero, one equation in (15) would fail.  If
all three were nonzero, (15) would give \(c=ab,b=ac,a=bc\), hence
\(a^2=b^2=c^2=1\), contradicting linear independence.  Thus every independent
support triple has at most one nonorthogonal pair.

If the support does not span \(\mathbb R^3\), it lies in a great circle.
Otherwise choose an independent support triple.  If exactly one pair, say
\(y,z\), is nonorthogonal, then \(x\perp y,z\), so
\(P=\operatorname{span}(y,z)=x^\perp\).  Write any other support vector as
\(w=\alpha x+p\), \(p\in P\).  If both \(\alpha\) and \(p\) are nonzero,
applying the preceding triple property to \((x,y,w)\) and \((x,z,w)\)
(using the remaining triple if \(p\) is parallel to one of \(y,z\)) gives a
contradiction.  Hence \(w\in P\) or \(w\parallel x\).

If the chosen triple is an orthonormal basis, the same property first shows
that every other support vector has at most two nonzero coordinates.  Thus
every non-axis support line lies in a coordinate plane.  Non-axis lines in two
different coordinate planes, together with a suitable shared basis axis,
form an independent triple with two nonorthogonal pairs, a contradiction.
So again the support lies in one plane plus its perpendicular axis.

The great-circle energy is
\(|\widehat\nu(2)|^2+|\widehat\nu(6)|^2+2/3\).  On a plane plus its pole, the
exact pole--equator identity is

\[
 E=6(w-\tfrac13)^2+(1-w)^2
 (|\widehat\nu(2)|^2+|\widehat\nu(6)|^2)\geq0.
\]

This proves the proposition. \(\square\)

There is also an exact bridge from (13) to orthonormal-basis potentials.  For
every orthonormal frame \((x,u,v)\),

\[
 K(x\cdot y)+K(u\cdot y)+K(v\cdot y)
 =96(x\cdot y)^2(u\cdot y)^2(v\cdot y)^2.                 \tag{16}
\]

Rotating \((u,v)\) around \(x\), the phase that realizes \(|B(x)|\) gives

\[
 A(x)-|B(x)|
 =8\min_{(u,v)}\int(x\cdot y)^2(u\cdot y)^2(v\cdot y)^2d\mu(y).
\]

At a support root, where \(U(x)=E\), this becomes

\[
 A(x)-|B(x)|
 =\frac1{12}\left(E+\min_{(u,v)}(U(u)+U(v))\right).       \tag{17}
\]

The immediate consequence of \(U\geq E\) is only
\(A(x)-|B(x)|\geq E/4\), which is vacuous when \(E<0\); nevertheless (17)
identifies the precise ONB-completion quantity whose saturation is
\(\mathcal F=0\).

For later use, the cone generated by ONB measures is itself an exactly solved
copositive cone.  If \(R=(r_i)\) and \(S=(s_j)\) are two orthonormal frames,
then

\[
 I(\nu_R,\nu_S)
 =\frac{32}{3}\sum_{j=1}^3\prod_{i=1}^3(r_i\cdot s_j)^2\geq0,             \tag{18}
\]

where \(\nu_R=\frac13\sum_i\delta_{r_i}\).  Indeed the ONB potential is
\(U_{\nu_R}(z)=32\prod_i(r_i\cdot z)^2\).  Integrating (18) over any two
probability distributions of frames proves that every mixture of ONB measures
has nonnegative energy.

Therefore an identity showing that a hypothetical negative global minimizer
has \(\mathcal F(\mu)=0\) would finish the entire problem.  Contact
factorization alone cannot supply that identity: Sections 3 and 4 exhibit,
respectively, finite and infinite exact first-order KKT examples with
\(\mathcal F>0\), both eliminated only after invoking second variation.

## 6. Remaining exact obstruction

The unresolved assertion in this route is:

> If \(q\geq0\), \(q=0\) on \(\operatorname{supp}\mu\), and the full mass and
> positional second variations of \(E\) are nonnegative, then
> \(\mathcal F(\mu)=0\).

Equivalently, in the finite branch one must prove that no zero set of a psd
ternary sextic with at most ten projective points supports positive weights
satisfying (2)--(3) with \(E<0\).  The general classification of psd ternary
sextics does not prove this: non-SOS sextics have a 17-dimensional family of
ten-zero configurations, not only the Robinson configuration.  In the
infinite branch, (5)--(6) leave quadratic and cubic curve factors, and the
latitude example shows why the contraction norm in (4) is essential.

## 7. The ONB-potential cone: exact polar and an exact obstruction

Let \(H_3\) be the seven-dimensional space of harmonic ternary cubics.  The
multiplication map

\[
 \Gamma:\operatorname{Sym}^2(H_3)\longrightarrow
 \mathbb R[x,y,z]_6,
 \qquad \Gamma(G)(z)=h(z)^{\mathsf T}G h(z),               \tag{19}
\]

is an isomorphism.  One quick proof is the multiplicity-free decomposition

\[
 \operatorname{Sym}^2(H_3)=H_0\oplus H_2\oplus H_4\oplus H_6
 =\mathbb R[x,y,z]_6;
\]

the multiplication map is equivariant and nonzero on each summand.  For an
orthonormal frame \(R=(r_1,r_2,r_3)\), put

\[
 v_R(z)=(r_1\cdot z)(r_2\cdot z)(r_3\cdot z)\in H_3.
\]

Thus the closed cone generated by ONB potentials is

\[
 \mathcal C_{\rm ONB}
 =\Gamma\!\left(\operatorname{cone}
       \{32v_Rv_R^{\mathsf T}:R\in SO(3)\}\right).        \tag{20}
\]

Using the pullback trace pairing
\(L_A(\Gamma(G))=\operatorname{tr}(AG)\), its polar has the exact
description

\[
 \boxed{\mathcal C_{\rm ONB}^*
   =\{A\in\operatorname{Sym}(H_3):
             v_R^{\mathsf T}Av_R\geq0\text{ for every }R\in SO(3)\}.}   \tag{21}
\]

This also exposes a fatal gap in an argument based only on the polynomial
facts \(q\geq0\), \(q\in H_0\oplus H_4\oplus H_6\), and \(E\leq0\).
Consider the Robinson--Schur sextic

\[
 R=x^6+y^6+z^6-\sum_{\rm sym}x^4y^2+3x^2y^2z^2.
\]

It is nonnegative by Schur's inequality applied to
\(a=x^2,b=y^2,c=z^2\).  Moreover

\[
 \Delta^2R=240(x^2+y^2+z^2),\qquad
 \int_{S^2}R\,d\sigma=\frac27,                            \tag{22}
\]

so it has no \(H_2\) part.  Consequently

\[
 q_*=\frac{16}{15}R                                      \tag{23}
\]

is nonnegative, lies in \(H_0\oplus H_4\oplus H_6\), and has constant
component \(32/105\), exactly the value corresponding to the formal label
\(E=0\).  Any \(cR\) with \(c>16/15\) has constant component greater than
\(32/105\), hence corresponds formally to \(E<0\).

Nevertheless none of these forms belongs to \(\mathcal C_{\rm ONB}\).
Indeed (20) consists of sums of squares of cubics, whereas \(R\) is not a
sum of squares.  Here is a completely elementary verification of the latter
fact.  Its ten projective zeros are the six edge lines

\[
 [1:\pm1:0],\quad[1:0:\pm1],\quad[0:1:\pm1]
\]

and the four body-diagonal lines \([1:\pm1:\pm1]\).  Evaluation of the ten
cubic monomials at these ten points gives a \(10\times10\) matrix of
determinant \(-128\).  Thus no nonzero cubic vanishes at all ten points.  If
\(R=\sum_j f_j^2\), every \(f_j\) would vanish at all ten zeros, forcing
every \(f_j=0\), a contradiction.

There is also an explicit rational polar separator.  In the harmonic-cubic
basis

\[
\begin{aligned}
 &(z(5z^2-3r^2),\ x(5z^2-r^2),\ y(5z^2-r^2),\ z(x^2-y^2),\\
 &\hspace{38mm}2xyz,\ x(x^2-3y^2),\ y(3x^2-y^2)),
\end{aligned}
\]

the unique matrix \(G_R=\Gamma^{-1}(R)\) is

\[
\begin{pmatrix}
1/4&0&0&0&0&0&0\\
0&1/8&0&0&0&-1/4&0\\
0&0&1/8&0&0&0&1/4\\
0&0&0&-1/4&0&0&0\\
0&0&0&0&-3&0&0\\
0&-1/4&0&0&0&3/8&0\\
0&0&1/4&0&0&0&3/8
\end{pmatrix}.                                                \tag{24}
\]

Taking \(A\) to be the rank-one projector onto the fifth basis coordinate
gives \(v_R^{\mathsf T}Av_R\geq0\) for every frame, so \(A\) belongs to the
polar (21), while \(L_A(R)=\operatorname{tr}(AG_R)=-3<0\).

Therefore

\[
 q\geq0,\quad q\in H_0\oplus H_4\oplus H_6,\quad
 q_0\geq\frac{32}{105}                                   \tag{25}
\]

does **not** imply \(q\in\mathcal C_{\rm ONB}\).  The scaled Robinson form
is not claimed to be \(U_\mu-E(\mu)\) for a probability measure; rather,
it proves exactly that any ONB-cone argument must use the missing moment
self-consistency/contact-support or second-variation information.  Positivity
and the sign of the constant term alone cannot establish cone membership.

For this obstruction the missing self-consistency can itself be seen exactly.
Suppose a probability measure supported on the ten zeros above had gap
\(q=cR\).  Writing arbitrary weights on the ten lines, expanding
\(U_\mu-E|z|^6=cR\), and using the forced constant relation

\[
 E=\frac{32}{105}-\frac{2c}{7},
\]

gives a linear system with the unique solution

\[
 c=\frac{64}{151},\qquad
 w_{\rm edge}=\frac{184}{2265},\qquad
 w_{\rm diagonal}=\frac{387}{3020}.
\]

It is exactly the Robinson KKT measure of Section 3 and has
\(E=416/2265>0\).  In particular the formally zero-energy scaling (23) is
not the gap of even a signed measure on its own contact set.  This calculation
pinpoints the extra information that a successful strengthening of the cone
route would have to exploit.
