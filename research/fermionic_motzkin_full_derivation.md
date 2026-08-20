# From copositivity to the finite-dimensional tangent--Plucker inequality

This note gives a self-contained derivation of the finite-dimensional
reduction for the energy

\[
E(\mu):=\iint K(x\cdot y)\,d\mu(x)d\mu(y),
\qquad
K(t)=32t^6-48t^4+20t^2-\frac43.
\]

Here \(\mu\) is a probability measure on the real projective plane. We
represent projective points by unit vectors \(x\in S^2\). Every expression
below is even in \(x\), so the choice of representative sign is irrelevant.

The final result is that the original measure problem reduces to one sharp
quadratic inequality on a 28-dimensional tangent--Veronese moment body.

## 1. A five-dimensional model

Let

\[
W:=\operatorname{Sym}^2_0(\mathbb R^3)
\]

be the space of real symmetric traceless \(3\times3\) matrices, equipped with
the Frobenius inner product

\[
\langle A,B\rangle=\operatorname{tr}(AB).
\]

This space has dimension five.

For \(x\in S^2\), define the two-dimensional subspace

\[
W_x:=\left\{\sqrt2\,x\odot v:v\perp x\right\}\subset W,
\]

where

\[
x\odot v=\frac{xv^T+vx^T}{2}.
\]

If \(v\perp x\) and \(\|v\|=1\), then \(\sqrt2\,x\odot v\) has Frobenius
norm one. Thus an orthonormal basis \(v_1,v_2\) of \(x^\perp\) gives an
orthonormal basis

\[
u_i(x)=\sqrt2\,x\odot v_i,\qquad i=1,2,
\]

of \(W_x\).

Let

\[
Q_x:W\longrightarrow W
\]

be the orthogonal projection onto \(W_x\). Hence \(Q_x\) is a rank-two
projection. Geometrically, the spaces \(W_x\) are the tangent planes of the
Veronese surface

\[
x\longmapsto xx^T-\frac13I
\]

inside \(W\).

## 2. Passing to Plucker vectors

The oriented two-plane \(W_x\) determines the unit bivector

\[
z_x:=u_1(x)\wedge u_2(x)\in\bigwedge^2W.
\]

Changing the orientation of \(v_1,v_2\) changes the sign of \(z_x\), but the
rank-one projection

\[
R_x:=|z_x\rangle\langle z_x|
\]

does not depend on that choice. Thus \(Q_x\) is a rank-two projection on
\(W\), while \(R_x\) is a rank-one projection on the ten-dimensional space
\(\bigwedge^2W\).

We next calculate their pairwise overlaps.

## 3. The two overlap polynomials

Set

\[
s=(x\cdot y)^2.
\]

By rotational invariance, take

\[
x=e_1,
\qquad
y=t e_1+\sqrt{1-t^2}\,e_2,
\qquad
s=t^2.
\]

Choose the tangent frames

\[
v_1=e_2,\quad v_2=e_3
\]

at \(x\), and

\[
w_1=-\sqrt{1-t^2}\,e_1+t e_2,
\quad
w_2=e_3
\]

at \(y\).

For unit orthogonal pairs \(a\perp b\) and \(c\perp d\),

\[
\left\langle
\sqrt2\,a\odot b,\sqrt2\,c\odot d
\right\rangle
=(a\cdot c)(b\cdot d)+(a\cdot d)(b\cdot c).
\]

Therefore the \(2\times2\) overlap matrix between the two tangent frames is

\[
\begin{pmatrix}
2t^2-1&0\\
0&t
\end{pmatrix}.
\]

It follows that

\[
\operatorname{tr}(Q_xQ_y)
=(2t^2-1)^2+t^2
=1-3s+4s^2. \tag{3.1}
\]

For the Plucker projections,

\[
\operatorname{tr}(R_xR_y)=|\langle z_x,z_y\rangle|^2.
\]

The inner product of two unit simple bivectors is the determinant of the
corresponding plane-overlap matrix. Hence

\[
\operatorname{tr}(R_xR_y)
=t^2(2t^2-1)^2
=s(2s-1)^2. \tag{3.2}
\]

Write

\[
q(s):=1-3s+4s^2,
\qquad
r(s):=s(2s-1)^2.
\]

A direct calculation gives

\[
4\left(2r(s)-q(s)+\frac23\right)
=32s^3-48s^2+20s-\frac43.
\]

Since \(s=(x\cdot y)^2\), this is exactly

\[
\boxed{
K(x\cdot y)
=4\left(
2\operatorname{tr}(R_xR_y)
-\operatorname{tr}(Q_xQ_y)
+\frac23
\right).} \tag{3.3}
\]

This is the fundamental identity.

## 4. The moment operators \(F\) and \(G\)

Define

\[
F:=\int Q_x\,d\mu(x)
\quad\text{on }W,
\]

and

\[
G:=\int R_x\,d\mu(x)
\quad\text{on }\bigwedge^2W.
\]

Because \(Q_x\) has rank two and \(R_x\) has rank one,

\[
F\succeq0,
\qquad
G\succeq0,
\]

and

\[
\operatorname{tr}F=2,
\qquad
\operatorname{tr}G=1.
\]

Moreover,

\[
\iint\operatorname{tr}(Q_xQ_y)\,d\mu(x)d\mu(y)
=\operatorname{tr}(F^2)=\|F\|_{HS}^2,
\]

and similarly

\[
\iint\operatorname{tr}(R_xR_y)\,d\mu(x)d\mu(y)
=\operatorname{tr}(G^2)=\|G\|_{HS}^2.
\]

Integrating (3.3) gives

\[
E(\mu)
=4\left(2\|G\|_{HS}^2-\|F\|_{HS}^2+\frac23\right).
\]

Equivalently,

\[
\boxed{
\frac{E(\mu)}8
=\|G\|_{HS}^2
-\frac12\|F\|_{HS}^2
+\frac13.} \tag{4.1}
\]

This turns the original kernel problem into a purity inequality between
\(G\) and its one-body contraction \(F\).

## 5. Why \(F\) is the one-body contraction of \(G\)

For an operator \(H\) on \(\bigwedge^2W\), define its contraction
\(\gamma(H)\) by

\[
\langle e_i,\gamma(H)e_k\rangle
=\sum_j
\left\langle
e_i\wedge e_j,
H(e_k\wedge e_j)
\right\rangle
\]

in any orthonormal basis \((e_i)\) of \(W\).

If

\[
H=|u\wedge v\rangle\langle u\wedge v|
\]

with \(u,v\) orthonormal, then

\[
\gamma(H)=|u\rangle\langle u|+|v\rangle\langle v|.
\]

Consequently,

\[
\gamma(R_x)=Q_x,
\]

and hence

\[
\boxed{F=\gamma(G).} \tag{5.1}
\]

Thus \(F\) and \(G\) are not independent positive operators.

## 6. Pinching in an eigenbasis of \(F\)

Choose an orthonormal eigenbasis

\[
e_0,\ldots,e_4
\]

of \(F\). Use the ten bivectors \(e_i\wedge e_j\), \(i<j\), as a basis of
\(\bigwedge^2W\).

Define

\[
p_{ij}
=\left\langle
e_i\wedge e_j,
G(e_i\wedge e_j)
\right\rangle.
\]

Because \(G\succeq0\),

\[
p_{ij}\ge0,
\qquad
\sum_{i<j}p_{ij}=\operatorname{tr}G=1.
\]

Let

\[
d_i:=\sum_{j\ne i}p_{ij}.
\]

By the contraction identity, \(d_i\) is precisely the eigenvalue of \(F\)
corresponding to \(e_i\). In particular,

\[
\sum_i d_i=2.
\]

Regard the ten quantities \(p_{ij}\) as weights on the edges of the complete
graph \(K_5\). Then \(d_i\) is the weighted degree of vertex \(i\).

Define

\[
\mathcal A(p)
:=\sum_{\substack{e<f\\e\cap f\ne\varnothing}}p_ep_f,
\]

the sum of products over pairs of distinct adjacent edges. Finally let

\[
R_{\mathrm{off}}
:=\sum_{e\ne f}|G_{ef}|^2,
\]

where the sum is over ordered distinct bivector basis elements. Then

\[
\|G\|^2=\sum_ep_e^2+R_{\mathrm{off}}. \tag{6.1}
\]

On the other hand,

\[
\sum_i d_i^2=2\sum_ep_e^2+2\mathcal A(p),
\]

because each edge is incident to two vertices, while two distinct edges
contribute exactly when they meet. Therefore

\[
\frac12\|F\|^2=\sum_ep_e^2+\mathcal A(p). \tag{6.2}
\]

Substituting (6.1)--(6.2) into (4.1) gives the exact graph identity

\[
\boxed{
\frac{E(\mu)}8
=\frac13-\mathcal A(p)+R_{\mathrm{off}}.} \tag{6.3}
\]

## 7. The capped Motzkin--Straus lemma

We now prove the following statement.

> If the nonnegative edge weights \(p_e\) on \(K_5\) satisfy
> \[
> \sum_ep_e=1,
> \qquad
> d_i\le\frac23\quad\text{for every vertex }i,
> \]
> then
> \[
> \mathcal A(p)\le\frac13.
> \]

### 7.1. Some degree equals \(2/3\)

Suppose \(d_0=2/3\). Write

\[
a_i=p_{0i},\qquad 1\le i\le4,
\]

for the four star weights at vertex zero. Then

\[
\sum_i a_i=\frac23.
\]

Let \(b_{ij}=p_{ij}\), \(1\le i<j\le4\), be the six remaining weights.
Their total mass is \(1/3\). Define their outer degrees

\[
r_i=\sum_{\substack{j=1\\j\ne i}}^4b_{ij},
\]

and let \(\mathcal A_b\) be the adjacent-pair sum among the six \(b\)-edges.

Splitting adjacent pairs into star--star, star--outer, and outer--outer
pairs gives

\[
\mathcal A(p)
=\frac29-\frac12\|a\|^2+a\cdot r+\mathcal A_b.
\]

Completing the square,

\[
\mathcal A(p)
=\frac29
-\frac12\|a-r\|^2
+\frac12\|r\|^2
+\mathcal A_b.
\]

For the outer \(K_4\),

\[
\frac12\|r\|^2+\mathcal A_b
=\sum b_{ij}^2+2\mathcal A_b.
\]

This is at most

\[
\left(\sum b_{ij}\right)^2=\frac19,
\]

because the square of the total mass contains the same terms plus the
nonnegative products of disjoint edges. Therefore

\[
\mathcal A(p)\le\frac29+\frac19=\frac13.
\]

### 7.2. No degree constraint is active

Take a maximizer of \(\mathcal A\) having minimal support. If two positive
edges are disjoint, transfer weight between them while preserving their sum.

Because disjoint edges do not occur together in \(\mathcal A\), its value is
affine along this transfer. Move in a nondecreasing direction until either
one edge weight becomes zero, contradicting minimal support, or some degree
becomes \(2/3\), reducing to the preceding case.

Therefore, unless the preceding case occurs, all positive edges must be
pairwise intersecting. A pairwise-intersecting family of edges in \(K_5\) is
either contained in a star or contained in a triangle.

A star is impossible because its central degree would equal the total mass
\(1>2/3\). On a triangle,

\[
\mathcal A(p)
=p_1p_2+p_1p_3+p_2p_3
\le\frac{(p_1+p_2+p_3)^2}{3}
=\frac13.
\]

This proves the lemma. Consequently, whenever

\[
\lambda_{\max}(F)\le\frac23,
\]

equation (6.3) gives

\[
E(\mu)\ge0.
\]

## 8. Why isotropy was sufficient

If \(\mu\) is isotropic, meaning

\[
\int xx^T\,d\mu(x)=\frac13I,
\]

then for every \(S\in W\),

\[
\langle S,Q_xS\rangle
=2\left(|Sx|^2-(x^TSx)^2\right).
\]

Therefore, for \(\|S\|=1\),

\[
\begin{aligned}
\langle S,FS\rangle
&=2\int\left(|Sx|^2-(x^TSx)^2\right)d\mu(x)\\
&\le2\int|Sx|^2\,d\mu(x)\\
&=2\operatorname{tr}\left(S^2\int xx^T\,d\mu(x)\right)\\
&=\frac23.
\end{aligned}
\]

Thus

\[
F\preceq\frac23I.
\]

The isotropic proof is therefore exactly the capped Motzkin--Straus case.
Removing isotropy means handling the possibility

\[
c:=\lambda_{\max}(F)>\frac23.
\]

## 9. Exact scalar reduction when \(c>2/3\)

Choose vertex zero corresponding to a top eigenvector of \(F\), so

\[
c=d_0>\frac23.
\]

Write

\[
a_i=p_{0i},\qquad 1\le i\le4,
\]

and retain \(b_{ij}=p_{ij}\) for the six outer edges. Set

\[
r_i=\sum_{j\ne i}b_{ij},
\]

and define the disjoint-edge expression

\[
D_b=b_{12}b_{34}+b_{13}b_{24}+b_{14}b_{23}.
\]

Let

\[
\Delta:=3c-2>0,
\qquad
u:=a-r.
\]

Since

\[
\sum_i a_i=c,
\qquad
\sum_i r_i=2(1-c),
\]

we have

\[
\sum_i u_i=3c-2=\Delta.
\]

Remove the mean by setting

\[
u_0=u-\frac{\Delta}{4}\mathbf1.
\]

Then

\[
\|u\|^2=\|u_0\|^2+\frac{\Delta^2}{4}. \tag{9.1}
\]

The adjacency sum decomposes as

\[
\mathcal A(p)
=\frac{c^2}{2}-\frac12\|a\|^2+a\cdot r+\mathcal A_b.
\]

Completing the square gives

\[
\mathcal A(p)
=\frac{c^2}{2}
-\frac12\|a-r\|^2
+\frac12\|r\|^2+\mathcal A_b.
\]

For the outer graph,

\[
\frac12\|r\|^2+\mathcal A_b
=(1-c)^2-2D_b.
\]

Hence

\[
\mathcal A(p)
=\frac{c^2}{2}
+(1-c)^2
-2D_b
-\frac12\|u\|^2.
\]

Substituting this into (6.3), then using (9.1), yields

\[
\boxed{
\frac{E(\mu)}8
=R_{\mathrm{off}}
+2D_b
+\frac12\|u_0\|^2
-\frac{\Delta^2}{24}.} \tag{9.2}
\]

The elementary scalar identity used here is

\[
\frac13-\frac{c^2}{2}-(1-c)^2+\frac{(3c-2)^2}{8}
=-\frac{(3c-2)^2}{24}.
\]

Consequently, the whole theorem reduces to proving

\[
\boxed{
(3c-2)^2
\le
24R_{\mathrm{off}}
+48D_b
+12\|u_0\|^2} \tag{9.3}
\]

for moment matrices \(G\) arising from the tangent planes \(W_x\).

## 10. A basis-free form of the remaining inequality

The graph formula depends on a full eigenbasis of \(F\). There is a cleaner
invariant version.

Let \(S\in W\) be a unit top eigenvector of \(F\):

\[
FS=cS.
\]

Write

\[
V=S^\perp,
\qquad
\dim V=4.
\]

Then

\[
\bigwedge^2W=(S\wedge V)\oplus\bigwedge^2V.
\]

Identify \(S\wedge V\) with \(V\), and write \(G\) in block form

\[
G=
\begin{pmatrix}
\mathsf A&\mathsf C\\
\mathsf C^T&\mathsf B
\end{pmatrix},
\]

where \(\mathsf A\) acts on \(V\), \(\mathsf B\) acts on
\(\bigwedge^2V\), and \(\mathsf C:\bigwedge^2V\to V\).

Let

\[
\mathsf R:=\gamma_V(\mathsf B)
\]

be the one-body contraction of \(\mathsf B\) on \(V\). Then

\[
\operatorname{tr}\mathsf A=c,
\qquad
\operatorname{tr}\mathsf B=1-c,
\qquad
\operatorname{tr}\mathsf R=2(1-c).
\]

Because \(S\) is an eigenvector of \(F\), the \(S\)-to-\(V\) block of \(F\)
vanishes, while

\[
F|_V=\mathsf A+\mathsf R.
\]

Thus

\[
\|F\|^2=c^2+\|\mathsf A+\mathsf R\|^2, \tag{10.1}
\]

and

\[
\|G\|^2
=\|\mathsf A\|^2+2\|\mathsf C\|^2+\|\mathsf B\|^2. \tag{10.2}
\]

## 11. The four-dimensional Hodge identity

Let

\[
\star:\bigwedge^2V\to\bigwedge^2V
\]

be the Hodge involution, and write

\[
\star\mathsf B\star=\star\circ\mathsf B\circ\star.
\]

A direct contraction calculation in four dimensions gives

\[
\boxed{
\|\mathsf R\|^2
=\|\mathsf B\|^2
+(\operatorname{tr}\mathsf B)^2
-\operatorname{tr}(\mathsf B\star\mathsf B\star).} \tag{11.1}
\]

For example, in an orthonormal basis \(v_1,\ldots,v_4\) of \(V\), write

\[
R_{ik}
=\sum_j
\left\langle
v_i\wedge v_j,
\mathsf B(v_k\wedge v_j)
\right\rangle.
\]

Expanding \(\sum_{i,k}R_{ik}^2\) and pairing complementary bivectors under
\(\star\) gives (11.1).

Since \(\mathsf B\succeq0\) and \(\star\mathsf B\star\succeq0\),

\[
\operatorname{tr}(\mathsf B\star\mathsf B\star)\ge0. \tag{11.2}
\]

## 12. Completing the square

Substitute (10.1)--(10.2) into (4.1):

\[
\begin{aligned}
\frac E8
&=\|\mathsf A\|^2
+2\|\mathsf C\|^2
+\|\mathsf B\|^2\\
&\qquad
-\frac12c^2
-\frac12\|\mathsf A+\mathsf R\|^2
+\frac13.
\end{aligned}
\]

Using

\[
\|\mathsf A\|^2
-\frac12\|\mathsf A+\mathsf R\|^2
=\frac12\|\mathsf A-\mathsf R\|^2-\|\mathsf R\|^2,
\]

we obtain

\[
\frac E8
=\frac12\|\mathsf A-\mathsf R\|^2
+2\|\mathsf C\|^2
+\|\mathsf B\|^2-\|\mathsf R\|^2
-\frac12c^2+\frac13.
\]

Apply the Hodge identity and use \(\operatorname{tr}\mathsf B=1-c\):

\[
\begin{aligned}
\frac E8
&=\frac12\|\mathsf A-\mathsf R\|^2
+2\|\mathsf C\|^2\\
&\quad
+\operatorname{tr}(\mathsf B\star\mathsf B\star)
-(1-c)^2-\frac12c^2+\frac13.
\end{aligned}
\]

Now

\[
\operatorname{tr}(\mathsf A-\mathsf R)
=c-2(1-c)=3c-2=\Delta.
\]

Define the traceless matrix

\[
\mathsf U
:=\mathsf A-\mathsf R-\frac{\Delta}{4}I_V.
\]

Then

\[
\|\mathsf A-\mathsf R\|^2
=\|\mathsf U\|^2+\frac{\Delta^2}{4}.
\]

After simplifying the scalar terms, we get

\[
\boxed{
\frac{E(\mu)}8
=\frac12\|\mathsf U\|^2
+2\|\mathsf C\|^2
+\operatorname{tr}(\mathsf B\star\mathsf B\star)
-\frac{\Delta^2}{24}.} \tag{12.1}
\]

Therefore \(E(\mu)\ge0\) is equivalent, in the remaining region
\(c>2/3\), to

\[
\boxed{
\Delta^2
\le
12\|\mathsf U\|^2
+48\|\mathsf C\|^2
+24\operatorname{tr}(\mathsf B\star\mathsf B\star),
\qquad
\Delta=3c-2.} \tag{12.2}
\]

This is the finite-dimensional tangent-star inequality.

## 13. Why the tangent condition is indispensable

Inequality (12.2) is false for a general positive operator \(G\) on
\(\bigwedge^2W\).

Indeed, fix an orthonormal basis \(S,e_1,\ldots,e_4\) of \(W\), and set

\[
G_\star
=\frac14\sum_{i=1}^4
|S\wedge e_i\rangle\langle S\wedge e_i|.
\]

Then

\[
F=\operatorname{diag}\left(1,\frac14,\frac14,\frac14,\frac14\right).
\]

Thus \(c=1\), \(\Delta=1\), while

\[
\mathsf A=\frac14I_4,
\qquad
\mathsf B=0,
\qquad
\mathsf C=0,
\qquad
\mathsf U=0.
\]

So the right-hand side of (12.2) is zero, even though its left-hand side is
one. Correspondingly,

\[
E(G_\star)=-\frac13.
\]

This \(G_\star\) is a legitimate mixture of arbitrary two-planes in \(W\),
but it is not a mixture of the special tangent planes \(W_x\). Therefore the
final proof must use the coherent tangent--Veronese structure, not merely
positivity of \(G\).

## 14. Why the reduction is genuinely finite-dimensional

Choose an orthonormal basis \(E_0,\ldots,E_4\) of \(W\). The Plucker
coordinates of \(z_x\) are

\[
z_{ij}(x)=2x\cdot(E_ix\times E_jx),
\qquad 0\le i<j\le4.
\]

Each \(z_{ij}(x)\) is a homogeneous cubic polynomial in the three coordinates
of \(x\). Hence

\[
G_{ij,kl}
=\int z_{ij}(x)z_{kl}(x)\,d\mu(x)
\]

depends only on the homogeneous degree-six moments of \(\mu\).

There are

\[
\binom{6+3-1}{3-1}=28
\]

monomials of degree six in three variables. Thus all admissible \(G\)'s lie
in a 28-dimensional linear moment space inside the 55-dimensional space of
symmetric \(10\times10\) matrices.

Equivalently, define the compact tangent moment body

\[
\mathcal T
=\operatorname{conv}\left\{
|z_x\rangle\langle z_x|:x\in S^2
\right\}.
\]

The original copositivity theorem is now reduced to the following.

> For every \(G\in\mathcal T\), let \(F=\gamma(G)\).
>
> If \(\lambda_{\max}(F)\le2/3\), capped Motzkin--Straus proves \(E\ge0\).
>
> If \(c=\lambda_{\max}(F)>2/3\), prove the finite-dimensional quadratic
> inequality (12.2).

By Caratheodory's theorem, every \(G\in\mathcal T\) can also be represented
using a uniformly bounded number of atoms--at most 29 without optimizing the
affine dimension. Thus the remaining statement can literally be written as
a polynomial inequality in finitely many point coordinates and weights.

The unresolved issue is not infinite dimensionality. It is proving the sharp
inequality (12.2) using the special algebraic relations satisfied by the
tangent Plucker vectors \(z_x\).
