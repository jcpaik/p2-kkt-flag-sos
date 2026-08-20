# Exact copositivity proof

Let

\[
K(t)=32t^6-48t^4+20t^2-\frac43,
\qquad
E(\mu)=\iint K(x\cdot y)\,d\mu(x)d\mu(y),
\]

where \(\mu\) is a probability measure on \(\mathbb {RP}^2\), represented
by unit vectors in \(S^2\).  This note proves

\[
\boxed{E(\mu)\geq0}
\]

for every \(\mu\), without an isotropy assumption.

## 1. Tangent--Pluecker reduction

Put \(W=\operatorname {Sym}^2_0(\mathbb R^3)\), with the Frobenius inner
product.  For a unit vector \(x\), let \(Q_x\) be the rank-two projection
onto

\[
W_x=\{\sqrt2\,x\odot v:v\perp x\}\subset W,
\]

and let \(R_x\) be the rank-one projection onto the Pluecker line
\(\bigwedge^2W_x\subset\bigwedge^2W\).  If \(s=(x\cdot y)^2\), direct
calculation in a tangent frame gives

\[
\operatorname {tr}(Q_xQ_y)=1-3s+4s^2,
\qquad
\operatorname {tr}(R_xR_y)=s(2s-1)^2.
\]

Consequently

\[
K(x\cdot y)=4\left(2\operatorname {tr}(R_xR_y)
-\operatorname {tr}(Q_xQ_y)+\frac23\right).
\]

Define

\[
F=\int Q_x\,d\mu(x),\qquad G=\int R_x\,d\mu(x).
\]

Then \(F\succeq0\), \(G\succeq0\), \(\operatorname {tr}F=2\),
\(\operatorname {tr}G=1\), and the one-body contraction of \(G\) is
\(F\).  Integrating the preceding identity yields

\[
\boxed{\frac{E(\mu)}8=\|G\|^2-\frac12\|F\|^2+\frac13.} \tag{1}
\]

## 2. The capped branch

Diagonalize \(F\) in an orthonormal basis \(e_0,\ldots,e_4\) of \(W\).
In the wedge basis put

\[
p_{ij}=\langle e_i\wedge e_j,G(e_i\wedge e_j)\rangle,
\qquad
d_i=\sum_{j\ne i}p_{ij}.
\]

The \(p_{ij}\) are nonnegative edge weights of total mass one on \(K_5\),
and \(d_i\) are the eigenvalues of \(F\).  If

\[
\mathcal A(p)=\sum_{e<f,\ e\cap f\ne\varnothing}p_ep_f
\]

and \(R_{\rm off}\) is the squared Hilbert--Schmidt norm of the off-diagonal
part of \(G\) in this wedge basis, then (1) becomes

\[
\frac{E(\mu)}8=\frac13-\mathcal A(p)+R_{\rm off}. \tag{2}
\]

The capped Motzkin--Straus lemma says that

\[
\max_i d_i\le\frac23\quad\Longrightarrow\quad\mathcal A(p)\le\frac13.
\]

For completeness, at a maximizer of minimal support, two disjoint positive
edges can be compressed until an edge disappears or a degree reaches
\(2/3\).  Pairwise-intersecting support is a triangle or a star; the star is
infeasible and the triangle gives \(\mathcal A\le1/3\).  If
\(d_0=2/3\), write the four incident weights as \(a_i\), the remaining
weights as \(b_{ij}\), and their outer degrees as \(r_i\).  Completing a
square gives

\[
\mathcal A=\frac29-\frac12\|a-r\|^2
+\frac12\|r\|^2+\mathcal A_b\le\frac29+\left(\sum b_{ij}\right)^2
=\frac13.
\]

Thus (2) proves \(E(\mu)\ge0\) whenever \(\lambda_{\max}(F)\le2/3\).

## 3. Exact remainder above the cap

It remains to consider a unit top eigenvector \(S\in W\),

\[
FS=cS,\qquad c>\frac23,
\qquad \Delta=3c-2>0.
\]

Set \(V=S^\perp\), and split

\[
\bigwedge^2W=(S\wedge V)\oplus\bigwedge^2V,
\qquad
G=\begin{pmatrix}A&C\\C^T&B\end{pmatrix}.
\]

Identify \(S\wedge V\) with \(V\), let \(R=\gamma_V(B)\) be the
one-body contraction of \(B\), and put

\[
U=A-R-\frac{\Delta}{4}I_V.
\]

If \(*\) is the Hodge involution on \(\bigwedge^2V\), the four-dimensional
contraction identity

\[
\|R\|^2=\|B\|^2+(\operatorname {tr}B)^2-
\operatorname {tr}(B*B*)
\]

and (1) give the exact decomposition

\[
\boxed{
\frac{E(\mu)}8=rac12\|U\|^2+2\|C\|^2
+\operatorname {tr}(B*B*)-\frac{\Delta^2}{24}.} \tag{3}
\]

Both \(B\) and \(*B*\) are positive semidefinite, so their trace product is
nonnegative.  It is therefore enough to prove

\[
2\|U\|^2+8\|C\|^2\ge\frac{\Delta^2}{6}. \tag{4}
\]

## 4. The spectral q-module certificate

After changing the sign of \(S\) and permuting the physical axes, every
unit traceless symmetric \(S\) has the form

\[
S=S_h=\frac{\sqrt3E_0-hE_1}{\sqrt{3+h^2}},\qquad0\le h\le1,
\]

where

\[
E_0=\frac{\operatorname {diag}(1,-1,0)}{\sqrt2},\qquad
E_1=\frac{\operatorname {diag}(1,1,-2)}{\sqrt6}.
\]

Complete this to the orthonormal basis

\[
S_h,\quad H_h=\frac{hE_0+\sqrt3E_1}{\sqrt{3+h^2}},
\quad T_{12},T_{13},T_{23}
\]

of \(W\), where \(T_{ij}=(e_ie_j^T+e_je_i^T)/\sqrt2\).
Average under the physical coordinate sign changes.  This leaves \(S_h\),
\(c\), and \(\Delta\) fixed.  Since the averaging is an orthogonal
projection,

\[
\|\bar U\|\le\|U\|,\qquad \|\bar C\|\le\|C\|. \tag{5}
\]

Write \(X=x_1^2,Y=x_2^2,Z=x_3^2\), so \(X,Y,Z\ge0\) and
\(X+Y+Z=1\).  Define

\[
\ell_h=(3-h)X-(3+h)Y+2hZ
\]

and

\[
\begin{aligned}
w_X&=\frac{-17h^2+32h+20}{70},\\
w_Y&=\frac{40h^2-73h+40}{140},\\
w_Z&=\frac{12-5h^2}{14},\\
g_h(X,Y,Z)&=(w_XX+w_YY+w_ZZ)\ell_h^2.
\end{aligned} \tag{6}
\]

All three weights are strictly positive on \([0,1]\): the numerator of
\(w_X\) is concave and positive at both endpoints, the numerator of
\(w_Y\) has negative discriminant, and \(12-5h^2\ge7\).  Hence

\[
g_h\ge0. \tag{7}
\]

Here is the exact algebraic certificate.  In the basis
\(V=(H_h,T_{12},T_{13},T_{23})\), let

\[
\mathbf u=(\bar U_{00},\bar U_{11},\bar U_{22},\bar U_{33})
\]

and, using the wedge order \(01,02,03,12,13,23\), let

\[
\mathbf c=(\bar C_{1,01},\bar C_{1,23},\bar C_{2,02},
\bar C_{2,13},\bar C_{3,03},\bar C_{3,12}).
\]

Direct expansion of the tangent Pluecker coordinates

\[
z_{ij}(x)=2x\cdot(E_ix\times E_jx)
\]

gives

\[
r_h:=\Delta+\int g_h\,d\mu=L(h)\cdot\mathbf u+M(h)\cdot\mathbf c. \tag{8}
\]

The four entries of \(L\) are

\[
\begin{aligned}
L_0&=-\frac{194h^4-435h^3+208h^2+81h+120}{280},\\
L_1&=-\frac{206h^4+435h^3-1168h^2-81h+440}{280},\\
L_2&= \frac{274h^4-101h^3+132h^2-753h+280}{280},\\
L_3&= \frac{126h^4+101h^3-1092h^2+753h+280}{280}.
\end{aligned} \tag{9}
\]

Put \(q=\sqrt{h^2+3}\).  The six entries of \(M\) are

\[
\begin{aligned}
M_0={}&-\frac{\sqrt3h(12h^5+512h^4-2779h^3-3829h^2+3847h+1741)}{840(h+1)},\\
M_1={}&-\frac{\sqrt3q(12h^5+364h^4-159h^3-90h^2+717h+960)}{1260(h+1)},\\
M_2={}& \frac{\sqrt3(12h^6-32h^5-983h^4-2373h^3-2405h^2+2645h+960)}{840(h+1)},\\
M_3={}&-\frac{\sqrt3q(6h^5+887h^4+582h^3-1155h^2+288h+480)}{1260(h+1)},\\
M_4={}& \frac{\sqrt3(12h^6+908h^5-545h^4-223h^3+4749h^2-85h-960)}{840(h+1)},\\
M_5={}& \frac{\sqrt3q(6h^5-523h^4-741h^3+1065h^2+429h+480)}{1260(h+1)}.
\end{aligned} \tag{10}
\]

This is a polynomial identity in the ten degree-three monomials in
\(X,Y,Z\); no inequality or numerical fitting enters (8).

The weighted dual cost of (8) is

\[
D(h):=\frac{\|L(h)\|^2}{2}+\frac{\|M(h)\|^2}{8}
=\frac{P(h)}{5644800(h+1)^2}, \tag{11}
\]

where

\[
\begin{aligned}
P(h)={}&1584h^{12}+117408h^{11}+10693280h^{10}+2137032h^9
+13971489h^8+45236046h^7\\
&-16424313h^6+4522956h^5+68504799h^4-88091010h^3
-11034711h^2+48384000h+24192000.
\end{aligned}
\]

Moreover \(D(h)<6\) on \([0,1]\).  Indeed,

\[
6-D(h)=\frac{R(h)}{5644800(h+1)^2},
\]

and the degree-twelve Bernstein coefficients of \(R\) on \([0,1]\) are

\[
\begin{gathered}
9676800,\ 11289600,\ \frac{298820637}{22},\ \frac{186523506}{11},
\ \frac{3576110762}{165},\ \frac{1838154161}{66},
\ \frac{10926245933}{308},\\
\frac{88743939}{2},\ \frac{2973164684}{55},
\ \frac{3483459192}{55},\ \frac{2287502816}{33},
\ \frac{194932736}{3},\ 33264640.
\end{gathered} \tag{12}
\]

Every number in (12) is positive, proving the claim exactly.

Weighted Cauchy--Schwarz applied to (8) now gives

\[
r_h^2\le D(h)\bigl(2\|\bar U\|^2+8\|\bar C\|^2\bigr)
\le6\bigl(2\|U\|^2+8\|C\|^2\bigr). \tag{13}
\]

By (7) and \(\Delta>0\), \(r_h\ge\Delta\).  Therefore (13) proves (4).
Substitution into (3) finishes the uncapped branch, while Section 2 handled
the capped branch.  Hence \(E(\mu)\ge0\) for every probability measure
\(\mu\) on \(\mathbb {RP}^2\).

The coefficient identity, the dual cost, and all Bernstein signs are
checked symbolically in `fermionic_spectral_qmodule_simple.py`.
