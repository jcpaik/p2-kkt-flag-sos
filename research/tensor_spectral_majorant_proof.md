# Exact tensor proof of copositivity

Let

\[
K(t)=32t^6-48t^4+20t^2-\frac43,
\qquad
E(\mu)=\iint K(x\cdot y)\,d\mu(x)d\mu(y).
\]

This note proves (E(\mu)\geq0) for every probability measure on
\(\mathbb {RP}^2\).  The long coefficient identity in the last step is
checked exactly, over the rational-function field, by
`tensor_spectral_majorant_exact.py`.

## 1. Tangent--Pluecker reduction

Put (W=\operatorname{Sym}^2_0(\mathbb R^3)), with the Frobenius inner
product.  For (x\in S^2), let (Q_x) be the rank-two orthogonal
projection onto

\[
W_x=\{\sqrt2\,x\odot v:v\perp x\}\subset W.
\]

Choose an oriented orthonormal basis of (W_x), let (z_x\in\Lambda^2W)
be its unit wedge, and put (R_x=|z_x\rangle\langle z_x|).  Directly
computing the overlap of tangent frames gives, with
\(s=(x\cdot y)^2\),

\[
\operatorname{tr}(Q_xQ_y)=1-3s+4s^2,
\qquad
\operatorname{tr}(R_xR_y)=s(2s-1)^2.
\]

Consequently

\[
K(x\cdot y)
=4\left(2\operatorname{tr}(R_xR_y)-\operatorname{tr}(Q_xQ_y)+\frac23\right).
\]

Define

\[
F=\int Q_x\,d\mu(x),\qquad G=\int R_x\,d\mu(x).
\]

Then (F=\gamma(G)), the one-particle contraction, and

\[
F\succeq0,quad G\succeq0,quad \operatorname{tr}F=2,quad
\operatorname{tr}G=1,
\]

while

\[
\boxed{\frac{E(\mu)}8=\|G\|_{HS}^2-\frac12\|F\|_{HS}^2+\frac13.}\tag{1}
\]

## 2. The capped case

Diagonalize (F) in an orthonormal basis (e_0,\ldots,e_4) of (W), and
write

\[
p_{ij}=\langle e_i\wedge e_j,G(e_i\wedge e_j)\rangle\geq0.
\]

The (p_{ij}) are edge weights of total mass one on (K_5), and their
weighted degrees are the eigenvalues (d_i) of (F).  If
\(R_{\rm off}\) is the squared Hilbert--Schmidt norm of the off-diagonal
part of (G) in this wedge basis, and \(\mathcal A(p)\) is the sum of
\(p_ep_f\) over distinct adjacent edges, then (1) is exactly

\[
\frac E8=\frac13-\mathcal A(p)+R_{\rm off}.\tag{2}
\]

We use the capped weighted Motzkin--Straus lemma

\[
\sum_ep_e=1,qquad d_i\leq\frac23
\quad\Longrightarrow\quad
\mathcal A(p)\leq\frac13.\tag{3}
\]

For completeness, if (d_0=2/3), split the weights into the four star
weights (a_i=p_{0i}) and the six outer weights (b_{ij}), whose total
mass is (1/3).  With outer degrees (r_i),

\[
\mathcal A(p)
=\frac29-\frac12\|a-r\|^2+\frac12\|r\|^2+\mathcal A(b)
\leq\frac29+\left(\sum b_{ij}\right)^2=\frac13.
\]

If no cap is active, take a maximizer with minimal support.  Weight can be
transferred between two disjoint positive edges along a direction on which
\(\mathcal A\) is affine, until an edge disappears or a cap becomes active.
Minimality therefore leaves a pairwise-intersecting edge family.  It is a
star or a triangle; a star violates the cap, and on a triangle the usual
three-variable inequality gives \(\mathcal A\leq1/3\).  This proves (3).

Thus (E\geq0) whenever \(\lambda_{\max}(F)\leq2/3\).

## 3. The uncapped block identity

Suppose now that (c=\lambda_{\max}(F)>2/3), and let (S\in W) be a unit
top eigenvector.  Put (V=S^\perp).  Relative to

\[
\Lambda^2W=(S\wedge V)\oplus\Lambda^2V
\]

and the identification (S\wedge V\simeq V), write

\[
G=\begin{pmatrix}A&C\\C^T&B\end{pmatrix},
\qquad R=\gamma_V(B),
\qquad \Delta=3c-2>0,
\]

and

\[
U=A-R-\frac\Delta4I_V.
\]

Since (FS=cS), the off-diagonal (S\)-to-(V) block of (F) vanishes
and (F|_V=A+R).  In dimension four, Hodge star satisfies

\[
\|R\|^2=\|B\|^2+(\operatorname{tr}B)^2
-\operatorname{tr}(B\star B\star).
\]

Substitution into (1), followed by completing the trace square, gives

\[
\boxed{
\frac E8=\frac12\|U\|^2+2\|C\|^2
+\operatorname{tr}(B\star B\star)-\frac{\Delta^2}{24}.}
\tag{4}
\]

Both (B) and \(\star B\star\) are positive semidefinite, so

\[
\operatorname{tr}(B\star B\star)\geq0.\tag{5}
\]

It remains to prove

\[
2\|U\|^2+8\|C\|^2\geq\frac{\Delta^2}{6}.\tag{6}
\]

## 4. A pointwise spectral majorant

Every line through a unit element of (W) is, after a rotation and a
permutation of coordinates, represented by

\[
S_h=\frac{\sqrt3E_0-hE_1}{\sqrt{3+h^2}},\qquad0\leq h\leq1,
\]

where

\[
E_0=\frac{\operatorname{diag}(1,-1,0)}{\sqrt2},\qquad
E_1=\frac{\operatorname{diag}(1,1,-2)}{\sqrt6}.
\]

Complete (S_h) to the orthonormal basis

\[
S_h,quad H_h=\frac{hE_0+\sqrt3E_1}{\sqrt{3+h^2}},quad
E_{xy},E_{xz},E_{yz}.
\]

For a unit vector (x=(x,y,z)), set (X=x^2,Y=y^2,Z=z^2), and define

\[
\ell_h=(3-h)X-(3+h)Y+2hZ.
\]

The three weights

\[
\begin{aligned}
w_X&=\frac{-17h^2+32h+20}{70},\\
w_Y&=\frac{40h^2-73h+40}{140},\\
w_Z&=\frac{12-5h^2}{14}
\end{aligned}\tag{7}
\]

are positive on (0\leq h\leq1).  Indeed (w_X) is concave and positive
at both endpoints,

\[
40h^2-73h+40=40\left(h-\frac{73}{80}\right)^2+\frac{1071}{160},
\]

and (12-5h^2\geq7).  Hence

\[
g_h(x)=(w_XX+w_YY+w_ZZ)\ell_h^2\geq0.\tag{8}
\]

Construct the pointwise Pluecker vector using the above orbital basis:

\[
z_{ij}(x)=2x\cdot(E_ix\times E_jx),\qquad i<j.
\]

Its squared norm is one.  Form the pointwise blocks (A_x,C_x,B_x), and
define (R_x,U_x,\Delta_x) by the same linear formulas as above.  There
are coefficient matrices (L_h\in\operatorname{Sym}(V)) and
\(M_h\in\operatorname{Hom}(\Lambda^2V,V)\) for which the following is a
polynomial identity in (x,y,z,h):

\[
\boxed{\Delta_x+g_h(x)=\langle L_h,U_x\rangle+\langle M_h,C_x\rangle.}
\tag{9}
\]

Here (L_h\) is diagonal.  Its entries, in the basis
\((H_h,E_{xy},E_{xz},E_{yz})\), are

\[
\begin{aligned}
l_0&=-\frac{194h^4-435h^3+208h^2+81h+120}{280},\\
l_1&=-\frac{206h^4+435h^3-1168h^2-81h+440}{280},\\
l_2&= \frac{274h^4-101h^3+132h^2-753h+280}{280},\\
l_3&= \frac{126h^4+101h^3-1092h^2+753h+280}{280}.
\end{aligned}\tag{10}
\]

Order the outer bivectors as

\[
(H\wedge xy,H\wedge xz,H\wedge yz,xy\wedge xz,xy\wedge yz,xz\wedge yz).
\]

The only nonzero entries of (M_h) are at

\[
(xy,0),(xy,5),(xz,1),(xz,4),(yz,2),(yz,3),
\]

and, in that order, equal

\[
\begin{aligned}
m_0={}&-\frac{\sqrt3h(12h^5+512h^4-2779h^3-3829h^2+3847h+1741)}{840(h+1)},\\
m_1={}&-\frac{\sqrt{3(h^2+3)}(12h^5+364h^4-159h^3-90h^2+717h+960)}{1260(h+1)},\\
m_2={}& \frac{\sqrt3(12h^6-32h^5-983h^4-2373h^3-2405h^2+2645h+960)}{840(h+1)},\\
m_3={}&-\frac{\sqrt{3(h^2+3)}(6h^5+887h^4+582h^3-1155h^2+288h+480)}{1260(h+1)},\\
m_4={}& \frac{\sqrt3(12h^6+908h^5-545h^4-223h^3+4749h^2-85h-960)}{840(h+1)},\\
m_5={}& \frac{\sqrt{3(h^2+3)}(6h^5-523h^4-741h^3+1065h^2+429h+480)}{1260(h+1)}.
\end{aligned}\tag{11}
\]

Identity (9) follows by substituting the ten explicit cubic Pluecker
coordinates and comparing the 28 sextic monomials.  Thus it is an exact
algebraic identity, not a numerical or optimization claim.

The weighted dual norm of these coefficients is

\[
D(h)=\frac{\|L_h\|^2}{2}+\frac{\|M_h\|^2}{8}
=\frac{N(h)}{5644800(h+1)^2},\tag{12}
\]

where

\[
\begin{aligned}
N(h)={}&1584h^{12}+117408h^{11}+10693280h^{10}+2137032h^9\\
&+13971489h^8+45236046h^7-16424313h^6+4522956h^5\\
&+68504799h^4-88091010h^3-11034711h^2+48384000h+24192000.
\end{aligned}
\]

Moreover (D(h)\leq6) on the whole interval.  After multiplying
\(6-D(h)\) by the positive denominator, its numerator has the following
degree-twelve Bernstein coefficients on ([0,1]):

\[
\begin{gathered}
9676800, 11289600, \frac{298820637}{22},\ \frac{186523506}{11},
\ \frac{3576110762}{165},\ \frac{1838154161}{66},\\
\frac{10926245933}{308},\ \frac{88743939}{2},\ \frac{2973164684}{55},
\ \frac{3483459192}{55},\ \frac{2287502816}{33},
\ \frac{194932736}{3},\ 33264640.
\end{gathered}
\]

They are all positive, proving the claim.

## 5. Completion of the proof

Integrate (9).  Linearity of all block constructions gives

\[
r:=\Delta+\int g_h\,d\mu
=\langle L_h,U\rangle+\langle M_h,C\rangle.
\]

By (8) and \(\Delta>0\), (r\geq\Delta>0).  Weighted
Cauchy--Schwarz and (12) now give

\[
r^2
\leq D(h)\bigl(2\|U\|^2+8\|C\|^2\bigr)
\leq6\bigl(2\|U\|^2+8\|C\|^2\bigr).
\]

Therefore

\[
2\|U\|^2+8\|C\|^2\geq\frac{r^2}{6}
\geq\frac{\Delta^2}{6},
\]

which is (6).  Equations (4)--(5) imply (E(\mu)\geq0) in the uncapped
case.  Together with Section 2, this proves copositivity for every
probability measure on \(\mathbb {RP}^2\).

