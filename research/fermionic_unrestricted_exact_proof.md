# Exact unrestricted fermionic proof of copositivity

This note closes the unrestricted gap left in
`fermionic_motzkin_unrestricted.md`.  It proves the original kernel energy

\[
 E(\mu)=\iint\left(32(x^Ty)^6-48(x^Ty)^4+20(x^Ty)^2-\frac43\right)
 \,d\mu(x)d\mu(y)
\]

is nonnegative for every probability measure on \(\mathbb {RP}^2\).
All coefficients below are exact.  The finite polynomial calculation is
verified by `fermionic_unrestricted_qmodule_verify.py`.

## 1. Fermionic reduction

Use the tangent-plane operators \(F,G\) from
`tensor_fermionic_isotropic_proof.md`.  Thus

\[
 \frac{E(\mu)}8=\|G\|^2-\frac12\|F\|^2+\frac13. \tag{1}
\]

If every eigenvalue of \(F\) is at most \(2/3\), the capped
Motzkin--Straus lemma proves (1) is nonnegative.  It remains to consider a
unit top eigenvector \(S\in H_2\) with

\[
 c=\langle S,FS\rangle>\frac23,
 \qquad \Delta=3c-2>0. \tag{2}
\]

Split

\[
 \bigwedge^2H_2=(S\wedge S^\perp)\oplus\bigwedge^2S^\perp,
 \qquad G=\begin{pmatrix}A&C\\C^T&B\end{pmatrix}.
\]

Let \(R\) be the one-particle contraction of \(B\), let \(*\) be the
four-dimensional Hodge involution on \(\bigwedge^2S^\perp\), and put

\[
 U=A-R-\frac{\Delta}{4}I_4.
\]

The exact invariant identity from the unrestricted reduction is

\[
 \frac{E(\mu)}8
 =\frac12\|U\|^2+2\|C\|^2+\operatorname{tr}(B*B*)
 -\frac{\Delta^2}{24}. \tag{3}
\]

Because \(B\succeq0\) and \(*B*\succeq0\), the Hodge term is
nonnegative.  It is therefore enough to prove

\[
 2\|U\|^2+8\|C\|^2\ge\frac{\Delta^2}{6}. \tag{4}
\]

The next two sections prove the stronger estimate (4) directly from the
tangent-Veronese structure.

## 2. Spectral chart and a nonnegative correction

After an orthogonal change of coordinates and, if necessary, replacing
\(S\) by \(-S\), every unit \(S\) has the form

\[
 S=\frac{\sqrt3E_0-hE_1}{\sqrt{3+h^2}},\qquad 0\le h\le1, \tag{5}
\]

where

\[
 E_0=\frac{\operatorname{diag}(1,-1,0)}{\sqrt2},\qquad
 E_1=\frac{\operatorname{diag}(1,1,-2)}{\sqrt6}.
\]

For completeness, choose the sign so that the eigenvalues of \(S\) are
\(a,b,-a-b\) with \(a\ge b\ge0\), and permute the axes to put the negative
eigenvalue in the middle slot.  If \(r=b/a\in[0,1]\), then the ratio of
the two nonnegative eigenvalues in (5) is \(2h/(3-h)\).  Thus
\(h=3r/(2+r)\in[0,1]\), and trace zero plus unit norm fixes the common
scale.  This proves that (5) covers every spectral type, not only the two
endpoint orbits.

Complete it to the orthonormal basis

\[
 (S,H,E_{xy},E_{xz},E_{yz}),\qquad
 H=\frac{hE_0+\sqrt3E_1}{\sqrt{3+h^2}}. \tag{6}
\]

For a unit root \(x=(x,y,z)\), write \(X=x^2,Y=y^2,Z=z^2\), and set

\[
 \ell_h=(3-h)X-(3+h)Y+2hZ. \tag{7}
\]

In fact \(\ell_h=\sqrt{6(3+h^2)}\,x^TSx\).  Define

\[
 g_h(x)=(w_XX+w_YY+w_ZZ)\ell_h^2, \tag{8}
\]

with

\[
 \begin{aligned}
 w_X&=\frac{-17h^2+32h+20}{70},\\
 w_Y&=\frac{40h^2-73h+40}{140},\\
 w_Z&=\frac{12-5h^2}{14}.
 \end{aligned} \tag{9}
\]

All three weights are positive on \([0,1]\): \(w_X\) is concave and
positive at both endpoints, \(w_Y\) has minimum
\(1071/22400\) at \(h=73/80\), and \(w_Z\ge1/2\).  Consequently

\[
 g_h(x)\ge0. \tag{10}
\]

## 3. Exact tangent identity and its cost

Use the outer-edge order

\[
 (H\wedge E_{xy},H\wedge E_{xz},H\wedge E_{yz},
 E_{xy}\wedge E_{xz},E_{xy}\wedge E_{yz},E_{xz}\wedge E_{yz})
\]

for the columns of \(C\), and index its rows by
\((H,E_{xy},E_{xz},E_{yz})\).  There are exact coefficient vectors
\(l=(l_0,\ldots,l_3)\) and \(m=(m_0,\ldots,m_5)\) such that

\[
 \boxed{
 \Delta+\int g_h\,d\mu
 =\sum_{i=0}^3l_iU_{ii}
 +m_0C_{1,0}+m_1C_{1,5}+m_2C_{2,1}
 +m_3C_{2,4}+m_4C_{3,2}+m_5C_{3,3}.} \tag{11}
\]

The four \(l_i\) are

\[
\begin{aligned}
l_0={}&-\frac{194h^4-435h^3+208h^2+81h+120}{280},\\
l_1={}&-\frac{206h^4+435h^3-1168h^2-81h+440}{280},\\
l_2={}& \frac{274h^4-101h^3+132h^2-753h+280}{280},\\
l_3={}& \frac{126h^4+101h^3-1092h^2+753h+280}{280}.
\end{aligned} \tag{12}
\]

The six \(m_i\) are

\[
\begin{aligned}
m_0={}&-\frac{\sqrt3h(12h^5+512h^4-2779h^3-3829h^2+3847h+1741)}{840(h+1)},\\
m_1={}&-\frac{\sqrt{3(3+h^2)}(12h^5+364h^4-159h^3-90h^2+717h+960)}{1260(h+1)},\\
m_2={}& \frac{\sqrt3(12h^6-32h^5-983h^4-2373h^3-2405h^2+2645h+960)}{840(h+1)},\\
m_3={}&-\frac{\sqrt{3(3+h^2)}(6h^5+887h^4+582h^3-1155h^2+288h+480)}{1260(h+1)},\\
m_4={}& \frac{\sqrt3(12h^6+908h^5-545h^4-223h^3+4749h^2-85h-960)}{840(h+1)},\\
m_5={}& \frac{\sqrt{3(3+h^2)}(6h^5-523h^4-741h^3+1065h^2+429h+480)}{1260(h+1)}.
\end{aligned} \tag{13}
\]

Identity (11) is just a homogeneous sextic polynomial identity before
integration.  With

\[
 z_{ij}(x)=2x\cdot(E_ix\times E_jx),
\]

substitute \(G_x=z(x)z(x)^T\), form \(U_x,C_x,\Delta_x\), and expand
both sides.  Every coefficient agrees exactly; the symbolic verifier does
this expansion independently.

Here \(\Delta_x=3\operatorname{tr}A_x-2\operatorname{tr}G_x\), so
\(G\mapsto(\Delta,U,C)\) is linear.  For unit roots
\(\operatorname{tr}G_x=1\), and integrating \(\Delta_x\) gives exactly
\(3c-2=\Delta\).  Thus no symmetry averaging or hidden isotropy
assumption is used in passing from the pointwise identity to (11).

The weighted squared dual norm of (11) is

\[
 D(h)=\frac12\sum_{i=0}^3l_i^2+\frac18\sum_{i=0}^5m_i^2
 =\frac{P(h)}{5644800(h+1)^2}, \tag{14}
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

Moreover

\[
 6-D(h)=\frac{R(h)}{5644800(h+1)^2}, \tag{15}
\]

and the degree-twelve Bernstein coefficients of \(R\) on \([0,1]\) are

\[
\begin{gathered}
9676800,\ 11289600,\ \frac{298820637}{22},\ \frac{186523506}{11},
\ \frac{3576110762}{165},\ \frac{1838154161}{66},\\
\frac{10926245933}{308},\ \frac{88743939}{2},
\ \frac{2973164684}{55},\ \frac{3483459192}{55},
\ \frac{2287502816}{33},\ \frac{194932736}{3},\ 33264640.
\end{gathered} \tag{16}
\]

They are all strictly positive, so

\[
 0<D(h)<6\qquad(0\le h\le1). \tag{17}
\]

## 4. Closing the uncapped branch

Put

\[
 r=\Delta+\int g_h\,d\mu.
\]

By (2) and (10), \(r\ge\Delta>0\).  Weighted Cauchy applied to (11)
and then inclusion of the omitted entries gives

\[
\begin{aligned}
r^2
&\le D(h)\left(
2\sum_{i=0}^3U_{ii}^2
+8\sum_{(i,j)\in\mathcal C}C_{ij}^2\right)\\
&\le D(h)\left(2\|U\|^2+8\|C\|^2\right),
\end{aligned} \tag{18}
\]

where \(\mathcal C\) is the six-entry set in (11).  Equations (17)--(18)
give

\[
 2\|U\|^2+8\|C\|^2
 \ge\frac{r^2}{D(h)}
 \ge\frac{\Delta^2}{6}. \tag{19}
\]

Thus (4) holds.  Substitution into (3) yields \(E(\mu)\ge0\) in the
uncapped branch.  The capped Motzkin--Straus lemma handles the other
branch, so the kernel is copositive for every probability measure on
\(\mathbb {RP}^2\).

The argument is exact and includes all spectral types, including the
rank-two case \(h=0\) and the axial double-eigenvalue case \(h=1\).
