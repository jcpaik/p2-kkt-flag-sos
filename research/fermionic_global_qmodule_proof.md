# Exact unrestricted fermionic certificate

This note closes the uncapped branch of the fermionic--Motzkin reduction
for every probability measure on (\mathbb{RP}^2).  Combined with the
capped weighted-(K_5) lemma in `fermionic_cap_purity.md`, it proves the
copositivity functional exactly:

\[
E(\mu)=32p_6-48p_4+20p_2-\frac43\ge0.
\]

## 1. Fermionic split

Let (W=\operatorname{Sym}^2_0(\mathbb R^3)).  The tangent plane at
(x\in S^2) gives a unit Pluecker vector (z_x\in\Lambda^2W), and set

\[
G=\int |z_x\rangle\langle z_x|\,d\mu(x),
\qquad F=2\operatorname{Tr}_2G.
\]

Then (G\succeq0), (\operatorname{tr}G=1),
(\operatorname{tr}F=2), and the exact overlap calculation gives

\[
\frac{E(\mu)}8=\|G\|_F^2-\frac12\|F\|_F^2+\frac13. \tag{1}
\]

Let (c=\lambda_{\max}(F)), and choose a unit top eigenvector
(S\in W).  Relative to

\[
\Lambda^2W=(S\wedge S^\perp)\oplus\Lambda^2S^\perp,
\qquad
G=\begin{pmatrix}A&C\\C^T&B\end{pmatrix},
\]

put

\[
\delta=3c-2,
\qquad
U=A-\gamma(B)-\frac\delta4I_4, \tag{2}
\]

where (\gamma(B)) is the one-particle contraction of (B).  The exact
Hodge identity on (\Lambda^2S^\perp\simeq\Lambda^2\mathbb R^4) yields

\[
\frac{E(\mu)}8
=\frac12\|U\|_F^2+2\|C\|_F^2
+\operatorname{tr}(B*B*)-\frac{\delta^2}{24}. \tag{3}
\]

Since (B\succeq0) and (*B*\succeq0), their trace pairing is
nonnegative.  Therefore the uncapped case (c>2/3) follows from

\[
2\|U\|_F^2+8\|C\|_F^2\ge\frac{\delta^2}{6}. \tag{4}
\]

If (c\le2/3), every eigenvalue of (F) is at most (2/3), and the
capped weighted-(K_5) lemma proves (1) directly.  Thus (4) is the only
remaining case.

## 2. Spectral chart for the top orbital

After a spatial rotation, a sign choice for (S), and a permutation of
the axes, every line spanned by a unit traceless symmetric matrix has the
form

\[
S_h=\frac{\sqrt3S_0-hH_0}{\sqrt{3+h^2}}
=\frac{\operatorname{diag}(3-h,-3-h,2h)}
       {\sqrt{6(3+h^2)}},
\qquad 0\le h\le1, \tag{5}
\]

where

\[
S_0=\frac{\operatorname{diag}(1,-1,0)}{\sqrt2},
\qquad
H_0=\frac{\operatorname{diag}(1,1,-2)}{\sqrt6}.
\]

Indeed choose the sign so that (S) has two nonnegative eigenvalues,
order them (a\ge b\ge0), and put (r=b/a\in[0,1]).  The ratio in (5)
is (2h/(3-h)), so (h=3r/(2+r)\in[0,1]).

Complete (S_h) to the orthonormal orbital basis

\[
S_h,quad
H_h=\frac{hS_0+\sqrt3H_0}{\sqrt{3+h^2}},quad
E_{xy},E_{xz},E_{yz}. \tag{6}
\]

All block coordinates below refer to this basis.

## 3. A nonnegative pointwise correction

Write (X=x^2,Y=y^2,Z=z^2), and define

\[
\ell_h=(3-h)X-(3+h)Y+2hZ. \tag{7}
\]

Thus
(\ell_h=\sqrt{6(3+h^2)}\,x^TS_hx).  Define

\[
g_h(x)=\bigl(w_X(h)X+w_Y(h)Y+w_Z(h)Z\bigr)\ell_h^2, \tag{8}
\]

with

\[
\begin{aligned}
w_X&=\frac{-17h^2+32h+20}{70},\\
w_Y&=\frac{40h^2-73h+40}{140},\\
w_Z&=\frac{12-5h^2}{14}.
\end{aligned} \tag{9}
\]

Their degree-two Bernstein coefficients on ([0,1]) are respectively

\[
\left(\frac27,\frac{18}{35},\frac12\right),\qquad
\left(\frac27,\frac1{40},\frac1{20}\right),\qquad
\left(\frac67,\frac67,\frac12\right). \tag{10}
\]

Hence every weight is strictly positive and (g_h(x)\ge0).

For a pure tangent projector (G_x=|z_x\rangle\langle z_x|), form
(U_x,C_x,\delta_x) linearly as in (2).  There is an exact pointwise
homogeneous-sextic identity

\[
\boxed{\quad
\delta_x+g_h(x)=\Lambda_h(U_x,C_x).
\quad} \tag{11}
\]

Here

\[
\Lambda_h(U,C)
=\sum_{i=0}^3\lambda_iU_{ii}
+\eta_0C_{1,0}+\eta_1C_{1,5}
+\eta_2C_{2,1}+\eta_3C_{2,4}
+\eta_4C_{3,2}+\eta_5C_{3,3}, \tag{12}
\]

where the (C)-indices are zero-based inside the (4\times6) block.
The exact coefficient functions are recorded in the Appendix below.
Expanding both sides of (11) in the ten even sextic monomials proves the
identity coefficient by coefficient.  The symbolic audit
`fermionic_global_qmodule_verify.py` performs precisely this expansion.

The coefficients were chosen so that their weighted dual norm is

\[
D(h):=\frac12\sum_{i=0}^3\lambda_i^2
+\frac18\sum_{j=0}^5\eta_j^2
=\frac{N(h)}{5644800(1+h)^2}, \tag{13}
\]

where

\[
\begin{aligned}
N(h)={}&1584h^{12}+117408h^{11}+10693280h^{10}+2137032h^9\\
&+13971489h^8+45236046h^7-16424313h^6+4522956h^5\\
&+68504799h^4-88091010h^3-11034711h^2\\
&+48384000h+24192000.
\end{aligned} \tag{14}
\]

Moreover (D(h)<6) throughout ([0,1]).  Indeed, if

\[
R(h)=5644800(1+h)^2(6-D(h)),
\]

then its degree-twelve Bernstein coefficients are

\[
\begin{gathered}
9676800, 11289600, \frac{298820637}{22},
\frac{186523506}{11},
\frac{3576110762}{165},
\frac{1838154161}{66},\\
\frac{10926245933}{308},
\frac{88743939}{2},
\frac{2973164684}{55},
\frac{3483459192}{55},
\frac{2287502816}{33},
\frac{194932736}{3},
33264640.
\end{gathered} \tag{15}
\]

Every entry is positive, so (R(h)>0) on ([0,1]).

## 4. Integration and weighted Cauchy

The maps (G\mapsto\delta,U,C) are linear.  Integrating (11) therefore
gives

\[
r:=\delta+\int g_h(x)\,d\mu(x)=\Lambda_h(U,C). \tag{16}
\]

In the uncapped branch (\delta>0), and (8)--(10) give (r\ge\delta>0).
Weighted Cauchy--Schwarz, (12), and (13) imply

\[
r^2\le D(h)\bigl(2\|U\|_F^2+8\|C\|_F^2\bigr)
\le6\bigl(2\|U\|_F^2+8\|C\|_F^2\bigr). \tag{17}
\]

Consequently (4) holds.  Finally, rewriting (3) as

\[
\frac{E(\mu)}8
=\frac14\left(2\|U\|_F^2+8\|C\|_F^2
-\frac{\delta^2}{6}\right)
+\operatorname{tr}(B*B*)
\]

proves (E(\mu)\ge0) in the uncapped branch.  Together with the capped
weighted-(K_5) lemma, this proves copositivity for every probability
measure, with no isotropy or finite-support assumption.

## Appendix: coefficients in the pointwise identity

\[
\begin{aligned}
\lambda_0&=-\frac{194h^4-435h^3+208h^2+81h+120}{280},\\
\lambda_1&=-\frac{206h^4+435h^3-1168h^2-81h+440}{280},\\
\lambda_2&=\frac{274h^4-101h^3+132h^2-753h+280}{280},\\
\lambda_3&=\frac{126h^4+101h^3-1092h^2+753h+280}{280}.
\end{aligned}
\]

Writing (d=\sqrt{h^2+3}),

\[
\begin{aligned}
\eta_0={}&-\frac{\sqrt3h}{840(h+1)}
(12h^5+512h^4-2779h^3-3829h^2+3847h+1741),\\
\eta_1={}&-\frac{\sqrt3d}{1260(h+1)}
(12h^5+364h^4-159h^3-90h^2+717h+960),\\
\eta_2={}&\frac{\sqrt3}{840(h+1)}
(12h^6-32h^5-983h^4-2373h^3-2405h^2+2645h+960),\\
\eta_3={}&-\frac{\sqrt3d}{1260(h+1)}
(6h^5+887h^4+582h^3-1155h^2+288h+480),\\
\eta_4={}&\frac{\sqrt3}{840(h+1)}
(12h^6+908h^5-545h^4-223h^3+4749h^2-85h-960),\\
\eta_5={}&\frac{\sqrt3d}{1260(h+1)}
(6h^5-523h^4-741h^3+1065h^2+429h+480).
\end{aligned}
\]
