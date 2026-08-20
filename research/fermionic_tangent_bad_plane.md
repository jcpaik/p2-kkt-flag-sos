# A sharp pointwise exclusion of the fermionic four-star

Let \(S\in\operatorname{Sym}^2_0(\mathbb R^3)\) have Frobenius norm one and
let

\[
 q_S(x)=\langle S,Q_xS\rangle
       =2\bigl(|Sx|^2-(x^TSx)^2\bigr).
\]

For orthonormal \(E,F\in H_2\), the normalized tangent Pluecker coordinate
is

\[
 z_{E,F}(x)=2x\mathbin\cdot(Ex\mathbin\times Fx).
\]

Diagonalize \(S\), ordering its eigenvalues as \(a\ge m\ge b\), and call
the corresponding coordinate directions \(e_1,e_3,e_2\), respectively
(so \(e_1,e_2\) are the two extreme eigendirections).  Inside \(S^\perp\),
let \(B_S\) be the two-plane spanned by

* the unit diagonal traceless matrix \(H\perp S\), and
* \(E_{12}=(e_1e_2^T+e_2e_1^T)/\sqrt2\).

Then the following sharp pointwise inequality holds:

\[
 \boxed{\quad z_{S,H}(x)^2+z_{S,E_{12}}(x)^2
       \le 4\bigl(1-q_S(x)\bigr).\quad} \tag{1}
\]

The constant four is sharp in the limit in which \(S\) has spectrum
\((1/\sqrt2,0,-1/\sqrt2)\) and \(x\) tends to either of the two roots for
which \(S\in W_x\).

## Exact coordinate proof

Write

\[
 a={-m+d\over2},\qquad b={-m-d\over2},\qquad d^2=2-3m^2.
\]

The ordering \(a\ge m\ge b\) says \(d\ge3|m|\).  Put

\[
 h={3m\over d}\in[-1,1],\qquad d^2={6\over3+h^2}.
\]

For \(X=x_1^2,Y=x_2^2,Z=x_3^2\), set

\[
 r=Z,\qquad t={X-Y\over X+Y}\in[-1,1]
\]

(with arbitrary \(t\) if \(X+Y=0\)).  Direct expansion gives

\[
 q_S={3\over3+h^2}\left[
 (1-r)^2(1-t^2)+r(1-r)(1+h^2-2ht)\right]
\]

and

\[
 z_{S,H}^2+z_{S,E_{12}}^2
 =3r(1-r)^2\left[(1-t^2)+{(1-ht)^2\over3+h^2}\right].
\]

Consequently
\[
 (3+h^2)\left[4(1-q_S)-
   \bigl(z_{S,H}^2+z_{S,E_{12}}^2\bigr)\right]=P(r,t,h),
\]
where

\[
\begin{aligned}
 P(r,t,h)={}&12r^2(2-r)
 +(4-15r+18r^2-3r^3)h^2\\
 &+(12-15r-6r^2+9r^3)t^2
 +(6r^3-36r^2+30r)ht. \tag{2}
\end{aligned}
\]

Here is a short exact check that \(P\ge0\) on the cube.  For fixed \(r\),
the trace of the quadratic-form matrix in \((h,t)\) is

\[
 2(3r^3+6r^2-15r+8)>0\qquad(0\le r\le1).
\]

For completeness, the cubic inside the parentheses has its interior
minimum at \(r=(-2+\sqrt{19})/3\), where its value is
\((178-38\sqrt{19})/9>0\).

Thus that matrix is either positive semidefinite, in which case (2) is
immediate, or has exactly one negative eigenvalue, in which case a minimum
on the square \([-1,1]^2\) occurs on its boundary.  The symmetry
\(P(r,-t,-h)=P(r,t,h)\) leaves only \(h=1\) and \(t=1\).

For \(h=1\), \(P\) is a quadratic in \(t\) with leading coefficient

\[
 3(1-r)^2(3r+4)
\]

and discriminant

\[
 192(1-r)^2(3r^4-6r^3-3r^2+3r-1)\le0.
\]

The last inequality is immediate from the degree-four Bernstein
coefficients of the negative quartic,

\[
 (1,\tfrac14,0,\tfrac74,4).
\]

For \(t=1\), the coefficient of \(h^2\) is

\[
 A(r)=-3r^3+18r^2-15r+4>0.
\]

Its interior minimum is \(22-14\sqrt{21}/3>0\).

Its discriminant is \(192g(r)\), where

\[
 g(r)=r^3-6r^2+5r-1.
\]

If \(g\le0\), the quadratic is nonnegative.  If \(g>0\), its vertex lies
to the left of \([-1,1]\), because its linear coefficient \(B\) satisfies

\[
 B-2A=4(3g+1)>0.
\]

The minimum on the interval is therefore at \(h=-1\), where

\[
 P(r,1,-1)=4(1-3g(r))>0.
\]

Indeed the exact maximum of \(g\) on \([0,1]\) is
\(-7+14\sqrt{21}/9<1/3\).  This proves (1).

## Integrated consequence

For a probability measure \(\mu\), fix \(S\) and let

\[
 V_S=\int v_S(x)v_S(x)^T\,d\mu(x)
\]

be the four-by-four star block of \(G=\int R_xd\mu(x)\), relative to
\(S\wedge S^\perp\).  Its trace is
\(c_S=\langle S,FS\rangle=\int q_Sd\mu\).  Integrating (1) yields the
exact linear constraint

\[
 \boxed{\quad \operatorname{tr}(P_{B_S}V_S)\le4(1-c_S).\quad} \tag{3}
\]

In particular, if \(c_S\ge4/5\), then

\[
 \operatorname{tr}(P_{B_S^\perp}V_S)\ge5c_S-4,
 \qquad
 \|V_S\|_F^2\ge{(5c_S-4)^2\over2}. \tag{4}
\]

Equation (3) is a quantitative tangent-Veronese obstruction to the
uniform four-edge star.  It is not a consequence of the generic
Grassmann/fermionic constraints.
