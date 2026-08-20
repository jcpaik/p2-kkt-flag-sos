# Exact Hankel proof on the boundary `det(I-A_2)=0`

Let \(L\) be a degree-six ternary functional, normalized by
\(L(r^6)=1\), whose middle catalecticant

\[
H_L(f,g)=L(fg),\qquad f,g\in\mathbb R[x,y,z]_3,
\]

is positive semidefinite.  Put \(\rho_x=2xx^{\mathsf T}-I\) on the
unit sphere and let

\[
A_2[S]=L(\rho_xS\rho_x),\qquad B=I-A_2
\]

on the five-dimensional space \(\operatorname{Sym}^2_0\mathbb R^3\).
Homogenizing the integrand to degree six gives, for every traceless
symmetric \(S\),

\[
\langle S,B[S]\rangle
=4L\!\left(r^2\bigl(r^2x^{\mathsf T}S^2x-(x^{\mathsf T}Sx)^2\bigr)\right)
=4\sum_{i,k}L\!\left(\bigl[x_k(x\mathbin\times Sx)_i\bigr]^2\right)\geq0. \tag{1}
\]

Thus \(B\succeq0\) on the full Hankel cone, not only for representing
measures.

Assume now that \(\det B=0\).  Choose \(0\ne S\in\ker B\).  Equality in
(1), together with \(H_L\succeq0\), implies

\[
x_k(x\mathbin\times Sx)_i\in\ker H_L\qquad(1\leq i,k\leq3). \tag{2}
\]

After an orthogonal change of coordinates, \(S\) is diagonal.

## Three distinct eigenvalues

If \(S=\operatorname{diag}(\lambda_1,\lambda_2,\lambda_3)\) has three
distinct eigenvalues, the components of \(x\times Sx\) are nonzero
multiples of \(yz,xz,xy\).  Equation (2) therefore puts every mixed
cubic

\[
x^2y,\ xy^2,\ x^2z,\ xz^2,\ y^2z,\ yz^2,\ xyz
\]

in \(\ker H_L\).  Only \(x^3,y^3,z^3\) remain.  Their cross pairings
also vanish by the Hankel identities, for example

\[
L(x^3y^3)=H_L(x^2y,xy^2)=0.
\]

Consequently

\[
L=a\,\operatorname{ev}_{(1,0,0)}
 +b\,\operatorname{ev}_{(0,1,0)}
 +c\,\operatorname{ev}_{(0,0,1)},\qquad a,b,c\geq0.
\]

After normalization this is a measure on at most three orthogonal
projective points, whose \(P_2\)-energy is nonnegative.

## One double eigenvalue

If \(S=\operatorname{diag}(\lambda,\lambda,\mu)\), \(\lambda\ne\mu\),
then the nonzero components of \(x\times Sx\) are multiples of \(xz\)
and \(yz\).  Equation (2) kills

\[
x^2z,\ xyz,\ xz^2,\ y^2z,\ yz^2.
\]

Hence \(H_L\) is block diagonal between
\(\mathbb R[x,y]_3\) and \(\mathbb Rz^3\).  The first block is a PSD
binary sextic Hankel functional and the second is a nonnegative
multiple of evaluation at the perpendicular axis.  Every nonnegative
binary sextic is a sum of squares of binary cubics; by duality, every
PSD binary middle-catalectic functional is a conic combination of real
point evaluations on \(\mathbb{RP}^1\).  Therefore \(L\) has a
representing measure supported on a projective great circle together
with its perpendicular pole.

For such a normalized measure, writing \(w\) for the pole mass and
\(\nu\) for the circle measure,

\[
E(L)=6\left(w-\frac13\right)^2
 +(1-w)^2\left(|\widehat\nu(2)|^2+|\widehat\nu(6)|^2\right)\geq0.
\]

The scalar case \(S=\lambda I\) is impossible because \(S\) is
traceless and nonzero.  This exhausts the boundary and proves

\[
\boxed{\det(I-A_2)=0\quad\Longrightarrow\quad Q(H_L)\geq0}
\]

on the entire PSD Hankel cone.

This is an exact boundary theorem.  The remaining quantitative interior
step for the proposed sharp strengthening is

\[
Q(H_L)\geq \frac{32}{105}
 \det\!\left(\frac54(I-A_2)\right)^2,
\]

whose right side vanishes precisely on the boundary just classified.
