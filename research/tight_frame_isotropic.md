# Exact tight-frame reductions for the kernel `E`

This note records exact results for the strategy that first tries to force a
minimizer to be isotropic, and then proves the desired inequality on the cone
of probabilistic unit-norm tight frames.  It also records several tempting but
false strengthenings, so that they are not used as hidden assumptions.

## 1. Tensor-moment formulation

Put `P_x=xx^T` and

\[
 S_k=\int P_x^{\otimes k}\,d\mu(x),\qquad
 p_{2k}=\operatorname{tr}(S_k^2)
        =\iint (x\mathbin\cdot y)^{2k}\,d\mu(x)d\mu(y).
\]

If `mu` is isotropic, then `S_1=I/3`, and hence

\[
 \frac{E(\mu)}{16}
 =q(\mu):=2\|S_3\|_{\rm HS}^2-3\|S_2\|_{\rm HS}^2+\frac13.
 \tag{1}
\]

Equivalently,

\[
q=\iint s(1-s)(1-2s)\,d\mu d\mu,
\qquad s=(x\mathbin\cdot y)^2.
\tag{2}
\]

Thus the isotropic subproblem is the sharp purity inequality

\[
3\|S_2\|_{\rm HS}^2\le 2\|S_3\|_{\rm HS}^2+\|S_1\|_{\rm HS}^2.
\tag{3}
\]

## 2. Exact symmetrized trace transport

On `Sym^3(R^3)`, define

\[
B=\Pi_{\rm sym}\frac{(S_2)_{12}+(S_2)_{13}+(S_2)_{23}}3\Pi_{\rm sym}.
\]

Direct index contraction, using `Tr_2 S_2=S_1=I/3`, gives

\[
\langle S_3,B\rangle=p_4,
\qquad
\|B\|_{\rm HS}^2=\frac79p_4+\frac4{27}.
\tag{4}
\]

Indeed, before restriction to the symmetric cube, the three diagonal
contractions contribute `3*3 p_4`, while the six crossed contractions reduce
to `tr(S_1^2)=1/3`; projection to the symmetric cube gives (4).  Equivalently,
this can be verified on the irreducible `H_4` and scalar summands, where the
transport factor on `H_4` is `7/9`.

The centered Cauchy--Schwarz consequence is

\[
p_6\ge \frac17+\frac97\left(p_4-\frac15\right)
=\frac97p_4-\frac4{35}.
\tag{5}
\]

This is valid but is not strong enough near the orthonormal-basis face.

### Chebyshev/projection form

There is a particularly economical equivalent form of the isotropic target.
Since

\[
T_3(t)=4t^3-3t,\qquad
E=2\iint\big(T_3(x\cdot y)^2+(x\cdot y)^2\big)\,d\mu d\mu-\frac43,
\]

isotropy gives

\[
E\ge0\quad\Longleftrightarrow\quad
\iint T_3(x\cdot y)^2\,d\mu d\mu\ge\frac13.
\tag{5a}
\]

For an atomic measure, let `A_ij=sqrt(w_i w_j) T_3(x_i dot x_j)`.
Let `P` have kernel `3t`; isotropy says that `P` is an orthogonal projection
of rank three.  Let `B` have kernel `P_3(t)=(5t^3-3t)/2`; then `B` is positive
semidefinite, has rank at most seven and trace one.  The exact identity

\[
A=\frac85B-\frac15P
\]

turns (5a) into

\[
12\|B\|_{\rm HS}^2-3\operatorname{tr}(BP)-1\ge0.
\tag{5b}
\]

If `v_3(x)` is the normalized degree-three harmonic feature and
`u(x)=sqrt(3)x`, then, with

\[
C=\int v_3v_3^T\,d\mu,\qquad M=\int v_3u^T\,d\mu,
\]

one has `||B||^2=||C||^2` and `tr(BP)=||M||^2`.  The generic covariance
Schur complement `C-MM^T>=0` is too weak; the missing input is precisely the
coherent-cubic relation between the `H_4` part of `C` and `M`.

## 3. Reduction-map inequalities

For each pair of tensor slots,

\[
R_{12|3}:=(S_2)_{12}\otimes I-S_3
=\int P_x\otimes P_x\otimes(I-P_x)\,d\mu(x)\succeq0,
\tag{6}
\]

and similarly for the two other placements.  Since the trace of a product of
positive semidefinite matrices is nonnegative,

\[
0\le\operatorname{tr}(R_{12|3}R_{13|2})
=p_6-2p_4+\frac13.
\tag{7}
\]

Here the crossed first term is exactly `tr(S_1^2)=1/3`, and each mixed term is
`p_4`.  Inequality (7) is exact, but (1) differs from it by `p_6-p_4<=0`, so it
does not prove the target.

The scalar sixth-harmonic positivity inequality is

\[
p_6\ge\frac{15}{11}p_4-\frac{10}{77}.
\tag{8}
\]

It proves (1) whenever `p_4<=17/63`; the unresolved isotropic range under this
estimate is `17/63<p_4<=1/3`.

## 4. The signed localizer, and an exact counterexample to its positivity

Define

\[
T_\mu=\Pi_{\rm sym}\left(\frac13I-3B+2S_3\right)\Pi_{\rm sym}.
\tag{9}
\]

It is the symmetrization of

\[
\int P_x\otimes(I-P_x)\otimes(I-2P_x)\,d\mu(x),
\]

and it satisfies the exact identity

\[
q(\mu)=\operatorname{tr}(S_3T_\mu).
\tag{10}
\]

It would suffice to prove `T_mu>=0`, but that statement is false.  For the
regular tetrahedral frame

\[
x_1=(1,1,1)/\sqrt3,\quad x_2=(1,-1,-1)/\sqrt3,\quad
x_3=(-1,1,-1)/\sqrt3,\quad x_4=(-1,-1,1)/\sqrt3
\]

with weights `1/4`, the exact spectrum of `T_mu` on `Sym^3(R^3)` is

\[
\left\{-\frac19,-\frac19,-\frac19,
0,0,0,\frac19,\frac5{27},\frac5{27},\frac5{27}\right\}.
\tag{11}
\]

Nevertheless `q(mu)=14/243>0`.  Thus positivity of the localizer is strictly
stronger than the desired self-pairing and cannot be assumed.

## 5. Volume-sampling form

For three independent samples put

\[
a=X\cdot Y,\quad b=X\cdot Z,\quad c=Y\cdot Z,
\qquad D=1+2abc-a^2-b^2-c^2.
\]

Isotropy gives, conditionally on `X,Y`,

\[
\mathbb E_ZD=\frac{1-a^2}{3}.
\]

Consequently

\[
q=\mathbb E\left[D\sum_{\rm cyc}(a^2-2a^4)\right].
\tag{12}
\]

The normalization is `E D=2/9`.  Thus (12) is an expectation under
three-point volume sampling.  The bracket in (12) is not pointwise
nonnegative; for example `a=b=c=4/5` gives a positive Gram determinant and a
negative bracket.  Any proof through volume sampling must therefore use more
than pointwise Gram positivity.

The weighted cross-product law is also isotropic: if `V=X cross Y`, then

\[
\mathbb E(VV^T)=\frac29I.
\tag{13}
\]

This follows by conditioning on `X` and using
`E_Y YY^T=I/3`.  In radial form, with `r=|X cross Y|^2=1-a^2`,

\[
q=\mathbb E\big[r(1-r)(2r-1)\big].
\tag{14}
\]

Equations (13)--(14) are exact, but isotropy of `V` alone does not imply (14).

## 6. Projective-deformation Hessian trace

Let `S` range through an orthonormal basis of traceless symmetric `3 x 3`
matrices, and let

\[
V_S(x)=Sx-(x^TSx)x
\]

be the infinitesimal projective deformation of the sphere.  For `t=x dot y`,
the common deformation of both arguments satisfies

\[
\sum_S(D_St)^2=2(1-t^2)^2,
\qquad
\sum_SD_S^2t=-2t(1-t^2).
\tag{15}
\]

Therefore the trace of the common projective Hessian of a zonal kernel is

\[
\mathcal LK(t)=2(1-t^2)\big((1-t^2)K''(t)-tK'(t)\big).
\tag{16}
\]

For the present kernel this is

\[
\mathcal LK(t)=16(t-1)(t+1)(2t^2-1)(72t^4-72t^2+5)
\]

or, expanded,

\[
2304t^8-5760t^6+4768t^4-1392t^2+80.
\tag{17}
\]

Hence every minimizer stable under common projective deformations obeys the
exact necessary condition

\[
2304p_8-5760p_6+4768p_4-1392p_2+80\ge0.
\tag{18}
\]

This condition by itself does not force isotropy, but it is an additional exact
constraint available in a KKT/stability proof.

## 6a. Clebsch--Gordan block data for the symmetric reduction

The positive symmetric reduction `R=B_transport-S_3` from Section 3 becomes
especially explicit under

\[
\operatorname{Sym}^3(\mathbb R^3)=H_3\oplus r^2H_1
\quad (\dim H_3=7,\ \dim H_1=3).
\]

In orthonormal tensor coordinates it has the exact block form

\[
R=\begin{pmatrix}
\mathcal A&X\\ X^T&\frac2{45}I_3
\end{pmatrix},
\qquad \mathcal A\succeq0,
\tag{18a}
\]

with

\[
\operatorname{tr}\mathcal A=\frac8{15},\qquad
\|X\|_{\rm HS}^2=\frac{p_4}{15}-\frac1{75}
=\frac{5p_4-1}{75}.
\tag{18b}
\]

The one-slot partial trace is

\[
\operatorname{Tr}_1R=\frac4{27}\Pi_{\operatorname{Sym}^2}
-\frac29S_2,
\tag{18c}
\]

and therefore

\[
\|\operatorname{Tr}_1R\|_{\rm HS}^2
=\frac{16}{243}+\frac4{81}p_4.
\tag{18d}
\]

The desired isotropic inequality is equivalent, in these blocks, to

\[
\|\mathcal A\|_{\rm HS}^2
\ge \frac7{225}+\frac{13}{6}\|X\|_{\rm HS}^2.
\tag{18e}
\]

More invariantly, `mathcal A` has only its scalar, `H_4`, and `H_6`
multipoles (the `H_2` multipole is killed by isotropy).  With
`u=||X||_HS^2` and with `mathcal A_6` denoting the orthogonal `H_6`
component, the Clebsch--Gordan norms give

\[
\|\mathcal A\|_{\rm HS}^2
=\frac{64}{1575}+\frac4{33}u+\|\mathcal A_6\|_{\rm HS}^2.
\tag{18f}
\]

Thus the entire isotropic problem is equivalently the single sharp coherent
multipole inequality

\[
\|\mathcal A_6\|_{\rm HS}^2+\frac1{105}
\ge\frac{45}{22}u.
\tag{18g}
\]

At an orthonormal-basis measure, `u=2/225` and (18g) is an equality.

The scalar block and (18b) follow by contracting the tensor indices.  The
Schur complement from (18a),
`mathcal A >= (45/2)XX^T`, is not by itself strong enough to imply (18e): in a
generic block-matrix optimization the lower eigenvalue constraints are
inactive throughout the admissible range `p_4<=1/3`.  Thus an exact closure
must also use the `H_4/H_6` coherent-state or PPT relations, not PSD of this
single block alone.

## 7. Perron shortcut that fails

At a negative stable atomic equilibrium, the shifted weighted kernel matrix
with entries

\[
A_{ij}=\sqrt{w_iw_j}\,(K(x_i\cdot x_j)+4/3)
\]

is positive semidefinite, has nonnegative entries, has trace `4`, and has
`sqrt(w)` as its Perron vector with eigenvalue `E+4/3`.  A rank bound of three
would finish the theorem, but such a rank bound is false even on the equality
face.  A pole of weight `1/3` together with a regular projective equatorial
`n`-gon of total weight `2/3` (for example `n=5`) gives a positive semidefinite
shifted matrix of rank six.  Thus a Perron proof needs a genuinely geometric
spectral bound, not merely `rank(A)<=3`.

## 8. Tangent two-plane / fermionic reformulation

Let `rho_x=2P_x-I` and let `pi_2(rho_x)` act on the five-dimensional
Euclidean space `H_2=Sym^2_0(R^3)` by conjugation.  Set

\[
 Q_x=\frac{I-\pi_2(\rho_x)}2.
\]

Then `Q_x` is the rank-two orthogonal projection onto

\[
 W_x=\{(xv^T+vx^T)/\sqrt2:v\perp x\}\subset H_2,
\]

the tangent two-plane of the Veronese surface at `P_x-I/3`.  If
`R_x=\bigwedge^2Q_x`, viewed as the rank-one Pluecker projection on
`Lambda^2 H_2`, and

\[
 F=\int Q_x\,d\mu(x),\qquad G=\int R_x\,d\mu(x),
\]

then the normalized spin-two gap operator is

\[
 \widehat B=\frac54(I-A_2)=\frac52F. \tag{19}
\]

For two points put `s=(x dot y)^2`.  The two principal squared cosines of
`W_x,W_y` are `s` and `(2s-1)^2`.  Consequently

\[
 \operatorname{tr}(Q_xQ_y)=1-3s+4s^2,
 \qquad
 \operatorname{tr}(R_xR_y)=s(2s-1)^2. \tag{20}
\]

The kernel itself therefore has the exact fusion-frame form

\[
 \frac{K(x\cdot y)}4
 =2\operatorname{tr}(R_xR_y)-\operatorname{tr}(Q_xQ_y)+\frac23,
\]

and hence

\[
 E=8\operatorname{tr}(G^2)-4\operatorname{tr}(F^2)+\frac83. \tag{21}
\]

Thus the normalized determinant conjecture is precisely

\[
 \|G\|_{\rm HS}^2-\frac12\|F\|_{\rm HS}^2+\frac13
 \ge \frac4{105}\det\!\left(\frac52F\right)^2. \tag{22}
\]

The special tangent-Veronese constraint is essential.  Equation (22), and
even its right-hand side zero version, is false for arbitrary probability
measures on `Gr(2,5)`: take the four planes
`span(e_1,e_i)`, `i=2,3,4,5`, with equal weights.  Then

\[
 F=\operatorname{diag}(1,1/4,1/4,1/4,1/4),\quad
 \operatorname{tr}F^2=5/4,
\]

while the four Pluecker vectors are orthonormal, so `tr G^2=1/4`.
The left side of (22) is `-1/12`, equivalently (21) gives `E=-1/3`.

Under `Lambda^2(H_2)=H_1 direct-sum H_3`, the squared norm of the `H_1`
component of every oriented Pluecker vector of `W_x` is exactly `1/5`
(and that of the `H_3` component is `4/5`).  This gives the fixed block
trace constraints

\[
 \operatorname{tr}G_{11}=1/5,\qquad \operatorname{tr}G_{33}=4/5. \tag{23}
\]

At Haar measure these force
`tr G^2 >= (1/5)^2/3+(4/5)^2/7=11/105`, with equality, while
`F=(2/5)I`.  This explains the sharp Haar constant in (22); a generic Haar
measure on all of `Gr(2,5)` instead has `G=I/10` and violates (22).

### Mixed-discriminant expansion

Let `D_5` denote the polarized determinant normalized by
`D_5(A,A,A,A,A)=det A`.  Multilinearity gives the exact Cauchy--Binet
formula

\[
 \det F=\mathbb E_{x_1,\ldots,x_5}D_5(Q_{x_1},\ldots,Q_{x_5}), \tag{24}
\]

and, if `(u_{i1},u_{i2})` is any orthonormal basis of `W_{x_i}`,

\[
 D_5(Q_{x_1},\ldots,Q_{x_5})
 =\frac1{5!}\sum_{a_1,\ldots,a_5\in\{1,2\}}
 \det[u_{1a_1}\ \cdots\ u_{5a_5}]^2\ge0. \tag{25}
\]

The square in (22) is therefore an exact ten-sample expectation of the
product of two nonnegative mixed discriminants.  A naive pointwise
domination by the average of the 25 cross kernels is false: random exact
plane tuples already give a negative average cross `K` while both mixed
discriminants are strictly positive.  Any use of (24)--(25) must retain
within-tuple terms or exploit identities special to the Veronese tangent
planes.

## 9. Spectral form under isotropy and a false axial extremality

On `H_2`, isotropy gives

\[
 \langle S,S_2S\rangle=\int(x^TSx)^2d\mu
 \leq\int x^TS^2x\,d\mu=\frac13\|S\|^2.
\]

If `lambda_1,...,lambda_5` are the eigenvalues of `S_2|_{H_2}`, then

\[
 0\leq\lambda_i\leq\frac13,\qquad \sum_i\lambda_i=\frac23,
 \qquad
 \widehat B=5\left(\frac13I-S_2|_{H_2}\right). \tag{26}
\]

Consequently

\[
 \det\widehat B=5^5\prod_i\left(\frac13-\lambda_i\right),
 \qquad p_4=\frac19+\sum_i\lambda_i^2. \tag{27}
\]

This is an exact spectral reduction of the determinant, but the spectrum of
`S_2` alone does not control `p_6`; the Veronese/separable extension to
`S_3` remains essential.

A tempting further reduction to axisymmetric measures is false.  Direct
optimization over Parseval frames gives an isotropic eight-atom frame with

\[
 A_4=0.1109713027,\quad A_6=0.0019819856,
 \quad\det\widehat B=0.7261482197,
\]

whose determinant gap `q-(2/105)det(Bhat)^2` is `0.0023608456`.  The minimum
over all axisymmetric moment laws having the same `A_4` is larger,
`0.0035934015`.  Thus even though axial families pass the conjecture, they
are not extremal at fixed fourth-harmonic norm.  A reproducible search and
the frame coordinates are in `research/isotropic_axial_envelope_opt.py`.
