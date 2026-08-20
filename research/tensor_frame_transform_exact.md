# The weighted ONB frame transform

Let \(\mu\) be a probability measure on \(\mathbb {RP}^{2}\).  For an
ordered nonparallel pair \((x,y)\), put

\[
 a=x\cdot y,\qquad s=1-a^{2},\qquad
 u={y-ax\over\sqrt{s}},\qquad n={x\times y\over\sqrt{s}},
\]

and let

\[
 R_{xy}=(x,u,n),\qquad
 \nu_{R_{xy}}={1\over3}(\delta_x+\delta_u+\delta_n),\qquad
 w_{xy}={a^{2}s^{2}\over4}.
\]

The definition is projectively well defined; choices of signs and orientation
do not change \(\nu_R\).  Define the positive (usually subprobability)
measure

\[
 \mathcal T\mu=24\iint w_{xy}\nu_{R_{xy}}\,d\mu(x)d\mu(y).       \tag{1}
\]

Every \(\nu_R\) is isotropic, hence so is \(\mathcal T\mu\).  Its mass is

\[
 \tau=\mathcal T\mu(1)
 =6\bigl(p_2-2p_4+p_6\bigr)\le {8\over9}.              \tag{2}
\]

The upper bound follows from
\(a^2(1-a^2)^2\le4/27\).  Notice that \(\mu-\mathcal T\mu\) need not be a
positive measure.  Only the first-axis part of (1) is dominated by \(\mu\):
its Radon--Nikodym density is
\(8\int w_{xy}\,d\mu(y)\le8/27\).  The other two axes are new.

## Exact cross identity

Write

\[
 J(t)=144t^6-216t^4+87t^2-5.
\]

For any orthonormal frame \(R=(r_1,r_2,r_3)\), isotropy of \(\nu_R\) and
the ONB identity give

\[
 I_J(\mu,\nu_R)
 =144\int\prod_{i=1}^3(r_i\cdot z)^2\,d\mu(z).          \tag{3}
\]

Since

\[
 w_{xy}\prod_{r\in R_{xy}}(r\cdot z)^2
 ={a^2\over4}(x\cdot z)^2
       (y\cdot z-a,x\cdot z)^2\det(x,y,z)^2,
\]

(1)--(3) imply the denominator-free potential formula

\[
 U_{J,\mathcal T\mu}(z)
 =864\iint a^2(x\cdot z)^2(y\cdot z-a,x\cdot z)^2
          \det(x,y,z)^2\,d\mu(x)d\mu(y).               \tag{4}
\]

Consequently, for

\[
 \mathcal F=8\iiint a^2b^2(c-ab)^2
 \det\operatorname{Gram}(x,y,z)\,d\mu^3,
\]

one has the exact identity

\[
 \boxed{I_J(\mu,\mathcal T\mu)=108\mathcal F}.          \tag{5}
\]

Thus the proposed four-sample strengthening is

\[
 J(\mu)-108\mathcal F-288\mathcal A_4
 =I_J(\mu,\mu-\mathcal T\mu)-288\mathcal A_4.          \tag{6}
\]

## Harmonic moments

The Legendre expansion is

\[
 J={48\over35}P_0+{22\over7}P_2-{1728\over385}P_4
       +{768\over77}P_6.                                \tag{7}
\]

For any consistently normalized harmonic evaluation vector \(Y_\ell\),

\[
 m_\ell(\mathcal T\mu)
 =8\iint w_{xy}\sum_{r\in R_{xy}}Y_\ell(r)\,d\mu(x)d\mu(y).
                                                               \tag{8}
\]

In particular \(m_2(\mathcal T\mu)=0\).  The zonal degree-four and
degree-six moments, which are sometimes more convenient than (8), are

\[
 \int P_4(r\cdot z)d\mathcal T\mu(r)
 ={7\over4}\iint a^2s^2
 \left(5\sum_{i=1}^3(r_i\cdot z)^4-3\right)d\mu^2,      \tag{9}
\]

and

\[
 \int P_6(r\cdot z)d\mathcal T\mu(r)
 ={1\over8}\iint a^2s^2
 \left(231\sum_{i=1}^3(r_i\cdot z)^6
 -315\sum_{i=1}^3(r_i\cdot z)^4+90\right)d\mu^2.       \tag{10}
\]

Formula (4) shows that the apparent denominators in (8)--(10) cancel after
the ordered-pair symmetrization and convolution with \(J\).

## Why a two-by-two Schur argument loses the circle term

The projected-circle term is not a Schur residual of the \(J\)-Gram matrix
on \(\operatorname{span}\{\mu,\mathcal T\mu\}\).  Indeed, take

\[
 \mu={1\over2}(\delta_x+\delta_y),\qquad
 q=(x\cdot y)^2\in(0,1).
\]

Both atoms lie on a great circle.  The only nonzero ordered-pair terms in
(1) are the two ONBs \(R_{xy}\) and \(R_{yx}\), and those ONBs share the
normal to the circle.  Hence \(\mathcal T\mu\) is a pole--equator measure
with zero \(J\)-energy.  Moreover its ONB potential vanishes on the original
circle, so

\[
 J(\mathcal T\mu)=I_J(\mu,\mathcal T\mu)=0.             \tag{11}
\]

On the other hand, direct circle Fourier calculation gives

\[
 \mathcal A_4
 =q^3(1-q)^6(8q^2-12q+5)>0.                            \tag{12}
\]

Thus the two-by-two Gram determinant is zero while the conditional circle
square is strictly positive.  Any successful Schur/iteration proof of (6)
must retain the latent plane (or its degree-two and degree-six circle
features); averaging it away into \(\mathcal T\mu\) irreversibly loses the
quantity in (12).

## Exact counterexample to the proposed coefficient 288

In fact, (6) itself is false.  Let

\[
 \mu_0={1\over3}\delta_{e_3}
 +{3\over100}(\delta_{e_1}+\delta_{e_2})
 +{91\over300}(\delta_{(3,4,0)/5}+\delta_{(-4,3,0)/5}).
\]

This is a convex mixture of two ONB measures sharing the pole.  Put

\[
 q={1\over901}(36,48,899),\qquad
 \mu=(1-10^{-8})\mu_0+10^{-8}\delta_q.
\]

The exact rational verifier `tensor_pair_circle_counterexample_exact.py`
gives

\[
 J=1.170211576571\ldots10^{-12},\quad
 \mathcal F=9.567977810729\ldots10^{-15},\quad
 \mathcal A_4=6.465091976191\ldots10^{-16},
\]

and

\[
 \boxed{J-108\mathcal F-288\mathcal A_4
 =-4.932467590170\ldots10^{-14}<0}.
\]

All coordinates and weights are rational and the sign is checked using
exact `Fraction` arithmetic.  Thus the coefficient 288 cannot be used in a
universal certificate; this is a first-order escape from a non-Haar
pole--equator zero face.
