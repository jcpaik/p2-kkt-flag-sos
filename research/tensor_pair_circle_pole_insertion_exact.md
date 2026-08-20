# Exact pole--Haar insertion audit for the projected-circle gap

Put

\[
 \mu_0={1\over3}\delta_{e_3}+{2\over3}\sigma_{e_3^\perp},\qquad
 q_z=(\sqrt z,0,\sqrt{1-z}),\qquad
 \mu_{\epsilon,z}=(1-\epsilon)\mu_0+\epsilon\delta_{q_z}.
\]

Let (F=8\mathbb E[a^2b^2(c-ab)^2\det\operatorname{Gram}(X,Y,Z)]), and let
(A_4\) be the pair-projected circle square with root weight
((x\cdot y)^4(1-(x\cdot y)^2)^6).  The exact Fourier-flow computation in
`tensor_pair_circle_pole_insertion_exact.py` gives

\[
 G(\epsilon,z):=J(\mu_{\epsilon,z})-108F(\mu_{\epsilon,z})
                 -288A_4(\mu_{\epsilon,z})
 =\epsilon H(\epsilon,z).
\]

In the tensor-product Bernstein basis
(B_i^3(\epsilon)B_j^{11}(z)), the coefficient rows of (H) are

\[
\begin{array}{c|rrrrrrrrrrrr}
i=0&0&0&27/110&269/330&509/330&529/231&4117/1386&192/55&409/110&229/66&27/11&0\\
i=1&10/3&10/3&532/165&2462/495&1901/396&31435/5544&121925/19712&928657/126720&33929/6336&72461/10560&21757/4224&1181/384\\
i=2&20/3&20/3&848/165&812/165&17/3&217/33&4293/616&853/132&22841/4224&21953/4224&29285/4224&2263/384\\
i=3&10&10&10&10&10&10&10&10&10&10&10&10.
\end{array}
\]

Every coefficient is nonnegative.  Hence (G(\epsilon,z)\ge0) exactly for
((\epsilon,z)\in[0,1]^2).  Its Gateaux coefficient at the zero face is

\[
 {\partial G\over\partial\epsilon}(0,z)
 ={z^2(1-z)\over6}
 (85z^5-271z^4+435z^3-327z^2+159z+81).
\]

The last factor has Bernstein coefficients

\[
 81,\quad {564\over5},\quad {1119\over10},\quad {609\over5},
 \quad {659\over5},\quad 162,
\]

so it is positive on ([0,1]).  In particular, (c=288) is not sharp on
this ordinary one-atom chord: the leading (A_4) ratio tends to (1152)
as (z\downarrow0).  If (288) is globally sharp, it comes from a more
degenerate/higher-order equality-face perturbation.

All entries above are produced from exact rational graph coefficients and
exact circle Fourier counts; no floating reconstruction is used.
