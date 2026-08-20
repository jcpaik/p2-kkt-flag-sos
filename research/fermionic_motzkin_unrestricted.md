# Unrestricted fermionic--Motzkin reduction

> **Resolved.**  The tangent estimate isolated below is proved for every
> spectral type in `fermionic_unrestricted_exact_proof.md`, with the exact
> symbolic certificate checked by `fermionic_unrestricted_qmodule_verify.py`.

This note starts from `tensor_fermionic_isotropic_proof.md` and isolates the
exact extra lemma needed to remove isotropy. It is a reduction, not yet a
proof of the final tangent lemma.

Let (W=\operatorname{Sym}^2_0(\mathbb R^3)), let

\[
G=\int R_x\,d\mu(x)\quad\text{on }\bigwedge^2W,
\qquad F=2\operatorname{Tr}_2G=\int Q_x\,d\mu(x).
\]

As in the isotropic proof,

\[
\frac{E(\mu)}8=\|G\|_{HS}^2-\frac12\|F\|_{HS}^2+\frac13. \tag{1}
\]

## 1. The pinched Motzkin identity

Choose an orthonormal eigenbasis (e_0,\ldots,e_4) of (F). In the wedge
basis (e_i\wedge e_j), write

\[
p_{ij}=G_{ij,ij},\qquad
d_i=\sum_{j\ne i}p_{ij},\qquad
R_{\rm off}=\sum_{(i,j)\ne(k,l)}|G_{ij,kl}|^2.
\]

Then the (d_i) are the eigenvalues of (F), 
(\sum p_{ij}=1), and, if

\[
A(p)=\sum_{e<f,\ e\cap f\ne\varnothing}p_ep_f,
\]

the exact pinching identity is

\[
\frac{E(\mu)}8=\frac13-A(p)+R_{\rm off}. \tag{2}
\]

When every (d_i\le2/3), the capped Motzkin--Straus lemma from the
isotropic proof gives (A(p)\le1/3). Thus only a top occupation
(c=d_0>2/3) remains.

Put

\[
a_i=p_{0i}\quad(1\le i\le4),\qquad
b_{ij}=p_{ij}\quad(1\le i<j\le4),
\]

let (r_i=\sum_{j\ne i}b_{ij}), let

\[
D_b=b_{12}b_{34}+b_{13}b_{24}+b_{14}b_{23},
\]

and set

\[
\Delta=3c-2,\qquad
u=a-r,\qquad
u_0=u-\frac{\Delta}{4}{\bf1}.
\]

Splitting adjacent pairs into star--star, star--outer, and outer--outer
pairs gives

\[
A(p)=\frac{c^2}{2}+\|r\|^2-\sum b_{ij}^2
      -\frac12\|a-r\|^2.
\]

Since

\[
\|r\|^2-\sum b_{ij}^2=(1-c)^2-2D_b,
\]

and (\|u\|^2=\Delta^2/4+\|u_0\|^2), this becomes

\[
\boxed{
\frac{E(\mu)}8
=R_{\rm off}+2D_b+\frac12\|u_0\|^2-\frac{\Delta^2}{24}.} \tag{3}
\]

Consequently the unrestricted theorem follows from the quantitative
tangent-star lemma

\[
(3c-2)^2\le
24R_{\rm off}+48D_b+12\|u_0\|^2. \tag{4}
\]

## 2. Basis-free matrix form of the same remainder

There is a cleaner invariant form of (3). Fix only the top eigenvector
(S=e_0\) of (F), and split

\[
\bigwedge^2W=(S\wedge S^\perp)\oplus\bigwedge^2S^\perp,
\qquad
G=\begin{pmatrix}A&C\\ C^T&B\end{pmatrix}.
\]

Thus (A) is (4\times4), (B) is (6\times6), and
(c=\operatorname{tr}A). Let (R=\operatorname{Tr}_1B) be the
one-particle contraction of (B) on (S^\perp); with this normalization,
(\operatorname{tr}R=2\operatorname{tr}B). Let (*) be the Hodge
involution on (\bigwedge^2S^\perp), and put

\[
\Delta=3c-2,
\qquad U=A-R-\frac{\Delta}{4}I_4.
\]

Because (S) is an eigenvector of (F), the (S\)-to-(S^\perp) block
of (F) vanishes, while its (S^\perp) block is (A+R). In four
dimensions the exact Hodge-contraction identity is

\[
\boxed{
\|R\|^2=\|B\|^2+(\operatorname{tr}B)^2
         -\operatorname{tr}(B*B*).} \tag{5}
\]

Expanding (1), using (\operatorname{tr}B=1-c), and completing the
trace part of (A-R) gives

\[
\boxed{
\frac{E(\mu)}8
=\frac12\|U\|^2+2\|C\|^2+\operatorname{tr}(B*B*)
 -\frac{\Delta^2}{24}.} \tag{6}
\]

Since (B\succeq0) and (*B*\succeq0),

\[
\operatorname{tr}(B*B*)\ge0.
\]

Thus (4) is equivalently the basis-free tangent estimate

\[
\boxed{
\Delta^2\le
12\|U\|^2+48\|C\|^2+24\operatorname{tr}(B*B*).} \tag{7}
\]

For a general two-fermion state (7) is false: the uniform four-edge star

\[
G_\star={1\over4}\sum_{i=1}^4
|e_0\wedge e_i\rangle\langle e_0\wedge e_i|
\]

has (\Delta=1) and all three terms on the right of (7) zero. It is not a
tangent-Veronese mixture. Indeed, two distinct tangent planes (W_x,W_y)
share a line only if (x\perp y), and a fixed (S\in W_x) can occur for at
most two orthogonal projective roots (x). Formula (7) is the required
quantitative version of that qualitative exclusion.

Equations (3) and (6) exhibit exactly the desired architecture:
Motzkin--Straus handles the capped region, while a special
SOS/Pluecker/Hodge estimate must exclude the four-star in the complement.
The remaining unproved statement is (7) for convex mixtures of the tangent
Veronese orbit.
