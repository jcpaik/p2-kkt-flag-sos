# Exact removal of the fermionic cap for a rank-two top orbital

This note proves the unrestricted purity inequality whenever a top
eigenvector of the one-particle marginal, regarded as a traceless symmetric
three-by-three matrix, has spectrum \((\lambda,-\lambda,0)\).  It does not
assume that the other four one-particle orbitals diagonalize in any preferred
basis.

Let

\[
 S={\operatorname{diag}(1,-1,0)\over\sqrt2}
\]

and complete it by

\[
 H={\operatorname{diag}(1,1,-2)\over\sqrt6},\quad
 E_{xy},\quad E_{xz},\quad E_{yz}.
\]

Relative to
\(\bigwedge^2H_2=(S\wedge S^\perp)\oplus\bigwedge^2S^\perp\), write

\[
 G=\begin{pmatrix}A&C\\C^T&B\end{pmatrix},\qquad
 R=\gamma(B),\qquad
 \delta=\operatorname{tr}A-2\operatorname{tr}B,
\]

and put

\[
 U=A-R-{\delta\over4}I_4.
\]

The basis-free uncapped Motzkin identity is

\[
 {E\over8}={1\over2}\|U\|^2+2\|C\|^2
 +\operatorname{tr}(B*B*)-{\delta^2\over24}.                 \tag{1}
\]

Here \(*\) is Hodge star on \(\bigwedge^2S^\perp\).  Since \(B\succeq0\)
and \(*B*\succeq0\), their trace pairing is nonnegative.

For a probability measure on the sphere define

\[
 m=\int z^2\,d\mu,\qquad
 \alpha=\int(x^2-y^2)^2\,d\mu.
\]

The Pluecker coordinates satisfy the pointwise identity

\[
 2z_{S,E_{xy}}+z_{E_{xz},E_{yz}}=z.
\]

Consequently

\[
 m=4A_{xy,xy}+B_{xz\wedge yz,xz\wedge yz}
   +4C_{xy,xz\wedge yz}.                                  \tag{2}
\]

More importantly, direct pointwise expansion gives the following linear
identity, hence it remains true after integration:

\[
 \boxed{
 1-3m=-2U_{xy,xy}+U_{xz,xz}+U_{yz,yz}
 +2\sqrt3\bigl(C_{xz,H\wedge xz}-C_{yz,H\wedge yz}\bigr).} \tag{3}
\]

Let \(L=\operatorname{diag}(0,-2,1,1)\) on
\((H,xy,xz,yz)\), and let \(M\) have only the two entries

\[
 M_{xz,H\wedge xz}=2\sqrt3,\qquad
 M_{yz,H\wedge yz}=-2\sqrt3.
\]

Then \(\|L\|^2=6\), \(\|M\|^2=24\), and (3) is
\(1-3m=\langle L,U\rangle+\langle M,C\rangle\).  Weighted
Cauchy--Schwarz therefore gives

\[
 (1-3m)^2
 \leq\left({\|L\|^2\over2}+{\|M\|^2\over8}\right)
       \left(2\|U\|^2+8\|C\|^2\right)
 =6\left(2\|U\|^2+8\|C\|^2\right).                       \tag{4}
\]

Adding the nonnegative Hodge term yields the exact scalar estimate

\[
 2\|U\|^2+8\|C\|^2+4\operatorname{tr}(B*B*)
 \geq{(1-3m)^2\over6}.                                    \tag{5}
\]

For the chosen \(S\),

\[
 \langle S,FS\rangle=1-m-\alpha,qquad
 \delta=1-3m-3\alpha.
\]

Combining (1) and (5) gives

\[
 \boxed{
 {E\over8}\geq {3\over8}\alpha
       \left({2\over3}-2m-\alpha\right).}                 \tag{6}
\]

If \(S\) is a top orbital with occupation greater than \(2/3\), then
\(m+\alpha<1/3\), and hence

\[
 {2\over3}-2m-\alpha
 =2\left({1\over3}-m-\alpha\right)+\alpha>0.
\]

Thus \(E\geq0\) throughout the previously uncapped region for every
measure whose top orbital has rank two.  If all occupations are at most
\(2/3\), the capped Motzkin--Straus lemma already gives \(E\geq0\).

The pointwise identities and constants are checked symbolically by
`fermionic_rank2_scalar_verify.py`.
