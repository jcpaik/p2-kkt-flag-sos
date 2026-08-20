# Exact failure of the global affine fermionic relaxation

Let $V=H_2(\mathbb R^3)$, with orthonormal basis

\[
\begin{aligned}
 E_0&=\frac1{\sqrt2}\operatorname{diag}(1,-1,0),&
 E_1&=\frac1{\sqrt6}\operatorname{diag}(1,1,-2),\\
 E_2&=\frac{e_1e_2^T+e_2e_1^T}{\sqrt2},&
 E_3&=\frac{e_1e_3^T+e_3e_1^T}{\sqrt2},&
 E_4&=\frac{e_2e_3^T+e_3e_2^T}{\sqrt2}.
\end{aligned}
\]

For $i<j$, write $e_{ij}=E_i\wedge E_j$.  In the pair ordering

\[
01,02,03,04,12,13,14,23,24,34,
\]

consider the diagonal density on $\bigwedge^2V$

\[
 G=\sum_{i<j}p_{ij}|e_{ij}\rangle\langle e_{ij}|,
\qquad
 p=\left(0,\frac2{43},0,0,\frac{14}{43},\frac2{43},
 \frac2{43},\frac{23}{86},\frac{23}{86},0\right).
 \tag{1}
\]

Clearly $G\succeq0$ and $\operatorname{tr}G=1$.  Its one-particle
contraction is diagonal:

\[
 F=\operatorname{Tr}_2G
 =\operatorname{diag}\left(
 \frac2{43},\frac{18}{43},\frac{39}{43},
 \frac{27}{86},\frac{27}{86}\right).
 \tag{2}
\]

## The exact affine Clebsch--Gordan relation

For a rotation generator $J_a:S\mapsto[L_a,S]$ on $V$, identify its
two-form with a vector in $\bigwedge^2V$, and normalize it by $1/\sqrt5$.
The three resulting orthonormal vectors span the $H_1$ summand in

\[
 \bigwedge^2H_2=H_1\oplus H_3.
\]

In the edge basis above, their rows are

\[
 W_1=
 \begin{pmatrix}
0&0&0\\
0&0&-2/\sqrt5\\
0&1/\sqrt5&0\\
1/\sqrt5&0&0\\
0&0&0\\
0&\sqrt{3/5}&0\\
-\sqrt{3/5}&0&0\\
-1/\sqrt5&0&0\\
0&1/\sqrt5&0\\
0&0&-1/\sqrt5
\end{pmatrix}.
\tag{3}
\]

Consequently

\[
 M:=5G_{11}=5W_1^T\operatorname{diag}(p)W_1
 =\operatorname{diag}\left(\frac{35}{86},\frac{35}{86},\frac8{43}\right).
 \tag{4}
\]

Thus $\operatorname{tr}G_{11}=1/5$, as for every tangent-Veronese
mixture.  For $A\in H_2$, put

\[
 T_A(S)=AS+SA-\frac23\operatorname{tr}(AS)I.
\]

The exact orbit identity relating the spin-two part of $F$ to $G_{11}$
is

\[
 \Pi_2(F)=\frac37T_{M-I/3}.
 \tag{5}
\]

For (1), both sides of (5), in the $E_i$ basis, are exactly

\[
 \operatorname{diag}\left(
 \frac{19}{301},-\frac{19}{301},\frac{19}{301},
 -\frac{19}{602},-\frac{19}{602}\right).
 \tag{6}
\]

Hence this state satisfies the full proposed global affine
Clebsch--Gordan relaxation.

## Negative purity gap

Directly,

\[
 \|G\|_F^2=\frac{945}{3698},\qquad
 \|F\|_F^2=\frac{4427}{3698}.
\]

Therefore

\[
 \boxed{
 \|G\|_F^2-\frac12\|F\|_F^2+\frac13=-\frac5{516}<0.
 }
\tag{7}
\]

So positivity, the fixed $H_1/H_3$ block traces, and the complete affine
spin-two relation do not prove the desired purity inequality.

There is a stronger pointwise-valid marginal constraint.  If

\[
 \langle S,P_M S\rangle=2\operatorname{tr}(S^2M),
\]

then a true tangent mixture obeys $P_M-F\succeq0$.  The counterexample
above violates it: the five eigenvalues of $P_M-F$ include
$-4/43$.  Thus (7) does not settle the relaxation with that additional
LMI.  Numerical non-diagonal counterexamples to the stronger relaxation are
recorded in `tensor_fermionic_general_relaxation_opt.py`, but they have not
yet been converted into an exact certificate.
