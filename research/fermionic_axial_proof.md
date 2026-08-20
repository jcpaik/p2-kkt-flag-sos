# Exact axial-top closure of the uncapped fermionic inequality

This note proves the fermionic purity inequality whenever a top orbital of
the one-particle marginal has a repeated eigenvalue as a traceless
quadratic form.  No symmetry of the measure is assumed.

## Setup

Let (W=H_2), let (G\succeq0) be the tangent two-fermion moment matrix on
(\Lambda^2W), and let (F=2\operatorname{Tr}_2G).  Normalize
(\operatorname{tr}G=1), hence (\operatorname{tr}F=2).  For a unit top
eigenvector (S\in W), write (c=\langle S,FS\rangle) and
(\delta=3c-2).  Relative to

\[
\Lambda^2W=(S\wedge S^\perp)\oplus\Lambda^2S^\perp,
\qquad
G=\begin{pmatrix}A&C\\ C^*&B\end{pmatrix},
\]

put

\[
U=A-\gamma(B)-\frac{\delta}{4}I_4.
\]

The exact Motzkin split is

\[
P:=\|G\|_F^2-\frac12\|F\|_F^2+\frac13
=\frac12\|U\|_F^2+2\|C\|_F^2
\operatorname{tr}(B*B*)-\frac{\delta^2}{24}. \tag{1}
\]

Because (B\succeq0) and (*B*\succeq0),
(\operatorname{tr}(B*B*)\ge0).

Assume now that (S) is axial.  After a spatial rotation,

\[
S=\frac1{\sqrt6}\operatorname{diag}(1,-2,1).
\]

The capped case (c\le2/3) is already covered by the exact fermionic
Motzkin--Straus lemma.  It remains to treat (c>2/3).

## Axial twirl

Average the measure over rotations in the (xz)-plane, which are exactly
the connected stabilizer of (S).  Denote averaged blocks by
\(\bar U,\bar C\).  The scalar (c), and hence (\delta), is unchanged.
The block construction is equivariant and the group acts orthogonally, so
Jensen gives

\[
\|U\|_F^2\ge\|\bar U\|_F^2,
\qquad
\|C\|_F^2\ge\|\bar C\|_F^2. \tag{2}
\]

For (x=(x,y,z)\in S^2), set (t=y^2), and define

\[
a=\int t(1-t)\,d\mu, \qquad
L=\int(3t-13t^2+12t^3)\,d\mu, \qquad
V=\int(t-4t^2+3t^3)\,d\mu. \tag{3}
\]

Use the orthonormal completion

\[
H=\frac1{\sqrt2}\operatorname{diag}(1,0,-1),
\quad E_{xy},E_{xz},E_{yz}
\]

of (S).  Direct averaging of the tangent Pluecker feature gives

\[
\bar U=\frac L4\operatorname{diag}(1,-1,1,-1), \tag{4}
\]

while (\bar C) has exactly four nonzero entries, each equal up to sign
to (\sqrt3V/4).  Consequently

\[
2\|\bar U\|_F^2+8\|\bar C\|_F^2
=\frac12L^2+6V^2. \tag{5}
\]

The same pointwise calculation gives

\[
c=3a,
\qquad
\delta=9a-2. \tag{6}
\]

There is a pointwise factorization on (0\le t\le1):

\[
-\bigl(3t-13t^2+12t^3\bigr)
-2\bigl(t-4t^2+3t^3\bigr)-9t(1-t)+2
=2(1-t)(3t-1)^2\ge0. \tag{7}
\]

After integration, (6), (7), and the uncapped assumption (\delta>0)
give

\[
-(L+2V)\ge\delta>0. \tag{8}
\]

Weighted Cauchy--Schwarz, with dual norm
(2\cdot1^2+2^2/6=8/3), gives

\[
\frac12L^2+6V^2\ge\frac38(L+2V)^2
\ge\frac38\delta^2\ge\frac16\delta^2. \tag{9}
\]

Combining (2), (5), and (9),

\[
2\|U\|_F^2+8\|C\|_F^2\ge\frac{\delta^2}{6}. \tag{10}
\]

Multiplying (1) by four and using
(4\operatorname{tr}(B*B*)\ge0) now yields (4P\ge0).  Thus the original
copositivity functional (E=8P) is nonnegative for every measure whose
top marginal orbital is axial.

The algebra behind (4)--(6) is checked exactly in
`fermionic_axial_verify.py`.
