# The root-determinant weight: exact expansion and zero-set classification

This note records exact facts about the three-sample functional

\[
 F(\mu)=8\,\mathbb E_{X,Y,Z}
 \left[(X\!\cdot\!Y)^2(X\!\cdot\!Z)^2
 \bigl(Y\!\cdot\!Z-(X\!\cdot\!Y)(X\!\cdot\!Z)\bigr)^2
 \det\operatorname{Gram}(X,Y,Z)\right].
\]

It does **not** claim the still-needed inequality `F E >= 0` or
`E >= 24 F`.

## 1. Pointwise square and the Lorentz determinant

Put

\[
 a=X\cdot Y,\qquad b=X\cdot Z,\qquad c=Y\cdot Z,
 \qquad D=1+2abc-a^2-b^2-c^2.
\]

Since `D = det(X,Y,Z)^2`, the integrand is the pointwise square

\[
 8\,[ab(c-ab)\det(X,Y,Z)]^2\ge0.
\]

Equivalently, choose an oriented orthonormal tangent frame `(u,v)` at `X`
and set

\[
 A_X=\int (X\cdot Y)^2(1-(X\cdot Y)^2)^2\,d\mu(Y),
 \qquad
 B_X=\int (X\cdot Y)^2((u+iv)\cdot Y)^4\,d\mu(Y).
\]

Then

\[
 F=\int(A_X^2-|B_X|^2)\,d\mu(X)\ge0.
\]

Indeed, for two leaves `Y,Z`, write `R=c-ab`.  The projected complex
inner product is `R+i det(X,Y,Z)`, and hence

\[
 \operatorname{Re}(R+i\sqrt D)^4=R^4-6R^2D+D^2,
 \qquad
 ((1-a^2)(1-b^2))^2=(R^2+D)^2.
\]

Their difference is `8 R^2 D`.

## 2. Exact triangle-label vector

For

\[
 \tau_{ijk}=\mathbb E[(X\cdot Y)^i(X\cdot Z)^j(Y\cdot Z)^k]
\]

with the exponents sorted canonically,

\[
\boxed{
\begin{aligned}
F={}&8\tau_{044}-16\tau_{046}-16\tau_{133}+32\tau_{135}
 +16\tau_{155}\\
 &+8\tau_{222}-24\tau_{224}-40\tau_{244}+32\tau_{333}.
\end{aligned}}
\]

Exact checks give `F(sigma_S2)=64/11025`, while `F=0` on the ONB and on
every pole--equator zero measure.

## 3. Symmetric kernel

Exchangeability gives the symmetric form

\[
 F=\mathbb E h(X,Y,Z),
\]

where

\[
 h=\frac83D\left[
 a^2b^2(c-ab)^2+a^2c^2(b-ac)^2+b^2c^2(a-bc)^2
 \right]\ge0.
\]

## 4. Exact classification of `F=0`

**Proposition.** If `F(mu)=0`, then the projective support of `mu` is
contained either

1. in a projective line (a great circle on `S^2`), or
2. in the union of one pole line and its orthogonal projective line.

Consequently every `F=0` measure has `E(mu)>=0`.

**Proof.** Since `h` is continuous and nonnegative, `F=0` implies `h=0`
on every triple from `supp(mu)`.  For a linearly independent triple `D>0`.
If none of `a,b,c` is zero, the three square summands force

\[
 c=ab,\qquad b=ac,\qquad a=bc.
\]

The first two equations imply `a^2=1`, contradicting independence.  Thus
one edge is zero.  If, say, `a=0`, the third square is `b^4c^4`, so `b=0`
or `c=0`.  Hence every independent triple contains one line orthogonal to
the other two.

If the support spans at most a plane, case 1 holds.  Otherwise choose an
independent triple and write it as `(p,y,z)` with `p` orthogonal to both
`y,z`; put `P=p^perp`.  If `y` and `z` are not orthogonal, every support
line outside `P` must be orthogonal to both `y,z`, hence must equal the pole
`p`.  This is case 2.

It remains to handle the case where the chosen triple is an ONB.  Relative
to that ONB, the double-orthogonality property first shows that every other
support vector has at most two nonzero coordinates.  If a non-axis support
line exists, it lies in one coordinate plane and supplies a nonorthogonal
pair spanning that plane; the preceding paragraph applies with the remaining
axis as pole.  If no such line exists, the support consists only of the three
axes and is contained in a pole--orthogonal-plane union as well.

Finally, a measure on one great circle satisfies

\[
 E=|\widehat\mu(2)|^2+|\widehat\mu(6)|^2+\frac23>0.
\]

For `mu=w delta_p+(1-w)nu`, with `nu` on `p^perp`, the exact identity is

\[
 E=6(w-\tfrac13)^2+(1-w)^2
 (|\widehat\nu(2)|^2+|\widehat\nu(6)|^2)\ge0.
\]

This proves the consequence.  QED.

## 5. What remains

The classification reduces the original conjecture to showing that a
hypothetical negative global minimizer must satisfy `F=0`.  Either of the
following exact inequalities would do:

\[
 FE\ge0
 \quad\text{or}\quad
 E\ge cF\quad(c>0).
\]

Numerics single out the sharp candidate `c=24`, but no numerical lower bound
or partial stratum argument is an exact certificate for either statement.
