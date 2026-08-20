# Exact fermionic proof of the isotropic inequality

This note proves the target exactly for every isotropic probability measure
on `RP^2`.  It is independent of `PLAN.md` and uses only the tangent-plane
exterior-square representation.

## 1. Tangent two-planes

Let `H_2=Sym^2_0(R^3)` with the Frobenius inner product.  For a unit vector
`x`, let

\[
 W_x=\{(xv^T+vx^T)/\sqrt2:v\perp x\}\subset H_2
\]

and denote by `Q_x` the orthogonal projection onto this two-plane.  If
`z_x` is either unit orientation of `W_x`, regarded as a unit vector of
`wedge^2 H_2`, put `R_x=|z_x><z_x|`.  Define

\[
 F=\int Q_x\,d\mu(x),\qquad G=\int R_x\,d\mu(x).
\]

Thus `F>=0`, `tr F=2`, `G>=0`, and `tr G=1`.  The one-particle contraction
of a normalized decomposable bivector satisfies

\[
 2\operatorname{Tr}_2 |u\wedge v\rangle\langle u\wedge v|
 =|u\rangle\langle u|+|v\rangle\langle v|,
\]

so, by averaging,

\[
 F=2\operatorname{Tr}_2G. \tag{1}
\]

If `s=(x^Ty)^2`, the two squared principal cosines between `W_x` and
`W_y` are `s` and `(2s-1)^2`.  Consequently

\[
 \operatorname{tr}(Q_xQ_y)=1-3s+4s^2,
 \qquad
 \operatorname{tr}(R_xR_y)=s(2s-1)^2.
\]

Direct substitution in the kernel gives

\[
 \frac{K(x^Ty)}4
 =2\operatorname{tr}(R_xR_y)-\operatorname{tr}(Q_xQ_y)+\frac23.
\]

After averaging twice,

\[
 E(\mu)=8\|G\|_{HS}^2-4\|F\|_{HS}^2+\frac83. \tag{2}
\]

## 2. Isotropy gives the occupation cap

For `S in H_2`, projection onto the displayed orthonormal tangent basis
gives

\[
 \langle S,Q_xS\rangle
 =2\bigl(|Sx|^2-(x^TSx)^2\bigr).
\]

If `mu` is isotropic, `int xx^T dmu=I/3`, and hence

\[
 \langle S,FS\rangle
 =\frac23\|S\|^2-2\int(x^TSx)^2d\mu(x)
 \le \frac23\|S\|^2.
\]

Therefore

\[
 0\preceq F\preceq\frac23I_5. \tag{3}
\]

## 3. A two-fermion purity lemma

**Lemma.**  Let `G>=0` be a real two-fermion density operator on
`wedge^2 R^5`, with `tr G=1`, and put `F=2 Tr_2 G`.  If
`F<=(2/3)I`, then

\[
 \|G\|_{HS}^2-\frac12\|F\|_{HS}^2+\frac13\ge0. \tag{4}
\]

**Proof.**  Orthogonally diagonalize `F`.  In that basis, average `G` over
all diagonal sign matrices `diag(epsilon_1,...,epsilon_5)`, acting through
their exterior square.  Distinct wedge-basis vectors `e_i wedge e_j` have
distinct sign characters, so the average is diagonal:

\[
 \bar G=\sum_{i<j}p_{ij}|e_i\wedge e_j\rangle
                    \langle e_i\wedge e_j|,
 \qquad p_{ij}\ge0,\quad\sum_{i<j}p_{ij}=1.
\]

The averaging fixes `F`, while orthogonal projection cannot increase the
Hilbert--Schmidt norm, so `||G||^2>=||bar G||^2`.  If

\[
 d_i=\sum_{j\ne i}p_{ij},
\]

then `F=diag(d_1,...,d_5)` and (3) is exactly `d_i<=2/3`.

View the ten `p_ij` as weights on the edges of `K_5`.  Let

\[
 A=\sum_{\substack{e<f\\e\cap f\ne\varnothing}}p_ep_f
\]

be the sum over distinct adjacent edge pairs.  Since every such pair has a
unique common endpoint,

\[
 \frac12\sum_i d_i^2=\sum_ep_e^2+A.
\]

It remains to prove `A<=1/3`; then

\[
 \|\bar G\|^2-\frac12\|F\|^2+\frac13
 =\sum_ep_e^2-\frac12\sum_id_i^2+\frac13
 =\frac13-A\ge0.
\]

We prove the edge inequality by the standard Motzkin--Straus compression.
First suppose some vertex, say `0`, is saturated: `d_0=2/3`.  Write
`a_i=p_{0i}` for `i=1,...,4`, and write `b_ij=p_ij` for the six inner
edges.  Then

\[
 \sum_i a_i=\frac23,\qquad \sum_{i<j}b_{ij}=\frac13.
\]

Let `r_i=sum_{j ne i}b_ij`, and let `A_b` be the adjacent-pair sum among
the inner edges.  Splitting adjacent pairs into star--star, star--inner,
and inner--inner pairs gives

\[
 A=\frac29-\frac12\|a\|^2+a\cdot r+A_b.
\]

Moreover `||r||^2=2 sum b_ij^2+2A_b`.  Completing the square therefore
gives the exact identity

\[
 A=\frac29+\sum_{i<j}b_{ij}^2+2A_b-\frac12\|a-r\|^2.
\]

The sum `sum b_ij^2+2A_b` omits only the nonnegative products of disjoint
inner edges from `(sum b_ij)^2=1/9`.  Hence `A<=2/9+1/9=1/3`.

Now suppose no vertex is saturated.  If two positive underlying edges are
disjoint, keep their combined weight fixed and transfer weight from one to
the other.  Because the two corresponding vertices of the line graph are
nonadjacent, `A` is affine in this transfer parameter.  Move in the
nondecreasing direction until either one weight becomes zero or one of the
four affected vertex degrees reaches `2/3`.  In the latter case the previous
paragraph applies; in the former case the support strictly shrinks.  Repeating
finitely many times either reaches a saturated vertex or leaves a
pairwise-intersecting family of two-subsets of a five-set.

Such a family is contained either in a star or in a triangle: if it contains
`{1,2}` and ` {1,3}`, every further member contains `1` or is `{2,3}`; in
the latter case all members lie in that triangle.  A star carrying total
weight one violates the degree cap at its center.  On a triangle the three
degree caps force every edge weight to be at least `1/3`, hence all three are
`1/3`, and `A=1/3`.  This proves the lemma.  QED

## 4. Conclusion

Apply the lemma to the operators `F,G` from Section 1.  Equation (3) supplies
its only nonautomatic hypothesis, and (2) yields

\[
 E(\mu)=8\left(\|G\|^2-\frac12\|F\|^2+\frac13\right)\ge0.
\]

Thus the kernel is copositive on the full isotropic slice, sharply; the
pole--equator zero measures attain equality.

This note does **not** assert an unrestricted isotropization reduction.
