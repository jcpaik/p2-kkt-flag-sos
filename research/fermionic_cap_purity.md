# A sharp capped-1RDM fermionic purity inequality

Let \(V=\mathbb R^5\), let \(G\succeq0\) be a trace-one operator on
\(\bigwedge^2V\), and let

\[
 F=\operatorname{Tr}_2G
\]

be its one-particle contraction, normalized so that
\(\operatorname{tr}F=2\). Then the following implication is exact.

**Theorem.** If \(F\preceq \frac23I\), then

\[
 \boxed{\quad
 \|G\|_{\rm HS}^2-\frac12\|F\|_{\rm HS}^2+\frac13\geq0.
 \quad}                                                     \tag{1}
\]

No decomposability or representing-measure assumption on \(G\) is needed.

## 1. Pinching reduces the theorem to a weighted graph

Choose an orthonormal eigenbasis \(e_1,\ldots,e_5\) of \(F\). Average \(G\)
under all diagonal sign changes

\[
 e_i\longmapsto\epsilon_i e_i,
 \qquad \epsilon_i\in\{\pm1\}.
\]

This is the Hilbert--Schmidt orthogonal projection onto the diagonal algebra
in the wedge basis \(e_i\wedge e_j\). Consequently the pinched state
\(G_0\) satisfies

\[
 \|G_0\|_{\rm HS}\leq\|G\|_{\rm HS},\qquad
 \operatorname{Tr}_2G_0=F,
\]

and has the form

\[
 G_0=\sum_{i<j}p_{ij}
 |e_i\wedge e_j\rangle\langle e_i\wedge e_j|,
 \qquad p_{ij}\geq0,\qquad \sum_{i<j}p_{ij}=1.             \tag{2}
\]

The eigenvalues of \(F\) are the weighted degrees

\[
 d_i=\sum_{j\ne i}p_{ij}\leq\frac23.                       \tag{3}
\]

Let

\[
 A(p)=\sum_{\substack{e<f\\e\cap f\ne\varnothing}}p_ep_f
\]

be the sum over distinct adjacent edges of \(K_5\). Expanding the five
degree squares gives

\[
 \sum_i d_i^2=2\sum_ep_e^2+2A(p).
\]

Thus the left side of (1), for \(G_0\), is exactly

\[
 \frac13-A(p).                                             \tag{4}
\]

It remains to prove \(A(p)\leq1/3\).

## 2. The capped weighted-\(K_5\) lemma

**Lemma.** If nonnegative weights of total one are put on the edges of
\(K_5\), and every weighted vertex degree is at most \(2/3\), then

\[
 \sum_{\substack{e<f\\e\cap f\ne\varnothing}}p_ep_f\leq\frac13. \tag{5}
\]

**Proof.** Maximize \(A\) on this compact polytope, and among maximizers
choose one with the smallest support. Suppose first that no degree
constraint is active. If two support edges \(e,f\) are disjoint, transfer
an amount \(t\) from \(e\) to \(f\). Because \(e,f\) are nonadjacent in the
line graph, \(A\) is affine in \(t\). At a maximum its slope is zero;
otherwise a sufficiently small transfer in one direction increases \(A\).
The transfer may therefore be continued without changing \(A\) until either
one weight vanishes or a degree constraint becomes active. The first event
contradicts minimal support, while the second reduces to the active-constraint
case below.

If there are no disjoint support edges, the support is a clique in the line
graph of \(K_5\). Such a clique is contained either in a triangle or in a
star. On a triangle, Motzkin--Straus (or
\(uv+uw+vw\leq(u+v+w)^2/3\)) gives \(A\leq1/3\). A star carrying total
weight one has center degree one and is infeasible. Hence it only remains
to handle a maximizer with an active degree constraint.

Relabel so that \(d_1=2/3\). Put

\[
 a_i=p_{1i}\quad(2\leq i\leq5),
 \qquad q_{ij}=p_{ij}\quad(2\leq i<j\leq5),
\]

and let \(r_i=\sum_{j\ne i,\,j\geq2}q_{ij}\) be the outer weighted degrees.
Then

\[
 \sum_i a_i=\frac23,\qquad \sum_{i<j}q_{ij}=\frac13.
\]

Writing \(A_{\rm out}\) for the adjacent-edge sum among the \(q\)'s,

\[
\begin{aligned}
 A
 &=\frac12\left(\frac49-\|a\|^2\right)
   +a\mathbin\cdot r+A_{\rm out},\\
 A_{\rm out}&=\frac12\|r\|^2-\sum q_{ij}^2.
\end{aligned}
\]

Completing the square yields

\[
 A=\frac29+\|r\|^2-\sum q_{ij}^2
       -\frac12\|a-r\|^2.                                 \tag{6}
\]

Now \(\|r\|^2-\sum q_{ij}^2\) is the total weight of ordered pairs of
outer edges which are equal or meet. It is at most the weight of all ordered
pairs, namely

\[
 \left(\sum q_{ij}\right)^2=\frac19.
\]

Equation (6) proves \(A\leq2/9+1/9=1/3\), completing the lemma and the
theorem. \(\square\)

## 3. Consequence for isotropic cubic Hankel states

For the tangent-Veronese reformulation, put

\[
 Q_x=\frac{I-\pi_2(\rho_x)}2,
 \qquad F=L(Q_x),
 \qquad G=L(\bigwedge^2Q_x).
\]

The kernel energy is

\[
 E=8\left(\|G\|_{\rm HS}^2-\frac12\|F\|_{\rm HS}^2+\frac13\right). \tag{7}
\]

Assume the first marginal of the normalized degree-six Hankel functional is
isotropic, \(R_1=I/3\). For every traceless symmetric \(S\) with
\(\|S\|_{\rm HS}=1\), direct expansion of \(Q_x\) gives

\[
\begin{aligned}
 \langle S,F[S]\rangle
 &=2L\!\left(r^4x^{\mathsf T}S^2x
       -r^2(x^{\mathsf T}Sx)^2\right)\\
 &=\frac23-2\sum_{k=1}^3
      L\!\left([x_k(x^{\mathsf T}Sx)]^2\right)
 \leq\frac23.                                             \tag{8}
\end{aligned}
\]

The last inequality uses only positivity of the middle cubic catalecticant.
Thus \(F\preceq(2/3)I\), and (1) and (7) prove

\[
 \boxed{E\geq0}
\]

on the full isotropic PSD Hankel cone, including pseudo-moment states with no
representing measure.
