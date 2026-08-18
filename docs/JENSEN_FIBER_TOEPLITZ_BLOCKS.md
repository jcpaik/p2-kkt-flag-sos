# The Jensen and fiber-Toeplitz blocks: definitions and validity

*Standalone reference (2026-08-18).  These are the two constraint
families that broke the weighted pole's inertness; full campaign log
with measurements in
[Toeplitz blocks note](TOEPLITZ_BLOCKS_NOTE.md), origin story in
[Cylindrical domination](CYLINDRICAL_DOMINATION.md) §6 (the
circle-pair theorem's Toeplitz coupling) and
[Unprojected escape note](UNPROJECTED_ESCAPE_NOTE.md).
Implementation: `toeplitz_blocks.py` (exact expansions + self-tests),
`toeplitz_export.py` (first-class exports).*

## 0. Setting

Roots $x_1,x_2$ and leaves $y,y'$ are sampled independently from an
antipodal probability measure $\mu$ on $S^2$.  Write

$$t_1=x_1\!\cdot\!y,\quad t_2=x_2\!\cdot\!y,\quad s=x_1\!\cdot\!x_2 .$$

The azimuth $\varphi$ of the leaf about the root-pair frame enters
polynomially through the complex combination

$$w:=(t_2-st_1)+i\det(x_1,x_2,y),\qquad
|w|^2=(1-t_1^2)(1-s^2),\qquad w=\sin\delta\,\sin\theta\,e^{i\varphi}$$

(frame $x_1=e_3$, $x_2$ at angle $\delta$ from $x_1$, leaf at polar
angle $\theta$).  $\operatorname{Re}w^k$ is a polynomial in
$(t_1,t_2,s)$; $\operatorname{Im}w^k$ carries one factor
$\det(x_1,x_2,y)$.  Expectations of leaf monomials reduce to exact
rational combinations of the moment labels: same-leaf products give
3-sample (`triangle`) labels, two-leaf products give 4-sample
(`graph_4`) labels, and $h_2$-multiplied copies give the disconnected
$p_2\times(\cdot)$ product labels.

## 1. The Jensen (averaging-contraction) block

**Data.**  A vector $\varphi=(\varphi_\alpha)$ of polynomial leaf
functions of the rooted triple — monomials $t_1^it_2^js^k$ in the
leaf-even parity sectors, or $\det(x_1,x_2,y)\cdot t_1^it_2^js^k$ in
the odd sectors — and a nonnegative root weight
$\rho\in\{1,\ 1-s^2\}$.

**The two matrices.**

$$T_{ab}=\mathbb E\big[\rho\,\varphi_a(y)\,\varphi_b(y)\big]
\quad(\text{same leaf}),\qquad
G_{ab}=\mathbb E\big[\rho\,\varphi_a(y)\,\varphi_b(y')\big]
\quad(\text{independent leaves}).$$

$G$ is the classical two-root flag Gram (already in the hierarchy);
$T$ is its same-leaf shadow (triangle sector).

**The block.**

$$\boxed{\;T-G\;=\;\mathbb E_{x_1,x_2}\big[\rho\cdot
\mathrm{Cov}_y(\varphi\mid x_1,x_2)\big]\;\succeq\;0\;}$$

**Validity (all measures).**  $\mathrm{Cov}_y(\varphi)\succeq0$
pointwise in the roots and $\rho\ge0$; take expectations.  No KKT
content; composes with the reduction lemma.  Equivalently, with the
mixed instant/averaged vector $f=(\varphi(y),\,\mathbb E_y\varphi)$,
$\mathbb E[ff^{\top}]=\begin{pmatrix}T&G\\G&G\end{pmatrix}\succeq0$,
which is exactly $G\succeq0$ *and* $T-G\succeq0$: the "unfolding"
inequality $[\![A]\!]^2\le[\![A^2]\!]$ of flag calculus, one new
matrix inequality per parity sector.

**Valid companions** (using $0\le h_2=\tfrac{3p_2-1}2\le1$ for every
measure):

$$h_2\,(T-G)\succeq0,\qquad (1-h_2)\,G\succeq0,\qquad
(1-h_2)\,(T-G)\succeq0 .$$

The $h_2$-localized copies load on the $p_2\times$(triangle/graph$_4$)
sectors — the measured escape's home; the $(1-h_2)$ "seesaw"
complements couple the plain and localized towers (the escape was
measured shuttling mass between them at a uniform ratio $\approx-4.9$,
a direction nothing else constrains).

**Unrooted (pair-sample) copy.**  Taking a whole pair $(X,Y)$ as the
sample with flags $(X\!\cdot\!Y)^d$:
$T_{ab}=p_{a+b}$ (pair labels), $G_{ab}=p_ap_b$ (product labels),
$T-G\succeq0$ — a Hankel-vs-product coupling not implied by any Gram.

**Why this is the circle-pair Toeplitz mechanism.**  The double-angle
coupling $1+\operatorname{Re}\zeta_{2k}=2\mathbb E\cos^2k\varphi\ge
2(\mathbb E\cos k\varphi)^2$ that closed the circle-pair theorem is a
diagonal entry of $T-G$ after the automatic double-angle reduction
($(\operatorname{Re}w^k)^2=\tfrac12[(|w|^2)^k+\operatorname{Re}w^{2k}]$
is a polynomial identity the moment reducer applies for free).  Grams
over modulated leaves with *independent* leaves are already spanned by
the two-root blocks — the mechanism is reachable only by mixing
same-leaf with averaged-leaf entries, which is what $T-G$ does.

**Sharpness / saturation.**  $\mathrm{Cov}_y(\varphi)=0$ exactly on
fiber-deterministic configurations (leaf angle a function of the
roots).  Hence the blocks vanish identically on coherent one-orbit
measures (ONB-face compatible), and, dually, a deep escape direction
can *saturate* them by concentrating fibers — the measured
re-steepening mechanism at degree 14.

## 2. The fiber-Toeplitz block

**Data.**  A radial depth $K\ge2$, a set of root polynomials
$g_a(t_1,t_2,s)$, and the frame variable $w$ of §0.

**The matrix.**  Index by pairs $(j,a)$, $j=0,\dots,K$:
$V_{(j,a)}=|w|^{K}e^{2ij\varphi}\,g_a$.  Then
$\mathbb E_y[VV^{H}]$ (conditional on the roots, then averaged over
them) is Hermitian PSD, and its real part has polynomial entries

$$\boxed{\;M_{(j,a),(k,b)}
=\mathbb E\big[(|w|^2)^{\,K-|j-k|}\,
\operatorname{Re}\!\big(w^{2(j-k)}\big)\,g_a\,g_b\big]\;\succeq\;0\;}$$

— the radially-weighted trigonometric moment matrix of the leaf's
azimuthal fiber distribution $(\sin\delta\sin\theta)^{2K}d\mu_y$
pushed to the circle.  Entries are `triangle` labels; $h_2$-localized
copies give $p_2\times$`triangle`.  Only the cosine
($\operatorname{Re}$) part is used — the imaginary parts carry single
$\det$ factors whose 3-point moments are outside the label algebra —
and $\operatorname{Re}$ of a Hermitian PSD matrix is PSD.

**Validity.**  Positivity is *moment-matrix* positivity of a
nonnegative fiber measure — not a sum of squares.  Machine-checked:
exact rational equality of the label expansion against direct
evaluation, and exact PSD at the ONB and 4-point-cross measures.

**Why no Gram spans it (the Fejér–Riesz obstruction).**  Every band
of $M$ carries the *constant* radial weight $|w|^{2K}$, while a
polynomial Gram $\mathbb E[v_j\bar v_k]$ with $v_j=w^{2j}q_j$ forces
band weights $|w|^{4\min(j,k)}q_jq_k$; matching them requires
$q_j=|w|^{K-2j}$ — non-polynomial for $j>K/2$.  For $K\le1$ the
content reduces to squares; for $K\ge2$ the block is strictly outside
the within-degree square cone.  This is the flag-algebra incarnation
of the pointwise $7\times7$ fiber-Toeplitz PSD identified as the
honest frontier of the cylindrical route.

**Pair-sector companions.**  The localized Hankel
$[\,p_{a+b}-p_{a+b+2}\,]\succeq0$ (moment matrix of $(1-t^2)d\nu$ on
pair moments) and the weighted pair-Jensen
$\mathrm{Cov}\big((1-t^2)t^a\big)\succeq0$, with $h_2$-localized
copies landing on $p_2\times$pair(-product) labels.

## 3. Measured effect (summary; details in the campaign log)

Stacked on the degree-14/16 weighted problems ($h_2E$ target), these
families produced, in three design iterations (v1: 11 families, v2:
29, v3: 43 incl. the $K\le4$ fiber-Toeplitz tower):

* the first nonzero movement of the weighted selector trace law ever
  measured, ending at exponent $\varepsilon^{-0.21}$ (from
  $\varepsilon^{-1.03}$) with no deep re-steepening at degree 16;
* dual-bound improvements at fixed degree:
  $-2.79\times10^{-5}\to-3.66\times10^{-6}$ (deg 14),
  $-1.31\times10^{-6}\to-5.25\times10^{-9}$ (deg 16, the current
  record wall);
* the $\varepsilon$-feasibility repair of the sharp-face-projected
  (we1) cone, validating the weighted-(E1) prediction that the missing
  objects were admissible-tower cuts.

Every family is exported in exact rational arithmetic with an
in-run verification (label expansion vs direct evaluation, exact PSD
at sharp measures); validity ledgers in the campaign log §3 and §8.
