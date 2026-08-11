# Degree-14 five-point command: complete mathematical reference

> **Historical note (2026-08-11).** This document describes the original
> *isotropic* formulation, in which every degree-two sampled vertex was
> contracted through \(\mathbb E[XX^{\mathsf T}]=\tfrac13I\). The current
> `sos_search.py` no longer performs that contraction: `("pair", 2)` and
> all degree-two-vertex graphs are independent labels, the target is
> \(T=-\tfrac14+\tfrac{15}{4}p_2-9p_4+6p_6\), and the isotropy deficit is
> carried by the `harmonic_2`, `harmonic_flag_*`, and `spin2_flag` blocks.
> Label counts, block lists, and the reported bound differ accordingly; see
> [RESULTS.md](../RESULTS.md) for the non-isotropic record. The derivations
> below remain correct for the isotropic branch they were written for.

## 1. Purpose and scope

This document explains every mathematical object enabled by the command

```sh
python3 sos_search.py \
  --dual --summary-only --scale-constraints --rank-relations \
  --degree 14 --no-pointwise-sos \
  --harmonics --three-point-flags --four-point-flags --two-root-flags \
  --max-flag-arity 5 --max-root-factor-degree 2 \
  --gradient --potential --potential-matrices \
  --hessian --four-point-hessian \
  --global-gap --global-tangent-gaps \
  --tolerance 1e-10
```

The command constructs a finite semidefinite relaxation for the isotropic KKT
branch of the \(P_2\)-kernel problem. It solves the **dual moment problem**:
among formal moment functionals satisfying the selected flag, harmonic, KKT,
and Gram-rank constraints, minimize the normalized energy.

Three qualifications are essential.

1. The formal moment variables in the dual need not come from an actual
   probability measure. A negative dual optimum is therefore a relaxation
   gap, not a counterexample.
2. The code uses the isotropy identity
   \[
   \mathbb E[XX^{\mathsf T}]=\frac13I.
   \]
   This command by itself does not justify reducing every possible global
   minimizer to the isotropic class.
3. MOSEK returns floating-point data. A value near zero is not an exact SOS
   certificate unless a primal identity with exact PSD matrices is separately
   recovered and verified.

The phrase **KKT-infused flag algebra** is descriptive rather than a standard
named theorem. In this project it means:

- ordinary flag-square and harmonic positivity;
- exact moment identities from symmetry, isotropy, and dimension three;
- first- and second-order KKT constraints satisfied by a hypothetical global
  minimizer;
- positive semidefinite multipliers for inequalities and unrestricted
  multipliers for equalities.

## 2. The kernel, energy, and normalized target

Let \(\mu\) be an antipodally symmetric probability measure on \(S^2\), or
equivalently a probability measure on
\(\mathbb {RP}^2=S^2/\{x\sim -x\}\). The kernel is

\[
K(t)=32t^6-48t^4+20t^2-\frac43.
\]

Its energy is

\[
E_K(\mu)
  =\iint K(x\cdot y)\,d\mu(x)d\mu(y).
\]

Write

\[
p_r=\mathbb E[(X\cdot Y)^r],
\qquad X,Y\stackrel{\mathrm{iid}}{\sim}\mu.
\]

Then

\[
E_K=32p_6-48p_4+20p_2-\frac43.
\]

The implementation assumes isotropy:

\[
\mathbb E[XX^{\mathsf T}]=\frac13I.
\]

Consequently,

\[
p_2
=\mathbb E_Y\!\left[
  Y^{\mathsf T}\mathbb E[XX^{\mathsf T}]Y
\right]
=\frac13.
\]

The energy therefore becomes

\[
E_K
=32p_6-48p_4+\frac{16}{3}
=\frac{16}{3}\left(1-9p_4+6p_6\right).
\]

The SDP target is

\[
\boxed{T=1-9p_4+6p_6=\frac{3}{16}E_K.}
\]

In the code this is the coefficient vector

```text
constant :  1
pair_4   : -9
pair_6   :  6
```

The dual objective is the formal evaluation of this vector.

## 3. Moment labels as multigraphs

### 3.1 Gram monomials

For independent samples \(X_0,\ldots,X_{n-1}\), a Gram monomial is

\[
m_{\mathbf e}(X_0,\ldots,X_{n-1})
=\prod_{0\le i<j<n}(X_i\cdot X_j)^{e_{ij}},
\]

where every \(e_{ij}\) is a nonnegative integer.

It is convenient to regard \(\mathbf e\) as a loopless multigraph:

- sampled points are vertices;
- the multiplicity of edge \(ij\) is \(e_{ij}\);
- the total Gram degree is \(\sum_{i<j}e_{ij}\).

The corresponding moment is

\[
M_{\mathbf e}=\mathbb E[m_{\mathbf e}].
\]

The program stores only a canonical label for each moment after applying the
following exact reductions.

### 3.2 Antipodal parity

Changing the sign of one independently sampled projective representative
\(X_i\) leaves its law unchanged. The monomial changes by

\[
(-1)^{d_i},
\qquad
d_i=\sum_{j\ne i}e_{ij}.
\]

Hence

\[
d_i\ \text{odd for some }i
\quad\Longrightarrow\quad
M_{\mathbf e}=0.
\]

In graph language, every surviving vertex has even multigraph degree.

### 3.3 Isotropic contraction

If a sampled vertex \(X_i\) has multigraph degree two, it can be integrated
out exactly. For example,

\[
\mathbb E_{X_i}
  [(X_i\cdot u)(X_i\cdot v)]
=u^{\mathsf T}\mathbb E[X_iX_i^{\mathsf T}]v
=\frac13(u\cdot v).
\]

Thus a degree-two path \(u-X_i-v\) contracts to an edge \(u-v\), with factor
\(1/3\). A doubled edge gives

\[
\mathbb E_{X_i}[(X_i\cdot u)^2]=\frac13.
\]

The reducer repeats these contractions until no degree-two vertex remains.

### 3.4 Isolated vertices and disconnected components

An isolated sampled vertex contributes \(1\) and is removed. If the graph
splits into connected components \(H_1,\ldots,H_s\), independence gives

\[
M_H=\prod_{\alpha=1}^s M_{H_\alpha}.
\]

The code represents this by a `product` label.

### 3.5 Exchangeability

The samples are iid, so relabeling the vertices does not change the moment.
Each connected graph is replaced by the lexicographically least edge-exponent
tuple among all vertex permutations.

The surviving labels have forms such as

```text
("constant",)
("pair", r)
("triangle", e01, e02, e12)
("graph_4", ...)
("graph_5", ...)
("product", label_1, label_2, ...)
```

The stated command produces 574 distinct canonical labels.

## 4. What a flag PSD constraint means

Suppose \(R=(X_0,\ldots,X_{r-1})\) is a tuple of shared **roots**, \(Y\) is an
independent leaf, and

\[
f_1(R,Y),\ldots,f_s(R,Y)
\]

are polynomial flag functions in their Gram entries. Define the conditional
flag vector

\[
F(R)
=\mathbb E_Y
  \begin{bmatrix}
  f_1(R,Y)\\
  \vdots\\
  f_s(R,Y)
  \end{bmatrix}.
\]

For every \(c\in\mathbb R^s\),

\[
\mathbb E_R[(c^{\mathsf T}F(R))^2]\ge0.
\]

Equivalently, the flag moment matrix

\[
\mathcal M
=\mathbb E_R[F(R)F(R)^{\mathsf T}]
\succeq0.
\]

Introducing a second independent leaf \(W\) expands its entries as

\[
\mathcal M_{ij}
=\mathbb E_{R,Y,W}[f_i(R,Y)f_j(R,W)].
\]

Each entry is reduced to the canonical moment labels of Section 3. This gives
coefficient matrices \(A_L\) satisfying

\[
\mathcal M(y)=\sum_L y_LA_L,
\]

where \(y_L\) is the dual variable assigned to label \(L\).

This gluing operation explains the arity terminology:

- one root plus two leaves gives a three-point block;
- two roots plus two leaves gives a four-point block;
- three roots plus two leaves gives a five-point block.

Projective sign characters split a flag space into sectors. Only flags with
the same parity character are placed in the same Gram block.

## 5. Primal identities and the dual moment relaxation

Every positive block family \(\beta\) has coefficient matrices
\(A_L^{(\beta)}\). Every equality block \(\gamma\) has matrices
\(B_L^{(\gamma)}\), and every scalar identity \(r\) has coefficients
\(R_L^{(r)}\).

An exact primal certificate would be a coefficient identity

\[
t_L
=\sum_\beta
  \left\langle A_L^{(\beta)},Q_\beta\right\rangle
 +\sum_\gamma
  \left\langle B_L^{(\gamma)},S_\gamma\right\rangle
 +\sum_r\lambda_rR_L^{(r)}
\quad\text{for every label }L,
\]

with

\[
Q_\beta\succeq0,
\]

while \(S_\gamma\) and \(\lambda_r\) are unrestricted because they multiply
identities equal to zero.

The command uses `--dual`, so it solves the conic dual:

\[
\begin{aligned}
\text{minimize}\quad&
  y_{\mathrm{constant}}-9y_{\mathrm{pair},4}
  +6y_{\mathrm{pair},6},\\
\text{subject to}\quad&
  y_{\mathrm{constant}}=1,\\
&\sum_Ly_LA_L^{(\beta)}\succeq0
  &&\text{for every positive block }\beta,\\
&\sum_Ly_LB_L^{(\gamma)}=0
  &&\text{for every matrix equality block }\gamma,\\
&\sum_Ly_LR_L^{(r)}=0
  &&\text{for every scalar identity }r.
\end{aligned}
\]

Every actual isotropic measure satisfying the encoded KKT conditions defines
a feasible vector \(y_L=M_L(\mu)\). The formal feasible set is larger, so the
dual optimum is a lower bound for the actual KKT branch.

## 6. `--degree 14` and `--no-pointwise-sos`

`--degree 14` limits the glued Gram monomials to total degree at most \(14\).
For an ordinary square this means that each half has degree at most \(7\).
For a localizing multiplier of degree \(d\), the square basis is shortened so
that

\[
2\deg(f)+d\le14.
\]

Some geometric kernels contain terms of mixed degree. The implementation uses
their conservative maximum degree when selecting a basis.

By default the script can add an ordinary polynomial SOS block

\[
v(a,b,c)^{\mathsf T}Qv(a,b,c),
\qquad Q\succeq0,
\]

in the three Gram variables of three points. The option
`--no-pointwise-sos` disables that block. The command also does not contain
`--gram-module`, so it does not add the generic three-point Gram-determinant
or principal-minor quadratic module.

This does **not** disable flag squares. All the conditional PSD blocks below
remain active.

## 7. `--three-point-flags`: one-root \(O(2)\) harmonic flags

Let \(X\) be the root and \(Y,Z\) the two leaves. Put

\[
a=X\cdot Y,\qquad b=Y\cdot Z,\qquad c=Z\cdot X.
\]

Project the leaves to the tangent plane at \(X\):

\[
P_XY=Y-aX,\qquad P_XZ=Z-cX.
\]

Then

\[
(P_XY)\cdot(P_XZ)=b-ac
\]

and

\[
\|P_XY\|^2\|P_XZ\|^2=(1-a^2)(1-c^2).
\]

The stabilizer of \(X\) is \(O(2)\). Its weight-\(k\) zonal kernel is the
polynomial

\[
R_k(a,b,c)
=\bigl((1-a^2)(1-c^2)\bigr)^{k/2}
  \cos(k(\phi_Y-\phi_Z)),
\]

where \(\phi_Y,\phi_Z\) are tangent-plane angles. Square roots cancel through
the recurrence

\[
\begin{aligned}
R_0&=1,\\
R_1&=b-ac,\\
R_{k+1}
&=2(b-ac)R_k
 -(1-a^2)(1-c^2)R_{k-1}.
\end{aligned}
\]

For leaf powers \(r,s\), the moment-matrix entry is

\[
\mathbb E[a^r c^sR_k(a,b,c)].
\]

The admissible leaf powers have parity \(r\equiv k\pmod2\), ensuring the
correct projective character. At degree \(14\), the active blocks are:

| block | \(k\) | leaf powers | size |
|---|---:|---|---:|
| `flag_0` | 0 | \(0,2,4,6\) | 4 |
| `flag_1` | 1 | \(1,3,5\) | 3 |
| `flag_2` | 2 | \(0,2,4\) | 3 |
| `flag_3` | 3 | \(1,3\) | 2 |
| `flag_4` | 4 | \(0,2\) | 2 |
| `flag_5` | 5 | \(1\) | 1 |
| `flag_6` | 6 | \(0\) | 1 |

For every \(k\), the dual imposes

\[
\left(
\mathbb E[a^r c^sR_k(a,b,c)]
\right)_{r,s}\succeq0.
\]

## 8. `--four-point-flags`: the empty-type pair block

For an even power \(r\), consider the unrooted two-sample flag average

\[
f_r=\mathbb E_{X,Y}[(X\cdot Y)^r]=p_r.
\]

Squaring a linear combination of these averages gives

\[
\mathbb E[
  (X_0\cdot X_1)^r
  (X_2\cdot X_3)^s
]
=p_rp_s.
\]

Thus

\[
\left(
\mathbb E[
  (X_0\cdot X_1)^r
  (X_2\cdot X_3)^s
]
\right)_{r,s}
\succeq0.
\]

At degree \(14\), the basis is

\[
r\in\{0,2,4,6\},
\]

so `empty_type_flag` is a \(4\times4\) PSD block. This block connects the
single-pair moments to disconnected `product` labels.

## 9. `--two-root-flags`: four-point conditional squares

Let \(X,Z\) be two shared roots and \(Y\) one leaf. A basis flag is

\[
f_{ijk}(X,Z;Y)
=(X\cdot Y)^i(Z\cdot Y)^j(X\cdot Z)^k,
\qquad i+j+k\le7.
\]

Gluing a second leaf \(W\) gives

\[
\begin{aligned}
&f_{ijk}(X,Z;Y)f_{i'j'k'}(X,Z;W)\\
&\quad=
(X\cdot Z)^{k+k'}
(X\cdot Y)^i(X\cdot W)^{i'}
(Z\cdot Y)^j(Z\cdot W)^{j'}.
\end{aligned}
\]

The three parity characters of a flag are

\[
\bigl(i+j,\ i+k,\ j+k\bigr)\pmod2,
\]

corresponding respectively to the leaf \(Y\), root \(X\), and root \(Z\).
The code uses four sectors.

### 9.1 Even-leaf sectors

The direct conditional squares are:

- `two_root_even_00`:
  \[
  i+j\equiv0,\quad i+k\equiv0,\quad j+k\equiv0\pmod2;
  \]
- `two_root_even_11`:
  \[
  i+j\equiv0,\quad i+k\equiv1,\quad j+k\equiv1\pmod2.
  \]

Each unlocalized degree-14 sector has 30 basis flags.

### 9.2 Odd-leaf orientation sectors

If \(i+j\) is odd, the raw conditional leaf integral vanishes by antipodal
symmetry. The code pairs it with the oriented tangent factor

\[
\det(X,Z,Y).
\]

After gluing, the multiplier is

\[
\det(X,Z,Y)\det(X,Z,W),
\]

which is independent of the choice of orientation after taking the product.
With

\[
\begin{aligned}
r&=X\cdot Z,&
p&=X\cdot Y,&
q&=X\cdot W,\\
s&=Z\cdot Y,&
t&=Z\cdot W,&
u&=Y\cdot W,
\end{aligned}
\]

the exact Gram polynomial is

\[
u-r^2u-st-pq+rpt+rqs.
\]

The two sectors are:

- `two_root_odd_01`:
  \[
  i+j\equiv1,\quad i+k\equiv0,\quad j+k\equiv1;
  \]
- `two_root_odd_10`:
  \[
  i+j\equiv1,\quad i+k\equiv1,\quad j+k\equiv0.
  \]

Each has 30 basis flags before localization.

### 9.3 Root-pair localizing blocks

Because \(X,Z\in S^2\),

\[
1-(X\cdot Z)^2\ge0.
\]

The command adds a second PSD block in every sector multiplied by this
nonnegative principal minor. The square basis now has degree at most \(6\).
The sizes are:

| block | size |
|---|---:|
| `two_root_even_00_minor` | 24 |
| `two_root_even_11_minor` | 20 |
| `two_root_odd_01_minor` | 20 |
| `two_root_odd_10_minor` | 20 |

For the odd sectors the localizing multiplier is

\[
\det(X,Z,Y)\det(X,Z,W)\bigl(1-(X\cdot Z)^2\bigr).
\]

## 10. Five-point root-weighted flags

The pair of options

```text
--max-flag-arity 5
--max-root-factor-degree 2
```

adds the systematic arity-five layer. There are three shared roots

\[
R=(X_0,X_1,X_2)
\]

and one leaf \(Y\) in each half of the square. A basis flag is indexed by:

- root-edge exponents
  \[
  \rho=(\rho_{01},\rho_{02},\rho_{12});
  \]
- root-to-leaf exponents
  \[
  \lambda=(\lambda_0,\lambda_1,\lambda_2).
  \]

It represents

\[
f_{\rho,\lambda}(R;Y)
=\prod_{0\le i<j\le2}(X_i\cdot X_j)^{\rho_{ij}}
 \prod_{i=0}^2(X_i\cdot Y)^{\lambda_i}.
\]

The restrictions are

\[
|\rho|\le2,
\qquad
|\rho|+|\lambda|\le7,
\qquad
|\lambda|\equiv0\pmod2.
\]

The first inequality is `--max-root-factor-degree 2`. The second ensures that
gluing two halves has degree at most \(14\). The last makes the leaf
projectively even.

The parity character at root \(X_i\) is

\[
\sigma_i
=\lambda_i+\sum_{j\ne i}\rho_{\min(i,j),\max(i,j)}
\pmod2.
\]

Only flags with equal \(\sigma\) can be paired. Because the total degree is
even, the active signatures and sizes are:

| block | root signature | size |
|---|---|---:|
| `weighted_flag_5_000` | \((0,0,0)\) | 92 |
| `weighted_flag_5_011` | \((0,1,1)\) | 80 |
| `weighted_flag_5_101` | \((1,0,1)\) | 80 |
| `weighted_flag_5_110` | \((1,1,0)\) | 80 |

For a row flag \((\rho,\lambda)\) and column flag
\((\rho',\lambda')\), gluing leaves \(Y,W\) produces

\[
\prod_{i<j}(X_i\cdot X_j)^{\rho_{ij}+\rho'_{ij}}
\prod_i(X_i\cdot Y)^{\lambda_i}
\prod_i(X_i\cdot W)^{\lambda'_i}.
\]

Reducing this five-vertex graph gives the corresponding matrix entry.

Since `--max-flag-arity` equals \(5\), no six-point block is added. Since the
command omits `--higher-rank-matrices`, it also does not add matrix-valued
Gram-determinant identities on the five-point flag spaces.

## 11. `--harmonics`: two-point spherical-harmonic positivity

The Legendre addition formula gives

\[
A_\ell
=\mathbb E[P_\ell(X\cdot Y)]
=\frac{4\pi}{2\ell+1}
  \sum_{m=-\ell}^{\ell}
  \left|
    \int Y_{\ell m}(x)\,d\mu(x)
  \right|^2
\ge0.
\]

The exact normalization depends on the convention for \(Y_{\ell m}\), but the
sign does not. Each \(A_\ell\) is therefore a \(1\times1\) PSD block.

At degree \(14\), the command adds

\[
\ell=4,6,8,10,12,14.
\]

Degree zero is the normalization \(A_0=1\). Under isotropy, the degree-two
harmonic vanishes, so it is not added as a separate positive block.

The polynomials used are:

\[
\begin{aligned}
P_4(t)
&=\frac{35t^4-30t^2+3}{8},\\
P_6(t)
&=\frac{231t^6-315t^4+105t^2-5}{16},\\
P_8(t)
&=\frac{
6435t^8-12012t^6+6930t^4-1260t^2+35
}{128},\\
P_{10}(t)
&=\frac{
46189t^{10}-109395t^8+90090t^6-30030t^4+3465t^2-63
}{256},\\
P_{12}(t)
&=\frac{
676039t^{12}-1939938t^{10}+2078505t^8
-1021020t^6+225225t^4-18018t^2+231
}{1024},\\
P_{14}(t)
&=\frac{
5014575t^{14}-16900975t^{12}+22309287t^{10}
-14549535t^8+4849845t^6-765765t^4+45045t^2-429
}{2048}.
\end{aligned}
\]

## 12. KKT conditions for a minimizing measure

Define the rooted potential

\[
U_\mu(x)=\int K(x\cdot y)\,d\mu(y).
\]

Then

\[
E_K(\mu)=\int U_\mu(x)\,d\mu(x).
\]

If \(\mu\) is a global minimizer, adding infinitesimal mass at \(z\) gives

\[
\left.\frac{d}{d\varepsilon}\right|_{\varepsilon=0}
E_K\bigl((1-\varepsilon)\mu+\varepsilon\delta_z\bigr)
=2(U_\mu(z)-E_K(\mu)).
\]

Therefore the global first-variation KKT conditions are

\[
\begin{aligned}
U_\mu(z)-E_K(\mu)&\ge0
&&\text{for every }z\in S^2,\\
U_\mu(x)-E_K(\mu)&=0
&&\text{for }\mu\text{-almost every }x.
\end{aligned}
\]

Moving a support point tangentially gives

\[
\nabla_{S^2}U_\mu(x)=0
\]

on the support, and second-order local minimality gives

\[
\operatorname{Hess}_{S^2}U_\mu(x)\succeq0
\]

on the tangent plane.

These KKT statements justify the next block families. They are not valid for
an arbitrary non-minimizing measure; the intended proof architecture applies
them to a hypothetical negative global minimizer.

The kernel derivatives used throughout are

\[
\begin{aligned}
K'(t)&=192t^5-192t^3+40t,\\
K''(t)&=960t^4-576t^2+40.
\end{aligned}
\]

## 13. `--potential`: scalar support-potential identities

For independent \(X,Y,Z,W\sim\mu\), support stationarity gives

\[
U_\mu(X)=E_K(\mu)
\qquad\mu\text{-almost surely}.
\]

Multiplying by \((X\cdot Z)^r\) and averaging gives

\[
\boxed{
\mathbb E[
  K(X\cdot Y)(X\cdot Z)^r
]
-
\mathbb E[K(X\cdot Y)]
\mathbb E[(Z\cdot W)^r]
=0.
}
\]

These are equality constraints, so their multipliers are unrestricted.

At degree \(14\), the code tests even powers

\[
r=0,2,4,6,8.
\]

After antipodal and isotropic reduction, the \(r=0,2\) relations are
identically zero. Three nonzero relations remain, for \(r=4,6,8\).

For example, the \(r=4\) relation reduces to

\[
\begin{aligned}
0={}&
-48M_{\triangle(0,4,4)}
+32M_{\triangle(0,4,6)}\\
&+48p_4^2-32p_4p_6.
\end{aligned}
\]

The \(r=6\) and \(r=8\) formulas are generated in the same way.

## 14. `--potential-matrices`: matrix support-potential identities

The scalar potential equalities are strengthened by multiplying
\(U_\mu(X)-E_K(\mu)\) by one-root flag Gram matrices.

Using the variables of Section 7, a typical entry is

\[
\begin{aligned}
0={}&
\mathbb E\!\left[
  K(X\cdot Y)
  (X\cdot Z)^r(X\cdot W)^s
  R_k(X;Z,W)
\right]\\
&-
E_K(\mu)\,
\mathbb E\!\left[
  (X\cdot Z)^r(X\cdot W)^s
  R_k(X;Z,W)
\right].
\end{aligned}
\]

In the isotropic branch,

\[
E_K=\frac{16}{3}-48p_4+32p_6,
\]

which is the expression used when expanding the second term.

Because these matrices equal zero rather than being nonnegative, their primal
matrix multipliers are unrestricted and the dual imposes matrix equalities.

The degree-14 bases are:

| free block | \(k\) | leaf powers | size |
|---|---:|---|---:|
| `potential_flag_0` | 0 | \(0,2,4\) | 3 |
| `potential_flag_1` | 1 | \(1,3\) | 2 |
| `potential_flag_2` | 2 | \(0,2\) | 2 |
| `potential_flag_3` | 3 | \(1\) | 1 |
| `potential_flag_4` | 4 | \(0\) | 1 |

The kernel consumes degree six, explaining why these bases are shorter than
the ordinary one-root flag bases.

## 15. `--gradient`: geometric first-order stationarity

The spherical gradient is

\[
G_\mu(X)
=\nabla_{S^2}U_\mu(X)
=\int
K'(X\cdot Y)
\bigl(Y-(X\cdot Y)X\bigr)\,d\mu(Y).
\]

For an auxiliary sample \(Z\), let

\[
a=X\cdot Y,\qquad
b=Y\cdot Z,\qquad
c=Z\cdot X.
\]

The tangent projection of \(Z\) is \(Z-cX\), so

\[
\bigl(Y-aX\bigr)\cdot(Z-cX)=b-ac.
\]

The scalar gradient kernel is therefore

\[
\boxed{g(a,b,c)=K'(a)(b-ac).}
\]

At a KKT support point \(X\), \(G_\mu(X)=0\). Multiplying its contraction with
the auxiliary tangent direction by \(c^r\) and averaging gives

\[
\mathbb E[c^rK'(a)(b-ac)]=0.
\]

The code considers \(r=0,\ldots,7\). Antipodal parity and isotropic
contraction make all but \(r=3,5,7\) vanish identically. The three retained
relations are unrestricted scalar equalities.

For instance, \(r=3\) reduces to

\[
\begin{aligned}
0={}&
-192M_{\triangle(0,4,6)}
+192M_{\triangle(1,3,5)}\\
&+192M_{\triangle(0,4,4)}
-192M_{\triangle(1,3,3)}.
\end{aligned}
\]

## 16. `--hessian`: scalarized spherical-Hessian positivity

For tangent vectors \(v,w\perp X\), the spherical Hessian is

\[
H_X(v,w)
=\int
\left[
K''(a)(v\cdot Y)(w\cdot Y)
-K'(a)a(v\cdot w)
\right]d\mu(Y),
\qquad a=X\cdot Y.
\]

At a locally minimizing support point,

\[
H_X\succeq0.
\]

Take the tangent vector

\[
v_\parallel=Z-cX,
\qquad c=X\cdot Z.
\]

Then

\[
v_\parallel\cdot Y=b-ac,
\qquad
\|v_\parallel\|^2=1-c^2,
\]

so the parallel Hessian kernel is

\[
\boxed{
h_\parallel(a,b,c)
=K''(a)(b-ac)^2-K'(a)a(1-c^2).
}
\]

For the perpendicular tangent vector

\[
v_\perp=X\times Z,
\]

one has

\[
(v_\perp\cdot Y)^2
=D(a,b,c),
\]

where

\[
D(a,b,c)
=1+2abc-a^2-b^2-c^2
=\det\operatorname{Gram}(X,Y,Z).
\]

Thus

\[
\boxed{
h_\perp(a,b,c)
=K''(a)D(a,b,c)-K'(a)a(1-c^2).
}
\]

The command integrates these nonnegative rooted Hessian quantities against a
nonnegative polynomial in \(c\). Every univariate polynomial nonnegative on
\([-1,1]\) can be represented in the form

\[
s(c)^{\mathsf T}Q_0s(c)
+(1-c^2)t(c)^{\mathsf T}Q_1t(c),
\qquad Q_0,Q_1\succeq0,
\]

at an appropriate degree.

For degree \(14\):

- \(s(c)=(1,c,c^2,c^3)\), producing \(4\times4\)
  `hessian_sos` and `perpendicular_hessian_sos` blocks;
- \(t(c)=(1,c,c^2)\), producing \(3\times3\)
  `hessian_minor` and `perpendicular_hessian_minor` blocks multiplied by
  \(1-c^2\).

All four matrices are PSD variables in the primal and PSD moment constraints
in the dual.

## 17. `--four-point-hessian`: bilinear tangent-field Hessian blocks

Scalarizing the Hessian along one tangent vector does not use its full matrix
structure. The four-point blocks test bilinear tangent vector fields.

Let

\[
\begin{aligned}
a&=X\cdot Y,&
c&=X\cdot Z,&
d&=X\cdot W,\\
b&=Y\cdot Z,&
e&=Y\cdot W,&
f&=Z\cdot W.
\end{aligned}
\]

For

\[
v_Z=Z-cX,\qquad v_W=W-dX,
\]

the bilinear Hessian kernel is

\[
\boxed{
H_\parallel
=K''(a)(b-ac)(e-ad)
-K'(a)a(f-cd).
}
\]

For the rotated tangent fields

\[
Jv_Z=X\times Z,\qquad Jv_W=X\times W,
\]

the identity

\[
(Jv_Z\cdot Y)(Jv_W\cdot Y)
=(f-cd)(1-a^2)-(b-ac)(e-ad)
\]

gives

\[
\boxed{
H_\perp
=K''(a)
  \left[(f-cd)(1-a^2)-(b-ac)(e-ad)\right]
-K'(a)a(f-cd).
}
\]

The vector-field basis uses odd auxiliary powers

\[
c^r v_Z,\qquad r\in\{1,3\}.
\]

Odd powers make the complete vector field invariant under the projective sign
change \(Z\mapsto -Z\). Consequently, the command creates two \(2\times2\)
PSD blocks:

```text
four_point_parallel_hessian
four_point_perpendicular_hessian
```

The command does not specify `--max-hessian-arity`, so no five-point Hessian
flag blocks are included.

## 18. `--global-gap`: the averaged global first-variation gap

The KKT inequality

\[
q_\mu(z)=U_\mu(z)-E_K(\mu)\ge0
\]

holds for every trial point \(z\) at a global minimizer.

Average \(z\) uniformly over \(S^2\). The degree-zero Legendre coefficient of
\(K\) is

\[
\int_{S^2}K(z\cdot y)\,d\sigma(z)=\frac{32}{105}.
\]

Therefore

\[
\frac{32}{105}-E_K(\mu)\ge0.
\]

Using isotropy,

\[
\frac{32}{105}-E_K
=-\frac{176}{35}+48p_4-32p_6.
\]

This is the \(1\times1\) PSD block `global_uniform_gap`.

## 19. `--global-tangent-gaps`: explicit trial-point gaps

These blocks evaluate the global KKT gap at two trial points constructed from
a support point \(X\) and an auxiliary sample \(Z\).

Use

\[
a=X\cdot Y,\qquad b=Y\cdot Z,\qquad c=Z\cdot X,
\qquad d=1-c^2.
\]

### 19.1 Parallel tangent trial

For \(d>0\), define

\[
z_\parallel=\frac{Z-cX}{\sqrt d}.
\]

Then

\[
z_\parallel\cdot Y=\frac{b-ac}{\sqrt d}.
\]

Because \(U_\mu(X)=E_K\) on the support, the cleared-denominator gap is

\[
\boxed{
G_\parallel(a,b,c)
=d^3
\left[
K\!\left(\frac{b-ac}{\sqrt d}\right)-K(a)
\right].
}
\]

Since \(K\) is even and of degree six, this expression is a polynomial:

\[
\begin{aligned}
G_\parallel
={}&32(b-ac)^6
-48(b-ac)^4d
+20(b-ac)^2d^2\\
&-\frac43d^3-d^3K(a).
\end{aligned}
\]

### 19.2 Perpendicular tangent trial

Define

\[
z_\perp=\frac{X\times Z}{\sqrt d}.
\]

The squared inner product with \(Y\) is

\[
(z_\perp\cdot Y)^2
=\frac{D(a,b,c)}{d},
\]

where \(D\) is the three-point Gram determinant. The polynomial gap is

\[
\boxed{
\begin{aligned}
G_\perp
={}&32D^3-48D^2d+20Dd^2\\
&-\frac43d^3-d^3K(a).
\end{aligned}
}
\]

Both gaps are nonnegative after averaging over \(Y\).

The base gap degree budget is 12. At total degree 14, the ordinary multiplier
basis is

\[
(1,c),
\]

giving a \(2\times2\) PSD block for each direction. A second localizing block
uses the multiplier \(1-c^2\) with constant basis, giving a \(1\times1\) block
for each direction:

```text
global_parallel_tangent_gap          size 2
global_parallel_tangent_gap_minor    size 1
global_perpendicular_tangent_gap     size 2
global_perpendicular_tangent_gap_minor size 1
```

The formulas extend continuously through \(d=0\) because the code uses the
cleared polynomial expressions.

## 20. `--rank-relations`: exact dimension-three Gram identities

Four vectors in \(\mathbb R^3\) are linearly dependent. Hence their
\(4\times4\) Gram determinant is identically zero:

\[
\det
\begin{pmatrix}
1&g_{01}&g_{02}&g_{03}\\
g_{01}&1&g_{12}&g_{13}\\
g_{02}&g_{12}&1&g_{23}\\
g_{03}&g_{13}&g_{23}&1
\end{pmatrix}
=0.
\]

The determinant has Gram degree at most four. To remain within total degree
14, the code multiplies it by canonical four-vertex Gram monomials of degree
at most

\[
14-4=10.
\]

Only multipliers with even degree at every vertex survive projective parity.
The results are:

1. reduced using antipodality and isotropy;
2. canonicalized under all four-vertex permutations;
3. normalized and deduplicated.

This produces 61 independent stored scalar relations for the stated command.
Their primal coefficients are unrestricted, and the dual requires every
formal pairing with the determinant relation to vanish.

These relations are crucial: the formal moment labels otherwise forget that
the sampled vectors live in dimension three.

## 21. Complete degree-14 block inventory

The command instantiates the following objects before any optional facial
reduction.

### 21.1 Positive semidefinite blocks

| family | number of blocks | sizes |
|---|---:|---|
| one-root harmonic flags | 7 | \(4,3,3,2,2,1,1\) |
| empty-type four-point flag | 1 | \(4\) |
| two-root flags | 8 | \(30,30,30,30,24,20,20,20\) |
| five-point root-weighted flags | 4 | \(92,80,80,80\) |
| Legendre harmonics | 6 | six \(1\times1\) blocks |
| uniform global gap | 1 | \(1\) |
| tangent global gaps | 4 | \(2,1,2,1\) |
| scalar Hessian localizers | 4 | \(4,4,3,3\) |
| bilinear four-point Hessians | 2 | \(2,2\) |

There are 37 PSD blocks in total.

### 21.2 Matrix equality blocks

The five potential-stationarity flag matrices have sizes

\[
3,2,2,1,1.
\]

They are equality blocks, not PSD blocks.

### 21.3 Scalar equalities

| family | count |
|---|---:|
| gradient stationarity | 3 |
| scalar potential stationarity | 3 |
| four-vector rank identities | 61 |

There are 67 scalar equalities.

### 21.4 Formal variables

After collecting all labels used by the target and the selected blocks, the
dual has 574 formal moment variables, one of which is fixed by

\[
y_{\mathrm{constant}}=1.
\]

## 22. Solver and output options

### 22.1 `--dual`

Builds the formal moment minimization problem rather than the primal
coefficient-decomposition problem. The dual is numerically more stable on the
sharp equality face, but it does not itself provide certificate matrices.

### 22.2 `--scale-constraints`

For each PSD or equality matrix, let

\[
s=\max_L\|A_L\|_{\max}.
\]

The code replaces the matrix by

\[
\frac1s\sum_Ly_LA_L.
\]

It similarly divides each scalar identity by its largest coefficient. Since
\(s>0\), this changes neither PSD signs nor equality solution sets. It only
changes numerical conditioning.

### 22.3 `--summary-only`

Suppresses the potentially large dictionary of nonzero formal moments. The
JSON output still reports:

- solver status;
- objective;
- number of labels;
- minimum PSD-block eigenvalue;
- maximum matrix-equality residual;
- maximum scalar-relation residual.

### 22.4 `--tolerance 1e-10`

Sets MOSEK's primal-feasibility, dual-feasibility, and relative-gap tolerances
to \(10^{-10}\). This is a requested numerical termination threshold, not an
exact-arithmetic guarantee.

## 23. How to interpret the reported objective

For the published run, MOSEK reports a value close to

\[
-4.8\times10^{-8}.
\]

The correct interpretation is:

- the selected finite formal-moment relaxation nearly proves \(T\ge0\);
- its small negative optimum may be numerical error, a relaxation gap, or a
  non-attainment/closure phenomenon;
- it is not an actual measure of negative energy;
- it is not an exact SOS proof.

If the dual optimum were rigorously nonnegative, then every actual isotropic
measure satisfying the encoded KKT conditions would have \(T\ge0\). If a
rational primal identity with exact PSD matrices were recovered, the identity
would give a machine-checkable certificate for that same isotropic KKT
branch.

To turn this branch result into unrestricted copositivity, one still needs
either:

1. a proof that a hypothetical negative global minimizer can be assumed
   isotropic; or
2. a certificate hierarchy that retains the full second-moment matrix rather
   than imposing isotropy.

## 24. Direct map from command options to implementation

| option | main implementation |
|---|---|
| `--degree 14` | degree budgets inside `solve` |
| `--no-pointwise-sos` | skips the `sos` module term |
| `--three-point-flags` | `tangent_harmonic_polynomials`, `flag_expectation_matrix` |
| `--four-point-flags` | `empty_type_flag_expectation_matrix` |
| `--two-root-flags` | `two_root_flag_expectation_matrix` |
| `--max-flag-arity 5` | `rooted_weighted_flag_sectors`, `rooted_weighted_flag_expectation_matrix` |
| `--max-root-factor-degree 2` | root-edge degree bound in the weighted flag basis |
| `--harmonics` | `harmonic_pair_vector` |
| `--potential` | `potential_stationarity_relation` |
| `--potential-matrices` | `potential_flag_relation_matrix` |
| `--gradient` | gradient polynomial from `kernel_polynomials` |
| `--hessian` | scalar kernels from `kernel_polynomials` |
| `--four-point-hessian` | `four_point_hessian_polynomials`, `four_point_hessian_expectation_matrix` |
| `--global-gap` | `global_uniform_gap` block in `solve` |
| `--global-tangent-gaps` | `global_tangent_gap_polynomials` |
| `--rank-relations` | `four_point_rank_relations` |
| `--dual` | formal moment branch in `solve` |
| `--scale-constraints` | positive block/relation rescaling in the dual branch |
| `--summary-only` | omits the formal moment dictionary from JSON |
| `--tolerance 1e-10` | MOSEK interior-point tolerance parameters |

The graph reducer used by all these families is
`reduce_graph_matrix`, with canonical cached entry point
`graph_expectation_label`.
