# The unprojected escape, fingerprinted

*2026-08-17.  PLAN Next-actions #1 executed on regenerated data (new $E$
normalization throughout).  Tools: `sdpa_runs/fingerprint_blocks.py`,
`sdpa_runs/fingerprint_expand.py`; data files under `sdpa_runs/`:
`sel_h2all_{1em3,1em4,5em5}.result`, `fingerprint_norms_e3e4.json`,
`fingerprint_D_e3e4.json`, blocks dump `blocks_deg14_h2w_h2all.json`.*

## 1. Setup and trace law

Minimal-trace selector certificates for $h_2E\ge-\varepsilon$ on the
degree-14 KKT-inclusive weighted problem (`deg14_h2w_h2all.dat-s`,
m = 790, bound wall $-2.7908758\times10^{-5}$, re-verified this date):

| $\varepsilon$ (new scale) | $\operatorname{tr}C_w$ | phase |
|---|---:|---|
| $10^{-3}$ | $764.99$ | pdFEAS |
| $10^{-4}$ | $8133.57$ | pdFEAS |

Growth $10.63\times$ per decade — local exponent $\varepsilon^{-1.03}$,
a clean simple pole on a decade far from the wall ($36\times$–$3.6\times$
above it).  The historical $16.7\times$ (legacy points) was
wall-inflated; the pole itself is simple, confirming a *single*
multiplier deficit rather than a higher-order one.

## 2. Where the growth lives (block level)

Frobenius growth $\varepsilon=10^{-3}\to10^{-4}$, absolute-mass
carriers:

| block | frob@$10^{-3}$ | frob@$10^{-4}$ | ratio | trace@$10^{-4}$ |
|---|---:|---:|---:|---:|
| `h2loc_empty_type_flag` | $5.1\times10^{-8}$ | $517.9$ | $10^{10}$ | $517.9$ |
| `h2loc_two_root_even_11` | $23.8$ | $1438.2$ | $60.5$ | $2849.6$ |
| `h2loc_two_root_odd_10/01` | $21.6$ | $391.0$ | $18.1$ | $745.9$ each |
| `h2loc_two_root_odd_*_minor` | $3.6$ | $126.4$ | $35.1$ | $186.0$ each |
| `h2loc_two_root_even_00` | $261.4$ | $1071.3$ | $4.1$ | $1769.5$ |
| `h2loc_two_root_even_11_minor` | $42.3$ | $214.6$ | $5.1$ | $407.1$ |

All other blocks are static or numerically negligible.  The escape is
carried by the **$h_2$-localized two-root modulated blocks** plus an
explosive activation of the **$h_2$-localized product block**
(`h2loc_empty_type_flag`, from machine zero).

## 3. Label expansion of the growth direction (the verdict)

$D=e(Y_{10^{-4}})-e(Y_{10^{-3}})$ via the gauge-invariant fingerprint
$e_L(Y)=\sum_b\langle A^b_L,Y_b\rangle$.  Sector masses of $|D|$
(total $5.0\times10^{4}$):

| sector | share |
|---|---:|
| $p_2\times$`graph_4` products | $57.7\%$ |
| `graph_4` connected | $19.2\%$ |
| $p_2\times$pair-products ($p_2p_4p_6$-type) | $10.8\%$ |
| $p_2\times$`triangle` products | $9.1\%$ |
| `triangle` connected | $3.2\%$ |
| **pair sector ($p_k$)** | $0.07\%$ |

Dominant individual labels: $p_2\times$(`graph_4` with edge
multiplicities up to 5), $p_2\times$(`triangle` $(1,3,5)$, $(0,4,6)$,
$(1,5,5)$, $(1,3,3)$), $p_2\times(p_2p_4p_6,\ p_2p_4p_4,\ p_2p_2p_4)$,
both signs.  $\langle\text{target},D\rangle=-3.64\times10^{4}<0$
(improving direction, as required).

**Verdict.**  The surviving escape of the weighted problem is the
**$h_2$-localized collision/theta tower**: $h_2\times$(high-modulation
two-root content), with essentially *zero* pair-sector component.
Consequences by the sign rule ([Gap cuts note](GAP_CUTS_NOTE.md)):

1. Every scalar invariant expanding in pair labels alone — powers of
   $h_2$, the $e_k(I-A_2)$ scalar-gap cuts, any $g_\ell$ combination —
   pairs $\approx0$ with $D$: dead as weights *and* as cuts against
   this escape.  This independently re-derives and extends the PLAN §4
   conclusion.
2. Matched objects must load on the product/collision sector.  The two
   implementable candidate classes, in order:
   (a) **$h_2$-localized theta atoms** — the diagonal window cuts of
   [Theta atom note](THETA_ATOM_NOTE.md) applied to the
   $h_2\times$two-root modulated family (`h2loc_two_root_*`), i.e.
   atoms $\sum_nq^{n^2}\,h_2\!\cdot\!\mathcal Q[\hat G_n]$ with the
   same rational tail majorants (validity: $h_2\ge0$ times a valid
   square family; the tail bound needs $h_2\le1$, which holds since
   $h_2=\frac{3p_2-1}2\le1$);
   (b) the **arity-2 weighted-(E1) projection** of the
   `h2loc_two_root` layer (in progress), which must annihilate the
   sharp-face component of exactly these blocks.
3. The $\approx20\%$ connected `graph_4` component is the natural
   4-sample shadow of the quadrupole-word tower
   ([Multi-weight program](MULTI_WEIGHT_PROGRAM.md)); whether the
   `graph_5` extension is required is decided by whether (a)+(b) kill
   the pole without it.

**Cross-check (measured).**  Third selector point
$\varepsilon=5.33\times10^{-5}$: $\operatorname{tr}=13637.2$ (pFEAS,
wall-softened as expected at $1.9\times$ the wall).  The growth
direction on $10^{-4}\to5.33\times10^{-5}$
(`fingerprint_D_e4e5.json`) has sector masses $64.5\%$
$p_2\times$`graph_4`, $21.4\%$ `graph_4`, $9.1\%$ $p_2\times$`triangle`,
$2.5\%$ $p_2\times$pair-products, $0.0\%$ pair — the same sectors and
the same dominant labels (top label
$p_2\times(\texttt{graph\_4},0,1,3,1,3,2)$ in both intervals).  The
fingerprint is stable across $\varepsilon$, hence not a wall artifact.
The degree-16 selector after the queued dual solve remains the fully
wall-free confirmation.

## 4. Rounds two and three: the Jensen blocks bend, then halve, the pole (2026-08-17, late)

Selector traces on the weighted degree-14 problem (new scale):

| ε | control | +Jensen v1 (11 families) | +v2 (29 families) |
|---|---:|---:|---:|
| 1e−3 | 764.99 | 733.47 | **437.51** |
| 1e−4 | 8133.57 | 7034.96 | **1544.18** |
| 5.33e−5 | 13637.2 | — | **1909.11** |
| growth/decade | 10.63× | 9.59× | **3.53×** |
| local exponent | 1.03 | 0.98 | **0.55** (→0.34 below 1e−4) |

The conditional-covariance (Jensen/unfolding) families are the first
objects to move this law at all; the v2 set halves the pole order in
one design iteration and keeps flattening at smaller ε.

**Residual escape (v2), fingerprinted** (`fingerprint_norms_toep2.json`,
`fingerprint_D_toep2.json`): sectors rotated — p2×graph_4 + graph_4
collapsed from 77% to 3.8%; the residual is 49.7% p2×triangle,
29.3% p2×pair-products, 17.1% triangle (top labels
p2×triangle(1,3,5), (0,4,6), (1,5,5); block carriers
h2loc_two_root_even_00, h2loc_flag_2/3).  The certificate now pumps
the same-leaf (T) side of the Jensen inequality; the matched v3
object is the **conditional-Toeplitz moment-matrix positivity of the
fiber** (PSD by Toeplitz positivity, not by squaring — the flag
incarnation of the 7×7 fiber PSD of the cylindrical analysis), plus
its p2-product/h2loc copies.  v3 in progress.
