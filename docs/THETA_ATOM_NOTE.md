# Theta-atom adjunction: entire-kernel generators against the $\varepsilon\to0$ pole

Companion to [PLAN §6](../PLAN.md) (open item "exact $\varepsilon=0$
feasibility over the projected cone + theta atoms"),
[Two-root generators](TWO_ROOT_GENERATORS.md) (the modulated families
$\hat C_n$, $\hat S_n$), [Wrapper lemmas](WRAPPER_LEMMAS.md) L7 (the
generator bounds), and [Gap cuts note](GAP_CUTS_NOTE.md) (the
ray-pairing sign rule).  Everything in §1–§3 is exact mathematics;
measurements are recorded in §4–§6.

Status: formalization (§1–§3) complete; pairing measurement (§4), A/B
(§5), $\varepsilon=0$ assembly (§6), verdict (§7), and the same-day
weighted-escape response (§8: $h_2$-localized atoms, sup-caps, GMP
export) recorded below.  All 2026-08-17, $E$ normalization.

## 0. Setting and notation

Labels $\ell$ and pseudo-moment vectors $y=(y_\ell)$ are as in
[Implementation](IMPLEMENTATION.md); for a genuine antipodal
probability measure $\mu$, $y_\ell(\mu)$ is the corresponding
multigraph moment.  For a polynomial $A(t_1,t_2,s)$ the **two-root
quadratic form** is

$$\mathcal Q[A](\mu)\;=\;\iint\Big(\int A(x_1\cdot y,\;x_2\cdot y,\;x_1\cdot x_2)\,d\mu(y)\Big)^{2}\,d\mu(x_1)\,d\mu(x_2)\;\ge\;0 ,$$

a $1\times1$ two-root Gram block.  Expanding the square gives the exact
rational label expansion
$\mathcal Q[A](\mu)=\sum_\ell q^{A}_\ell\,y_\ell(\mu)$ with
$q^A\in\mathbb Q^{(\text{4-point labels})}$: with leaves $y,y'$ and
roots $x_1,x_2$,

$$\mathcal Q[A]=E\big[A(x_1\!\cdot\! y,\,x_2\!\cdot\! y,\,x_1\!\cdot\! x_2)\;A(x_1\!\cdot\! y',\,x_2\!\cdot\! y',\,x_1\!\cdot\! x_2)\big],$$

a linear combination of four-vertex multigraph labels of degree
$2\deg A$.  The modulated generators are those of
[Two-root generators](TWO_ROOT_GENERATORS.md):

$$\hat C_n=C_n+\tfrac13T_{|n|}(s),\qquad \hat S_n=S_n+\tfrac13T_{|n|}(s),
\qquad n\in\mathbb Z,$$

$\deg\hat C_n=|n|+1$ ($n\ne0$), $\deg\hat S_n=|n|+5$ ($n\ne0$), with
the machine-checked closed forms and Gram-body bounds
([Wrapper lemmas](WRAPPER_LEMMAS.md) L7)

$$|\hat C_n|\le\tfrac43+4|n|,\qquad |\hat S_n|\le\tfrac43+12|n|
\qquad\text{on }|t_1|,|t_2|,|s|\le1 .$$

Write $\hat G^{(2)}_n=\hat C_n$, $\hat G^{(6)}_n=\hat S_n$, and
$c_2=4$, $c_6=12$ for the slope constants, $a=\tfrac43$.

## 1. What object is adjoined, and on which side

**Definition (theta atom, diagonal/incoherent form).**  For a family
$f\in\{2,6\}$ and a rational $q\in(0,1)$, the *theta-atom moment* of a
measure $\mu$ is the number

$$\tau_{f,q}(\mu)\;:=\;\sum_{n\in\mathbb Z}q^{n^2}\,\mathcal Q[\hat G^{(f)}_n](\mu)\;\in\;[0,\;T^{(f)}_q(0)+\mathcal Q[\hat G^{(f)}_0](\mu)] ,$$

absolutely convergent because
$0\le\mathcal Q[\hat G^{(f)}_n]\le(a+c_f|n|)^2$ (the leaf integral of a
probability measure is bounded by the sup of the leaf, L7).

**Where it enters.**  On the **moment side** (the relaxation the solver
actually solves), each atom adjoins one new scalar pseudo-moment
variable $y_\tau$ — a new label
$(\texttt{theta},f,q)$ — linked to the polynomial labels by one-sided
linear cuts (§2).  Dually, on the **certificate side**, the cut
multipliers assemble into a new certificate generator: the resummed
series $\sum_n q^{n^2}\mathcal Q[\hat G^{(f)}_n]$ may appear in a
certificate with *either sign* on its high-modulation part, the
negative part being paid for by an exactly rational constant (the tail
majorant).  No closed-form value of $\tau_{f,q}$ is ever needed: the
atom exists in the SDP only through its cuts, exactly as the kernel
labels $k_G$ of [Exact zero program](EXACT_ZERO_PROGRAM.md) §4.2.

**Why the diagonal (not the coherent square).**  The alternative
"coherent" atom — the single leaf
$\Theta_q=\sum_nq^{n^2}\hat G_n$ used inside one square
$\mathcal Q[\Theta_q]$ — expands as
$\sum_{m,n}q^{m^2+n^2}\,\langle\hat G_m,\hat G_n\rangle_\mu$ with
cross terms of both signs; its truncation error is controlled only by
Cauchy–Schwarz windows of width $\sqrt{T}$, which are not rational and
not one-sided.  The diagonal atom has *exactly rational, one-sided*
truncation control (§2), which is what the exact-certificate program
requires.  No strength is lost: by Cauchy–Schwarz on the leaf
integrals, at every measure
$$\mathcal Q\Big[\textstyle\sum_nq^{n^2}\hat G_n\Big]\;\le\;
\Big(\textstyle\sum_nq^{n^2}\Big)\cdot\sum_nq^{n^2}\,\mathcal Q[\hat G_n]
\;=\;\theta_3(q)\,\tau_{f,q}(\mu),$$
so every coherent-square certificate term is dominated by the diagonal
atom up to the explicit rational-majorizable constant
$\theta_3(q)=\sum_nq^{n^2}$, and the escape diagnostics (PLAN §6)
concern exactly the diagonal profile $q^{n^2}$.
(The coherent object stays available as a *finite* leaf
$\sum_{|n|\le N}q^{n^2}\hat G_n$ inside the ordinary projected two-root
blocks — it is polynomial, hence already in the hierarchy.)

## 2. The cuts, and their validity for every measure

Fix $f$, $q$, and abbreviate
$L_n(y)=\sum_\ell q^{\hat G_n}_\ell\,y_\ell$ (the exact label expansion
of $\mathcal Q[\hat G^{(f)}_n]$), so $L_n(y(\mu))=\mathcal Q[\hat
G^{(f)}_n](\mu)$.  The adjunction consists of, for a finite set of
truncation orders $N$:

* **(G$_n$) generator blocks** ($1\times1$, valid squares):
  $L_n(y)\ \ge\ 0$ for each carried $|n|\le N_{\max}$.  For $|n|$ inside
  the (E1)-projected two-root basis these are implied by the existing
  projected Gram blocks (each $\hat G_n$ is an admissible basis
  direction); beyond the solver degree they are new valid $1\times1$
  blocks (their labels are new but are still genuine moments).

* **(T$^-_N$) lower truncation cuts**:
  $$y_\tau\;-\;\sum_{|n|\le N}q^{n^2}L_n(y)\;\ge\;0 .$$
  Valid at every measure because the tail
  $\sum_{|n|>N}q^{n^2}\mathcal Q[\hat G_n](\mu)$ is a sum of
  nonnegative terms — "each tail is a sum of admissible PSD blocks"
  (PLAN §6).

* **(T$^+_N$) upper truncation cuts**:
  $$\sum_{|n|\le N}q^{n^2}L_n(y)\;+\;T^{(f)}_q(N)\;-\;y_\tau\;\ge\;0 ,$$
  with the exactly rational tail majorant $T^{(f)}_q(N)$ of §2.1.
  Valid at every measure because
  $\mathcal Q[\hat G_n](\mu)\le(a+c_f|n|)^2$ pointwise (L7 bound +
  Cauchy–Schwarz for the probability leaf integral).

**Validity lemma.**  *For every antipodal probability measure $\mu$,
the assignment $y=y(\mu)$, $y_\tau=\tau_{f,q}(\mu)$ satisfies (G$_n$),
(T$^-_N$), (T$^+_N$) for all $n$, $N$.  Hence the extended problem is
still a relaxation of $\inf_\mu E$, its optimum is still a valid lower
bound, and — since every constraint is an identity or a valid
inequality for* all *measures (no KKT content) — the adjunction
composes with the reduction lemma and with the all-measures certificate
track (L5/L6(i)).*

Proof: (G$_n$) is $\mathcal Q[\hat G_n](\mu)\ge0$, a square.
(T$^-_N$): the omitted tail is a sum of squares times positive
$q^{n^2}$.  (T$^+_N$): each omitted term is at most
$q^{n^2}(a+c_f|n|)^2$ and §2.1 majorizes the sum.  $\square$

All data are rational: $q^{n^2}\in\mathbb Q$, the $q^{\hat G_n}_\ell$
are exact rationals (Chebyshev recurrence with rational seeds), and
$T^{(f)}_q(N)\in\mathbb Q$.  The construction is certification-ready.

### 2.1 The exactly rational tail majorants

For rational $q\in(0,1)$, integers $N\ge0$, and the geometric
domination $q^{n^2}\le q^{N^2}\,q^{(2N+1)(n-N)}$ for $n>N$ (because
$n^2-N^2=(n-N)(n+N)\ge(n-N)(2N+1)$ when $n-N\ge1$), put
$r=q^{2N+1}\in(0,1)$, $\alpha=a+c_fN$, $m=n-N$.  Using
$\sum_{m\ge1}r^m=\frac r{1-r}$,
$\sum_{m\ge1}mr^m=\frac r{(1-r)^2}$,
$\sum_{m\ge1}m^2r^m=\frac{r(1+r)}{(1-r)^3}$:

**First-moment tail** (for linear atom variants and diagnostics):

$$\sum_{n>N}q^{n^2}(a+c_fn)\;\le\;
q^{N^2}\Big[\alpha\,\frac{r}{1-r}+c_f\,\frac{r}{(1-r)^2}\Big]
\;=:\;t^{(f)}_q(N)\;\in\mathbb Q .$$

**Second-moment tail** (the one the atom cuts use):

$$\sum_{n>N}q^{n^2}(a+c_fn)^2\;\le\;
q^{N^2}\Big[\alpha^2\,\frac{r}{1-r}
+2\alpha c_f\,\frac{r}{(1-r)^2}
+c_f^2\,\frac{r(1+r)}{(1-r)^3}\Big] ,$$

and the two-sided-in-$n$ majorant used in (T$^+_N$) is

$$T^{(f)}_q(N)\;=\;2\,q^{N^2}\Big[\alpha^2\,\frac{r}{1-r}
+2\alpha c_f\,\frac{r}{(1-r)^2}
+c_f^2\,\frac{r(1+r)}{(1-r)^3}\Big]\;\in\;\mathbb Q ,$$

valid because $|\hat G^{(f)}_{-n}|\le a+c_f n$ as well (the closed
form $C_n=C_0T_n(s)+2t_1(t_2-st_1)U_{n-1}(s)$ extends to $n<0$ with
$T_{-n}=T_n$, $U_{-n}=-U_{n-2}$, so the same triangle bound applies;
machine-checked in `theta_atoms.py --self-test` alongside the exact
recurrence/closed-form agreement).

Example values (`theta_atoms.py --majorants` prints exact fractions;
machine-checked against 400-term 60-digit partial sums): at
$q=\tfrac12$, $T^{(2)}_q(N)=238,\ 14.4,\ 0.733,\ 9.28\times10^{-3},\
2.72\times10^{-5},\ 1.87\times10^{-8}$ for $N=0..5$ and
$T^{(6)}_q(N)=1860,\ 108,\ 5.76,\ 7.52\times10^{-2},\
2.25\times10^{-4},\ 1.57\times10^{-7}$ — the constants collapse
super-exponentially in $N$, which is what makes deep-window atom cuts
numerically tight.

## 3. What the adjunction actually adds (elimination of $y_\tau$)

$y_\tau$ occurs only in (T$^\pm$).  Eliminating it, the extended
relaxation's shadow on the polynomial labels is: for every carried pair
$N'<N$,

$$\boxed{\;\sum_{N'<|n|\le N}q^{n^2}\,L_n(y)\;\le\;T^{(f)}_q(N')\;}
\tag{W$_{N',N}$}$$

— take (T$^+_{N'}$) and (T$^-_{N}$).  (Pairs $N\le N'$ give
consequences of (G$_n$) alone, i.e. nothing new; and any $y$
satisfying all (W$_{N',N}$) with the (G$_n$) lifts back to a feasible
$y_\tau$, e.g. $y_\tau=\max_N\sum_{|n|\le N}q^{n^2}L_n(y)$, so the
boxed window inequalities are *exactly* the adjunction's content.)

(W$_{N',N}$) is a genuinely new valid inequality: a uniform,
$q$-graded bound on the total high-modulation diagonal two-root mass —
the polynomial relaxation at any finite degree contains nothing that
couples the modulation axis to a constant.  This is the
finite-dimensional shadow of "the certificate is an entire function of
theta type in the modulation index".

**Recession-cone consequence (the pole mechanism).**  A recession
direction $r$ of the extended problem must satisfy the homogenized
cuts (constants drop):

$$L_n(r)\ \ge\ 0\quad(|n|\le N_{\max}),\qquad
\sum_{N'<|n|\le N}q^{n^2}\,L_n(r)\ \le\ 0
\;\Longrightarrow\;
L_n(r)=0\quad\text{for all }N'<|n|\le N .$$

So the atom kills every escape ray that carries *any* strictly positive
diagonal $\hat C_n$/$\hat S_n$ mass in the window $N'<|n|\le N$ — and
only those.  By the sign rule of [Gap cuts note](GAP_CUTS_NOTE.md),
the atom-induced constraint (W$_{N',N}$) pairs against a ray $r$ as
$-\sum_{N'<|n|\le N}q^{n^2}L_n(r)$: **negative pairing (= positive
window mass) kills the ray; zero pairing for every $q$, $f$, window
means the theta ansatz misses the surviving escape.**  This is the
cheap decisive test of §4.

**Effective window at solver degree $d$.**  The top-degree
cancellation in the closed form works only forward:
$\deg\hat G_n=n+\mathrm{off}_f$ for $n\ge1$ but $|n|+\mathrm{off}_f+1$
for $n\le-1$, $\mathrm{off}_2=1$, $\mathrm{off}_6=5$ (verified exactly
in `theta_atoms.py --self-test`; the docs' symmetric statement holds
"up to seed offsets").  $L_n$ involves labels of degree $2\deg\hat
G_n$, so within the degree-$d$ label algebra the atom windows reach
$-(d/2-\mathrm{off}_f-1)\le n\le d/2-\mathrm{off}_f$; at $d=14$:
$n\in[-5,6]$ (family 2) and $n\in[-1,2]$ (family 6).  Larger $|n|$ may
be carried too — their labels are new but still genuine moments — but
their coordinates are constrained only by (G$_n$) and the cuts, so a
solver can zero them on a ray; the *load-bearing* window is the
within-degree one.

**Certificate-side reading.**  A dual solution assigns multipliers
$\lambda_N\ge0$ to (T$^-_N$) and $\rho_{N'}\ge0$ to (T$^+_{N'}$) with
$\sum\lambda_N=\sum\rho_{N'}$ (the $y_\tau$ coefficient must cancel),
i.e. certificates gain the terms

$$\sum_{n}q^{n^2}\Big[\underbrace{\textstyle\sum_{N'\ge|n|}\rho_{N'}-\sum_{N\ge|n|}\lambda_N}_{=:c_n,\ \text{may be negative}}\Big]\,\mathcal Q[\hat G_n](\mu)
\;+\;\sum_{N'}\rho_{N'}T^{(f)}_q(N') ,$$

nonnegative for every $\mu$ by construction.  The atom therefore lets a
certificate *subtract* high-modulation squares at rate $q^{n^2}$,
paying an exact rational constant — precisely the resummation the
escaping $\varepsilon$-certificates were trying to build term by term.

## 4. Ray-pairing measurement (C2)

Measured 2026-08-17 with `theta_atoms.py --pair` against the
regenerated recession ray of the (E1)-projected degree-14 dual
(`sdpa_runs/ray_nocuts_deg14.json`; pruned nine-toggle base +
`--e1-project e1_projection_deg14.json --scale-constraints --find-ray`,
CLARABEL `optimal_inaccurate`).  The regenerated ray *is* the
documented one, now in the $E$ normalization: squared norm
$1.9586=55.7\times(3/16)^2$, $\langle E,r\rangle=-1.0000$,
$g_4\cdot r=1.0031$ ($\times16/3=5.3499$, documented $5.3498$),
$g_2\cdot r=4.4\times10^{-4}$, $g_{6..8}\cdot r=0$.  All eight
`theta_atoms.py --self-test` checks pass (closed form $=$ recurrence
on $n\in[-8,8]$, pole-value lemma, sup bounds, majorant domination,
label expansion $=$ direct integration to $5\times10^{-11}$, sandwich
validity at random measures).

**Per-generator diagonal pairings**
$\ell_n(r)=\langle\mathcal Q[\hat G_n],r\rangle$ (all labels of every
expansion are present in the 642-label problem; the only coordinate
outside the stored support is the constant, which is $0$ by the ray
normalization):

| family | $n$ | $\ell_n(r)$ |
|---|---:|---:|
| 2 | $-5\dots6$ | $2.0\times10^{-4}\dots6.0\times10^{-4}$ (all $\ge0$; at the ray's noise floor $\approx g_2\cdot r$) |
| 6 | $-1$ | $\mathbf{0.0201}$ |
| 6 | $0$ | $\mathbf{0.0570}$ |
| 6 | $+1$ | $\mathbf{0.0200}$ |
| 6 | $+2$ | $\mathbf{0.0318}$ |

Every $\ell_n(r)\ge0$, as forced by the projected two-root blocks
(§3).  The decisive fact: **the escape ray carries strictly positive
mode-6 diagonal theta content**, two orders of magnitude above the
solve's noise floor, while the mode-2 content is at noise level.

**Window-cut pairings** (recession pairing
$-\sum_{N'<|n|\le N}q^{n^2}\ell_n(r)$; negative $=$ ray killed):

| family | window | $q=1/4$ | $q=1/2$ | $q=3/4$ | $q=9/10$ |
|---|---|---:|---:|---:|---:|
| 6 | $(0,1]$ | $-0.0100$ | $-0.0201$ | $-0.0301$ | $-0.0361$ |
| 6 | $(0,2]$ | $-0.0102$ | $-0.0221$ | $-0.0402$ | $-0.0570$ |
| 6 | $(1,2]$ | $-0.00012$ | $-0.0020$ | $-0.0101$ | $-0.0209$ |
| 2 | $(0,6]$ | $-0.00017$ | $-0.00055$ | $-0.0012$ | $-0.0016$ |

(Full table: `sdpa_runs/ray_nocuts_deg14_theta_pairings.json`.)

**Verdict (sign rule).**  *Every* tested $q$ kills the ray through the
mode-6 family, with pairings $-0.010\dots-0.057$; the $n=0$ generator
$\hat S_0=T_6(t_1)+\tfrac13$ alone carries the largest share
($\ell_0=0.057$), i.e. even the innermost truncation cut
$(0,N]$ suffices.  The mode-2 family pairs negatively too but at the
ray's noise floor, so no claim is made for it.  The theta ansatz sees
the surviving escape — the working conjecture of PLAN §6 passes its
first falsifiable test.  (Note the $g_4$ paradox is only apparent: the
ray is a *pure $g_4$* escape in the pair sector, but its two-root
diagonal content loads on the mode-6 generators — the $\hat S_n$
squares expand over exactly the triangle/graph_4 labels where the
ray's mass sits.)

## 5. A/B feasibility (C3)

Driver: `theta_ab.py` — captures the exact `--find-ray` problem built
by `sos_search.solve` (solve-call interception; no re-implementation)
and re-solves it with the recession-form atom constraints appended:
(G) rows $L_{n,\alpha}(r)\ge0$ and lumped window cuts
$\sum_{|n|>N'}q^{n^2}L_{n,\alpha}(r)\le0$.  Label order verified
against `--dump-blocks`; $q=\tfrac12$.

**Round 1 — plain diagonal atoms** (leaf $=\hat G_n$ itself;
`sdpa_runs/theta_ab_ray_nocuts_q1_2.json`,
`theta_ab_ray_fingerprint_q1_2.json`):

| variant | added rows | status | squared norm |
|---|---:|---|---:|
| control | 0 | `optimal_inaccurate` | 1.9586 |
| f6 (mode-6 atom) | 6 | **`optimal` — ray survives** | 1.9963 |
| f2 (mode-2 atom) | 18 | `optimal_inaccurate` — survives | 1.9613 |
| both | 24 | `optimal` — survives | 1.9963 |

**The ray does not die.**  The atoms do cut the documented least-norm
ray (its window pairings are $-0.02\dots-0.057$, §4), but the escape
*cone* is fat: the solver returns a nearby recession direction with
every within-degree theta diagonal forced to exact zero at a $+1.9\%$
norm cost.  The surviving ray is *the same escape*: pure $g_4$
($g_4\cdot r=1.0026$, all other $g_\ell=0$, $h_2$-orthogonal), mass
triangle 45.8% / pair 29.1% / graph$_4$ 22.1% — barely rotated from
the control (40.5/27.7/28.6).  Contrast: the three scalar gap cuts
$e_{2,3,4}(I-A_2)\ge0$ make the same problem **infeasible**
(docs/GAP_CUTS_NOTE.md §5).  The diagonal theta tower is a
codimension-few slice of the escape cone, not its dual description.

**Round 2 — localized escalation.**  The plain-atom rows
$\mathcal Q[\hat G_n]\ge0$ are already implied by the projected
two-root blocks; the genuinely new content was only the window upper
cuts.  The escalated atoms use leaves $\hat G_n\cdot m_\alpha$ over
*all* monomial multipliers $m_\alpha(t_1,t_2,s)$ within the degree cap
(per-multiplier atoms $\tau_{f,q,m}=\sum_nq^{n^2}\mathcal Q[\hat
G_nm]$, same L7 majorants since $|m_\alpha|\le1$ on the cube).  These
squares are **not** (E1)-admissible, so they reach the (E1)-dead
directions the projection removed — new content for both the G rows
and the windows.  Row counts after dropping odd-leaf-parity multipliers
(their $\mathcal Q$ vanishes identically on antipodal measures):
138 (family 2) + 6 (family 6) leaves.

**Round 2 results**
(`sdpa_runs/theta_ab_ray_localized_q1_2.json`,
`theta_ab_bound_nocuts_q1_2.json`):

| experiment | status | value |
|---|---|---:|
| ray, f6loc (+8 rows) | `optimal` — survives | norm 1.9963 |
| ray, bothloc (+146 rows) | `user_limit` (CLARABEL) | norm $\ge8.5\times10^9$ at cutoff |
| bound, control | **`unbounded`** (MOSEK, clean) | $-\infty$ |
| bound, f6loc | `unbounded_inaccurate` | $-\infty$ |
| bound, bothloc (+207 rows) | `optimal_inaccurate` (SCS only) | $\approx-1851$ |

Reading: the full localized package inflates the least-norm ray by
$\ge4.5$ orders of magnitude and turns the projected dual from cleanly
unbounded into "bounded at $\approx-1.9\times10^3$" (double-precision
SCS, unverified) — the atoms remove (or nearly remove) exact recession
**without restoring any usable bound**: the value is six orders of
magnitude below the unprojected degree-14 bound $-4.49\times10^{-3}$,
never mind $0$.  Mode-6 alone does nothing beyond zeroing its own
diagonals.

**Where the theta-orthogonal escape lives**
(`theta_ray_blocks.py`, block fingerprint of the surviving `both` ray;
Frobenius norms of $A_b(r)$):

| block | $\|A_b(r)\|_F$ |
|---|---:|
| `hessian_sos` / `hessian_minor` | 11.6 / 9.0 |
| `perpendicular_hessian_sos` / `_minor` | 4.6 / 4.6 |
| `global_perpendicular_tangent_gap` (+minor) | 3.2 / 2.4 |
| `harmonic_4` (the $g_4$ mode) | 1.0 |
| `global_parallel_tangent_gap` (+minor) | 0.52 / 0.33 |
| `two_root_odd_01` / `odd_10` (orientation sectors, rank 3) | 0.21 / 0.21 |
| all even two-root sectors | $\le2\times10^{-4}$ (zeroed by the cuts) |

The escape routes through (a) the KKT-only multiplier blocks
(Hessian and tangent-gap families), (b) the $g_4$ harmonic block, and
(c) the **orientation-odd two-root sectors** — precisely the sectors
the theta families $\hat C_n,\hat S_n$ do not touch ((E1) leaves the
odd sector "nearly free", [Two-root generators](TWO_ROOT_GENERATORS.md)
§8, so the projection keeps it large).  The $q^{n^2}$ diagonal profile
of PLAN §6 was a correct description of the *even-sector* content of
the escape, but that content is removable at $+1.9\%$ norm; the
load-bearing content is elsewhere.

**Round 3 — the ray dies on the all-measures projected cone**
(`sdpa_runs/theta_ab_ray_allm_q1_2.json`; same driver, cone `allm` =
pruned base minus the four KKT toggles, still fully (E1)-projected):

| variant | status | norm |
|---|---|---:|
| control | `optimal` (MOSEK) — ray exists | 0.7474 |
| f6loc | `optimal` — survives | 0.7548 |
| **bothloc (+146 rows)** | **`infeasible` — no recession direction survives** | — |

The proof-carrying (all-measures) projected cone *is* re-bounded by
the localized theta package.  Moreover the recession cone of the
all-measures problem **contains** that of the KKT-inclusive problem
(the KKT cone only adds PSD blocks and equality relations, both of
which shrink the recession cone), so the Round-2 `user_limit` /
"$-1851$" outcomes are explained: the KKT-cone-with-atoms problem is
*also* recession-free, and the double-precision solvers simply failed
to certify it through the near-degenerate geometry.  Corrected
headline: **the full localized theta-atom package kills every
recession direction of the (E1)-projected degree-14 dual on both
cones** — directly certified on the all-measures cone, implied by
inclusion on the KKT cone.  The scalar (non-localized) atoms of PLAN
§6 do *not* suffice; the localization over polynomial multipliers
$m_\alpha$ (reaching the (E1)-dead directions) is what closes the
cone.

## 6. $\varepsilon=0$ assembly (C4)

**Boundedness and margin, all-measures projected cone, degree 14**
(double precision; `sdpa_runs/theta_ab_bound_allm_q1_2.json`,
`theta_ab_ray_allm_gapcuts_q1_2.json`,
`theta_ab_bound_allm_gapcuts_q1_2.json`):

| assembly | ray | $\varepsilon=0$ dual bound |
|---|---|---:|
| projected cone alone | exists (norm 0.747) | `unbounded` (MOSEK, clean) |
| + localized theta atoms (146 rows, $q=\tfrac12$) | **infeasible** | $\approx-5807$ (`optimal_inaccurate`, CLARABEL) |
| + gap scalar cuts alone (3 rows) | **infeasible** (MOSEK, clean) | $\mathbf{-3.18615}$ (`optimal`, MOSEK) |
| + gap cuts + atoms (the full C4 stack) | — | $-3.18615$ (identical to 6 digits) |

Margin report: the assembled $\varepsilon=0$ problem over the
projected all-measures cone + gap cuts + atoms is *feasible and
bounded* with double-precision margin $-3.186$ — but the atoms
contribute nothing on top of the gap cuts ($\Delta\approx2\times
10^{-6}$, solver noise), exactly the slack-at-optimum pattern the gap
cuts themselves show on the unprojected weighted problem (PLAN §4).
Two independent mechanisms now certify that the projected-cone
recession structure is repairable — three scalar spectral cuts, or 146
theta rows — and both repairs leave the same message: **the
finite-degree value of the projected sharp-ansatz cone is
catastrophically weak** ($-3.19$ against $-6.1\times10^{-3}$ for the
unprojected all-measures cone).  The certificate strength at degree 14
lives in the (E1)-*inadmissible* squares; the sharp ansatz plus any
boundedness repair cannot reach $\varepsilon=0$ at this degree.  The
atom majorant constants ($T^{(2)}_{1/2}(N)=238,\,14.4,\,0.73,\,9.3
\times10^{-3},\dots$; $T^{(6)}_{1/2}(N)=1860,\,108,\,5.8,\,7.5\times
10^{-2},\dots$ — exact rationals in the export sidecar) collapse
super-exponentially in $N$, so the atoms' inertness is not a
constant-size artifact: the deep windows are numerically tight and
still inactive.

**$\varepsilon$-trace behavior with atoms** (certificate-side selector
A/B, `theta_trace.py`, degree 14, pruned KKT-inclusive cone, unscaled,
double precision; `sdpa_runs/theta_trace_ab.json`).  The atom window
cuts enter the selector as extra certificate generators
$\lambda_c\,(T_c-\mathrm{row}_c(y))$, $\lambda_c\ge0$:

| $\varepsilon$ | control $\operatorname{tr}C$ | with atoms | atom multipliers |
|---:|---:|---:|---|
| 0.3 | 6.4977 (MOSEK) | **6.4977 — identical to 4 decimals** | 3 active, $\lambda_{\max}=1.60$ |
| 0.1 | 245.50 (CLARABEL, inaccurate) | solver failure | — |
| 0.05 | solver failure | solver failure | — |

The multipliers *activate* but produce zero trace reduction: the atoms
give the optimum alternative degenerate representations, not cheaper
ones.  Below $\varepsilon\approx0.1$ the degree-14 selector is beyond
double precision (documented); the pole-domain
($\varepsilon\le10^{-4}$) re-measurement with atoms is a GMP job — the
export below is the artifact for it.

**GMP-ready export** (`theta_export.py`):
`sdpa_runs/deg14_allm_gapcuts_theta.dat-s` — the $\varepsilon=0$
all-measures **unprojected** pruned cone (harmonics + three/four-point
+ two-root flags + rank relations, no KKT families) + the three gap
scalar cuts + the full localized theta window-cut package at
$q=\tfrac12$, appended as one diagonal block.  Everything new is exact:
generator expansions by rational Chebyshev recurrence, window weights
$q^{n^2}$, majorants $T^{(f)}_q(N)$ (exact fractions in the sidecar
`.map.json`), translated to the exported $z$-coordinates through the
exact `base_point`/`directions` data of the base export — no float
round trip.  Every cut row is verified exactly at the ONB measure
before writing.  On the *unprojected* cone the lower (G) rows are
implied by the full two-root Gram blocks, so only the window upper
cuts are appended.  Why not the projected cone: `--e1-project`
conjugates block matrices in floating point
(`conjugate_label_matrices`), so a projected export cannot currently
be certification-grade without an exact-conjugation pass — recorded
here as the missing piece if a projected GMP measurement is ever
wanted.

## 7. Verdict on the theta route

1. **The working conjecture of PLAN §6 fails as stated.**  The
   diagonal theta atoms $\Theta^{(2)}_q,\Theta^{(6)}_q$ — scalar
   resummations of the admissible generators — do *not* make the
   $\varepsilon=0$ certificate finite at degree 14, and do not even
   kill the projected escape ray: the escape's even-sector
   $q^{n^2}$-diagonal content (§4) is real but removable at $+1.9\%$
   norm.  The load-bearing escape lives in the orientation-odd
   two-root sectors, the $g_4$ harmonic direction, and (on the KKT
   cone) the Hessian/tangent-gap multiplier blocks — none of which the
   $\hat C_n/\hat S_n$ families touch.
2. **The localized theta package is a genuine boundedness repair.**
   Adjoining the per-multiplier atoms $\tau_{f,q,m}$ (leaves $\hat
   G_nm_\alpha$, 146 within-degree rows) kills *every* recession
   direction of the (E1)-projected all-measures dual (`infeasible`,
   clean), and by cone inclusion also of the KKT-inclusive dual.  This
   is the first mechanism *derived from the (E1) solution itself* that
   re-bounds the projected problem — but the gap scalar cuts achieve
   the same with 3 rows, and both leave the projected bound at
   $\approx-3.19$: recession repair, not certificate strength.
3. **At the finite-degree optimum the atoms are inert** — bound
   unchanged to 6 digits on top of the gap cuts, selector trace
   unchanged to 4 decimals at $\varepsilon=0.3$ — the same
   slack-at-optimum pattern as every valid cut measured so far (PLAN
   §4).  Their content is the $\varepsilon\to0$ attainment structure;
   whether they damp the *pole* ($\varepsilon\le10^{-4}$) is exactly
   the question the GMP-ready export is built to answer.
4. **What the theta ansatz misses** (the constructive output of the
   negative result): the escape's odd-sector content.  (E1) leaves the
   orientation-odd sectors nearly free
   ([Two-root generators](TWO_ROOT_GENERATORS.md) §8), the projection
   therefore keeps them large, and the measured escape uses them.  If
   a resummed completion of the sharp ansatz exists, it needs
   *odd-sector modulated families* (orientation leaves
   $\det(x_1,x_2,y)B(t_1,t_2,s)$ with their own $q^{n^2}$ towers) —
   objects (E1) currently constrains only by two slice conditions —
   and/or the kernel family $DQ-\tfrac13(1-s^2)Q(0,0,s)$.  Classifying
   the odd-sector analogue of $\hat C_n$ is the concrete next
   derivation this note motivates.
5. **Practical guidance**: for bound-chasing at finite degree the
   theta atoms are dominated by the gap cuts (3 rows, same effect,
   cleanly solvable).  The atoms' unique asset — exact rational
   control of an *infinite* generator series — matters only where the
   escaping series itself matters, i.e. in the attainment/pole regime.
   Spend GMP time there (the export), not on more double-precision
   bound sweeps.

## 8. The weighted escape (same-day update): $h_2$-localized atoms and sup-caps

[Unprojected escape note](UNPROJECTED_ESCAPE_NOTE.md) (written in
parallel) fingerprints the surviving escape $D$ of the *weighted*
problem (`deg14_h2w_h2all`, $\varepsilon=10^{-3}\to10^{-4}$ selector
certificates, simple pole $10.6\times$/decade): 78% of $|D|$ in
$p_2\times(\cdot)$ product labels, block carriers `h2loc_two_root_*` —
i.e. the **$h_2$-localized theta tower**.  Response, per the
coordinating analysis:

**$h_2$-localized atoms (formalization).**  For $f,q$ as before,
$$\tau^{\mathrm{loc}}_{f,q}(\mu)=\sum_nq^{n^2}\,h_2(\mu)\,
\mathcal Q[\hat G^{(f)}_n](\mu),$$
with label rows $L^{\mathrm{loc}}_n=$ the exact $h_2$-shift of $L_n$:
$h_2\cdot F$ expands as $\tfrac32\,(p_2\times\text{labels})-\tfrac12\,
\text{labels}$ (`theta_atoms.h2_localize`, exact).  Validity: each
term is $h_2\ge0$ times a valid square, so lower cuts hold; and
$h_2=\tfrac{3p_2-1}2\in[0,1]$ on probability measures, so
$h_2\mathcal Q[\hat G_n]\le(a+c_f|n|)^2$ pointwise and **the same
rational tail majorants $T^{(f)}_q(N)$ apply verbatim** to the upper
cuts.  All-measures valid; composes with the reduction lemma.

**Sup-cap cuts.**  The same L7 bound gives, for every single $n$, the
valid rows
$$\mathcal Q[\hat G_n](y)\ \le\ (a+c_f|n|)^2,\qquad
h_2\mathcal Q[\hat G_n](y)\ \le\ (a+c_f|n|)^2$$
(exact rational constants) — the "$N'=n$ slice" of the atom package,
previously subsumed in the majorants, now carried individually.

**Pairing against $D$** (`theta_fingerprint_pair.py`, plain label-space
dot products as directed;
`sdpa_runs/fingerprint_D_e3e4_theta_pairings.json`).  The
per-generator pairings $\ell_n(D)$ are *sign-alternating in $n$* with
magnitudes growing to $10^5$–$10^6$ at the degree edge — the collision
boundary ($s=\pm1$, $T_n(\pm1)=(\pm1)^n$) signature, consistent with
the flat $p_j$ tail of the old ray.  Consequences:

* the one-sided **window** functionals pair with mixed signs
  ($q$-dependent cancellations); kill-sign windows exist (plain
  family 6: $-1.8\times10^6$ at $q=\tfrac34$, low 1; plain family 2:
  $-1.1\times10^3$ at $q=\tfrac12$) but are not uniform in $q$ —
  the "every-$q$ kills" pattern of the projected ray (§4) does not
  recur here;
* the **per-$n$ sup-caps** are the sharp instruments: 13 of 28 pair
  with the kill sign, led by
  $\operatorname{cap}[\text{plain},6,n{=}2]:-5.56\times10^6$,
  $\operatorname{cap}[h_2\text{loc},6,n{=}0]:-5.53\times10^6$,
  $\operatorname{cap}[h_2\text{loc},2,n{=}-3]:-5.3\times10^5$,
  $\operatorname{cap}[\text{plain},2,n{=}-4]:-4.1\times10^5$.  The
  weighted escape *grows* the plain $\hat S_2$-square and the
  $h_2$-localized $\hat S_0$-square (among others) without bound —
  exactly what a pointwise sup bound forbids.

Honest caveat: $D$ is a certificate-side fingerprint
($e(Y_{10^{-4}})-e(Y_{10^{-3}})$), not a feasible moment direction (no
constraint forces $\ell_n(D)\ge0$, and indeed the signs alternate), so
these pairings are the operational matching heuristic of
[Gap cuts note](GAP_CUTS_NOTE.md) practice, not a recession-cone
theorem as in §3.  The decisive test is the GMP selector re-run on the
augmented problem.

**Deliverable**: `sdpa_runs/deg14_h2w_h2all_theta.dat-s` — the
weighted degree-14 problem (`deg14_h2w_h2all` base + its exact
`.map.json`) augmented with, all exactly rational and
ONB-verified: plain windows ($q=\tfrac12$), $h_2$-localized windows,
and the full per-$n$ sup-cap family (plain + localized).  Queue for
GMP: (i) the bound re-solve (expected: unchanged, cuts slack at the
$\varepsilon$-away optimum), (ii) the $\varepsilon$-trace selector at
$10^{-3}/10^{-4}$ against the control values $765.0/8133.6$ — if the
cap generators let the selector shed the $h_2$-localized tower mass,
the growth ratio drops below $10.6\times$/decade; that is the theta
route's remaining live claim.

## Files

| artifact | content |
|---|---|
| `theta_atoms.py` | exact generators, expansions, majorants, pairing CLI, 8 self-tests |
| `theta_ab.py` | ray/bound A/B driver (capture-based, three cones, gap cuts, localized atoms) |
| `theta_ray_blocks.py` | per-block fingerprint of a ray |
| `theta_trace.py` | certificate-side $\varepsilon$-trace selector A/B with atom generators |
| `theta_export.py` | exact GMP-ready exports (windows, localized windows, sup-caps) |
| `theta_fingerprint_pair.py` | pairing of atom functionals against a label-space fingerprint |
| `sdpa_runs/ray_nocuts_deg14.json` | regenerated projected-dual escape ray ($E$ normalization) |
| `sdpa_runs/ray_nocuts_deg14_theta_pairings.json` | full C2 pairing table |
| `sdpa_runs/theta_ab_*.json` | all A/B runs (ray + bound, three cones) |
| `sdpa_runs/theta_trace_ab.json` | $\varepsilon$-trace A/B |
| `sdpa_runs/deg14_allm_gapcuts_theta.dat-s` (+`.map.json`, `.base.dat-s`) | GMP-ready $\varepsilon=0$ artifact (all-measures cone) |
| `sdpa_runs/fingerprint_D_e3e4_theta_pairings.json` | §8 pairing table against the weighted escape $D$ |
| `sdpa_runs/deg14_h2w_h2all_theta.dat-s` (+`.map.json`) | GMP-ready weighted problem + 136 theta cut rows |
