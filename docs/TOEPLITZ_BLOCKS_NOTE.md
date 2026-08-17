# Toeplitz blocks: the circle-pair mechanism as flag-algebra constraints

*2026-08-17, Agent T.  Translates the Toeplitz/double-angle coupling
that closed the circle-pair case
([Cylindrical domination](CYLINDRICAL_DOMINATION.md) §6, item D5) into
valid all-measures flag blocks aimed at the measured weighted escape
([Unprojected escape note](UNPROJECTED_ESCAPE_NOTE.md)).  Code:
`toeplitz_blocks.py` (families, self-tests, pairing CLI),
`toeplitz_export.py` (first-class GMP export), `toeplitz_ab.py`
(double-precision selector A/B).  Artifacts under `sdpa_runs/`.*

**Verdicts.**

1. **T1(a): the per-pair Toeplitz Gram is already spanned.**  Every
   conditional Gram over modulated leaves $(1,\hat C_k,\hat C_{2k})$
   — in the flag reading, two independent leaves — is a principal
   restriction of the existing `two_root_even_00` $\oplus$
   `two_root_even_11` monomial blocks at degree 14
   (machine-checked: each within-degree $\hat C_n,\hat S_n$ lies in
   one parity-sector monomial basis, cap 7; cross-parity entries
   vanish identically for antipodal measures).  The pointwise
   double-angle identity $\cos^2 k\varphi=(1+\cos2k\varphi)/2$ is a
   polynomial identity in $(t_1,t_2,s)$ which the moment reducer
   applies automatically.  Nothing new there — consistent with the
   theta-cut inertness.
2. **T1(b): what is genuinely new is the *unfolding* coupling** —
   the load-bearing inequality of the circle-pair proof is the fiber
   Cauchy–Schwarz $(\mathbb E\cos k\varphi)^2\le\mathbb E\cos^2
   k\varphi$, which compares a **two-leaf** square (graph$_4$ labels)
   with a **same-leaf** average (triangle labels).  Its exact flag
   form is one matrix inequality per two-root sector,
   $$T-G=\mathbb E_{x_1,x_2}\big[\rho\cdot\mathrm{Cov}_y(\varphi)\big]\ \succeq\ 0,$$
   imposed by no family of the current nine-toggle cone (§2.1).  The
   h2-localized copies land on $p_2\times$triangle $+\,p_2\times$
   graph$_4$ — the sectors carrying $58\%+9\%$ of the escape.  The
   literal cross-pair candidate (two disjoint root pairs) either
   factorizes or produces graph$_5$/graph$_6$ labels, where the
   measured escape has $0.0\%$ mass — deprioritized on data (§2.4).
3. **A second new family found in the data: the $h_2$-complement
   (seesaw) blocks** $(1-h_2)\cdot G\succeq0$.  The measured escape
   moves plain and $h_2$-localized copies of the *same* generator in
   opposite directions (pairing ratio $\approx-5$); the complement
   blocks are exactly the valid couplings that forbid this seesaw,
   and their label support lies entirely inside the current problem
   (zero missing labels).
4. **T2: every family sees the escape at matrix level.**
   $M(D)=\sum_LD_L(T_L-G_L)$ has 8–16 negative eigenvalues per
   family, negative mass up to $-9.7\times10^4$ on
   $\|D\|_1=5\times10^4$, stable across both measured decades
   ($10^{-3}\!\to\!10^{-4}$ and $10^{-4}\!\to\!5.3\times10^{-5}$).
   Honest caveat and baselines in §4.
5. **T3: first-class GMP export done** — not a post-hoc file append:
   the 337 dropped dependent directions of the base export make naive
   appending unsound; the families are handed to
   `sos_search.export_sdpa_problem` itself so the exact elimination /
   image-selection machinery redoes the bookkeeping (§5).
   `sdpa_runs/deg14_h2w_h2all_toep.dat-s` + selectors at
   $\varepsilon=10^{-3},10^{-4}$ are READY for the orchestrator's GMP
   queue.
6. **First cut family measured *active* at a selector optimum.**  In
   the double-precision A/B (all solves `optimal`,
   $\varepsilon=10^{-1},3\times10^{-2},10^{-2}$) the augmented
   problem's minimal trace drops by $1.3$–$2.0\%$ and the
   `jensen_even_*` blocks carry trace $0.15$–$0.19$ in the optimal
   certificates — in contrast to the gap cuts and theta atoms, which
   were slack/inert everywhere they were measured.  This range is
   pre-pole; the pole-decade verdict is the queued GMP measurement
   (§6).

## 0. Notation

Roots $x_1,x_2$ and leaf $y$ sampled independently from $\mu$;
$t_1=x_1\!\cdot\!y$, $t_2=x_2\!\cdot\!y$, $s=x_1\!\cdot\!x_2$.  The
azimuth $\varphi$ of the leaf about the root-pair frame enters
polynomially through
$$w:=(t_2-st_1)+i\det(x_1,x_2,y),\qquad
|w|^2=(1-t_1^2)(1-s^2),\qquad
w=\sin\delta\,\sin\theta\,e^{i\varphi}$$
(frame $x_1=e_3$, $x_2$ at angle $\delta$ in the $x_1$–$x_2$ plane,
leaf at polar angle $\theta$): $\operatorname{Re}w^k$ is a polynomial
in $(t_1,t_2,s)$, $\operatorname{Im}w^k=\det\cdot(\text{polynomial})$.
Two-root sectors as in `sos_search.py`: leaf monomials
$t_1^it_2^js^k$ split by the parities of $(i{+}j,\,i{+}k,\,j{+}k)$
into `even_00`, `even_11` (leaf-even), `odd_01`, `odd_10` (leaf-odd,
orientation leaves $\det\cdot B$); `_minor` carries the root weight
$1-s^2$; `h2loc_` copies are multiplied by
$h_2=\tfrac{3p_2-1}2\in[0,1]$.

## 1. T1(a): what the existing degree-14 cone already contains

**Fact 1 (spanning).**  Let $v=(1,A_1,\dots,A_r)$ be any polynomial
leaf functions with $\deg A_i\le7$ whose monomials lie in a single
parity sector, and let
$M_{ab}=\mathbb E\big[(\mathbb E_yv_a)(\mathbb E_{y'}v_b)\big]$ be the
two-root flag Gram.  Then $M=P^{\top}G_{\text{sector}}P$ where $P$ is
the coefficient matrix of the $v_a$ over the sector's monomial basis
and $G_{\text{sector}}$ is the full `two_root_*` block: every dual
multiplier for $M$ lifts to the full block as $P\Lambda P^{\top}$.
The modulated generators fit: $\deg\hat C_n=|n|{+}1$,
$\deg\hat S_n=|n|{+}5$, with $\hat C_{2k},\hat S_{2k}\in$ `even_00`
and $\hat C_{2k+1},\hat S_{2k+1}\in$ `even_11` (machine-checked for
all within-degree $n$, `toeplitz_blocks.py --self-test`).  Mixed
triples such as $(1,\hat C_3,\hat C_6)$ split block-diagonally:
cross-parity entries have odd degree at a root vertex and vanish
identically for antipodal measures.

**Fact 2 (double angle is free).**  $(\operatorname{Re}w^k)^2=
\tfrac12\big((1-t_1^2)^k(1-s^2)^k+\operatorname{Re}w^{2k}\big)$ is a
polynomial identity, applied automatically by
`graph_expectation_label` whenever a same-leaf product is expanded.

So the theta tower, the $(1,\hat C_k,\hat C_{2k})$ Grams, and every
"Toeplitz matrix over modulated leaves with independent leaves" are
inside the current cone.  This sharpens the theta verdict: the
mechanism cannot be reached by any *two-leaf* object at fixed degree.

## 2. T1(b): the genuinely new object

### 2.1 The averaging-contraction (conditional-Jensen) blocks

For a vector $\varphi=(\varphi_\alpha)$ of leaf functions of the
rooted triple $(x_1,x_2,y)$ — monomials $t_1^it_2^js^k$ (even
sectors) or $\det(x_1,x_2,y)\cdot t_1^it_2^js^k$ (odd sectors) — and
a nonnegative root weight $\rho\in\{1,\,1-s^2\}$ define

$$T_{ab}=\mathbb E\big[\rho\,\varphi_a(y)\varphi_b(y)\big]
\quad(\text{same leaf; triangle labels}),\qquad
G_{ab}=\mathbb E\big[\rho\,\varphi_a(y)\varphi_b(y')\big]
\quad(\text{independent leaves; the existing two-root Gram}).$$

**Theorem (validity, all measures).**
$T-G=\mathbb E_{x_1,x_2}\big[\rho\cdot\mathrm{Cov}_y(\varphi)\big]
\succeq0$ for every probability measure $\mu$ on $S^2$; the
$h_2$-localized copy $h_2\cdot(T-G)\succeq0$ and the complements
$(1-h_2)\cdot G\succeq0$, $(1-h_2)(T-G)\succeq0$ are likewise valid
because $0\le h_2\le1$ for every measure
($p_2=\operatorname{tr}\Sigma^2\in[1/3,1]$).  No KKT content; the
constraints compose with the reduction lemma.
*Proof.*  $\mathrm{Cov}_y(\varphi)\succeq0$ pointwise in the roots;
$\rho\ge0$; take expectations.  The label expansions are exact:
$T$ entries reduce through 3-point monomials $(i,j,k)\mapsto$
`expectation_label`, $G$ entries through the standard 4-point
two-root reduction (no $y$–$y'$ edge), $h_2$-shifts through
$\tfrac32\,p_2\times L-\tfrac12L$.  $\square$

Equivalently, with $f=(\varphi(y),\,\mathbb E_y\varphi)$ the mixed
instant/averaged flag vector on the 3-root type,
$U=\mathbb E[ff^{\top}]=\begin{pmatrix}T&G\\G&G\end{pmatrix}\succeq0$,
and since $G(I-GG^{+})=0$ automatically, $U\succeq0\iff G\succeq0$
and $T-G\succeq0$: **one new matrix inequality per sector**, the
"unfolding" $[\![A]\!]^2\le[\![A^2]\!]$ of flag calculus.  No family
of the current nine-toggle cone imposes it: the two-root Grams
constrain $G$ alone, the three-point `flag_k` Grams constrain
*different* triangle-sector matrices (one root, two averaged leaves),
`--no-pointwise-sos` removes the same-leaf moment matrix entirely,
and no existing block mixes same-leaf with averaged-leaf entries.

### 2.2 Why this *is* the Toeplitz mechanism

D5's measure-level matrix
$[[\sigma,\mu_k,\mu_{2k}],[\ast,\sigma,\mu_k],[\ast,\ast,\sigma]]
\succeq0$ tested against latitude functions $(g_0,g_1,g_2)$ is
$\int|g_0+g_1e^{ik\varphi}+g_2e^{2ik\varphi}|^2d\mu\ge0$: squares of
functions *pointwise on one sample*, referencing the frame.  With the
frame sampled as a root pair, those are exactly the same-leaf squares
building $T$; the mode averages $\zeta_k$ are the averaged flags
building $G$; and the proof of D5,
$1+\operatorname{Re}\zeta_{2k}=2\mathbb E\cos^2k\varphi\ge
2(\mathbb E\cos k\varphi)^2$, is the centered vector
$\varphi=\cos k\varphi$-rep $-$ its conditional mean — i.e. a
diagonal entry of $T-G$ after the automatic double-angle reduction of
Fact 2.  Sharpness transfers: at a coherent one-orbit measure
$\{\pm u\}$ (the flag image of $|\zeta_k|=1$) every even-sector
$T-G$ vanishes identically — machine-checked exactly, alongside
strict positivity of the $C_1$-direction gap at generic circle
measures (`--self-test`).  Odd fiber modes are not lost: for
antipodal measures $\mu_{2k+1}$ is odd across latitudes, and
mode-$(2k{+}1)$ content enters the leaf-even sectors through
latitude-odd radials $t_1^{\text{odd}}\operatorname{Re}w^{2k+1}$ and
the odd sectors through $\det$-pairing — all inside the monomial
bases.

### 2.3 The unrooted (pair-sample) copy

Taking the whole pair $(X,Y)$ as the sample, flags $(X\!\cdot\!Y)^d$,
$d\in\{0,2,4,6\}$:
$$T_{ab}=p_{a+b}\ (\text{pair labels}),\qquad
G_{ab}=p_a\,p_b\ (\text{product labels}),\qquad T-G\succeq0,$$
the covariance matrix of $((X\!\cdot\!Y)^d)_d$ — a Hankel-vs-product
coupling between the pair sector and the pair-product sector that the
`empty_type_flag` Gram alone does not imply.  Its $h_2$-localized
copy loads on $p_2\times$pair-products ($10.8\%$ of the escape) and
pairs against the explosively-activating `h2loc_empty_type_flag`.

### 2.4 The cross-pair candidate, settled by the data

The mission's candidate (i) — two disjoint root pairs — was
evaluated first.  A product of two pair-flags factorizes
($\mathbb E[F(\text{pair}_1)]\,\mathbb E[G(\text{pair}_2)]$: rank
one, no content).  Non-factorizing versions must either share the
leaf ($\mathbb E_y[\Phi_a\Phi_b]$, $\Phi_a(y)=\mathbb E_{\text{pair}}
[W_a]$: graph$_5$ entries) or contract a harmonic index between the
pairs (`spin2_flag` pattern: graph$_6$ entries).  The measured escape
has $0.0\%$ mass on graph$_5$/graph$_6$ (both decades).  The
"$4$-to-$6$-sample" sector where the escape actually lives is
graph$_4$ + $p_2\times$graph$_4$ — which is exactly the label support
of $T-G$ (cross corner) and of its $h_2$-localization ($p_2\times$
labels are the 6-sample disconnected ones).  So candidate (ii)
absorbs candidate (i)'s intent; the literal 4-root blocks are
deprioritized on data, not on principle.

### 2.5 The seesaw and the $h_2$-complement family

Measured per-generator Jensen pairings $\langle t_n-g_n,D\rangle$
(1×1 diagonal of $T-G$ on $\hat C_n/\hat S_n$; `--pair-d`) alternate
in sign with $|n|$ and **anti-correlate between plain and h2loc
copies of the same generator**: the ratio
h2loc/plain lies in $[-6.3,-4.5]$ (mean $\approx-4.9$) for *all* 14
generators with pairing above $10^3$ — e.g. $\hat C_{-4}$:
$-1.44\times10^6$ plain vs $+7.04\times10^6$ h2loc; $\hat S_0$:
$+9.7\times10^5$ plain vs $-4.8\times10^6$ h2loc.  Since the h2loc
row is $\tfrac32(p_2\times L)-\tfrac12L$, a uniform ratio $\approx-4.9$
says the escape's $p_2\times L$ content is $\approx-2.9\times$ its
plain-$L$ content, coherently across the whole modulated tower
(maximally non-measure-like: a genuine measure near ONB has ratio
$+p_2\approx+\tfrac13$).  The escape shuttles two-root mass between
the plain and $h_2$-localized towers — a direction nothing in the
problem constrains, because the plain and h2loc blocks are separate
PSD constraints.  The valid
couplings that forbid it are the complements
$$G_{\text{plain}}-G_{\text{h2loc}}=(1-h_2)\,G\succeq0
\qquad(\text{and }(1-h_2)(T-G)\succeq0),$$
whose label support is exactly the union of the existing plain and
h2loc block labels — **zero labels outside the current problem**.

## 3. Validity ledger

| item | statement | status |
|---|---|---|
| J1 | $T-G=\mathbb E[\rho\,\mathrm{Cov}_y(\varphi)]\succeq0$, all sectors, $\rho\in\{1,1-s^2\}$ | theorem (§2.1); exact label expansions machine-checked against direct covariance on random atomic antipodal measures, rel. gap $<10^{-14}$, and exactly (rational equality) at ONB and the 4-point cross |
| J2 | pair-sample copy $[p_{a+b}]-[p_ap_b]\succeq0$ | same theorem, sample = pair |
| J3 | $h_2\ge0$, $h_2\le1$ for every measure | $p_2=\operatorname{tr}\Sigma^2\in[1/3,1]$ |
| J4 | $h_2$-localized and $(1-h_2)$-complement copies PSD | J1–J3; scalar $\ge0$ times PSD |
| J5 | exact PSD-ness at ONB and cross measures | rational pivoted-elimination PSD check, every exported family |
| J6 | $h_2$-shift convention matches `--h2-localized-all` | exact match (gap $0.0$) against `blocks_deg14_h2w_h2all.json`, plain and h2loc, incl. minors |
| J7 | spanning Fact 1 / degree-parity of $\hat C_n,\hat S_n$ | machine-checked (`--self-test`) |
| J8 | coherence sharpness: even-sector $T-G\equiv0$ at $\{\pm u\}$ | machine-checked exactly |

KKT-only assumptions: none.  All families are all-measures valid and
compose with the reduction lemma.  (They are **not** members of the
sharp (E1) ansatz: $\mathrm{Cov}_y\ne0$ at the zero measures, so
these are $\varepsilon$-regime/pole instruments, same status as the
unprojected cone, not sharp-certificate generators.)

## 4. T2: pairing against the measured escape

$D=e(Y_{10^{-4}})-e(Y_{10^{-3}})$ (`fingerprint_D_e3e4.json`, 1283
labels, $\|D\|_1=5.0\times10^4$); test
$M(D)=\sum_LD_L(T_L-G_L)$ per family (mission criterion: a valid PSD
family with $M(D)\not\succeq0$ can cut the escape).  Full tables:
`sdpa_runs/toeplitz_T2_e3e4.json`, `toeplitz_T2_e4e5.json`.

| family (size) | min eig | #neg | neg sum | labels outside $D$ |
|---|---:|---:|---:|---:|
| jensen_pair (4) | $-2.2\times10^2$ | 2 | $-3.3\times10^2$ | 0 |
| jensen_even_00 (30) | $-2.4\times10^3$ | 14 | $-9.0\times10^3$ | 4 |
| jensen_even_11 (30) | $-1.8\times10^3$ | 12 | $-7.4\times10^3$ | 4 |
| h2loc_jensen_even_00 (30) | $-9.3\times10^3$ | 12 | $-3.5\times10^4$ | 8 |
| h2loc_jensen_even_11 (30) | $-1.1\times10^4$ | 14 | $-4.0\times10^4$ | 8 |
| h2loc_jensen_pair (4) | $-3.6\times10^3$ | 1 | $-3.6\times10^3$ | 0 |
| **h2comp_gram_even_00 (30)** | $-9.2\times10^3$ | 14 | $-4.0\times10^4$ | **0** |
| **h2comp_gram_even_11 (30)** | $-1.0\times10^4$ | 14 | $-4.0\times10^4$ | **0** |
| **h2comp_gram_odd_01/10 (30)** | $-1.4\times10^4$ | 15 | $-7.0\times10^4$ | **0** |
| h2comp_gram_pair (4) | $-4.7\times10^3$ | 1 | $-4.7\times10^3$ | 0 |
| jensen/h2loc odd, minors, h2comp_cov | $-2.1\times10^3\ldots-2.8\times10^4$ | 7–16 | to $-9.7\times10^4$ | 4–54 |

Same pattern, larger magnitudes, on the second decade
($10^{-4}\!\to\!5.3\times10^{-5}$): h2comp_gram families reach neg
sum $-1.7\times10^5$; the fingerprint is stable.

**Honest caveats.**  (i) $D$ is a certificate-side fingerprint, not a
moment direction; nothing forces $M(D)\succeq0$ even for *existing*
families, and indeed the corner baselines $T(D)$, $G(D)$ alone are
similarly indefinite.  The eigen test is the necessary-condition
filter of the mission, not a trace-law theorem.  (ii) The theta
sup-caps also had kill-sign pairings and were inert on the GMP
selector trace (765.0/8133.57 unchanged) — the lesson being that
codimension-few *scalar* slices reroute.  What is different here:
these are full PSD **matrix** families (8–16 independent negative
directions each), they tie the growing graph$_4$/$p_2\times$graph$_4$
sectors to the nearly-static triangle/pair sectors (Jensen) and the
plain tower to the h2loc tower (complements), and the h2comp_gram
family lives entirely on labels the problem already couples.  Whether
that converts into trace reduction is exactly what the T3 artifacts
measure.

Per-generator Jensen rows (the interpretable 1×1 slices; full list in
the JSONs): sign-alternating in $n$, anti-correlated plain vs h2loc
(§2.5), magnitudes $10^2$–$7\times10^6$.

## 5. T3: artifacts (READY for the GMP queue)

**Why a first-class export.**  The base export drops 337 directions
whose images on the base blocks are linearly dependent
(`dropped_dependent_directions: 337`); a block appended to the
finished `.dat-s` would silently pin those coordinates — an *invalid
strengthening*.  `toeplitz_export.py` therefore patches
`sos_search.export_sdpa_problem` and passes the families into the
original exact pipeline, which redoes the equality elimination and
image selection with the new blocks included (m grows as previously
collapsed directions become distinguishable).  Family entries are
integers/halves; the float hand-off is verified round-trip-exact
entry by entry, and every family is re-verified (exact equality +
exact PSD at ONB and cross) inside the run.

Appended families (11 blocks): `jensen_pair`, `h2loc_jensen_pair`,
`h2comp_gram_pair` (4×4 each); `jensen_even_00/11`,
`h2loc_jensen_even_00/11` (cap 7, trimmed 30→25/24 to the run's label
set); `h2comp_gram_even_00/11`, `h2comp_gram_odd_01/10` (30×30,
untrimmed).

**Export outcome** (2026-08-17, `deg14_h2w_h2all_toep.export.log`):
m = **794** (base 790: only 4 of the 337 dropped directions become
distinguishable through the new blocks — the GMP cost stays at the
base level), 85 blocks, 1,002,498 entries, 66.8 MB,
`objective_shift_exact = 2/3`, `equality_rank` 184 unchanged.  The
near-coincidence 333 vs 337 also quantifies the soundness point:
a naive file-append would have silently pinned exactly those 4
coordinates.

| artifact | content |
|---|---|
| `sdpa_runs/deg14_h2w_h2all_toep.dat-s` (+`.map.json`, `.families.json`, `.export.log`) | the augmented problem, m=794, 85 blocks; **bound = objValPrimal + 2/3** |
| `sdpa_runs/sel_toep_1em3.dat-s` | selector, m=795, 86 blocks (+slack); `--bound=-6.676666666666666666666666666666666666667E-1` |
| `sdpa_runs/sel_toep_1em4.dat-s` | selector, m=795, 86 blocks (+slack); `--bound=-6.667666666666666666666666666666666666667E-1` |
| `sdpa_runs/blocks_deg14_h2w_h2all_toep.json` | 85-block dump (base dump + the 11 families) for `fingerprint_blocks.py` / `fingerprint_expand.py` on `sel_toep_*` results |

GMP queue (orchestrator; one at a time, 128-bit for selectors):
compare `sel_toep_*` traces against the controls **765.0 / 8133.57**;
a growth ratio below $10.6\times$/decade is the success signal.
(Adding certificate generators can only lower the minimal trace; the
question is whether the drop is material and whether the *ratio*
moves.)

**Double-precision A/B** (`toeplitz_ab.py`): rebuilds the selector in
label space from the captured model and compares control vs +Jensen
with MOSEK.  At $\varepsilon\le10^{-3}$ double precision fails
(documented wall); on the tractable ladder
$\varepsilon\in\{10^{-1},3\times10^{-2},10^{-2}\}$ the blocks are
*active* and reduce the trace at every point (§6).

## 6. Measurements on this date

- Self-tests: **6/6** (`toeplitz_blocks.py --self-test`): exact
  expansions $=$ direct covariance (rel.\ gap $<10^{-14}$), PSD at
  random atomic measures, coherence sharpness at one-orbit measures,
  strict $C_1$ gap at generic circle measures, $h_2\ge0$, generator
  spans.  Dump cross-check (J6): **exact, gap 0.0**, eight blocks
  including h2loc and minors.
- T2 eigen tables (§4): both decades;
  `sdpa_runs/toeplitz_T2_e3e4.json`, `toeplitz_T2_e4e5.json`.
- Export: succeeded, m = 794, 11 families all verified in-run (exact
  PSD at ONB + cross, exact float round-trip); selectors and dump
  built (§5).
- Double-precision selector A/B (`toeplitz_ab.py`, MOSEK, all six
  solves `optimal`; label-space form of the selector, KKT-inclusive
  cone):

  | $\varepsilon$ | control $\operatorname{tr}$ | +Jensen $\operatorname{tr}$ | $\Delta$ |
  |---:|---:|---:|---:|
  | $10^{-1}$ | 5.6038 | **5.4920** | $-2.0\%$ |
  | $3\times10^{-2}$ | 7.8830 | **7.7691** | $-1.4\%$ |
  | $10^{-2}$ | 12.6538 | **12.4915** | $-1.3\%$ |

  **The Jensen blocks are *active* at the optimum** —
  `jensen_even_11` carries trace $0.15$–$0.19$ and `jensen_even_00`
  up to $0.16$ in the optimal certificates (top-10 blocks) — in sharp
  contrast to every previously measured valid cut (gap cuts, theta
  windows/caps: multipliers at zero effect, traces unchanged to 4+
  digits).  Honest limits: the reduction is small and this
  $\varepsilon$-range is pre-pole (growth only $2.26\times$/decade
  here vs $10.6\times$ in the pole decade); at $\varepsilon=10^{-3}$
  MOSEK fails on control and augmented problem alike (the documented
  double-precision wall), so whether the *pole coefficient* drops is
  decided by the GMP `sel_toep_*` solves.  Data:
  `sdpa_runs/toeplitz_ab_selector.json` (per-block traces included).

## Files

| file | content |
|---|---|
| `toeplitz_blocks.py` | exact families ($T$, $G$, $T-G$, complements), self-tests, `--pair-d` eigen tables, `--coverage`, `--check-dump` |
| `toeplitz_export.py` | first-class export driver + exact verification + selector builder |
| `toeplitz_ab.py` | label-space double-precision selector A/B |
| `toeplitz_dump_blocks.py` | augmented 85-block dump writer for the fingerprint tools |
| `sdpa_runs/toeplitz_T2_e3e4.json`, `toeplitz_T2_e4e5.json` | T2 eigenvalue tables and generator rows |
| `sdpa_runs/deg14_h2w_h2all_toep.dat-s` (+sidecars), `sel_toep_1em3/1em4.dat-s` | T3 GMP-ready artifacts |
| `sdpa_runs/toeplitz_families_exact.json` | exact rational label matrices of the 11 exported families (3703 label matrices; reproducibility) |
