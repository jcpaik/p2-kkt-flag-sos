# The equal-potential defect square in the all-measures cone

*Agent D2, 2026-08-20.  Mission: port the KKT tower's dominant object —
the equal-potential defect square
$h_2\,[K(X\!\cdot\!Y)-K(Z\!\cdot\!Y)]^2$ of
docs/LIMIT_EXTRACTION_NOTE.md §4a — into the proof-carrying
(all-measures) cone as valid blocks, gate every candidate through the
escape pre-test, and stage the v4 exports.  Code:
`toeplitz_blocks.py` (`build_defect_family`, `build_krow_family`,
`build_leafpoly_family`, `s4`/`s6` root weights, self-tests),
`toeplitz_export.py` (`--version v4`).  Data:
`sdpa_runs/fingerprint_D_am16_toep3_deep.json` (new),
`sdpa_runs/defect_pretest_deg{14,16}.json`.*

Notation.  $p(t)=8t^6-12t^4+5t^2=\tfrac14\big(K(t)+\tfrac43\big)$
(coefficient vector $(5,-12,8)$ on $t^2,t^4,t^6$),
$\Delta p = p(t_1)-p(t_2)$ with $t_1=X\!\cdot\!Y$, $t_2=Z\!\cdot\!Y$,
$s=X\!\cdot\!Z$; $U_p(x)=\int p(x\cdot y)\,d\mu(y)$, and the **defect**
$D(X,Z)=U_p(X)-U_p(Z)=\mathbb E_Y[\Delta p\mid X,Z]$.  Since constants
cancel in differences, $K(X\!\cdot\!Y)-K(Z\!\cdot\!Y)=4\,\Delta p$.

## 1. The objects and their validity (D2-1)

### 1.1 (a) The plain CE-square is the variance form, and both are already spanned

The two-root CE-square with leaf $K(t_1)-K(t_2)$ is, by definition of
the two-leaf Gram corner ($Y,Y'$ independent leaves over the root pair),

$$Q[\Delta K]\;=\;\mathbb E\big[(K(X\!\cdot\!Y)-K(Z\!\cdot\!Y))
(K(X\!\cdot\!Y')-K(Z\!\cdot\!Y'))\big]
\;=\;\mathbb E_{X,Z}\big[(U(X)-U(Z))^2\big]
\;=\;2\,\mathrm{Var}_\mu(U),$$

the last step because $X,Z$ are independent.  **So candidate (a) and
the scalar variance row (b) are the same object up to the factor 2**:
$e(\mathrm{Var}) = \tfrac12\,v_K^{\mathsf T} G\, v_K$ with $v_K$ the
coefficient vector of $\Delta K$ in the `even_00` leaf basis.  The
label identity (machine-checked in `toeplitz_blocks.py --self-test`,
"defect Gram (0,0) == 2 Var(U_p)", gap $2\times10^{-13}$ float /
exact-rational in `verify_family`) is

$$\mathrm{Var}(U_p)=\sum_{i,j}k_ik_j\Big[
\big(\text{'triangle'},0,2i,2j\big)
-\big(\text{'pair'},2i\big)\big(\text{'pair'},2j\big)\Big],
\qquad (k_1,k_2,k_3)=(5,-12,8),$$

triangle labels from $\mathbb E[U_p^2]$ (the $Y$–$Z$ edge carries
exponent 0), disconnected products from $(\mathbb E U_p)^2$.

**Spanning check**: $\Delta K$ has leaf degree 6, all monomials in the
`even_00` parity sector, so it lies in the base `two_root_even_00`
Gram block for every cap $\ge6$ — i.e. at degrees 14/16/18 in both
cones (caps 7/8/9).  The $h_2$-localized copy sits in
`h2loc_two_root_even_00` (`--h2-localized-all`), the
$(1-h_2)$-localized copy in `h2comp_gram_even_00` (v1 family, kept
untrimmed 45/45 at deg-16), the $(1-s^2)$-weighted copy in the base
minor block.  **Conclusion (a): already inside the exported two-root
blocks at every degree.**  The reason the am tower does not exploit it
the way the KKT tower does is therefore exactly the missing KKT
relation: in the KKT-inclusive cone the gradient rows make
$h_2\mathcal Q[\Delta K]$ an *exact null* (LIMIT_EXTRACTION §4a,
pairing $3.3\times10^{-20}$), so the dual adds it with a free
multiplier $\sim10^{10}$; in the am cone it is not null (pairing
$4.3\times10^{-3}$) and no valid *inequality* can substitute for an
*identity* the feasible span does not know.

### 1.2 (b) The variance family and the section lemma

All the "richer" variance objects are covariances and therefore valid
for every measure:

* **defect Gram** ($G$-corner on leaves $s^m\Delta p$):
  $G_{mm'}=\mathbb E[s^{m+m'}D^2]=:m_{(m+m')/2}$ (odd $s$-powers vanish
  by antipodality) — the Hankel moment matrix in $s^2$ of the
  nonnegative measure $D(X,Z)^2\,d\mu(X)d\mu(Z)$ pushed to $s$; PSD as
  a second-moment matrix of averaged flags.  $m_0=2\mathrm{Var}(U_p)$.
* **defect covariance** ($A=T-G$ on the same leaves):
  $\mathbb E_{X,Z}[\mathrm{Cov}_Y(s^m\Delta p)]\succeq0$ — the
  root-conditional covariance (conditional Jensen), same proof as the
  v1 `jensen_*` families.
* **$h_2$- and $(1-h_2)$-localized copies**: nonnegative scalar
  ($0\le h_2\le1$) times a PSD family.
* **$(1-s^2)$- and $s^{2m}$-weighted copies**: pointwise-nonnegative
  root weights.
* **kernel-row (unconjugated) versions** `krow_*`: leaves
  $s^m(t_1^{2i}-t_2^{2i})$, $i=1,2,3$ — contain the conjugated
  direction $(5,-12,8)$ per modulation block.

All verified in-run: exact label expansion $=$ exact-rational direct
integration, exact PSD at the ONB measure and the 4-point cross
(`toeplitz_export.verify_family`), plus float agreement; and at the
ONB measure (equal-potential face) every defect **Gram** entry
vanishes identically — the family is precisely the all-measures shadow
of the KKT equal-potential identity (self-test "defect Gram vanishes
at the equal-potential ONB measure", exactly 0).

**Section lemma (the structural fact of this mission).**  For any
root-measurable weight $w(s)\ge0$ and leaf vector $\phi$,

$$\mathbb E[w(s)\,\mathrm{Cov}_Y(\phi)] \;=\;
\mathbb E[\mathrm{Cov}_Y(\sqrt{w}\,\phi)]
\quad\text{(bilinearity of $\mathrm{Cov}_Y$ over root factors)},$$

so every weighted defect/variance object with polynomial
$\sqrt{w}=s^m$ is an **exact congruence section** $P^{\mathsf T}MP$ of
the plain family at leaf cap increased by $m$; and every congruence
section of a PSD block already in the export is implied by it — the
dual multiplier $P\Lambda P^{\mathsf T}$ is available inside the
parent.  Consequence, checked row-by-row against the deg-16 am export:

| object | parent in v3 export | status at deg-16 |
|---|---|---|
| $Q[\Delta K]$, $\mathrm{Var}(U)$ row | base `two_root_even_00` | **implied** (leaf deg 6 $\le$ 8) |
| $h_2\cdot$, $(1-h_2)\cdot$, $(1-s^2)\cdot$ copies | `h2loc_two_root_*`, `h2comp_gram_*`, base minor | **implied** |
| $m{=}1$ copy $s\Delta p$ (cov) | `jensen_even_11_d8` (rows kept) | **implied** |
| $m{=}2$ copy $s^2\Delta p$ (gram) | base `two_root_even_00` (cap 8) | **implied** |
| $m{=}2$ copy $s^2\Delta p$ (cov) | `jensen_even_00_d8` rows **(6,0,2),(0,6,2) trimmed** | **out of vocabulary** (its own entries too — trimmed away) |
| $m{=}3$ copy $s^3\Delta p$ (gram) | leaf degree 9 $>$ base cap 8 — **no parent exists** | **genuinely new** (in-vocabulary; the shifted Hankel row) |
| $(1-s^2)s^{2m}\Delta p$ (gram), $m\le2$ | base minor block cap 7; $m{=}2$ leaf degree 8 $>7$ | $m\le1$ implied; **$m{=}2$ genuinely new** |
| $s^4$/$s^6$-weighted Jensen, caps 6/5 | see the sharpened accounting below | mostly sections; **a few genuinely new rows** |

So **within the two-root polynomial algebra, every in-vocabulary
defect/variance object of leaf degree $\le$ the base cap is exactly
redundant** — the am cone already contains the equal-potential shadow
as inequalities; what it lacks and cannot be given at finite degree is
the *identity*.  The genuinely new finite-degree content is the
**trim/cap-boundary re-import**, and the section lemma cuts both ways:
a weighted row whose $s^k$-image row was dropped from the parent *for
out-of-vocabulary entries* has the **same** entries and is dropped
from the weighted family too.  The exact inventory of kept rows that
are *not* congruence sections of any exported family
(machine-computed against the kept bases; deg-16):

* `jensen_even_00_s6_d5` (+ `h2loc` copy): row $(3,1,1)$ — the
  $s^3$-shifted leaf is $(3,1,4)$, dropped from `jensen_even_11_d8`,
  but *this* family's restricted pairing keeps its entries
  in-vocabulary;
* `jensen_even_11_s6_d5` (+ `h2loc`): row $(4,0,1)$ — image $(4,0,4)$;
* `h2comp_cov_*_s4_d6` / `_s6_d5`: all rows whose image has leaf
  degree 8 (6+3 / 3+6 rows per sector) — the $(1-h_2)$-complement
  parents exist only at cap 7, so the whole shell-8 complement
  coupling is new;
* `defect_gram_d9` / `h2loc_` / `comp_` copies: the $m{=}3$ row
  $s^3\Delta p$ (leaf degree 9, beyond every cap) — extending the
  defect Hankel to $[m_1\,m_2;m_2\,m_3]\succeq0$;
* `defect_gram_minor_d8`: the $m{=}2$ row (beyond the minor cap 7).

At deg-14 the analogues: `jensen_even_00_s6_d4` row $(2,0,2)$,
`jensen_even_11_s4_d5` row $(2,0,3)$ (both hitting leaf $(2,0,5)$,
dropped from `jensen_even_11_d7`), the shell-7 `h2comp_cov` rows, and
the defect-Gram $m{=}2$ row (leaf degree 8 $>$ cap 7 — at deg-14 even
the *plain* modulated Gram tower outruns the base block).  These
targets are exactly the measured escape carriers: the deg-18 am
escape's top swap-antisymmetric directions are $(2,4,2)-(4,2,2)$ and
$(5,1,1)-(1,5,1)$, and its $X\!\cdot\!Z$-modulation index climbs
$k\approx1\to2\to3$ (LIMIT_EXTRACTION §5) — $k=2,3$ are the
$s^4$/$s^6$ weights and the $m=2,3$ defect modulations.

### 1.3 (c) The modulated copies $(X\!\cdot\!Z)^{2m}$, $m=1,2$

Per-degree availability of the conjugated tower $s^m\Delta p$ (entry
degree $2(6+m)\le d$ for the Gram corner):

| leaf | deg-14 | deg-16 | deg-18 |
|---|---|---|---|
| $\Delta p$ (cov / gram) | in / in | in / in | in / in |
| $s\,\Delta p$ | **out** / in | in / in | in / in |
| $s^2\Delta p$ | out / in* | **out** / in | in / in |
| $s^3\Delta p$ | out / out | out / **in*** | (expected in) |

("in" = exported after trimming; "out" = at least one entry label
outside the vocabulary; * = **beyond every existing block's cap, i.e.
genuinely new there** — the deg-14 vocabulary happens to contain the
degree-16 star labels `('graph_4',0,0,4,0,6,6)`-type that the $m=2$
Gram entries need, and the deg-16 vocabulary the corresponding $m=3$
stars.)  The trimmed families are what v4 ships; the $m$-tail that
never fits at any fixed degree is the object of the resummation
program (§4).

## 2. The pre-test (D2-2)

**New escape fingerprint.**  The am-tower deep growth direction at
degree 16 was computed from the proof-carrying selector ladder,
$D = e(Y_{5\mathrm{e}{-}6}) - e(Y_{1\mathrm{e}{-}5})$
(`sel16_am_toep3_1em5/5em6.result`, both pdFEAS; label matrices from
the exporter capture, all 110 blocks) —
`sdpa_runs/fingerprint_D_am16_toep3_deep.json`, 1944 labels.  Its top
components are $p_2\times$`graph_4` products in swap-paired
sign-opposite pairs, confirming the boundary-layer picture.

**M(D) eigenvalue tables** (`family_pairing`, T2 protocol: negative
eigenvalue = the escape violates the candidate = cut).  Full tables in
`sdpa_runs/defect_pretest_deg{14,16}.json`; deg-16 summary against
the am16 deep escape (calibration rows [P] are families already in the
export):

| family (trimmed) | size | min eig | max eig | verdict |
|---|---:|---:|---:|---|
| `defect_cov` | 2 | $-2.22\times10^3$ | $-4.15\times10^2$ | cut ($M(D)\prec0$ entirely) |
| `h2loc_defect_cov` | 2 | $+1.99\times10^3$ | $+1.52\times10^4$ | **weight-signed — excluded** |
| `defect_gram` | 3 | $-1.01\times10^3$ | $+9.29\times10^2$ | cut |
| `h2loc_defect_gram` | 3 | $-8.68\times10^3$ | $+4.40\times10^3$ | cut |
| `defect_comp_cov` ($(1-h_2)$) | 2 | $-1.75\times10^4$ | $-2.41\times10^3$ | **strongest cut** |
| `defect_comp_gram` | 3 | $-5.38\times10^3$ | $+9.58\times10^3$ | cut |
| `defect_cov_minor`, `defect_gram_minor` | 1, 2 | $-1.8\times10^3$, $-1.6\times10^2$ | | cut |
| `krow_cov` / `krow_gram` (+h2loc) | 8–9 | $-9.7\ldots-88$ | | cut (mild; exactly redundant — not exported) |
| `jensen_even_00_s4_d6` (+11, h2loc, comp) | 17–20 | $-9.9\ldots-62$ | | cut, **new content** |
| `jensen_*_s6_d5` (+h2loc, comp) | 14 | $-3.2\ldots-24$ | | cut, **new content** |
| [P] `jensen_even_00_d8` | 38 | $-29.1$ | $+32.9$ | (present) |
| [P] `h2loc_jensen_even_00_d8` | 38 | $-164$ | $+145$ | (present) |
| [P] `jensen_*_s2_d7` | 26–27 | $-21$, $-20$ | | (present) |

Scale caveat: the defect families' large raw eigenvalues carry the
integer conjugation norm $|v|^2=466$; per unit leaf norm their
violation ($\approx-4.7$) is *below* the present parents' min
eigenvalues — consistent with §1.2's redundancy finding.  What is
qualitatively new in the table is (i) the **whole-matrix negativity**
of `defect_cov` / `defect_comp_cov` (every direction of the family
decreases along the escape — the escape is assembling the defect
square itself, as §4a predicted), and (ii) the s4/s6 families cutting
at the same level as the present families while carrying rows the
export does not have.  Against the deg-14 KKT escapes
(`fingerprint_D_toep2/toep3_deep`) all of the above are also
kill-signed (except `h2loc_defect_cov` at v2), 30–100× stronger — the
KKT dual's appetite for the defect direction.

**Scalar sign-rule rows (§6.3 deviatoric test)** against the am16 deep
escape:

| functional | $\langle\cdot,D\rangle$ | verdict |
|---|---:|---|
| $w_{\mathrm{dev}}=\iint(t^2-\tfrac13)^2 = p_4-\tfrac23p_2+\tfrac19$ | $+0.985$ | **weight** (positive: can absorb the pole; matches §4b's surviving deviatoric leaf) |
| $h_2$ | $-0.566$ | (escape drives $h_2\downarrow0$: wall adversary $p_2\to1/3^+$) |
| $e_2(I-A_2)=6+6p_2-8p_4$ | $-7.95$ | cut-signed; not exported (historically inert at GMP, ENRICHMENTS §2.4) |
| target $h_2E$ | $+457$ | (trace growth, expected) |

**Gate decision.**  Exported in v4: the defect suite with the Gram
towers extended to $m\le3$ (`defect_cov_d9`, `defect_gram_d9`,
`h2loc_defect_gram_d9`, `defect_comp_cov_d9`, `defect_comp_gram_d9`,
`defect_cov_minor_d7`, `defect_gram_minor_d8`; ≤4×4 each; the
$m{=}3$-extended Grams re-tested against the am16 escape: min eigs
$-1.11\times10^3$ / $-8.68\times10^3$ / $-6.71\times10^3$ /
$-6.80\times10^2$, all cut-signed; `h2loc_defect_cov` excluded by the
sign rule) and the twelve s4/s6 families
(`jensen_{even_00,even_11}_{s4,s6}` + `h2loc_` + `h2comp_cov_`
copies — carrying the genuinely new rows of §1.2).  Excluded:
`krow_*` (exactly redundant, no unique rows), $e_2$ scalar cut,
deviatoric weight (a *weight*, i.e. a future multiplier-dictionary
entry $w_{\mathrm{dev}}E$ or $(h_2+\kappa w_{\mathrm{dev}})E$, not a
block).

## 3. Staged exports (D2-3)

`toeplitz_export.py --version v4` = v3 (43 families) + the 19 new
families of §2, built through the same interception driver: exact
elimination redone with the new images (no post-hoc file append),
every family trimmed to the run's own label set and verified in-run
(exact rational equality against direct integration + exact PSD at
the ONB and cross measures + float-exactness of every entry), and the
export accounting checked (`dropped_objective_inconsistent = 0`).

### 3.1 Staged files

* **`sdpa_runs/deg14_am_toep4.dat-s`** (staged 2026-08-20): m = 812
  (v3: 810), 57 families (18 new; `defect_cov_minor` trimmed away at
  deg-14), base blocks byte-identical to `deg14_am_toep3` (58 names +
  sizes checked), `dropped_objective_inconsistent = 0`
  (`dropped_dependent_directions = 331`, same elimination behavior as
  v3).  Deg-14 nuance: this file was produced with the $m\le2$ Gram
  towers; the later $m\le3$ extension trims back to $m\le2$ at deg-14,
  so the file is already the fixed point.  Wall solve (200-bit,
  ~20 min):

      cd sdpa_runs && ./sdpa_gmp -ds deg14_am_toep4.dat-s \
          -o deg14_am_toep4.result -p param_200bit.sdpa

  Reported bound = objValPrimal + 2/3.  v3 reference: −1.2081e−5.

* **`sdpa_runs/deg16_am_toep4.dat-s`** (staged 2026-08-20): m = 1267
  (v3: 1265), 128 blocks, 62 families — all 19 new families landed
  (`defect_gram_d9` / `h2loc_` / `comp_` kept the full $m\le3$ tower,
  4×4; `defect_gram_minor_d8` kept $m\le2$; `defect_cov_d9` /
  `comp_cov` trimmed to $m\le1$; s4 families 20+17 rows, s6 families
  14+14 untrimmed) — base blocks identical to `deg16_am_toep3`
  (66 names + sizes), `dropped_objective_inconsistent = 0`
  (`dropped_dependent_directions = 445`).  Selectors staged:
  `sel16_am_toep4_{1em4,1em5,5em6}.dat-s` (m = 1268, 129 blocks +
  slack, exact 40-digit bounds as before).  Wall solve (200-bit) then
  the selector ladder (128-bit):

      cd sdpa_runs && ./sdpa_gmp -ds deg16_am_toep4.dat-s \
          -o deg16_am_toep4.result -p param_200bit.sdpa
      ./sdpa_gmp -ds sel16_am_toep4_1em4.dat-s -o sel16_am_toep4_1em4.result -p param_128bit.sdpa
      ./sdpa_gmp -ds sel16_am_toep4_1em5.dat-s -o sel16_am_toep4_1em5.result -p param_128bit.sdpa
      ./sdpa_gmp -ds sel16_am_toep4_5em6.dat-s -o sel16_am_toep4_5em6.result -p param_128bit.sdpa

  v3 references: wall −1.2252e−8; selector traces 1808.99 / 3055.87 /
  3949.27 at 1e−4/1e−5/5e−6.  Success criterion (§5): wall well below
  −1.2252e−8 and the geometric 3-point law breaking.

## 4. Boundary-layer resummation: the generating-function spec (D2-4)

Derivation only; implementation next round.  This is §6.2 of
LIMIT_EXTRACTION_NOTE made precise, fused with the defect square.

**Why the escape has the $(u^m-w^m)$ structure.**  Write $u=t_1^2$,
$w=t_2^2$.  The measured escape directions are swap-antisymmetric
even-sector monomial differences, e.g.
$(2,4,2)-(4,2,2)=-\,uw\,s^2\,(u-w)$ and
$(5,1,1)-(1,5,1)=t_1t_2s\,(u^2-w^2)$: every one factors through
$u^m-w^m$ for some $m\le3$, i.e. through the defect factor
$(u-w)=t_1^2-t_2^2$.  The kernel-conjugated defect leaf is itself the
combination $\Delta p=8(u^3-w^3)-12(u^2-w^2)+5(u-w)$.  So the
boundary layer is, in each shell, a low-$X\!\cdot\!Z$-modulated
multiple of the *same* difference structure the defect square
carries — the escape is assembling $\Delta$-potential cross terms
from boundary monomials, one shell per degree.  A finite family
tower can chase it (v3 → v4 moves the boundary out by one shell); only
a resummation closes it.

**The dictionary.**  Admissible spin-0 kernels of the (E1)
classification: $f\in\{T_2+\tfrac13,\;T_6+\tfrac13\}$
(docs/SHARP_STRUCTURE.md II.1); note $\Delta f$ kills the $+\tfrac13$,
so $\Delta(T_2+\tfrac13)=T_2(t_1)-T_2(t_2)=2(u-w)$ and
$\Delta(T_6+\tfrac13)=32(u^3-w^3)-48(u^2-w^2)+18(u-w)$.  The kernel
itself decomposes **exactly** over the admissible dictionary:

$$K \;=\; \big(T_2+\tfrac13\big)\;+\;\big(T_6+\tfrac13\big),
\qquad\text{so}\qquad
\Delta p=\tfrac14\,\Delta K=\tfrac14\big(\Delta T_2+\Delta T_6\big),$$

i.e. the defect leaf is the equal-weight sum of the two admissible
generators' defects — the resummation dictionary of §6.2 contains the
defect square as its $(m{=}0)$ member by construction.

**Generating labels.**  For $f$ as above and rational $q\in(0,1)$,
define the modulation generating function of the defect square

$$k_G^{(f)}(q)\;:=\;\sum_{m\ge0} q^m\, \mathbb E\big[s^{2m} D_f^2\big]
\;=\;\mathbb E\Big[\frac{D_f(X,Z)^2}{1-q\,s^2}\Big],
\qquad D_f=U_f(X)-U_f(Z),$$

and its truncations $k_G^{(f)}(q;N)=\sum_{m\le N}q^m\,m^{(f)}_m$ with
$m^{(f)}_m=\mathbb E[s^{2m}D_f^2]$ — the entries of §1.2's Hankel
family ($m^{(f)}_m$ has label degree $2\deg f+2m$, so degree $d$
carries the window $m\le(d-2\deg f)/2$: $m\le1$ (deg 14), $m\le2$
(16), $m\le3$ (18) for $f=T_6+\tfrac13$ — exactly the measured
escape's $k\le3$).  Carried as one new coordinate $y_\tau$ per
$(f,q)$ with the **one-sided truncation cuts** (the theta-atom scheme
of ENRICHMENTS §4.5 transplanted from the azimuthal $q^{n^2}$ axis to
the modulation $q^m$ axis):

$$(\mathrm T^+_{N'}):\quad y_\tau\;\le\;k_G^{(f)}(q;N')+R^{(f)}_q(N'),
\qquad
(\mathrm T^-_N):\quad y_\tau\;\ge\;k_G^{(f)}(q;N),$$

with the exactly-rational tail majorant, from
$|D_f|\le2\sup_{[-1,1]}|f|=\tfrac83$ and $s^2\le1$ pointwise:

$$R^{(f)}_q(N')\;=\;\frac{64}{9}\cdot\frac{q^{N'+1}}{1-q}
\;\ge\;\sum_{m>N'}q^m m^{(f)}_m
\qquad(\text{rational for rational }q).$$

(For the kernel-conjugated tower itself — $f$ replaced by $p$,
$\sup_{[-1,1]}|p|=1$ exactly, attained at $t=\pm1$ — the constant
sharpens to $R^{(p)}_q(N')=4\,q^{N'+1}/(1-q)$; the v4 defect-Gram
rows $m^{(p)}_m\ge0$-with-Hankel are the lower half of this package,
already exported for $m\le3$.)

Eliminating $y_\tau$, the polynomial-label shadow is the window
inequality, for every $N'<N$ within degree:

$$\boxed{\ \sum_{N'<m\le N}q^m\,m^{(f)}_m(y)\ \le\ R^{(f)}_q(N')\ }$$

— a uniform $q$-graded bound on the total high-modulation defect mass
that **no finite-degree polynomial relaxation contains** (nothing in
the algebra couples the modulation axis to a constant).  Homogenized
on a recession ray $r$: $m^{(f)}_m(r)\ge0$ (the v4 Hankel rows) plus
$\sum q^m m_m(r)\le0$ force $m^{(f)}_m(r)=0$ for all $N'<m\le N$ —
the atom kills every escape ray carrying strictly positive modulated
defect-square mass in the window.

**Localized escalation** (what §6.2's "$(X\!\cdot\!Z)^2$ products"
buys): the same construction with $D_f$ replaced by the vector
$D_f\cdot g_\alpha$, $g_\alpha$ monomials in $(t_1,t_2,s)$ within a
reduced cap (per-multiplier generating labels
$k_{G,\alpha}$, majorant $\sup|g_\alpha|\le1$ on the cube so
$R^{(f)}_q$ applies verbatim), and the $h_2$-localized copies
$h_2\,m^{(f)}_m$ (valid: $h_2\in[0,1]$, same majorants) — mirroring
ENRICHMENTS §4.6.  These reach the swap-antisymmetric even-sector
escape *directly*: the defect factor supplies $(u-w)$, the
multipliers supply the $uw s^2$-type dressing measured in §5 of
LIMIT_EXTRACTION_NOTE.

**Exactness discipline.**  $m^{(f)}_m$ label matrices are exact
rationals (same `same_leaf_entry`/`two_leaf_entry` machinery, already
in `build_leafpoly_family`); $q$ and $R^{(f)}_q(N')$ exact rationals;
the cuts are $1\times1$ rows plus one free variable per $(f,q)$ —
implementable in `sdpa_selector`-style post-processing or as
first-class rows in the exporter.  Sharpness knob: replace the crude
$\sup$ bound $\tfrac{64}{9}$ by the exact
$m^{(f)}_m\le\mathbb E[D_f^2]\,\sup s^{2m}\le m^{(f)}_0$ chain, or by
Chebyshev-graded weights $T_m(s)$ (the $|T_m|\le1$ bound is tight at
the collision boundary $s=\pm1$ where the escape's
$T_n(\pm1)=(\pm1)^n$ signature lives).

## 5. Verdict and consequences

1. The equal-potential defect square **is representable and valid**
   in the all-measures cone at every degree — but as an *inequality
   family* it is (in-vocabulary) exactly redundant: the am cone's
   deficit relative to the KKT cone is the missing *identity*, which
   no valid all-measures block can supply at finite degree.  The
   648× KKT advantage is not portable by blocks alone.
2. The finite-degree payload of the defect program is the
   **trim-boundary re-import** (s4/s6 families) — new rows exactly
   where the escape lives — plus dedicated defect blocks that give
   the dual concentrated multipliers along the direction it is
   measurably assembling (whole-matrix negative pairings).
   v4 success criterion (per the wall-law analysis): the deg-16 am
   wall moves well below $-1.2252\times10^{-8}$ and the 3-point
   geometric fit's $c_\infty<0$ breaks; if the wall barely moves, the
   redundancy analysis of §1.2 is confirmed empirically and the
   resummation (§4) is the only remaining route — it is the
   degree-uniform closure of exactly this construction.
3. The deviatoric square $w_{\mathrm{dev}}=(t^2-\tfrac13)^2$ is a
   **weight** (pairing $+0.985$ with the current am escape): queue it
   as a multiplier-dictionary entry (weighted target
   $(h_2+\kappa\,w_{\mathrm{dev}})E$ or a dedicated
   $w_{\mathrm{dev}}$-weighted family), not as a cut.
