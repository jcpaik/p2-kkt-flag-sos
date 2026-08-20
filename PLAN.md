# Plan: exact certificate for the $P_2$ kernel

New agent?  Read [CLAUDE.md](CLAUDE.md) first (environment, hardware
rules, update discipline), then this page top-to-bottom until the
detailed record, then pick the first open item in **Next actions**.

## Status (dashboard)

Goal: an exact certificate of $E(\mu)\ge0$.
Where things stand (details and provenance in §3–§6 below):

| fact | value | where |
|---|---|---|
| best unweighted bound, degree 14 / 18 | $-4.4856\times10^{-3}$ / $-7.5595\times10^{-5}$ | RESULTS.md |
| weighted target $h_2E$ + full localized module, degree 14 | $\mathbf{-2.7909\times10^{-5}}$ ($161\times$ gain; re-verified 2026-08-17 post-regeneration, $-2.7908758\times10^{-5}$) | §5 Stage 2 |
| adding the $e_k(I-A_2)$ cuts to that bound | identical to 12 digits (cuts slack at optimum) | §4 |
| $\varepsilon$-trace law, weighted, deg 14 (new scale, 2026-08-17) | $765\to8134\to13637$ at $\varepsilon=10^{-3},10^{-4},5.3\times10^{-5}$ — $10.6\times$/decade, simple pole ($\varepsilon^{-1.03}$) | docs/SUBCASES_AND_RECORD.md |
| projected-dual escape ray | pure $g_4$; **killed** by the scalar cuts (`--find-ray` infeasible A/B; re-verified 2026-08-17) | §4, docs/ENRICHMENTS.md |
| **the surviving escape, fingerprinted (2026-08-17)** | the $h_2$-**localized theta tower**: $78\%$ $p_2\times$(graph_4/triangle/pair-products), $22\%$ connected graph_4/triangle, $0\%$ pair sector — every scalar $p_k$-invariant is blind to it | docs/SUBCASES_AND_RECORD.md |
| second weight candidate | $e_5(I-A_2)=\det(I-A_2)$: valid, dense-positive, vanishes on the **whole pole–equator stratum**; needs `graph_5` labels | docs/ENRICHMENTS.md |
| theta-atom adjunction (measured 2026-08-17) | localized atoms re-bound the projected dual (`infeasible` A/B — first (E1)-derived mechanism); but **inert on the finite-$\varepsilon$ pole**: theta-augmented selector traces $765.0/8133.57$ = controls to 6 digits | docs/ENRICHMENTS.md |
| **two new exactly-solved sub-cases (2026-08-17)** | axisymmetric measures (rational Handelman-8 certificate, 6/6); antipodal circle-pairs, any latitude (Toeplitz coupling, 27/27, sharp at ONB) | §1, docs/SUBCASES_AND_RECORD.md |
| scalar mode-domination | **falsified exactly** (min $-128159/27456$; cross shadow $-1/192$): the missing constraint is measure-valued Toeplitz joint positivity | docs/SUBCASES_AND_RECORD.md |
| weighted-(E1) arity 2 (2026-08-17) | classified (69/69): admissible leaves = root-equivariant traceless contractions; $\hat C_n$ tower survives, $\hat S_n$ dies; third-order law forces unweighted-(E1) bases on the $h_2$loc layer | docs/SHARP_STRUCTURE.md |
| **split verdict on the pole (2026-08-17)** | the selector growth direction is $80.6\to95.7\%$ **inside** the sharp-face admissible span as $\varepsilon\downarrow$: the pole is non-attainment *inside* the sharp cone — projections cannot block it, only cuts on the admissible tower | docs/SHARP_STRUCTURE.md |
| degree-18, sharp-face projected | **reachable**: `deg18_h2w_h2all_am_we1.dat-s`, m = 1015 (< unprojected deg-16's 1243); deg16 le1 m = 1173 | sdpa_runs/, docs/SHARP_STRUCTURE.md |
| weighted deg-16 at 128-bit | spurious pdINF after 37 min — **use 200-bit** for weighted bound solves at degree ≥ 16 | sdpa_runs/REGEN_NOTES.md |
| **weighted degree-16 bound (2026-08-17, 200-bit, 71 min)** | $\mathbf{-1.3140\times10^{-6}}$ — $21\times$ better than weighted deg-14, $57\times$ better than *unweighted deg-18*; the super-geometric weighted decay is confirmed; deg-16 wall now $\approx1.3\times10^{-6}$ | `deg16_h2w_h2all_200.result` |
| **Jensen/covariance cut tower vs the pole (2026-08-17/18)** | v1/v2/v3 traces at $10^{-4}$: $7035/1544/1473$ vs control $8134$; window exponent $1.03\to0.55$; **dual bound $-2.79\times10^{-5}\to-3.6605\times10^{-6}$ ($7.6\times$) at degree 14**; below $5\times10^{-5}$ the law re-steepens — the escape *saturates* the covariance inequalities (deterministic-fiber kernel) | docs/SUBCASES_AND_RECORD.md §4–5, docs/ENRICHMENTS.md |
| $e_5(I-A_2)$ implemented (2026-08-18) | exact 35-label `graph_5` expansion; $0\le e_5\le(4/5)^5$ two-sided validity; **first scalar invariant to move a bound** ($-1.920\to-1.832\times10^{-4}$ double); deep-escape pairing cut-signed ($-1.26\times10^5$); grows positive-quadratically off the stratum | docs/ENRICHMENTS.md |
| degree-16 selector law | $4073.7\to14068.7$ at $10^{-4}\to10^{-5}$ ($3.45\times$/decade); +Jensen v1 stack: $3752.6/8151.9$ ($2.17\times$, exp 0.34) | sdpa_runs/sel16_* |
| **deg-16 + v3 stack (43 families, 2026-08-18)** | trace law $959.5/1546.0/1805.3$ at $10^{-4}/10^{-5}/5\times10^{-6}$ — **exponent 0.21, no re-steepening**; dual wall $\mathbf{-5.2504\times10^{-9}}$ ($250\times$ over plain deg-16; gap $1.8\times10^{-14}$) | docs/SUBCASES_AND_RECORD.md §7–8 |
| bare $e_5$ cut at GMP | inert (bound and $10^{-4}$ trace = baseline); $e_5$'s route is the localized-module weight | docs/ENRICHMENTS.md, REGEN_NOTES |
| all-measures (proof-carrying) selector law | $2366.8\to13873.3$ ($5.86\times$/decade) | sdpa_runs/sel_am_* |
| **proof-carrying deg-16 + v3 (2026-08-18/19)** | wall $\mathbf{-1.2252\times10^{-8}}$ (2.3× the KKT copy — KKT toggles nearly redundant); ladder $1809/3056/3949$ (exp 0.23/0.37) | docs/UNPROJECTED_ESCAPE_NOTE.md §9 |
| **FIRST EXACT RATIONAL CERTIFICATE (2026-08-19)** | $h_2E\ge-1.4824\times10^{-8}$ for every antipodal measure — exact, self-contained checker `verify_h2E_bound.py` (stdlib-only), PASS 104 s; the exactification pipeline is proven end-to-end | certificates/, docs/EXACT_CERTIFICATE_NOTE.md |
| sharp-face (we1) route | **closed**: we1+v3 selector pUNBD at $10^{-5}$ (cannot certify at any trace) | REGEN_NOTES |
| **degree-18 + v3 (2026-08-19)** | wall $\mathbf{-2.21142\times10^{-14}}$ (pdOPT; cascade accelerating: $700\times\to240{,}000\times$ per two degrees); trace law exponent 0.32 at depth — **no plateau: zero is a limit, not attained** | docs/UNPROJECTED_ESCAPE_NOTE.md §10 |
| **exact certificate #2 (2026-08-19)** | $h_2E\ge-2.745\times10^{-14}$ (KKT-inclusive relaxation; gradient rows load-bearing), checker PASS 139 s; pipeline upgraded (iterative-refinement $\lambda$, fraction-free exact PSD) | certificates/, docs/EXACT_CERTIFICATE_NOTE.md §6 |
| **proof-carrying degree-18 (2026-08-20)** | wall $\mathbf{-1.43335\times10^{-11}}$ (pdOPT); KKT/am gap widens with degree ($2.3\times$ at 16, $648\times$ at 18) | sdpa_runs/deg18_am_toep3.* |
| **CAPSTONE EXACT CERTIFICATE (2026-08-20)** | $h_2E\ge-1.7356\times10^{-11}$ for **every** antipodal probability measure — unconditional, zero KKT content, reduction lemma composes; checker PASS 83 s; $854\times$ beyond the deg-16 record | certificates/h2E_geq_deg18_am.json.gz |
| (E1) complementary-slackness equations | solved in closed form, one- and two-root layers, 54 exact checks | §6, docs/SHARP_STRUCTURE.md |
| wrapper lemmas for the final proof | drafted, machine-checked | docs/FOUNDATIONS.md |
| conjecture status in the literature (swept 2026-08-17) | open; kernel untouched; nearest relative (BMN 2409.16508 §8) also open | §7 below |

**Normalization note (2026-08-17).**  The solver target was changed
from the legacy $(3/16)E$ to $E$ itself
($E=-\tfrac43+20p_2-48p_4+32p_6$; weighted target
$h_2E=\tfrac23-12p_2+24p_4-16p_6+30p_2p_2-72p_2p_4+48p_2p_6$).  Bounds
in this file are converted (legacy $\times16/3$).  New exports carry
shifts $-4/3$ ($E$) and $2/3$ ($h_2E$).

**Regeneration note (2026-08-17, later the same day).**  The gitignored
`sdpa_runs/` and `.venv` were found wiped and were fully regenerated
from code in the new normalization (params from the
`~/opt/sdpa-gmp/param.sdpa` template; venv from requirements.txt; all
exports re-made; 22/22 tests pass).  The degree-14 weighted anchor
re-solved to $-2.7908758\times10^{-5}$, matching the pre-wipe record to
8 digits, so regenerated data is trustworthy.  All legacy-scale
artifacts (selector traces $330/5529$, legacy ray files) are gone;
the regenerated trace law and fingerprints (new scale) supersede them.
Keep small critical artifacts committed (docs notes, and
`sdpa_runs/REGEN_NOTES.md` describes the reproduction commands).

Standing conclusions a newcomer should not re-derive: no finite-degree
certificate for bare $E$ exists (non-attainment, $1.07/\varepsilon$ (legacy scale)
trace pole); the $h_2$-weighted reformulation is lossless (reduction
lemma, §5) and strictly stronger numerically; proof-carrying runs must
use the all-measures cone (§5 composition caveat); scalar functions of
invariants ($\sqrt{h_2}$-adjunction, $h_2$ powers, $e_k(I-A_2)$ as
weights) are measured dead — see the sign rule in
docs/ENRICHMENTS.md.

## Next actions (work queue, in order)

*(2026-08-20: the Jensen/fiber-Toeplitz campaign and the degree-14/16/18
ladder are complete — see the dashboard, docs/UNPROJECTED_ESCAPE_NOTE.md
and docs/LIMIT_EXTRACTION_NOTE.md.  The proof-carrying tower decays
purely geometrically (986×/855× per two degrees), so degree escalation
alone cannot land on zero: the closing construction below is the queue.)*

1. **Adjoin the equal-potential defect as an all-measures object**
   (LIMIT_EXTRACTION_NOTE §6.1): the KKT tower's dominant certificate
   object is $h_2\mathrm{loc}\,[K(X\!\cdot\!Y)-K(Z\!\cdot\!Y)]^2$
   (+ its $(X\!\cdot\!Z)^{2m}$-modulated copies, m = 1, 2).  Build the
   valid all-measures analogue (variance-form labels via the one-sided
   truncation principle), pair-test against the stored escape data,
   export, solve.
2. **Resum the even-sector boundary layer** (§6.2): the residual escape
   is a swap-antisymmetric $q$-profile in leaf degree at fiber
   modulation $k\le3$; the finite dictionary is the generating-function
   labels $k_G$ over the two admissible spin-0 kernels
   ($T_2+\tfrac13$, $T_6+\tfrac13$) and their $(X\!\cdot\!Z)^2$
   products, with exactly-rational one-sided truncation cuts.
3. **Cheap pre-test** (§6.3): the deviatoric-square weight
   $w=(t^2-\tfrac13)^2$ — pair against the escape fingerprints before
   any solve.
4. **e5-localized module** (docs/E5_WEIGHT_NOTE.md): the
   $(h_2+\kappa e_5)E$ multiplier target still needs its localization
   layer; deg-16-scale build.
5. **Exactification pipeline on every new wall** (now sub-minute per
   solve, docs/EXACT_CERTIFICATE_NOTE.md): keep the certificates/
   ladder current; the moment a wall lands at ≥ 0 in the all-measures
   cone, the reduction lemma closes the conjecture.

## Map of the repository

Code: `sos_search.py` (hierarchy, exports, `--find-ray`, all toggles —
see docs/FOUNDATIONS.md), `solve_e1.py` (exact (E1) solver +
projections), `sdpa_selector.py` (min-trace selector transform),
`sdpa_extract.py` (parse GMP results), `verify_certificate.py` /
`verify_exact_structure.py` (independent exact checkers).

Documents: strategy — docs/ENRICHMENTS.md; session notes with
full derivations — docs/ENRICHMENTS.md (cuts/sign rule),
docs/ENRICHMENTS.md (the $161\times$ mechanism); theory
reference — docs/FOUNDATIONS.md, docs/FOUNDATIONS.md,
docs/FOUNDATIONS.md, docs/SHARP_STRUCTURE.md,
docs/SHARP_STRUCTURE.md, docs/SHARP_STRUCTURE.md,
docs/FOUNDATIONS.md; measurement log — RESULTS.md (committed
record) and the checklists below (live).

Data: `sdpa_runs/` (gitignored) — solver binary, parameter files, all
exported problems, results, ray JSONs, exact expansions.

---

Everything in §0–§2 is verified in exact arithmetic by
`python3 verify_exact_structure.py` (25/25 checks).  §3–§6 are the
detailed live record; the checklists there are the provenance for the
dashboard above.

## 0. Clean reformulations

The kernel is, with $t=\cos\theta$,

$$K(t)=32t^6-48t^4+20t^2-\tfrac43 \;=\; \cos 6\theta+\cos 2\theta+\tfrac23 .$$

Legendre expansion $K=\tfrac{32}{105}P_0+\tfrac87P_2-\tfrac{384}{385}P_4+\tfrac{512}{231}P_6$:
**exactly one negative coefficient** ($\ell=4$).

Equivalent forms of the conjecture $E(\mu)\ge0$:

1. **Squared form.** $E=2\iint W(x\cdot y)-\tfrac43$, where
   $$W(t)=T_3(t)^2+t^2=\cos^2 3\theta+\cos^2\theta=16t^6-24t^4+10t^2 .$$
   Claim: $\iint W\ge\tfrac23$.
2. **Frame form.** $\iint T_3(x\cdot y)^2\,d\mu\,d\mu+\big\|\int xx^{\mathsf T}d\mu\big\|_F^2\ \ge\ \tfrac23 .$
3. **Isotropic form.** If $A_2=0$ then $E\ge0\iff 11-36A_4+80A_6\ge0$.
4. **Double-angle form.** $W=1-u+2u^3=(1+u)\big(u^2+(1-u)^2\big)$, $u=\cos2\theta=2t^2-1$.

Neither summand of (2) suffices alone: $\min\iint T_3^2=0.2678884\ldots<\tfrac13$
in general, and $=\tfrac13$ exactly under isotropy — but the isotropic version of
(2) is algebraically *identical* to the target, so there is no free split.

## 1. Exactly solved sub-cases

* **Measures on a great circle.** $E(\nu)=|\hat\nu(2)|^2+|\hat\nu(6)|^2+\tfrac23\ \ge\ \tfrac23>0$.
  Immediate from $K=\cos6\theta+\cos2\theta+\tfrac23$; the planar Fourier basis has
  no negative coefficient. This is the only place the problem is easy.
* **Axisymmetric measures (2026-08-17, new).**  $E\ge0$ for every
  rotation-invariant $\mu$ — equivalently the kernel
  $c_0(s,t)$ ($s=a^2$; azimuthal average of $K$) is copositive on
  $[0,1]$.  Exact rational certificate: symmetric-Handelman degree 8
  (the feasibility threshold; degrees 3–7 infeasible) + a PSD Gram
  part whose form vanishes at the zero measure
  $\sigma^*=\tfrac23\delta_0+\tfrac13\delta_1$, every Handelman term
  corner-vanishing.  Verified independently in pure rational
  arithmetic: `cylinder_cert.py m1 --verify
  sdpa_runs/cylinder/m1_certificate.json` (6/6).  See
  [Cylindrical domination](docs/SUBCASES_AND_RECORD.md).
* **Antipodal circle-pair measures (2026-08-17, new).**  $E\ge0$ for
  every antipodal measure supported on a pair of antipodal circles of
  *any* latitude with *arbitrary* fibers, with equality exactly on the
  ONB orbit.  Proof = scalar mode domination + the Toeplitz/Herglotz
  coupling $|\zeta_{2k}|\ge2|\zeta_k|^2-1$ for the mode pairs
  $(2,4)$ and $(3,6)$ (double-angle trigonometry of the fiber), seven
  Sturm-count facts, and the sharp factorization
  $4v_6+v_3=4(1-s)^3(3s-1)^2(6s+1)$.  Machine-checked:
  `cylinder_cert.py pair-theorem` (27/27).  Notably *pure scalar*
  domination is provably insufficient even here (cross shadow
  $-1/192$ at $s=\tfrac12$, ONB shadow $-64/729$ at $s=\tfrac13$):
  the mode-$2k$/mode-$k$ coupling is essential — the first exact
  confirmation that the missing information is Toeplitz-type joint
  positivity, not scalar bounds.
* **Pole–equator measures.** For $\mu=w\,\delta_{\pm e}+(1-w)\nu$ with $\nu$ on $e^\perp$:
  $$\boxed{\;E(\mu)=6\Big(w-\tfrac13\Big)^2+(1-w)^2\Big(|\hat\nu(2)|^2+|\hat\nu(6)|^2\Big)\;}$$
  an exact SOS identity. Zero iff $w=\tfrac13$ and $\hat\nu(2)=\hat\nu(6)=0$.
* **ONB KKT certificate.** For $\mu_{\mathrm{ONB}}$: $E=0$ and the rooted potential is
  $U(z)=32\,z_1^2z_2^2z_3^2\ \ge\ 0$ exactly.

## 2. The minimizer set (the central obstruction)

Numerically, $\min E=0$ over 3–30 atoms, and **every** converged global
minimizer (100+ runs) is of pole–equator form. So the minimizer set is
conjecturally the $SO(3)$-orbit of the codimension-4 convex family above —
infinite-dimensional. This matches Bilyk–Matzke–Nathe (arXiv:2409.16508, §8),
who call these "pole–equator measures".

Consequences, both verified:

* **No 2-point (Delsarte/Yudin) LP certificate exists at any degree.**
  Pole + Haar-equator has $E=0$ and a distance set equal to all of $[-1,1]$,
  forcing $h\equiv K$, contradicting $\hat K_4<0$. Imposing isotropy does not help.
* Facial reduction is **already complete**: pole–equator faces of regular order
  $m$ annihilate the target for every $m\ne1,3$ (checked exactly), but adding
  orders $4,5,7$ to the published ONB + continuous reduction leaves the reduced
  dimension unchanged at 480. So the published audit was not missing faces.

**Cylindrical decomposition.** Writing $a=\cos\psi_x,\ b=\cos\psi_y$ and
$p^2=(1-a^2)(1-b^2)$,
$$K(x\cdot y)=\sum_{k=0}^{6}c_k(a,b)\cos\big(k(\varphi_x-\varphi_y)\big),$$
$$c_6=p^6,\quad c_5=12abp^5,\quad c_4=6p^4(11a^2b^2-a^2-b^2),$$
$$c_3=4abp^3(55a^2b^2-15a^2-15b^2+3),\ \ c_2,\ c_1,\ c_0 \text{ as computed}.$$
$c_5,c_6$ are rank-one positive kernels; $c_0,c_2,c_3,c_4$ are **not** PSD
(determinants $-34992$, $-6480$, $-60$, $-1$). So this decomposition alone does
not certify — as it must not, since $K$ is not a positive-definite kernel.
*(Correction 2026-08-17: $c_1\not\equiv0$ — the odd mode is present;
exact formulas for all $c_k$ are hard-coded and 16/16-verified in
`cylinder_cert.py`.)*

## 3. State of the SDP search

* Baseline degree-14 five-point dual reproduced exactly: $-2.548701208\times10^{-7}$.
* **All** stronger settings plateau in $[-2.45\times10^{-5},-2.4\times10^{-7}]$ with
  *no* monotone improvement — root-factor 3 ($-5.60\mathrm e{-6}$), arity 6
  ($-2.41\mathrm e{-7}$), degree 16 ($-1.66\mathrm e{-5}$), `--gram-module`
  ($-1.08\mathrm e{-5}$), `--higher-rank-matrices` ($-6.19\mathrm e{-7}$).
  Adding blocks can only *raise* a dual minimum, so these are solver noise:
  **the dual optimum is 0 to solver accuracy.**
* Fixes made to `sos_search.py` while diagnosing:
  - `--eliminate-free` **breaks the primal**. Without it the primal solves
    (degree 10, $\varepsilon=0.1$: trace $0.989$); with it MOSEK fails even on
    trivially feasible instances. The rational RREF quotient rows are too
    ill-conditioned in floating point. The published exactification audit used
    this flag.
  - the primal branch had no block (column) rescaling, only row rescaling; added.
  - added `--max-margin CAP`: maximise $\lambda$ s.t. $Q_\beta\succeq\lambda I$
    and $\sum\operatorname{tr}Q_\beta\le\text{CAP}$. A strictly interior Gram
    solution is the prerequisite for rational rounding; trace minimisation is not.

## 4. Roadmap to a proof

What separates the current state from a complete proof is one exact
certificate object plus routine wrapper lemmas.  The proof does **not**
require the minimizer classification (§2): the zero-set knowledge only
guides the search; a found certificate is verified self-contained.

Decision order:

- [x] **Attainment test: measured, pole survives** (degree 14, GMP,
      both cones): $\operatorname{tr}C_w$ grows $13$–$17\times$ per
      $\varepsilon$-decade ($330\to5529$ KKT cone, $704\to9044$
      all-measures cone at $10^{-4}\to10^{-5}$); both degree-14 weighted
      bounds certified in $[-5.3\times10^{-5},0)$.  No finite-degree-14 weighted
      certificate; the escape is the pure $g_4$ mode (below).
- [ ] **If attained: exactification.**  Max-margin primal (interior Gram
      point) on the reduced face → rational rounding → independent exact
      verification (`verify_certificate.py`).  With the reduction lemma
      (§5, proved) this closes the conjecture.
- [ ] **Weighted (E1)** (highest-leverage parallel step): the zero set of
      $h_2E$ contains *all isotropic measures*, so every square of a sharp
      weighted certificate must vanish there — the winning
      `--h2-localized-all` module has exactly this form.  Derive the full
      weighted admissible cone as in §6, project, and push the weighted
      hierarchy to degree 18–20 at the reduced size.
- [ ] **The pole survived $h_2$ (measured): the primary route is now the
      spin-2 operator gap.**  The projected-dual escape ray is a *pure
      $g_4$ direction* ($g_4\cdot\text{ray}=5.35$, all other $g_\ell$ and
      $E[\det\mathrm{Gram}_3]$ pairings $=0$; $-\frac{72}{385}\times5.35=-1.000$
      reproduces the normalization exactly) — the $\hat K_4<0$ mode is the
      whole escape.  The constraint that kills it while staying
      sharpness-compatible is the SO(3) operator bound of
      [Structure](docs/FOUNDATIONS.md) §4: $I-A_2\succeq0$ with
      $A_2=\int\pi_2(\rho_x)\,d\mu$.  Verified exactly (`solve_e1.py`
      Part E): the axial quadrupole $\mathrm{diag}(1,1,-2)$ is an
      eigenvalue-1 eigenvector of $A_2(\mu^*)$ **identically across the
      whole zero family**, and $A_2(\mathrm{ONB})$ has spectrum
      $\{1,1,-\frac13,-\frac13,-\frac13\}$ — the gap is active at
      sharpness (unlike the rejected $|M|\le1$ box cuts).  Next: implement
      operator-gap localizing blocks
      $G_{\alpha\beta}=E[v_\alpha^{\mathsf T}(I-A_2)v_\beta]$ over
      equivariant spin-2 test features $v_\alpha$ (entries are existing
      multigraph labels), re-run `--find-ray` (the $g_4$ ray must die),
      then the bound sweep and the $\varepsilon$-trace law over
      cone + operator gaps.  The theta atoms (§6) remain the fallback for
      whatever escape survives.  *First measurements (implemented as
      `--spin2-operator-gap`, multiplier $4ab(c-ab)$): valid, ONB-audit
      clean with the minus block PSD-with-kernel; but the depth-1 feature
      localization does not yet pin the ray ($g_4$-escape persists, norm
      $60\to156$) — escalate feature depth / $A_2$-word blocks, and
      re-measure in GMP.*

      **Update: the scalar shadow of the gap already kills the escape.**
      The elementary invariants $e_k(I-A_2)\ge0$ expand exactly into
      $\le k$-sample labels via
      $\operatorname{tr}A_2^k=E[\chi_2(\rho_{x_1}\cdots\rho_{x_k})]$,
      $\chi_2(R)=(\operatorname{tr}R)^2-\operatorname{tr}R-1$
      (`sdpa_runs/gap_invariant_expansions.json`;
      $\operatorname{tr}(I-A_2)=4$ identically).  Ray pairings:
      $e_2=-9.79$, $e_3=-19.26$, $e_4=-5.85$ (calibration $g_2=0.005$,
      $E[P_4]=5.353$ both reproduced), so each cut individually forbids
      the $g_4$ ray; none qualifies as a second *weight* (face values
      nonzero on the continuous strata except $e_4$ at the ONB alone,
      and every pairing has the cut sign, not the weight sign) — see
      [Exact zero program](docs/ENRICHMENTS.md) §5.  Implemented
      as `--gap-scalar-cuts` (three $1\times1$ blocks;
      $e_2=6+6p_2-8p_4$, i.e. $p_4\le\tfrac34(1+p_2)$).  **A/B
      `--find-ray` on the (E1)-projected problem: without cuts
      `optimal_inaccurate`, norm $55.7$,
      $E[P_4]\cdot\text{ray}=5.3498$ (control reproduced,
      `sdpa_runs/ray_nocuts_deg14.json`); with cuts `infeasible`
      (`sdpa_runs/ray_gapcuts_deg14.json`) — no improving recession
      direction of any kind survives, so the projected dual is bounded
      again.**  First stacked measurement (degree 14, 200-bit,
      `deg14_h2w_h2all_cuts`): the weighted bound with cuts is
      $-2.7908757\times10^{-5}$ — **identical to the no-cuts value to 12
      digits**, i.e. the cuts are inactive (slack) at the finite-degree
      optimum.  As expected from the ray analysis: their content is the
      $\varepsilon\to0$ attainment structure, not the
      $\varepsilon$-away-from-zero bound.  $\varepsilon$-trace law on
      the cut problem, measured (`sel_cuts_1em5`, 128-bit, terminated
      pFEAS at 37% duality gap):
      $\operatorname{tr}C(10^{-5})\in[3430,4980]$ versus $\le330$ at
      $10^{-4}$ (cuts enlarge the certificate cone, so the no-cuts
      $330$ is an upper bound there) — still $\ge10\times$ per decade.
      **The cuts kill the projected recession ray but not the
      finite-$\varepsilon$ pole**: the $\varepsilon$-certificates grow
      through directions of the *unprojected* cone invisible to the
      (E1)-projection analysis.  Next diagnostic: fingerprint the trace
      mass of the $\varepsilon$-certificates themselves (parse the
      selector `yMat` blocks via `sdpa_extract.py`, compare block norms
      across $\varepsilon$) to identify the unprojected escape, then
      match it with a weight or a cut by the sign rule of
      [Gap cuts note](docs/ENRICHMENTS.md).  Queued behind it: the
      deg16 weighted solve (memory-bound on this machine).
- [x] **Wrapper lemmas drafted** — [Wrapper lemmas](docs/FOUNDATIONS.md):
      antipodal reduction, minimizer existence, first/second-order
      conditions, identity validity, certificate-validity lemma, the
      weighted-composition requirement with both fixes, and the theta-atom
      bounds ($|\hat C_n|\le\frac43+4n$, $|\hat S_n|\le\frac43+12n$ via the
      closed form $C_n=C_0T_n(s)+2t_1(t_2-st_1)U_{n-1}(s)$, machine-checked)
      with valid one-sided truncation cuts.  Remaining: paper-grade
      write-up polish and the $W$-KKT encoding (L6(ii)) if route (ii) is
      ever needed.
- [ ] Superseded items: the isotropy caveat is closed by the §5 reduction
      lemma; the degree-14 max-margin go/no-go is subsumed by the
      weighted attainment test above.

## 5. Enriched-algebra program (exact zero without degree escalation)

Worked out in [Exact zero program](docs/ENRICHMENTS.md).  The
$\operatorname{tr}C(\varepsilon)\sim1.07/\varepsilon$ simple pole says the
sharp certificate escapes along a recession direction $Y_0$ — classical SOS
non-attainment, repaired by a positive multiplier vanishing to first order
along the leak, i.e. $h_2=(3p_2-1)/2$.

- [x] Reduction lemma: $h_2E\ge0$ for all measures $\Rightarrow$ $E\ge0$ for
      all measures (anisotropic density + weak-\* continuity).  No isotropic
      stratum remains.
- [x] **Composition caveat (audited 2026-08-14).**  The reduction lemma
      needs $h_2E\ge0$ for *all* measures, but the Stage-1/2 runs include
      the KKT-only families (`--gradient --potential --hessian
      --global-tangent-gaps` and their $p_2$-shifted copies), valid only at
      critical measures.  A sharp certificate over that cone proves only
      "any counterexample minimizer is isotropic" — vacuous, because the
      zero family itself is isotropic ($h_2=0$ there).  Proof-carrying
      numbers must come from (i) the **all-measures cone** (drop the four
      KKT toggles; every remaining family is an identity or a valid
      square), or (ii) KKT re-encoded for $W=h_2E$ itself (first variation
      $\Phi^W_\mu=h_2(\mu)\,U_\mu+E(\mu)\,u^{h_2}_\mu$, encodable in the
      label algebra).  The KKT-inclusive numbers remain useful as
      attainment diagnostics only.
- [x] All-measures weighted cone measured (MOSEK double, degree 14,
      flags+harmonics+rank only): plain $-1.07\times10^{-2}$; with
      `--h2-localized-all` $\mathbf{-3.8\times10^{-4}}$ — stronger than
      the Stage-1 KKT-inclusive numbers.  Dropping the KKT families costs
      little; the proof-carrying cone is the strong one.
- [x] Weighted (E1) at arity 0/1 (`solve_e1.py` Part D, exact, 54 checks
      total): the zero set of $h_2E$ contains every isotropic measure, so
      each *pure* square of a sharp all-measures certificate must have
      leaves that are pure degree-2 spherical harmonics in the leaf
      variable, while $h_2\times$block squares are exempt.  Surviving
      pure-square layer: spin 0 $\{t^2-\tfrac13\}$, spin 1 $\{t\}$, spin 2
      $\{1\}$ — exactly the deviatoric-tensor flag — and pair layer
      $\{3p_2-1\}$.  `--h2-localized-all` is the *forced* certificate
      structure, not a heuristic.  Weighted arity-2 classification: open.
- [x] `--h2-weighted-target` implemented: $h_2E$ is polynomial in the existing
      label algebra (products $p_2p_j$ are four-point disconnected labels).
- [x] Stage 1: degree-14 SDPA-GMP bounds for $h_2E$ (200-bit,
      `epsilonStar 1e-25`, pruned nine-toggle base):

      | module | m | bound for $h_2E$ |
      |---|---:|---:|
      | none (products only via `empty_type_flag`) | 398 | $-2.0839555\times10^{-3}$ |
      | + `--h2-localized-flags` (one-root $h_2\times$flag) | 414 | $-1.2031830\times10^{-3}$ |
      | + $p_2\times$(all families) | 790 | $-2.0839555\times10^{-3}$ (identical to none) |
      
      Unweighted baseline $-4.4856\times10^{-3}$.  Lesson: pure
      $p_2$-multiplied PSD blocks touch only fresh product labels and are
      **vacuous** for the dual bound; localized blocks must carry the
      $-\tfrac12$ constant part of $h_2=\tfrac{3p_2-1}{2}$ so that product
      and base labels appear in one PSD constraint.  Implemented as
      `--h2-localized-all` ($h_2\times$every PSD family,
      $p_2\times$every equality family).
- [x] Stage 2: `--h2-localized-all` at degree 14 (m = 790, 60 iterations,
      41 min): bound $-2.7909\times10^{-5}$ — a $161\times$ improvement over
      the unweighted degree-14 bound at the same degree and arity, and
      $2.7\times$ better than the *unweighted degree-18* value
      ($-7.5595\times10^{-5}$, m = 1089).  The full localized module is the
      strongest per-variable configuration measured so far in the
      non-isotropic hierarchy.
- [x] Stage 3 (both points measured, SDPA-GMP 128-bit, degree 14, full
      localized module): $\operatorname{tr}C_w(10^{-4})=330.125$,
      $\operatorname{tr}C_w(10^{-5})=5529.29$ (pdFEAS, gap $1.5\times10^{-9}$).
      **Verdict: the pole survives one $h_2$** — growth $16.7\times$ over the
      decade (local exponent $\varepsilon^{-1.2}$), against $\approx1\times$
      for attainment.  The $h_2$ weight only collapses the pole coefficient
      ($19\times$ below the unweighted law at $10^{-5}$).  Caveat: the
      degree-14 weighted wall is $-2.79\times10^{-5}$, so the $10^{-5}$
      point is partially edge-dominated; the clean confirmation is a
      degree-16 selector (wall pushed to $\sim10^{-7}$) after the queued
      deg16 dual solve.  Consistent with the projected-dual recession ray
      (§6): the residual escape is $h_2$-orthogonal with a collision-type
      tail, so no power of $h_2$ alone can repair it — the next denominator
      must be collision-graded, or the certificate needs the theta atoms.
- [x] ~~Queued: `sel_h2all_1em5` selector~~ — measured (Stage 3 above:
      $5529.29$).  The deg-16 dual solve remains queued as **Next
      actions #2** at the top of this file.
- [ ] Stage 3 completion: if any trace law plateaus, run §4 max-margin on
      the weighted problem and start rational rounding
      (`verify_certificate.py` pipeline unchanged) — **Next actions #5**.
- [ ] ~~Next denominator factor ($h_2^2$, $h_2+\beta h_4+\gamma h_6$)~~ —
      superseded: the residual escape is $h_2$-orthogonal and
      collision-graded (§6 ray data), so no power of $h_2$ helps, and the
      measured $e_k(I-A_2)$ candidates all carry the cut sign, not the
      weight sign (§4 update; docs/ENRICHMENTS.md).  A new weight, if
      any, must come out of the **Next actions #1** fingerprint.
- [ ] Fallback: entire-kernel dictionary (heat kernel at rational
      temperatures, Gegenbauer $(1-2rt+r^2)^{-1/2}$), projected onto the
      (E1) circle-mode subspace; one-sided truncation cuts keep the SDP
      exactly rational.  Square-root adjunction of scalars ($s=\sqrt{h_2}$)
      is *dominated* (it reproduces the Cauchy–Schwarz closure already in
      the spin-2 block) — see the negative result in the program doc §3.

## 6. The (E1) equations solved (Route 3: solve the pattern, don't guess it)

Full derivation and 54 exact checks: [E1 admissible](docs/SHARP_STRUCTURE.md)
and [Two-root generators](docs/SHARP_STRUCTURE.md),
`python3 solve_e1.py`.  Instead of extrapolating solver coefficients, the
complementary-slackness equations (E1) are solved in closed form; the
finite-degree pattern *is* the solution.

- [x] **One-root leaves classified, all spins, all degrees**: spin 0
      $\{T_2+\frac13,\,T_6+\frac13\}$, spin 1 $\{t,U_5\}$, spin 2
      $\{1,(4t^2-1)^2\}$, spin 3 $\{t^3\}$, spin $\ge4$ empty — a
      7-dimensional layer independent of degree, with
      $K=(T_2+\frac13)+(T_6+\frac13)$ and spin-1 = $\partial_\theta$(spin-0).
- [x] **Unrooted families**: pair flags survive on $\{3p_2-1,\ T\}$ only;
      `harmonic_flag_l` survive exactly on the spin-0 radials; $g_2\equiv0$
      on the family, $g_{\ge4}$ dead as sharp-certificate terms (kept as
      valid inequalities), uniform global gap dead by (E3).
- [x] **Two-root leaves solved structurally**: circle modes $\{0,\pm2,\pm6\}$
      with two free $\delta$-profiles; canonical modulated generators
      $\hat C_n=C_n+\frac13T_n(s)$, $\hat S_n=S_n+\frac13T_n(s)$ admissible
      for every $n$; profile-map kernel $=\{DQ-\frac13(1-s^2)Q(0,0,s)\}$
      (dimension identity $52=26+26$ at degree 8); orientation sector nearly
      free (closed-form codimension).  Exact dimension tables per solver
      sector (`solve_e1.py --table`).
- [x] **`--e1-project` in `sos_search.py`** (+ `--e1-project-families`),
      exact rational bases via `solve_e1.py --export-projection`; ONB audit:
      every projected flag block vanishes at the ONB.
- [x] **Measurement (degree 14, scaled)**: baseline $-4.5\times10^{-3}$;
      one-root projection alone costs only $\sim2\times10^{-4}$ (the
      hierarchy never used the dead one-root directions); two-root alone
      $-2.56\times10^{-2}$; one-root+two-root $-0.869$; full projection
      **unbounded** at degrees 14 and 16 — the $\varepsilon>0$ escape runs
      through (E1)-dead directions of both layers jointly.
- [x] Recession ray extracted (`--find-ray`, CLARABEL; data in
      `sdpa_runs/ray_projected_deg14.json`): the escape is $h_2$-orthogonal
      ($p_2\approx0.003$) with a flat $p_4..p_{14}$ collision-type tail,
      mass 40% triangle / 29% graph_4 / 28% pair — high-modulation two-root
      content, i.e. the truncated theta tower.  (Qualitative: double
      precision, `optimal_inaccurate`.)
- [ ] Exact $\varepsilon=0$ feasibility over the projected cone + theta
      atoms $\Theta_q^{(2)}=\sum_n q^{n^2}\hat C_n$,
      $\Theta_q^{(6)}=\sum_n q^{n^2}\hat S_n$ at rational $q$, with
      one-sided truncation cuts (each tail is a sum of admissible PSD
      blocks) — the dictionary route restricted to the solved ansatz.
- [ ] Combine with §5: the $h_2$-weighted target's certificate is also
      subject to (E1) with the weighted zero set; re-derive the admissible
      leaves for $h_2E$ and project the winning `--h2-localized-all` module.
- [ ] (E3)/(E4) analogues for the multiplier layers ($\rho$ supported on
      $\{U=E\}$, $B$ against the support Hessian).

## 7. Literature position (swept 2026-08-17)

Full annotated sweep in the session record; verdict: **the conjecture
is open** and this exact kernel appears in no indexed work; the parent
problem (Bilyk–Matzke–Nathe arXiv:2409.16508 §8, geodesic projective
energies) is likewise unresolved, with the same pole–equator family
appearing there empirically.  Methodological anchors worth citing in
any write-up:

* Reznick (1995) / Scheiderer — positive multipliers repairing SOS
  non-attainment: the classical frame for the $h_2$-weighted route
  (docs/ENRICHMENTS.md §1–2 describes exactly this mechanism).
* Klep–Magron–Volčič, *SOS certificates for polynomial moment
  inequalities* (arXiv:2306.05761, FoCM 2025) — certificates for
  polynomial expressions in measure moments; structurally our setting
  ($E$ polynomial in $p_2,p_4,p_6$), including the
  perturbation-vs-attainment dichotomy we measure.
* Dostert–de Laat–Moustrou (arXiv:2001.00256) and follow-up
  (arXiv:2403.16874) — the exact-rational-rounding pipeline for SDP
  bounds (our `verify_certificate.py` discipline is the same
  paradigm).
* de Laat (arXiv:1610.04905) — moment hierarchies for continuous
  energy on $S^2$: numerically sharp, never exact — consistent with
  the non-attainment this program works around.
