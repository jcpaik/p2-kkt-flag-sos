# Cylindrical mode domination

Track D: condition on the measure's own principal axis, expand the kernel in
azimuthal Fourier modes, and dominate the modes by explicit penalties in the
height marginal.  Everything here is independent of the flag-algebra
hierarchy.  Code: [`cylinder_cert.py`](../cylinder_cert.py); tests:
`test_cylinder_cert.py`; artifacts: `sdpa_runs/cylinder/`.

**Verdicts (2026-08-17).**

1. **M0.** The mode decomposition, parities, sharpness audit and every
   penalty ingredient are verified in exact arithmetic
   (`python cylinder_cert.py verify`, 16/16; independent slow re-derivation
   by symbolic integration: `derive`).
2. **M1 — new exactly solved sub-case.**  The conjecture restricted to
   *axisymmetric* measures is equivalent to copositivity of an explicit
   bicubic kernel `C(s,t)` on probability measures on `[0,1]`, and this is
   now **proved** by an exact rational PSD + Handelman certificate with the
   correct facial structure at the zero measure
   (`sdpa_runs/cylinder/m1_certificate.json`, re-verifiable in pure rational
   arithmetic by `python cylinder_cert.py m1 --verify PATH`; 6/6 checks).
3. **M2 — the designed master bound is false, and repairably so.**  The
   box-QP master bound `B(sigma)` fails badly at mid-latitude profiles
   (min ≈ −4.68).  Sharper: *no* per-mode scalar-domination scheme can work —
   the best conceivable one, `B_tight`, is exactly `−1/192 < 0` at the
   profile of the orthonormal 4-point cross and exactly `−64/729 < 0` at the
   profile of the ONB minimizer.  What scalar domination provably loses is
   the fiber moment coupling `|zeta_{2k}| >= 2|zeta_k|^2 − 1` between mode
   `k` and mode `2k`.
4. **M3 (partially executed) — second new exactly solved sub-case.**  Adding
   only the `(k, 2k)` Toeplitz couplings for `k = 2, 3` closes the
   single-profile case completely: **`E(mu) >= 0` for every antipodal
   measure supported on a pair of antipodal circles (any latitude, arbitrary
   fiber measures), with equality exactly on the ONB orbit** — proved in
   exact arithmetic (`verify_circle_pair_theorem`, 27/27 checks).  This
   strictly extends the great-circle sub-case of PLAN §1 and contains the
   ONB zero in its interior, sharply.
5. One claim from the track brief is **corrected**: the mode-2 principal-axis
   "gauge" `int (1-a^2) dmu_2 = 0` is *invalid* (only its imaginary part can
   be normalized away; the real part equals `lambda_2 − lambda_3 >= 0`, and
   the 4-point cross has `lambda_2 − lambda_3 = 1/2 != 0`).  The mode-1
   constraint `int a sqrt(1-a^2) dmu_1 = 0` *is* valid.  See §7.

---

## 1. Setup and exact mode decomposition

Fix a unit vector `e`.  For `x` on the sphere write `a = x·e in [−1,1]` and
`phi` for the azimuth around `e`.  For two points,

    x·y = a b + p cos(phi_x − phi_y),        p = sqrt((1−a^2)(1−b^2)).

With `K(t) = 32t^6 − 48t^4 + 20t^2 − 4/3 = cos 6θ + cos 2θ + 2/3`,

    K(x·y) = sum_{k=0}^{6} c_k(a,b) cos(k (phi_x − phi_y)).

**Exact mode kernels.**  With `s = a^2, t = b^2`, `eps = k mod 2`, each
`c_k(a,b) = sum_{ij} M_k[i,j] f^k_i(a) f^k_j(b)` on the basis
`f^k_i(a) = a^{eps+2i} (1−a^2)^{k/2}`:

| k | basis functions (per variable) | matrix `M_k` |
|---|---|---|
| 0 | `1, a^2, a^4, a^6` | `[[2/3,−4,12,−10],[−4,84,−270,210],[12,−270,840,−630],[−10,210,−630,462]]` |
| 1 | `a, a^3, a^5` times `sqrt(1−a^2)` | `[[16,−96,120],[−96,624,−720],[120,−720,792]]` |
| 2 | `1, a^2, a^4` times `(1−a^2)` | `[[1,−6,15],[−6,132,−270],[15,−270,495]]` |
| 3 | `a, a^3` times `(1−a^2)^{3/2}` | `[[12,−60],[−60,220]]` |
| 4 | `1, a^2` times `(1−a^2)^2` | `[[0,−6],[−6,66]]` |
| 5 | `a (1−a^2)^{5/2}` | `[[12]]` |
| 6 | `(1−a^2)^3` | `[[1]]` |

Notes.  `c_1` is *not* zero (PLAN §2 lists only `c_0, c_2..c_6`; the missing
`c_1 = 8ab·p·(99s^2t^2−90s^2t−90st^2+15s^2+15t^2+78st−12s−12t+2)` is supplied
here).  The determinant fingerprints of PLAN §2 are matched in the
integer-normalized bases: `det((3/2)M_0) = −34992`, `det M_2 = −6480`,
`det(M_3/4) = −60`, `det(M_4/6) = −1`.  Modes 5, 6 are rank-one *positive*;
`c_2, c_3, c_4` (and `c_0`) are indefinite.

**Verification** (exact, fast): the polynomial identity
`K(ab + pz) = sum_k c_k(a,b) T_k(z)` holds in
`Q[a,b,z][p]/(p^2 − (1−a^2)(1−b^2))` — `T_k(cos ψ) = cos kψ` plus uniqueness
of Fourier coefficients makes this equivalent to the mode decomposition.
An independent slow derivation by symbolic Fourier integrals is
`python cylinder_cert.py derive` (~2 min).

**Mode measures.**  For a finite positive measure `mu` let
`dmu_k(a) = int_fiber e^{ik phi} dmu` (a complex measure on `[−1,1]`) and
`sigma = mu_0` (the height marginal, `sigma >= 0`, mass 1 for probability
measures).  Then, with `Q_k := int int c_k(a,b) dmu_k(a) d(conj mu_k)(b)`
(real, no factor 2, because `mu_{−k} = conj(mu_k)` and `c_k` is symmetric):

    E(mu) = sum_{k=0}^{6} Q_k,    Q_5, Q_6 >= 0  (rank-one positive kernels).

**Parity.**  Antipodality `(a,phi) -> (−a, phi+pi)` gives
`dmu_k(−a) = (−1)^k dmu_k(a)`; `sigma` is even.  This matches the parity of
the basis functions, and for an even kernel `K` antipodality is free:
`E(mu-tilde) = E(mu)` for the antipodal symmetrization, so sub-cases proved
for all measures imply the antipodal statement and conversely.

**Principal-axis normalization.**  Take `e` = top eigenvector of
`Sigma = int x x^T dmu`.  Then

* `lambda_1 = int a^2 dsigma >= 1/3` (top eigenvalue of a trace-1 PSD matrix);
* `int a sqrt(1−a^2) dmu_1 = int x_3(x_1 + i x_2) dmu = 0` (eigenvector
  property) — the first mode-1 moment vanishes: `z_1 = 0`;
* after rotating the frame about `e`, `Im int (1−a^2) dmu_2 = 2 int x_1 x_2 = 0`
  and `Re int (1−a^2) dmu_2 = lambda_2 − lambda_3 in [0, 3 lambda_1 − 1]`.

Realizability: every even `sigma` with `int a^2 dsigma >= 1/3` occurs (take
the axisymmetric measure with marginal `sigma`), so the constraint set for
the master bound is exactly right.

## 2. Domination toolbox (each item is a theorem; ledger in §7)

**(D1) Total-variation domination.**  `|mu_k| <= sigma` as measures: for
Borel `A`, `|mu_k|(A) = sup_{|g|<=1} |int_A g dmu_k| <= mu(A × fiber)
= sigma(A)`.

**(D2) Moment box.**  `z_i = int f_i dmu_k` satisfies
`|z_i| <= int |f_i| d|mu_k| <= int |f_i| dsigma =: m_i(sigma)`.

**(D3) Signed-abs box-QP penalty.**  For real symmetric `M` and complex `z`
with `|z_i| <= m_i`:

    z^H M z >= min { u^T M~ u : 0 <= u_i <= m_i },
    M~_ii = M_ii,  M~_ij = −|M_ij| (i != j),

because `z^H M z >= sum_i M_ii |z_i|^2 − sum_{i!=j} |M_ij| |z_i||z_j|` at
`u_i = |z_i|`.  The box minimum of the homogeneous quadratic is computed
*exactly* by enumerating the `3^n` KKT activity patterns
(`box_qp_min`, exact rational arithmetic; singular free blocks are safely
skipped since face minima recur on their boundaries).

**(D4) Kernel-split domination** (stated for completeness; used implicitly).
If `c = P − D` with `P` a PSD kernel and `D >= 0` pointwise, then for
`|nu| <= sigma`: `int int c dnu d(conj nu) >= − int int D d|nu| d|nu|
>= − int int D dsigma dsigma`.

**(D5) Toeplitz mode coupling** (the repair; classical).  For a probability
measure `rho` on the circle with moments `zeta_k = int e^{ik phi} drho`:

    |zeta_{2k}| >= 2 |zeta_k|^2 − 1.

Proof: rotate so `zeta_k >= 0`; then `1 + Re zeta_{2k} = 2 int cos^2(k phi)
>= 2 (int cos k phi)^2 = 2 zeta_k^2` by Cauchy–Schwarz.  (Equivalently: the
3×3 Toeplitz Gram of `(1, e^{ik phi}, e^{2ik phi})` is PSD.)  Measure-level
form for the master program: the matrix measure
`[[sigma, mu_k, mu_{2k}], [conj mu_k, sigma, mu_k], [conj mu_{2k}, conj mu_k,
sigma]] >= 0`.

## 3. The master bound `B(sigma)` and its sharpness audit

For discrete even `sigma = sum_j w_j (delta_{alpha_j} + delta_{−alpha_j})/2`
define (implemented in `master_bound`, exact or float):

    B(sigma) = int int c_0 dsigma dsigma − sum_{k=1}^{4} pen_k(sigma),

where `pen_k` is the (D3) box-QP penalty for `M_k` with the dominated
moments of `|f^k_i|`, mode 1 restricted to the `(a^3, a^5)` block by the
axis constraint `z_1 = 0`, and modes 5, 6 dropped (PSD).  By construction
**`E(mu) >= B(sigma(mu))` is a theorem** for every antipodal probability
measure with principal axis `e`.

**Sharpness audit (exact).**  On the pole–equator family
`sigma_w = (1−w) delta_0 + w (delta_1 + delta_{−1})/2`: every dominated
moment that any penalty uses vanishes (`a^2(1−a^2)^2`, `a(1−a^2)^{3/2}`,
`a^{2i}(1−a^2)` for `i >= 1`, `a^{2i−1}sqrt(1−a^2)` all vanish on
`{0, ±1}`), the surviving mode-2 box has positive diagonal (`M_2[0,0] = 1`),
so `pen_k(sigma_w) = 0` for all k, and

    B(sigma_w) = 6 (w − 1/3)^2                 (exact; `family` CLI).

The mode ledger on the family reproduces PLAN §1's identity exactly:
equator mass `1−w` enters mode 2 through `f^2_0(0) = 1` and mode 6 through
`f^6_0(0) = 1`, giving `E = 6(w−1/3)^2 + (1−w)^2(|nu-hat(2)|^2 +
|nu-hat(6)|^2)`.

## 4. M1: the axisymmetric sub-case, solved exactly

**Reduction.**  For axisymmetric `mu` all `mu_k = 0` (k >= 1), so
`E = int int c_0 dsigma dsigma`; conversely any probability `sigma` on
`[−1,1]` occurs.  By evenness of `c_0`, WLOG `sigma` even; pushing forward to
`s = a^2` turns `c_0` into the bicubic

    C(s,t) = sum_{i,j<=3} M_0[i,j] s^i t^j     on [0,1]^2,

with anchors `C(0,0) = 2/3`, `C(1,t) = K(sqrt t)`, `C(1,1) = 8/3`.  The
axisymmetric case of the conjecture `<=>` `int int C dtau dtau >= 0` for
every probability `tau` on `[0,1]`.  (Not among PLAN §1's solved sub-cases;
numerically `min = 0` attained at `tau* = (2/3) delta_0 + (1/3) delta_1`,
i.e. `w = 1/3`.)

**Certificate.**  With `Phi = (1, s, s^2, s^3)`:

    C(s,t) = Phi(s)^T G Phi(t) + sum_beta lambda_beta B_beta(s,t),

* `G = P W P^T` PSD, `W` a rational 3×3 PSD matrix, `P` a rational basis of
  the orthogonal complement of the zero-measure moment vector
  `v* = (1, 1/3, 1/3, 1/3)` — so `G v* = 0` (facial sharpness);
* `B_beta` are symmetrized Handelman products
  `(s^i (1−s)^j t^k (1−t)^l + sym)/2` that vanish on the corner set
  `{0,1}^2` (facial sharpness at `supp tau* × supp tau*`), degree <= 8;
* all `lambda_beta >= 0`; only **8** products are needed, every one
  containing the factor `x(1−x)^2` whose `tau*`-average generates the
  first-order stationarity function `F(t) = (2/3)C(0,t) + (1/3)C(1,t)
  = 4t(1−t)^2` (exact).

Any such datum proves `int int C dtau dtau = int int Gram >= 0 (+) int int
H >= 0` for every probability `tau`.  Both parts vanish at `tau*`, matching
the zero.

**Status.**  Handelman degree 8 is the threshold for this formulation
(max-margin values: −0.855, −0.798, −0.383, −0.166, −0.107 at degrees
3..7; **+0.164** at 8); no multiplier was needed.  The exact rational
certificate (interior margin 0.067 before rounding; RREF rounding with
small denominators) is stored at `sdpa_runs/cylinder/m1_certificate.json`
and re-verified from the file in pure rational arithmetic:

    .venv/bin/python cylinder_cert.py m1 --verify sdpa_runs/cylinder/m1_certificate.json
    # identity, W PSD, G PSD, G v* = 0, lambda >= 0, corner-vanishing: 6/6

**Lemma (axisymmetric case).**  *For every probability measure `tau` on
`[0,1]`, `int int C dtau dtau >= 0`.  Hence `E(mu) >= 0` for every
axisymmetric measure `mu` on the sphere.*  (Machine-checked exact witness;
the reduction is the elementary argument above.)

## 5. M2: falsification of the master bound, with exact witnesses

**The designed bound fails.**  Structured scans and global optimization
(`m2` CLI; artifacts `m2_scan.json`, `m2_optimize.json`) put the minimum of
`B` at single mid-latitude pairs — with 2, 3, or 4 free atoms the optimizer
*collapses all atoms onto one latitude* `alpha ≈ 0.7151`
(`min B ≈ −4.6758`), so the worst case of the designed bound is exactly the
single-orbit family analyzed below.  Exact witnesses:

    sigma = pair(1/sqrt 2):  B = −128159/27456 ≈ −4.6678   (exact box-QPs)
       [c0 energy 25/96; pen_1..4 = 63/104, 257/64, 3/11, 3/88]
    sigma = pair(1/sqrt 3):  B = −106552/34749 ≈ −3.0663

The mode-2 box-QP dominates the loss: the box relaxation lets the three
dominated moments decouple, while realizable moment vectors from one
dominated measure are collinear.

**No per-mode scalar domination can work (exact impossibility).**  For the
single-pair profile `sigma_s = (delta_alpha + delta_{−alpha})/2`,
`s = alpha^2`, every antipodal measure with that height profile has

    E = C(s,s) + sum_{k=1}^{6} |rho-hat(k)|^2 v_k(s),
    v_k(s) = f^k(alpha)^T M_k f^k(alpha)   (exact polynomials; `rho` = upper fiber):

    v_1 = 8s(1−s) q_1,  q_1 = 99s^4−180s^3+108s^2−24s+2 > 0 on [0,1]
    v_2 = (1−s)^2 q_2,  q_2 = 495s^4−540s^3+162s^2−12s+1   (< 0 on ≈ (0.408, 0.622))
    v_3 = 4s(1−s)^3 g,  g = 55s^2−30s+3                    (< 0 on (3/11 ± 2 sqrt15/55))
    v_4 = 6s(1−s)^4(11s−2),                                (< 0 on (0, 2/11))
    v_5 = 12s(1−s)^5 >= 0,   v_6 = (1−s)^6 >= 0.

A scheme that knows only `|mu_k| <= sigma` (+ parity + the valid axis
constraints) must allow the *coherent* dominated measure
`mu_k = zeta · (fiber)` with `|zeta| = 1` for each single mode; hence its
bound is at most

    B_tight(s) = C(s,s) + sum_{k=2}^{4} min(0, v_k(s)),

and exactly:

    B_tight(1/2) = 25/96 − 17/64            = **−1/192**   (mode-2 leak)
    B_tight(1/3) = 64/243 − 256/729         = **−64/729**  (mode-3 leak)

Both witnesses are shadows of the two canonical zero configurations:

* `s = 1/2` is the profile of the **orthonormal 4-point cross**
  `{±u, ±v}, u·v = 0`, whose exact mode ledger is
  `E = 25/96 − 17/64 + 21/32 + 1/64 = 2/3`: scalar domination keeps the
  `−17/64` (mode 2) but cannot see the forced `+21/32` (mode 4) payback.
* `s = 1/3` is the profile of the **ONB minimizer** (±basis vectors seen
  along the cube diagonal: two 3-point fibers at `a = ±1/sqrt 3`), with the
  exact ledger `E = 64/243 − 256/729 + 64/729 = 0`: the dropped mode-6
  credit `+64/729` is precisely the undershoot.

In both cases the missing fact is the fiber moment coupling (D5):
`|zeta_2| = 1` forces `|zeta_4| = 1`, `|zeta_3| = 1` forces `|zeta_6| = 1`.
True fiber minima (numeric diagnostics): `min E(pair(1/sqrt2)) ≈ 0.1262 > 0`
and `min E(pair(1/sqrt3)) ≈ 0.0000` — the optimizer rediscovers the ONB zero.

**Repairs within the design freedom.**

* *h-gauge in mode 2*: **invalid as specified** — see §7 (the real part of
  `int (1−a^2) dmu_2` does not vanish; it equals `lambda_2 − lambda_3`,
  e.g. `= 1/2` for the cross).  The frame rotation only kills the imaginary
  part; the surviving real one-sided constraint
  `z_1^{(2)} = Re z_1^{(2)} in [0, 3 lambda_1 − 1]` is recorded as a usable
  cut, but it does *not* exclude the coherent leak (`zeta = 1` satisfies it).
* *Keeping PSD credits / tighter joint QPs*: already inside `B_tight`;
  insufficient by the exact witnesses above.
* *Toeplitz coupling (M3)*: sufficient on single-orbit profiles — next
  section — and exactly sharp at the ONB shadow.

## 6. M3: the coupling repair and the circle-pair theorem (exact)

Add to the scalar tools the single coupling (D5) for the pairs `(2,4)` and
`(3,6)`.  For the single-pair profile this yields the bound

    E >= C2(s) + min { v2 y2 + v4 y4 + v3 y3 + v6 y6 :
                       y in [0,1]^4,  y4 >= ((2 y2 − 1)^+)^2,
                       y6 >= ((2 y3 − 1)^+)^2 }        (y_k = |rho-hat(k)|^2,
                       v1, v5 credits dropped; v4 loaded worst-case where v4 < 0).

**Theorem (circle-pair case).**  *For every antipodal measure supported on
`{a = ±alpha}` (any `alpha`, arbitrary fiber measures), `E(mu) >= 0`;
equality holds exactly for `s = alpha^2 = 1/3` with uniform 3-point fibers —
the SO(3) orbit of the ONB configuration.*

Proof scheme, fully machine-checked in rational arithmetic
(`verify_circle_pair_theorem`, 27/27):

* the `(3,6)` coupling is *always* in its smooth branch:
  `4 v6 + v3 = 4(1−s)^3 (3s−1)^2 (6s+1) >= 0` — an exact identity, with the
  ONB latitude appearing as the double root `s = 1/3`;
* on the `(3,6)` window, the coupled minimum equals
  `C2 + v3/2 − s^2 g^2 = −(3s−1)^2 (891s^4−216s^3−81s^2−6s−2)/3 >= 0`
  because the quartic factor is negative there (PB) — *sharp at `s = 1/3`*;
* the `(2,4)` coupled minimum clears denominators to the polynomial facts
  (P2b), (PC), (PD); the `v4 < 0` region is (PA1)–(PA2); everything else is
  `C2 > 0` (P0) and `q_1 > 0` (P1);
* each polynomial fact is an exact Sturm root-count on an interval with
  rational endpoints, with rational bracketing certificates for the
  algebraic window endpoints.

Equality analysis: all facts are strict except (PB) at `s = 1/3`; there the
minimizing `y` forces `|zeta_3| = |zeta_6| = 1` and `zeta_1 = zeta_2 = 0`
(since `v_1, v_2 > 0` at `1/3`), i.e. uniform 3-point fibers: exactly ONB.

Quantified gain at the two critical profiles:

| profile | designed `B` | best scalar `B_tight` | + (k,2k) coupling | true fiber min |
|---|---|---|---|---|
| pair(1/sqrt2) | −4.6678 | −1/192 ≈ −0.0052 | **+1733/14336 ≈ +0.1209** | ≈ 0.1262 |
| pair(1/sqrt3) | −3.0663 | −64/729 ≈ −0.0878 | **0 (exact, sharp)** | 0 (ONB) |

**What remains for the full conjecture (the honest frontier).**  For general
`sigma` the exact reformulation is: `E = int int c_0 dsigma dsigma + sum_k
Q_k(zeta_k sigma)` where `zeta_k(a)` are the conditional fiber moments,
constrained pointwise by the 7×7 Toeplitz PSD condition (and parity, and the
axis normalizations).  Scalar domination keeps only `|zeta_k| <= 1` — proved
insufficient; the `(k,2k)` couplings close all single-atom profiles — proved
above; multi-atom profiles couple different latitudes through the indefinite
kernels `M_k`, and the convex tool for that level is the measure-valued
matrix domination of (D5):

    [[sigma, mu_2, mu_4], [*, sigma, mu_2], [*, *, sigma]] >= 0,
    [[sigma, mu_3, mu_6], [*, sigma, mu_3], [*, *, sigma]] >= 0,

i.e. matrix-kernel completions of the scalar penalties.  This is the
concrete next implementation step (not started here); the single-orbit
results pin down both the necessity and the expected sharpness structure
(pole–equator face for mode 0/2, ONB orbit for modes 3/6).

## 7. Validity ledger

Inequalities used, with status:

| item | statement | status |
|---|---|---|
| L1 | mode decomposition `K = sum c_k cos(k dphi)`, tables §1 | exact identity, machine-checked (fast Chebyshev route + slow integral route) |
| L2 | `E = sum Q_k`, no factor 2 | identity (`mu_{−k} = conj mu_k`, `c_k` symmetric) |
| L3 | `Q_5, Q_6 >= 0` | rank-one positive kernels (table §1) |
| L4 | `|mu_k| <= sigma`; moment boxes (D2) | measure theory, proof in §2 |
| L5 | box-QP penalty (D3) + exact KKT enumeration | proof in §2; exact arithmetic |
| L6 | parity of `mu_k` under antipodality | direct computation |
| L7 | `lambda_1 >= 1/3`; mode-1 axis constraint `z_1 = 0` | eigenvector property of the principal axis |
| L8 | mode-2 frame normalization: `Im int (1−a^2) dmu_2 = 0`, `Re ... = lambda_2 − lambda_3 in [0, 3 lambda_1 −1]` | valid one-sided version; the two-sided "gauge" of the track brief is **false** (cross counterexample) — do not use |
| L9 | Toeplitz coupling (D5) | Cauchy–Schwarz, proof in §2 |
| L10 | M1 certificate | exact rational witness, 6/6 independent checks |
| L11 | circle-pair theorem | exact, 27/27 checks (Sturm counts, rational endpoints) |
| — | `B(sigma) >= 0`? | **FALSE** — exact witnesses §5; superseded by the coupled route |

KKT-only assumptions used: none.  Every bound above holds for all
(antipodal) measures; the principal-axis constraints are normalizations, not
KKT conditions.

## 8. Artifacts and reproduction

All under `sdpa_runs/cylinder/` (gitignored, like the rest of `sdpa_runs/`):

* `m1_certificate.json` — exact rational M1 certificate (8 Handelman terms,
  denominators < 10^8).  Verify:
  `.venv/bin/python cylinder_cert.py m1 --verify sdpa_runs/cylinder/m1_certificate.json`.
  A committed copy lives at `certificates/m1_axisymmetric.json` (same
  file; verified 6/6 from that path too) so the theorem's witness
  survives any `sdpa_runs/` wipe.
* `m2_scan.json` — structured violation scan of the designed `B` (float).
* `m2_optimize.json` — multi-atom differential-evolution minima of `B`.
* `m2_verdict.json` — the exact numbers of §5–§6 in one place.

CLI entry points: `verify` (M0, 16 exact checks), `derive` (slow independent
re-derivation), `family` (sharpness audit), `m1 --solve/--verify`,
`m2 [--quick]`, `pair-theorem` (the 27 exact checks of §6).
Tests: `python -m pytest test_cylinder_cert.py -q`.
