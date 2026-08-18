# Exactly solved sub-cases, and the measurement record

Companion to [FOUNDATIONS.md](FOUNDATIONS.md) (definitions, validity
discipline, the reduction lemma), [SHARP_STRUCTURE.md](SHARP_STRUCTURE.md)
(the (E1)/(E2) sharp-face analysis), and [ENRICHMENTS.md](ENRICHMENTS.md)
(the cut/weight families: Jensen, fiber-Toeplitz, theta atoms, $e_5$).

* **Part A** collects the theorems that are *proved* — exact,
  independently re-verifiable sub-cases of the conjecture
  $E(\mu)\ge0$, $K(t)=32t^6-48t^4+20t^2-\tfrac43$ — together with the
  exact falsification results that delimit what scalar domination can
  do.  Code: `cylinder_cert.py`; tests `test_cylinder_cert.py`.
* **Part B** is the consolidated measurement record of the certificate
  program: baselines, escape fingerprints, selector trace laws, the
  dual-bound cascade, and the standing verdicts.  Data files under
  `sdpa_runs/`.

---

# PART A — Exactly solved sub-cases

Track: condition on the measure's own principal axis, expand the kernel
in azimuthal Fourier modes, and dominate the modes by explicit
penalties in the height marginal.  Everything here is independent of
the flag-algebra hierarchy.  The mode decomposition, parities,
sharpness audit and every penalty ingredient are verified in exact
arithmetic (`python cylinder_cert.py verify`, **16/16**; independent
slow re-derivation by symbolic integration: `derive`).

## A.1 Setup and exact mode decomposition

Fix a unit vector `e`.  For `x` on the sphere write `a = x·e in [−1,1]`
and `phi` for the azimuth around `e`.  For two points,

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

Notes.  `c_1` is *not* zero (older mode lists give only
`c_0, c_2..c_6`; the missing correction is
`c_1 = 8ab·p·(99s^2t^2−90s^2t−90st^2+15s^2+15t^2+78st−12s−12t+2)`).
Determinant fingerprints in the integer-normalized bases:
`det((3/2)M_0) = −34992`, `det M_2 = −6480`, `det(M_3/4) = −60`,
`det(M_4/6) = −1`.  Modes 5, 6 are rank-one *positive*; `c_2, c_3, c_4`
(and `c_0`) are indefinite.

**Verification** (exact, fast): the polynomial identity
`K(ab + pz) = sum_k c_k(a,b) T_k(z)` holds in
`Q[a,b,z][p]/(p^2 − (1−a^2)(1−b^2))` — `T_k(cos ψ) = cos kψ` plus
uniqueness of Fourier coefficients makes this equivalent to the mode
decomposition.  An independent slow derivation by symbolic Fourier
integrals is `python cylinder_cert.py derive`.

**Mode measures.**  For a finite positive measure `mu` let
`dmu_k(a) = int_fiber e^{ik phi} dmu` (a complex measure on `[−1,1]`)
and `sigma = mu_0` (the height marginal, `sigma >= 0`, mass 1 for
probability measures).  Then, with
`Q_k := int int c_k(a,b) dmu_k(a) d(conj mu_k)(b)` (real, no factor 2,
because `mu_{−k} = conj(mu_k)` and `c_k` is symmetric):

    E(mu) = sum_{k=0}^{6} Q_k,    Q_5, Q_6 >= 0  (rank-one positive kernels).

**Parity.**  Antipodality `(a,phi) -> (−a, phi+pi)` gives
`dmu_k(−a) = (−1)^k dmu_k(a)`; `sigma` is even.  This matches the
parity of the basis functions, and for an even kernel `K` antipodality
is free: `E(mu-tilde) = E(mu)` for the antipodal symmetrization, so
sub-cases proved for all measures imply the antipodal statement and
conversely.

**Principal-axis normalization.**  Take `e` = top eigenvector of
`Sigma = int x x^T dmu`.  Then

* `lambda_1 = int a^2 dsigma >= 1/3` (top eigenvalue of a trace-1 PSD
  matrix);
* `int a sqrt(1−a^2) dmu_1 = int x_3(x_1 + i x_2) dmu = 0` (eigenvector
  property) — the first mode-1 moment vanishes: `z_1 = 0`;
* after rotating the frame about `e`,
  `Im int (1−a^2) dmu_2 = 2 int x_1 x_2 = 0` and
  `Re int (1−a^2) dmu_2 = lambda_2 − lambda_3 in [0, 3 lambda_1 − 1]`.

Realizability: every even `sigma` with `int a^2 dsigma >= 1/3` occurs
(take the axisymmetric measure with marginal `sigma`), so the
constraint set for the master bound is exactly right.

## A.2 Domination toolbox (each item is a theorem; ledger in A.9)

**(D1) Total-variation domination.**  `|mu_k| <= sigma` as measures:
for Borel `A`, `|mu_k|(A) = sup_{|g|<=1} |int_A g dmu_k| <=
mu(A × fiber) = sigma(A)`.

**(D2) Moment box.**  `z_i = int f_i dmu_k` satisfies
`|z_i| <= int |f_i| d|mu_k| <= int |f_i| dsigma =: m_i(sigma)`.

**(D3) Signed-abs box-QP penalty.**  For real symmetric `M` and complex
`z` with `|z_i| <= m_i`:

    z^H M z >= min { u^T M~ u : 0 <= u_i <= m_i },
    M~_ii = M_ii,  M~_ij = −|M_ij| (i != j),

because `z^H M z >= sum_i M_ii |z_i|^2 − sum_{i!=j} |M_ij| |z_i||z_j|`
at `u_i = |z_i|`.  The box minimum of the homogeneous quadratic is
computed *exactly* by enumerating the `3^n` KKT activity patterns
(`box_qp_min`, exact rational arithmetic; singular free blocks are
safely skipped since face minima recur on their boundaries).

**(D4) Kernel-split domination.**  If `c = P − D` with `P` a PSD kernel
and `D >= 0` pointwise, then for `|nu| <= sigma`:
`int int c dnu d(conj nu) >= − int int D d|nu| d|nu| >= − int int D
dsigma dsigma`.

**(D5) Toeplitz mode coupling** (the repair; classical).  For a
probability measure `rho` on the circle with moments
`zeta_k = int e^{ik phi} drho`:

    |zeta_{2k}| >= 2 |zeta_k|^2 − 1.

Proof: rotate so `zeta_k >= 0`; then
`1 + Re zeta_{2k} = 2 int cos^2(k phi) >= 2 (int cos k phi)^2 = 2
zeta_k^2` by Cauchy–Schwarz.  (Equivalently: the 3×3 Toeplitz Gram of
`(1, e^{ik phi}, e^{2ik phi})` is PSD.)  Measure-level form for the
master program: the matrix measure
`[[sigma, mu_k, mu_{2k}], [conj mu_k, sigma, mu_k], [conj mu_{2k},
conj mu_k, sigma]] >= 0`.

## A.3 The master bound `B(sigma)` and its sharpness audit

For discrete even `sigma = sum_j w_j (delta_{alpha_j} +
delta_{−alpha_j})/2` define (implemented in `master_bound`, exact or
float):

    B(sigma) = int int c_0 dsigma dsigma − sum_{k=1}^{4} pen_k(sigma),

where `pen_k` is the (D3) box-QP penalty for `M_k` with the dominated
moments of `|f^k_i|`, mode 1 restricted to the `(a^3, a^5)` block by
the axis constraint `z_1 = 0`, and modes 5, 6 dropped (PSD).  By
construction **`E(mu) >= B(sigma(mu))` is a theorem** for every
antipodal probability measure with principal axis `e`.

**Sharpness audit (exact).**  On the pole–equator family
`sigma_w = (1−w) delta_0 + w (delta_1 + delta_{−1})/2`: every dominated
moment that any penalty uses vanishes (`a^2(1−a^2)^2`,
`a(1−a^2)^{3/2}`, `a^{2i}(1−a^2)` for `i >= 1`,
`a^{2i−1}sqrt(1−a^2)` all vanish on `{0, ±1}`), the surviving mode-2
box has positive diagonal (`M_2[0,0] = 1`), so `pen_k(sigma_w) = 0`
for all k, and

    B(sigma_w) = 6 (w − 1/3)^2                 (exact; `family` CLI).

The mode ledger on the family reproduces the known pole–equator
identity exactly: equator mass `1−w` enters mode 2 through
`f^2_0(0) = 1` and mode 6 through `f^6_0(0) = 1`, giving
`E = 6(w−1/3)^2 + (1−w)^2(|nu-hat(2)|^2 + |nu-hat(6)|^2)`.

## A.4 Theorem (axisymmetric case), solved exactly

**Reduction.**  For axisymmetric `mu` all `mu_k = 0` (k >= 1), so
`E = int int c_0 dsigma dsigma`; conversely any probability `sigma` on
`[−1,1]` occurs.  By evenness of `c_0`, WLOG `sigma` even; pushing
forward to `s = a^2` turns `c_0` into the bicubic

    C(s,t) = sum_{i,j<=3} M_0[i,j] s^i t^j     on [0,1]^2,

with anchors `C(0,0) = 2/3`, `C(1,t) = K(sqrt t)`, `C(1,1) = 8/3`.
The axisymmetric case of the conjecture `<=>`
`int int C dtau dtau >= 0` for every probability `tau` on `[0,1]`.
(Numerically `min = 0` attained at
`tau* = (2/3) delta_0 + (1/3) delta_1`, i.e. `w = 1/3`.)

**Certificate.**  With `Phi = (1, s, s^2, s^3)`:

    C(s,t) = Phi(s)^T G Phi(t) + sum_beta lambda_beta B_beta(s,t),

* `G = P W P^T` PSD, `W` a rational 3×3 PSD matrix, `P` a rational
  basis of the orthogonal complement of the zero-measure moment vector
  `v* = (1, 1/3, 1/3, 1/3)` — so `G v* = 0` (facial sharpness);
* `B_beta` are symmetrized Handelman products
  `(s^i (1−s)^j t^k (1−t)^l + sym)/2` that vanish on the corner set
  `{0,1}^2` (facial sharpness at `supp tau* × supp tau*`), degree <= 8;
* all `lambda_beta >= 0`; only **8** products are needed, every one
  containing the factor `x(1−x)^2` whose `tau*`-average generates the
  first-order stationarity function
  `F(t) = (2/3)C(0,t) + (1/3)C(1,t) = 4t(1−t)^2` (exact).

Any such datum proves `int int C dtau dtau = int int Gram >= 0 (+)
int int H >= 0` for every probability `tau`.  Both parts vanish at
`tau*`, matching the zero.

**Degree-8 feasibility threshold.**  Handelman degree 8 is the
threshold for this formulation (max-margin values at degrees 3..8):

| Handelman degree | 3 | 4 | 5 | 6 | 7 | **8** |
|---|---:|---:|---:|---:|---:|---:|
| max margin | −0.855 | −0.798 | −0.383 | −0.166 | −0.107 | **+0.164** |

No multiplier was needed.  The exact rational certificate (interior
margin 0.067 before rounding; RREF rounding with small denominators)
is stored at `sdpa_runs/cylinder/m1_certificate.json`; a **committed
copy lives at `certificates/m1_axisymmetric.json`** (same file,
verified 6/6 from that path too), so the theorem's witness survives
any `sdpa_runs/` wipe.  Independent verification, pure rational
arithmetic:

    .venv/bin/python cylinder_cert.py m1 --verify sdpa_runs/cylinder/m1_certificate.json
    # identity, W PSD, G PSD, G v* = 0, lambda >= 0, corner-vanishing: 6/6

**Lemma (axisymmetric case).**  *For every probability measure `tau`
on `[0,1]`, `int int C dtau dtau >= 0`.  Hence `E(mu) >= 0` for every
axisymmetric measure `mu` on the sphere.*  (Machine-checked exact
witness; the reduction is the elementary argument above.)

## A.5 Exact falsification of scalar mode-domination

**The designed bound fails.**  Structured scans and global optimization
(`m2` CLI; artifacts `m2_scan.json`, `m2_optimize.json`) put the
minimum of `B` at single mid-latitude pairs — with 2, 3, or 4 free
atoms the optimizer *collapses all atoms onto one latitude*
`alpha ≈ 0.7151` (`min B ≈ −4.6758`), so the worst case of the
designed bound is exactly the single-orbit family analyzed below.
Exact witnesses:

    sigma = pair(1/sqrt 2):  B = −128159/27456 ≈ −4.6678   (exact box-QPs)
       [c0 energy 25/96; pen_1..4 = 63/104, 257/64, 3/11, 3/88]
    sigma = pair(1/sqrt 3):  B = −106552/34749 ≈ −3.0663

The mode-2 box-QP dominates the loss: the box relaxation lets the
three dominated moments decouple, while realizable moment vectors from
one dominated measure are collinear.

**No per-mode scalar domination can work (exact impossibility).**  For
the single-pair profile `sigma_s = (delta_alpha + delta_{−alpha})/2`,
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
`mu_k = zeta · (fiber)` with `|zeta| = 1` for each single mode; hence
its bound is at most

    B_tight(s) = C(s,s) + sum_{k=2}^{4} min(0, v_k(s)),

and exactly:

    B_tight(1/2) = 25/96 − 17/64            = **−1/192**   (mode-2 leak)
    B_tight(1/3) = 64/243 − 256/729         = **−64/729**  (mode-3 leak)

Both witnesses are shadows of the two canonical zero configurations:

* `s = 1/2` is the profile of the **orthonormal 4-point cross**
  `{±u, ±v}, u·v = 0`, whose exact mode ledger is
  `E = 25/96 − 17/64 + 21/32 + 1/64 = 2/3`: scalar domination keeps
  the `−17/64` (mode 2) but cannot see the forced `+21/32` (mode 4)
  payback.
* `s = 1/3` is the profile of the **ONB minimizer** (±basis vectors
  seen along the cube diagonal: two 3-point fibers at
  `a = ±1/sqrt 3`), with the exact ledger
  `E = 64/243 − 256/729 + 64/729 = 0`: the dropped mode-6 credit
  `+64/729` is precisely the undershoot.

In both cases the missing fact is the fiber moment coupling (D5):
`|zeta_2| = 1` forces `|zeta_4| = 1`, `|zeta_3| = 1` forces
`|zeta_6| = 1`.  True fiber minima (numeric diagnostics):
`min E(pair(1/sqrt2)) ≈ 0.1262 > 0` and
`min E(pair(1/sqrt3)) ≈ 0.0000` — the optimizer rediscovers the ONB
zero.

**Repairs within the design freedom.**

* *h-gauge in mode 2*: **invalid as specified** — see A.7 (the real
  part of `int (1−a^2) dmu_2` does not vanish; it equals
  `lambda_2 − lambda_3`, e.g. `= 1/2` for the cross).  The frame
  rotation only kills the imaginary part; the surviving real one-sided
  constraint `z_1^{(2)} = Re z_1^{(2)} in [0, 3 lambda_1 − 1]` is
  recorded as a usable cut, but it does *not* exclude the coherent
  leak (`zeta = 1` satisfies it).
* *Keeping PSD credits / tighter joint QPs*: already inside `B_tight`;
  insufficient by the exact witnesses above.
* *Toeplitz coupling*: sufficient on single-orbit profiles — next
  section — and exactly sharp at the ONB shadow.

## A.6 Theorem (circle-pair case): the coupling repair, exact

Add to the scalar tools the single coupling (D5) for the pairs `(2,4)`
and `(3,6)`.  For the single-pair profile this yields the bound

    E >= C2(s) + min { v2 y2 + v4 y4 + v3 y3 + v6 y6 :
                       y in [0,1]^4,  y4 >= ((2 y2 − 1)^+)^2,
                       y6 >= ((2 y3 − 1)^+)^2 }        (y_k = |rho-hat(k)|^2,
                       v1, v5 credits dropped; v4 loaded worst-case where v4 < 0).

**Theorem (circle-pair case).**  *For every antipodal measure
supported on `{a = ±alpha}` (any `alpha`, arbitrary fiber measures),
`E(mu) >= 0`; equality holds exactly for `s = alpha^2 = 1/3` with
uniform 3-point fibers — the SO(3) orbit of the ONB configuration.*

This strictly extends the great-circle sub-case and contains the ONB
zero in its interior, sharply.  Proof scheme, fully machine-checked in
rational arithmetic (`verify_circle_pair_theorem`, **27/27**):

* the `(3,6)` coupling is *always* in its smooth branch:
  `4 v6 + v3 = 4(1−s)^3 (3s−1)^2 (6s+1) >= 0` — an exact identity,
  with the ONB latitude appearing as the double root `s = 1/3`;
* on the `(3,6)` window, the coupled minimum equals
  `C2 + v3/2 − s^2 g^2 = −(3s−1)^2 (891s^4−216s^3−81s^2−6s−2)/3 >= 0`
  because the quartic factor is negative there (PB) — *sharp at
  `s = 1/3`*;
* the `(2,4)` coupled minimum clears denominators to the polynomial
  facts (P2b), (PC), (PD); the `v4 < 0` region is (PA1)–(PA2);
  everything else is `C2 > 0` (P0) and `q_1 > 0` (P1);
* each of the seven polynomial facts is an exact Sturm root-count on
  an interval with rational endpoints, with rational bracketing
  certificates for the algebraic window endpoints.

Equality analysis: all facts are strict except (PB) at `s = 1/3`;
there the minimizing `y` forces `|zeta_3| = |zeta_6| = 1` and
`zeta_1 = zeta_2 = 0` (since `v_1, v_2 > 0` at `1/3`), i.e. uniform
3-point fibers: exactly ONB.

Quantified gain at the two critical profiles:

| profile | designed `B` | best scalar `B_tight` | + (k,2k) coupling | true fiber min |
|---|---|---|---|---|
| pair(1/sqrt2) | −4.6678 | −1/192 ≈ −0.0052 | **+1733/14336 ≈ +0.1209** | ≈ 0.1262 |
| pair(1/sqrt3) | −3.0663 | −64/729 ≈ −0.0878 | **0 (exact, sharp)** | 0 (ONB) |

## A.7 The invalid gauge, corrected

One claim from the original track brief is **corrected**: the mode-2
principal-axis "gauge" `int (1−a^2) dmu_2 = 0` is *invalid* — only its
imaginary part can be normalized away.  Exactly:

    Im int (1−a^2) dmu_2 = 0   (after frame rotation),
    Re int (1−a^2) dmu_2 = lambda_2 − lambda_3 in [0, 3 lambda_1 − 1],

and the 4-point cross has `lambda_2 − lambda_3 = 1/2 != 0`.  The
two-sided gauge must not be used; the one-sided real constraint is a
valid cut.  The mode-1 constraint `int a sqrt(1−a^2) dmu_1 = 0` *is*
valid (eigenvector property).

## A.8 The honest frontier (what remains for the full conjecture)

For general `sigma` the exact reformulation is:
`E = int int c_0 dsigma dsigma + sum_k Q_k(zeta_k sigma)` where
`zeta_k(a)` are the conditional fiber moments, constrained pointwise
by the 7×7 Toeplitz PSD condition (and parity, and the axis
normalizations).  Scalar domination keeps only `|zeta_k| <= 1` —
proved insufficient (A.5); the `(k,2k)` couplings close all
single-atom profiles — proved above (A.6); multi-atom profiles couple
different latitudes through the indefinite kernels `M_k`, and the
convex tool for that level is the measure-valued matrix domination of
(D5):

    [[sigma, mu_2, mu_4], [*, sigma, mu_2], [*, *, sigma]] >= 0,
    [[sigma, mu_3, mu_6], [*, sigma, mu_3], [*, *, sigma]] >= 0,

i.e. matrix-kernel completions of the scalar penalties.  The
single-orbit results pin down both the necessity and the expected
sharpness structure (pole–equator face for mode 0/2, ONB orbit for
modes 3/6).  The flag-algebra incarnation of this frontier is the
Jensen/fiber-Toeplitz family program of
[ENRICHMENTS.md](ENRICHMENTS.md), measured in Part B.

## A.9 Validity ledger

| item | statement | status |
|---|---|---|
| L1 | mode decomposition `K = sum c_k cos(k dphi)`, tables A.1 | exact identity, machine-checked (fast Chebyshev route + slow integral route) |
| L2 | `E = sum Q_k`, no factor 2 | identity (`mu_{−k} = conj mu_k`, `c_k` symmetric) |
| L3 | `Q_5, Q_6 >= 0` | rank-one positive kernels (table A.1) |
| L4 | `|mu_k| <= sigma`; moment boxes (D2) | measure theory, proof in A.2 |
| L5 | box-QP penalty (D3) + exact KKT enumeration | proof in A.2; exact arithmetic |
| L6 | parity of `mu_k` under antipodality | direct computation |
| L7 | `lambda_1 >= 1/3`; mode-1 axis constraint `z_1 = 0` | eigenvector property of the principal axis |
| L8 | mode-2 frame normalization: `Im int (1−a^2) dmu_2 = 0`, `Re ... = lambda_2 − lambda_3 in [0, 3 lambda_1 −1]` | valid one-sided version; the two-sided "gauge" is **false** (cross counterexample) — do not use |
| L9 | Toeplitz coupling (D5) | Cauchy–Schwarz, proof in A.2 |
| L10 | axisymmetric certificate | exact rational witness, 6/6 independent checks |
| L11 | circle-pair theorem | exact, 27/27 checks (Sturm counts, rational endpoints) |
| — | `B(sigma) >= 0`? | **FALSE** — exact witnesses A.5; superseded by the coupled route |

KKT-only assumptions used: none.  Every bound above holds for all
(antipodal) measures; the principal-axis constraints are
normalizations, not KKT conditions.

## A.10 Artifacts and reproduction

Under `sdpa_runs/cylinder/` (gitignored): `m1_certificate.json` (exact
rational axisymmetric certificate, 8 Handelman terms, denominators
< 10^8), `m2_scan.json`, `m2_optimize.json`, `m2_verdict.json`.
Committed witness: `certificates/m1_axisymmetric.json`.

CLI entry points of `cylinder_cert.py`: `verify` (16 exact checks),
`derive` (slow independent re-derivation), `family` (sharpness audit),
`m1 --solve/--verify`, `m2 [--quick]`, `pair-theorem` (the 27 exact
checks of A.6).  Tests: `python -m pytest test_cylinder_cert.py -q`.

---

# PART B — The measurement record

All numbers below are in the **$E$ normalization** (target $E$, export
shift $-4/3$; target $h_2E$, shift $+2/3$) unless marked legacy.
Selector convention: minimal-trace certificates for
$h_2E \ge -\varepsilon$, i.e. SDPA bound line
$t_0 = -(2/3+\varepsilon)$ (e.g. `sel_toep_1em3`:
`--bound=-6.6766…E-1`; `sel_toep_1em4`: `-6.6676…E-1`; the
$5.33\times10^{-5}$ point uses the exact `-6.6672E-1`, the
$2\times10^{-5}$ point `-6.66686…E-1`).  The reported bound of any
export = objValPrimal + objective_shift.  Trace-law exponents are
local: $\operatorname{tr}\propto\varepsilon^{-\gamma}$ on the quoted
interval.

## B.1 Regenerated baseline

* Degree-14 unweighted pruned nine-toggle base (m = 398):
  $-4.48560259\times10^{-3}$ (200-bit, feasibility errors
  $\sim10^{-26}$).
* Degree-14 KKT-inclusive weighted problem `deg14_h2w_h2all.dat-s`
  (m = 790): anchor bound $\mathbf{-2.7908758\times10^{-5}}$,
  re-verified after the full regeneration.
* Stage table of the $h_2$-weighted reformulation (degree 14, same
  label algebra throughout; the $161\times$ gain):

| # | configuration | m | bound |
|---|---|---:|---:|
| 0 | unweighted $E$ (baseline) | 398 | $-4.4856\times10^{-3}$ |
| 1 | $h_2E$, same blocks | 398 | $-2.0840\times10^{-3}$ |
| 2 | $h_2E$ + `--h2-localized-flags` | 414 | $-1.2032\times10^{-3}$ |
| 3 | $h_2E$ + $p_2\times$(all families) — negative control | 790 | $-2.0840\times10^{-3}$ (identical to run 1 to eleven digits) |
| 4 | $h_2E$ + $h_2\times$(all families) (`--h2-localized-all`) | 790 | $\mathbf{-2.7909\times10^{-5}}$ |

The run-3 control proves the gain is carried entirely by the
$-\tfrac12$ base-label coupling of the $h_2$-localizer (a
$p_2$-multiplied block touches only fresh product labels and is
vacuous).  Legacy attainment point: minimal trace for legacy-scaled
$h_2E\ge-10^{-4}$ was 330.1 vs the unweighted pole law
$1.07/\varepsilon$ prediction $\sim10{,}700$ — a $32\times$ pole
collapse (ratios scale-free).

## B.2 Escape fingerprints

**Control trace law** (weighted deg-14, selectors
`sel_h2all_{1em3,1em4,5em5}`): growth $10.63\times$/decade, local
exponent $\varepsilon^{-1.03}$ — a clean simple pole (the historical
$16.7\times$ was wall-inflated), confirming a *single* multiplier
deficit.

**Block-level carriers** (Frobenius growth
$\varepsilon=10^{-3}\to10^{-4}$):

| block | frob@$10^{-3}$ | frob@$10^{-4}$ | ratio | trace@$10^{-4}$ |
|---|---:|---:|---:|---:|
| `h2loc_empty_type_flag` | $5.1\times10^{-8}$ | $517.9$ | $10^{10}$ | $517.9$ |
| `h2loc_two_root_even_11` | $23.8$ | $1438.2$ | $60.5$ | $2849.6$ |
| `h2loc_two_root_odd_10/01` | $21.6$ | $391.0$ | $18.1$ | $745.9$ each |
| `h2loc_two_root_odd_*_minor` | $3.6$ | $126.4$ | $35.1$ | $186.0$ each |
| `h2loc_two_root_even_00` | $261.4$ | $1071.3$ | $4.1$ | $1769.5$ |
| `h2loc_two_root_even_11_minor` | $42.3$ | $214.6$ | $5.1$ | $407.1$ |

All other blocks static or negligible: the escape is the
**$h_2$-localized collision/theta tower**.

**Sector masses of the growth direction**
$D=e(Y_{10^{-4}})-e(Y_{10^{-3}})$ ($\|D\|_1=5.0\times10^4$;
$\langle\text{target},D\rangle=-3.64\times10^4<0$):

| sector | $10^{-3}\!\to\!10^{-4}$ | $10^{-4}\!\to\!5.33\times10^{-5}$ |
|---|---:|---:|
| $p_2\times$`graph_4` products | $57.7\%$ | $64.5\%$ |
| `graph_4` connected | $19.2\%$ | $21.4\%$ |
| $p_2\times$pair-products ($p_2p_4p_6$-type) | $10.8\%$ | $2.5\%$ |
| $p_2\times$`triangle` products | $9.1\%$ | $9.1\%$ |
| `triangle` connected | $3.2\%$ | — |
| pair sector ($p_k$) | $0.07\%$ | $0.0\%$ |

Top labels: $p_2\times$(`graph_4`, edge multiplicities up to 5) —
top label $p_2\times(\texttt{graph\_4},0,1,3,1,3,2)$ in both
intervals — $p_2\times$`triangle` $(1,3,5)$, $(0,4,6)$, $(1,5,5)$,
$(1,3,3)$, $p_2\times(p_2p_4p_6,\,p_2p_4p_4,\,p_2p_2p_4)$, both signs.
Stable across $\varepsilon$: not a wall artifact.  Consequence (sign
rule): every scalar invariant expanding in pair labels alone (powers
of $h_2$, $e_k(I-A_2)$ scalar cuts, $g_\ell$ combinations) pairs
$\approx0$ with $D$ — dead as weights *and* cuts.

**v2 residual rotation** (`fingerprint_D_toep2.json`,
$\|D\|_1=1.33\times10^4$, down from $5.0\times10^4$):
$p_2\times$graph$_4$ + graph$_4$ collapsed $77\%\to3.8\%$; residual
$49.7\%$ $p_2\times$triangle, $29.3\%$ $p_2\times$pair-products,
$17.1\%$ triangle; top labels $p_2\times$triangle$(1,3,5)$, $(0,4,6)$,
$(1,5,5)$; carriers `h2loc_two_root_even_00`, `h2loc_flag_2/3`.  The
certificate pumps the same-leaf (T) side of the Jensen inequality —
the matched v3 object is the fiber-Toeplitz moment matrix.

**Deep saturation return** (deg-14, below $\sim5\times10^{-5}$): the
law re-steepens (local exponent 1.17 on
$5.33\times10^{-5}\to2\times10^{-5}$, at $5.5\times$ the toep wall)
and the deep escape returns through `h2loc_two_root_even_11`
($10.9\to753$) — the sector the Jensen blocks suppress in the window.

**Weighted-(E1) admissible split** (selector growth per h2loc block,
`fingerprint_admissible_split_e1.json`):

| interval | $\|\Delta Y\|_F$ | admissible | dead | dead share |
|---|---:|---:|---:|---:|
| $10^{-3}\to10^{-4}$ | 1864.4 | 1829.2 | 360.8 | 19.4% |
| $10^{-4}\to5.3\times10^{-5}$ | 3974.7 | 3971.0 | 170.5 | **4.3%** |

The pole is carried by content *inside* the sharp admissible cone: no
facial reduction can remove it.

## B.3 Selector trace laws (all measured configurations)

Weighted target $h_2E$, KKT-inclusive cone unless noted; GMP
(128-bit selectors, 200-bit bounds).

| configuration | $10^{-3}$ | $10^{-4}$ | $5.33\times10^{-5}$ | $2\times10^{-5}$ | law |
|---|---:|---:|---:|---:|---|
| deg-14 control | 764.99 | 8133.57 | 13637.2 | — | $10.63\times$/dec, exp 1.03 |
| deg-14 + theta windows/caps | 765.0 | 8133.57 | — | — | **inert** (identical to 6 digits) |
| deg-14 + e5 cut | — | 8133.566 | (unconverged) | (unconverged) | **inert at GMP**; double-precision $-1.832\times10^{-4}$ was noise |
| deg-14 + Jensen v1 (11 families) | 733.47 | 7034.96 | — | — | $9.59\times$/dec, exp 0.98 |
| deg-14 + v2 (29 families) | **437.51** | **1544.18** | **1909.11** | — | $3.53\times$/dec, exp 0.55 ($\to$0.34 below $10^{-4}$) |
| deg-14 + v3 (39 families) | **359.83** | **1472.80** | **1883.73** | **5934.81** | re-steepens deep: exp 1.17 on the last interval (saturation) |
| all-measures cone (proof-carrying) | 2366.76 | 13873.30 | 48232.6 | — | $5.86\times$/dec (exp 0.77) |

Degree 16 (grid $\varepsilon=10^{-4}/10^{-5}/5\times10^{-6}$):

| configuration | $10^{-4}$ | $10^{-5}$ | $5\times10^{-6}$ | law |
|---|---:|---:|---:|---|
| deg-16 plain | 4073.66 | 14068.66 | — | $3.45\times$/dec, exp 0.54 |
| deg-16 + Jensen v1 | 3752.58 | 8151.85 | — | $2.17\times$/dec, exp 0.34, **no deep re-steepening** |
| deg-16 + v3 (43 families, K=4) | **959.52** | **1546.05** | **1805.27** | $1.61\times$/dec, exp **0.207** (sub-decade 0.224), no re-steepening down to $\approx3.8\times$ the wall |

**Exponent ladder** (deepest clean interval per configuration):
control 1.03 → deg-16 plain 0.54 → +Jensen v1 0.34 → **+v2/v3 0.21**.
Each (degree × family-design) iteration roughly halves the residual
singularity; the certificate trace stays below 2000 approaching the
wall.

Deep e5-cut selector points ($5.33\times10^{-5}$, $2\times10^{-5}$) at
128-bit reached only pFEAS — unconverged, not counted as measurements.
Related precision note: the weighted deg-16 *bound* solve at 128-bit
produced a spurious pdINF; weighted bound solves at degree $\ge16$
require 200-bit.

## B.4 Dual-bound cascade

All 200-bit GMP, primal = dual to $\ge13$ digits unless noted:

| stage | bound | factor |
|---|---:|---|
| deg-14 unweighted | $-4.4856\times10^{-3}$ | 1 |
| deg-14 weighted ($h_2E$, `--h2-localized-all`) | $-2.7908758\times10^{-5}$ | $161\times$ |
| deg-14 + v2 Jensen (toep2 = toep3, pdOPT) | $-3.66051\times10^{-6}$ | $7.6\times$ further |
| deg-16 weighted | $-1.31\times10^{-6}$ | |
| **deg-16 + v3** | $\mathbf{-5.2504467\times10^{-9}}$ | $250\times$ over plain deg-16; $\sim850{,}000\times$ over the program start |

The deg-16+v3 wall is genuine, not a convergence artifact:
primal/dual agree to $1.8\times10^{-14}$ relative, feasibility errors
$5\times10^{-26}$ / $1.5\times10^{-34}$.  Active at the wall: the
Jensen d7/d8 towers (plain + h2loc) and `h2loc_two_root_even_00`.

Side branches (closed):

* deg-16 am `le1` (h2loc-layer-projected all-measures cone) dual:
  $-4.5762\times10^{-3}$ — the projection route is closed as a bound
  route (projection weakens the moment side; the deg-14 pure-layer
  weighted-(E1) projection is dual-unbounded at *every*
  $\varepsilon$).
* deg-14 e5-cut dual: $-2.7908758\times10^{-5}$, identical to
  baseline — bare cut inert at GMP precision.

## B.5 Wall adversary (deg-16+v3)

Pseudo-moments of the optimal dual at the wall (`sdpa_extract`):

    p2 = 1/3 + 9.7e−5   (h2 ≈ 1.5e−4)
    p4 = 0.29043        (between ONB 1/3 and pole+Haar 5/18 — partial equator mode-4 mass)
    p6 = 0.26892

— a perturbed zero-family configuration on which the cone still
tolerates $E \approx -3.6\times10^{-5}$.

## B.6 Cut-family verdicts in detail

### Theta atoms (windows, caps, localized package)

Ray pairings (regenerated (E1)-projected deg-14 escape ray; squared
norm 1.9586, $\langle E,r\rangle=-1.0000$, $g_4\cdot r = 1.0031$):
mode-6 diagonals $\ell_{-1},\ell_0,\ell_1,\ell_2 =
0.0201/0.0570/0.0200/0.0318$ (mode-2 at noise floor
$\sim g_2\cdot r = 4.4\times10^{-4}$); window-cut pairings all
kill-signed, $-0.010\dots-0.057$ across $q\in\{1/4,1/2,3/4,9/10\}$.
Yet A/B (q = 1/2):

| experiment | outcome |
|---|---|
| plain diagonal atoms (f6/f2/both), ray | ray survives at $+1.9\%$ norm (1.9586 → 1.9963), rotated escape: pure $g_4$, theta diagonals forced to zero |
| localized atoms (146 rows), all-measures projected cone, ray | **infeasible — every recession direction killed** (implies the same for the KKT cone by inclusion) |
| $\varepsilon=0$ bound, projected cone + gap cuts | $-3.18615$ (`optimal`, MOSEK) |
| $\varepsilon=0$ bound, + atoms on top of gap cuts | $-3.18615$ (identical to 6 digits; $\Delta\approx2\times10^{-6}$ noise) |
| selector trace at $\varepsilon=0.3$, control vs atoms | 6.4977 vs 6.4977 (identical to 4 decimals; 3 multipliers active, $\lambda_{\max}=1.60$, zero trace effect) |
| GMP pole-decade selector, `deg14_h2w_h2all_theta` | 765.0 / 8133.57 — **inert** |

Against the weighted escape $D$: 13 of 28 per-$n$ sup-caps pair
kill-signed (led by $\operatorname{cap}[\text{plain},6,n{=}2]:
-5.56\times10^6$, $\operatorname{cap}[h_2\text{loc},6,n{=}0]:
-5.53\times10^6$), but the GMP selector was unchanged: scalar
codimension-few slices reroute.  Tail majorants (exact rationals): at
$q=1/2$, $T^{(2)}_q(N)=238,\ 14.4,\ 0.733,\ 9.28\times10^{-3},\dots$
and $T^{(6)}_q(N)=1860,\ 108,\ 5.76,\ 7.52\times10^{-2},\dots$ for
$N=0..3$ — super-exponential collapse, so inertness is not a
constant-size artifact.

### Jensen / fiber-Toeplitz families (the movers)

Structural facts (machine-checked; `toeplitz_blocks.py --self-test`
6/6, dump cross-check exact): every conditional Toeplitz Gram over
modulated leaves with *independent* leaves is already spanned by the
`two_root_*` blocks; the new content is the unfolding inequality
$T-G=\mathbb E[\rho\,\mathrm{Cov}_y(\varphi)]\succeq0$ per sector,
its $h_2$-localized and $(1-h_2)$-complement copies, the pair-sample
copy $[p_{a+b}]-[p_ap_b]\succeq0$, and (v3) the fiber-Toeplitz moment
matrices — PSD by Toeplitz positivity, provably outside every
polynomial Gram span for $K\ge2$.  Escape seesaw: h2loc/plain pairing
ratio $\approx-4.9$ uniformly over 14 generators (a genuine measure
near ONB has $+p_2\approx+\tfrac13$) — the complement blocks forbid
it with zero labels outside the problem.

T2 eigen gate (necessary-condition filter, $M(D)=\sum_LD_L(T_L-G_L)$):
v1/v2 families 8–16 negative eigenvalues each, negative mass to
$-9.7\times10^4$ on $\|D\|_1=5\times10^4$, stable across both decades;
v3 fiber-Toeplitz families dwarf them against the rotated residual —
`h2loc_ftoep2_even_11_r3` neg sum $-9.8\times10^5$ on
$\|D\|_1=1.33\times10^4$ (up to $74\times$ the residual mass).

Double-precision selector A/B (MOSEK, all solves `optimal`):

| $\varepsilon$ | control | v1 | v2 | v3 |
|---:|---:|---:|---:|---:|
| $10^{-1}$ | 5.6038 | 5.4920 ($-2.0\%$) | 5.4408 ($-2.9\%$) | **5.1176** ($-8.7\%$) |
| $3\times10^{-2}$ | 7.8830 | 7.7691 ($-1.4\%$) | 7.7393 ($-1.8\%$) | **7.4630** ($-5.3\%$) |
| $10^{-2}$ | 12.6538 | 12.4915 ($-1.3\%$) | 12.2363 ($-3.3\%$) | **11.6718** ($-7.8\%$) |

First families measured *active* at a selector optimum
(`jensen_even_11` trace 0.15–0.19; `ftoep2_even_11_r3` top-6 block,
trace 0.31/0.20) — in contrast to gap cuts and theta atoms (slack
everywhere).  Deg-16 v3 A/B at $\varepsilon=10^{-2}$: control 11.8252
→ **10.7443** ($-9.1\%$).  GMP verdicts: B.3/B.4.

Export bookkeeping (first-class, exact): naive file-appends are
unsound (the base export drops 337 image-dependent directions);
`toeplitz_export.py` re-runs `sos_search.export_sdpa_problem` with the
families included.  m stays 794 (deg-14 v1 = v2 = v3) and 1245
(deg-16 v1 = v3); `dropped_objective_inconsistent = 0` on all.

Deg-18 sharp-face stack (`deg18_we1_toep3.dat-s`, base m = 1015, both
layers projected): at $\varepsilon=10^{-2}$ the pure we1 selector
fails (MOSEK error, the predicted $\varepsilon$-pathology of the
unstacked sharp-face cone) while the v3-stacked selector is optimal,
trace 20.4883.  GMP selectors pending.

### The $e_5(I-A_2)$ weight/cut

Exact expansion `gap_elementary_vector(5)`: 35 labels (11 `graph_5`),
identically 0 on the pole–equator stratum (symbolic), $(4/5)^5$ at
uniform, first vanishing invariant ($e_4$ does not vanish).
Double-precision measurements (all-measures cone, degree 14, baseline
$-1.9198\times10^{-4}$ reproduced):

| configuration | bound (double) |
|---|---:|
| + e5 cut module, target $h_2E$ | $-1.83199\times10^{-4}$ ($-4.6\%$; **later shown to be double-precision noise — GMP inert**) |
| `--e5-weight 1/4` (1×1 coverage only) | unbounded |
| `--e5-weight 1/4` (Hankel + AM–GM module) | $-1.7649\times10^{-2}$ |
| `--e5-weight 1` | $-7.0237\times10^{-2}$ |
| + arity-5 weighted flags, $\kappa=1/4$ | $-1.37225\times10^{-2}$ ($-22\%$: the `graph_5` sector works, but degree 14 cannot tame $\kappa e_5E$) |

Escape pairings of the visible $\le4$-sample shadow:
$\langle e_5^{\le4},D_{e3e4}\rangle=-6.37\times10^5$,
$\langle\cdot,D_{e4e5}\rangle=+2.36\times10^4$,
toep2 residual $-1.18\times10^4$, **deep toep3
($5.33\times10^{-5}\to2\times10^{-5}$): $-1.26\times10^5$ —
cut-signed on the deep escape**.  Face geometry (exact series):
$e_5$ grows positive-quadratically transverse to the pole–equator
stratum in every model family (F1 $\tfrac{256}{81}u^2$, F2
$\tfrac{256}{81}(1-r^2)u^2$, quartic $\tfrac{1024}{81}u^4$ at the
mode-4 corner, F3 $\tfrac{64}{81}t^2$, F4 $\tfrac{128}{81}u^2$; flat
on-stratum).  The weight role awaits the $e_5$-localized module
(degree-16-scale build).

### Weighted-(E1) projection (structure, not a bound route)

Classification at arity 2 (`solve_e1.py` Part F, 69 exact checks;
`audit_weighted_e1.py`, 11 test measures + 3 negative controls): the
weighted-admissible two-root pure layer shrinks 30→9/9/7/7 (majors);
the mode-2 tower $\hat C_n$ survives, mode-6 $\hat S_n$ dies; the
h2loc layer must satisfy the *unweighted* (E1) conditions
(second-order slackness).  Measured (degree 14, all-measures cone,
MOSEK dual):

| configuration | m | bound for $h_2E$ |
|---|---:|---:|
| baseline (no projection) | 806 | $-1.9198\times10^{-4}$ |
| pure layer → weighted-(E1) | 456 | **unbounded** (recession ray; $g_2$-mode, 55% $p_2\times$products) |
| h2loc layer → unweighted-(E1) | 747 | $-4.8667\times10^{-3}$ |
| both layers | 397 | unbounded |

Projected exports: `deg14_h2w_h2all_am_we1.dat-s` (m = 397),
`deg16_h2w_h2all_am_we1.dat-s` (m = 607),
`deg16_h2w_h2all_am_le1.dat-s` (m = 1173, GMP dual
$-4.5762\times10^{-3}$), `deg18_h2w_h2all_am_we1.dat-s` (m = 1015).
The projection's operative use is structural (pins the sharp limit
object, certifies sharpness-compatibility of added families), not
bound acceleration.

## B.7 Saturation phenomenology (conclusion)

Each convex cut family flattens the trace law over a window and is
then *saturated* on its equality manifold: the covariance (Jensen)
inequalities have kernel Cov = 0 (deterministic fibers), and the sharp
escape asymptotically saturates them — the same Cauchy–Schwarz-closure
behaviour predicted for the $\sqrt{h_2}$ leak.  At degree 14 the
saturation returned below $\sim5\times10^{-5}$ through
`h2loc_two_root_even_11`; at degree 16 the richer cone + families show
no re-steepening in the measured window, i.e. the saturation face
recedes with degree.  Convex families get saturated one by one; the
**irreducible repair is a multiplier vanishing on the saturation
face** (fiber-deterministic, quadrupole-invariant measures).  Proven
candidates: $e_5(I-A_2)$ with an $e_5$-localized module, the W-KKT
re-encoding, deeper-K fiber-Toeplitz at higher degree.  Zero remains
unattained at every finite structure measured: every wall is $<0$,
decaying super-geometrically under (degree × families); the endgames
are (A) a configuration whose wall lands at 0 with a flat law
(attainment), or (B) the multiplier target $(h_2+\kappa e_5)E$ with
zero attained by construction.

## B.8 Standing verdicts

Measured **dead** (do not re-try without new structure):

* scalar $p_k$-invariants as weights or cuts — powers of $h_2$,
  $g_\ell$ combinations, $e_k(I-A_2)$ scalar cuts: pair $\approx0$
  with the escape (pair-sector mass 0.07%);
* the bare $e_5$ cut at GMP precision (dual and selector identical to
  baseline; the double-precision gain was noise);
* theta windows and sup-caps at finite $\varepsilon$ (selector traces
  identical to controls; bound unchanged to 6 digits) — their one
  genuine effect is recession repair of the projected cone (146
  localized rows: `infeasible`), matched by 3 gap-cut rows;
* $\sqrt{h_2}$-adjunction (dominated: reproduces the Cauchy–Schwarz
  closure already in the spin-2 block);
* pure facial projections as bound routes (weakens the moment side;
  pure-layer weighted-(E1) cone dual-unbounded at every
  $\varepsilon$).

Measured **movers**:

* the $h_2$ weight with its localized module ($161\times$; lossless by
  the reduction lemma);
* degree (the strongest single lever: deg-14+v2 → deg-16+v3 moved the
  wall $\sim700\times$);
* the Jensen/fiber-Toeplitz family tower (first families active at
  selector optima; exponent 1.03 → 0.21, wall to
  $-5.25\times10^{-9}$).
