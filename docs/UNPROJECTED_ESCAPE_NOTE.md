
## 9. Proof-carrying copy and the closing of we1 (2026-08-18)

* `deg16_am_toep3` (all-measures cone + 43 v3 families, m = 1265;
  **no KKT content** — composes with the reduction lemma): dual wall
  **−1.2252e−8**, only 2.3× the KKT-inclusive −5.25e−9.  The v3
  families subsume nearly all of the KKT families' strength: the
  mechanism is proof-grade.  Selector ladder in progress.
* `deg18_we1_toep3` (sharp-face projected + v3): selector at 1e−5
  **pUNBD** — that cone cannot certify h2E ≥ −1e−5 at any trace; the
  we1 route is closed for bounds (as the split verdict predicted).
* Unprojected degree-18: base exported (m = 2184); v3 stack building;
  extrapolated wall ≈ −1e−11 if the deg14→16 rate (~700×/2 degrees)
  sustains.

Proof-carrying ladder complete (2026-08-19): `sel16_am_toep3` traces
1808.99 / 3055.87 / 3949.27 at ε = 1e−4/1e−5/5e−6 (decade exponent
0.228, sub-decade 0.370), wall −1.2252e−8 — the same nearly-flat
behavior as the KKT copy, in the cone that composes with the reduction
lemma.  Degree-18 + v3 dual (m = 2192) launched.

## 10. The degree-18 + v3 wall: −2.21142e−14 (2026-08-19)

`deg18_h2w_h2all_toep3` (m = 2192, 1.04 GB, 200-bit, **pdOPT**):
bound = **−2.21142e−14** (primal; dual −2.21133e−14 — agreement to 5
digits, genuinely resolved, not a solver floor).  Wall adversary:
p2 = 1/3 + 2.9e−7 (h2 ≈ 4.3e−7), p4 = 0.28693, p6 = 0.26373.

The (degree × v3-families) cascade is **accelerating**:
deg-14+v2 −3.66e−6 → deg-16+v3 −5.25e−9 (700×) →
**deg-18+v3 −2.21e−14 (240,000×)**.  Steady non-attainment predicts
geometric decay; this acceleration is the signature expected if the
cone attains the sharp certificate at some finite degree.  The
discriminating test is the trace law approaching this wall (selectors
at 1e−5 … 1e−12): a plateau = attainment = the exactification
pipeline closes the weighted conjecture.  Ladder launched.

**Deg-18 stack trace law (2026-08-19, measured):** tr(1e−6) = 4860.95
(pFEAS), tr(1e−9) = 43838.0 (pFEAS): 2.08×/decade, exponent 0.318 —
**no attainment plateau**; the ~ε^−0.3 law persists at degree 18.
Verdict: the (degree × families) cascade drives the wall to zero
super-geometrically (−3.7e−6 → −5.3e−9 → −2.2e−14) but zero is a
limit, not attained; the residual fractional-order (~ε^{−0.3})
singularity is the surviving leak at every degree.  Endgame paths
unchanged: (i) limit-certificate extraction across the deg-14/16/18
certificate sequence; (ii) the e5-localized / multi-weight target
attaining zero by construction; (iii) the am-cone deg-18 solve for the
proof-carrying exact certificate at the ~1e−13 scale (staged, m=2216).

**Proof-carrying degree-18 wall (2026-08-20, pdOPT):**
`deg18_am_toep3` (all-measures cone + 43 v3 families, m = 2216):
bound = **−1.43335e−11** (primal/dual agree to 7 digits).  The
KKT-inclusive/all-measures gap widens with degree (2.3× at deg-16,
648× at deg-18 — the gradient rows carry real weight at 18).  Exact
rounding to the all-measures theorem h2E ≥ −3e−11 dispatched — the
cycle's capstone certificate (~500× beyond the committed −1.48e−8).
