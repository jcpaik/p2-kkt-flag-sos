# Agent onboarding

**Mission**: prove $E(\mu)\ge0$ for every antipodal probability
measure on $S^2$, kernel $K(t)=32t^6-48t^4+20t^2-\tfrac43$, by finding
one exact, independently verifiable certificate.  Everything else in
this repository exists to produce, sharpen, or verify that object.

**Start here, in order**:

1. [PLAN.md](PLAN.md) — "Status" and "Next actions" at the top are the
   live state and the work queue.  Do the first unclaimed queue item.
2. [docs/EXACT_ZERO_PROGRAM.md](docs/EXACT_ZERO_PROGRAM.md) — the
   strategy document (why finite degree cannot land on zero, the
   weighted/rational-certificate program, the weight-vs-cut sign rule).
3. [docs/MATHEMATICAL_BACKGROUND.md](docs/MATHEMATICAL_BACKGROUND.md)
   and [docs/IMPLEMENTATION.md](docs/IMPLEMENTATION.md) — definitions
   (labels, flags, blocks, degree/arity) and the map from math to
   `sos_search.py` options.  Read on demand, not up front.

**Update discipline**: every measurement, positive or negative, gets a
dated entry in the relevant PLAN.md checklist item, with the data file
saved under `sdpa_runs/`.  Negative results (vacuous blocks, failed
candidates, dead rays) are recorded with the same care as gains — they
prune the search space for the next agent.  New derivations go in a
`docs/*.md` note cross-linked from PLAN.md.

## Environment facts

- Python: use `.venv/bin/python` (has cvxpy, numpy, sympy, mosek;
  MOSEK license at `~/mosek/mosek.lic`).  Tests: `python -m pytest -q`.
- High-precision SDP solver: `sdpa_runs/sdpa_gmp` (SDPA-GMP 7.1.3,
  built from nakatamaho/sdpa-gmp; rebuild recipe in
  [docs/H2_WEIGHTED_EXPERIMENT.md](docs/H2_WEIGHTED_EXPERIMENT.md) §1).
  Invoke with option-style flags only:
  `./sdpa_gmp -ds PROBLEM.dat-s -o OUT.result -p param_200bit.sdpa`
  (positional arguments are misparsed).  Parameter files
  `param_200bit.sdpa` (200-bit, `epsilonStar 1e-25`, for bounds) and
  `param_128bit.sdpa` (128-bit, `1e-16`, for large/selector solves)
  are in `sdpa_runs/`.
- Export pipeline: `sos_search.py --export-sdpa` writes exact rational
  SDPA files; **the reported bound = objValPrimal + objective_shift**,
  where objective_shift is printed in the export JSON (currently
  $-4/3$ for target $E$, $+2/3$ for target $h_2E$; legacy exports made before 2026-08-17 used the $(3/16)$-scaled target with shifts $-1/4$ and $1/8$).
- `sdpa_runs/` is gitignored and holds problems, results, parameter
  files, ray data, and exact expansion JSONs.  Keep it that way.

## Hardware rules (8 GB laptop, frequently swap-thrashed)

- Run **one** GMP solve at a time.  Parallel solves get OOM-killed.
- Before launching, check `sysctl vm.swapusage`; above ~6 GB used,
  expect kills or ~2 iterations/hour thrashing — defer the solve.
- Typical durations when healthy: m≈400 → ~6 min; m≈790 → ~40 min;
  m≈1250 (degree 16) → hours; budget accordingly and run solves in the
  background, reading results only from the `.result` file after the
  process exits.
- Never trust a mid-run status over the file contents: verify
  completion (`phase.value`, feasibility errors) in the result file.

## Non-negotiables

- Exact arithmetic at the boundaries: exports are rational; a
  numerical bound near zero is **never** a certificate.  Only
  rationally rounded and independently verified identities
  (`verify_certificate.py`, `verify_exact_structure.py`) count.
- Validity discipline: every new block or relation must be provably
  valid — for all measures (identities, squares, operator bounds), or
  clearly marked KKT-only (valid at minimizers).  The reduction lemma
  composes only with the all-measures cone; see the composition caveat
  in PLAN.md §5.
- When adding a candidate constraint or weight, first run the cheap
  pairing test against the current escape data (sign rule:
  positive ray pairing → weight, negative → cut; see
  [docs/GAP_CUTS_NOTE.md](docs/GAP_CUTS_NOTE.md)) before any big solve.
