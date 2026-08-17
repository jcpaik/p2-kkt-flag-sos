#!/usr/bin/env python3
"""Certificate-side epsilon-trace A/B: do theta atoms damp the pole?

The minimal-trace selector min sum tr Q_b s.t. e(Q) = E + eps is the
attainment diagnostic (docs/EXACT_ZERO_PROGRAM.md section 2.4): a
1/eps trace law means the sharp certificate escapes.  The theta-atom
adjunction gives certificates the extra generators

    lambda_c * ( T_c - sum_labels row_c[l] y_l ),   lambda_c >= 0,

(the window cuts of docs/THETA_ATOM_NOTE.md, valid at every measure),
which on the certificate side means each label equality gains
lambda_c * (T_c 1[l = const] - row_c[l]).  This script captures the
exact selector problem built by sos_search.solve (unscaled), rebuilds
the per-label equalities with the atom terms, and sweeps epsilon with
and without the atoms.  Double precision, structural comparison only.

Usage:
  .venv/bin/python theta_trace.py --epsilons 0.3,0.1,0.03,0.01 \
      --out sdpa_runs/theta_trace_ab.json
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

import cvxpy as cp
import numpy as np

import sos_search
import theta_atoms
from theta_export import atom_cut_rows

BASE_TOGGLES = [
    "--degree", "14", "--no-pointwise-sos", "--harmonics",
    "--three-point-flags", "--four-point-flags", "--two-root-flags",
    "--gradient", "--potential", "--hessian", "--global-tangent-gaps",
    "--rank-relations",
]


class _CapturedProblem(Exception):
    def __init__(self, problem: cp.Problem):
        self.problem = problem


def build_args(extra: list[str]) -> argparse.Namespace:
    saved = sys.argv
    try:
        sys.argv = ["sos_search.py"] + BASE_TOGGLES + extra
        return sos_search.parse_args()
    finally:
        sys.argv = saved


def dump_labels(tag: str) -> list[str]:
    path = Path("sdpa_runs") / f"theta_trace_dump_{tag}.json"
    sos_search.solve(build_args(["--dump-blocks", str(path)]))
    with open(path) as handle:
        return json.load(handle)["labels"]


def capture_selector(epsilon: float) -> cp.Problem:
    args = build_args(["--target-epsilon", str(epsilon)])
    original = cp.Problem.solve

    def interceptor(self, *positional, **keyword):
        raise _CapturedProblem(self)

    cp.Problem.solve = interceptor
    try:
        sos_search.solve(args)
    except _CapturedProblem as captured:
        return captured.problem
    finally:
        cp.Problem.solve = original
    raise RuntimeError("selector problem was never constructed")


def with_atoms(problem: cp.Problem, labels: list[str], cuts):
    """Rebuild the selector with atom multipliers lambda_c >= 0."""
    count = len(labels)
    equalities = list(problem.constraints[-count:])
    head = list(problem.constraints[:-count])
    lam = cp.Variable(len(cuts), nonneg=True, name="theta_multipliers")
    # per-label delta: sum_c lam_c * (T_c 1[l = const] - row_c[l])
    coefficient_matrix = np.zeros((count, len(cuts)))
    index_of = {label: position for position, label in enumerate(labels)}
    constant_index = index_of[str(("constant",))]
    for cut_index, (description, row, majorant) in enumerate(cuts):
        coefficient_matrix[constant_index, cut_index] += float(majorant)
        for label, value in row.items():
            coefficient_matrix[index_of[str(label)], cut_index] -= float(
                value
            )
    rebuilt = []
    for position, constraint in enumerate(equalities):
        delta = coefficient_matrix[position, :] @ lam
        rebuilt.append(
            constraint.args[0] + delta == constraint.args[1]
        )
    return cp.Problem(problem.objective, head + rebuilt), lam


def solve_selector(problem: cp.Problem):
    try:
        problem.solve(
            solver="MOSEK",
            mosek_params={
                "MSK_DPAR_INTPNT_CO_TOL_PFEAS": 1e-9,
                "MSK_DPAR_INTPNT_CO_TOL_DFEAS": 1e-9,
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": 1e-9,
            },
        )
        return problem.status, problem.value, "MOSEK"
    except cp.error.SolverError:
        pass
    try:
        problem.solve(solver="CLARABEL")
        return problem.status, problem.value, "CLARABEL"
    except cp.error.SolverError:
        return "solver_failure", None, "none"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epsilons", default="0.3,0.1,0.03,0.01")
    parser.add_argument("--q", default="1/2")
    parser.add_argument("--out", default="sdpa_runs/theta_trace_ab.json")
    args = parser.parse_args()

    epsilons = [float(text) for text in args.epsilons.split(",")]
    q = Fraction(args.q)
    labels = dump_labels("main")
    print(f"{len(labels)} labels; computing atom cut rows (q = {q}) ...")
    cuts = atom_cut_rows(q)
    print(f"{len(cuts)} window cuts")

    results = []
    for epsilon in epsilons:
        problem = capture_selector(epsilon)
        count = len(labels)
        if len(problem.constraints) < count:
            raise SystemExit("constraint/label count mismatch")
        status0, value0, solver0 = solve_selector(problem)
        atom_problem, lam = with_atoms(problem, labels, cuts)
        status1, value1, solver1 = solve_selector(atom_problem)
        active = None
        if lam.value is not None:
            active = int(np.count_nonzero(lam.value > 1e-6))
        entry = {
            "epsilon": epsilon,
            "control": {
                "status": status0,
                "trace": None if value0 is None else float(value0),
                "solver": solver0,
            },
            "atoms": {
                "status": status1,
                "trace": None if value1 is None else float(value1),
                "solver": solver1,
                "active_multipliers": active,
                "max_multiplier": (
                    None if lam.value is None else float(lam.value.max())
                ),
            },
        }
        results.append(entry)
        print(json.dumps(entry))

    Path(args.out).write_text(
        json.dumps(
            {"q": str(q), "toggles": BASE_TOGGLES, "sweep": results},
            indent=1,
        )
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
