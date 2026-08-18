#!/usr/bin/env python3
"""A/B test: does the theta-atom adjunction kill the projected escape ray?

Captures the exact `--find-ray` problem built by sos_search.solve (no
re-implementation, no sos_search surgery), then re-solves it with the
recession-form theta-atom constraints appended:

  (G_n)      L_n(r) >= 0                    for within-degree n, both families
  (W_{N',N}) sum_{N'<|n|<=N} q^{n^2} L_n(r) <= 0   (homogenized window cuts)

where L_n is the exact label expansion of the 1x1 two-root block
Q[Ghat_n] (theta_atoms.py, docs/ENRICHMENTS.md sections 2-3).  The
window cuts are the recession form of the one-sided truncation sandwich;
constants T^f_q(N') drop out along a ray.  A/B statuses:

  control            feasible ray (documented escape)         -> baseline
  + atom windows     infeasible                               -> ray dead

Also rebuilds the moment-bound problem (min <E, y> over the same cone,
y_const = 1) from the captured constraints for the boundedness A/B.

Usage:
  .venv/bin/python theta_ab.py --mode ray --variants control,f6,f2,both
  .venv/bin/python theta_ab.py --mode bound --variants control,both [--gap-cuts]
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

RAY_TOGGLES = [
    "--degree", "14", "--no-pointwise-sos", "--harmonics",
    "--three-point-flags", "--four-point-flags", "--two-root-flags",
    "--gradient", "--potential", "--hessian", "--global-tangent-gaps",
    "--rank-relations", "--e1-project", "e1_projection_deg14.json",
    "--scale-constraints",
]

UNPROJECTED_TOGGLES = [
    toggle
    for toggle in RAY_TOGGLES
    if toggle not in ("--e1-project", "e1_projection_deg14.json")
]

KKT_TOGGLES = (
    "--gradient", "--potential", "--hessian", "--global-tangent-gaps",
)

ALL_MEASURES_TOGGLES = [
    toggle for toggle in RAY_TOGGLES if toggle not in KKT_TOGGLES
]

CONES = {
    "kkt": RAY_TOGGLES,
    "allm": ALL_MEASURES_TOGGLES,
    "unprojected": UNPROJECTED_TOGGLES,
}


class _CapturedProblem(Exception):
    def __init__(self, problem: cp.Problem):
        self.problem = problem


ACTIVE_CONE = "kkt"


def build_args(extra: list[str]) -> argparse.Namespace:
    saved = sys.argv
    try:
        sys.argv = ["sos_search.py"] + CONES[ACTIVE_CONE] + extra
        return sos_search.parse_args()
    finally:
        sys.argv = saved


def dump_labels(gap_cuts: bool, tag: str) -> list[str]:
    path = Path("sdpa_runs") / f"theta_ab_dump_{tag}.json"
    extra = ["--dump-blocks", str(path)]
    if gap_cuts:
        extra.append("--gap-scalar-cuts")
    sos_search.solve(build_args(extra))
    with open(path) as handle:
        payload = json.load(handle)
    return payload["labels"]


def capture_ray_problem(gap_cuts: bool) -> cp.Problem:
    extra = ["--find-ray"]
    if gap_cuts:
        extra.append("--gap-scalar-cuts")
    args = build_args(extra)
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
    raise RuntimeError("find-ray problem was never constructed")


def atom_rows(labels: list[str], degree: int = 14, localized: bool = False):
    """Linear functionals l_{n,a} as index/coefficient arrays.

    localized=False: one row per generator (leaf Ghat_n itself).
    localized=True: one row per leaf Ghat_n * m_a over all monomial
    multipliers m_a fitting the degree cap (these squares are not
    (E1)-admissible, so they are new content on the projected cone).
    """
    index_of = {label: position for position, label in enumerate(labels)}
    rows = {}
    for family in (2, 6):
        n_min, n_max = theta_atoms.family_window(family, degree)
        for n in range(n_min, n_max + 1):
            if localized:
                leaves = theta_atoms.localized_leaves(family, n, degree)
            else:
                leaves = [((0, 0, 0), theta_atoms.hat_generator(family, n))]
            for multiplier, leaf in leaves:
                if sum(multiplier[:2]) % 2:
                    # odd leaf parity: the y-integral vanishes for every
                    # antipodal measure, so Q[leaf] is identically zero.
                    continue
                expansion = theta_atoms.expand_q_block(leaf)
                indices, values, missing = [], [], []
                for label, coefficient in expansion.items():
                    key = str(label)
                    if key in index_of:
                        indices.append(index_of[key])
                        values.append(float(coefficient))
                    else:
                        missing.append(key)
                if missing:
                    raise SystemExit(
                        f"labels missing from problem for family {family} "
                        f"n={n} m={multiplier}: {missing}"
                    )
                if not indices:
                    continue
                rows[(family, n, multiplier)] = (
                    np.array(indices, dtype=int),
                    np.array(values),
                )
    return rows


def linear_expression(variable, row):
    indices, values = row
    return values @ variable[indices]


def atom_constraints(variable, rows, families, q, degree: int = 14):
    """Recession-form (G) rows and window cuts for the requested families.

    Rows are keyed (family, n, alpha).  Windows are lumped over alpha;
    with the G rows present this forces the same recession zeros as the
    per-multiplier atoms (any positive combination does).
    """
    constraints = []
    described = []
    for family in families:
        n_min, n_max = theta_atoms.family_window(family, degree)
        keys = sorted(key for key in rows if key[0] == family)
        for key in keys:
            constraints.append(linear_expression(variable, rows[key]) >= 0)
        described.append(f"G[{family}] x{len(keys)}")
        cap = max(n_max, -n_min)
        for low in range(0, cap):
            window_keys = [key for key in keys if abs(key[1]) > low]
            if not window_keys:
                continue
            expression = sum(
                float(q) ** (key[1] * key[1])
                * linear_expression(variable, rows[key])
                for key in window_keys
            )
            constraints.append(expression <= 0)
            described.append(
                f"W[{family},q={q},({low},{cap}]] x{len(window_keys)}"
            )
    return constraints, described


def solve_with_fallback(problem: cp.Problem, tolerance: float = 1e-9):
    try:
        problem.solve(
            solver="MOSEK",
            mosek_params={
                "MSK_DPAR_INTPNT_CO_TOL_PFEAS": tolerance,
                "MSK_DPAR_INTPNT_CO_TOL_DFEAS": tolerance,
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": tolerance,
            },
        )
        return problem.status, problem.value, "MOSEK"
    except cp.error.SolverError:
        pass
    try:
        problem.solve(solver="CLARABEL")
        return problem.status, problem.value, "CLARABEL"
    except cp.error.SolverError:
        pass
    try:
        problem.solve(solver="SCS", max_iters=200000)
        return problem.status, problem.value, "SCS"
    except cp.error.SolverError:
        return "solver_failure", None, "none"


def get_ray_variable(problem: cp.Problem):
    for variable in problem.variables():
        if variable.name() == "ray_moments":
            return variable
    raise RuntimeError("ray_moments variable not found")


FAMILY_SETS = {
    "control": ((), False),
    "f2": ((2,), False),
    "f6": ((6,), False),
    "both": ((2, 6), False),
    "f2loc": ((2,), True),
    "f6loc": ((6,), True),
    "bothloc": ((2, 6), True),
}


def run_ray_ab(variants, q, gap_cuts: bool, out_path: Path):
    tag = ACTIVE_CONE + ("_gapcuts" if gap_cuts else "_nocuts")
    labels = dump_labels(gap_cuts, tag)
    problem = capture_ray_problem(gap_cuts)
    variable = get_ray_variable(problem)
    if variable.shape[0] != len(labels):
        raise SystemExit(
            f"label count mismatch: dump {len(labels)}, "
            f"variable {variable.shape[0]}"
        )
    plain_rows = atom_rows(labels)
    localized_rows = (
        atom_rows(labels, localized=True)
        if any(FAMILY_SETS[name][1] for name in variants)
        else None
    )
    results = {}
    for variant in variants:
        families, localized = FAMILY_SETS[variant]
        rows = localized_rows if localized else plain_rows
        added, described = atom_constraints(variable, rows, families, q)
        candidate = cp.Problem(
            problem.objective, list(problem.constraints) + added
        )
        status, value, solver = solve_with_fallback(candidate)
        entry = {
            "status": status,
            "squared_norm": None if value is None else float(value),
            "solver": solver,
            "added_constraints": len(added),
            "atom_families": list(families),
            "q": str(q),
        }
        if (
            variable.value is not None
            and status in ("optimal", "optimal_inaccurate")
        ):
            ray = variable.value
            entry["theta_diagonal"] = {
                str(key): float(np.dot(values, ray[indices]))
                for key, (indices, values) in sorted(rows.items())
            }
            entry["support"] = int(np.count_nonzero(np.abs(ray) > 1e-7))
            entry["ray"] = {
                labels[index]: float(ray[index])
                for index in np.argsort(-np.abs(ray))
                if abs(ray[index]) > 1e-7
            }
        results[variant] = entry
        print(f"[ray/{tag}] {variant}: {status} "
              f"norm={entry['squared_norm']} solver={solver} "
              f"(+{len(added)} atom rows)")
    payload = {
        "mode": "ray",
        "gap_cuts": gap_cuts,
        "q": str(q),
        "labels": len(labels),
        "toggles": CONES[ACTIVE_CONE],
        "results": results,
    }
    out_path.write_text(json.dumps(payload, indent=1))
    print(f"wrote {out_path}")


def rebuild_bound_problem(problem: cp.Problem, labels: list[str]):
    """Moment-bound problem from the captured ray constraints.

    The find-ray constraint list ends with [constant == 0, target == -1];
    replace them by constant == 1 and objective min <E, y>.
    """
    variable = get_ray_variable(problem)
    core = list(problem.constraints[:-2])
    constant_index = labels.index(str(("constant",)))
    target = {
        str(("constant",)): -4.0 / 3.0,
        str(("pair", 2)): 20.0,
        str(("pair", 4)): -48.0,
        str(("pair", 6)): 32.0,
    }
    objective_vector = np.zeros(len(labels))
    for key, value in target.items():
        objective_vector[labels.index(key)] = value
    problem = cp.Problem(
        cp.Minimize(objective_vector @ variable),
        core + [variable[constant_index] == 1],
    )
    return problem, variable


def inhomogeneous_atom_constraints(variable, rows, families, q):
    """Full sandwich content in eliminated form: G rows plus, for each
    multiplier m, the nested window cuts of the per-m atom with their
    exact rational tail majorants (docs/ENRICHMENTS.md section 3).
    """
    added = []
    for family in families:
        keys = sorted(key for key in rows if key[0] == family)
        for key in keys:
            added.append(linear_expression(variable, rows[key]) >= 0)
        groups: dict[tuple, list] = {}
        for key in keys:
            groups.setdefault(key[2], []).append(key)
        for _, group_keys in sorted(groups.items()):
            orders = [key[1] for key in group_keys]
            cap = max(max(orders), -min(orders))
            for low in range(0, cap):
                window = [key for key in group_keys if abs(key[1]) > low]
                if not window:
                    continue
                expression = sum(
                    float(q) ** (key[1] * key[1])
                    * linear_expression(variable, rows[key])
                    for key in window
                )
                majorant = float(theta_atoms.tail_majorant(family, q, low))
                added.append(expression <= majorant)
    return added


def run_bound_ab(variants, q, gap_cuts: bool, out_path: Path):
    tag = ACTIVE_CONE + ("_gapcuts" if gap_cuts else "_nocuts")
    labels = dump_labels(gap_cuts, tag)
    captured = capture_ray_problem(gap_cuts)
    base, variable = rebuild_bound_problem(captured, labels)
    plain_rows = atom_rows(labels)
    localized_rows = (
        atom_rows(labels, localized=True)
        if any(FAMILY_SETS[name][1] for name in variants)
        else None
    )
    results = {}
    for variant in variants:
        families, localized = FAMILY_SETS[variant]
        rows = localized_rows if localized else plain_rows
        added = (
            inhomogeneous_atom_constraints(variable, rows, families, q)
            if families
            else []
        )
        candidate = cp.Problem(
            base.objective, list(base.constraints) + added
        )
        status, value, solver = solve_with_fallback(candidate)
        entry = {
            "status": status,
            "bound": None if value is None else float(value),
            "solver": solver,
            "added_constraints": len(added),
            "atom_families": list(families),
            "q": str(q),
        }
        results[variant] = entry
        print(f"[bound/{tag}] {variant}: {status} bound={entry['bound']} "
              f"solver={solver} (+{len(added)} atom rows)")
    payload = {
        "mode": "bound",
        "gap_cuts": gap_cuts,
        "q": str(q),
        "labels": len(labels),
        "toggles": CONES[ACTIVE_CONE],
        "results": results,
    }
    out_path.write_text(json.dumps(payload, indent=1))
    print(f"wrote {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("ray", "bound"), default="ray")
    parser.add_argument("--variants", default="control,f6,f2,both")
    parser.add_argument("--q", default="1/2")
    parser.add_argument("--gap-cuts", action="store_true")
    parser.add_argument("--cone", choices=tuple(CONES), default="kkt")
    parser.add_argument("--out")
    args = parser.parse_args()

    global ACTIVE_CONE
    ACTIVE_CONE = args.cone
    variants = [name.strip() for name in args.variants.split(",")]
    q = Fraction(args.q)
    tag = args.cone + ("_gapcuts" if args.gap_cuts else "_nocuts")
    default_out = (
        Path("sdpa_runs")
        / f"theta_ab_{args.mode}_{tag}_q{q.numerator}_{q.denominator}.json"
    )
    out_path = Path(args.out) if args.out else default_out
    if args.mode == "ray":
        run_ray_ab(variants, q, args.gap_cuts, out_path)
    else:
        run_bound_ab(variants, q, args.gap_cuts, out_path)


if __name__ == "__main__":
    main()
