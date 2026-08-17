#!/usr/bin/env python3
"""Double-precision selector A/B for the Jensen/Toeplitz blocks.

Rebuilds, directly in label space, the minimal-trace selector of the
KKT-inclusive weighted degree-14 problem from the capture written by
toeplitz_export.py:

    minimize   sum_b tr(Lambda_b)
    subject to target_L + eps [L = constant]
                 = sum_b <A^b_L, Lambda_b> + sum_f <A^f_L, F_f>
                   + sum_r lambda_r r_L          for every label L,
               Lambda_b >= 0  (PSD certificate blocks),
               F_f symmetric free, lambda_r free,

which is the label-space form of the GMP selector (sdpa_selector.py on
the exported problem; control values tr = 765.0 at eps = 1e-3 and
8133.6 at eps = 1e-4, docs/UNPROJECTED_ESCAPE_NOTE.md).  The A/B
compares the control block list against control + Jensen/Toeplitz
families.  Double precision only -- the decisive measurement stays the
GMP queue; this is the cheap trend check.

Usage:
  .venv/bin/python toeplitz_ab.py [--eps 1e-3,1e-4] [--out PATH]
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np

CAPTURE_PATH = Path(
    "/private/tmp/claude-501/-Users-jcpaik-Documents-research-"
    "p2-kkt-flag-sos/2d5da291-4f4d-44a8-bb7b-d3ca80702a32/scratchpad/"
    "toeplitz_capture_deg14_h2w.pkl"
)


def solve_selector(
    target: dict,
    psd_blocks: list,
    free_label_matrices: list,
    relations: list,
    epsilon: float,
    solver: str,
) -> dict:
    import cvxpy as cp
    import scipy.sparse as sp

    labels = set(target)
    for _, label_matrices in psd_blocks:
        labels.update(label_matrices)
    for label_matrices in free_label_matrices:
        labels.update(label_matrices)
    for relation in relations:
        labels.update(relation)
    ordered = sorted(labels, key=str)
    index = {label: i for i, label in enumerate(ordered)}
    constant = ("constant",)
    n_labels = len(ordered)

    def coupling_matrix(label_matrices, size):
        rows, cols, vals = [], [], []
        for label, matrix in label_matrices.items():
            symmetric = 0.5 * (matrix + matrix.T)
            flat = np.asarray(symmetric).reshape(-1)
            nz = np.nonzero(flat)[0]
            rows.extend([index[label]] * len(nz))
            cols.extend(nz.tolist())
            vals.extend(flat[nz].tolist())
        return sp.csr_matrix(
            (vals, (rows, cols)), shape=(n_labels, size * size)
        )

    psd_vars = []
    total = 0
    for name, label_matrices in psd_blocks:
        size = next(iter(label_matrices.values())).shape[0]
        variable = cp.Variable((size, size), symmetric=True, name=name)
        psd_vars.append((name, variable))
        total = total + coupling_matrix(label_matrices, size) @ cp.vec(
            variable, order="C"
        )
    for count, label_matrices in enumerate(free_label_matrices):
        size = next(iter(label_matrices.values())).shape[0]
        variable = cp.Variable(
            (size, size), symmetric=True, name=f"free_{count}"
        )
        total = total + coupling_matrix(label_matrices, size) @ cp.vec(
            variable, order="C"
        )
    if relations:
        rows, cols, vals = [], [], []
        for r, relation in enumerate(relations):
            for label, coefficient in relation.items():
                rows.append(index[label])
                cols.append(r)
                vals.append(float(coefficient))
        relation_matrix = sp.csr_matrix(
            (vals, (rows, cols)), shape=(n_labels, len(relations))
        )
        relation_vars = cp.Variable(
            len(relations), name="relation_multipliers"
        )
        total = total + relation_matrix @ relation_vars

    rhs = np.zeros(n_labels)
    for label, value in target.items():
        rhs[index[label]] = float(value)
    rhs[index[constant]] += epsilon

    constraints = [variable >> 0 for _, variable in psd_vars]
    constraints.append(total == rhs)

    objective = cp.Minimize(
        sum(cp.trace(variable) for _, variable in psd_vars)
    )
    problem = cp.Problem(objective, constraints)
    try:
        if solver == "MOSEK":
            problem.solve(solver=cp.MOSEK, verbose=False)
        else:
            problem.solve(solver=cp.CLARABEL, verbose=False)
    except Exception as error:  # noqa: BLE001
        return {"status": f"solver error: {error}", "trace": None}
    trace = None
    if problem.value is not None and np.isfinite(problem.value):
        trace = float(problem.value)
    per_block = {}
    if trace is not None:
        for name, variable in psd_vars:
            if variable.value is not None:
                per_block[name] = float(np.trace(variable.value))
    return {
        "status": problem.status,
        "trace": trace,
        "per_block_trace_top": dict(
            sorted(
                per_block.items(), key=lambda kv: -kv[1]
            )[:12]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eps", default="1e-3,1e-4")
    parser.add_argument("--solver", default="MOSEK")
    parser.add_argument(
        "--out", default="sdpa_runs/toeplitz_ab_selector.json"
    )
    args = parser.parse_args()

    with open(CAPTURE_PATH, "rb") as handle:
        capture = pickle.load(handle)
    target = {
        label: float(value)
        for label, value in capture["target"].items()
    }
    psd_blocks = capture["psd_blocks"]
    extra_blocks = capture["extra_blocks"]
    free_label_matrices = capture["free_label_matrices"]
    relations = capture["relations"]
    print(
        f"capture: {len(psd_blocks)} base blocks, "
        f"{len(extra_blocks)} Jensen blocks, "
        f"{len(free_label_matrices)} free families, "
        f"{len(relations)} relations"
    )

    results = []
    for epsilon_text in args.eps.split(","):
        epsilon = float(epsilon_text)
        for variant, block_list in (
            ("control", psd_blocks),
            ("toeplitz", psd_blocks + extra_blocks),
        ):
            outcome = solve_selector(
                target,
                block_list,
                free_label_matrices,
                relations,
                epsilon,
                args.solver,
            )
            trace_text = (
                f"{outcome['trace']:.4f}"
                if outcome["trace"] is not None
                else "-"
            )
            print(
                f"eps={epsilon:g} {variant:>9}: "
                f"status={outcome['status']} trace={trace_text}"
            )
            results.append(
                {
                    "eps": epsilon,
                    "variant": variant,
                    **outcome,
                }
            )
    Path(args.out).write_text(json.dumps(results, indent=1))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
