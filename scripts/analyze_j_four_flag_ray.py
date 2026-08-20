"""Analyze the exact degree-six four-sample dual ray for the J target.

This is a research helper, not part of the certificate search CLI.  It reads
an unreduced block dump, repeats the exact ONB and pole--equator facial
reductions, and solves for a moment-side Farkas ray.
"""

from __future__ import annotations

import argparse
import ast
import json
from fractions import Fraction
from pathlib import Path

import cvxpy as cp
import numpy as np
import sympy as sp
from scipy.linalg import qr

import sos_search as ss


def load_blocks(path: Path):
    payload = json.loads(path.read_text())
    labels = [ast.literal_eval(value) for value in payload["labels"]]
    target = {
        ast.literal_eval(label): Fraction(str(value))
        for label, value in payload["target"].items()
    }
    blocks = {
        name: {
            ast.literal_eval(label): np.asarray(matrix, dtype=float)
            for label, matrix in matrices.items()
        }
        for name, matrices in payload["blocks"].items()
    }
    return labels, target, blocks


def reduce_on_face(blocks, label_value, *, onb: bool = False):
    reduced = {}
    for name, matrices in blocks.items():
        moment = (
            ss.exact_onb_moment_matrix(matrices)
            if onb
            else ss.exact_moment_matrix(matrices, label_value)
        )
        kernel_exact = ss.exact_nullspace(moment)
        if kernel_exact.cols == 0:
            continue
        kernel = np.asarray(kernel_exact, dtype=float)
        new_matrices = {
            label: kernel.T @ matrix @ kernel
            for label, matrix in matrices.items()
        }
        new_matrices = {
            label: matrix
            for label, matrix in new_matrices.items()
            if np.max(np.abs(matrix)) > 1e-13
        }
        if new_matrices:
            reduced[name] = new_matrices
    return reduced


def exact_matrix(matrix: np.ndarray) -> sp.Matrix:
    return sp.Matrix(
        [
            [
                sp.Rational(value.numerator, value.denominator)
                for value in (
                    ss.rationalize_float(float(entry)) for entry in row
                )
            ]
            for row in matrix
        ]
    )


def rationalize_ray(
    labels,
    target,
    blocks,
    values,
    denominator: int,
    rank_degree: int,
    print_moments: bool,
):
    indices = {label: index for index, label in enumerate(labels)}
    relations = ss.four_point_rank_relations(rank_degree)
    rows = []
    right_hand_side = []

    constant_row = [Fraction(0) for _ in labels]
    constant_row[indices[("constant",)]] = Fraction(1)
    rows.append(constant_row)
    right_hand_side.append(Fraction(0))

    target_row = [Fraction(0) for _ in labels]
    for label, coefficient in target.items():
        target_row[indices[label]] = coefficient
    rows.append(target_row)
    right_hand_side.append(Fraction(-1))

    for relation in relations:
        row = [Fraction(0) for _ in labels]
        for label, coefficient in relation.items():
            row[indices[label]] = ss.rationalize_float(float(coefficient))
        rows.append(row)
        right_hand_side.append(Fraction(0))

    equality = sp.Matrix(
        [[sp.Rational(v.numerator, v.denominator) for v in row]
         for row in rows]
    )
    rhs = sp.Matrix(
        [sp.Rational(v.numerator, v.denominator) for v in right_hand_side]
    )
    _, _, numerical_pivots = qr(
        np.asarray(equality, dtype=float),
        mode="economic",
        pivoting=True,
    )
    pivot_columns = list(map(int, numerical_pivots[: equality.rows]))
    free_columns = [index for index in range(len(labels))
                    if index not in pivot_columns]
    rational_values = [
        sp.Rational(round(float(value) * denominator), denominator)
        for value in values
    ]
    free_vector = sp.Matrix([rational_values[index] for index in free_columns])
    pivot_matrix = equality[:, pivot_columns]
    free_matrix = equality[:, free_columns]
    pivot_values = pivot_matrix.inv() * (rhs - free_matrix * free_vector)
    for index, value in zip(pivot_columns, pivot_values, strict=True):
        rational_values[index] = value

    vector = sp.Matrix(rational_values)
    assert equality * vector == rhs
    print("rational pivots", [repr(labels[index]) for index in pivot_columns])
    print("rational denominator", denominator)
    print("rational norm", float(sp.sqrt((vector.T * vector)[0])))
    for name, coefficient_matrices in blocks.items():
        matrix = sp.zeros(next(iter(coefficient_matrices.values())).shape[0])
        for label, coefficient_matrix in coefficient_matrices.items():
            matrix += vector[indices[label]] * exact_matrix(coefficient_matrix)
        leading_minors = [
            sp.factor(matrix[:order, :order].det())
            for order in range(1, matrix.rows + 1)
        ]
        if not all(value > 0 for value in leading_minors):
            raise AssertionError((name, leading_minors))
        print(
            "exact-pd",
            name,
            "leading-minor-signs",
            [sp.sign(value) for value in leading_minors],
            "min-eig",
            float(min(np.linalg.eigvalsh(np.asarray(matrix, dtype=float)))),
        )
    if print_moments:
        print("rational moments")
        for label, value in zip(labels, rational_values, strict=True):
            if value:
                print(repr(label), str(value))
    return rational_values


def solve_ray(
    labels,
    target,
    blocks,
    *,
    solver: str,
    max_margin: bool,
    rational_denominator: int | None,
    rank_degree: int,
    print_moments: bool,
):
    indices = {label: index for index, label in enumerate(labels)}
    y = cp.Variable(len(labels))
    constraints = []
    matrices = {}
    margin = cp.Variable() if max_margin else None
    for name, coefficient_matrices in blocks.items():
        matrix = sum(
            y[indices[label]] * coefficient
            for label, coefficient in coefficient_matrices.items()
        )
        matrices[name] = matrix
        if max_margin:
            constraints.append(matrix - margin * np.eye(matrix.shape[0]) >> 0)
        else:
            constraints.append(matrix >> 0)
    relations = ss.four_point_rank_relations(rank_degree)
    constraints.extend(
        sum(Fraction(coefficient) * y[indices[label]]
            for label, coefficient in relation.items()) == 0
        for relation in relations
    )
    constraints.extend(
        [
            y[indices[("constant",)]] == 0,
            sum(coefficient * y[indices[label]]
                for label, coefficient in target.items()) == -1,
        ]
    )
    if max_margin:
        constraints.append(cp.norm(y, 2) <= 1)
        problem = cp.Problem(cp.Maximize(margin), constraints)
    else:
        problem = cp.Problem(cp.Minimize(cp.sum_squares(y)), constraints)
    kwargs = {}
    if solver == "MOSEK":
        kwargs["mosek_params"] = {
            "MSK_DPAR_INTPNT_CO_TOL_PFEAS": 1e-12,
            "MSK_DPAR_INTPNT_CO_TOL_DFEAS": 1e-12,
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": 1e-12,
        }
    value = problem.solve(solver=solver, **kwargs)
    print("status", problem.status, "objective", value,
          "norm", np.linalg.norm(y.value))
    for name, matrix in matrices.items():
        eigenvalues = np.linalg.eigvalsh(matrix.value)
        print(
            name,
            matrix.shape,
            "eig",
            np.array2string(eigenvalues, precision=12, suppress_small=False),
        )
    if print_moments:
        print("moments")
        for label, value in sorted(
            zip(labels, y.value, strict=True), key=lambda item: -abs(item[1])
        ):
            if abs(value) > 1e-9:
                print(repr(label), repr(float(value)))
    if rational_denominator is not None:
        rationalize_ray(
            labels,
            target,
            blocks,
            y.value,
            rational_denominator,
            rank_degree,
            print_moments,
        )
    return y.value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dump", type=Path)
    parser.add_argument("--solver", default="CLARABEL")
    parser.add_argument("--max-margin", action="store_true")
    parser.add_argument("--rational-denominator", type=int)
    parser.add_argument("--rank-degree", type=int, default=2)
    parser.add_argument("--print-moments", action="store_true")
    args = parser.parse_args()

    labels, target, blocks = load_blocks(args.dump)
    blocks = reduce_on_face(blocks, ss.onb_label_value, onb=True)
    blocks = reduce_on_face(
        blocks,
        lambda label: ss.pole_equator_label_value(label, 0),
    )
    print("blocks", {name: next(iter(matrices.values())).shape[0]
                     for name, matrices in blocks.items()})
    solve_ray(
        labels,
        target,
        blocks,
        solver=args.solver,
        max_margin=args.max_margin,
        rational_denominator=args.rational_denominator,
        rank_degree=args.rank_degree,
        print_moments=args.print_moments,
    )


if __name__ == "__main__":
    main()
