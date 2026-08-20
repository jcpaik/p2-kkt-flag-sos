"""Search a low-complexity exact certificate for the axial SO(2) zero mode."""

from __future__ import annotations

import cvxpy as cp
import numpy as np
import sympy as sp

import fermionic_axial_zero_mode_exact as axial


VARIABLES = sp.Matrix(axial.s)


def matrix(polynomial: sp.Expr) -> np.ndarray:
    return np.array(sp.hessian(sp.expand(polynomial), axial.s) / 2, dtype=float)


TARGET_EXACT = sp.hessian(axial.data()[-1], axial.s) / 2
TARGET = np.array(TARGET_EXACT, dtype=float)


def generators():
    s0, s1, s2, s3 = axial.s
    moment_determinants = {
        "D0": s0 * s2 - s1**2,
        "D1": s1 * s3 - s2**2,
        "Dbar": (s0 - s1) * (s2 - s3) - (s1 - s2) ** 2,
    }
    top = {
        "ka": -s0 + 6 * s1 - 5 * s2,
        "kb": -s0 + 9 * s1 - 10 * s2,
        "delta": -2 * s0 + 9 * s1 - 9 * s2,
    }
    bernstein = {
        "b0": s0 - 3 * s1 + 3 * s2 - s3,
        "b1": s1 - 2 * s2 + s3,
        "b2": s2 - s3,
        "b3": s3,
    }
    out = dict(moment_determinants)
    for top_name, top_form in top.items():
        for moment_name, moment_form in bernstein.items():
            out[f"{top_name}_{moment_name}"] = top_form * moment_form
    for first, (name1, form1) in enumerate(top.items()):
        for name2, form2 in list(top.items())[first:]:
            out[f"{name1}_{name2}"] = form1 * form2
    return out


GENERATORS = generators()


def solve():
    sos = cp.Variable((4, 4), symmetric=True)
    weights = cp.Variable(len(GENERATORS), nonneg=True)
    constraints = [sos >> 0]
    expression = sos
    for index, polynomial in enumerate(GENERATORS.values()):
        expression = expression + weights[index] * matrix(polynomial)
    constraints.append(expression == TARGET)
    problem = cp.Problem(cp.Minimize(cp.trace(sos) + cp.sum(weights)), constraints)
    value = problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-10,
        tol_feas=1e-10,
        tol_gap_rel=1e-10,
        max_iter=3000,
    )
    print("status/value", problem.status, value)
    print("sos eig", np.linalg.eigvalsh(sos.value))
    for name, coefficient in zip(GENERATORS, weights.value):
        if coefficient > 1e-7:
            print(name, coefficient)
    print("sos", sos.value)
    return problem.status, sos.value, weights.value


if __name__ == "__main__":
    solve()
