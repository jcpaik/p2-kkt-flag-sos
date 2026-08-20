"""Audit E - (4/9) h2 on the Hankel and 40D fermionic relaxations.

Here h2=(3 p2-1)/2.  Thus the strengthened target is the quadratic
pair-moment functional

    W = 32 p6 - 48 p4 + (58/3) p2 - 10/9.

The script uses exact linear constraints (up to the numerical CG bases used
by tensor_fermionic_general_relaxation_opt.py) and CCCP/Frank--Wolfe only as
counterexample searches.  A negative value is an obstruction for the
corresponding relaxation, not for genuine measures unless the state is
separable/tangent.
"""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np

import tensor_fermionic_general_relaxation_opt as rel
import tensor_hankel_fw as hankel


def fermionic_weighted_value(X):
    F = (rel.F_LINEAR @ X.reshape(-1)).reshape(5, 5)
    Q = np.sum(X * X) - np.sum(F * F) / 2 + 1 / 3
    A = X[:3, :3]
    h2 = (75 / 2) * np.sum((A - np.eye(3) / 15) ** 2)
    return 8 * Q - (4 / 9) * h2


def fermionic_cccp(initial, iterations=100):
    X = cp.Variable((10, 10), symmetric=True)
    xvec = cp.reshape(X, (100,), order="C")
    F = cp.reshape(rel.F_LINEAR @ xvec, (5, 5), order="C")
    constraints = [
        X >> 0,
        cp.trace(X) == 1,
        cp.trace(X[:3, :3]) == 1 / 5,
        rel.LINEAR_CONSTRAINT_REDUCED @ xvec == 0,
        rel.L4_REDUCED @ xvec == 0,
    ]
    current = initial.copy()
    for _ in range(iterations):
        F0 = (rel.F_LINEAR @ current.reshape(-1)).reshape(5, 5)
        A0 = current[:3, :3]
        # W/8 = ||X||^2 - ||F||^2/2 + 1/3
        #       - (25/12)||A-I/15||^2.
        # Linearize both concave quadratic terms.
        flinear = (rel.F_LINEAR.T @ F0.reshape(-1)).reshape(10, 10)
        center = (flinear + flinear.T) / 4
        alinear = np.zeros((10, 10))
        # The tangent of -c||A-A_*||^2 has linear coefficient
        # -2c<A0-A_*,A>, hence contributes c(A0-A_*) to the square center.
        alinear[:3, :3] = (25 / 12) * (A0 - np.eye(3) / 15)
        center += alinear
        problem = cp.Problem(cp.Minimize(cp.sum_squares(X - center)), constraints)
        try:
            problem.solve(
                solver="CLARABEL",
                tol_gap_abs=1e-10,
                tol_feas=1e-10,
                tol_gap_rel=1e-10,
                max_iter=1000,
            )
        except cp.error.SolverError:
            problem.solve(solver="SCS", eps=1e-8, max_iters=300000)
        current = X.value
    return fermionic_weighted_value(current), current


def fermionic_extreme_restarts(restarts=100, seed=20260820):
    rng = np.random.default_rng(seed)
    X = cp.Variable((10, 10), symmetric=True)
    xvec = cp.reshape(X, (100,), order="C")
    direction = cp.Parameter((10, 10), symmetric=True)
    constraints = [
        X >> 0,
        cp.trace(X) == 1,
        cp.trace(X[:3, :3]) == 1 / 5,
        rel.LINEAR_CONSTRAINT_REDUCED @ xvec == 0,
        rel.L4_REDUCED @ xvec == 0,
    ]
    extreme_problem = cp.Problem(
        cp.Minimize(cp.sum(cp.multiply(direction, X))), constraints
    )
    best = None
    for restart in range(restarts):
        raw = rng.normal(size=(10, 10))
        direction.value = (raw + raw.T) / 2
        extreme_problem.solve(solver="CLARABEL")
        value, state = fermionic_cccp(X.value)
        if best is None or value < best[0]:
            best = value, state
            F = (rel.F_LINEAR @ state.reshape(-1)).reshape(5, 5)
            print(
                "fermionic",
                restart,
                value,
                "rank spectrum",
                np.linalg.eigvalsh(state),
                "F spectrum",
                np.linalg.eigvalsh(F),
                flush=True,
            )
    return best


def hankel_weighted_value(y):
    rho = hankel.matrices(y)
    return (
        32 * np.sum(rho[3] ** 2)
        - 48 * np.sum(rho[2] ** 2)
        + (58 / 3) * np.sum(rho[1] ** 2)
        - 10 / 9
    )


def hankel_weighted_gradient(y):
    rho = hankel.matrices(y)
    return (
        64 * np.einsum("ab,abk->k", rho[3], hankel.maps[3])
        - 96 * np.einsum("ab,abk->k", rho[2], hankel.maps[2])
        + (116 / 3) * np.einsum("ab,abk->k", rho[1], hankel.maps[1])
    )


def hankel_restarts(restarts=50, iterations=500, seed=20260821):
    rng = np.random.default_rng(seed)
    starts = [np.array([hankel.uniform_moment(a) for a in hankel.deg6])]
    for _ in range(restarts - 1):
        hankel.direction.value = rng.normal(size=len(hankel.deg6))
        hankel.problem.solve(solver="CLARABEL")
        starts.append(hankel.y_variable.value.copy())
    best = None
    for number, y in enumerate(starts):
        for iteration in range(iterations):
            hankel.direction.value = hankel_weighted_gradient(y)
            hankel.problem.solve(solver="CLARABEL")
            vertex = hankel.y_variable.value.copy()
            delta = vertex - y
            e0 = hankel_weighted_value(y)
            e1 = hankel_weighted_value(y + delta)
            em = hankel_weighted_value(y + delta / 2)
            curvature = 2 * (e1 + e0 - 2 * em)
            slope = 4 * em - e1 - 3 * e0
            if curvature > 1e-14:
                step = np.clip(-slope / (2 * curvature), 0, 1)
            else:
                step = 1 if e1 < e0 else 0
            y += step * delta
            if step < 1e-11:
                break
        value = hankel_weighted_value(y)
        if best is None or value < best[0]:
            best = value, y.copy()
            print(
                "hankel",
                number,
                value,
                "rho3 spectrum",
                np.linalg.eigvalsh(hankel.matrices(y)[3]),
                flush=True,
            )
    return best


if __name__ == "__main__":
    hankel_restarts()
    fermionic_extreme_restarts()
