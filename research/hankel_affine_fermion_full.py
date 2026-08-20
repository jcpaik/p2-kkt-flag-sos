"""CCCP search in the exact affine-plus-upper fermionic relaxation.

Unlike tensor_fermionic_general_relaxation_opt.py, this uses a concrete
physical basis of H2 and the induced edge basis of wedge^2 H2.  This makes
the output suitable for recognizing exact algebraic counterexamples.
"""

import itertools
import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np

import hankel_affine_fermion_diag as base


pairs = base.pairs
W1 = base.W1

# F = Gamma vec(G), with row-major vectorization.
GAMMA = np.zeros((25, 100))
for a, (i, j) in enumerate(pairs):
    for b, (k, ell) in enumerate(pairs):
        col = 10 * a + b
        GAMMA[5 * i + k, col] += j == ell
        GAMMA[5 * i + ell, col] -= j == k
        GAMMA[5 * j + k, col] -= i == ell
        GAMMA[5 * j + ell, col] += i == k


def pi2_linear():
    out = np.zeros((25, 25))
    for i in range(5):
        for j in range(5):
            X = np.zeros((5, 5))
            X[i, j] = 1
            out[:, 5 * i + j] = base.pi2(X).reshape(-1)
    return out


PI2 = pi2_linear()


def t_linear():
    out = np.zeros((25, 9))
    for i in range(3):
        for j in range(3):
            A = np.zeros((3, 3))
            A[i, j] = 1
            out[:, 3 * i + j] = base.T_matrix(A).reshape(-1)
    return out


TMAP = t_linear()


def p_linear():
    out = np.zeros((25, 9))
    for a in range(3):
        for b in range(3):
            M = np.zeros((3, 3))
            M[a, b] = 1
            P = np.array(
                [
                    [
                        np.trace(M @ (base.E[i] @ base.E[j] + base.E[j] @ base.E[i]))
                        for j in range(5)
                    ]
                    for i in range(5)
                ]
            )
            out[:, 3 * a + b] = P.reshape(-1)
    return out


PMAP = p_linear()


def tangent_z(x):
    seed = np.eye(3)[np.argmin(np.abs(x))]
    u = np.cross(x, seed)
    u /= np.linalg.norm(u)
    v = np.cross(x, u)
    cols = []
    for y in (u, v):
        S = (np.outer(x, y) + np.outer(y, x)) / np.sqrt(2)
        cols.append(np.array([np.trace(A @ S) for A in base.E]))
    frame = np.stack(cols, axis=1)
    z = np.array(
        [frame[i, 0] * frame[j, 1] - frame[i, 1] * frame[j, 0] for i, j in pairs]
    )
    return z / np.linalg.norm(z)


def orbit_mixture(points, weights):
    return sum(w * np.outer(tangent_z(x), tangent_z(x)) for x, w in zip(points, weights))


def contractions(X):
    F = (GAMMA @ X.reshape(-1)).reshape(5, 5)
    M = 5 * W1.T @ X @ W1
    P = (PMAP @ M.reshape(-1)).reshape(5, 5)
    return F, M, P


def gap(X):
    F, _, _ = contractions(X)
    return np.sum(X * X) - 0.5 * np.sum(F * F) + 1 / 3


def cccp(initial, iterations=40, verbose=False, linear_factor=1.0):
    X = cp.Variable((10, 10), symmetric=True)
    xvec = cp.reshape(X, (100,), order="C")
    fvec = GAMMA @ xvec
    F = cp.reshape(fvec, (5, 5), order="C")
    M = 5 * W1.T @ X @ W1
    mvec = cp.reshape(M, (9,), order="C")
    P = cp.reshape(PMAP @ mvec, (5, 5), order="C")
    relation = PI2 @ fvec - (3 / 7) * TMAP @ cp.reshape(
        M - np.eye(3) / 3, (9,), order="C"
    )
    constraints = [
        X >> 0,
        cp.trace(X) == 1,
        cp.trace(W1.T @ X @ W1) == 1 / 5,
        relation == 0,
        P - F >> 0,
    ]
    current = initial.copy()
    for it in range(iterations):
        F0, _, _ = contractions(current)
        coefficient = (GAMMA.T @ F0.reshape(-1)).reshape(10, 10)
        objective = cp.Minimize(
            cp.sum_squares(X)
            - linear_factor * cp.sum(cp.multiply(coefficient, X))
        )
        problem = cp.Problem(objective, constraints)
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=1e-10,
            tol_feas=1e-10,
            tol_gap_rel=1e-10,
            max_iter=1000,
        )
        current = X.value
        if verbose:
            F0, M0, P0 = contractions(current)
            print(it, gap(current), np.linalg.eigvalsh(current), np.linalg.eigvalsh(P0 - F0))
    return gap(current), current


def main(restarts=30, iterations=40):
    rng = np.random.default_rng(20260820)
    best = None
    for restart in range(restarts):
        points = rng.normal(size=(20, 3))
        points /= np.linalg.norm(points, axis=1)[:, None]
        initial = orbit_mixture(points, rng.dirichlet(np.ones(len(points))))
        # The exact tangent majorization uses factor one.  The half-gradient
        # variant is also useful as a continuation heuristic and was the
        # convention in the earlier block-basis search.
        result = cccp(initial, iterations=iterations, linear_factor=0.5)
        if best is None or result[0] < best[0]:
            best = result
            F, M, P = contractions(result[1])
            print("best", restart, result[0])
            print("G", repr(result[1]))
            print("F", repr(F))
            print("M", repr(M))
            print("eig G", np.linalg.eigvalsh(result[1]))
            print("eig upper", np.linalg.eigvalsh(P - F), flush=True)


if __name__ == "__main__":
    main()
