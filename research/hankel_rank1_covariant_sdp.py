"""Search a covariant sum-of-H-squares proof for the rank-one slack branch.

For the reconstructed middle catalecticant H_f=I+s C_f in the spherical
L2-orthonormal cubic basis, a PSD covariance Z on End(V) would prove the
branch immediately if

    sum_a ||B_a f||^2       = (32/105)||f||^2,
    sum_a C_f(B_a f,B_a f) = -||f||^4,

where Z=sum_a vec(B_a)vec(B_a)^T.  Then

    sum_a H_f(B_a f,B_a f)
      = ||f||^2(32/105-s||f||^2) = ||f||^2 E.

This script restricts Z to the exact SO(3) commutant (30 symmetric
parameters) and tests coefficient equality by generic samples.
"""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

import hankel_rank1_slack_invariant as data


def action(S):
    """Raw-monomial matrix of f -> (Sx).grad f."""
    O = np.zeros((10, 10))
    for j, alpha in enumerate(data.C3):
        for i in range(3):
            if alpha[i] == 0:
                continue
            for k in range(3):
                beta = list(alpha)
                beta[i] -= 1
                beta[k] += 1
                O[data.C3.index(tuple(beta)), j] += alpha[i] * S[i, k]
    return O


ROTATIONS = []
for i, j in ((0, 1), (0, 2), (1, 2)):
    S = np.zeros((3, 3))
    S[i, j], S[j, i] = 1, -1
    ROTATIONS.append(S)

# Whiten the spherical L2 metric on cubics.
eigenvalues, eigenvectors = np.linalg.eigh(data.GRAM3)
SQRT_G = eigenvectors @ np.diag(np.sqrt(eigenvalues)) @ eigenvectors.T
INV_SQRT_G = eigenvectors @ np.diag(1 / np.sqrt(eigenvalues)) @ eigenvectors.T
GENERATORS_V = [SQRT_G @ action(S) @ INV_SQRT_G for S in ROTATIONS]


def end_generator(generator):
    out = np.zeros((100, 100))
    for column in range(100):
        B = np.zeros((10, 10))
        B.reshape(-1)[column] = 1
        out[:, column] = (generator @ B - B @ generator).reshape(-1)
    return out


GENERATORS = [end_generator(g) for g in GENERATORS_V]


def symmetric_basis(n):
    out = []
    for i in range(n):
        E = np.zeros((n, n))
        E[i, i] = 1
        out.append(E)
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n))
            E[i, j] = E[j, i] = 1 / np.sqrt(2)
            out.append(E)
    return out


def invariant_basis():
    casimir = -sum(g @ g for g in GENERATORS)
    values, vectors = np.linalg.eigh(casimir)
    full = []
    labels = []
    for ell in range(7):
        E = vectors[:, np.abs(values - ell * (ell + 1)) < 2e-7]
        local_generators = [E.T @ g @ E for g in GENERATORS]
        symmetric = symmetric_basis(E.shape[1])
        equations = np.stack(
            [
                np.concatenate([(g @ B - B @ g).reshape(-1) for g in local_generators])
                for B in symmetric
            ],
            axis=1,
        )
        kernel = (
            np.eye(len(symmetric))
            if np.linalg.norm(equations) < 1e-7
            else null_space(equations, rcond=1e-8)
        )
        labels.append((ell, E.shape[1], kernel.shape[1]))
        for column in range(kernel.shape[1]):
            B = sum(kernel[row, column] * symmetric[row] for row in range(len(symmetric)))
            full.append(E @ B @ E.T)
    print("commutant", labels, "total", len(full))
    return full


BASIS = invariant_basis()


def c_matrix(z):
    raw = INV_SQRT_G @ z
    square = data.square_coefficients(raw)
    C = np.array(
        [
            [data.D_MATRIX[data.PRODUCT_INDEX[i, j]] @ square for j in range(10)]
            for i in range(10)
        ]
    )
    return INV_SQRT_G.T @ C @ INV_SQRT_G


def sample_rows(z):
    # T maps vec_C(B) to Bz.
    T = np.kron(np.eye(10), z.reshape(1, 10))
    K0 = T.T @ T
    K1 = T.T @ c_matrix(z) @ T
    return (
        np.array([np.sum(B * K0) for B in BASIS]),
        np.array([np.sum(B * K1) for B in BASIS]),
        (32 / 105) * (z @ z),
        -(z @ z) ** 2,
    )


def solve(samples=100, seed=441):
    rng = np.random.default_rng(seed)
    rows = [sample_rows(rng.normal(size=10)) for _ in range(samples)]
    A0 = np.stack([row[0] for row in rows])
    A1 = np.stack([row[1] for row in rows])
    b0 = np.array([row[2] for row in rows])
    b1 = np.array([row[3] for row in rows])
    coefficients = cp.Variable(len(BASIS))
    Z = sum(coefficients[i] * BASIS[i] for i in range(len(BASIS)))
    constraints = [Z >> 0, A0 @ coefficients == b0, A1 @ coefficients == b1]
    problem = cp.Problem(cp.Minimize(cp.trace(Z)), constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                result = problem.solve(
                    solver=solver,
                    tol_gap_abs=1e-9,
                    tol_feas=1e-9,
                    tol_gap_rel=1e-9,
                    max_iter=2000,
                )
            else:
                result = problem.solve(solver=solver, eps=1e-7, max_iters=300000)
            print(solver, problem.status, result)
            if coefficients.value is not None:
                Zv = sum(coefficients.value[i] * BASIS[i] for i in range(len(BASIS)))
                print("Z eigenvalues", np.linalg.eigvalsh(Zv))
                errors = []
                for _ in range(100):
                    r0, r1, t0, t1 = sample_rows(rng.normal(size=10))
                    errors.append((r0 @ coefficients.value - t0, r1 @ coefficients.value - t1))
                print("validation", np.max(np.abs(errors), axis=0))
            if problem.status in ("optimal", "optimal_inaccurate"):
                break
        except cp.error.SolverError as error:
            print(solver, error)


PAIRS2 = [(i, j) for i in range(10) for j in range(i, 10)]


def quadratic_features(z):
    return np.array([z[i] * z[j] for i, j in PAIRS2])


def solve_with_quartic_sos(samples=900, seed=1441):
    """Relax the failed identity to F1 <= -||f||^4 via a quartic SOS."""
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(samples):
        z = rng.normal(size=10)
        r0, r1, t0, t1 = sample_rows(z)
        q = quadratic_features(z)
        rows.append((r0, r1, t0, t1, np.outer(q, q)))
    A0 = np.stack([row[0] for row in rows])
    A1 = np.stack([row[1] for row in rows])
    b0 = np.array([row[2] for row in rows])
    b1 = np.array([row[3] for row in rows])
    gram_rows = np.stack([row[4].reshape(-1) for row in rows])

    coefficients = cp.Variable(len(BASIS))
    Z = sum(coefficients[i] * BASIS[i] for i in range(len(BASIS)))
    W = cp.Variable((len(PAIRS2), len(PAIRS2)), symmetric=True)
    constraints = [
        Z >> 0,
        W >> 0,
        A0 @ coefficients == b0,
        gram_rows @ cp.reshape(W, (len(PAIRS2) ** 2,), order="C") == b1 - A1 @ coefficients,
    ]
    problem = cp.Problem(cp.Minimize(cp.trace(Z) + cp.trace(W)), constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                result = problem.solve(
                    solver=solver,
                    tol_gap_abs=2e-8,
                    tol_feas=2e-8,
                    tol_gap_rel=2e-8,
                    max_iter=3000,
                )
            else:
                result = problem.solve(solver=solver, eps=2e-6, max_iters=500000)
            print("quartic", solver, problem.status, result)
            if coefficients.value is not None and W.value is not None:
                Zv = sum(coefficients.value[i] * BASIS[i] for i in range(len(BASIS)))
                print("eig Z", np.linalg.eigvalsh(Zv)[:20])
                print("eig W", np.linalg.eigvalsh(W.value)[:20])
                errors = []
                for _ in range(200):
                    z = rng.normal(size=10)
                    r0, r1, t0, t1 = sample_rows(z)
                    q = quadratic_features(z)
                    errors.append(
                        (r0 @ coefficients.value - t0, r1 @ coefficients.value + q @ W.value @ q - t1)
                    )
                print("quartic validation", np.max(np.abs(errors), axis=0))
            if problem.status in ("optimal", "optimal_inaccurate"):
                break
        except cp.error.SolverError as error:
            print("quartic", solver, error)


def critical_gauge_basis():
    """Cross-covariances encoding arbitrary invariant linear multipliers of Pf=0."""
    raw_h1 = np.zeros((10, 3))
    for k in range(3):
        for i in range(3):
            alpha = [0, 0, 0]
            alpha[k] += 1
            alpha[i] += 2
            raw_h1[data.C3.index(tuple(alpha)), k] += 1
    q1, _ = np.linalg.qr(SQRT_G @ raw_h1)
    p1 = q1 @ q1.T
    p3 = np.eye(10) - p1
    identity = np.eye(10).reshape(-1)
    out = []
    for multiplier in (p1, p3):
        vector = multiplier.reshape(-1)
        out.append(np.outer(identity, vector) + np.outer(vector, identity))
    return out


GAUGE = critical_gauge_basis()


def solve_mod_criticality(samples=120, seed=2441):
    """Search the exact identity modulo the stationarity equation P_f f=0."""
    rng = np.random.default_rng(seed)
    rows = [sample_rows(rng.normal(size=10)) for _ in range(samples)]
    A0 = np.stack([row[0] for row in rows])
    A1 = np.stack([row[1] for row in rows])
    G0 = []
    G1 = []
    for _ in range(samples):
        pass
    # Re-evaluate the same sample sequence so the gauge columns use exactly
    # the K0,K1 underlying A0,A1.
    rng = np.random.default_rng(seed)
    for _ in range(samples):
        z = rng.normal(size=10)
        T = np.kron(np.eye(10), z.reshape(1, 10))
        K0 = T.T @ T
        K1 = T.T @ c_matrix(z) @ T
        G0.append([np.sum(B * K0) for B in GAUGE])
        G1.append([np.sum(B * K1) for B in GAUGE])
    G0 = np.asarray(G0)
    G1 = np.asarray(G1)
    b0 = np.array([row[2] for row in rows])
    b1 = np.array([row[3] for row in rows])

    coefficients = cp.Variable(len(BASIS))
    gauge = cp.Variable(len(GAUGE))
    Z = sum(coefficients[i] * BASIS[i] for i in range(len(BASIS)))
    # Z-GAUGE has the target coefficients.  Z itself is a Gram covariance,
    # while GAUGE vanishes after contraction with P_f f=0.
    constraints = [
        Z >> 0,
        A0 @ coefficients - G0 @ gauge == b0,
        A1 @ coefficients - G1 @ gauge == b1,
    ]
    problem = cp.Problem(cp.Minimize(cp.trace(Z)), constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                result = problem.solve(
                    solver=solver,
                    tol_gap_abs=1e-9,
                    tol_feas=1e-9,
                    tol_gap_rel=1e-9,
                    max_iter=3000,
                )
            else:
                result = problem.solve(solver=solver, eps=1e-7, max_iters=500000)
            print("critical gauge", solver, problem.status, result)
            if coefficients.value is not None:
                Zv = sum(coefficients.value[i] * BASIS[i] for i in range(len(BASIS)))
                print("gauge", gauge.value)
                print("gauge eig Z", np.linalg.eigvalsh(Zv))
                errors = []
                for _ in range(100):
                    z = rng.normal(size=10)
                    r0, r1, t0, t1 = sample_rows(z)
                    T = np.kron(np.eye(10), z.reshape(1, 10))
                    K0 = T.T @ T
                    K1 = T.T @ c_matrix(z) @ T
                    q0 = np.array([np.sum(B * K0) for B in GAUGE])
                    q1 = np.array([np.sum(B * K1) for B in GAUGE])
                    errors.append(
                        (r0 @ coefficients.value - q0 @ gauge.value - t0,
                         r1 @ coefficients.value - q1 @ gauge.value - t1)
                    )
                print("gauge validation", np.max(np.abs(errors), axis=0))
            if problem.status in ("optimal", "optimal_inaccurate"):
                break
        except cp.error.SolverError as error:
            print("critical gauge", solver, error)


if __name__ == "__main__":
    solve()
    solve_with_quartic_sos()
    solve_mod_criticality()
