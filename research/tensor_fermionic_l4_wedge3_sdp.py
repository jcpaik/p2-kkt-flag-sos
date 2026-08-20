"""Sample-screen a higher-compound certificate on the 40D l2+l4 slice.

Search

    tr(G) Q(G) = tr(G) <W2,C2(G)> + <W3,C3(G)>,

with W2,W3 PSD.  This is a homogeneous cubic certificate valid for G PSD.
Random coefficient samples are followed by independent validation; a feasible
screen is only a candidate, while infeasibility already rejects this ansatz.
"""

import itertools
import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import block_diag, null_space

import tensor_fermionic_l4_compound_sdp as data
import tensor_wedge_relation as wedge


BASIS = np.stack(data.BASIS)
Q = data.Q
PAIRS = np.array(list(itertools.combinations(range(10), 2)))
TRIPLES = np.array(list(itertools.combinations(range(10), 3)))


def compound(G, indices):
    submatrices = G[indices[:, None, :, None], indices[None, :, None, :]]
    return np.linalg.det(submatrices)


def values(z):
    G = np.einsum("a,aij->ij", z, BASIS)
    return (
        G,
        compound(G, PAIRS),
        compound(G, TRIPLES),
        np.trace(G) * (z @ Q @ z),
    )


def solve(samples=700, seed=818):
    rng = np.random.default_rng(seed)
    rows = [values(rng.normal(size=len(BASIS))) for _ in range(samples)]
    W2 = cp.Variable((45, 45), symmetric=True)
    W3 = cp.Variable((120, 120), symmetric=True)
    constraints = [W2 >> 0, W3 >> 0]
    for G, C2, C3, target in rows:
        constraints.append(
            np.trace(G) * cp.sum(cp.multiply(W2, C2))
            + cp.sum(cp.multiply(W3, C3))
            == target
        )
    problem = cp.Problem(cp.Minimize(cp.trace(W2) + cp.trace(W3)), constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                value = problem.solve(
                    solver=solver,
                    tol_gap_abs=2e-8,
                    tol_feas=2e-8,
                    tol_gap_rel=2e-8,
                    max_iter=2000,
                )
            else:
                value = problem.solve(solver=solver, eps=2e-6, max_iters=300000)
            print(solver, problem.status, value)
            if W2.value is not None:
                print("eig W2", np.linalg.eigvalsh(W2.value)[:10])
                print("eig W3", np.linalg.eigvalsh(W3.value)[:20])
                errors = []
                for _ in range(200):
                    G, C2, C3, target = values(rng.normal(size=len(BASIS)))
                    errors.append(
                        np.trace(G) * np.sum(W2.value * C2)
                        + np.sum(W3.value * C3)
                        - target
                    )
                print("validation", np.max(np.abs(errors)))
                np.savez(
                    "research/tensor_fermionic_l4_wedge3_solution.npz",
                    W2=W2.value,
                    W3=W3.value,
                )
            if problem.status in ("optimal", "optimal_inaccurate"):
                break
        except cp.error.SolverError as error:
            print(solver, error)


def exterior_generator(generator, degree):
    indices = list(itertools.combinations(range(10), degree))
    index = {entry: i for i, entry in enumerate(indices)}
    out = np.zeros((len(indices), len(indices)))
    for column, source_tuple in enumerate(indices):
        for position, source in enumerate(source_tuple):
            for target in range(10):
                value = generator[target, source]
                if value == 0:
                    continue
                target_tuple = list(source_tuple)
                target_tuple[position] = target
                if len(set(target_tuple)) < degree:
                    continue
                inversions = sum(
                    target_tuple[a] > target_tuple[b]
                    for a in range(degree)
                    for b in range(a + 1, degree)
                )
                out[index[tuple(sorted(target_tuple))], column] += (-1) ** inversions * value
    return out


def symmetric_basis(dimension):
    out = []
    for i in range(dimension):
        E = np.zeros((dimension, dimension))
        E[i, i] = 1
        out.append(E)
    for i in range(dimension):
        for j in range(i + 1, dimension):
            E = np.zeros((dimension, dimension))
            E[i, j] = E[j, i] = 1 / np.sqrt(2)
            out.append(E)
    return out


GENERATORS = [block_diag(a, b) for a, b in zip(wedge.gw1, wedge.gw3)]


def invariant_exterior_data(degree):
    generators = [exterior_generator(generator, degree) for generator in GENERATORS]
    casimir = -sum(generator @ generator for generator in generators)
    eigenvalues, eigenvectors = np.linalg.eigh(casimir)
    full_basis = []
    local_blocks = []
    labels = []
    for ell in range(16):
        E = eigenvectors[:, np.abs(eigenvalues - ell * (ell + 1)) < 1e-7]
        if E.shape[1] == 0:
            continue
        restricted = [E.T @ generator @ E for generator in generators]
        symmetric = symmetric_basis(E.shape[1])
        equations = np.stack(
            [
                np.concatenate(
                    [(generator @ B - B @ generator).reshape(-1) for generator in restricted]
                )
                for B in symmetric
            ],
            axis=1,
        )
        kernel = (
            np.eye(len(symmetric))
            if np.linalg.norm(equations) < 1e-9
            else null_space(equations, rcond=1e-9)
        )
        local = []
        for column in range(kernel.shape[1]):
            B = sum(kernel[row, column] * symmetric[row] for row in range(len(symmetric)))
            B = (B + B.T) / 2
            local.append(B)
            full_basis.append(E @ B @ E.T)
        local_blocks.append(local)
        labels.append((ell, E.shape[1], len(local)))
    print("exterior", degree, labels)
    return full_basis, local_blocks


def solve_invariant(samples=180, seed=1818):
    basis2, blocks2 = invariant_exterior_data(2)
    basis3, blocks3 = invariant_exterior_data(3)
    rng = np.random.default_rng(seed)
    rows2 = []
    rows3 = []
    targets = []
    for _ in range(samples):
        G, C2, C3, target = values(rng.normal(size=len(BASIS)))
        rows2.append(np.trace(G) * np.array([np.sum(B * C2) for B in basis2]))
        rows3.append(np.array([np.sum(B * C3) for B in basis3]))
        targets.append(target)
    rows2 = np.stack(rows2)
    rows3 = np.stack(rows3)
    targets = np.array(targets)
    raw = np.concatenate([rows2, rows3], axis=1)
    coefficient, *_ = np.linalg.lstsq(raw, targets, rcond=None)
    print(
        "linear span rank/residual",
        np.linalg.matrix_rank(raw),
        np.linalg.norm(raw @ coefficient - targets),
        "variables",
        raw.shape[1],
    )

    variables2 = [cp.Variable(len(block)) for block in blocks2]
    variables3 = [cp.Variable(len(block)) for block in blocks3]
    constraints = []
    for variable, block in zip(variables2, blocks2):
        constraints.append(sum(variable[k] * block[k] for k in range(len(block))) >> 0)
    for variable, block in zip(variables3, blocks3):
        constraints.append(sum(variable[k] * block[k] for k in range(len(block))) >> 0)
    vector2 = cp.hstack(variables2)
    vector3 = cp.hstack(variables3)
    constraints.append(rows2 @ vector2 + rows3 @ vector3 == targets)
    objective = 0
    for variable, block in zip(variables2, blocks2):
        objective += cp.trace(sum(variable[k] * block[k] for k in range(len(block))))
    for variable, block in zip(variables3, blocks3):
        objective += cp.trace(sum(variable[k] * block[k] for k in range(len(block))))
    problem = cp.Problem(cp.Minimize(objective), constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                result = problem.solve(solver=solver, tol_gap_abs=1e-9, tol_feas=1e-9)
            else:
                result = problem.solve(solver=solver, eps=2e-7, max_iters=300000)
            print("invariant", solver, problem.status, result)
            if problem.status in ("optimal", "optimal_inaccurate"):
                break
        except cp.error.SolverError as error:
            print("invariant", solver, error)


if __name__ == "__main__":
    solve_invariant()
