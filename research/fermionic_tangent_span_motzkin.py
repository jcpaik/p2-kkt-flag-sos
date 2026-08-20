"""Diagnose the Motzkin defect after imposing the exact tangent moment span.

This is deliberately a numerical discovery script.  It fixes the natural
orthonormal basis of ``Sym^2_0(R^3)`` with

    e0 = diag(1,-1,0)/sqrt(2)

and constructs the 28-dimensional linear span of the tangent Pluecker
projectors.  After additionally imposing that the one-particle contraction
is diagonal, it minimizes the exact pinched Motzkin functional over the
remaining linear coordinates.  A later exact verifier should replace the
sample/SVD construction if this exposes a usable identity.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize


def h2_basis() -> list[np.ndarray]:
    out = [
        np.diag([1.0, -1.0, 0.0]) / np.sqrt(2),
        np.diag([1.0, 1.0, -2.0]) / np.sqrt(6),
    ]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        matrix = np.zeros((3, 3))
        matrix[i, j] = matrix[j, i] = 1 / np.sqrt(2)
        out.append(matrix)
    return out


H2 = h2_basis()
EDGES = list(combinations(range(5), 2))
SYMMETRIC_PAIRS = [(i, j) for i in range(10) for j in range(i, 10)]
DIAGONAL_INDICES = [SYMMETRIC_PAIRS.index((i, i)) for i in range(10)]
OFF_INDICES = [i for i in range(55) if i not in DIAGONAL_INDICES]


def pluecker(x: np.ndarray, basis: list[np.ndarray] | None = None) -> np.ndarray:
    if basis is None:
        basis = H2
    return np.array(
        [2 * x @ np.cross(basis[i] @ x, basis[j] @ x) for i, j in EDGES]
    )


def symmetric_coordinates(matrix: np.ndarray) -> np.ndarray:
    return np.array(
        [matrix[i, j] * (np.sqrt(2) if i != j else 1) for i, j in SYMMETRIC_PAIRS]
    )


def symmetric_matrix(coordinates: np.ndarray) -> np.ndarray:
    matrix = np.zeros((10, 10))
    for value, (i, j) in zip(coordinates, SYMMETRIC_PAIRS):
        matrix[i, j] = matrix[j, i] = value / (np.sqrt(2) if i != j else 1)
    return matrix


def contraction_matrix() -> np.ndarray:
    edge_index = {edge: index for index, edge in enumerate(EDGES)}

    def oriented_edge(i: int, j: int) -> tuple[int | None, int]:
        if i < j:
            return edge_index[i, j], 1
        if j < i:
            return edge_index[j, i], -1
        return None, 0

    gamma = np.zeros((25, 55))
    for column in range(55):
        state = symmetric_matrix(np.eye(55)[column])
        marginal = np.zeros((5, 5))
        for i in range(5):
            for k in range(5):
                for j in range(5):
                    ij, left_sign = oriented_edge(i, j)
                    kj, right_sign = oriented_edge(k, j)
                    if left_sign and right_sign:
                        marginal[i, k] += left_sign * right_sign * state[ij, kj]
        gamma[:, column] = marginal.reshape(-1)
    return gamma


def tangent_span(
    samples: int = 300,
    seed: int = 82026,
    basis: list[np.ndarray] | None = None,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    columns = []
    for _ in range(samples):
        x = rng.normal(size=3)
        x /= np.linalg.norm(x)
        z = pluecker(x, basis)
        columns.append(symmetric_coordinates(np.outer(z, z)))
    left, singular, _ = np.linalg.svd(np.stack(columns, axis=1), full_matrices=False)
    rank = int(np.sum(singular > 1e-10))
    assert rank == 28
    return left[:, :rank]


def reduced_maps(
    basis: list[np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    span = tangent_span(basis=basis)
    gamma = contraction_matrix()
    off_marginal = [5 * i + j for i in range(5) for j in range(i + 1, 5)]
    eigen_kernel = null_space(gamma[off_marginal] @ span)
    feasible_span = span @ eigen_kernel
    diagonal = feasible_span[DIAGONAL_INDICES]
    off = feasible_span[OFF_INDICES]
    left_null = null_space(diagonal.T)
    return diagonal, off, left_null


def optimize_span(
    restarts: int = 100,
    seed: int = 1820,
    basis: list[np.ndarray] | None = None,
    verbose: bool = True,
):
    diagonal, off, left_null = reduced_maps(basis)
    rng = np.random.default_rng(seed)
    adjacency = np.array(
        [
            [bool(set(edge) & set(other)) and edge != other for other in EDGES]
            for edge in EDGES
        ],
        dtype=float,
    )

    def occupations(p: np.ndarray) -> np.ndarray:
        return np.array(
            [sum(p[k] for k, edge in enumerate(EDGES) if i in edge) for i in range(5)]
        )

    def value(coordinates: np.ndarray) -> float:
        p = diagonal @ coordinates
        off_coordinates = off @ coordinates
        adjacent = 0.5 * p @ adjacency @ p
        return off_coordinates @ off_coordinates + 1 / 3 - adjacent

    constraints = [
        {"type": "eq", "fun": lambda coordinates: np.sum(diagonal @ coordinates) - 1},
        {
            "type": "ineq",
            "fun": lambda coordinates: occupations(diagonal @ coordinates)[0] - 2 / 3,
        },
    ]
    for i in range(1, 5):
        constraints.append(
            {
                "type": "ineq",
                "fun": lambda coordinates, i=i: occupations(diagonal @ coordinates)[0]
                - occupations(diagonal @ coordinates)[i],
            }
        )
    for edge in range(10):
        constraints.append(
            {
                "type": "ineq",
                "fun": lambda coordinates, edge=edge: (diagonal @ coordinates)[edge],
            }
        )

    # A least-squares preimage of the uniform ONB diagonal is feasible and
    # lies on the cap boundary; perturb it for independent starts.
    target = np.zeros(10)
    target[EDGES.index((0, 1))] = 1 / 3
    target[EDGES.index((0, 2))] = 1 / 3
    target[EDGES.index((1, 2))] = 1 / 3
    base, *_ = np.linalg.lstsq(diagonal, target, rcond=None)
    best = None
    for _ in range(restarts):
        initial = base + 0.02 * rng.normal(size=len(base))
        result = minimize(
            value,
            initial,
            method="SLSQP",
            constraints=constraints,
            options={"maxiter": 5000, "ftol": 1e-13},
        )
        if result.success and (best is None or result.fun < best.fun):
            best = result
    assert best is not None
    p = diagonal @ best.x
    if verbose:
        print("left-null diagonal relations", left_null.T)
        print("minimum", best.fun)
        print("p", p)
        print("occupations", occupations(p))
        print("off norm squared", np.sum((off @ best.x) ** 2))
    return best


def random_completions(count: int = 20, restarts: int = 20, seed: int = 991):
    """Keep e0 fixed and rotate its four-dimensional orthogonal complement."""
    rng = np.random.default_rng(seed)
    best = None
    for sample in range(count):
        raw = rng.normal(size=(4, 4))
        orthogonal, _ = np.linalg.qr(raw)
        basis = [H2[0]]
        basis.extend(
            [sum(orthogonal[j, i] * H2[j + 1] for j in range(4)) for i in range(4)]
        )
        result = optimize_span(
            restarts=restarts,
            seed=seed + sample,
            basis=basis,
            verbose=False,
        )
        if best is None or result.fun < best[0]:
            best = result.fun, orthogonal
            print("completion best", sample, best[0], flush=True)
    return best


if __name__ == "__main__":
    optimize_span()
