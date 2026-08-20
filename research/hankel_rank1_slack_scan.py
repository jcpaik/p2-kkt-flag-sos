"""Numerically classify rank-one SOS slacks for the Hankel KKT problem.

If the KKT slack is ``q=s*f^2`` for a ternary cubic ``f``, self-consistency
and complementarity reduce to criticality of an explicit SO(3)-invariant
quartic Rayleigh quotient.  This script builds that quotient directly from
the exact spherical monomial inner product and searches its critical points.

This is a diagnostic, not a proof.  All structural formulas used to build the
quotient are exact; scipy null spaces and the optimization are floating point.
"""

from __future__ import annotations

import math

import numpy as np
import torch
from scipy.linalg import null_space


torch.set_default_dtype(torch.float64)


def compositions(degree: int) -> list[tuple[int, int, int]]:
    return [
        (i, j, degree - i - j)
        for i in range(degree + 1)
        for j in range(degree - i + 1)
    ]


def sphere_average(alpha: tuple[int, int, int]) -> float:
    if any(entry % 2 for entry in alpha):
        return 0.0
    half = [entry // 2 for entry in alpha]
    numerator = 1
    for entry in half:
        numerator *= math.prod(range(1, 2 * entry, 2)) if entry else 1
    denominator = math.prod(range(1, 2 * sum(half) + 2, 2))
    return numerator / denominator


MONOMIALS = {degree: compositions(degree) for degree in range(7)}
INDEX = {
    degree: {alpha: index for index, alpha in enumerate(MONOMIALS[degree])}
    for degree in MONOMIALS
}


def harmonic_embedding(ell: int) -> np.ndarray:
    """Columns span r^(6-ell) H_ell in the sextic monomial basis."""
    if ell == 0:
        harmonic = np.eye(1)
    else:
        laplacian = np.zeros((len(MONOMIALS[ell - 2]), len(MONOMIALS[ell])))
        for column, alpha in enumerate(MONOMIALS[ell]):
            for coordinate in range(3):
                if alpha[coordinate] >= 2:
                    beta = list(alpha)
                    beta[coordinate] -= 2
                    laplacian[INDEX[ell - 2][tuple(beta)], column] += (
                        alpha[coordinate] * (alpha[coordinate] - 1)
                    )
        harmonic = null_space(laplacian)

    embedded = harmonic
    degree = ell
    for _ in range((6 - ell) // 2):
        multiply_r2 = np.zeros(
            (len(MONOMIALS[degree + 2]), len(MONOMIALS[degree]))
        )
        for column, alpha in enumerate(MONOMIALS[degree]):
            for coordinate in range(3):
                beta = list(alpha)
                beta[coordinate] += 2
                multiply_r2[INDEX[degree + 2][tuple(beta)], column] += 1
        embedded = multiply_r2 @ embedded
        degree += 2
    return embedded


EMBEDDINGS = {ell: harmonic_embedding(ell) for ell in (0, 2, 4, 6)}
CHANGE = np.hstack([EMBEDDINGS[ell] for ell in (0, 2, 4, 6)])
CHANGE_INVERSE = np.linalg.inv(CHANGE)

PROJECTORS: dict[int, np.ndarray] = {}
offset = 0
for ell in (0, 2, 4, 6):
    dimension = 2 * ell + 1
    selector = np.zeros((28, 28))
    selector[offset : offset + dimension, offset : offset + dimension] = np.eye(
        dimension
    )
    PROJECTORS[ell] = CHANGE @ selector @ CHANGE_INVERSE
    offset += dimension

GRAM6 = np.array(
    [
        [
            sphere_average(
                tuple(
                    MONOMIALS[6][row][coordinate]
                    + MONOMIALS[6][column][coordinate]
                    for coordinate in range(3)
                )
            )
            for column in range(28)
        ]
        for row in range(28)
    ]
)
GRAM3 = np.array(
    [
        [
            sphere_average(
                tuple(
                    MONOMIALS[3][row][coordinate]
                    + MONOMIALS[3][column][coordinate]
                    for coordinate in range(3)
                )
            )
            for column in range(10)
        ]
        for row in range(10)
    ]
)

SQUARE_MAP = np.zeros((28, 10, 10))
for row, alpha in enumerate(MONOMIALS[3]):
    for column, beta in enumerate(MONOMIALS[3]):
        total = tuple(alpha[k] + beta[k] for k in range(3))
        SQUARE_MAP[INDEX[6][total], row, column] += 1

LEGENDRE_COEFFICIENT = {
    2: 8 / 7,
    4: -384 / 385,
    6: 512 / 231,
}
INVERSE_KERNEL_WEIGHT = {
    ell: (2 * ell + 1) / LEGENDRE_COEFFICIENT[ell]
    for ell in LEGENDRE_COEFFICIENT
}
QUARTIC_MATRIX = sum(
    INVERSE_KERNEL_WEIGHT[ell]
    * (PROJECTORS[ell].T @ GRAM6 @ PROJECTORS[ell])
    for ell in (2, 4, 6)
)
SPHERE_INTEGRAL6 = np.array([sphere_average(alpha) for alpha in MONOMIALS[6]])


def slack_hankel(cubic: np.ndarray) -> tuple[float, float, np.ndarray]:
    cubic = cubic / np.sqrt(cubic @ GRAM3 @ cubic)
    square = np.einsum("kij,i,j->k", SQUARE_MAP, cubic, cubic)
    rayleigh = float(square @ QUARTIC_MATRIX @ square)
    scale = -1 / rayleigh
    functional = SPHERE_INTEGRAL6 + scale * GRAM6 @ sum(
        INVERSE_KERNEL_WEIGHT[ell] * (PROJECTORS[ell] @ square)
        for ell in (2, 4, 6)
    )
    hankel = np.zeros((10, 10))
    for row, alpha in enumerate(MONOMIALS[3]):
        for column, beta in enumerate(MONOMIALS[3]):
            total = tuple(alpha[k] + beta[k] for k in range(3))
            hankel[row, column] = functional[INDEX[6][total]]
    return rayleigh, scale, hankel


def scan(restarts: int = 120, steps: int = 3000) -> None:
    quartic = torch.tensor(QUARTIC_MATRIX)
    gram3 = torch.tensor(GRAM3)
    square_map = torch.tensor(SQUARE_MAP)
    records = []

    for seed in range(restarts):
        torch.manual_seed(seed)
        raw = torch.randn(10, requires_grad=True)
        optimizer = torch.optim.Adam([raw], lr=0.02)
        stationarity = None
        for _ in range(steps):
            optimizer.zero_grad()
            cubic = raw / torch.sqrt(raw @ gram3 @ raw)
            square = torch.einsum("kij,i,j->k", square_map, cubic, cubic)
            rayleigh = square @ quartic @ square
            gradient = torch.autograd.grad(rayleigh, cubic, create_graph=True)[0]
            natural_gradient = torch.linalg.solve(gram3, gradient)
            tangent_gradient = natural_gradient - (gradient @ cubic) * cubic
            stationarity = tangent_gradient @ gram3 @ tangent_gradient
            stationarity.backward()
            optimizer.step()

        with torch.no_grad():
            cubic = raw / torch.sqrt(raw @ gram3 @ raw)
            residual = float(stationarity)
            rayleigh, scale, hankel = slack_hankel(cubic.numpy())
            if residual < 1e-10:
                records.append(
                    (
                        rayleigh,
                        float(np.linalg.eigvalsh(hankel)[0]),
                        residual,
                        scale,
                    )
                )

    clusters: list[list[tuple[float, float, float, float]]] = []
    for record in sorted(records):
        if not clusters or abs(record[0] - clusters[-1][0][0]) > 1e-5:
            clusters.append([record])
        else:
            clusters[-1].append(record)

    print("converged", len(records), "clusters", len(clusters))
    for cluster in clusters:
        print(
            "R",
            cluster[0][0],
            "count",
            len(cluster),
            "minimum Hankel eigenvalue range",
            min(record[1] for record in cluster),
            max(record[1] for record in cluster),
        )


if __name__ == "__main__":
    scan()
