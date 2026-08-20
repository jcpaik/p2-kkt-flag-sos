"""Numerical audit of the scalar fermionic bridge for every spectrum of S.

Write

    S=(sqrt(3) S0-h H0)/sqrt(3+h^2),
    H=(h S0+sqrt(3) H0)/sqrt(3+h^2),   -1 <= h <= 1.

On the exact (sampled/SVD) tangent moment span, impose only FS=cS and test

    R0 - r_h^2/6 >= 0,

where

    r_h = 3 q0 (mass-(1-h^2) E[z^2]) - 2 mass,
    q0=3/(3+h^2).

This is a discovery script; an exact parametric verifier must follow any
positive numerical result.
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.linalg import null_space

import tensor_star_canonical_exact as canonical


PAIRS = canonical.PAIRS


def orbitals(h: float) -> list[np.ndarray]:
    old = [np.array(matrix, dtype=float) for matrix in canonical.ORB]
    d = np.sqrt(3 + h * h)
    return [
        (np.sqrt(3) * old[0] - h * old[1]) / d,
        (h * old[0] + np.sqrt(3) * old[1]) / d,
        *old[2:],
    ]


def pluecker(point: np.ndarray, basis: list[np.ndarray]) -> np.ndarray:
    return np.array(
        [2 * point @ np.cross(basis[i] @ point, basis[j] @ point) for i, j in PAIRS]
    )


def contraction(G: np.ndarray) -> np.ndarray:
    out = np.zeros((5, 5))
    edge = {pair: index for index, pair in enumerate(PAIRS)}
    for i in range(5):
        for k in range(5):
            for j in range(5):
                if j in (i, k):
                    continue
                pi, si = ((i, j), 1) if i < j else ((j, i), -1)
                pk, sk = ((k, j), 1) if k < j else ((j, k), -1)
                out[i, k] += si * sk * G[edge[pi], edge[pk]]
    return (out + out.T) / 2


HODGE = np.zeros((6, 6))
for i, j, sign in ((0, 5, 1), (1, 4, -1), (2, 3, 1)):
    HODGE[i, j] = HODGE[j, i] = sign


def outer_contraction(B: np.ndarray) -> np.ndarray:
    out = np.zeros((4, 4))
    pairs = list(itertools.combinations(range(4), 2))
    edge = {pair: index for index, pair in enumerate(pairs)}
    for i in range(4):
        for k in range(4):
            for j in range(4):
                if j in (i, k):
                    continue
                pi, si = ((i, j), 1) if i < j else ((j, i), -1)
                pk, sk = ((k, j), 1) if k < j else ((j, k), -1)
                out[i, k] += si * sk * B[edge[pi], edge[pk]]
    return (out + out.T) / 2


def tangent_slice(h: float, samples: int = 300, seed: int = 82026):
    rng = np.random.default_rng(seed)
    basis = orbitals(h)
    rows = []
    for _ in range(samples):
        point = rng.normal(size=3)
        point /= np.linalg.norm(point)
        z = pluecker(point, basis)
        rows.append(np.outer(z, z).reshape(-1))
    _, singular, vh = np.linalg.svd(np.stack(rows), full_matrices=False)
    rank = int(np.sum(singular > 1e-10))
    assert rank == 28
    raw = [row.reshape(10, 10) for row in vh[:rank]]
    equations = np.array([[contraction(G)[0, i] for G in raw] for i in range(1, 5)])
    kernel = null_space(equations)
    sliced = [sum(kernel[i, j] * raw[i] for i in range(rank)) for j in range(kernel.shape[1])]
    assert len(sliced) == 24
    return sliced


def target_matrix(h: float, adjusted: bool = False):
    basis = tangent_slice(h)
    us, cs, bs, masses, z2s = [], [], [], [], []
    d = np.sqrt(3 + h * h)
    w = np.zeros(10)
    w[1] = 2 * np.sqrt(3) / d
    w[4] = 2 * h / d
    w[9] = 1
    for G in basis:
        A, C, B = G[:4, :4], G[:4, 4:], G[4:, 4:]
        F = contraction(G)
        mass = np.trace(G)
        delta = 3 * F[0, 0] - 2 * mass
        U = A - outer_contraction(B) - delta * np.eye(4) / 4
        us.append(U.reshape(-1))
        cs.append(C.reshape(-1))
        bs.append(B)
        masses.append(mass)
        z2s.append(w @ G @ w)
    us = np.stack(us)
    cs = np.stack(cs)
    masses = np.array(masses)
    z2s = np.array(z2s)
    q0 = 3 / (3 + h * h)
    zcoefficient = 1 - h * h
    if adjusted:
        zcoefficient = 1 - 2 * h * h / 3
    r = 3 * q0 * (masses - zcoefficient * z2s) - 2 * masses
    matrix = 2 * (us @ us.T) + 8 * (cs @ cs.T) - np.outer(r, r) / 6
    for i, Bi in enumerate(bs):
        for j, Bj in enumerate(bs):
            matrix[i, j] += 4 * np.trace(Bi @ HODGE @ Bj @ HODGE)
    return (matrix + matrix.T) / 2


if __name__ == "__main__":
    for h in np.linspace(-1, 1, 17):
        original = np.linalg.eigvalsh(target_matrix(float(h)))
        adjusted = np.linalg.eigvalsh(target_matrix(float(h), adjusted=True))
        print(
            f"h={h:+.3f}",
            "orig min/rank", original[0], np.sum(original > 1e-8),
            "adj min/rank", adjusted[0], np.sum(adjusted > 1e-8),
        )
