"""Random audit of the D2-twirled direct square bound for general spectrum."""

from __future__ import annotations

import itertools

import numpy as np

import fermionic_general_spectrum_scalar_audit as general


SIGNS = np.array(((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)))


def d2_atom(squares: np.ndarray, basis: list[np.ndarray]) -> np.ndarray:
    point = np.sqrt(squares)
    out = np.zeros((10, 10))
    for signs in SIGNS:
        zeta = general.pluecker(signs * point, basis)
        out += np.outer(zeta, zeta) / 4
    return out


def outer_contraction(B: np.ndarray) -> np.ndarray:
    out = np.zeros((4, 4))
    edges = list(itertools.combinations(range(4), 2))
    for a, (i, j) in enumerate(edges):
        for b, (k, ell) in enumerate(edges):
            value = B[a, b]
            if j == ell:
                out[i, k] += value
            if j == k:
                out[i, ell] -= value
            if i == ell:
                out[j, k] -= value
            if i == k:
                out[j, ell] += value
    return out


def direct_data(G: np.ndarray):
    F = general.contraction(G)
    A, C, B = G[:4, :4], G[:4, 4:], G[4:, 4:]
    delta = 3 * F[0, 0] - 2
    U = A - outer_contraction(B) - delta * np.eye(4) / 4
    core = 2 * np.sum(U * U) + 8 * np.sum(C * C)
    hodge = np.trace(B @ general.HODGE @ B @ general.HODGE)
    return F, delta, core, hodge


def search(h: float, trials: int = 300_000, seed: int = 20260820):
    rng = np.random.default_rng(seed)
    basis = general.orbitals(h)
    best = (np.inf, None)
    accepted = 0
    for _ in range(trials):
        number = rng.integers(1, 7)
        locations = rng.dirichlet(np.ones(3), size=number)
        weights = rng.dirichlet(np.ones(number))
        G = sum(
            weight * d2_atom(location, basis)
            for weight, location in zip(weights, locations)
        )
        F, delta, core, hodge = direct_data(G)
        if delta <= 0 or F[0, 0] + 1e-11 < np.max(np.diag(F)):
            continue
        accepted += 1
        gap = core - delta * delta / 6
        if gap < best[0]:
            best = (gap, (locations, weights, F, delta, core, hodge))
    return accepted, best


if __name__ == "__main__":
    for h in np.linspace(0, 1, 11):
        accepted, best = search(float(h), trials=20_000, seed=20260820)
        print(f"h={h:.2f} accepted={accepted} best_core_gap={best[0]:+.10g}")
        if best[1] is not None and best[0] < -1e-8:
            print("  F", np.diag(best[1][2]), "delta", best[1][3], "hodge", best[1][5])
