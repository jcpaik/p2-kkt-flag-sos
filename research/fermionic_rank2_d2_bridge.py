"""Audit the exact D2-sector bridge for a rank-two top orbital.

For S=diag(1,-1,0)/sqrt(2), decompose a tangent Pluecker moment G into
the four characters of the stabilizer D2 of S.  The nontrivial character
sectors obey

    Q_chi + 3/22 ||K_chi||_F^2 >= 0.

The invariant-sector certificate in ``tensor_star_canonical_exact.py``
has an explicit nonnegative product remainder L.  This script tests the
candidate final bridge

    L >= 3/22 ||K_nontrivial||_F^2

on genuine random measures and reports exact ingredients for any failure.
"""

from __future__ import annotations

import itertools

import numpy as np
import sympy as sp

import tensor_star_canonical_exact as canonical
import tensor_star_rank2_d2_exact as d2


EDGES = canonical.PAIRS


def pluecker(point: np.ndarray) -> np.ndarray:
    x = np.asarray(point, dtype=float)
    orbitals = [np.array(matrix, dtype=float) for matrix in canonical.ORB]
    return np.array(
        [2 * x @ np.cross(orbitals[i] @ x, orbitals[j] @ x) for i, j in EDGES]
    )


def contraction(G: np.ndarray) -> np.ndarray:
    out = np.zeros((5, 5))
    edge_index = {edge: index for index, edge in enumerate(EDGES)}
    for i in range(5):
        for k in range(5):
            for j in range(5):
                if j in (i, k):
                    continue
                left, sl = ((i, j), 1) if i < j else ((j, i), -1)
                right, sr = ((k, j), 1) if k < j else ((j, k), -1)
                out[i, k] += sl * sr * G[edge_index[left], edge_index[right]]
    return (out + out.T) / 2


def invariant_product_remainder(G: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Return L, diagonal top slack, and its off-diagonal part."""
    p = np.diag(G)
    F = contraction(G)
    c = F[0, 0]
    caps = c - np.diag(F)[1:]

    def edge(i: int, j: int) -> float:
        return p[EDGES.index((min(i, j), max(i, j)))]

    L = (
        4 * edge(1, 2) * edge(3, 4) / 3
        + edge(1, 2) * caps[0] / 3
        + (edge(2, 3) * caps[2] + edge(2, 4) * caps[3]) / 3
        + (edge(2, 3) * caps[3] + edge(2, 4) * caps[2]) / 12
    )
    K = c * np.eye(4) - F[1:, 1:]
    Koff = K - np.diag(np.diag(K))
    return L, caps, Koff


def random_audit(samples: int = 20000, atoms: int = 8, seed: int = 82026):
    rng = np.random.default_rng(seed)
    worst = None
    for _ in range(samples):
        points = rng.normal(size=(atoms, 3))
        points /= np.linalg.norm(points, axis=1)[:, None]
        weights = rng.dirichlet(np.ones(atoms))
        vectors = np.stack([pluecker(point) for point in points])
        G = np.einsum("i,ij,ik->jk", weights, vectors, vectors)
        F = contraction(G)
        if np.linalg.norm(F[0, 1:]) > 1e-8:
            # A random fixed chart does not make S an eigenvector.  Rotate
            # the measure is not enough to turn an arbitrary top orbital
            # into the rank-two orbit, so this routine is only a raw stress
            # test of the proposed algebraic bridge.
            continue
        L, caps, Koff = invariant_product_remainder(G)
        gap = L - 3 * np.sum(Koff * Koff) / 22
        if worst is None or gap < worst[0]:
            worst = (gap, L, np.sum(Koff * Koff), caps, G, F)
    return worst


def exact_sector_penalties():
    """Eliminate each five-dimensional character sector at fixed K off block.

    Returns the exact two-by-two quadratic H such that

        Q_chi >= - k_chi.T H k_chi,

    where k_chi consists of the two independent upper-triangular entries of
    the 4x4 top-slack matrix K in that character.  This is the sharp fiber
    minimum, rather than the coarser 3/22 Frobenius-norm completion.
    """
    gm, fm, km = d2.exact_rank_two_slice()
    q = d2.purity_form(gm, fm)
    out = {}
    for character in d2.NONTRIVIAL:
        basis = d2.sector_basis(gm, character)
        qs = sp.simplify(basis.T * q * basis)
        physical_k = sp.simplify(km * basis)
        entries = []
        for i in range(4):
            for j in range(i + 1, 4):
                row = physical_k[4 * i + j, :]
                if row != sp.zeros(1, 5):
                    entries.append(row)
        assert len(entries) == 2
        L = sp.Matrix.vstack(*entries)

        # Coordinates x=N y + R k, where N spans ker(L) and L R=I.
        N = sp.Matrix.hstack(*L.nullspace())
        R = sp.simplify(L.T * (L * L.T).inv())
        assert L * R == sp.eye(2)
        A = sp.simplify(N.T * qs * N)
        B = sp.simplify(N.T * qs * R)
        C = sp.simplify(R.T * qs * R)
        assert A.is_positive_definite
        effective = sp.simplify(C - B.T * A.inv() * B)
        H = -effective
        out[character] = (L, H)
    return out


def invariant_quadratic_data():
    """Return exact Q0, cap, mass, and z^2-moment maps in nine weights."""
    _, forced = canonical.exact_least_norm_form()
    indices = [0, 2, 3, 4, 5, 6, 7, 8, 9]
    substitution = sp.zeros(10, 9)
    for column, index in enumerate(indices):
        substitution[index, column] = 1
    substitution[1, 0] = substitution[1, 3] = sp.Rational(1, 3)
    q0 = sp.simplify(substitution.T * (forced + canonical.graph_form()) * substitution)
    center = sp.Matrix([[int(i == 0) for i, _ in EDGES]]) * substitution
    caps = sp.zeros(4, 9)
    for orbital in range(1, 5):
        degree = sp.Matrix(
            [[int(i == orbital or j == orbital) for i, j in EDGES]]
        ) * substitution
        caps[orbital - 1, :] = center - degree
    mass = sp.ones(1, 10) * substitution

    x, y, z = canonical.R
    coefficients = sp.symbols("ell0:10")
    polynomial = sp.expand(
        sum(coefficients[i] * canonical.PLUECKER[i] ** 2 for i in range(10))
        - z**2 * (x**2 + y**2 + z**2) ** 2
    )
    equations = [coefficient for coefficient in sp.Poly(polynomial, x, y, z).coeffs()]
    solution = next(iter(sp.linsolve(equations, coefficients)))
    # Pick the zero value for any free diagonal-Hankel gauge parameters.
    free = sorted(set().union(*(value.free_symbols for value in solution)), key=str)
    gauge = {symbol: 0 for symbol in free if symbol not in set(coefficients)}
    # linsolve normally parametrizes with some of the original ell symbols.
    for symbol in coefficients:
        if symbol in set().union(*(value.free_symbols for value in solution)):
            gauge[symbol] = 0
    z2_diagonal = sp.Matrix([[sp.simplify(value.subs(gauge)) for value in solution]])
    m = sp.simplify(z2_diagonal * substitution)
    assert sp.expand(
        sum(z2_diagonal[i] * canonical.PLUECKER[i] ** 2 for i in range(10))
        - z**2 * (x**2 + y**2 + z**2) ** 2
    ) == 0
    return q0, caps, mass, m


def bad_k_form(K: np.ndarray) -> float:
    """Sharp sum of the three eliminated nontrivial-sector defects."""
    k12, k13, k14 = K[0, 1], K[0, 2], K[0, 3]
    k23, k24, k34 = K[1, 2], K[1, 3], K[2, 3]
    return (
        -5 * k12**2
        - 8 * np.sqrt(3) * k12 * k34
        - 3 * k34**2
        + k13**2
        + 4 * np.sqrt(3) * k13 * k24
        - 3 * k24**2
        + k14**2
        + 4 * np.sqrt(3) * k14 * k23
        - 3 * k23**2
    ) / 11


def abstract_k_audit(samples: int = 200000, seed: int = 20260820):
    """Stress-test rank-two T0 >= 4H(K) in the generic PSD-K relaxation."""
    q0_exact, caps_exact, mass_exact, m_exact = invariant_quadratic_data()
    q0 = np.array(q0_exact, dtype=float)
    cap_map = np.array(caps_exact, dtype=float)
    mass_map = np.array(mass_exact, dtype=float).reshape(-1)
    m_map = np.array(m_exact, dtype=float).reshape(-1)
    rng = np.random.default_rng(seed)
    worst = None
    accepted = 0
    for _ in range(samples):
        p = rng.exponential(size=9)
        caps = cap_map @ p
        if np.min(caps) < 0:
            continue
        accepted += 1
        raw = rng.normal(size=(4, 4))
        gram = raw @ raw.T
        scale = np.sqrt(np.diag(gram))
        correlation = gram / np.outer(scale, scale)
        K = np.sqrt(caps[:, None] * caps[None, :]) * correlation
        mass = mass_map @ p
        m = m_map @ p
        c = mass - np.sum(0.0 * p) - 0.0  # documented below via alpha
        # alpha=t-m-c, while c is the star occupation F_00.  In the
        # invariant diagonal chart it is the sum of the first four edges.
        full_p = np.array(
            [p[0], (p[0] + p[3]) / 3, p[1], p[2], *p[3:]], dtype=float
        )
        center = np.sum(full_p[:4])
        alpha = mass - m - center
        invariant_t = (
            4 * (p @ q0 @ p)
            - (mass - 3 * m) * alpha
            + 1.5 * alpha**2
        )
        gap = invariant_t - 4 * bad_k_form(K)
        normalized = gap / max(np.sum(p) ** 2, 1e-30)
        if worst is None or normalized < worst[0]:
            worst = (
                normalized,
                gap,
                p,
                caps,
                K,
                invariant_t,
                4 * bad_k_form(K),
                mass,
                m,
                alpha,
            )
    return accepted, worst


if __name__ == "__main__":
    for character, (linear_map, penalty) in exact_sector_penalties().items():
        print("character", character)
        print("K entry map =", linear_map)
        print("sharp penalty H =", penalty)
    accepted, worst = abstract_k_audit(samples=50000)
    print("abstract accepted", accepted)
    print("abstract worst", worst)
