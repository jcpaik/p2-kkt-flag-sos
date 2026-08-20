"""Exactify the all-spectrum D2 correction certificate.

For

    S_h=(sqrt(3) S0-h H0)/sqrt(3+h^2),  0 <= h <= 1,

restrict tangent Pluecker moments to the D2-invariant even sextic sector,
with coordinates X=x^2, Y=y^2, Z=z^2.  This file constructs exactly the
quadratic core

    2 ||U||^2 + 8 ||C||^2

and searches symbolically for weights w_X,w_Y,w_Z such that

    r = delta + E[(w_X X+w_Y Y+w_Z Z) ell_h^2]

lies in the range of the core, where

    ell_h=(3-h)X-(3+h)Y+2hZ.

If the weights are nonnegative and the minimum dual cost D is at most 6,
then the pointwise correction is nonnegative and weighted Cauchy proves
the desired direct bridge core >= delta^2/6 after D2 twirling.
"""

from __future__ import annotations

import itertools

import sympy as sp


h = sp.symbols("h", real=True)
X, Y, Z = sp.symbols("X Y Z", nonnegative=True)
x, y, z = sp.symbols("x y z", real=True)
MONOMIALS = [(a, b, 3 - a - b) for a in range(4) for b in range(4 - a)]
PAIRS = list(itertools.combinations(range(5), 2))


def coefficient_row(polynomial: sp.Expr) -> sp.Matrix:
    poly = sp.Poly(sp.expand(polynomial), x, y, z)
    return sp.Matrix([[
        poly.coeff_monomial(x ** (2 * a) * y ** (2 * b) * z ** (2 * c))
        for a, b, c in MONOMIALS
    ]])


def symbolic_maps():
    s0 = sp.diag(1, -1, 0) / sp.sqrt(2)
    h0 = sp.diag(1, 1, -2) / sp.sqrt(6)
    denominator = sp.sqrt(3 + h**2)
    sh = (sp.sqrt(3) * s0 - h * h0) / denominator
    hh = (h * s0 + sp.sqrt(3) * h0) / denominator
    orbitals = [sh, hh]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        matrix = sp.zeros(3)
        matrix[i, j] = matrix[j, i] = 1 / sp.sqrt(2)
        orbitals.append(matrix)

    point = sp.Matrix((x, y, z))
    pluecker = [
        sp.expand(2 * point.dot((orbitals[i] * point).cross(orbitals[j] * point)))
        for i, j in PAIRS
    ]
    gm = sp.zeros(100, 10)
    for i in range(10):
        for j in range(10):
            gm[10 * i + j, :] = coefficient_row(pluecker[i] * pluecker[j])

    def grow(i, j):
        return gm[10 * i + j, :]

    tr_a = sum((grow(i, i) for i in range(4)), sp.zeros(1, 10))
    mass = sum((grow(i, i) for i in range(10)), sp.zeros(1, 10))
    delta = sp.simplify(3 * tr_a - 2 * mass)

    outer_pairs = [(i - 1, j - 1) for i, j in PAIRS[4:]]
    contracted = sp.zeros(16, 10)
    for a, (i, j) in enumerate(outer_pairs):
        for b, (k, ell) in enumerate(outer_pairs):
            value = grow(a + 4, b + 4)
            if j == ell:
                contracted[4 * i + k, :] += value
            if j == k:
                contracted[4 * i + ell, :] -= value
            if i == ell:
                contracted[4 * j + k, :] -= value
            if i == k:
                contracted[4 * j + ell, :] += value

    urows = []
    for i in range(4):
        for j in range(4):
            urows.append(sp.simplify(
                grow(i, j)
                - contracted[4 * i + j, :]
                - (delta / 4 if i == j else sp.zeros(1, 10))
            ))
    crows = [grow(i, j + 4) for i in range(4) for j in range(6)]
    core = sp.zeros(10)
    for row in urows:
        core += 2 * row.T * row
    for row in crows:
        core += 8 * row.T * row
    return core, delta, urows, crows


def xyz_coefficient_row(polynomial: sp.Expr) -> sp.Matrix:
    poly = sp.Poly(sp.expand(polynomial), X, Y, Z)
    return sp.Matrix([[
        poly.coeff_monomial(X**a * Y**b * Z**c)
        for a, b, c in MONOMIALS
    ]])


def correction_row(wx, wy, wz):
    ell = (3 - h) * X - (3 + h) * Y + 2 * h * Z
    return xyz_coefficient_row((wx * X + wy * Y + wz * Z) * ell**2)


def derive(verbose: bool = True):
    core, delta, urows, crows = symbolic_maps()
    print("constructed core") if verbose else None
    raw_rows = [row for row in (*urows, *crows) if any(value != 0 for value in row)]
    normalized_rows = []
    seen = set()
    for row in raw_rows:
        pivot = next(value for value in row if value != 0)
        normalized = sp.Matrix([[
            sp.factor(sp.cancel(value / pivot)) for value in row
        ]])
        key = tuple(normalized)
        if key not in seen:
            normalized_rows.append(normalized)
            seen.add(key)
    feature_map = sp.Matrix.vstack(*normalized_rows)
    print("nonzero/normalized feature rows", len(raw_rows), len(normalized_rows)) if verbose else None
    _, pivots = feature_map.rref(simplify=False)
    print("generic rank", len(pivots)) if verbose else None
    null = feature_map.nullspace(simplify=False)
    print("nullity", len(null)) if verbose else None
    wx, wy, wz = sp.symbols("w_X w_Y w_Z", real=True)
    corrected = sp.simplify(delta + correction_row(wx, wy, wz))
    constraints = [sp.factor((vector.T * corrected.T)[0]) for vector in null]
    print("range constraints", constraints) if verbose else None
    solution = sp.solve(constraints, (wx, wy), dict=True, simplify=True)
    print("weight solutions", solution) if verbose else None
    assert constraints == [0, 0]

    # Every one of delta, X ell^2, Y ell^2, Z ell^2 is in the feature
    # range.  Compute their exact dual Gram matrix on any independent set
    # of eight moment coordinates.
    source_rows = [delta]
    source_rows.extend(
        correction_row(*weights)
        for weights in ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    )
    source = sp.Matrix.vstack(*[row[:, list(pivots)] for row in source_rows]).T
    core_sub = core.extract(pivots, pivots)
    solved = core_sub.inv(method="DM") * source
    dual_gram = sp.simplify(source.T * solved)
    print("dual Gram", dual_gram.applyfunc(sp.factor)) if verbose else None
    wvector = sp.Matrix((wx, wy, wz))
    optimum = sp.simplify(-dual_gram[1:, 1:].inv(method="DM") * dual_gram[1:, 0])
    optimum = optimum.applyfunc(lambda value: sp.factor(sp.cancel(value)))
    optimum_cost = sp.factor(sp.cancel(
        dual_gram[0, 0]
        + 2 * (dual_gram[0, 1:] * optimum)[0]
        + (optimum.T * dual_gram[1:, 1:] * optimum)[0]
    ))
    print("optimal weights", optimum) if verbose else None
    print("optimal dual cost", optimum_cost) if verbose else None
    return {
        "core": core,
        "delta": delta,
        "corrected": corrected,
        "null": null,
        "constraints": constraints,
        "pivots": pivots,
        "dual_gram": dual_gram,
        "optimal_weights": optimum,
        "optimal_cost": optimum_cost,
    }


if __name__ == "__main__":
    derive()
