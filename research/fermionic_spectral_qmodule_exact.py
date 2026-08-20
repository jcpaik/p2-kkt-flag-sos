"""Exact spectral q-module certificate for the unrestricted tangent bridge.

Let S be a unit traceless symmetric matrix, diagonalized and parametrized by

    S=(sqrt(3) E0-h E1)/sqrt(3+h^2),     0 <= h <= 1.

The D2 stabilizer of this diagonal chart acts by coordinate sign changes.
After D2 averaging, the entries of G are linear in the ten cubic moments of
X=x^2, Y=y^2, Z=z^2.  This script derives exact weights w_X,w_Y,w_Z for

    ell=(3-h)X-(3+h)Y+2hZ,
    g_h=(w_X X+w_Y Y+w_Z Z) ell^2 >= 0,

such that r_h=delta+E[g_h] belongs to the linear span of the U and C
entries.  It then computes the least weighted dual norm D(h).  If
w_X,w_Y,w_Z>=0 and D(h)<=6 on [0,1], weighted Cauchy gives

    2||U||^2+8||C||^2 >= r_h^2/D(h) >= delta^2/6

whenever delta>0.  Jensen under the D2 average and the nonnegative Hodge
term then prove the full unrestricted fermionic estimate.
"""

from __future__ import annotations

import itertools

import sympy as sp


h = sp.symbols("h", real=True)
x, y, z = sp.symbols("x y z", real=True)
X, Y, Z = sp.symbols("X Y Z", nonnegative=True)
MONOMIALS = [(a, b, 3 - a - b) for a in range(4) for b in range(4 - a)]
PAIRS = list(itertools.combinations(range(5), 2))


def orbitals():
    E0 = sp.diag(1, -1, 0) / sp.sqrt(2)
    E1 = sp.diag(1, 1, -2) / sp.sqrt(6)
    S = (sp.sqrt(3) * E0 - h * E1) / sp.sqrt(3 + h**2)
    H = (h * E0 + sp.sqrt(3) * E1) / sp.sqrt(3 + h**2)
    out = [S, H]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        matrix = sp.zeros(3)
        matrix[i, j] = matrix[j, i] = 1 / sp.sqrt(2)
        out.append(matrix)
    return out


def pluecker():
    point = sp.Matrix([x, y, z])
    basis = orbitals()
    return [
        sp.expand(2 * point.dot((basis[i] * point).cross(basis[j] * point)))
        for i, j in PAIRS
    ]


def even_coefficient_row(polynomial: sp.Expr) -> sp.Matrix:
    poly = sp.Poly(sp.expand(polynomial), x, y, z)
    return sp.Matrix(
        [[
            poly.coeff_monomial(x ** (2 * a) * y ** (2 * b) * z ** (2 * c))
            for a, b, c in MONOMIALS
        ]]
    )


def contraction_rows(G: list[list[sp.Matrix]]) -> list[list[sp.Matrix]]:
    out = [[sp.zeros(1, 10) for _ in range(5)] for _ in range(5)]
    for a, (i, j) in enumerate(PAIRS):
        for b, (k, ell) in enumerate(PAIRS):
            row = G[a][b]
            if j == ell:
                out[i][k] += row
            if j == k:
                out[i][ell] -= row
            if i == ell:
                out[j][k] -= row
            if i == k:
                out[j][ell] += row
    return out


def outer_contraction_rows(B: list[list[sp.Matrix]]) -> list[list[sp.Matrix]]:
    pairs = list(itertools.combinations(range(4), 2))
    out = [[sp.zeros(1, 10) for _ in range(4)] for _ in range(4)]
    for a, (i, j) in enumerate(pairs):
        for b, (k, ell) in enumerate(pairs):
            row = B[a][b]
            if j == ell:
                out[i][k] += row
            if j == k:
                out[i][ell] -= row
            if i == ell:
                out[j][k] -= row
            if i == k:
                out[j][ell] += row
    return out


def maps():
    zeta = pluecker()
    G = [[even_coefficient_row(zeta[i] * zeta[j]) for j in range(10)] for i in range(10)]
    F = contraction_rows(G)
    mass = sum((G[i][i] for i in range(10)), sp.zeros(1, 10))
    top = sum((G[i][i] for i in range(4)), sp.zeros(1, 10))
    delta = sp.simplify(3 * top - 2 * mass)
    B = [[G[i + 4][j + 4] for j in range(6)] for i in range(6)]
    R = outer_contraction_rows(B)
    U = []
    for i in range(4):
        for j in range(4):
            row = G[i][j] - R[i][j]
            if i == j:
                row -= delta / 4
            if row != sp.zeros(1, 10):
                U.append(sp.simplify(row))
    C = []
    for i in range(4):
        for j in range(6):
            row = G[i][j + 4]
            if row != sp.zeros(1, 10):
                C.append(sp.simplify(row))
    return G, F, mass, delta, U, C


def polynomial_row(polynomial: sp.Expr) -> sp.Matrix:
    poly = sp.Poly(sp.expand(polynomial), X, Y, Z)
    return sp.Matrix(
        [[poly.coeff_monomial(X**a * Y**b * Z**c) for a, b, c in MONOMIALS]]
    )


def derive(verbose: bool = True):
    _, _, mass, delta, U, C = maps()
    feature = sp.Matrix.vstack(*U, *C)
    weights = sp.diag(*([2] * len(U) + [8] * len(C)))
    norm = sp.simplify(feature.T * weights * feature)
    # Row-normalization removes harmless sqrt(3+h^2) factors and makes the
    # rational-function nullspace calculation dramatically smaller.
    normalized_rows = []
    for index in range(feature.rows):
        row = feature.row(index)
        first = next(value for value in row if value != 0)
        normalized_rows.append(
            sp.Matrix([[sp.cancel(value / first) for value in row]])
        )
    range_feature = sp.Matrix.vstack(*normalized_rows)

    ell = (3 - h) * X - (3 + h) * Y + 2 * h * Z
    grows = [polynomial_row(variable * ell**2) for variable in (X, Y, Z)]
    wx, wy, wz = sp.symbols("w_X w_Y w_Z", real=True)
    r = sp.simplify(delta + wx * grows[0] + wy * grows[1] + wz * grows[2])

    null = range_feature.nullspace(simplify=False)
    assert len(null) == 2
    range_equations = [sp.factor((r * vector)[0]) for vector in null]
    range_solution = next(iter(sp.linsolve(range_equations, (wx, wy))))
    range_substitution = {wx: range_solution[0], wy: range_solution[1]}
    r = sp.simplify(r.subs(range_substitution))

    dual_variables = sp.symbols("d0:10")
    dual = sp.Matrix(dual_variables)
    solution = next(iter(sp.linsolve(list(norm * dual - r.T), dual_variables)))
    dual_substitution = dict(zip(dual_variables, solution))
    dual = sp.Matrix([dual_substitution[value] for value in dual_variables])
    D = sp.factor((r * dual)[0])
    # D is independent of the two nullspace gauges in the solution.
    free_weights = sorted(D.free_symbols & {wx, wy, wz}, key=str)
    gauges = sorted(D.free_symbols - {h, *free_weights}, key=str)
    assert not gauges, gauges
    optimum = next(
        iter(sp.linsolve([sp.diff(D, variable) for variable in free_weights], free_weights))
    )
    optimal_substitution = dict(zip(free_weights, optimum))
    ranged_weights = [
        sp.simplify(variable.subs(range_substitution)) for variable in (wx, wy, wz)
    ]
    final_weights = [
        sp.factor(value.subs(optimal_substitution)) for value in ranged_weights
    ]
    final_D = sp.factor(D.subs(optimal_substitution))
    final_r = sp.simplify(r.subs(optimal_substitution))

    if verbose:
        print("feature counts U/C", len(U), len(C), "rank", feature.rank())
        print("range solution wx,wy =", [sp.factor(value) for value in range_solution])
        print("optimal weights =", final_weights)
        print("D =", final_D)
        print("6-D =", sp.factor(6 - final_D))
        print("r row =", final_r)
    return {
        "mass": mass,
        "delta": delta,
        "U": U,
        "C": C,
        "feature": feature,
        "norm": norm,
        "ell": ell,
        "grows": grows,
        "weights": final_weights,
        "r": final_r,
        "D": final_D,
    }


if __name__ == "__main__":
    derive()
