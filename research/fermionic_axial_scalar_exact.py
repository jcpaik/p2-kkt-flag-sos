"""Exact axial (double-eigenvalue) tangent-star scalar audit.

Fix the unit axial orbital

    S = diag(1,-2,1)/sqrt(6)

and use the orthonormal completion

    H = diag(1,0,-1)/sqrt(2), xy, xz, yz.

For a tangent Pluecker moment matrix G, impose F S = c S, where F is
the one-particle contraction.  This script constructs the exact
24-dimensional linear slice and tests the desired endpoint form

    R0 - tr(G)^2/96,

where R0 = 2||U||^2 + 8||C||^2 + 4 tr(B *B*).

The file is a discovery/verifier aid.  Its printed inertia determines
whether the endpoint follows from an ordinary quadratic SOS on the
linear tangent/eigen slice, before using G >= 0.
"""

from __future__ import annotations

import itertools

import sympy as sp


rt = sp.sqrt
RVAR = sp.Matrix(sp.symbols("x y z"))


def orbital_basis() -> list[sp.Matrix]:
    out = [
        sp.diag(1, -2, 1) / rt(6),
        sp.diag(1, 0, -1) / rt(2),
    ]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        matrix = sp.zeros(3)
        matrix[i, j] = matrix[j, i] = 1 / rt(2)
        out.append(matrix)
    gram = sp.Matrix([[sp.trace(a * b) for b in out] for a in out])
    assert gram == sp.eye(5)
    return out


ORB = orbital_basis()
EDGES = list(itertools.combinations(range(5), 2))
PLUECKER = [
    sp.expand(2 * RVAR.dot((ORB[i] * RVAR).cross(ORB[j] * RVAR)))
    for i, j in EDGES
]


def tangent_span() -> tuple[list[tuple[int, int]], sp.Matrix]:
    monomials = [
        (a, b, 6 - a - b) for a in range(7) for b in range(7 - a)
    ]
    coordinates = [(a, b) for a in range(10) for b in range(a, 10)]
    span = sp.zeros(55, 28)
    x, y, z = RVAR
    for row, (a, b) in enumerate(coordinates):
        polynomial = sp.Poly(
            sp.expand(PLUECKER[a] * PLUECKER[b]
                      * (rt(2) if a < b else 1)),
            x, y, z,
        )
        for column, (i, j, k) in enumerate(monomials):
            span[row, column] = polynomial.coeff_monomial(x**i * y**j * z**k)
    assert span.rank() == 28
    return coordinates, span


def contraction_matrix() -> sp.Matrix:
    out = sp.zeros(25, 100)
    for a, (i, j) in enumerate(EDGES):
        for b, (k, ell) in enumerate(EDGES):
            column = 10 * a + b
            if j == ell:
                out[5 * i + k, column] += 1
            if j == k:
                out[5 * i + ell, column] -= 1
            if i == ell:
                out[5 * j + k, column] -= 1
            if i == k:
                out[5 * j + ell, column] += 1
    return out


def exact_axial_slice() -> tuple[sp.Matrix, sp.Matrix]:
    coordinates, span = tangent_span()
    gmap = sp.zeros(100, 28)
    for row, (a, b) in enumerate(coordinates):
        value = span[row, :] / (rt(2) if a < b else 1)
        gmap[10 * a + b, :] = value
        gmap[10 * b + a, :] = value
    fmap = contraction_matrix() * gmap
    eigen_constraint = fmap[[5 * i for i in range(1, 5)], :]
    null = sp.Matrix.hstack(*eigen_constraint.nullspace())
    assert null.shape == (28, 24)
    return sp.simplify(gmap * null), sp.simplify(fmap * null)


def outer_contraction_map(gm: sp.Matrix) -> sp.Matrix:
    """Map the outer B block to its 4x4 contraction gamma(B)."""
    out = sp.zeros(16, gm.cols)
    outer_edges = [(i - 1, j - 1) for i, j in EDGES[4:]]
    for a, (i, j) in enumerate(outer_edges):
        for b, (k, ell) in enumerate(outer_edges):
            value = gm[10 * (a + 4) + (b + 4), :]
            if j == ell:
                out[4 * i + k, :] += value
            if j == k:
                out[4 * i + ell, :] -= value
            if i == ell:
                out[4 * j + k, :] -= value
            if i == k:
                out[4 * j + ell, :] += value
    return out


def hodge_form(gm: sp.Matrix) -> sp.Matrix:
    hodge = sp.zeros(6)
    for i, j, sign in ((0, 5, 1), (1, 4, -1), (2, 3, 1)):
        hodge[i, j] = hodge[j, i] = sign
    rows = [10 * (a + 4) + (b + 4) for a in range(6) for b in range(6)]
    bmap = gm[rows, :]
    return sp.simplify(
        bmap.T * sp.kronecker_product(hodge, hodge) * bmap
    )


def endpoint_form():
    gm, fm = exact_axial_slice()
    amap = gm[[10 * a + b for a in range(4) for b in range(4)], :]
    cmap = gm[[10 * a + (b + 4) for a in range(4) for b in range(6)], :]
    rmap = outer_contraction_map(gm)
    tr_a = sum((gm[11 * a, :] for a in range(4)), sp.zeros(1, gm.cols))
    tr_b = sum((gm[11 * a, :] for a in range(4, 10)), sp.zeros(1, gm.cols))
    delta = tr_a - 2 * tr_b
    imap = sp.zeros(16, gm.cols)
    for i in range(4):
        imap[5 * i, :] = delta / 4
    umap = sp.simplify(amap - rmap - imap)
    mass = tr_a + tr_b
    r0 = sp.simplify(
        2 * umap.T * umap
        + 8 * cmap.T * cmap
        + 4 * hodge_form(gm)
    )
    endpoint = sp.simplify(r0 - mass.T * mass / 96)
    return gm, fm, mass, r0, endpoint


def circle_average_polynomial(polynomial: sp.Expr, t: sp.Symbol) -> sp.Expr:
    """Average a homogeneous sextic over y^2=t and the xz circle."""
    x, y, z = RVAR
    poly = sp.Poly(sp.expand(polynomial), x, y, z)
    out = 0
    for (px, py, pz), coefficient in poly.terms():
        if px % 2 or py % 2 or pz % 2:
            continue
        a, b, c = px // 2, py // 2, pz // 2
        circle = (
            sp.factorial(2 * a) * sp.factorial(2 * c)
            / (
                4 ** (a + c)
                * sp.factorial(a)
                * sp.factorial(c)
                * sp.factorial(a + c)
            )
        )
        out += coefficient * t**b * (1 - t) ** (a + c) * circle
    return sp.factor(out)


def axial_average_projector(t: sp.Symbol) -> sp.Matrix:
    zeta = sp.Matrix(PLUECKER)
    raw = zeta * zeta.T
    return raw.applyfunc(lambda value: circle_average_polynomial(value, t))


def split_value(matrix: sp.Matrix):
    """Return F, delta, U, C, B, Hodge and R0 for an actual 10x10 G."""
    flat = sp.Matrix(100, 1, list(matrix))
    f = sp.Matrix(5, 5, contraction_matrix() * flat)
    a = matrix[:4, :4]
    b = matrix[4:, 4:]
    cblock = matrix[:4, 4:]
    outer_edges = [(i - 1, j - 1) for i, j in EDGES[4:]]
    contracted = sp.zeros(4)
    for edge_a, (i, j) in enumerate(outer_edges):
        for edge_b, (k, ell) in enumerate(outer_edges):
            value = b[edge_a, edge_b]
            if j == ell:
                contracted[i, k] += value
            if j == k:
                contracted[i, ell] -= value
            if i == ell:
                contracted[j, k] -= value
            if i == k:
                contracted[j, ell] += value
    delta = sp.trace(a) - 2 * sp.trace(b)
    u = sp.simplify(a - contracted - delta * sp.eye(4) / 4)
    hodge = sp.zeros(6)
    for i, j, sign in ((0, 5, 1), (1, 4, -1), (2, 3, 1)):
        hodge[i, j] = hodge[j, i] = sign
    hodge_value = sp.simplify(sp.trace(b * hodge * b * hodge))
    squares = sp.simplify(
        2 * sp.trace(u.T * u) + 8 * sp.trace(cblock.T * cblock)
    )
    return f, delta, u, cblock, b, hodge_value, sp.simplify(squares + 4 * hodge_value)


def invariant_moment_form():
    """Give the axial-average square term in the moments E[t],E[t^2],E[t^3]."""
    t = sp.symbols("t")
    m1, m2, m3 = sp.symbols("m1 m2 m3", real=True)
    ring = axial_average_projector(t)

    def momentize(value):
        polynomial = sp.Poly(sp.expand(value), t)
        return sp.expand(
            polynomial.coeff_monomial(1)
            + polynomial.coeff_monomial(t) * m1
            + polynomial.coeff_monomial(t**2) * m2
            + polynomial.coeff_monomial(t**3) * m3
        )

    averaged = ring.applyfunc(momentize)
    f, delta, u, cblock, b, hodge, r0 = split_value(averaged)
    squares = sp.factor(r0 - 4 * hodge)
    return (m1, m2, m3), averaged, f, sp.factor(delta), sp.factor(squares)


def inertia_numeric(matrix: sp.Matrix, tol: float = 1e-9):
    import numpy as np

    values = np.linalg.eigvalsh(np.array(matrix.evalf(), dtype=float))
    return (
        int((values < -tol).sum()),
        int((abs(values) <= tol).sum()),
        int((values > tol).sum()),
        values,
    )


def axial_generator() -> tuple[sp.Matrix, sp.Matrix]:
    """Infinitesimal SO(2) generators on H2 and Lambda^2 H2."""
    j3 = sp.Matrix([[0, 0, 1], [0, 0, 0], [-1, 0, 0]])
    d = sp.zeros(5)
    for column, orbital in enumerate(ORB):
        derivative = j3 * orbital - orbital * j3
        for row, test in enumerate(ORB):
            d[row, column] = sp.trace(test * derivative)
    assert d.T == -d
    wedge = sp.zeros(10)
    edge_index = {edge: i for i, edge in enumerate(EDGES)}

    def add_wedge(row_i, row_j, value, column):
        if row_i == row_j or value == 0:
            return
        if row_i < row_j:
            wedge[edge_index[(row_i, row_j)], column] += value
        else:
            wedge[edge_index[(row_j, row_i)], column] -= value

    for column, (i, j) in enumerate(EDGES):
        for row in range(5):
            add_wedge(row, j, d[row, i], column)
            add_wedge(i, row, d[row, j], column)
    assert wedge.T == -wedge
    return d, wedge


def frequency_blocks(gm: sp.Matrix, form: sp.Matrix):
    """Return exact bases and restrictions for SO(2) frequencies 0,...,6."""
    _, wedge = axial_generator()
    # vec(WG-GW), with row-major vectorization.
    action100 = (
        sp.kronecker_product(wedge, sp.eye(10))
        + sp.kronecker_product(sp.eye(10), wedge)
    )
    image = sp.simplify(action100 * gm)
    # Recover the induced action on the injective 24-coordinate map.
    _, pivot_rows = gm.T.rref()
    rows = list(pivot_rows)
    square = gm[rows, :]
    assert square.det() != 0
    action = sp.simplify(square.inv() * image[rows, :])
    assert gm * action == image
    out = {}
    action2 = sp.simplify(action * action)
    for frequency in range(7):
        basis_vectors = (action2 + frequency**2 * sp.eye(24)).nullspace()
        if not basis_vectors:
            continue
        basis = sp.Matrix.hstack(*basis_vectors)
        block = sp.simplify(basis.T * form * basis)
        out[frequency] = (basis, block)
    assert sum(basis.cols for basis, _ in out.values()) == 24
    return action, out


def verify():
    mass_poly = sp.factor(sum(value * value for value in PLUECKER))
    gm, fm, mass, r0, endpoint = endpoint_form()
    print("Pluecker mass:", mass_poly)
    print("slice dimension:", gm.cols)
    for name, matrix in (("R0", r0), ("R0-mass^2/96", endpoint)):
        negative, zero, positive, values = inertia_numeric(matrix)
        print(name, "inertia", (negative, zero, positive))
        print(name, "smallest eigenvalues", values[:8])
        print(name, "exact rank", matrix.rank())
    _, blocks = frequency_blocks(gm, endpoint)
    for frequency, (basis, block) in blocks.items():
        negative, zero, positive, values = inertia_numeric(block)
        print(
            "endpoint frequency", frequency, "dimension/inertia",
            basis.cols, (negative, zero, positive),
            "smallest", values[: min(4, len(values))],
        )
    moments, _, fbar, deltabar, squares = invariant_moment_form()
    print("invariant F diagonal", [sp.factor(fbar[i, i]) for i in range(5)])
    print("invariant delta", deltabar)
    print("invariant 2||U||^2+8||C||^2", sp.factor(squares))
    print("invariant direct gap", sp.factor(squares - deltabar**2 / 6))
    return gm, fm, mass, r0, endpoint


if __name__ == "__main__":
    verify()
