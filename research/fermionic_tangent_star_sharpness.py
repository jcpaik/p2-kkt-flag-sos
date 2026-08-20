"""Exact audit of the sharp tangent--Motzkin constant.

The five-point family below is a split orthonormal-basis zero.  It proves
that the constant 7/120 in

    T = R_off + 2 D_b + ||u_0||^2/2 >= (7/120) Delta^2

cannot be increased.  Equivalently, since E/8 = T-Delta^2/24, the
coefficient 1/60 in E/8 >= Delta^2/60 is sharp.
"""

from __future__ import annotations

import sympy as sp


u = sp.symbols("u", nonnegative=True)
sqrt = sp.sqrt

H2 = [
    sp.diag(1, -1, 0) / sqrt(2),
    sp.diag(1, 1, -2) / sqrt(6),
]
for i, j in ((0, 1), (0, 2), (1, 2)):
    matrix = sp.zeros(3)
    matrix[i, j] = matrix[j, i] = 1 / sqrt(2)
    H2.append(matrix)

points = [
    sp.Matrix([sqrt(1 - u), sqrt(u), 0]),
    sp.Matrix([sqrt(1 - u), -sqrt(u), 0]),
    sp.Matrix([0, 1, 0]),
    sp.Matrix([0, sqrt(u), sqrt(1 - u)]),
    sp.Matrix([0, -sqrt(u), sqrt(1 - u)]),
]
weights = [
    sp.Rational(1, 6) + u,
    sp.Rational(1, 6) + u,
    sp.Rational(1, 3) - 4 * u,
    sp.Rational(1, 6) + u,
    sp.Rational(1, 6) + u,
]


def tangent_projector(x: sp.Matrix) -> sp.Matrix:
    values = [(x.T * matrix * x)[0] for matrix in H2]
    images = [matrix * x for matrix in H2]
    return sp.Matrix(
        5,
        5,
        lambda a, b: 2
        * ((images[a].T * images[b])[0] - values[a] * values[b]),
    )


F = sp.simplify(
    sum(
        (weight * tangent_projector(x) for weight, x in zip(weights, points)),
        sp.zeros(5),
    )
)

# The fourth basis vector is the top orbital throughout 0 <= u <= 1/12.
c = sp.factor(F[3, 3])
assert sp.expand(c - (sp.Rational(2, 3) + sp.Rational(10, 3) * u - 4 * u**2)) == 0
Delta = sp.factor(3 * c - 2)
assert sp.expand(Delta - 2 * u * (5 - 6 * u)) == 0

# The pair kernel equals E/8.
purity = 0
for wi, xi in zip(weights, points):
    for wj, xj in zip(weights, points):
        s = sp.expand((xi.dot(xj)) ** 2)
        kernel = 4 * s**3 - 6 * s**2 + sp.Rational(5, 2) * s - sp.Rational(1, 6)
        purity += wi * wj * kernel
purity = sp.factor(purity)

P6 = (
    3168 * u**6
    - 8160 * u**5
    + 7240 * u**4
    - 2608 * u**3
    + 368 * u**2
    - 12 * u
    + 5
)
assert purity == u**2 * P6 / 3

gap = sp.factor(purity - Delta**2 / 60)
positive_quartic = 3960 * u**4 - 10200 * u**3 + 9050 * u**2 - 3260 * u + 451
assert gap == 4 * u**4 * positive_quartic / 15

# On u=v/12, the quartic has the following strictly positive Bernstein
# coefficients (degree four), proving the displayed gap is nonnegative on
# the whole probability interval 0 <= u <= 1/12.
bernstein_coefficients = [
    sp.Integer(451),
    sp.Rational(4597, 12),
    sp.Rational(140677, 432),
    sp.Rational(26611, 96),
    sp.Rational(7567, 32),
]
v = sp.symbols("v")
bernstein = sum(
    coefficient * sp.binomial(4, k) * v**k * (1 - v) ** (4 - k)
    for k, coefficient in enumerate(bernstein_coefficients)
)
assert sp.expand(bernstein - positive_quartic.subs(u, v / 12)) == 0

# T = E/8 + Delta^2/24, so Delta^2/(24 T) tends to 5/7.
T = sp.factor(purity + Delta**2 / 24)
assert sp.limit(Delta**2 / (24 * T), u, 0, dir="+") == sp.Rational(5, 7)

print("c =", c)
print("Delta =", Delta)
print("E/8 =", purity)
print("E/8 - Delta^2/60 =", gap)
print("lim Delta^2/(24T) =", sp.limit(Delta**2 / (24 * T), u, 0, dir="+"))
