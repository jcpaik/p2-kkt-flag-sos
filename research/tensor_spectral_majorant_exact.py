"""Exact all-spectrum majorant closing the tangent-star inequality.

For the spectral chart

    S_h = (sqrt(3) E0-h E1)/sqrt(3+h^2),   0 <= h <= 1,

this verifier constructs the pointwise tangent Pluecker blocks U,C and
checks the identity

    delta_x + g_h(x) = <L_h,U_x> + <M_h,C_x>,

where g_h is manifestly nonnegative.  It also checks the exact weighted
dual norm D(h) and the positive Bernstein certificate for D(h) <= 6.

Together with the capped Motzkin--Straus lemma and the block/Hodge identity
in ``fermionic_motzkin_full_derivation.md``, this proves E(mu)>=0 for every
probability measure on RP^2.  Unlike the earlier D2 discovery calculation,
the identity verified here is pointwise and therefore does not require a
twirl.
"""

from __future__ import annotations

import itertools

import sympy as sp


h = sp.symbols("h", real=True)
x, y, z = sp.symbols("x y z", real=True)
point = sp.Matrix([x, y, z])
rt = sp.sqrt

E0 = sp.diag(1, -1, 0) / rt(2)
E1 = sp.diag(1, 1, -2) / rt(6)
S = (rt(3) * E0 - h * E1) / rt(3 + h**2)
H = (h * E0 + rt(3) * E1) / rt(3 + h**2)
ORBITALS = [S, H]
for i, j in ((0, 1), (0, 2), (1, 2)):
    matrix = sp.zeros(3)
    matrix[i, j] = matrix[j, i] = 1 / rt(2)
    ORBITALS.append(matrix)

PAIRS = list(itertools.combinations(range(5), 2))
PLUECKER = [
    sp.expand(2 * point.dot((ORBITALS[i] * point).cross(ORBITALS[j] * point)))
    for i, j in PAIRS
]
G = sp.Matrix(PLUECKER) * sp.Matrix(PLUECKER).T
A = G[:4, :4]
C = G[:4, 4:]
B = G[4:, 4:]


def contraction(matrix: sp.Matrix, pairs, dimension: int) -> sp.Matrix:
    out = sp.zeros(dimension)
    for a, (i, j) in enumerate(pairs):
        for b, (k, ell) in enumerate(pairs):
            value = matrix[a, b]
            if j == ell:
                out[i, k] += value
            if j == k:
                out[i, ell] -= value
            if i == ell:
                out[j, k] -= value
            if i == k:
                out[j, ell] += value
    return out


OUTER_PAIRS = list(itertools.combinations(range(4), 2))
R = contraction(B, OUTER_PAIRS, 4)
MASS = sp.expand(sum(value**2 for value in PLUECKER))
DELTA = sp.expand(3 * sp.trace(A) - 2 * MASS)
U = sp.simplify(A - R - DELTA * sp.eye(4) / 4)

# Manifestly nonnegative pointwise correction.
ELL = (3 - h) * x**2 - (3 + h) * y**2 + 2 * h * z**2
WX = (-17 * h**2 + 32 * h + 20) / 70
WY = (40 * h**2 - 73 * h + 40) / 140
WZ = (12 - 5 * h**2) / 14
CORRECTION = sp.expand((WX * x**2 + WY * y**2 + WZ * z**2) * ELL**2)

# Coefficient matrices in the exact pointwise identity.  L is diagonal.
LCOEFF = [
    -(194 * h**4 - 435 * h**3 + 208 * h**2 + 81 * h + 120) / 280,
    -(206 * h**4 + 435 * h**3 - 1168 * h**2 - 81 * h + 440) / 280,
    (274 * h**4 - 101 * h**3 + 132 * h**2 - 753 * h + 280) / 280,
    (126 * h**4 + 101 * h**3 - 1092 * h**2 + 753 * h + 280) / 280,
]

# The six nonzero entries of M are listed as (row of C, outer-edge column).
C_LOCATIONS = [(1, 0), (1, 5), (2, 1), (2, 4), (3, 2), (3, 3)]
MCOEFF = [
    -rt(3) * h
    * (12 * h**5 + 512 * h**4 - 2779 * h**3 - 3829 * h**2 + 3847 * h + 1741)
    / (840 * (h + 1)),
    -rt(3) * rt(h**2 + 3)
    * (12 * h**5 + 364 * h**4 - 159 * h**3 - 90 * h**2 + 717 * h + 960)
    / (1260 * (h + 1)),
    rt(3)
    * (12 * h**6 - 32 * h**5 - 983 * h**4 - 2373 * h**3 - 2405 * h**2 + 2645 * h + 960)
    / (840 * (h + 1)),
    -rt(3) * rt(h**2 + 3)
    * (6 * h**5 + 887 * h**4 + 582 * h**3 - 1155 * h**2 + 288 * h + 480)
    / (1260 * (h + 1)),
    rt(3)
    * (12 * h**6 + 908 * h**5 - 545 * h**4 - 223 * h**3 + 4749 * h**2 - 85 * h - 960)
    / (840 * (h + 1)),
    rt(3) * rt(h**2 + 3)
    * (6 * h**5 - 523 * h**4 - 741 * h**3 + 1065 * h**2 + 429 * h + 480)
    / (1260 * (h + 1)),
]

EXPECTED_D = (
    1584 * h**12
    + 117408 * h**11
    + 10693280 * h**10
    + 2137032 * h**9
    + 13971489 * h**8
    + 45236046 * h**7
    - 16424313 * h**6
    + 4522956 * h**5
    + 68504799 * h**4
    - 88091010 * h**3
    - 11034711 * h**2
    + 48384000 * h
    + 24192000
) / (5644800 * (h + 1) ** 2)

SIX_MINUS_NUMERATOR = -(
    1584 * h**12
    + 117408 * h**11
    + 10693280 * h**10
    + 2137032 * h**9
    + 13971489 * h**8
    + 45236046 * h**7
    - 16424313 * h**6
    + 4522956 * h**5
    + 68504799 * h**4
    - 88091010 * h**3
    - 44903511 * h**2
    - 19353600 * h
    - 9676800
)


def bernstein_coefficients(polynomial: sp.Expr, degree: int):
    poly = sp.Poly(sp.expand(polynomial), h)
    power = [poly.coeff_monomial(h**k) for k in range(degree + 1)]
    return [
        sp.factor(
            sum(
                power[k] * sp.binomial(j, k) / sp.binomial(degree, k)
                for k in range(j + 1)
            )
        )
        for j in range(degree + 1)
    ]


def verify():
    assert sp.factor(MASS - (x**2 + y**2 + z**2) ** 3) == 0
    right = sum(LCOEFF[i] * U[i, i] for i in range(4))
    right += sum(
        MCOEFF[k] * C[i, j]
        for k, (i, j) in enumerate(C_LOCATIONS)
    )
    assert sp.factor(sp.together(right - DELTA - CORRECTION)) == 0

    dual_cost = sp.factor(
        sum(value**2 for value in LCOEFF) / 2
        + sum(value**2 for value in MCOEFF) / 8
    )
    assert sp.factor(dual_cost - EXPECTED_D) == 0
    assert sp.factor(
        6 - EXPECTED_D
        - SIX_MINUS_NUMERATOR / (5644800 * (h + 1) ** 2)
    ) == 0

    expected_bernstein = [
        9676800,
        11289600,
        sp.Rational(298820637, 22),
        sp.Rational(186523506, 11),
        sp.Rational(3576110762, 165),
        sp.Rational(1838154161, 66),
        sp.Rational(10926245933, 308),
        sp.Rational(88743939, 2),
        sp.Rational(2973164684, 55),
        sp.Rational(3483459192, 55),
        sp.Rational(2287502816, 33),
        sp.Rational(194932736, 3),
        33264640,
    ]
    bernstein = bernstein_coefficients(SIX_MINUS_NUMERATOR, 12)
    assert bernstein == expected_bernstein
    assert all(value > 0 for value in bernstein)

    print("pointwise identity: verified")
    print("weights:", WX, WY, WZ)
    print("dual cost D(h):", sp.factor(EXPECTED_D))
    print("positive Bernstein coefficients of 6-D numerator:", bernstein)


if __name__ == "__main__":
    verify()
