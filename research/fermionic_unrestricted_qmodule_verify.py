"""Exact verifier for the unrestricted fermionic q-module certificate.

This is the finite symbolic check behind
``fermionic_unrestricted_exact_proof.md``.  It verifies, for the complete
spectral parameter interval 0 <= h <= 1, the polynomial identity

    sum_i l_i U_ii + sum_k m_k C_k
      = delta + (w_X x^2+w_Y y^2+w_Z z^2) ell_h^2,

and the exact weighted dual cost D(h)<=6.  No floating-point calculation
or SDP output is used by this verifier.
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
A, C, B = G[:4, :4], G[:4, 4:], G[4:, 4:]


def outer_contraction(B: sp.Matrix) -> sp.Matrix:
    outer_pairs = [(i - 1, j - 1) for i, j in PAIRS[4:]]
    out = sp.zeros(4)
    for a, (i, j) in enumerate(outer_pairs):
        for b, (k, ell) in enumerate(outer_pairs):
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


MASS = sp.expand(sum(PLUECKER[i] ** 2 for i in range(10)))
DELTA = sp.expand(3 * sp.trace(A) - 2 * MASS)
U = sp.simplify(A - outer_contraction(B) - DELTA * sp.eye(4) / 4)

# The six D2-invariant C entries, in the outer-edge order
# (H^xy,H^xz,H^yz,xy^xz,xy^yz,xz^yz).
C_LOCATIONS = ((1, 0), (1, 5), (2, 1), (2, 4), (3, 2), (3, 3))

WEIGHTS = (
    (-17 * h**2 + 32 * h + 20) / 70,
    (40 * h**2 - 73 * h + 40) / 140,
    (12 - 5 * h**2) / 14,
)
ELL = (3 - h) * x**2 - (3 + h) * y**2 + 2 * h * z**2
G_NONNEGATIVE = (WEIGHTS[0] * x**2 + WEIGHTS[1] * y**2 + WEIGHTS[2] * z**2) * ELL**2

L_COEFFICIENTS = (
    -(194 * h**4 - 435 * h**3 + 208 * h**2 + 81 * h + 120) / 280,
    -(206 * h**4 + 435 * h**3 - 1168 * h**2 - 81 * h + 440) / 280,
    (274 * h**4 - 101 * h**3 + 132 * h**2 - 753 * h + 280) / 280,
    (126 * h**4 + 101 * h**3 - 1092 * h**2 + 753 * h + 280) / 280,
)

M_COEFFICIENTS = (
    -rt(3) * h * (12 * h**5 + 512 * h**4 - 2779 * h**3 - 3829 * h**2 + 3847 * h + 1741) / (840 * (h + 1)),
    -rt(3) * rt(h**2 + 3) * (12 * h**5 + 364 * h**4 - 159 * h**3 - 90 * h**2 + 717 * h + 960) / (1260 * (h + 1)),
    rt(3) * (12 * h**6 - 32 * h**5 - 983 * h**4 - 2373 * h**3 - 2405 * h**2 + 2645 * h + 960) / (840 * (h + 1)),
    -rt(3) * rt(h**2 + 3) * (6 * h**5 + 887 * h**4 + 582 * h**3 - 1155 * h**2 + 288 * h + 480) / (1260 * (h + 1)),
    rt(3) * (12 * h**6 + 908 * h**5 - 545 * h**4 - 223 * h**3 + 4749 * h**2 - 85 * h - 960) / (840 * (h + 1)),
    rt(3) * rt(h**2 + 3) * (6 * h**5 - 523 * h**4 - 741 * h**3 + 1065 * h**2 + 429 * h + 480) / (1260 * (h + 1)),
)

P_COST = (
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
)
D = P_COST / (5644800 * (h + 1) ** 2)

GAP_NUMERATOR = (
    -1584 * h**12
    - 117408 * h**11
    - 10693280 * h**10
    - 2137032 * h**9
    - 13971489 * h**8
    - 45236046 * h**7
    + 16424313 * h**6
    - 4522956 * h**5
    - 68504799 * h**4
    + 88091010 * h**3
    + 44903511 * h**2
    + 19353600 * h
    + 9676800
)

EXPECTED_GAP_BERNSTEIN = (
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
)


def bernstein_coefficients(polynomial: sp.Expr, degree: int):
    poly = sp.Poly(sp.expand(polynomial), h)
    return tuple(
        sp.factor(
            sum(
                poly.coeff_monomial(h**k)
                * sp.binomial(j, k)
                / sp.binomial(degree, k)
                for k in range(j + 1)
            )
        )
        for j in range(degree + 1)
    )


def verify(verbose: bool = True):
    representation = sum(L_COEFFICIENTS[i] * U[i, i] for i in range(4))
    representation += sum(
        M_COEFFICIENTS[k] * C[i, j]
        for k, (i, j) in enumerate(C_LOCATIONS)
    )
    identity_residual = sp.factor(sp.together(representation - DELTA - G_NONNEGATIVE))
    assert identity_residual == 0

    computed_cost = sp.factor(
        sum(value**2 for value in L_COEFFICIENTS) / 2
        + sum(value**2 for value in M_COEFFICIENTS) / 8
    )
    assert sp.factor(computed_cost - D) == 0
    assert sp.factor(6 - D - GAP_NUMERATOR / (5644800 * (h + 1) ** 2)) == 0

    gap_bernstein = bernstein_coefficients(GAP_NUMERATOR, 12)
    assert gap_bernstein == EXPECTED_GAP_BERNSTEIN
    assert all(value > 0 for value in gap_bernstein)

    # Positivity of the three q-module weights on 0<=h<=1:
    # w_X is concave and positive at both endpoints; w_Y has its exact
    # global minimum at h=73/80; w_Z decreases to 1/2.
    assert WEIGHTS[0].subs(h, 0) > 0 and WEIGHTS[0].subs(h, 1) > 0
    assert sp.factor(WEIGHTS[1].subs(h, sp.Rational(73, 80))) == sp.Rational(1071, 22400)
    assert WEIGHTS[2].subs(h, 1) == sp.Rational(1, 2)

    if verbose:
        print("exact polynomial identity: OK")
        print("weights", WEIGHTS)
        print("D(h) =", D)
        print("6-D Bernstein coefficients", gap_bernstein)
    return True


if __name__ == "__main__":
    verify()
