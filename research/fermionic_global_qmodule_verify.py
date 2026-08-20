"""Exact verifier for the unrestricted fermionic q-module certificate.

For every spectral parameter 0 <= h <= 1 this checks the pointwise
identity

    delta_x + g_h(x) = sum_i l_i(h) U_ii(x)
                       + sum_j m_j(h) C_j(x),

where g_h is manifestly nonnegative.  It also checks the exact dual cost
D(h), positivity of the three correction weights, and D(h) <= 6 through
strictly positive Bernstein coefficients on [0,1].
"""

from __future__ import annotations

import itertools

import sympy as sp


h = sp.symbols("h", real=True)
x, y, z = sp.symbols("x y z", real=True)
rt = sp.sqrt


def orbital_data():
    e0 = sp.diag(1, -1, 0) / rt(2)
    e1 = sp.diag(1, 1, -2) / rt(6)
    sh = (rt(3) * e0 - h * e1) / rt(3 + h**2)
    hh = (h * e0 + rt(3) * e1) / rt(3 + h**2)
    orbitals = [sh, hh]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        matrix = sp.zeros(3)
        matrix[i, j] = matrix[j, i] = 1 / rt(2)
        orbitals.append(matrix)
    assert sp.simplify(sp.Matrix([
        [sp.trace(a * b) for b in orbitals] for a in orbitals
    ]) - sp.eye(5)) == sp.zeros(5)
    return orbitals


def tangent_blocks():
    orbitals = orbital_data()
    pairs = list(itertools.combinations(range(5), 2))
    point = sp.Matrix((x, y, z))
    zeta = sp.Matrix([
        sp.expand(2 * point.dot((orbitals[i] * point).cross(orbitals[j] * point)))
        for i, j in pairs
    ])
    gram = zeta * zeta.T
    a = gram[:4, :4]
    cblock = gram[:4, 4:]
    b = gram[4:, 4:]
    outer_pairs = [(i - 1, j - 1) for i, j in pairs[4:]]
    contracted = sp.zeros(4)
    for edge_a, (i, j) in enumerate(outer_pairs):
        for edge_b, (k, ell) in enumerate(outer_pairs):
            value = b[edge_a, edge_b]
            if j == ell:
                contracted[i, k] += value
            if j == k:
                contracted[i, ell] -= value
            if i == ell:
                contracted[j, k] -= value
            if i == k:
                contracted[j, ell] += value
    mass = sp.expand(sum(value**2 for value in zeta))
    delta = sp.expand(3 * sp.trace(a) - 2 * mass)
    u = sp.simplify(a - contracted - delta * sp.eye(4) / 4)
    return zeta, mass, delta, u, cblock


def certificate_coefficients():
    weights = (
        (-17 * h**2 + 32 * h + 20) / 70,
        (40 * h**2 - 73 * h + 40) / 140,
        (12 - 5 * h**2) / 14,
    )
    ell = (3 - h) * x**2 - (3 + h) * y**2 + 2 * h * z**2
    correction = (weights[0] * x**2 + weights[1] * y**2 + weights[2] * z**2) * ell**2

    l = (
        -(194*h**4 - 435*h**3 + 208*h**2 + 81*h + 120) / 280,
        -(206*h**4 + 435*h**3 - 1168*h**2 - 81*h + 440) / 280,
        (274*h**4 - 101*h**3 + 132*h**2 - 753*h + 280) / 280,
        (126*h**4 + 101*h**3 - 1092*h**2 + 753*h + 280) / 280,
    )
    m = (
        -rt(3)*h*(12*h**5 + 512*h**4 - 2779*h**3 - 3829*h**2 + 3847*h + 1741)
        / (840*(h + 1)),
        -rt(3)*rt(h**2 + 3)*(12*h**5 + 364*h**4 - 159*h**3 - 90*h**2 + 717*h + 960)
        / (1260*(h + 1)),
        rt(3)*(12*h**6 - 32*h**5 - 983*h**4 - 2373*h**3 - 2405*h**2 + 2645*h + 960)
        / (840*(h + 1)),
        -rt(3)*rt(h**2 + 3)*(6*h**5 + 887*h**4 + 582*h**3 - 1155*h**2 + 288*h + 480)
        / (1260*(h + 1)),
        rt(3)*(12*h**6 + 908*h**5 - 545*h**4 - 223*h**3 + 4749*h**2 - 85*h - 960)
        / (840*(h + 1)),
        rt(3)*rt(h**2 + 3)*(6*h**5 - 523*h**4 - 741*h**3 + 1065*h**2 + 429*h + 480)
        / (1260*(h + 1)),
    )
    return weights, ell, correction, l, m


def bernstein_coefficients(polynomial: sp.Expr, variable: sp.Symbol, degree: int):
    poly = sp.Poly(sp.expand(polynomial), variable)
    power = [poly.coeff_monomial(variable**j) for j in range(degree + 1)]
    return [
        sp.factor(sum(
            power[j] * sp.binomial(k, j) / sp.binomial(degree, j)
            for j in range(k + 1)
        ))
        for k in range(degree + 1)
    ]


def verify(verbose: bool = True):
    zeta, mass, delta, u, cblock = tangent_blocks()
    weights, ell, correction, l, m = certificate_coefficients()
    assert sp.factor(mass - (x**2 + y**2 + z**2)**3) == 0

    # These are the six D2-invariant cross-block entries used by the
    # certificate, indexed inside the 4x6 C block.
    locations = ((1, 0), (1, 5), (2, 1), (2, 4), (3, 2), (3, 3))
    linear = sum(l[i] * u[i, i] for i in range(4))
    linear += sum(m[k] * cblock[i, j] for k, (i, j) in enumerate(locations))
    identity_numerator = sp.factor(sp.together(linear - delta - correction))
    assert identity_numerator == 0

    cost = sp.factor(sum(value**2 for value in l) / 2 + sum(value**2 for value in m) / 8)
    expected_cost = (
        1584*h**12 + 117408*h**11 + 10693280*h**10 + 2137032*h**9
        + 13971489*h**8 + 45236046*h**7 - 16424313*h**6
        + 4522956*h**5 + 68504799*h**4 - 88091010*h**3
        - 11034711*h**2 + 48384000*h + 24192000
    ) / (5644800 * (h + 1)**2)
    assert sp.factor(cost - expected_cost) == 0

    # Positive degree-two Bernstein coefficients prove w_i >= 0 on [0,1].
    weight_bernstein = [bernstein_coefficients(value, h, 2) for value in weights]
    assert weight_bernstein == [
        [sp.Rational(2, 7), sp.Rational(18, 35), sp.Rational(1, 2)],
        [sp.Rational(2, 7), sp.Rational(1, 40), sp.Rational(1, 20)],
        [sp.Rational(6, 7), sp.Rational(6, 7), sp.Rational(1, 2)],
    ]

    # Clear the positive denominator in 6-D.  Every degree-12 Bernstein
    # coefficient of the numerator is strictly positive.
    numerator = sp.factor(5644800 * (h + 1)**2 * (6 - cost))
    cost_bernstein = bernstein_coefficients(numerator, h, 12)
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
    assert cost_bernstein == expected_bernstein
    assert all(value > 0 for value in cost_bernstein)

    if verbose:
        print("Pluecker mass =", mass)
        print("weights =", weights)
        print("weight Bernstein coefficients =", weight_bernstein)
        print("dual cost D(h) =", sp.factor(cost))
        print("Bernstein coefficients of cleared 6-D =", cost_bernstein)
        print("pointwise identity verified")
    return {
        "weights": weights,
        "ell": ell,
        "cost": cost,
        "cost_bernstein": cost_bernstein,
    }


if __name__ == "__main__":
    verify()
