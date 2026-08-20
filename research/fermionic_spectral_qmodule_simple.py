"""Exact verifier for the simple all-spectrum q-module certificate.

This is the final non-isotropic bridge.  For 0 <= h <= 1 it verifies

    delta + E[(w_X X+w_Y Y+w_Z Z) ell_h^2]
        = <L_h,U> + <M_h,C>,

with nonnegative weights and weighted dual cost at most six.  No floating
point or SDP calculation is used.
"""

from __future__ import annotations

import sympy as sp

from fermionic_spectral_qmodule_exact import X, Y, Z, h, maps, polynomial_row


def bernstein_coefficients(poly: sp.Expr, variable: sp.Symbol, degree: int):
    """Power-to-Bernstein conversion on [0,1]."""
    expanded = sp.Poly(sp.expand(poly), variable)
    power = [expanded.nth(j) for j in range(degree + 1)]
    return [
        sp.factor(
            sum(
                power[j] * sp.binomial(k, j) / sp.binomial(degree, j)
                for j in range(k + 1)
            )
        )
        for k in range(degree + 1)
    ]


def certificate():
    _, _, _, delta, u_rows, c_rows = maps()
    assert len(u_rows) == 4
    assert len(c_rows) == 6

    ell = (3 - h) * X - (3 + h) * Y + 2 * h * Z
    wx = (-17 * h**2 + 32 * h + 20) / 70
    wy = (40 * h**2 - 73 * h + 40) / 140
    wz = (12 - 5 * h**2) / 14
    correction = polynomial_row((wx * X + wy * Y + wz * Z) * ell**2)

    l = [
        -(194 * h**4 - 435 * h**3 + 208 * h**2 + 81 * h + 120) / 280,
        -(206 * h**4 + 435 * h**3 - 1168 * h**2 - 81 * h + 440) / 280,
        (274 * h**4 - 101 * h**3 + 132 * h**2 - 753 * h + 280) / 280,
        (126 * h**4 + 101 * h**3 - 1092 * h**2 + 753 * h + 280) / 280,
    ]
    root = sp.sqrt(h**2 + 3)
    m = [
        -sp.sqrt(3) * h
        * (12 * h**5 + 512 * h**4 - 2779 * h**3 - 3829 * h**2 + 3847 * h + 1741)
        / (840 * (h + 1)),
        -sp.sqrt(3) * root
        * (12 * h**5 + 364 * h**4 - 159 * h**3 - 90 * h**2 + 717 * h + 960)
        / (1260 * (h + 1)),
        sp.sqrt(3)
        * (12 * h**6 - 32 * h**5 - 983 * h**4 - 2373 * h**3 - 2405 * h**2 + 2645 * h + 960)
        / (840 * (h + 1)),
        -sp.sqrt(3) * root
        * (6 * h**5 + 887 * h**4 + 582 * h**3 - 1155 * h**2 + 288 * h + 480)
        / (1260 * (h + 1)),
        sp.sqrt(3)
        * (12 * h**6 + 908 * h**5 - 545 * h**4 - 223 * h**3 + 4749 * h**2 - 85 * h - 960)
        / (840 * (h + 1)),
        sp.sqrt(3) * root
        * (6 * h**5 - 523 * h**4 - 741 * h**3 + 1065 * h**2 + 429 * h + 480)
        / (1260 * (h + 1)),
    ]

    represented = sum(
        (coefficient * row for coefficient, row in zip(l, u_rows)),
        sp.zeros(1, 10),
    )
    represented += sum(
        (coefficient * row for coefficient, row in zip(m, c_rows)),
        sp.zeros(1, 10),
    )
    assert all(
        sp.factor(value) == 0
        for value in (represented - delta - correction)
    )

    dual_cost = sp.factor(sum(value**2 for value in l) / 2 + sum(value**2 for value in m) / 8)
    numerator = (
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
    claimed_cost = numerator / (5644800 * (h + 1) ** 2)
    assert sp.factor(dual_cost - claimed_cost) == 0

    gap_numerator = sp.factor(5644800 * (h + 1) ** 2 * (6 - claimed_cost))
    expected_gap = (
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
    assert sp.expand(gap_numerator - expected_gap) == 0

    # Positivity on [0,1] is exact because all Bernstein coefficients are
    # positive.  The correction weights use degree two; the gap uses degree
    # twelve.
    weight_polynomials = [70 * wx, 140 * wy, 14 * wz]
    for polynomial in weight_polynomials:
        assert all(value > 0 for value in bernstein_coefficients(polynomial, h, 2))
    gap_bernstein = bernstein_coefficients(expected_gap, h, 12)
    assert all(value > 0 for value in gap_bernstein)

    return {
        "weights": (wx, wy, wz),
        "ell": ell,
        "L": l,
        "M": m,
        "cost": claimed_cost,
        "gap_bernstein": gap_bernstein,
    }


if __name__ == "__main__":
    result = certificate()
    print("simple all-spectrum q-module certificate verified exactly")
    print("weights:", result["weights"])
    print("dual cost:", result["cost"])
    print("Bernstein coefficients of 6-cost:", result["gap_bernstein"])
