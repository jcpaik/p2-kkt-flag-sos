"""Collective transport second-variation constraints for the P2 energy.

This module adds the part of the positional Hessian that is missed by the
pointwise support Hessian condition.  At a minimizing measure, push every
point x along a tangent field V(x).  If X_t is the spherical geodesic with
initial velocity V(X), then

    d^2/dt^2 E[K(X_t . Y_t)]|_{t=0} >= 0.

For the vector-field basis

    V_r(x) = E_Z[(x.Z)^r P_{x^perp} Z],

with odd r (so the lifted field is well defined on RP^2), this yields a PSD
matrix of four-point Gram moments.  The returned matrix represents one half of
the full second variation, which has the same PSD condition.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import sympy as sp

from sos_search import GraphPolynomial, Label, graph_expectation_label


def _convert(poly: sp.Poly) -> GraphPolynomial:
    return [
        (
            Fraction(int(coefficient)),
            tuple(int(value) for value in exponent),
        )
        for exponent, coefficient in poly.terms()
    ]


def collective_transport_polynomials() -> tuple[GraphPolynomial, GraphPolynomial]:
    """Return the local and cross kernels in the collective second variation.

    Vertices are ordered (X,Y,Z,W), with Gram variables

        XY=a, XZ=c, XW=d, YZ=b, YW=e, ZW=f.

    Z generates a tangent vector at X and W generates either a second tangent
    vector at X (the local term) or a tangent vector at Y (the cross term).

    If v=P_{X^perp} Z and u=P_{X^perp} W, the local bilinear term is

        K''(a) (v.Y)(u.Y) - a K'(a) (v.u).

    If v=P_{X^perp} Z and w=P_{Y^perp} W, the off-diagonal transport term is

        K''(a) (v.Y)(X.w) + K'(a) (v.w).

    The full second variation divided by two is local + cross.
    """

    a, c, d, b, e, f = sp.symbols("a c d b e f")
    kernel_prime = 192 * a**5 - 192 * a**3 + 40 * a
    kernel_second = 960 * a**4 - 576 * a**2 + 40

    xz_y = b - a * c
    xw_y = e - a * d
    xz_xw = f - c * d

    # P_X Z paired with P_Y W.
    x_yw = d - a * e
    xz_yw = f - e * b - c * d + a * c * e

    local = sp.Poly(
        sp.expand(
            kernel_second * xz_y * xw_y
            - a * kernel_prime * xz_xw
        ),
        a,
        c,
        d,
        b,
        e,
        f,
    )
    cross = sp.Poly(
        sp.expand(
            kernel_second * xz_y * x_yw
            + kernel_prime * xz_yw
        ),
        a,
        c,
        d,
        b,
        e,
        f,
    )
    return _convert(local), _convert(cross)


def collective_transport_expectation_matrix(
    auxiliary_degrees: list[int],
) -> dict[Label, np.ndarray]:
    """Moment matrix for collective tangent transports.

    ``auxiliary_degrees`` should normally contain positive odd integers.  The
    row indexed by r corresponds to

        V_r(x) = E_Z[(x.Z)^r P_{x^perp} Z].

    For row r and column s:

    * the local term is weighted by (X.Z)^r (X.W)^s;
    * the cross term is weighted by (X.Z)^r (Y.W)^s.

    The matrix is symmetrized explicitly after canonical moment reduction.
    Only its symmetric part contributes to c^T M c, and the exact transport
    bilinear form is symmetric by exchangeability of the iid samples.
    """

    if any(degree <= 0 or degree % 2 == 0 for degree in auxiliary_degrees):
        raise ValueError("collective transport degrees must be positive odd integers")

    local, cross = collective_transport_polynomials()
    size = len(auxiliary_degrees)
    matrices: dict[Label, np.ndarray] = {}

    def add_terms(
        row: int,
        column: int,
        polynomial: GraphPolynomial,
        shift_edges: tuple[tuple[int, int], tuple[int, int]],
    ) -> None:
        for coefficient, exponent in polynomial:
            shifted = list(exponent)
            shifted[shift_edges[0][0]] += shift_edges[0][1]
            shifted[shift_edges[1][0]] += shift_edges[1][1]
            label, reduction_coefficient = graph_expectation_label(4, tuple(shifted))
            if label is None or reduction_coefficient == 0:
                continue
            matrix = matrices.setdefault(label, np.zeros((size, size)))
            matrix[row, column] += float(coefficient * reduction_coefficient)

    # Edge order for four vertices is XY, XZ, XW, YZ, YW, ZW.
    for row, left_degree in enumerate(auxiliary_degrees):
        for column, right_degree in enumerate(auxiliary_degrees):
            add_terms(
                row,
                column,
                local,
                ((1, left_degree), (2, right_degree)),
            )
            add_terms(
                row,
                column,
                cross,
                ((1, left_degree), (4, right_degree)),
            )

    result: dict[Label, np.ndarray] = {}
    for label, matrix in matrices.items():
        symmetric = (matrix + matrix.T) / 2.0
        if np.max(np.abs(symmetric)) > 1e-13:
            result[label] = symmetric
    return result


def maximum_collective_degree(total_degree: int) -> int:
    """Largest vector-field exponent allowed by a total Gram-degree cap."""

    return (total_degree - 8) // 2


def collective_transport_degrees(total_degree: int) -> list[int]:
    """Default odd vector-field powers fitting within ``total_degree``."""

    maximum = maximum_collective_degree(total_degree)
    return list(range(1, maximum + 1, 2))
