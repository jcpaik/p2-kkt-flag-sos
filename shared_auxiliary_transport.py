"""Collective-transport Hessian with one shared random auxiliary line.

For every fixed projective line ``z`` and every odd positive integer ``r``,

    V_{r,z}(x) = (x.z)^r P_{x^perp} z

is a well-defined tangent vector field on RP^2.  At a minimizing measure the
collective second variation is nonnegative for every linear combination of
these fields.  Averaging that PSD quadratic form over an independent
``z ~ mu`` produces a three-sample (rather than four-sample) moment block.

This block is tailored to the root-determinant weight F: both live entirely
in the triangle-label sector, and the shared auxiliary retains the
orientation/Gram-determinant coupling lost when the two vector-field leaves
are averaged independently.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import sympy as sp

from sos_search import Label, Polynomial, expectation_label


def _convert(poly: sp.Poly) -> Polynomial:
    return [
        (
            Fraction(int(coefficient)),
            tuple(int(value) for value in exponent),
        )
        for exponent, coefficient in poly.terms()
    ]


def shared_auxiliary_transport_polynomials() -> tuple[Polynomial, Polynomial]:
    """Base local/cross kernels in variables ``(a,b,c)=(X.Y,Y.Z,Z.X)``."""

    a, b, c = sp.symbols("a b c")
    kernel_prime = 192 * a**5 - 192 * a**3 + 40 * a
    kernel_second = 960 * a**4 - 576 * a**2 + 40

    pxz_y = b - a * c
    pyz_x = c - a * b
    pxz_norm_squared = 1 - c**2
    pxz_pyz = 1 - b**2 - c**2 + a * b * c

    local = sp.Poly(
        sp.expand(
            kernel_second * pxz_y**2
            - a * kernel_prime * pxz_norm_squared
        ),
        a,
        b,
        c,
    )
    cross = sp.Poly(
        sp.expand(
            kernel_second * pxz_y * pyz_x
            + kernel_prime * pxz_pyz
        ),
        a,
        b,
        c,
    )
    return _convert(local), _convert(cross)


def shared_auxiliary_transport_expectation_matrix(
    auxiliary_degrees: list[int],
) -> dict[Label, np.ndarray]:
    """Return the exact-label PSD block indexed by odd field powers."""

    if any(degree <= 0 or degree % 2 == 0 for degree in auxiliary_degrees):
        raise ValueError("shared-auxiliary degrees must be positive odd integers")
    local, cross = shared_auxiliary_transport_polynomials()
    size = len(auxiliary_degrees)
    matrices: dict[Label, np.ndarray] = {}

    def accumulate(
        row: int,
        column: int,
        polynomial: Polynomial,
        shift: tuple[int, int, int],
    ) -> None:
        for coefficient, exponent in polynomial:
            shifted = tuple(
                exponent[index] + shift[index] for index in range(3)
            )
            label, reduction = expectation_label(shifted)  # type: ignore[arg-type]
            if label is None or not reduction:
                continue
            matrix = matrices.setdefault(label, np.zeros((size, size)))
            matrix[row, column] += float(coefficient * reduction)

    for row, left_degree in enumerate(auxiliary_degrees):
        for column, right_degree in enumerate(auxiliary_degrees):
            # Both local fields live at X, hence powers c^(r+s).
            accumulate(
                row,
                column,
                local,
                (0, 0, left_degree + right_degree),
            )
            # The cross fields live at X and Y, hence c^r b^s.
            accumulate(
                row,
                column,
                cross,
                (0, right_degree, left_degree),
            )

    return {
        label: (matrix + matrix.T) / 2.0
        for label, matrix in matrices.items()
        if np.max(np.abs((matrix + matrix.T) / 2.0)) > 1e-13
    }


def shared_auxiliary_transport_degrees(total_degree: int) -> list[int]:
    """Odd powers whose pairwise entries fit under the Gram-degree cap."""

    maximum = (total_degree - 8) // 2
    return list(range(1, maximum + 1, 2))

