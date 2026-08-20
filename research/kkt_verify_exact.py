#!/usr/bin/env python3
"""Exact checks for research/kkt_classification_exact.md."""

from __future__ import annotations

from itertools import product

import sympy as sp


def kernel_from_squared_dot(s: sp.Expr) -> sp.Expr:
    return 32 * s**3 - 48 * s**2 + 20 * s - sp.Rational(4, 3)


def bernstein_coefficients(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    left: sp.Rational,
    right: sp.Rational,
) -> list[sp.Expr]:
    parameter = sp.symbols("parameter")
    transformed = sp.Poly(
        sp.expand(polynomial.subs(variable, left + (right - left) * parameter)),
        parameter,
    )
    degree = transformed.degree()
    power = [transformed.nth(index) for index in range(degree + 1)]
    return [
        sp.factor(
            sum(
                power[index]
                * sp.binomial(column, index)
                / sp.binomial(degree, index)
                for index in range(column + 1)
            )
        )
        for column in range(degree + 1)
    ]


def check_robinson() -> None:
    points: list[tuple[list[int], int]] = []
    for zero_coordinate in range(3):
        other = [index for index in range(3) if index != zero_coordinate]
        for sign in (1, -1):
            vector = [0, 0, 0]
            vector[other[0]] = 1
            vector[other[1]] = sign
            points.append((vector, 2))
    for first, second in product((1, -1), repeat=2):
        points.append(([1, first, second], 3))

    matrix = sp.zeros(10)
    for row, (x, norm_x) in enumerate(points):
        for column, (y, norm_y) in enumerate(points):
            dot = sum(x[index] * y[index] for index in range(3))
            matrix[row, column] = kernel_from_squared_dot(
                sp.Rational(dot**2, norm_x * norm_y)
            )

    weights = sp.Matrix(
        [sp.Rational(184, 2265)] * 6 + [sp.Rational(387, 3020)] * 4
    )
    energy = sp.Rational(416, 2265)
    assert sum(weights) == 1
    assert matrix * weights == energy * sp.ones(10, 1)
    assert (weights.T * matrix * weights)[0] == energy

    negative_directions = (
        sp.Matrix([-1, -1, 1, 1, 0, 0, 0, 0, 0, 0]),
        sp.Matrix([-1, -1, 0, 0, 1, 1, 0, 0, 0, 0]),
    )
    for direction in negative_directions:
        assert sum(direction) == 0
        assert matrix * direction == -direction

    x, y, z = sp.symbols("x y z")
    radius_squared = x**2 + y**2 + z**2
    coordinates = (x, y, z)
    potential = 0
    for (point, norm), weight in zip(points, weights, strict=True):
        dot = sum(point[index] * coordinates[index] for index in range(3))
        potential += weight * (
            32 * dot**6 / norm**3
            - 48 * radius_squared * dot**4 / norm**2
            + 20 * radius_squared**2 * dot**2 / norm
            - sp.Rational(4, 3) * radius_squared**3
        )
    robinson = (
        x**6
        + y**6
        + z**6
        - x**4 * y**2
        - x**4 * z**2
        - y**4 * x**2
        - y**4 * z**2
        - z**4 * x**2
        - z**4 * y**2
        + 3 * x**2 * y**2 * z**2
    )
    gap = sp.expand(potential - energy * radius_squared**3)
    assert sp.expand(gap - sp.Rational(64, 151) * robinson) == 0


def check_latitude() -> None:
    a, v = sp.symbols("a v")
    potential = sp.Rational(2, 3) * (
        693 * a**3 * v**3
        - 945 * a**3 * v**2
        + 315 * a**3 * v
        - 15 * a**3
        - 945 * a**2 * v**3
        + 1260 * a**2 * v**2
        - 405 * a**2 * v
        + 18 * a**2
        + 315 * a * v**3
        - 405 * a * v**2
        + 126 * a * v
        - 6 * a
        - 15 * v**3
        + 18 * v**2
        - 6 * v
        + 1
    )
    stationarity = (
        693 * a**5
        - 1575 * a**4
        + 1260 * a**3
        - 420 * a**2
        + 54 * a
        - 2
    )
    assert sp.diff(potential, v).subs(v, a) == 2 * stationarity

    coefficient = 2 * (231 * a**3 - 315 * a**2 + 105 * a - 5)
    constant = 2 * (462 * a**4 - 945 * a**3 + 630 * a**2 - 145 * a + 6)
    claimed_factorization = (v - a) ** 2 * (coefficient * v + constant)
    gap = sp.expand(potential - potential.subs(v, a))
    assert sp.expand(gap - claimed_factorization - 2 * (v - a) * stationarity) == 0

    left = sp.Rational(421, 1000)
    right = sp.Rational(8, 19)
    assert stationarity.subs(a, left) < 0 < stationarity.subs(a, right)
    assert sp.count_roots(stationarity, left, right) == 1

    energy = sp.factor(potential.subs(v, a))
    mode_two_inner = 495 * a**4 - 540 * a**3 + 162 * a**2 - 12 * a + 1
    positive_polynomials = (coefficient, constant - coefficient, energy)
    for polynomial in positive_polynomials:
        assert all(
            value > 0
            for value in bernstein_coefficients(polynomial, a, left, right)
        )
    assert all(
        value < 0
        for value in bernstein_coefficients(mode_two_inner, a, left, right)
    )


def check_spin_four_and_onb_identities() -> None:
    real, imaginary = sp.symbols("real imaginary", real=True)
    complex_value = real + sp.I * imaginary
    assert sp.expand(
        (real**2 + imaginary**2) ** 2 - sp.re(sp.expand(complex_value**4))
    ) == 8 * real**2 * imaginary**2

    first, second, third = sp.symbols("first second third", real=True)
    onb_sum = sum(
        kernel_from_squared_dot(coordinate**2)
        for coordinate in (first, second, third)
    )
    sphere_relation = first**2 + second**2 + third**2 - 1
    quotient, remainder = sp.div(
        sp.Poly(sp.expand(onb_sum - 96 * first**2 * second**2 * third**2), first),
        sp.Poly(sphere_relation, first),
    )
    assert remainder.as_expr() == 0
    assert quotient.as_expr() != 0


def check_onb_cone_obstruction() -> None:
    x, y, z = sp.symbols("x y z")
    radius_squared = x**2 + y**2 + z**2
    robinson = (
        x**6
        + y**6
        + z**6
        - x**4 * y**2
        - x**4 * z**2
        - y**4 * x**2
        - y**4 * z**2
        - z**4 * x**2
        - z**4 * y**2
        + 3 * x**2 * y**2 * z**2
    )

    def laplacian(polynomial: sp.Expr) -> sp.Expr:
        return sum(sp.diff(polynomial, variable, 2) for variable in (x, y, z))

    assert sp.expand(laplacian(laplacian(robinson)) - 240 * radius_squared) == 0

    # The spherical average follows from E[x^6]=1/7,
    # E[x^4 y^2]=1/35, and E[x^2 y^2 z^2]=1/105.
    spherical_average = (
        3 * sp.Rational(1, 7)
        - 6 * sp.Rational(1, 35)
        + 3 * sp.Rational(1, 105)
    )
    assert spherical_average == sp.Rational(2, 7)
    assert sp.Rational(16, 15) * spherical_average == sp.Rational(32, 105)

    points: list[list[int]] = []
    for zero_coordinate in range(3):
        other = [index for index in range(3) if index != zero_coordinate]
        for sign in (1, -1):
            vector = [0, 0, 0]
            vector[other[0]] = 1
            vector[other[1]] = sign
            points.append(vector)
    for first, second in product((1, -1), repeat=2):
        points.append([1, first, second])

    cubic_exponents = [
        (first, second, 3 - first - second)
        for first in range(4)
        for second in range(4 - first)
    ]
    evaluation = sp.Matrix(
        [
            [
                point[0] ** first * point[1] ** second * point[2] ** third
                for first, second, third in cubic_exponents
            ]
            for point in points
        ]
    )
    assert evaluation.det() == -128

    harmonic_cubics = sp.Matrix(
        [
            z * (5 * z**2 - 3 * radius_squared),
            x * (5 * z**2 - radius_squared),
            y * (5 * z**2 - radius_squared),
            z * (x**2 - y**2),
            2 * x * y * z,
            x * (x**2 - 3 * y**2),
            y * (3 * x**2 - y**2),
        ]
    )
    gram_variables = sp.symbols("robinson_gram_0:28")
    gram = sp.zeros(7)
    variable_index = 0
    for row in range(7):
        for column in range(row, 7):
            gram[row, column] = gram[column, row] = gram_variables[variable_index]
            variable_index += 1
    gram_equations = sp.Poly(
        sp.expand((harmonic_cubics.T * gram * harmonic_cubics)[0] - robinson),
        x,
        y,
        z,
    ).coeffs()
    gram_solution = sp.solve(gram_equations, gram_variables, dict=True)
    assert len(gram_solution) == 1
    robinson_gram = gram.subs(gram_solution[0])
    assert robinson_gram[4, 4] == -3

    weights = sp.symbols("cone_weight_0:10")
    scale = sp.symbols("cone_scale")
    energy = sp.Rational(32, 105) - sp.Rational(2, 7) * scale
    potential = 0
    for index, (point, norm) in enumerate(
        [(point, 2) for point in points[:6]]
        + [(point, 3) for point in points[6:]]
    ):
        dot = point[0] * x + point[1] * y + point[2] * z
        potential += weights[index] * (
            32 * dot**6 / norm**3
            - 48 * radius_squared * dot**4 / norm**2
            + 20 * radius_squared**2 * dot**2 / norm
            - sp.Rational(4, 3) * radius_squared**3
        )
    coefficient_equations = [
        coefficient
        for _, coefficient in sp.Poly(
            sp.expand(
                potential - energy * radius_squared**3 - scale * robinson
            ),
            x,
            y,
            z,
        ).terms()
    ]
    solutions = sp.solve(
        coefficient_equations + [sum(weights) - 1],
        [*weights, scale],
        dict=True,
    )
    assert solutions == [
        {
            scale: sp.Rational(64, 151),
            **{
                weights[index]: sp.Rational(184, 2265)
                for index in range(6)
            },
            **{
                weights[index]: sp.Rational(387, 3020)
                for index in range(6, 10)
            },
        }
    ]


def main() -> None:
    check_robinson()
    check_latitude()
    check_spin_four_and_onb_identities()
    check_onb_cone_obstruction()
    print("all exact KKT classification checks passed")


if __name__ == "__main__":
    main()
