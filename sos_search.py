#!/usr/bin/env python3
"""MOSEK search for KKT-infused flag/SOS certificates.

Independent samples from an antipodally symmetric measure on S^2 are
represented by multigraph monomials in their pairwise Gram entries.  No
isotropy assumption is made: every expectation with even vertex degrees is
retained as an independent canonical moment label, including second-moment
quantities such as E[(X.Y)^2].  The deficit from isotropy is expressed inside
the hierarchy by the degree-two harmonic flag square E[P_2(X.Y)] >= 0.  The
code therefore makes no unrecorded moment-closure assumption.

The hierarchy combines ordinary rooted flag squares of several arities with
the measure KKT conditions: global potential gaps, support stationarity, and
positive semidefiniteness of the rooted spherical Hessian.  Four-vector Gram
rank identities and their matrix-valued higher-arity forms are included as
exact equality relations.

All solver output is numerical.  A proof requires exact rational PSD matrices
and exact coefficient verification; this script currently finds only a
near-certificate.
"""

from __future__ import annotations

import argparse
import json
import itertools
import math
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import cvxpy as cp
import numpy as np
import sympy as sp

Exponent = tuple[int, int, int]
Label = tuple[object, ...]
Polynomial = list[tuple[Fraction, Exponent]]
GraphExponent = tuple[int, ...]
GraphPolynomial = list[tuple[Fraction, GraphExponent]]


def monomials(max_degree: int) -> list[Exponent]:
    if max_degree < 0:
        return []
    return [
        (i, j, total - i - j)
        for total in range(max_degree + 1)
        for i in range(total + 1)
        for j in range(total - i + 1)
    ]


def pair_label(power: int, factor: Fraction = Fraction(1)) -> tuple[Label, Fraction]:
    if power == 0:
        return ("constant",), factor
    return ("pair", power), factor


def expectation_label(exponent: Exponent) -> tuple[Label | None, Fraction]:
    """Reduce one triangle monomial using antipodality alone.

    The degree at X, Y, Z is respectively i+k, i+j, j+k.  Odd vertex degree
    has zero expectation.  Vertex degree zero reduces to a pair moment.  All
    remaining triangle moments are retained as independent labels, modulo
    S_3 symmetry.  No isotropy contraction is applied.
    """

    i, j, k = exponent
    vertex_degrees = i + k, i + j, j + k
    if any(degree % 2 for degree in vertex_degrees):
        return None, Fraction(0)

    if i + k == 0:
        return pair_label(j)
    if i + j == 0:
        return pair_label(k)
    if j + k == 0:
        return pair_label(i)

    return ("triangle",) + tuple(sorted(exponent)), Fraction(1)


def graph_edges(vertex_count: int) -> list[tuple[int, int]]:
    return [
        (left, right)
        for left in range(vertex_count)
        for right in range(left + 1, vertex_count)
    ]


def graph_matrix(
    vertex_count: int,
    edge_exponents: GraphExponent,
) -> list[list[int]]:
    matrix = [[0 for _ in range(vertex_count)] for _ in range(vertex_count)]
    for exponent, (left, right) in zip(
        edge_exponents,
        graph_edges(vertex_count),
        strict=True,
    ):
        matrix[left][right] = exponent
        matrix[right][left] = exponent
    return matrix


def matrix_edge_exponents(matrix: list[list[int]]) -> GraphExponent:
    return tuple(
        matrix[left][right]
        for left, right in graph_edges(len(matrix))
    )


def induced_matrix(
    matrix: list[list[int]],
    vertices: list[int],
) -> list[list[int]]:
    return [
        [matrix[left][right] for right in vertices]
        for left in vertices
    ]


def connected_components(matrix: list[list[int]]) -> list[list[int]]:
    remaining = set(range(len(matrix)))
    components: list[list[int]] = []
    while remaining:
        start = min(remaining)
        stack = [start]
        component: list[int] = []
        remaining.remove(start)
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            neighbors = {
                other
                for other, exponent in enumerate(matrix[vertex])
                if exponent and other in remaining
            }
            remaining.difference_update(neighbors)
            stack.extend(neighbors)
        components.append(sorted(component))
    return components


def canonical_connected_label(matrix: list[list[int]]) -> Label:
    vertex_count = len(matrix)
    canonical = min(
        tuple(
            matrix[permutation[left]][permutation[right]]
            for left, right in graph_edges(vertex_count)
        )
        for permutation in itertools.permutations(range(vertex_count))
    )
    if vertex_count == 2:
        return ("pair", canonical[0])
    if vertex_count == 3:
        return ("triangle",) + canonical
    return (f"graph_{vertex_count}",) + canonical


def canonical_graph_exponent(edge_exponents: GraphExponent) -> GraphExponent:
    matrix = graph_matrix(4, edge_exponents)
    return min(
        tuple(
            matrix[permutation[left]][permutation[right]]
            for left, right in graph_edges(4)
        )
        for permutation in itertools.permutations(range(4))
    )


def combine_component_labels(labels: list[Label]) -> Label:
    nonconstant = [label for label in labels if label != ("constant",)]
    if not nonconstant:
        return ("constant",)
    if len(nonconstant) == 1:
        return nonconstant[0]
    return ("product",) + tuple(sorted(nonconstant, key=str))


def multiply_labels(left: Label, right: Label) -> Label:
    def factors(label: Label) -> list[Label]:
        if label == ("constant",):
            return []
        if label and label[0] == "product":
            return list(label[1:])  # type: ignore[list-item]
        return [label]

    return combine_component_labels(factors(left) + factors(right))


def onb_label_value(label: Label) -> Fraction:
    """Evaluate a formal moment label on the uniform three-axis measure."""

    if label == ("constant",):
        return Fraction(1)
    if label[0] == "pair":
        return Fraction(1, 3)
    if label[0] == "product":
        value = Fraction(1)
        for factor in label[1:]:
            value *= onb_label_value(factor)  # type: ignore[arg-type]
        return value
    if label[0] == "triangle":
        vertex_count = 3
    elif isinstance(label[0], str) and label[0].startswith("graph_"):
        vertex_count = int(label[0].split("_")[1])
    else:
        raise ValueError(f"Unsupported ONB label: {label}")

    matrix = graph_matrix(
        vertex_count,
        tuple(int(value) for value in label[1:]),
    )
    components = connected_components(
        [
            [int(exponent > 0) for exponent in row]
            for row in matrix
        ]
    )
    return Fraction(1, 3 ** (vertex_count - len(components)))


@lru_cache(maxsize=None)
def equatorial_graph_moment(
    edge_exponents: GraphExponent,
    vertex_count: int,
    regular_order: int = 0,
) -> Fraction:
    """Exact moment of a graph monomial on a projective equator.

    ``regular_order == 0`` means Haar-uniform angle.  A positive order means
    the uniform measure on angles ``j*pi/regular_order``.  Expanding every
    cosine into its two Fourier characters reduces the integral to an integer
    circulation count.
    """

    flows: dict[tuple[int, ...], int] = {
        (0,) * vertex_count: 1
    }
    total_degree = 0
    for exponent, (left, right) in zip(
        edge_exponents,
        graph_edges(vertex_count),
        strict=True,
    ):
        if exponent == 0:
            continue
        total_degree += exponent
        next_flows: dict[tuple[int, ...], int] = {}
        for flow, weight in flows.items():
            for choice in range(exponent + 1):
                frequency = exponent - 2 * choice
                shifted = list(flow)
                shifted[left] += frequency
                shifted[right] -= frequency
                shifted_tuple = tuple(shifted)
                next_flows[shifted_tuple] = (
                    next_flows.get(shifted_tuple, 0)
                    + weight * math.comb(exponent, choice)
                )
        flows = next_flows

    if regular_order:
        modulus = 2 * regular_order
        numerator = sum(
            weight
            for flow, weight in flows.items()
            if all(value % modulus == 0 for value in flow)
        )
    else:
        numerator = flows.get((0,) * vertex_count, 0)
    return Fraction(numerator, 2**total_degree)


@lru_cache(maxsize=None)
def pole_equator_label_value(
    label: Label,
    regular_order: int = 0,
) -> Fraction:
    """Evaluate a moment on the 1/3-pole, 2/3-equator equality measure."""

    if label == ("constant",):
        return Fraction(1)
    if label[0] == "product":
        value = Fraction(1)
        for factor in label[1:]:
            value *= pole_equator_label_value(
                factor,  # type: ignore[arg-type]
                regular_order,
            )
        return value
    if label[0] == "pair":
        vertex_count = 2
    elif label[0] == "triangle":
        vertex_count = 3
    elif isinstance(label[0], str) and label[0].startswith("graph_"):
        vertex_count = int(label[0].split("_")[1])
    else:
        raise ValueError(f"Unsupported pole-equator label: {label}")

    edge_exponents = tuple(int(value) for value in label[1:])
    equator_moment = equatorial_graph_moment(
        edge_exponents,
        vertex_count,
        regular_order,
    )
    return (
        Fraction(1, 3**vertex_count)
        + Fraction(2, 3) ** vertex_count * equator_moment
    )


def rationalize_float(
    value: float,
    maximum_denominator: int = 10**9,
    tolerance: float = 1e-10,
) -> Fraction:
    """Recover the exact small rational used to construct a float."""

    rational = Fraction(float(value)).limit_denominator(maximum_denominator)
    if abs(float(rational) - float(value)) > tolerance:
        raise ValueError(f"Could not rationalize coefficient {value}")
    return rational


def exact_moment_matrix(
    label_matrices: dict[Label, np.ndarray],
    label_value,
) -> sp.Matrix:
    """Reconstruct a block moment matrix over the rationals."""

    size = next(iter(label_matrices.values())).shape[0]
    entries = [
        [Fraction(0) for _ in range(size)]
        for _ in range(size)
    ]
    for label, coefficient_matrix in label_matrices.items():
        moment = label_value(label)
        rows, columns = np.nonzero(np.abs(coefficient_matrix) > 1e-13)
        for row, column in zip(rows, columns, strict=True):
            entries[int(row)][int(column)] += (
                moment
                * rationalize_float(coefficient_matrix[row, column])
            )
    return sp.Matrix(
        [
            [
                sp.Rational(value.numerator, value.denominator)
                for value in row
            ]
            for row in entries
        ]
    )


def exact_onb_moment_matrix(
    label_matrices: dict[Label, np.ndarray],
) -> sp.Matrix:
    """Reconstruct a block's ONB moment matrix over the rationals."""

    return exact_moment_matrix(label_matrices, onb_label_value)


def exact_nullspace(moment_matrix: sp.Matrix) -> sp.Matrix:
    """Return a bounded rational basis for an exact moment-matrix kernel."""

    if moment_matrix.is_zero_matrix:
        return sp.eye(moment_matrix.rows)

    vectors = moment_matrix.nullspace(simplify=False)
    if not vectors:
        return sp.zeros(moment_matrix.rows, 0)

    primitive_vectors: list[sp.Matrix] = []
    for vector in vectors:
        denominators = [
            int(sp.denom(entry))
            for entry in vector
            if entry
        ]
        common_denominator = math.lcm(*denominators) if denominators else 1
        integers = [
            int(entry * common_denominator)
            for entry in vector
        ]
        divisor = math.gcd(*(abs(value) for value in integers if value)) or 1
        integers = [value // divisor for value in integers]
        scale = max(abs(value) for value in integers) or 1
        primitive_vectors.append(
            sp.Matrix(
                [
                    sp.Rational(value, scale)
                    for value in integers
                ]
            )
        )

    basis = sp.Matrix.hstack(*primitive_vectors)
    if moment_matrix * basis != sp.zeros(moment_matrix.rows, basis.cols):
        raise ValueError("Exact nullspace reconstruction failed")
    return basis


def exact_onb_nullspace(
    label_matrices: dict[Label, np.ndarray],
) -> sp.Matrix:
    """Return a bounded rational basis for the exact ONB kernel."""

    moment_matrix = exact_onb_moment_matrix(label_matrices)
    return exact_nullspace(moment_matrix)


def symmetric_matrix_generators(
    label_matrices: dict[Label, np.ndarray],
) -> list[dict[Label, Fraction]]:
    """Flatten a free symmetric matrix multiplier into exact columns.

    ``sum(M * X)`` counts an off-diagonal entry twice when both ``M`` and
    ``X`` are symmetric, so each independent upper-triangular coordinate of
    ``X`` has coefficient ``M_ij + M_ji``.
    """

    size = next(iter(label_matrices.values())).shape[0]
    generators: list[dict[Label, Fraction]] = []
    for row in range(size):
        for column in range(row, size):
            generator: dict[Label, Fraction] = {}
            for label, matrix in label_matrices.items():
                value = matrix[row, column]
                if row != column:
                    value += matrix[column, row]
                if abs(value) > 1e-13:
                    generator[label] = rationalize_float(float(value))
            if generator:
                generators.append(generator)
    return generators


def exact_equality_quotient_rows(
    ordered_labels: list[Label],
    free_label_matrices: Iterable[dict[Label, np.ndarray]],
    relations: Iterable[dict[Label, float]],
) -> tuple[list[dict[Label, Fraction]], int, int]:
    """Return an exact sparse basis annihilating all free equality terms.

    If ``F`` has the unrestricted KKT/rank generators as columns, a
    coefficient residual ``r`` can be absorbed by those generators exactly
    when every returned row ``q`` satisfies ``q^T r = 0``.  The rows form a
    basis of ``ker(F^T)`` and are obtained by rational RREF, without a
    floating-point rank decision.
    """

    generators: list[dict[Label, Fraction]] = []
    for label_matrices in free_label_matrices:
        generators.extend(symmetric_matrix_generators(label_matrices))
    generators.extend(
        {
            label: rationalize_float(float(coefficient))
            for label, coefficient in relation.items()
            if abs(coefficient) > 1e-13
        }
        for relation in relations
    )
    generators = [generator for generator in generators if generator]

    label_indices = {
        label: index for index, label in enumerate(ordered_labels)
    }
    transpose = sp.MutableSparseMatrix(
        len(generators),
        len(ordered_labels),
        {},
    )
    for row, generator in enumerate(generators):
        for label, coefficient in generator.items():
            transpose[row, label_indices[label]] = sp.Rational(
                coefficient.numerator,
                coefficient.denominator,
            )

    reduced, pivot_columns = transpose.rref(simplify=False)
    pivot_set = set(pivot_columns)
    quotient_rows: list[dict[Label, Fraction]] = []
    for free_column in range(len(ordered_labels)):
        if free_column in pivot_set:
            continue
        row: dict[Label, Fraction] = {
            ordered_labels[free_column]: Fraction(1)
        }
        for pivot_row, pivot_column in enumerate(pivot_columns):
            value = -reduced[pivot_row, free_column]
            if value:
                row[ordered_labels[pivot_column]] = Fraction(
                    int(sp.numer(value)),
                    int(sp.denom(value)),
                )
        quotient_rows.append(row)

    # An exact internal check protects the certificate search from an
    # accidentally omitted factor of two in a free matrix block.
    for row in quotient_rows:
        for generator in generators:
            pairing = sum(
                coefficient * generator.get(label, Fraction(0))
                for label, coefficient in row.items()
            )
            if pairing:
                raise ValueError("Equality quotient construction failed")

    return quotient_rows, len(generators), len(pivot_columns)


def reduce_graph_matrix(matrix: list[list[int]]) -> tuple[Label | None, Fraction]:
    """Evaluate all consequences of antipodality for a graph moment.

    Second moments are NOT contracted: a degree-two sampled vertex is kept in
    the canonical label, so no isotropy assumption enters the reduction.
    """

    if not matrix:
        return ("constant",), Fraction(1)

    degrees = [sum(row) for row in matrix]
    if any(degree % 2 for degree in degrees):
        return None, Fraction(0)

    active = [vertex for vertex, degree in enumerate(degrees) if degree]
    if len(active) < len(matrix):
        if not active:
            return ("constant",), Fraction(1)
        return reduce_graph_matrix(induced_matrix(matrix, active))

    components = connected_components(matrix)
    if len(components) > 1:
        labels: list[Label] = []
        coefficient = Fraction(1)
        for component in components:
            label, component_coefficient = reduce_graph_matrix(
                induced_matrix(matrix, component)
            )
            if label is None:
                return None, Fraction(0)
            labels.append(label)
            coefficient *= component_coefficient
        return combine_component_labels(labels), coefficient

    return canonical_connected_label(matrix), Fraction(1)


@lru_cache(maxsize=None)
def graph_expectation_label(
    vertex_count: int,
    edge_exponents: GraphExponent,
) -> tuple[Label | None, Fraction]:
    return reduce_graph_matrix(graph_matrix(vertex_count, edge_exponents))


def graph_expectation_vector(
    vertex_count: int,
    polynomial: GraphPolynomial,
) -> dict[Label, float]:
    vector: dict[Label, float] = {}
    for coefficient, edge_exponents in polynomial:
        label, reduction_coefficient = graph_expectation_label(
            vertex_count,
            edge_exponents,
        )
        if label is None or reduction_coefficient == 0:
            continue
        vector[label] = vector.get(label, 0.0) + float(
            coefficient * reduction_coefficient
        )
    return {
        label: value
        for label, value in vector.items()
        if abs(value) > 1e-13
    }


def add_exponents(left: Exponent, right: Exponent) -> Exponent:
    return tuple(left[index] + right[index] for index in range(3))  # type: ignore[return-value]


def multiply_polynomials(left: Polynomial, right: Polynomial) -> Polynomial:
    return [
        (left_coefficient * right_coefficient, add_exponents(left_exponent, right_exponent))
        for left_coefficient, left_exponent in left
        for right_coefficient, right_exponent in right
    ]


def shifted_polynomial(polynomial: Polynomial, shift: Exponent) -> Polynomial:
    return [
        (coefficient, add_exponents(exponent, shift))
        for coefficient, exponent in polynomial
    ]


def expectation_matrix(basis: list[Exponent], multiplier: Polynomial) -> dict[Label, np.ndarray]:
    size = len(basis)
    matrices: dict[Label, np.ndarray] = {}
    for row, left in enumerate(basis):
        for column, right in enumerate(basis):
            product_exponent = add_exponents(left, right)
            for multiplier_coefficient, multiplier_exponent in multiplier:
                label, reduction_coefficient = expectation_label(
                    add_exponents(product_exponent, multiplier_exponent)
                )
                if label is None or reduction_coefficient == 0:
                    continue
                matrix = matrices.setdefault(label, np.zeros((size, size)))
                matrix[row, column] += float(
                    multiplier_coefficient * reduction_coefficient
                )
    return matrices


def expectation_vector(polynomial: Polynomial) -> dict[Label, float]:
    vector: dict[Label, float] = {}
    for coefficient, exponent in polynomial:
        label, reduction_coefficient = expectation_label(exponent)
        if label is None or reduction_coefficient == 0:
            continue
        vector[label] = vector.get(label, 0.0) + float(
            coefficient * reduction_coefficient
        )
    return {label: value for label, value in vector.items() if abs(value) > 1e-13}


def add_polynomials(left: Polynomial, right: Polynomial) -> Polynomial:
    coefficients: dict[Exponent, Fraction] = {}
    for coefficient, exponent in left + right:
        coefficients[exponent] = coefficients.get(exponent, Fraction(0)) + coefficient
    return [
        (coefficient, exponent)
        for exponent, coefficient in coefficients.items()
        if coefficient
    ]


def scale_polynomial(polynomial: Polynomial, factor: Fraction) -> Polynomial:
    return [(factor * coefficient, exponent) for coefficient, exponent in polynomial]


def sympy_polynomial(poly: sp.Poly) -> Polynomial:
    return [
        (Fraction(int(coefficient)), tuple(int(value) for value in exponent))
        for exponent, coefficient in poly.terms()
    ]


def kernel_polynomials() -> tuple[Polynomial, Polynomial, Polynomial]:
    a, b, c = sp.symbols("a b c")
    kernel_prime = 192 * a**5 - 192 * a**3 + 40 * a
    kernel_second = 960 * a**4 - 576 * a**2 + 40
    tangent_inner_product = b - a * c

    gradient = sp.Poly(
        sp.expand(kernel_prime * tangent_inner_product),
        a,
        b,
        c,
    )
    hessian = sp.Poly(
        sp.expand(
            kernel_second * tangent_inner_product**2
            - kernel_prime * a * (1 - c**2)
        ),
        a,
        b,
        c,
    )
    gram_determinant = (
        1 + 2 * a * b * c - a**2 - b**2 - c**2
    )
    perpendicular_hessian = sp.Poly(
        sp.expand(
            kernel_second * gram_determinant
            - kernel_prime * a * (1 - c**2)
        ),
        a,
        b,
        c,
    )
    return (
        sympy_polynomial(gradient),
        sympy_polynomial(hessian),
        sympy_polynomial(perpendicular_hessian),
    )


def four_point_hessian_polynomials() -> tuple[GraphPolynomial, GraphPolynomial]:
    """Bilinear rooted Hessian kernels for two auxiliary tangent vectors.

    The six edge variables are ordered as

      XY=a, XZ=c, XW=d, YZ=b, YW=e, ZW=f.
    """

    a, c, d, b, e, f = sp.symbols("a c d b e f")
    kernel_prime = 192 * a**5 - 192 * a**3 + 40 * a
    kernel_second = 960 * a**4 - 576 * a**2 + 40
    yz_tangent = b - a * c
    yw_tangent = e - a * d
    zw_tangent = f - c * d

    parallel = sp.Poly(
        sp.expand(
            kernel_second * yz_tangent * yw_tangent
            - kernel_prime * a * zw_tangent
        ),
        a,
        c,
        d,
        b,
        e,
        f,
    )
    perpendicular_y_product = (
        zw_tangent * (1 - a**2) - yz_tangent * yw_tangent
    )
    perpendicular = sp.Poly(
        sp.expand(
            kernel_second * perpendicular_y_product
            - kernel_prime * a * zw_tangent
        ),
        a,
        c,
        d,
        b,
        e,
        f,
    )

    def convert(poly: sp.Poly) -> GraphPolynomial:
        return [
            (
                Fraction(int(coefficient)),
                tuple(int(value) for value in exponent),
            )
            for exponent, coefficient in poly.terms()
        ]

    return convert(parallel), convert(perpendicular)


def global_tangent_gap_polynomials() -> tuple[Polynomial, Polynomial]:
    """Cleared-denominator KKT gaps at two tangent trial points."""

    a, b, c = sp.symbols("a b c")
    tangent_norm_squared = 1 - c**2
    parallel_numerator = b - a * c
    gram_determinant = 1 + 2 * a * b * c - a**2 - b**2 - c**2

    def kernel_from_squared_inner_product(value: sp.Expr) -> sp.Expr:
        return 32 * value**3 - 48 * value**2 + 20 * value - sp.Rational(4, 3)

    root_kernel = kernel_from_squared_inner_product(a**2)
    parallel = sp.Poly(
        sp.expand(
            32 * parallel_numerator**6
            - 48 * parallel_numerator**4 * tangent_norm_squared
            + 20 * parallel_numerator**2 * tangent_norm_squared**2
            - sp.Rational(4, 3) * tangent_norm_squared**3
            - tangent_norm_squared**3 * root_kernel
        ),
        a,
        b,
        c,
    )
    perpendicular = sp.Poly(
        sp.expand(
            32 * gram_determinant**3
            - 48 * gram_determinant**2 * tangent_norm_squared
            + 20 * gram_determinant * tangent_norm_squared**2
            - sp.Rational(4, 3) * tangent_norm_squared**3
            - tangent_norm_squared**3 * root_kernel
        ),
        a,
        b,
        c,
    )
    return sympy_polynomial(parallel), sympy_polynomial(perpendicular)


def four_point_rank_relations(maximum_multiplier_degree: int) -> list[dict[Label, float]]:
    variables = sp.symbols("g01 g02 g03 g12 g13 g23")
    gram = sp.Matrix(
        [
            [1, variables[0], variables[1], variables[2]],
            [variables[0], 1, variables[3], variables[4]],
            [variables[1], variables[3], 1, variables[5]],
            [variables[2], variables[4], variables[5], 1],
        ]
    )
    determinant = sp.Poly(sp.expand(gram.det()), *variables)
    determinant_terms: GraphPolynomial = [
        (
            Fraction(int(coefficient)),
            tuple(int(value) for value in exponent),
        )
        for exponent, coefficient in determinant.terms()
    ]

    representatives: set[GraphExponent] = set()
    for total_degree in range(maximum_multiplier_degree + 1):
        for first in range(total_degree + 1):
            for second in range(total_degree - first + 1):
                for third in range(total_degree - first - second + 1):
                    for fourth in range(
                        total_degree - first - second - third + 1
                    ):
                        for fifth in range(
                            total_degree
                            - first
                            - second
                            - third
                            - fourth
                            + 1
                        ):
                            sixth = (
                                total_degree
                                - first
                                - second
                                - third
                                - fourth
                                - fifth
                            )
                            exponent = (
                                first,
                                second,
                                third,
                                fourth,
                                fifth,
                                sixth,
                            )
                            degrees = [
                                first + second + third,
                                first + fourth + fifth,
                                second + fourth + sixth,
                                third + fifth + sixth,
                            ]
                            if any(degree % 2 for degree in degrees):
                                continue
                            representatives.add(canonical_graph_exponent(exponent))

    relations: list[dict[Label, float]] = []
    normalized_relations: set[tuple[tuple[str, float], ...]] = set()
    for multiplier in sorted(representatives):
        polynomial = [
            (
                coefficient,
                tuple(
                    exponent[index] + multiplier[index]
                    for index in range(6)
                ),
            )
            for coefficient, exponent in determinant_terms
        ]
        relation = graph_expectation_vector(4, polynomial)
        if not relation:
            continue
        first_label = min(relation, key=str)
        scale = relation[first_label]
        normalized = tuple(
            sorted(
                (
                    (str(label), round(value / scale, 11))
                    for label, value in relation.items()
                )
            )
        )
        if normalized in normalized_relations:
            continue
        normalized_relations.add(normalized)
        relations.append(relation)
    return relations


def four_point_hessian_expectation_matrix(
    auxiliary_degrees: list[int],
    hessian_polynomial: GraphPolynomial,
) -> dict[Label, np.ndarray]:
    size = len(auxiliary_degrees)
    matrices: dict[Label, np.ndarray] = {}
    for row, left_degree in enumerate(auxiliary_degrees):
        for column, right_degree in enumerate(auxiliary_degrees):
            for coefficient, exponent in hessian_polynomial:
                shifted = list(exponent)
                shifted[1] += left_degree
                shifted[2] += right_degree
                label, reduction_coefficient = graph_expectation_label(
                    4,
                    tuple(shifted),
                )
                if label is None or reduction_coefficient == 0:
                    continue
                matrix = matrices.setdefault(label, np.zeros((size, size)))
                matrix[row, column] += float(coefficient * reduction_coefficient)
    return matrices


def potential_stationarity_relation(power: int) -> dict[Label, float]:
    """Return E[K(X.Y)(X.Z)^r] - E[K(X.Y)]E[(Z.W)^r]."""

    kernel_terms = [
        (Fraction(-4, 3), 0),
        (Fraction(20), 2),
        (Fraction(-48), 4),
        (Fraction(32), 6),
    ]
    triangle: GraphPolynomial = [
        (coefficient, (kernel_power, power, 0))
        for coefficient, kernel_power in kernel_terms
    ]
    disconnected_pairs: GraphPolynomial = [
        (coefficient, (kernel_power, 0, 0, 0, 0, power))
        for coefficient, kernel_power in kernel_terms
    ]
    relation = graph_expectation_vector(3, triangle)
    for label, value in graph_expectation_vector(4, disconnected_pairs).items():
        relation[label] = relation.get(label, 0.0) - value
    return {
        label: value
        for label, value in relation.items()
        if abs(value) > 1e-12
    }


def empty_type_flag_expectation_matrix(
    pair_degrees: list[int],
) -> dict[Label, np.ndarray]:
    """Gram block for squares of unrooted two-sample flag averages."""

    size = len(pair_degrees)
    matrices: dict[Label, np.ndarray] = {}
    for row, left_degree in enumerate(pair_degrees):
        for column, right_degree in enumerate(pair_degrees):
            label, reduction_coefficient = graph_expectation_label(
                4,
                (left_degree, 0, 0, 0, 0, right_degree),
            )
            if label is None or reduction_coefficient == 0:
                continue
            matrix = matrices.setdefault(label, np.zeros((size, size)))
            matrix[row, column] += float(reduction_coefficient)
    return matrices


def harmonic_flag_expectation_matrix(
    order: int,
    weight_degrees: list[int],
) -> dict[Label, np.ndarray]:
    """Spin-``order`` Gram block of harmonic-weighted unrooted pair flags.

    For v_a = int int Y_{lm}(x) (x.y)^a dmu(x) dmu(y), summing the Gram
    matrix over m gives

        G_ab = E[P_l(X.Z) (X.Y)^a (Z.W)^b] >= 0 (as a matrix),

    a four-point unrooted flag square.  For l = 2 this block carries the
    deviatoric second-moment correlations that an isotropy assumption would
    otherwise fix, so it subsumes the isotropy reduction into the flag
    decomposition.
    """

    variable = sp.symbols("t")
    legendre_terms = [
        (Fraction(int(coefficient.p), int(coefficient.q)), int(power))
        for (power,), coefficient in sp.Poly(
            sp.legendre(order, variable),
            variable,
        ).terms()
    ]
    size = len(weight_degrees)
    matrices: dict[Label, np.ndarray] = {}
    # Vertices: 0 = X, 1 = Z, 2 = Y, 3 = W.  Edges in graph_edges(4) order:
    # (0,1), (0,2), (0,3), (1,2), (1,3), (2,3).
    for row, left_degree in enumerate(weight_degrees):
        for column, right_degree in enumerate(weight_degrees):
            for coefficient, power in legendre_terms:
                label, reduction_coefficient = graph_expectation_label(
                    4,
                    (power, left_degree, 0, 0, right_degree, 0),
                )
                if label is None or reduction_coefficient == 0:
                    continue
                matrix = matrices.setdefault(label, np.zeros((size, size)))
                matrix[row, column] += float(
                    coefficient * reduction_coefficient
                )
    return {
        label: matrix
        for label, matrix in matrices.items()
        if np.max(np.abs(matrix)) > 1e-13
    }


def spin2_flag_basis(
    maximum_flag_degree: int,
) -> list[tuple[tuple[int, ...], int]]:
    """Unrooted spin-2 flags: a harmonic vertex with up to two leaves.

    A flag ``((a_1, ..., a_k), b)`` denotes

        v_m = int Y_2m(x) prod_i (x.y_i)^{a_i} (y_1.y_2)^b dmu^{1+k},

    with ``b = 0`` unless ``k = 2``.  Antipodal parity requires every leaf
    degree ``a_i + b`` to be even and ``sum a_i`` to be even.  The empty flag
    is the deviatoric second moment D itself, whose squared norm is a
    positive multiple of the harmonic energy h_2.
    """

    flags: list[tuple[tuple[int, ...], int]] = [((), 0)]
    for a in range(2, maximum_flag_degree + 1, 2):
        flags.append(((a,), 0))
    for a1 in range(1, maximum_flag_degree + 1):
        for a2 in range(a1, maximum_flag_degree + 1):
            if (a1 + a2) % 2:
                continue
            for b in range(maximum_flag_degree + 1):
                if (a1 + b) % 2 or (a2 + b) % 2:
                    continue
                if a1 + a2 + b > maximum_flag_degree:
                    continue
                flags.append(((a1, a2), b))
    return flags


def spin2_flag_expectation_matrix(
    flag_basis: list[tuple[tuple[int, ...], int]],
) -> dict[Label, np.ndarray]:
    """Gram block of unrooted spin-2 flags.

    Summing v_m(A) conj(v_m(B)) over m gives (up to a positive constant)

        G_AB = E[P_2(X.X') flag_A(X;Y) flag_B(X';Y')],

    a fully unrooted flag square on independent samples.  At a feasible
    point with h_2 = 0 positive semidefiniteness forces the entire D-row to
    vanish, which restores every isotropic contraction identity whose
    residual flag lies in the basis.  This subsumes the isotropy reduction
    into the flag decomposition.
    """

    legendre_terms = [
        (Fraction(3, 2), 2),
        (Fraction(-1, 2), 0),
    ]
    size = len(flag_basis)
    matrices: dict[Label, np.ndarray] = {}
    for row, (left_leaves, left_pair) in enumerate(flag_basis):
        for column, (right_leaves, right_pair) in enumerate(flag_basis):
            left_count = len(left_leaves)
            right_count = len(right_leaves)
            vertex_count = 2 + left_count + right_count
            left_root = 0
            right_root = 1
            left_offset = 2
            right_offset = 2 + left_count
            edge_indices = {
                edge: index
                for index, edge in enumerate(graph_edges(vertex_count))
            }
            base = [0] * len(edge_indices)
            for index, power in enumerate(left_leaves):
                base[edge_indices[(left_root, left_offset + index)]] += power
            if left_pair:
                base[
                    edge_indices[(left_offset, left_offset + 1)]
                ] += left_pair
            for index, power in enumerate(right_leaves):
                base[
                    edge_indices[(right_root, right_offset + index)]
                ] += power
            if right_pair:
                base[
                    edge_indices[(right_offset, right_offset + 1)]
                ] += right_pair
            for coefficient, power in legendre_terms:
                exponent = list(base)
                exponent[edge_indices[(left_root, right_root)]] += power
                label, reduction_coefficient = graph_expectation_label(
                    vertex_count,
                    tuple(exponent),
                )
                if label is None or reduction_coefficient == 0:
                    continue
                matrix = matrices.setdefault(label, np.zeros((size, size)))
                matrix[row, column] += float(
                    coefficient * reduction_coefficient
                )
    return {
        label: matrix
        for label, matrix in matrices.items()
        if np.max(np.abs(matrix)) > 1e-13
    }


def h2_localized_flag_expectation_matrix(
    leaf_degrees: list[int],
    tangent_harmonic: Polynomial,
) -> dict[Label, np.ndarray]:
    """One-root flag Gram block localized by the harmonic energy h_2.

    Multiplying a conditional flag square by the nonnegative moment
    quantity h_2 = E[P_2(W.W')] over two fresh independent samples yields

        h_2 * E[flag square] = E[P_2(W.W') flag(X;Y,Z)^2] >= 0,

    a five-vertex block whose entries are product labels coupling the
    spin-2 deficit to three-point flag moments.  These are the
    second-order h_2 couplings exploited by the sqrt(h_2) pseudo-moment
    leak, so this localization is the direct certificate-side response.
    """

    legendre_terms = [
        (Fraction(3, 2), 2),
        (Fraction(-1, 2), 0),
    ]
    edge_indices = {
        edge: index for index, edge in enumerate(graph_edges(5))
    }
    size = len(leaf_degrees)
    matrices: dict[Label, np.ndarray] = {}
    for row, left_degree in enumerate(leaf_degrees):
        for column, right_degree in enumerate(leaf_degrees):
            for coefficient, (i, j, k) in tangent_harmonic:
                shifted = (i + left_degree, j, k + right_degree)
                for legendre_coefficient, legendre_power in legendre_terms:
                    exponent = [0] * len(edge_indices)
                    exponent[edge_indices[(0, 1)]] = shifted[0]
                    exponent[edge_indices[(1, 2)]] = shifted[1]
                    exponent[edge_indices[(0, 2)]] = shifted[2]
                    exponent[edge_indices[(3, 4)]] = legendre_power
                    label, reduction_coefficient = graph_expectation_label(
                        5,
                        tuple(exponent),
                    )
                    if label is None or reduction_coefficient == 0:
                        continue
                    matrix = matrices.setdefault(
                        label,
                        np.zeros((size, size)),
                    )
                    matrix[row, column] += float(
                        coefficient
                        * legendre_coefficient
                        * reduction_coefficient
                    )
    return {
        label: matrix
        for label, matrix in matrices.items()
        if np.max(np.abs(matrix)) > 1e-13
    }


def harmonic_pair_vector(degree: int) -> dict[Label, np.ndarray]:
    variable = sp.symbols("t")
    polynomial = sp.Poly(sp.legendre(degree, variable), variable)
    vector: dict[Label, float] = {}
    for (power,), coefficient in polynomial.terms():
        label, reduction_coefficient = pair_label(power)
        vector[label] = vector.get(label, 0.0) + float(
            Fraction(int(coefficient.p), int(coefficient.q))
            * reduction_coefficient
        )
    return {
        label: np.array([[value]])
        for label, value in vector.items()
        if abs(value) > 1e-13
    }


def potential_flag_relation_matrix(
    leaf_degrees: list[int],
    tangent_harmonic: Polynomial,
) -> dict[Label, np.ndarray]:
    """Matrix form of (U(X)-E) times a one-root flag Gram matrix."""

    size = len(leaf_degrees)
    matrices: dict[Label, np.ndarray] = {}
    kernel_terms = [
        (Fraction(-4, 3), 0),
        (Fraction(20), 2),
        (Fraction(-48), 4),
        (Fraction(32), 6),
    ]
    energy_terms = [
        (Fraction(-4, 3), ("constant",)),
        (Fraction(20), ("pair", 2)),
        (Fraction(-48), ("pair", 4)),
        (Fraction(32), ("pair", 6)),
    ]

    def add_entry(label: Label, row: int, column: int, value: Fraction) -> None:
        matrix = matrices.setdefault(label, np.zeros((size, size)))
        matrix[row, column] += float(value)

    for row, left_degree in enumerate(leaf_degrees):
        for column, right_degree in enumerate(leaf_degrees):
            # E[K(X.Y) flag(X;Z,W)]
            for kernel_coefficient, kernel_power in kernel_terms:
                for flag_coefficient, (xz_power, zw_power, xw_power) in tangent_harmonic:
                    label, reduction_coefficient = graph_expectation_label(
                        4,
                        (
                            kernel_power,
                            xz_power + left_degree,
                            xw_power + right_degree,
                            0,
                            0,
                            zw_power,
                        ),
                    )
                    if label is not None and reduction_coefficient:
                        add_entry(
                            label,
                            row,
                            column,
                            kernel_coefficient
                            * flag_coefficient
                            * reduction_coefficient,
                        )

            # -E[K] E[flag(X;Z,W)]
            flag_vector: dict[Label, Fraction] = {}
            for flag_coefficient, (xz_power, zw_power, xw_power) in tangent_harmonic:
                label, reduction_coefficient = expectation_label(
                    (
                        xz_power + left_degree,
                        zw_power,
                        xw_power + right_degree,
                    )
                )
                if label is not None and reduction_coefficient:
                    flag_vector[label] = flag_vector.get(label, Fraction(0)) + (
                        flag_coefficient * reduction_coefficient
                    )
            for flag_label, flag_coefficient in flag_vector.items():
                for energy_coefficient, energy_label in energy_terms:
                    add_entry(
                        multiply_labels(energy_label, flag_label),
                        row,
                        column,
                        -energy_coefficient * flag_coefficient,
                    )

    return {
        label: matrix
        for label, matrix in matrices.items()
        if np.max(np.abs(matrix)) > 1e-13
    }


def two_root_flag_expectation_matrix(
    flag_basis: list[Exponent],
    leaf_pairing_multiplier: GraphPolynomial | None = None,
) -> dict[Label, np.ndarray]:
    """Gram block for squares after conditioning on two sampled roots.

    A basis exponent (i,j,k) denotes

      (X.Y)^i (Z.Y)^j (X.Z)^k,

    where X,Z are the roots and Y is the unlabeled leaf.
    """

    if leaf_pairing_multiplier is None:
        leaf_pairing_multiplier = [(Fraction(1), (0, 0, 0, 0, 0, 0))]

    size = len(flag_basis)
    matrices: dict[Label, np.ndarray] = {}
    for row, (xy_power, zy_power, xz_power) in enumerate(flag_basis):
        for column, (xw_power, zw_power, xz_power_right) in enumerate(flag_basis):
            base = (
                xz_power + xz_power_right,
                xy_power,
                xw_power,
                zy_power,
                zw_power,
                0,
            )
            for coefficient, shift in leaf_pairing_multiplier:
                label, reduction_coefficient = graph_expectation_label(
                    4,
                    tuple(
                        base[index] + shift[index]
                        for index in range(6)
                    ),
                )
                if label is None or reduction_coefficient == 0:
                    continue
                matrix = matrices.setdefault(label, np.zeros((size, size)))
                matrix[row, column] += float(
                    coefficient * reduction_coefficient
                )
    return matrices


@lru_cache(maxsize=None)
def exponent_tuples_up_to(
    length: int,
    maximum_total_degree: int,
) -> list[tuple[int, ...]]:
    """All nonnegative exponent tuples of a bounded total degree."""

    if length == 0:
        return [()] if maximum_total_degree >= 0 else []
    if maximum_total_degree < 0:
        return []
    tuples: list[tuple[int, ...]] = []
    for first in range(maximum_total_degree + 1):
        for tail in exponent_tuples_up_to(
            length - 1,
            maximum_total_degree - first,
        ):
            tuples.append((first,) + tail)
    return tuples


def rooted_star_flag_sectors(
    root_count: int,
    maximum_leaf_degree: int,
    leaf_parity: int,
) -> dict[tuple[int, ...], list[tuple[int, ...]]]:
    """Projective character sectors for a leaf joined to sampled roots.

    An exponent tuple ``e`` represents

        product_i (X_i.Y)^e_i.

    The parity at the leaf is ``sum(e) mod 2``.  Flags with the same root
    parity character may be paired in a conditional Gram block.
    """

    sectors: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for exponent in exponent_tuples_up_to(root_count, maximum_leaf_degree):
        if sum(exponent) % 2 != leaf_parity:
            continue
        signature = tuple(value % 2 for value in exponent)
        sectors.setdefault(signature, []).append(exponent)
    return sectors


def rooted_weighted_flag_sectors(
    root_count: int,
    maximum_total_degree: int,
    maximum_root_degree: int,
) -> dict[
    tuple[int, ...],
    list[tuple[tuple[int, ...], tuple[int, ...]]],
]:
    """Flag sectors allowing bounded polynomial factors among the roots."""

    root_edges = graph_edges(root_count)
    sectors: dict[
        tuple[int, ...],
        list[tuple[tuple[int, ...], tuple[int, ...]]],
    ] = {}
    for root_exponent in exponent_tuples_up_to(
        len(root_edges),
        maximum_root_degree,
    ):
        root_degree = sum(root_exponent)
        if root_degree > maximum_total_degree:
            continue
        for leaf_exponent in exponent_tuples_up_to(
            root_count,
            maximum_total_degree - root_degree,
        ):
            if sum(leaf_exponent) % 2:
                continue
            vertex_degrees = list(leaf_exponent)
            for power, (left, right) in zip(
                root_exponent,
                root_edges,
                strict=True,
            ):
                vertex_degrees[left] += power
                vertex_degrees[right] += power
            signature = tuple(degree % 2 for degree in vertex_degrees)
            sectors.setdefault(signature, []).append(
                (root_exponent, leaf_exponent)
            )
    return sectors


def rooted_weighted_flag_expectation_matrix(
    root_count: int,
    flag_basis: list[tuple[tuple[int, ...], tuple[int, ...]]],
) -> dict[Label, np.ndarray]:
    """Conditional squares with both root-root and root-leaf factors."""

    vertex_count = root_count + 2
    left_leaf = root_count
    right_leaf = root_count + 1
    target_edge_indices = {
        edge: index
        for index, edge in enumerate(graph_edges(vertex_count))
    }
    root_edges = graph_edges(root_count)
    size = len(flag_basis)
    matrices: dict[Label, np.ndarray] = {}
    for row, (left_root, left_leaf_exponent) in enumerate(flag_basis):
        for column, (right_root, right_leaf_exponent) in enumerate(flag_basis):
            exponent = [0] * len(target_edge_indices)
            for edge_index, edge in enumerate(root_edges):
                exponent[target_edge_indices[edge]] += (
                    left_root[edge_index] + right_root[edge_index]
                )
            for root, power in enumerate(left_leaf_exponent):
                exponent[target_edge_indices[(root, left_leaf)]] += power
            for root, power in enumerate(right_leaf_exponent):
                exponent[target_edge_indices[(root, right_leaf)]] += power
            label, reduction_coefficient = graph_expectation_label(
                vertex_count,
                tuple(exponent),
            )
            if label is None or reduction_coefficient == 0:
                continue
            matrix = matrices.setdefault(label, np.zeros((size, size)))
            matrix[row, column] += float(reduction_coefficient)
    return {
        label: matrix
        for label, matrix in matrices.items()
        if np.max(np.abs(matrix)) > 1e-13
    }


def rooted_star_flag_expectation_matrix(
    root_count: int,
    flag_basis: list[tuple[int, ...]],
    multiplier: GraphPolynomial | None = None,
) -> dict[Label, np.ndarray]:
    """Conditional flag-square block with ``root_count`` shared roots.

    Squaring a flag with one integrated leaf produces two independent leaves,
    hence a graph on ``root_count + 2`` sampled vertices.  This is the
    systematic five-/six-point extension of the selected two-root blocks.
    """

    vertex_count = root_count + 2
    if multiplier is None:
        multiplier = [
            (
                Fraction(1),
                tuple(0 for _ in graph_edges(vertex_count)),
            )
        ]
    left_leaf = root_count
    right_leaf = root_count + 1
    edge_indices = {
        edge: index
        for index, edge in enumerate(graph_edges(vertex_count))
    }
    size = len(flag_basis)
    matrices: dict[Label, np.ndarray] = {}
    for row, left_exponent in enumerate(flag_basis):
        for column, right_exponent in enumerate(flag_basis):
            graph_exponent = [0] * len(edge_indices)
            for root, power in enumerate(left_exponent):
                graph_exponent[edge_indices[(root, left_leaf)]] += power
            for root, power in enumerate(right_exponent):
                graph_exponent[edge_indices[(root, right_leaf)]] += power
            for coefficient, multiplier_exponent in multiplier:
                exponent = tuple(
                    graph_exponent[index] + multiplier_exponent[index]
                    for index in range(len(edge_indices))
                )
                label, reduction_coefficient = graph_expectation_label(
                    vertex_count,
                    exponent,
                )
                if label is None or reduction_coefficient == 0:
                    continue
                matrix = matrices.setdefault(label, np.zeros((size, size)))
                matrix[row, column] += float(
                    coefficient * reduction_coefficient
                )
    return {
        label: matrix
        for label, matrix in matrices.items()
        if np.max(np.abs(matrix)) > 1e-13
    }


def four_vertex_gram_determinant() -> GraphPolynomial:
    variables = sp.symbols("g01 g02 g03 g12 g13 g23")
    gram = sp.Matrix(
        [
            [1, variables[0], variables[1], variables[2]],
            [variables[0], 1, variables[3], variables[4]],
            [variables[1], variables[3], 1, variables[5]],
            [variables[2], variables[4], variables[5], 1],
        ]
    )
    determinant = sp.Poly(sp.expand(gram.det()), *variables)
    return [
        (
            Fraction(int(coefficient)),
            tuple(int(value) for value in exponent),
        )
        for exponent, coefficient in determinant.terms()
    ]


def lift_graph_polynomial(
    polynomial: GraphPolynomial,
    source_vertex_count: int,
    target_vertex_count: int,
    target_vertices: tuple[int, ...],
) -> GraphPolynomial:
    """Embed a graph polynomial into a larger labelled graph."""

    if len(target_vertices) != source_vertex_count:
        raise ValueError("Wrong number of target vertices")
    target_edge_indices = {
        edge: index
        for index, edge in enumerate(graph_edges(target_vertex_count))
    }
    result: GraphPolynomial = []
    for coefficient, source_exponent in polynomial:
        target_exponent = [0] * len(target_edge_indices)
        for power, (source_left, source_right) in zip(
            source_exponent,
            graph_edges(source_vertex_count),
            strict=True,
        ):
            target_edge = tuple(
                sorted(
                    (
                        target_vertices[source_left],
                        target_vertices[source_right],
                    )
                )
            )
            target_exponent[target_edge_indices[target_edge]] += power
        result.append((coefficient, tuple(target_exponent)))
    return result


def lifted_hessian_expectation_matrix(
    conditioning_root_count: int,
    flag_basis: list[tuple[int, ...]],
    hessian_polynomial: GraphPolynomial,
) -> dict[Label, np.ndarray]:
    """KKT Hessian block with extra shared conditioning roots.

    Root 0 is the Hessian base point.  The remaining roots are additional
    shared samples.  The Hessian integration point and the two vector-field
    leaves are independent, so the resulting arity is
    ``conditioning_root_count + 3``.

    A flag exponent ``e`` represents a tangent vector-field integrand

        product_i (X_i.Z)^e_i P_{X_0-perp} Z.

    Its total leaf degree must be odd so the integrand is projectively even.
    """

    vertex_count = conditioning_root_count + 3
    hessian_leaf = conditioning_root_count
    left_leaf = conditioning_root_count + 1
    right_leaf = conditioning_root_count + 2
    edge_indices = {
        edge: index
        for index, edge in enumerate(graph_edges(vertex_count))
    }

    # The four-point polynomial uses vertices (X,Y,Z,W).
    old_vertices = (0, hessian_leaf, left_leaf, right_leaf)
    lifted_terms: GraphPolynomial = []
    for coefficient, old_exponent in hessian_polynomial:
        lifted = [0] * len(edge_indices)
        for power, (old_left, old_right) in zip(
            old_exponent,
            graph_edges(4),
            strict=True,
        ):
            new_edge = tuple(
                sorted(
                    (
                        old_vertices[old_left],
                        old_vertices[old_right],
                    )
                )
            )
            lifted[edge_indices[new_edge]] += power
        lifted_terms.append((coefficient, tuple(lifted)))

    size = len(flag_basis)
    matrices: dict[Label, np.ndarray] = {}
    for row, left_exponent in enumerate(flag_basis):
        for column, right_exponent in enumerate(flag_basis):
            flag_shift = [0] * len(edge_indices)
            for root, power in enumerate(left_exponent):
                flag_shift[edge_indices[(root, left_leaf)]] += power
            for root, power in enumerate(right_exponent):
                flag_shift[edge_indices[(root, right_leaf)]] += power
            for coefficient, hessian_exponent in lifted_terms:
                exponent = tuple(
                    hessian_exponent[index] + flag_shift[index]
                    for index in range(len(edge_indices))
                )
                label, reduction_coefficient = graph_expectation_label(
                    vertex_count,
                    exponent,
                )
                if label is None or reduction_coefficient == 0:
                    continue
                matrix = matrices.setdefault(label, np.zeros((size, size)))
                matrix[row, column] += float(
                    coefficient * reduction_coefficient
                )
    return {
        label: matrix
        for label, matrix in matrices.items()
        if np.max(np.abs(matrix)) > 1e-13
    }


ORIENTATION_PAIRING: GraphPolynomial = [
    # det(X,Z,Y) det(X,Z,W)
    (Fraction(1), (0, 0, 0, 0, 0, 1)),
    (Fraction(-1), (2, 0, 0, 0, 0, 1)),
    (Fraction(-1), (0, 0, 0, 1, 1, 0)),
    (Fraction(-1), (0, 1, 1, 0, 0, 0)),
    (Fraction(1), (1, 1, 0, 0, 1, 0)),
    (Fraction(1), (1, 0, 1, 1, 0, 0)),
]
ROOT_PAIR_MINOR: GraphPolynomial = [
    (Fraction(1), (0, 0, 0, 0, 0, 0)),
    (Fraction(-1), (2, 0, 0, 0, 0, 0)),
]


def multiply_graph_polynomials(
    left: GraphPolynomial,
    right: GraphPolynomial,
) -> GraphPolynomial:
    return [
        (
            left_coefficient * right_coefficient,
            tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(len(left_exponent))
            ),
        )
        for left_coefficient, left_exponent in left
        for right_coefficient, right_exponent in right
    ]


def tangent_harmonic_polynomials(max_order: int) -> list[Polynomial]:
    """Return the O(2)-zonal kernels around one spherical root.

    If a=X.Y, b=Y.Z, c=Z.X, then R_k is

      ((1-a^2)(1-c^2))^(k/2) cos(k(phi_Y-phi_Z)).

    The square roots cancel through the Chebyshev recurrence, leaving a
    polynomial in the three Gram variables.  Positive semidefinite blocks
    against R_k are precisely one-root conditional flag squares in the
    weight-k representation of the root stabilizer O(2).
    """

    if max_order < 0:
        return []
    one = [(Fraction(1), (0, 0, 0))]
    if max_order == 0:
        return [one]

    tangent_dot = [
        (Fraction(1), (0, 1, 0)),
        (Fraction(-1), (1, 0, 1)),
    ]
    tangent_norm_product = [
        (Fraction(1), (0, 0, 0)),
        (Fraction(-1), (2, 0, 0)),
        (Fraction(-1), (0, 0, 2)),
        (Fraction(1), (2, 0, 2)),
    ]
    result = [one, tangent_dot]
    for order in range(1, max_order):
        next_polynomial = add_polynomials(
            scale_polynomial(
                multiply_polynomials(tangent_dot, result[order]),
                Fraction(2),
            ),
            scale_polynomial(
                multiply_polynomials(tangent_norm_product, result[order - 1]),
                Fraction(-1),
            ),
        )
        result.append(next_polynomial)
    return result


def flag_expectation_matrix(
    leaf_degrees: list[int],
    tangent_harmonic: Polynomial,
) -> dict[Label, np.ndarray]:
    """Expectation map for a one-root conditional flag-square block."""

    size = len(leaf_degrees)
    matrices: dict[Label, np.ndarray] = {}
    for row, left_degree in enumerate(leaf_degrees):
        for column, right_degree in enumerate(leaf_degrees):
            for coefficient, exponent in tangent_harmonic:
                shifted = (
                    exponent[0] + left_degree,
                    exponent[1],
                    exponent[2] + right_degree,
                )
                label, reduction_coefficient = expectation_label(shifted)
                if label is None or reduction_coefficient == 0:
                    continue
                matrix = matrices.setdefault(label, np.zeros((size, size)))
                matrix[row, column] += float(coefficient * reduction_coefficient)
    return matrices


ONE: Polynomial = [(Fraction(1), (0, 0, 0))]
GRAM_DETERMINANT: Polynomial = [
    (Fraction(1), (0, 0, 0)),
    (Fraction(-1), (2, 0, 0)),
    (Fraction(-1), (0, 2, 0)),
    (Fraction(-1), (0, 0, 2)),
    (Fraction(2), (1, 1, 1)),
]
PRINCIPAL_MINORS: list[Polynomial] = [
    [(Fraction(1), (0, 0, 0)), (Fraction(-1), (2, 0, 0))],
    [(Fraction(1), (0, 0, 0)), (Fraction(-1), (0, 2, 0))],
    [(Fraction(1), (0, 0, 0)), (Fraction(-1), (0, 0, 2))],
]


def make_psd_block(
    name: str,
    basis: list[Exponent],
    multiplier: Polynomial,
) -> tuple[cp.Variable, dict[Label, np.ndarray]]:
    variable = cp.Variable((len(basis), len(basis)), symmetric=True, name=name)
    return variable, expectation_matrix(basis, multiplier)


def fraction_decimal(value: Fraction, digits: int) -> str:
    """Decimal string of a rational with the requested significant digits."""

    import decimal

    with decimal.localcontext() as context:
        context.prec = digits
        quotient = (
            decimal.Decimal(value.numerator)
            / decimal.Decimal(value.denominator)
        )
        return format(quotient, "E")


def export_sdpa_problem(
    path: Path,
    digits: int,
    target: dict[Label, float],
    ordered_labels: list[Label],
    psd_blocks: list[tuple[str, dict[Label, np.ndarray]]],
    free_label_matrices: list[dict[Label, np.ndarray]],
    relations: list[dict[Label, float]],
) -> dict[str, object]:
    """Write the dual moment problem in exact SDPA sparse format.

    Every equality constraint (KKT relations, rank relations, free matrix
    blocks) is eliminated over the rationals: the moment vector is
    parameterized as y = y0 + sum_j z_j p_j with y0, p_j exact, y0 the
    normalized solution of y_constant = 1, and p_j spanning the equality
    kernel with vanishing constant coordinate.  The exported problem

        minimize  sum_j c_j z_j  subject to  M0 + sum_j z_j M_j >= 0

    is a pure semidefinite program in free variables, so an interior
    exists whenever the original relaxation has one.  All data are written
    as exact decimals so an extended-precision solver is not contaminated
    by binary64 rounding.  The reported bound is objValue + objective_shift.
    """

    constant = ("constant",)
    quotient_rows, generator_count, equality_rank = (
        exact_equality_quotient_rows(
            ordered_labels,
            free_label_matrices,
            relations,
        )
    )

    pivot_index = next(
        (
            index
            for index, row in enumerate(quotient_rows)
            if row.get(constant)
        ),
        None,
    )
    if pivot_index is None:
        raise ValueError("Equality kernel forces the constant moment to 0")
    pivot_row = quotient_rows[pivot_index]
    pivot_value = pivot_row[constant]
    base_point = {
        label: coefficient / pivot_value
        for label, coefficient in pivot_row.items()
    }
    directions: list[dict[Label, Fraction]] = []
    for index, row in enumerate(quotient_rows):
        if index == pivot_index:
            continue
        weight = row.get(constant, Fraction(0)) / pivot_value
        direction = dict(row)
        if weight:
            for label, coefficient in pivot_row.items():
                updated = (
                    direction.get(label, Fraction(0)) - weight * coefficient
                )
                if updated:
                    direction[label] = updated
                else:
                    direction.pop(label, None)
        if direction:
            directions.append(direction)

    rational_target = {
        label: rationalize_float(float(coefficient))
        for label, coefficient in target.items()
    }

    def pair_with_target(vector: dict[Label, Fraction]) -> Fraction:
        return sum(
            (
                coefficient * rational_target[label]
                for label, coefficient in vector.items()
                if label in rational_target
            ),
            Fraction(0),
        )

    objective_shift = pair_with_target(base_point)
    objective = [pair_with_target(direction) for direction in directions]

    # Sparse exact block data: for each block, label -> upper-tri entries.
    sparse_blocks: list[dict[Label, list[tuple[int, int, Fraction]]]] = []
    block_sizes: list[int] = []
    for _, label_matrices in psd_blocks:
        size = next(iter(label_matrices.values())).shape[0]
        block_sizes.append(size)
        sparse: dict[Label, list[tuple[int, int, Fraction]]] = {}
        for label, matrix in label_matrices.items():
            entries = []
            for row in range(size):
                for column in range(row, size):
                    value = 0.5 * (
                        matrix[row, column] + matrix[column, row]
                    )
                    if abs(value) > 1e-13:
                        entries.append(
                            (row, column, rationalize_float(float(value)))
                        )
            if entries:
                sparse[label] = entries
        sparse_blocks.append(sparse)

    def assemble(
        vector: dict[Label, Fraction],
    ) -> list[dict[tuple[int, int], Fraction]]:
        matrices: list[dict[tuple[int, int], Fraction]] = []
        for sparse in sparse_blocks:
            matrix: dict[tuple[int, int], Fraction] = {}
            for label, coefficient in vector.items():
                for row, column, value in sparse.get(label, ()):
                    updated = (
                        matrix.get((row, column), Fraction(0))
                        + coefficient * value
                    )
                    if updated:
                        matrix[(row, column)] = updated
                    else:
                        matrix.pop((row, column), None)
            matrices.append(matrix)
        return matrices

    # SDPA requires linearly independent constraint matrices.  Directions
    # supported on labels absent from every PSD block have zero image, and
    # further image-space dependencies are possible; select an image basis
    # by pivoted Cholesky on the Gram matrix of the flattened images.
    direction_images = [assemble(direction) for direction in directions]

    def flatten(
        matrices: list[dict[tuple[int, int], Fraction]],
    ) -> dict[tuple[int, int, int], Fraction]:
        flat: dict[tuple[int, int, int], Fraction] = {}
        for block_index, matrix in enumerate(matrices):
            for (row, column), value in matrix.items():
                flat[(block_index, row, column)] = value
        return flat

    flattened = [flatten(matrices) for matrices in direction_images]
    count = len(flattened)

    # Exact sparse Gaussian elimination over Q selects a maximal subset of
    # directions with linearly independent images.  Float rank decisions
    # are unusable here: dropping a nearly-parallel but independent
    # direction perturbs the feasible set, which the extended-precision
    # solver then correctly reports as infeasible.
    echelon: dict[tuple[int, int, int], dict[tuple[int, int, int], Fraction]]
    echelon = {}
    selected: list[int] = []
    for index in range(count):
        vector = dict(flattened[index])
        while vector:
            pivot_coordinate = min(
                vector,
                key=lambda coordinate: (
                    coordinate not in echelon,
                    coordinate,
                ),
            )
            row = echelon.get(pivot_coordinate)
            if row is None:
                pivot_value = vector[pivot_coordinate]
                echelon[pivot_coordinate] = {
                    coordinate: value / pivot_value
                    for coordinate, value in vector.items()
                }
                selected.append(index)
                break
            factor = vector[pivot_coordinate]
            for coordinate, value in row.items():
                updated = (
                    vector.get(coordinate, Fraction(0)) - factor * value
                )
                if updated:
                    vector[coordinate] = updated
                else:
                    vector.pop(coordinate, None)
    selected.sort()
    dropped = count - len(selected)

    # Normalize each kept direction exactly so the exported variables are
    # O(1): divide by the largest absolute image coefficient.
    normalizers: dict[int, Fraction] = {}
    for index in selected:
        largest = max(
            (
                abs(value)
                for matrix in direction_images[index]
                for value in matrix.values()
            ),
            default=Fraction(0),
        )
        normalizer = largest if largest else Fraction(1)
        normalizers[index] = normalizer
        direction_images[index] = [
            {
                coordinates: value / normalizer
                for coordinates, value in matrix.items()
            }
            for matrix in direction_images[index]
        ]
        objective[index] = objective[index] / normalizer

    variable_count = len(selected)
    lines: list[str] = [
        f"{variable_count} = mDIM",
        f"{len(block_sizes)} = nBLOCK",
        "(" + ", ".join(str(size) for size in block_sizes) + ") = bLOCKsTRUCT",
    ]
    lines.append(
        "{"
        + ", ".join(
            fraction_decimal(objective[index], digits)
            for index in selected
        )
        + "}"
    )
    entry_count = 0
    # F_0 = -M0 so that sum_j z_j M_j - F_0 = M0 + sum_j z_j M_j.
    for block_index, matrix in enumerate(assemble(base_point)):
        for (row, column), value in sorted(matrix.items()):
            lines.append(
                f"0 {block_index + 1} {row + 1} {column + 1} "
                f"{fraction_decimal(-value, digits)}"
            )
            entry_count += 1
    for variable_index, direction_index in enumerate(selected):
        for block_index, matrix in enumerate(
            direction_images[direction_index]
        ):
            for (row, column), value in sorted(matrix.items()):
                lines.append(
                    f"{variable_index + 1} {block_index + 1} "
                    f"{row + 1} {column + 1} "
                    f"{fraction_decimal(value, digits)}"
                )
                entry_count += 1
    path.write_text("\n".join(lines) + "\n")

    # Sidecar map for solution extraction: reconstruct the moment vector
    # y = y0 + sum_j z_j q_j / N_j from the solver's x, and identify the
    # certificate blocks (SDPA's dual yMat) with the named flag blocks.
    map_path = path.with_name(path.name + ".map.json")
    map_path.write_text(
        json.dumps(
            {
                "blocks": [
                    {"name": name, "size": size}
                    for (name, _), size in zip(
                        psd_blocks,
                        block_sizes,
                        strict=True,
                    )
                ],
                "objective_shift": str(objective_shift),
                "objective": [
                    str(objective[index]) for index in selected
                ],
                "base_point": {
                    str(label): str(value)
                    for label, value in base_point.items()
                },
                "directions": [
                    {
                        "normalizer": str(normalizers[index]),
                        "coefficients": {
                            str(label): str(value)
                            for label, value in directions[index].items()
                        },
                    }
                    for index in selected
                ],
            }
        )
    )

    return {
        "export": str(path),
        "map": str(map_path),
        "variables": variable_count,
        "dropped_dependent_directions": dropped,
        "labels": len(ordered_labels),
        "equality_generators": generator_count,
        "equality_rank": equality_rank,
        "block_sizes": {
            name: size
            for (name, _), size in zip(psd_blocks, block_sizes, strict=True)
        },
        "entries": entry_count,
        "objective_shift": float(objective_shift),
        "objective_shift_exact": str(objective_shift),
        "note": "bound = objValPrimal + objective_shift",
    }


def solve(args: argparse.Namespace) -> dict[str, object]:
    # T = (3/16) E = -1/4 + (15/4) p2 - 9 p4 + 6 p6, with no isotropy
    # substitution: p2 = E[(X.Y)^2] stays a genuine moment variable.
    target = {
        ("constant",): -0.25 + args.target_epsilon,
        ("pair", 2): 3.75,
        ("pair", 4): -9.0,
        ("pair", 6): 6.0,
    }
    gradient, hessian, perpendicular_hessian = kernel_polynomials()
    (
        four_point_parallel_hessian,
        four_point_perpendicular_hessian,
    ) = four_point_hessian_polynomials()
    (
        global_parallel_tangent_gap,
        global_perpendicular_tangent_gap,
    ) = global_tangent_gap_polynomials()

    blocks: list[tuple[str, cp.Variable, dict[Label, np.ndarray]]] = []
    free_blocks: list[tuple[str, cp.Variable, dict[Label, np.ndarray]]] = []
    constraints: list[cp.Constraint] = []

    module_terms: list[tuple[str, Polynomial, int]] = []
    if not args.no_pointwise_sos:
        module_terms.append(("sos", ONE, 0))
    if args.gram_module:
        module_terms.append(("det", GRAM_DETERMINANT, 3))
        module_terms.extend(
            (f"minor_{index}", polynomial, 2)
            for index, polynomial in enumerate(PRINCIPAL_MINORS)
        )

    for name, multiplier, multiplier_degree in module_terms:
        basis_degree = (args.degree - multiplier_degree) // 2
        basis = monomials(basis_degree)
        variable, label_matrices = make_psd_block(name, basis, multiplier)
        blocks.append((name, variable, label_matrices))
        constraints.append(variable >> 0)

    if args.three_point_flags:
        harmonics = tangent_harmonic_polynomials(args.degree // 2)
        for order, tangent_harmonic in enumerate(harmonics):
            # Every entry has degree at most degree.  Leaf parity must match
            # the O(2) weight for antipodally invariant (projective) flags.
            maximum_leaf_degree = (args.degree - 2 * order) // 2
            leaf_degrees = list(
                range(order % 2, maximum_leaf_degree + 1, 2)
            )
            if not leaf_degrees:
                continue
            variable = cp.Variable(
                (len(leaf_degrees), len(leaf_degrees)),
                symmetric=True,
                name=f"flag_{order}",
            )
            blocks.append(
                (
                    f"flag_{order}",
                    variable,
                    flag_expectation_matrix(leaf_degrees, tangent_harmonic),
                )
            )
            constraints.append(variable >> 0)

    if args.h2_localized_flags:
        localized_harmonics = tangent_harmonic_polynomials(
            max(0, (args.degree - 2) // 2)
        )
        for order, tangent_harmonic in enumerate(localized_harmonics):
            maximum_leaf_degree = (args.degree - 2 - 2 * order) // 2
            leaf_degrees = list(
                range(order % 2, maximum_leaf_degree + 1, 2)
            )
            if not leaf_degrees:
                continue
            label_matrices = h2_localized_flag_expectation_matrix(
                leaf_degrees,
                tangent_harmonic,
            )
            if not label_matrices:
                continue
            variable = cp.Variable(
                (len(leaf_degrees), len(leaf_degrees)),
                symmetric=True,
                name=f"h2_flag_{order}",
            )
            blocks.append((f"h2_flag_{order}", variable, label_matrices))
            constraints.append(variable >> 0)

    if args.four_point_flags:
        maximum_pair_degree = args.degree // 2
        pair_degrees = list(range(0, maximum_pair_degree + 1, 2))
        variable = cp.Variable(
            (len(pair_degrees), len(pair_degrees)),
            symmetric=True,
            name="empty_type_flag",
        )
        blocks.append(
            (
                "empty_type_flag",
                variable,
                empty_type_flag_expectation_matrix(pair_degrees),
            )
        )
        constraints.append(variable >> 0)

    if args.two_root_flags:
        flag_degree = args.degree // 2
        all_flags = monomials(flag_degree)
        sectors = [
            (
                "two_root_even_00",
                [
                    exponent
                    for exponent in all_flags
                    if (exponent[0] + exponent[1]) % 2 == 0
                    and (exponent[0] + exponent[2]) % 2 == 0
                    and (exponent[1] + exponent[2]) % 2 == 0
                ],
                None,
            ),
            (
                "two_root_even_11",
                [
                    exponent
                    for exponent in all_flags
                    if (exponent[0] + exponent[1]) % 2 == 0
                    and (exponent[0] + exponent[2]) % 2 == 1
                    and (exponent[1] + exponent[2]) % 2 == 1
                ],
                None,
            ),
            (
                "two_root_odd_01",
                [
                    exponent
                    for exponent in all_flags
                    if (exponent[0] + exponent[1]) % 2 == 1
                    and (exponent[0] + exponent[2]) % 2 == 0
                    and (exponent[1] + exponent[2]) % 2 == 1
                ],
                ORIENTATION_PAIRING,
            ),
            (
                "two_root_odd_10",
                [
                    exponent
                    for exponent in all_flags
                    if (exponent[0] + exponent[1]) % 2 == 1
                    and (exponent[0] + exponent[2]) % 2 == 1
                    and (exponent[1] + exponent[2]) % 2 == 0
                ],
                ORIENTATION_PAIRING,
            ),
        ]
        for name, flag_basis, multiplier in sectors:
            if not flag_basis:
                continue
            variable = cp.Variable(
                (len(flag_basis), len(flag_basis)),
                symmetric=True,
                name=name,
            )
            blocks.append(
                (
                    name,
                    variable,
                    two_root_flag_expectation_matrix(
                        flag_basis,
                        multiplier,
                    ),
                )
            )
            constraints.append(variable >> 0)

        localizing_degree = (args.degree - 2) // 2
        localizing_flags = monomials(localizing_degree)
        localizing_sectors = [
            (
                "two_root_even_00_minor",
                [
                    exponent
                    for exponent in localizing_flags
                    if (exponent[0] + exponent[1]) % 2 == 0
                    and (exponent[0] + exponent[2]) % 2 == 0
                    and (exponent[1] + exponent[2]) % 2 == 0
                ],
                ROOT_PAIR_MINOR,
            ),
            (
                "two_root_even_11_minor",
                [
                    exponent
                    for exponent in localizing_flags
                    if (exponent[0] + exponent[1]) % 2 == 0
                    and (exponent[0] + exponent[2]) % 2 == 1
                    and (exponent[1] + exponent[2]) % 2 == 1
                ],
                ROOT_PAIR_MINOR,
            ),
            (
                "two_root_odd_01_minor",
                [
                    exponent
                    for exponent in localizing_flags
                    if (exponent[0] + exponent[1]) % 2 == 1
                    and (exponent[0] + exponent[2]) % 2 == 0
                    and (exponent[1] + exponent[2]) % 2 == 1
                ],
                multiply_graph_polynomials(
                    ORIENTATION_PAIRING,
                    ROOT_PAIR_MINOR,
                ),
            ),
            (
                "two_root_odd_10_minor",
                [
                    exponent
                    for exponent in localizing_flags
                    if (exponent[0] + exponent[1]) % 2 == 1
                    and (exponent[0] + exponent[2]) % 2 == 1
                    and (exponent[1] + exponent[2]) % 2 == 0
                ],
                multiply_graph_polynomials(
                    ORIENTATION_PAIRING,
                    ROOT_PAIR_MINOR,
                ),
            ),
        ]
        for name, flag_basis, multiplier in localizing_sectors:
            if not flag_basis:
                continue
            variable = cp.Variable(
                (len(flag_basis), len(flag_basis)),
                symmetric=True,
                name=name,
            )
            blocks.append(
                (
                    name,
                    variable,
                    two_root_flag_expectation_matrix(
                        flag_basis,
                        multiplier,
                    ),
                )
            )
            constraints.append(variable >> 0)

    if args.max_flag_arity >= 5:
        maximum_leaf_degree = args.degree // 2
        for arity in range(5, args.max_flag_arity + 1):
            root_count = arity - 2
            if args.max_root_factor_degree > 0:
                sectors = rooted_weighted_flag_sectors(
                    root_count,
                    maximum_leaf_degree,
                    args.max_root_factor_degree,
                )
            else:
                sectors = rooted_star_flag_sectors(
                    root_count,
                    maximum_leaf_degree,
                    leaf_parity=0,
                )
            for signature, flag_basis in sorted(sectors.items()):
                if not flag_basis:
                    continue
                signature_name = "".join(str(value) for value in signature)
                if args.max_root_factor_degree > 0:
                    name = f"weighted_flag_{arity}_{signature_name}"
                    label_matrices = rooted_weighted_flag_expectation_matrix(
                        root_count,
                        flag_basis,
                    )
                else:
                    name = f"star_flag_{arity}_{signature_name}"
                    label_matrices = rooted_star_flag_expectation_matrix(
                        root_count,
                        flag_basis,
                    )
                if not label_matrices:
                    continue
                variable = cp.Variable(
                    (len(flag_basis), len(flag_basis)),
                    symmetric=True,
                    name=name,
                )
                blocks.append(
                    (
                        name,
                        variable,
                        label_matrices,
                    )
                )
                constraints.append(variable >> 0)

    if args.higher_rank_matrices and args.max_flag_arity >= 5:
        determinant = four_vertex_gram_determinant()
        maximum_leaf_degree = (args.degree - 4) // 2
        for arity in range(5, args.max_flag_arity + 1):
            root_count = arity - 2
            sectors = rooted_star_flag_sectors(
                root_count,
                maximum_leaf_degree,
                leaf_parity=0,
            )
            for vertices in itertools.combinations(range(arity), 4):
                lifted_determinant = lift_graph_polynomial(
                    determinant,
                    4,
                    arity,
                    vertices,
                )
                subset_name = "".join(str(vertex) for vertex in vertices)
                for signature, flag_basis in sorted(sectors.items()):
                    if not flag_basis:
                        continue
                    label_matrices = rooted_star_flag_expectation_matrix(
                        root_count,
                        flag_basis,
                        lifted_determinant,
                    )
                    if not label_matrices:
                        continue
                    signature_name = "".join(
                        str(value) for value in signature
                    )
                    name = (
                        f"rank_flag_{arity}_{subset_name}_"
                        f"{signature_name}"
                    )
                    variable = cp.Variable(
                        (len(flag_basis), len(flag_basis)),
                        symmetric=True,
                        name=name,
                    )
                    free_blocks.append(
                        (
                            name,
                            variable,
                            label_matrices,
                        )
                    )

    if args.harmonics:
        # Degree 2 encodes the isotropy deficit as a flag square:
        # E[P_2(X.Y)] = (3 p2 - 1)/2 = sum_m |mu-hat(2,m)|^2 >= 0.
        for degree in range(2, args.degree + 1, 2):
            variable = cp.Variable(
                (1, 1),
                symmetric=True,
                name=f"harmonic_{degree}",
            )
            blocks.append(
                (
                    f"harmonic_{degree}",
                    variable,
                    harmonic_pair_vector(degree),
                )
            )
            constraints.append(variable >> 0)

        # Rich spin-2 Gram block: harmonic vertex with up to two leaves.
        # Contains D itself, so h_2 -> 0 forces the contraction identities.
        spin2_degree = (args.degree - 2) // 2
        spin2_basis = spin2_flag_basis(spin2_degree)
        spin2_matrices = spin2_flag_expectation_matrix(spin2_basis)
        if spin2_matrices:
            variable = cp.Variable(
                (len(spin2_basis), len(spin2_basis)),
                symmetric=True,
                name="spin2_flag",
            )
            blocks.append(("spin2_flag", variable, spin2_matrices))
            constraints.append(variable >> 0)

        # Spin-l Gram blocks of harmonic-weighted unrooted pair flags.
        # The l = 2 block replaces the dropped isotropy contraction by
        # genuine flag squares of the deviatoric second-moment tensor.
        for order in range(2, args.degree + 1, 2):
            maximum_weight_degree = (args.degree - order) // 2
            weight_degrees = list(range(0, maximum_weight_degree + 1, 2))
            if len(weight_degrees) < 2:
                continue
            label_matrices = harmonic_flag_expectation_matrix(
                order,
                weight_degrees,
            )
            if not label_matrices:
                continue
            variable = cp.Variable(
                (len(weight_degrees), len(weight_degrees)),
                symmetric=True,
                name=f"harmonic_flag_{order}",
            )
            blocks.append(
                (
                    f"harmonic_flag_{order}",
                    variable,
                    label_matrices,
                )
            )
            constraints.append(variable >> 0)

    if args.global_gap:
        variable = cp.Variable((1, 1), symmetric=True, name="global_uniform_gap")
        blocks.append(
            (
                "global_uniform_gap",
                variable,
                {
                    # int U dsigma - E = 172/105 - 20 p2 + 48 p4 - 32 p6
                    ("constant",): np.array([[172.0 / 105.0]]),
                    ("pair", 2): np.array([[-20.0]]),
                    ("pair", 4): np.array([[48.0]]),
                    ("pair", 6): np.array([[-32.0]]),
                },
            )
        )
        constraints.append(variable >> 0)

    if args.global_tangent_gaps:
        gap_basis_degree = (args.degree - 12) // 2
        if gap_basis_degree >= 0:
            basis = [(0, 0, power) for power in range(gap_basis_degree + 1)]
            for name, polynomial in (
                ("global_parallel_tangent_gap", global_parallel_tangent_gap),
                (
                    "global_perpendicular_tangent_gap",
                    global_perpendicular_tangent_gap,
                ),
            ):
                variable, label_matrices = make_psd_block(
                    name,
                    basis,
                    polynomial,
                )
                blocks.append((name, variable, label_matrices))
                constraints.append(variable >> 0)

        gap_minor_degree = gap_basis_degree - 1
        if gap_minor_degree >= 0:
            basis = [(0, 0, power) for power in range(gap_minor_degree + 1)]
            for name, polynomial in (
                ("global_parallel_tangent_gap_minor", global_parallel_tangent_gap),
                (
                    "global_perpendicular_tangent_gap_minor",
                    global_perpendicular_tangent_gap,
                ),
            ):
                variable, label_matrices = make_psd_block(
                    name,
                    basis,
                    multiply_polynomials(polynomial, PRINCIPAL_MINORS[2]),
                )
                blocks.append((name, variable, label_matrices))
                constraints.append(variable >> 0)

    if args.potential_matrices:
        harmonics = tangent_harmonic_polynomials(max(0, (args.degree - 6) // 2))
        for order, tangent_harmonic in enumerate(harmonics):
            maximum_leaf_degree = (args.degree - 6 - 2 * order) // 2
            leaf_degrees = list(
                range(order % 2, maximum_leaf_degree + 1, 2)
            )
            if not leaf_degrees:
                continue
            variable = cp.Variable(
                (len(leaf_degrees), len(leaf_degrees)),
                symmetric=True,
                name=f"potential_flag_{order}",
            )
            free_blocks.append(
                (
                    f"potential_flag_{order}",
                    variable,
                    potential_flag_relation_matrix(
                        leaf_degrees,
                        tangent_harmonic,
                    ),
                )
            )

    if args.hessian:
        # A nonnegative polynomial r(c) is represented as
        # s(c)^T H_0 s(c) + (1-c^2)t(c)^T H_1 t(c).
        hessian_basis_degree = (args.degree - 8) // 2
        if hessian_basis_degree >= 0:
            basis = [(0, 0, power) for power in range(hessian_basis_degree + 1)]
            variable, label_matrices = make_psd_block(
                "hessian_sos",
                basis,
                hessian,
            )
            blocks.append(("hessian_sos", variable, label_matrices))
            constraints.append(variable >> 0)

            perpendicular_variable, perpendicular_label_matrices = make_psd_block(
                "perpendicular_hessian_sos",
                basis,
                perpendicular_hessian,
            )
            blocks.append(
                (
                    "perpendicular_hessian_sos",
                    perpendicular_variable,
                    perpendicular_label_matrices,
                )
            )
            constraints.append(perpendicular_variable >> 0)

        hessian_minor_degree = hessian_basis_degree - 1
        if hessian_minor_degree >= 0:
            basis = [(0, 0, power) for power in range(hessian_minor_degree + 1)]
            hessian_minor = multiply_polynomials(hessian, PRINCIPAL_MINORS[2])
            variable, label_matrices = make_psd_block(
                "hessian_minor",
                basis,
                hessian_minor,
            )
            blocks.append(("hessian_minor", variable, label_matrices))
            constraints.append(variable >> 0)

            perpendicular_hessian_minor = multiply_polynomials(
                perpendicular_hessian,
                PRINCIPAL_MINORS[2],
            )
            (
                perpendicular_variable,
                perpendicular_label_matrices,
            ) = make_psd_block(
                "perpendicular_hessian_minor",
                basis,
                perpendicular_hessian_minor,
            )
            blocks.append(
                (
                    "perpendicular_hessian_minor",
                    perpendicular_variable,
                    perpendicular_label_matrices,
                )
            )
            constraints.append(perpendicular_variable >> 0)

    if args.four_point_hessian:
        maximum_auxiliary_degree = (args.degree - 8) // 2
        auxiliary_degrees = list(
            range(1, maximum_auxiliary_degree + 1, 2)
        )
        if auxiliary_degrees:
            for name, polynomial in (
                ("four_point_parallel_hessian", four_point_parallel_hessian),
                (
                    "four_point_perpendicular_hessian",
                    four_point_perpendicular_hessian,
                ),
            ):
                variable = cp.Variable(
                    (len(auxiliary_degrees), len(auxiliary_degrees)),
                    symmetric=True,
                    name=name,
                )
                blocks.append(
                    (
                        name,
                        variable,
                        four_point_hessian_expectation_matrix(
                            auxiliary_degrees,
                            polynomial,
                        ),
                    )
                )
                constraints.append(variable >> 0)

    if args.max_hessian_arity >= 5:
        maximum_vector_flag_degree = (args.degree - 8) // 2
        if maximum_vector_flag_degree >= 1:
            for arity in range(5, args.max_hessian_arity + 1):
                conditioning_root_count = arity - 3
                sectors = rooted_star_flag_sectors(
                    conditioning_root_count,
                    maximum_vector_flag_degree,
                    leaf_parity=1,
                )
                for signature, flag_basis in sorted(sectors.items()):
                    if not flag_basis:
                        continue
                    signature_name = "".join(
                        str(value) for value in signature
                    )
                    for direction, polynomial in (
                        ("parallel", four_point_parallel_hessian),
                        ("perpendicular", four_point_perpendicular_hessian),
                    ):
                        name = (
                            f"hessian_flag_{arity}_{direction}_"
                            f"{signature_name}"
                        )
                        label_matrices = lifted_hessian_expectation_matrix(
                            conditioning_root_count,
                            flag_basis,
                            polynomial,
                        )
                        if not label_matrices:
                            continue
                        variable = cp.Variable(
                            (len(flag_basis), len(flag_basis)),
                            symmetric=True,
                            name=name,
                        )
                        blocks.append(
                            (
                                name,
                                variable,
                                label_matrices,
                            )
                        )
                        constraints.append(variable >> 0)

    if args.jacobi_scale_blocks:
        # Exact congruence rescaling A_L -> D A_L D with diagonal D chosen
        # from the largest diagonal coefficient per basis element.  This
        # preserves positive semidefiniteness and every certificate
        # identity (the PSD variable absorbs D^{-1} on both sides) while
        # compressing the coefficient range MOSEK must handle.
        def jacobi_rescale(
            label_matrices: dict[Label, np.ndarray],
        ) -> dict[Label, np.ndarray]:
            size = next(iter(label_matrices.values())).shape[0]
            diagonal_scale = np.ones(size)
            for index in range(size):
                largest = max(
                    float(
                        max(
                            np.max(np.abs(matrix[index, :])),
                            np.max(np.abs(matrix[:, index])),
                        )
                    )
                    for matrix in label_matrices.values()
                )
                if largest > 0:
                    diagonal_scale[index] = 1.0 / math.sqrt(largest)
            return {
                label: diagonal_scale[:, None]
                * matrix
                * diagonal_scale[None, :]
                for label, matrix in label_matrices.items()
            }

        blocks = [
            (name, variable, jacobi_rescale(label_matrices))
            for name, variable, label_matrices in blocks
            if label_matrices
        ]
        free_blocks = [
            (name, variable, jacobi_rescale(label_matrices))
            for name, variable, label_matrices in free_blocks
            if label_matrices
        ]

    gradient_relations: list[tuple[int, dict[Label, float]]] = []
    if args.gradient:
        for power in range(max(0, args.degree - 7) + 1):
            relation = expectation_vector(
                shifted_polynomial(gradient, (0, 0, power))
            )
            if relation:
                gradient_relations.append((power, relation))
    gradient_coefficients = (
        cp.Variable(len(gradient_relations), name="gradient_coefficients")
        if gradient_relations
        else None
    )

    potential_relations: list[tuple[int, dict[Label, float]]] = []
    if args.potential:
        for power in range(0, max(0, args.degree - 6) + 1, 2):
            relation = potential_stationarity_relation(power)
            if relation:
                potential_relations.append((power, relation))
    potential_coefficients = (
        cp.Variable(len(potential_relations), name="potential_coefficients")
        if potential_relations
        else None
    )

    rank_relations: list[tuple[int, dict[Label, float]]] = []
    if args.rank_relations:
        rank_relations = list(
            enumerate(
                four_point_rank_relations(max(0, args.degree - 4))
            )
        )
    rank_coefficients = (
        cp.Variable(len(rank_relations), name="rank_coefficients")
        if rank_relations
        else None
    )

    labels = set(target)
    for _, _, label_matrices in blocks:
        labels.update(label_matrices)
    for _, _, label_matrices in free_blocks:
        labels.update(label_matrices)
    for _, relation in gradient_relations:
        labels.update(relation)
    for _, relation in potential_relations:
        labels.update(relation)
    for _, relation in rank_relations:
        labels.update(relation)

    ordered_labels = sorted(labels, key=str)

    if args.export_sdpa:
        if args.jacobi_scale_blocks:
            raise ValueError(
                "--export-sdpa needs unscaled rational blocks"
            )
        return export_sdpa_problem(
            Path(args.export_sdpa),
            args.sdpa_digits,
            target,
            ordered_labels,
            [
                (name, label_matrices)
                for name, _, label_matrices in blocks
                if label_matrices
            ],
            [label_matrices for _, _, label_matrices in free_blocks],
            [
                relation
                for _, relation in (
                    gradient_relations
                    + potential_relations
                    + rank_relations
                )
            ],
        )

    if args.facial_reduce_onb:
        if args.target_epsilon != 0:
            raise ValueError(
                "ONB facial reduction is valid only for the sharp target"
            )
        reduced_blocks: list[
            tuple[str, cp.Variable, dict[Label, np.ndarray]]
        ] = []
        constraints = []
        for name, _, label_matrices in blocks:
            if args.exact_onb_face:
                exact_onb_kernel = exact_onb_nullspace(label_matrices)
                nullspace = np.array(exact_onb_kernel, dtype=float)
            else:
                onb_matrix = sum(
                    float(onb_label_value(label)) * coefficient_matrix
                    for label, coefficient_matrix in label_matrices.items()
                )
                eigenvalues, eigenvectors = np.linalg.eigh(onb_matrix)
                if eigenvalues[0] < -1e-8:
                    raise ValueError(
                        f"Block {name} is not PSD on the ONB: "
                        f"{eigenvalues[0]}"
                    )
                nullspace = eigenvectors[:, eigenvalues <= 1e-8]
            if nullspace.shape[1] == 0:
                continue
            reduced_matrices = {
                label: nullspace.T @ matrix @ nullspace
                for label, matrix in label_matrices.items()
            }
            reduced_matrices = {
                label: matrix
                for label, matrix in reduced_matrices.items()
                if np.max(np.abs(matrix)) > 1e-13
            }
            if not reduced_matrices:
                continue
            reduced_variable = cp.Variable(
                (nullspace.shape[1], nullspace.shape[1]),
                symmetric=True,
                name=f"{name}_onb_face",
            )
            reduced_blocks.append(
                (
                    name,
                    reduced_variable,
                    reduced_matrices,
                )
            )
            constraints.append(reduced_variable >> 0)
        blocks = reduced_blocks

    pole_equator_faces: list[int] = []
    if args.pole_equator_faces:
        for token in args.pole_equator_faces.split(","):
            token = token.strip().lower()
            pole_equator_faces.append(
                0 if token in {"continuous", "haar"} else int(token)
            )
    for regular_order in pole_equator_faces:
        if args.target_epsilon != 0:
            raise ValueError(
                "Pole-equator facial reduction is valid only for the "
                "sharp target"
            )
        target_value = sum(
            Fraction(str(coefficient))
            * pole_equator_label_value(label, regular_order)
            for label, coefficient in target.items()
        )
        if target_value:
            raise ValueError(
                "Requested pole-equator measure does not annihilate target"
            )

        reduced_blocks = []
        constraints = []
        face_name = (
            "continuous"
            if regular_order == 0
            else f"regular_{regular_order}"
        )
        for name, _, label_matrices in blocks:
            moment_matrix = exact_moment_matrix(
                label_matrices,
                lambda label, order=regular_order: (
                    pole_equator_label_value(label, order)
                ),
            )
            numerical_moment_matrix = np.array(moment_matrix, dtype=float)
            if (
                numerical_moment_matrix.size
                and np.linalg.eigvalsh(numerical_moment_matrix)[0] < -1e-8
            ):
                raise ValueError(
                    f"Block {name} is not PSD on pole-equator "
                    f"face {face_name}"
                )
            exact_kernel = exact_nullspace(moment_matrix)
            if exact_kernel.shape[1] == 0:
                continue
            kernel = np.array(exact_kernel, dtype=float)
            reduced_matrices = {
                label: kernel.T @ matrix @ kernel
                for label, matrix in label_matrices.items()
            }
            reduced_matrices = {
                label: matrix
                for label, matrix in reduced_matrices.items()
                if np.max(np.abs(matrix)) > 1e-13
            }
            if not reduced_matrices:
                continue
            reduced_variable = cp.Variable(
                (kernel.shape[1], kernel.shape[1]),
                symmetric=True,
                name=f"{name}_{face_name}_face",
            )
            reduced_blocks.append(
                (
                    name,
                    reduced_variable,
                    reduced_matrices,
                )
            )
            constraints.append(reduced_variable >> 0)
        blocks = reduced_blocks

    numerical_face_paths = [
        Path(token.strip())
        for token in args.numerical_faces.split(",")
        if token.strip()
    ] if args.numerical_faces else []
    for face_index, face_path in enumerate(numerical_face_paths):
        data = np.load(face_path)
        stored_labels = list(data["labels"])
        expected_labels = [repr(label) for label in ordered_labels]
        if stored_labels != expected_labels:
            raise ValueError(
                f"Numerical face {face_path} uses a different label basis"
            )
        face_moments = np.array(data["moments"], dtype=float)
        reduced_blocks = []
        constraints = []
        for name, _, label_matrices in blocks:
            moment_matrix = sum(
                face_moments[ordered_labels.index(label)] * matrix
                for label, matrix in label_matrices.items()
            )
            eigenvalues, eigenvectors = np.linalg.eigh(moment_matrix)
            if eigenvalues[0] < -100 * args.face_threshold:
                raise ValueError(
                    f"Numerical face {face_path} is not PSD on block "
                    f"{name}: {eigenvalues[0]}"
                )
            kernel = eigenvectors[
                :,
                eigenvalues <= args.face_threshold,
            ]
            if kernel.shape[1] == 0:
                continue
            reduced_matrices = {
                label: kernel.T @ matrix @ kernel
                for label, matrix in label_matrices.items()
            }
            reduced_matrices = {
                label: matrix
                for label, matrix in reduced_matrices.items()
                if np.max(np.abs(matrix)) > 1e-13
            }
            if not reduced_matrices:
                continue
            reduced_variable = cp.Variable(
                (kernel.shape[1], kernel.shape[1]),
                symmetric=True,
                name=f"{name}_numerical_face_{face_index}",
            )
            reduced_blocks.append(
                (
                    name,
                    reduced_variable,
                    reduced_matrices,
                )
            )
            constraints.append(reduced_variable >> 0)
        blocks = reduced_blocks

    if args.check_onb:
        quotient_count = 0
        equality_generator_count = 0
        equality_rank = 0
        if args.eliminate_free:
            quotient_rows, equality_generator_count, equality_rank = (
                exact_equality_quotient_rows(
                    ordered_labels,
                    (
                        label_matrices
                        for _, _, label_matrices in free_blocks
                    ),
                    (
                        relation
                        for _, relation in (
                            gradient_relations
                            + potential_relations
                            + rank_relations
                        )
                    ),
                )
            )
            quotient_count = len(quotient_rows)
        onb_values = {
            label: float(onb_label_value(label))
            for label in ordered_labels
        }
        block_checks = {}
        worst_block_eigenvalue = float("inf")
        for name, _, label_matrices in blocks:
            matrix = sum(
                onb_values[label] * coefficient_matrix
                for label, coefficient_matrix in label_matrices.items()
            )
            minimum_eigenvalue = float(np.linalg.eigvalsh(matrix)[0])
            worst_block_eigenvalue = min(
                worst_block_eigenvalue,
                minimum_eigenvalue,
            )
            block_checks[name] = minimum_eigenvalue
        free_checks = {}
        worst_free_residual = 0.0
        for name, _, label_matrices in free_blocks:
            matrix = sum(
                onb_values[label] * coefficient_matrix
                for label, coefficient_matrix in label_matrices.items()
            )
            residual = float(np.max(np.abs(matrix)))
            worst_free_residual = max(worst_free_residual, residual)
            free_checks[name] = residual
        relation_checks = {}
        worst_relation_residual = 0.0
        for family, relations in (
            ("gradient", gradient_relations),
            ("potential", potential_relations),
            ("rank", rank_relations),
        ):
            residuals = [
                abs(
                    sum(
                        coefficient * onb_values[label]
                        for label, coefficient in relation.items()
                    )
                )
                for _, relation in relations
            ]
            family_residual = max(residuals, default=0.0)
            relation_checks[family] = family_residual
            worst_relation_residual = max(
                worst_relation_residual,
                family_residual,
            )
        return {
            "target": sum(
                coefficient * onb_values[label]
                for label, coefficient in target.items()
            ),
            "labels": len(ordered_labels),
            "block_sizes": {
                name: variable.shape[0]
                for name, variable, _ in blocks
            },
            "free_block_sizes": {
                name: variable.shape[0]
                for name, variable, _ in free_blocks
            },
            "relation_counts": {
                "gradient": len(gradient_relations),
                "potential": len(potential_relations),
                "rank": len(rank_relations),
            },
            "equality_generators": equality_generator_count,
            "equality_rank": equality_rank,
            "quotient_constraints": quotient_count,
            "worst_block_eigenvalue": worst_block_eigenvalue,
            "worst_free_residual": worst_free_residual,
            "worst_relation_residual": worst_relation_residual,
            "block_checks": block_checks,
            "free_checks": free_checks,
            "relation_checks": relation_checks,
        }

    if args.find_face:
        label_indices = {
            label: index for index, label in enumerate(ordered_labels)
        }
        face_moments = cp.Variable(
            len(ordered_labels),
            name="face_moments",
        )
        face_constraints: list[cp.Constraint] = []
        trace_terms: list[cp.Expression] = []
        for _, _, label_matrices in blocks:
            if not label_matrices:
                continue
            raw_block_scale = max(
                np.max(np.abs(matrix))
                for matrix in label_matrices.values()
            )
            block_scale = (
                raw_block_scale if args.scale_constraints else 1.0
            )
            matrix = sum(
                face_moments[label_indices[label]]
                * (coefficient_matrix / block_scale)
                for label, coefficient_matrix in label_matrices.items()
            )
            face_constraints.append(matrix >> 0)
            trace_terms.append(cp.trace(matrix))
        for _, _, label_matrices in free_blocks:
            raw_block_scale = max(
                np.max(np.abs(matrix))
                for matrix in label_matrices.values()
            )
            block_scale = (
                raw_block_scale if args.scale_constraints else 1.0
            )
            matrix = sum(
                face_moments[label_indices[label]]
                * (coefficient_matrix / block_scale)
                for label, coefficient_matrix in label_matrices.items()
            )
            face_constraints.append(matrix == 0)
        for _, relation in (
            gradient_relations + potential_relations + rank_relations
        ):
            relation_scale = (
                max(abs(coefficient) for coefficient in relation.values())
                if args.scale_constraints
                else 1.0
            )
            face_constraints.append(
                sum(
                    (coefficient / relation_scale)
                    * face_moments[label_indices[label]]
                    for label, coefficient in relation.items()
                )
                == 0
            )
        face_constraints.extend(
            [
                sum(
                    coefficient * face_moments[label_indices[label]]
                    for label, coefficient in target.items()
                )
                == 0,
                sum(trace_terms) == 1,
            ]
        )
        face_problem = cp.Problem(
            cp.Minimize(cp.sum_squares(face_moments)),
            face_constraints,
        )
        face_value = face_problem.solve(
            solver="MOSEK",
            verbose=args.verbose,
            mosek_params={
                "MSK_DPAR_INTPNT_CO_TOL_PFEAS": args.tolerance,
                "MSK_DPAR_INTPNT_CO_TOL_DFEAS": args.tolerance,
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": args.tolerance,
            },
        )
        moment_values = {}
        block_ranks = {}
        if face_moments.value is not None:
            moment_values = {
                str(label): float(face_moments.value[index])
                for index, label in enumerate(ordered_labels)
                if abs(face_moments.value[index]) > 1e-9
            }
            for name, _, label_matrices in blocks:
                matrix = sum(
                    face_moments.value[label_indices[label]]
                    * coefficient_matrix
                    for label, coefficient_matrix in label_matrices.items()
                )
                eigenvalues = np.linalg.eigvalsh(matrix)
                block_ranks[name] = {
                    "size": len(eigenvalues),
                    "minimum_eigenvalue": float(eigenvalues[0]),
                    "rank_at_1e-7": int(
                        np.count_nonzero(eigenvalues > 1e-7)
                    ),
                }
            if args.output:
                np.savez(
                    Path(args.output),
                    moments=face_moments.value,
                    labels=np.array(
                        [repr(label) for label in ordered_labels],
                        dtype=str,
                    ),
                )
        return {
            "status": face_problem.status,
            "objective": (
                None if face_value is None else float(face_value)
            ),
            "labels": len(ordered_labels),
            "block_ranks": block_ranks,
            "moments": {} if args.summary_only else moment_values,
            "output": args.output,
        }

    if args.dual:
        label_indices = {
            label: index for index, label in enumerate(ordered_labels)
        }
        moments = cp.Variable(len(ordered_labels), name="moments")
        dual_constraints: list[cp.Constraint] = [
            moments[label_indices[("constant",)]] == 1
        ]
        if args.box_bounds:
            for label, index in label_indices.items():
                if label == ("constant",):
                    continue
                dual_constraints.extend(
                    [
                        moments[index] >= -1,
                        moments[index] <= 1,
                    ]
                )
        for _, _, label_matrices in blocks:
            if not label_matrices:
                continue
            raw_block_scale = max(
                np.max(np.abs(coefficient_matrix))
                for coefficient_matrix in label_matrices.values()
            )
            if raw_block_scale < 1.0e-12:
                continue
            block_scale = (
                raw_block_scale
                if args.scale_constraints
                else 1.0
            )
            matrix = sum(
                moments[label_indices[label]]
                * (coefficient_matrix / block_scale)
                for label, coefficient_matrix in label_matrices.items()
            )
            dual_constraints.append(matrix >> 0)
        for _, _, label_matrices in free_blocks:
            if not label_matrices:
                continue
            raw_block_scale = max(
                np.max(np.abs(coefficient_matrix))
                for coefficient_matrix in label_matrices.values()
            )
            if raw_block_scale < 1.0e-12:
                continue
            block_scale = (
                raw_block_scale
                if args.scale_constraints
                else 1.0
            )
            matrix = sum(
                moments[label_indices[label]]
                * (coefficient_matrix / block_scale)
                for label, coefficient_matrix in label_matrices.items()
            )
            dual_constraints.append(matrix == 0)
        for _, relation in gradient_relations + potential_relations + rank_relations:
            relation_scale = (
                max(abs(coefficient) for coefficient in relation.values())
                if args.scale_constraints
                else 1.0
            )
            dual_constraints.append(
                sum(
                    (coefficient / relation_scale)
                    * moments[label_indices[label]]
                    for label, coefficient in relation.items()
                )
                == 0
            )
        dual_objective = sum(
            coefficient * moments[label_indices[label]]
            for label, coefficient in target.items()
        )
        dual_problem = cp.Problem(
            cp.Minimize(dual_objective),
            dual_constraints,
        )
        dual_mosek_params: dict[str, object] = {
            "MSK_DPAR_INTPNT_CO_TOL_PFEAS": args.tolerance,
            "MSK_DPAR_INTPNT_CO_TOL_DFEAS": args.tolerance,
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": args.tolerance,
        }
        if args.mosek_solve_form != "free":
            dual_mosek_params["MSK_IPAR_INTPNT_SOLVE_FORM"] = (
                f"MSK_SOLVE_{args.mosek_solve_form.upper()}"
            )
        dual_value = dual_problem.solve(
            solver="MOSEK",
            verbose=args.verbose,
            mosek_params=dual_mosek_params,
        )
        moment_values = {}
        minimum_block_eigenvalue = None
        maximum_free_residual = None
        maximum_relation_residual = None
        if moments.value is not None:
            moment_values = {
                str(label): float(moments.value[index])
                for index, label in enumerate(ordered_labels)
                if abs(moments.value[index]) > 1e-8
            }
            block_eigenvalues = []
            for _, _, label_matrices in blocks:
                if not label_matrices:
                    continue
                matrix = sum(
                    moments.value[label_indices[label]] * coefficient_matrix
                    for label, coefficient_matrix in label_matrices.items()
                )
                block_eigenvalues.append(float(np.linalg.eigvalsh(matrix)[0]))
            minimum_block_eigenvalue = min(block_eigenvalues, default=0.0)

            free_residuals = []
            for _, _, label_matrices in free_blocks:
                if not label_matrices:
                    continue
                matrix = sum(
                    moments.value[label_indices[label]] * coefficient_matrix
                    for label, coefficient_matrix in label_matrices.items()
                )
                free_residuals.append(float(np.max(np.abs(matrix))))
            maximum_free_residual = max(free_residuals, default=0.0)

            relation_residuals = []
            for _, relation in (
                gradient_relations + potential_relations + rank_relations
            ):
                relation_residuals.append(
                    abs(
                        sum(
                            coefficient
                            * moments.value[label_indices[label]]
                            for label, coefficient in relation.items()
                        )
                    )
                )
            maximum_relation_residual = max(
                relation_residuals,
                default=0.0,
            )
        return {
            "status": dual_problem.status,
            "objective": None if dual_value is None else float(dual_value),
            "degree": args.degree,
            "target_epsilon": args.target_epsilon,
            "labels": len(ordered_labels),
            "scale_constraints": args.scale_constraints,
            "minimum_block_eigenvalue": minimum_block_eigenvalue,
            "maximum_free_residual": maximum_free_residual,
            "maximum_relation_residual": maximum_relation_residual,
            "moments": {} if args.summary_only else moment_values,
        }

    bounded_labels = [
        label for label in ordered_labels if label != ("constant",)
    ]
    box_plus = (
        cp.Variable(len(bounded_labels), nonneg=True, name="box_plus")
        if args.box_bounds
        else None
    )
    box_minus = (
        cp.Variable(len(bounded_labels), nonneg=True, name="box_minus")
        if args.box_bounds
        else None
    )
    bounded_label_indices = {
        label: index for index, label in enumerate(bounded_labels)
    }

    quotient_rows: list[dict[Label, Fraction]] | None = None
    equality_generator_count = 0
    equality_rank = 0
    if args.eliminate_free:
        if args.box_bounds:
            raise ValueError(
                "Exact equality elimination is incompatible with box bounds"
            )
        quotient_rows, equality_generator_count, equality_rank = (
            exact_equality_quotient_rows(
                ordered_labels,
                (
                    label_matrices
                    for _, _, label_matrices in free_blocks
                ),
                (
                    relation
                    for _, relation in (
                        gradient_relations
                        + potential_relations
                        + rank_relations
                    )
                ),
            )
        )
        for quotient_row in quotient_rows:
            terms: list[cp.Expression] = []
            projected_target = sum(
                float(coefficient) * target.get(label, 0.0)
                for label, coefficient in quotient_row.items()
            )
            equation_scale = abs(projected_target)
            for _, variable, label_matrices in blocks:
                projected_matrix: np.ndarray | None = None
                for label, coefficient in quotient_row.items():
                    matrix = label_matrices.get(label)
                    if matrix is None:
                        continue
                    term = float(coefficient) * matrix
                    projected_matrix = (
                        term
                        if projected_matrix is None
                        else projected_matrix + term
                    )
                if projected_matrix is None:
                    continue
                matrix_scale = float(np.max(np.abs(projected_matrix)))
                if matrix_scale <= 1e-13:
                    continue
                terms.append(
                    cp.sum(cp.multiply(projected_matrix, variable))
                )
                equation_scale = max(equation_scale, matrix_scale)
            if args.scale_constraints and equation_scale > 0:
                constraints.append(
                    sum(terms) / equation_scale
                    == projected_target / equation_scale
                )
            else:
                constraints.append(sum(terms) == projected_target)
    else:
        for label in ordered_labels:
            terms = []
            equation_scale = abs(target.get(label, 0.0))
            for _, variable, label_matrices in blocks:
                matrix = label_matrices.get(label)
                if matrix is not None:
                    terms.append(cp.sum(cp.multiply(matrix, variable)))
                    equation_scale = max(
                        equation_scale,
                        float(np.max(np.abs(matrix))),
                    )
            for _, variable, label_matrices in free_blocks:
                matrix = label_matrices.get(label)
                if matrix is not None:
                    terms.append(cp.sum(cp.multiply(matrix, variable)))
                    equation_scale = max(
                        equation_scale,
                        float(np.max(np.abs(matrix))),
                    )
            if gradient_coefficients is not None:
                for index, (_, relation) in enumerate(gradient_relations):
                    coefficient = relation.get(label)
                    if coefficient is not None:
                        terms.append(
                            coefficient * gradient_coefficients[index]
                        )
                        equation_scale = max(
                            equation_scale,
                            abs(coefficient),
                        )
            if potential_coefficients is not None:
                for index, (_, relation) in enumerate(potential_relations):
                    coefficient = relation.get(label)
                    if coefficient is not None:
                        terms.append(
                            coefficient * potential_coefficients[index]
                        )
                        equation_scale = max(
                            equation_scale,
                            abs(coefficient),
                        )
            # Rank-three Gram identities are unrestricted averaging-null
            # terms.
            if rank_coefficients is not None:
                for index, (_, relation) in enumerate(rank_relations):
                    coefficient = relation.get(label)
                    if coefficient is not None:
                        terms.append(
                            coefficient * rank_coefficients[index]
                        )
                        equation_scale = max(
                            equation_scale,
                            abs(coefficient),
                        )
            if box_plus is not None and box_minus is not None:
                if label == ("constant",):
                    terms.append(cp.sum(box_plus + box_minus))
                    equation_scale = max(equation_scale, 1.0)
                else:
                    index = bounded_label_indices[label]
                    terms.append(box_plus[index] - box_minus[index])
                    equation_scale = max(equation_scale, 1.0)
            if args.scale_constraints and equation_scale > 0:
                constraints.append(
                    sum(terms) / equation_scale
                    == target.get(label, 0.0) / equation_scale
                )
            else:
                constraints.append(
                    sum(terms) == target.get(label, 0.0)
                )

    objective_terms: list[cp.Expression] = [
        cp.trace(variable) for _, variable, _ in blocks
    ]
    if box_plus is not None and box_minus is not None:
        objective_terms.append(cp.sum(box_plus + box_minus))
    objective = cp.Minimize(sum(objective_terms))
    problem = cp.Problem(objective, constraints)
    mosek_params = {
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS": args.tolerance,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS": args.tolerance,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": args.tolerance,
    }
    if args.mosek_solve_form != "free":
        mosek_params["MSK_IPAR_INTPNT_SOLVE_FORM"] = (
            f"MSK_SOLVE_{args.mosek_solve_form.upper()}"
        )
    value = problem.solve(
        solver="MOSEK",
        verbose=args.verbose,
        mosek_params=mosek_params,
    )

    result: dict[str, object] = {
        "status": problem.status,
        "objective": None if value is None else float(value),
        "degree": args.degree,
        "target_epsilon": args.target_epsilon,
        "gram_module": args.gram_module,
        "pointwise_sos": not args.no_pointwise_sos,
        "gradient": args.gradient,
        "hessian": args.hessian,
        "three_point_flags": args.three_point_flags,
        "four_point_flags": args.four_point_flags,
        "two_root_flags": args.two_root_flags,
        "max_flag_arity": args.max_flag_arity,
        "max_root_factor_degree": args.max_root_factor_degree,
        "higher_rank_matrices": args.higher_rank_matrices,
        "harmonics": args.harmonics,
        "four_point_hessian": args.four_point_hessian,
        "max_hessian_arity": args.max_hessian_arity,
        "potential": args.potential,
        "potential_matrices": args.potential_matrices,
        "global_gap": args.global_gap,
        "global_tangent_gaps": args.global_tangent_gaps,
        "labels": len(labels),
        "gradient_powers": [power for power, _ in gradient_relations],
        "potential_powers": [power for power, _ in potential_relations],
        "rank_relations": len(rank_relations),
        "pole_equator_faces": pole_equator_faces,
        "numerical_faces": [str(path) for path in numerical_face_paths],
        "eliminate_free": args.eliminate_free,
        "equality_generators": equality_generator_count,
        "equality_rank": equality_rank,
        "quotient_constraints": (
            0 if quotient_rows is None else len(quotient_rows)
        ),
        "box_bounds": args.box_bounds,
        "blocks": {},
    }
    block_results: dict[str, object] = {}
    for name, variable, _ in blocks:
        if variable.value is None:
            block_results[name] = {"size": variable.shape[0]}
            continue
        eigenvalues = np.linalg.eigvalsh(variable.value)
        block_results[name] = {
            "size": variable.shape[0],
            "minimum_eigenvalue": float(eigenvalues[0]),
            "maximum_eigenvalue": float(eigenvalues[-1]),
            "rank_at_1e-7": int(np.count_nonzero(eigenvalues > 1e-7)),
        }
    result["blocks"] = block_results
    result["free_blocks"] = {
        name: {"size": variable.shape[0]}
        for name, variable, _ in free_blocks
    }
    if gradient_coefficients is not None and gradient_coefficients.value is not None:
        result["gradient_coefficients"] = gradient_coefficients.value.tolist()
    if potential_coefficients is not None and potential_coefficients.value is not None:
        result["potential_coefficients"] = potential_coefficients.value.tolist()
    if box_plus is not None and box_plus.value is not None:
        active_box_terms = {}
        for index, label in enumerate(bounded_labels):
            plus_value = float(box_plus.value[index])
            minus_value = float(box_minus.value[index])
            if plus_value > 1e-8 or minus_value > 1e-8:
                active_box_terms[str(label)] = {
                    "plus": plus_value,
                    "minus": minus_value,
                }
        result["active_box_terms"] = active_box_terms

    if args.output and problem.status in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        output_path = Path(args.output)
        arrays = {
            name: variable.value
            for name, variable, _ in blocks
            if variable.value is not None
        }
        arrays.update(
            {
                name: variable.value
                for name, variable, _ in free_blocks
                if variable.value is not None
            }
        )
        if gradient_coefficients is not None and gradient_coefficients.value is not None:
            arrays["gradient_coefficients"] = gradient_coefficients.value
        if potential_coefficients is not None and potential_coefficients.value is not None:
            arrays["potential_coefficients"] = potential_coefficients.value
        if box_plus is not None and box_plus.value is not None:
            arrays["box_plus"] = box_plus.value
            arrays["box_minus"] = box_minus.value
        np.savez(output_path, **arrays)
        result["output"] = str(output_path)

    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--gram-module", action="store_true")
    parser.add_argument("--no-pointwise-sos", action="store_true")
    parser.add_argument("--gradient", action="store_true")
    parser.add_argument("--hessian", action="store_true")
    parser.add_argument("--three-point-flags", action="store_true")
    parser.add_argument("--h2-localized-flags", action="store_true")
    parser.add_argument("--four-point-flags", action="store_true")
    parser.add_argument("--two-root-flags", action="store_true")
    parser.add_argument("--max-flag-arity", type=int, default=0)
    parser.add_argument("--max-root-factor-degree", type=int, default=0)
    parser.add_argument("--higher-rank-matrices", action="store_true")
    parser.add_argument("--harmonics", action="store_true")
    parser.add_argument("--four-point-hessian", action="store_true")
    parser.add_argument("--max-hessian-arity", type=int, default=0)
    parser.add_argument("--potential", action="store_true")
    parser.add_argument("--potential-matrices", action="store_true")
    parser.add_argument("--global-gap", action="store_true")
    parser.add_argument("--global-tangent-gaps", action="store_true")
    parser.add_argument("--rank-relations", action="store_true")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument(
        "--mosek-solve-form",
        choices=("free", "primal", "dual"),
        default="free",
    )
    parser.add_argument("--target-epsilon", type=float, default=0.0)
    parser.add_argument("--output")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--dual", action="store_true")
    parser.add_argument("--find-face", action="store_true")
    parser.add_argument("--box-bounds", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--scale-constraints", action="store_true")
    parser.add_argument("--jacobi-scale-blocks", action="store_true")
    parser.add_argument("--export-sdpa")
    parser.add_argument("--sdpa-digits", type=int, default=50)
    parser.add_argument("--check-onb", action="store_true")
    parser.add_argument("--facial-reduce-onb", action="store_true")
    parser.add_argument("--exact-onb-face", action="store_true")
    parser.add_argument(
        "--pole-equator-faces",
        help=(
            "Comma-separated exact equality faces: continuous and/or "
            "regular equator orders such as 4,5"
        ),
    )
    parser.add_argument(
        "--numerical-faces",
        help="Comma-separated NPZ exposing directions from --find-face",
    )
    parser.add_argument("--face-threshold", type=float, default=1e-7)
    parser.add_argument("--eliminate-free", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(solve(parse_args()), indent=2, sort_keys=True))
