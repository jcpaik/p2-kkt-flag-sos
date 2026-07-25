#!/usr/bin/env python3
"""MOSEK search for isotropic KKT-infused flag/SOS certificates.

Independent samples from an antipodally symmetric isotropic measure on S^2
are represented by multigraph monomials in their pairwise Gram entries.  A
degree-two sampled vertex is contracted using E[XX^T] = I/3; every remaining
expectation is retained as an independent canonical moment label.  The code
therefore makes no unrecorded moment-closure assumption.

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
    if power == 2:
        return ("constant",), factor / 3
    return ("pair", power), factor


def expectation_label(exponent: Exponent) -> tuple[Label | None, Fraction]:
    """Reduce one triangle monomial using antipodality and isotropy.

    The degree at X, Y, Z is respectively i+k, i+j, j+k.  Odd vertex degree
    has zero expectation.  Vertex degree zero reduces to a pair moment, while
    vertex degree two can be integrated using E[XX^T] = I/3.  Remaining
    triangle moments are retained as independent labels, modulo S_3 symmetry.
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

    if i + k == 2:
        power = j + 1 if (i, k) == (1, 1) else j
        return pair_label(power, Fraction(1, 3))
    if i + j == 2:
        power = k + 1 if (i, j) == (1, 1) else k
        return pair_label(power, Fraction(1, 3))
    if j + k == 2:
        power = i + 1 if (j, k) == (1, 1) else i
        return pair_label(power, Fraction(1, 3))

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


def reduce_graph_matrix(matrix: list[list[int]]) -> tuple[Label | None, Fraction]:
    """Evaluate all consequences of antipodality and E[XX^T]=I/3."""

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

    degree_two_vertices = [
        vertex for vertex, degree in enumerate(degrees) if degree == 2
    ]
    if degree_two_vertices:
        vertex = degree_two_vertices[0]
        incident = [
            (other, exponent)
            for other, exponent in enumerate(matrix[vertex])
            if exponent
        ]
        retained = [other for other in range(len(matrix)) if other != vertex]
        reduced = induced_matrix(matrix, retained)
        old_to_new = {
            old: new for new, old in enumerate(retained)
        }
        if len(incident) == 2:
            (left, left_exponent), (right, right_exponent) = incident
            if left_exponent != 1 or right_exponent != 1:
                raise ValueError("Unexpected degree-two contraction")
            reduced_left = old_to_new[left]
            reduced_right = old_to_new[right]
            reduced[reduced_left][reduced_right] += 1
            reduced[reduced_right][reduced_left] += 1
        elif len(incident) != 1 or incident[0][1] != 2:
            raise ValueError("Unexpected degree-two contraction")
        label, coefficient = reduce_graph_matrix(reduced)
        return label, coefficient / 3

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
    isotropic_energy_terms = [
        (Fraction(16, 3), ("constant",)),
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
                for energy_coefficient, energy_label in isotropic_energy_terms:
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


def solve(args: argparse.Namespace) -> dict[str, object]:
    target = {
        ("constant",): 1.0 + args.target_epsilon,
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
        for degree in range(4, args.degree + 1, 2):
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

    if args.global_gap:
        variable = cp.Variable((1, 1), symmetric=True, name="global_uniform_gap")
        blocks.append(
            (
                "global_uniform_gap",
                variable,
                {
                    ("constant",): np.array([[-176.0 / 35.0]]),
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
            onb_matrix = sum(
                float(onb_label_value(label)) * coefficient_matrix
                for label, coefficient_matrix in label_matrices.items()
            )
            eigenvalues, eigenvectors = np.linalg.eigh(onb_matrix)
            if eigenvalues[0] < -1e-8:
                raise ValueError(
                    f"Block {name} is not PSD on the ONB: {eigenvalues[0]}"
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

    if args.check_onb:
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
            "worst_block_eigenvalue": worst_block_eigenvalue,
            "worst_free_residual": worst_free_residual,
            "worst_relation_residual": worst_relation_residual,
            "block_checks": block_checks,
            "free_checks": free_checks,
            "relation_checks": relation_checks,
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

    for label in ordered_labels:
        terms: list[cp.Expression] = []
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
                    terms.append(coefficient * gradient_coefficients[index])
                    equation_scale = max(equation_scale, abs(coefficient))
        if potential_coefficients is not None:
            for index, (_, relation) in enumerate(potential_relations):
                coefficient = relation.get(label)
                if coefficient is not None:
                    terms.append(coefficient * potential_coefficients[index])
                    equation_scale = max(equation_scale, abs(coefficient))
        # Rank-three Gram identities are unrestricted averaging-null terms.
        if rank_coefficients is not None:
            for index, (_, relation) in enumerate(rank_relations):
                coefficient = relation.get(label)
                if coefficient is not None:
                    terms.append(coefficient * rank_coefficients[index])
                    equation_scale = max(equation_scale, abs(coefficient))
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
            constraints.append(sum(terms) == target.get(label, 0.0))

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
    parser.add_argument("--box-bounds", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--scale-constraints", action="store_true")
    parser.add_argument("--check-onb", action="store_true")
    parser.add_argument("--facial-reduce-onb", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(solve(parse_args()), indent=2, sort_keys=True))
