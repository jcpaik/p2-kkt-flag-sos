#!/usr/bin/env python3
"""Averaging-contraction (conditional-Jensen / Toeplitz) flag blocks.

Implements the mathematical program of docs/TOEPLITZ_BLOCKS_NOTE.md
(mission "Agent T": translate the circle-pair Toeplitz mechanism of
docs/CYLINDRICAL_DOMINATION.md section 6 into valid flag-algebra blocks).

The object.  For a vector phi = (phi_alpha) of leaf functions of the
rooted three-sample configuration (x1, x2, y) -- polynomials
B(t1, t2, s) with t1 = x1.y, t2 = x2.y, s = x1.x2, optionally times the
orientation det(x1, x2, y) -- and a nonnegative root weight
rho(x1, x2) in {1, 1 - s^2}, define

    T_ab = E[ rho * phi_a(y) phi_b(y) ]        (same leaf; 3-point labels)
    G_ab = E[ rho * phi_a(y) phi_b(y') ]       (independent leaves; the
                                                existing two_root Gram)

Then, writing Cov_y for the leaf-conditional covariance given the roots,

    T - G = E_{x1,x2}[ rho * Cov_y(phi) ]  >=  0   (PSD, every measure).

This single matrix inequality is the exact flag transcription of the
cylindrical Toeplitz coupling (D5): restricted to azimuthal-mode vectors
it contains  |zeta_{2k}| >= 2 |zeta_k|^2 - 1  (the double-angle identity
cos^2(k phi) = (1 + cos 2k phi)/2 is a polynomial identity in
(t1, t2, s) and is applied automatically by the moment reducer inside
the same-leaf entries of T), and with polynomial test multipliers it is
the measure-valued matrix domination

    [[sigma, mu_k, mu_2k], [., sigma, mu_k], [., ., sigma]] >= 0

conditioned on the sampled pair frame instead of a fixed axis.
Equivalently: the mixed Gram over instant + averaged flags
U = [[T, G], [G, G]] = E[f f^T], f = (phi(y), E_y phi), is PSD, and
U >= 0  <=>  G >= 0 and T - G >= 0 (Schur; G(I - G G^+) = 0 holds
automatically).  G >= 0 is already in the hierarchy; T - G is new.

The unrooted (empty-type) analogue takes the whole pair (X, Y) as the
sample: T_ab = p_{a+b} (pair labels), G_ab = p_a p_b (product labels),
so T - G is the covariance matrix of the random variable (X.Y)^d --
new coupling between the pair sector and the pair-product sector.

h2-localized copies multiply by the scalar h2(mu) = (3 p2 - 1)/2 >= 0
(valid for every measure since p2 = tr Sigma^2 >= 1/3), shifting labels
by p2-products: exactly the p2 x (...) sectors where the measured
weighted escape lives (docs/UNPROJECTED_ESCAPE_NOTE.md).

Everything exact-rational; self-tests compare the label expansions with
direct numerical integration on random atomic antipodal measures
(cylinder_cert.py verify style) and check PSD-ness of T - G there.

CLI:
  --self-test          exact + numeric validation (run before trusting)
  --pair-d PATH        M(D) eigenvalue tables against a fingerprint JSON
                       (default: sdpa_runs/fingerprint_D_e3e4.json)
  --coverage           label-coverage report against a blocks dump
  --dump-json PATH     exact expansions of all family matrices
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction

import numpy as np

import theta_atoms
from theta_atoms import (
    Poly3,
    hat_generator,
    label_moment,
    poly3_add,
    poly3_degree,
    poly3_eval,
    poly3_mul,
    random_antipodal_measure,
)
from sos_search import (
    ORIENTATION_PAIRING,
    ROOT_PAIR_MINOR,
    expectation_label,
    graph_expectation_label,
    monomials,
    multiply_graph_polynomials,
    multiply_labels,
)

Label = tuple
GraphPolynomial = list

ONE_POLY3: Poly3 = {(0, 0, 0): Fraction(1)}
ONE_MINUS_S2: Poly3 = {(0, 0, 0): Fraction(1), (0, 0, 2): Fraction(-1)}
# det(x1, x2, y)^2 = Gram determinant of the triple.
GRAM_DET: Poly3 = {
    (0, 0, 0): Fraction(1),
    (1, 1, 1): Fraction(2),
    (2, 0, 0): Fraction(-1),
    (0, 2, 0): Fraction(-1),
    (0, 0, 2): Fraction(-1),
}
ONE_PAIRING: GraphPolynomial = [(Fraction(1), (0, 0, 0, 0, 0, 0))]

PARITY = {
    "even_00": lambda i, j, k: (i + j) % 2 == 0
    and (i + k) % 2 == 0
    and (j + k) % 2 == 0,
    "even_11": lambda i, j, k: (i + j) % 2 == 0
    and (i + k) % 2 == 1
    and (j + k) % 2 == 1,
    "odd_01": lambda i, j, k: (i + j) % 2 == 1
    and (i + k) % 2 == 0
    and (j + k) % 2 == 1,
    "odd_10": lambda i, j, k: (i + j) % 2 == 1
    and (i + k) % 2 == 1
    and (j + k) % 2 == 0,
}


def sector_basis(sector: str, max_degree: int) -> list[tuple[int, int, int]]:
    """Leaf monomial basis of a two-root parity sector, in solver order."""
    return [e for e in monomials(max_degree) if PARITY[sector](*e)]


# ---------------------------------------------------------------------------
# Exact entries
# ---------------------------------------------------------------------------

def same_leaf_entry(
    left: Poly3, right: Poly3, weight: Poly3
) -> dict[Label, Fraction]:
    """E[weight(t1,t2,s) * left * right] with ONE leaf: 3-point labels.

    Variable map (X = x1, Y = y, Z = x2): a = X.Y = t1 (exponent i),
    b = Y.Z = t2 (exponent j), c = Z.X = s (exponent k), so a Poly3 key
    (i, j, k) feeds sos_search.expectation_label((i, j, k)) directly.
    """
    out: dict[Label, Fraction] = {}
    product = poly3_mul(poly3_mul(left, right), weight)
    for key, coefficient in product.items():
        label, reduction = expectation_label(key)
        if label is None or reduction == 0:
            continue
        value = out.get(label, Fraction(0)) + coefficient * reduction
        if value:
            out[label] = value
        else:
            out.pop(label, None)
    return out


def two_leaf_entry(
    left: Poly3, right: Poly3, pairing: GraphPolynomial
) -> dict[Label, Fraction]:
    """E[left(t1,t2,s) right(t1',t2',s) * pairing]: 4-point labels.

    Vertex order (X=0, Z=1, Y=2, W=3), edges (XZ, XY, XW, ZY, ZW, YW);
    identical to sos_search.two_root_flag_expectation_matrix and
    theta_atoms.expand_q_block, but exact-rational.
    """
    out: dict[Label, Fraction] = {}
    for (i, j, k), a in left.items():
        for (p, r, t), b in right.items():
            base = (k + t, i, p, j, r, 0)
            for coefficient, shift in pairing:
                exponent = tuple(
                    base[index] + shift[index] for index in range(6)
                )
                label, reduction = graph_expectation_label(4, exponent)
                if label is None or reduction == 0:
                    continue
                value = (
                    out.get(label, Fraction(0))
                    + a * b * coefficient * reduction
                )
                if value:
                    out[label] = value
                else:
                    out.pop(label, None)
    return out


def h2_shift_matrices(
    matrices: dict[Label, list[list[Fraction]]],
) -> dict[Label, list[list[Fraction]]]:
    """Label matrices of h2 * F given those of F (exact h2 = (3p2-1)/2)."""
    p2 = ("pair", 2)
    out: dict[Label, list[list[Fraction]]] = {}
    for label, matrix in matrices.items():
        size = len(matrix)
        for factor, shifted in (
            (Fraction(3, 2), multiply_labels(p2, label)),
            (Fraction(-1, 2), label),
        ):
            target = out.setdefault(
                shifted,
                [[Fraction(0)] * size for _ in range(size)],
            )
            for row in range(size):
                for column in range(size):
                    target[row][column] += factor * matrix[row][column]
    return {
        label: matrix
        for label, matrix in out.items()
        if any(any(value for value in row) for row in matrix)
    }


# ---------------------------------------------------------------------------
# Families
# ---------------------------------------------------------------------------

S2_POLY3: Poly3 = {(0, 0, 2): Fraction(1)}
S2_PAIRING: GraphPolynomial = [(Fraction(1), (2, 0, 0, 0, 0, 0))]


def build_two_root_family(
    sector: str,
    degree_cap: int,
    minor: bool = False,
    h2loc: bool = False,
    s2: bool = False,
) -> dict:
    """Exact T, G, and A = T - G label matrices for one sector family.

    ``s2`` multiplies the root weight by s^2 >= 0 (pointwise), the
    complementary localization to the ``minor`` weight 1 - s^2: it
    concentrates the conditional-Jensen coupling at nearly-collinear
    root pairs -- the collision boundary s = +-1 where the measured
    escape's sign-alternating T_n(+-1) signature lives.
    """
    basis = sector_basis(sector, degree_cap)
    odd = sector.startswith("odd")
    weight = GRAM_DET if odd else ONE_POLY3
    pairing = ORIENTATION_PAIRING if odd else ONE_PAIRING
    if minor:
        weight = poly3_mul(weight, ONE_MINUS_S2)
        pairing = multiply_graph_polynomials(pairing, ROOT_PAIR_MINOR)
    if s2:
        weight = poly3_mul(weight, S2_POLY3)
        pairing = multiply_graph_polynomials(pairing, S2_PAIRING)
    size = len(basis)

    T: dict[Label, list[list[Fraction]]] = {}
    G: dict[Label, list[list[Fraction]]] = {}
    for row, left_exp in enumerate(basis):
        left = {left_exp: Fraction(1)}
        for column, right_exp in enumerate(basis):
            right = {right_exp: Fraction(1)}
            for target, entry in (
                (T, same_leaf_entry(left, right, weight)),
                (G, two_leaf_entry(left, right, pairing)),
            ):
                for label, value in entry.items():
                    matrix = target.setdefault(
                        label,
                        [[Fraction(0)] * size for _ in range(size)],
                    )
                    matrix[row][column] += value
    if h2loc:
        T = h2_shift_matrices(T)
        G = h2_shift_matrices(G)
    A: dict[Label, list[list[Fraction]]] = {}
    for label in set(T) | set(G):
        t_matrix = T.get(label)
        g_matrix = G.get(label)
        matrix = [
            [
                (t_matrix[row][column] if t_matrix else Fraction(0))
                - (g_matrix[row][column] if g_matrix else Fraction(0))
                for column in range(size)
            ]
            for row in range(size)
        ]
        if any(any(value for value in row) for row in matrix):
            A[label] = matrix
    return {
        "kind": "two_root_jensen",
        "sector": sector,
        "degree_cap": degree_cap,
        "minor": minor,
        "h2loc": h2loc,
        "s2": s2,
        "basis": basis,
        "size": size,
        "T": T,
        "G": G,
        "A": A,
    }


def build_pair_family(max_pair_degree: int, h2loc: bool = False) -> dict:
    """Unrooted pair-sample covariance: T_ab = p_{a+b}, G_ab = p_a p_b.

    Sample = one whole pair (X, Y); flags (X.Y)^d, d even.  T - G is the
    covariance matrix of the vector ((X.Y)^d)_d: PSD for every measure.
    Matches empty_type_flag's flag list (pair_degrees = 0, 2, ..).
    """
    degrees = list(range(0, max_pair_degree + 1, 2))
    size = len(degrees)

    def pair_label_of(power: int) -> Label:
        return ("constant",) if power == 0 else ("pair", power)

    T: dict[Label, list[list[Fraction]]] = {}
    G: dict[Label, list[list[Fraction]]] = {}
    for row, a in enumerate(degrees):
        for column, b in enumerate(degrees):
            label = pair_label_of(a + b)
            T.setdefault(
                label, [[Fraction(0)] * size for _ in range(size)]
            )[row][column] += 1
            label, reduction = graph_expectation_label(
                4, (a, 0, 0, 0, 0, b)
            )
            if label is not None and reduction != 0:
                G.setdefault(
                    label, [[Fraction(0)] * size for _ in range(size)]
                )[row][column] += reduction
    if h2loc:
        T = h2_shift_matrices(T)
        G = h2_shift_matrices(G)
    A: dict[Label, list[list[Fraction]]] = {}
    for label in set(T) | set(G):
        t_matrix = T.get(label)
        g_matrix = G.get(label)
        matrix = [
            [
                (t_matrix[row][column] if t_matrix else Fraction(0))
                - (g_matrix[row][column] if g_matrix else Fraction(0))
                for column in range(size)
            ]
            for row in range(size)
        ]
        if any(any(value for value in row) for row in matrix):
            A[label] = matrix
    return {
        "kind": "pair_jensen",
        "sector": "empty_type",
        "degree_cap": max_pair_degree,
        "minor": False,
        "h2loc": h2loc,
        "basis": degrees,
        "size": size,
        "T": T,
        "G": G,
        "A": A,
    }


def subtract_matrices(
    left: dict[Label, list[list[Fraction]]],
    right: dict[Label, list[list[Fraction]]],
    size: int,
) -> dict[Label, list[list[Fraction]]]:
    out: dict[Label, list[list[Fraction]]] = {}
    for label in set(left) | set(right):
        l_matrix = left.get(label)
        r_matrix = right.get(label)
        matrix = [
            [
                (l_matrix[row][column] if l_matrix else Fraction(0))
                - (r_matrix[row][column] if r_matrix else Fraction(0))
                for column in range(size)
            ]
            for row in range(size)
        ]
        if any(any(value for value in row) for row in matrix):
            out[label] = matrix
    return out


def build_h2_complement_family(
    sector: str,
    degree_cap: int,
    minor: bool = False,
    which: str = "G",
) -> dict:
    """The (1 - h2)-localized family: plain minus h2-localized copy.

    Valid because h2 = (3 p2 - 1)/2 <= 1 (p2 <= 1), so 1 - h2 >= 0 and
    (1 - h2) * E[..] >= 0 for the PSD kernel averages G (two-root Gram)
    or T - G (conditional covariance).  ``which`` in {"G", "A"}.
    Directly targets the measured plain <-> h2loc seesaw of the escape.
    """
    plain = build_two_root_family(sector, degree_cap, minor=minor)
    localized = build_two_root_family(
        sector, degree_cap, minor=minor, h2loc=True
    )
    size = plain["size"]
    return {
        "kind": "h2_complement",
        "which": which,
        "sector": sector,
        "degree_cap": degree_cap,
        "minor": minor,
        "h2loc": False,
        "basis": plain["basis"],
        "size": size,
        "T": plain["T"],
        "G": plain["G"],
        "A": subtract_matrices(plain[which], localized[which], size),
    }


def build_pair_complement_family(
    max_pair_degree: int, which: str = "G"
) -> dict:
    plain = build_pair_family(max_pair_degree)
    localized = build_pair_family(max_pair_degree, h2loc=True)
    size = plain["size"]
    return {
        "kind": "h2_complement_pair",
        "which": which,
        "sector": "empty_type",
        "degree_cap": max_pair_degree,
        "minor": False,
        "h2loc": False,
        "basis": plain["basis"],
        "size": size,
        "T": plain["T"],
        "G": plain["G"],
        "A": subtract_matrices(plain[which], localized[which], size),
    }


# ---------------------------------------------------------------------------
# v3: fiber-Toeplitz blocks and pair-sector moment families
# ---------------------------------------------------------------------------
# w = (t2 - s t1) + i det(x1,x2,y) = sin(delta) sin(theta) e^{i phi}:
# |w|^2 = (1-t1^2)(1-s^2) and w^{2m} = A_m + i det B_m with A_m, B_m
# polynomial in (t1,t2,s) (det^2 = GRAM_DET).

ABS_W2: Poly3 = {
    (0, 0, 0): Fraction(1),
    (2, 0, 0): Fraction(-1),
    (0, 0, 2): Fraction(-1),
    (2, 0, 2): Fraction(1),
}
RE_W2: Poly3 = poly3_add(
    poly3_mul(
        {(0, 1, 0): Fraction(1), (1, 0, 1): Fraction(-1)},
        {(0, 1, 0): Fraction(1), (1, 0, 1): Fraction(-1)},
    ),
    GRAM_DET,
    Fraction(-1),
)
IM_W2_OVER_DET: Poly3 = {(0, 1, 0): Fraction(2), (1, 0, 1): Fraction(-2)}


def w_even_power_re(maximum: int) -> list[Poly3]:
    """[A_0 .. A_K] with w^{2m} = A_m + i det B_m (A_m = Re w^{2m})."""
    a_list, b_list = [dict(ONE_POLY3)], [dict()]
    for m in range(maximum):
        a_next = poly3_add(
            poly3_mul(a_list[m], RE_W2),
            poly3_mul(poly3_mul(b_list[m], IM_W2_OVER_DET), GRAM_DET),
            Fraction(-1),
        )
        b_next = poly3_add(
            poly3_mul(a_list[m], IM_W2_OVER_DET),
            poly3_mul(b_list[m], RE_W2),
        )
        a_list.append(a_next)
        b_list.append(b_next)
    return a_list


def build_fiber_toeplitz_family(
    order: int,
    radial_cap: int,
    sector: str = "even_00",
    h2loc: bool = False,
) -> dict:
    """Fiber-Toeplitz block: Re of the Hermitian Toeplitz(x)radial Gram.

    With V_{(j,a)} = |w|^K e^{2 i j phi} g_a(t1,t2,s), j = 0..K, the
    matrix E[V V^H] is Hermitian PSD for every measure (single-leaf
    average of a rank-one PSD, conditional on the roots, then averaged);
    its real part is PSD and has POLYNOMIAL entries

        M_{(j,a),(k,b)} = E[ (|w|^2)^{K-|j-k|} Re(w^{2(j-k)}) g_a g_b ],

    i.e. the radially-weighted trigonometric moment matrix of the
    leaf's azimuthal distribution about the pair frame -- Toeplitz
    positivity of the fiber measure, NOT a polynomial square: the
    Fejer-Riesz factors |w|^K e^{2ij phi} are non-polynomial for
    j > K/2, so no within-degree polynomial Gram reproduces the
    constant-diagonal structure (see docs/TOEPLITZ_BLOCKS_NOTE.md).
    Entries are triangle labels; h2loc copies give p2 x triangle.
    ``sector`` selects the per-vertex parity class of the radial
    factors g_a (pairs from one class keep entries leaf-even).
    """
    g_basis = sector_basis(sector, radial_cap)
    a_polys = w_even_power_re(order)
    absw2_pow = [dict(ONE_POLY3)]
    for _ in range(order):
        absw2_pow.append(poly3_mul(absw2_pow[-1], ABS_W2))
    indices = [
        (j, g_index)
        for j in range(order + 1)
        for g_index in range(len(g_basis))
    ]
    size = len(indices)
    A: dict[Label, list[list[Fraction]]] = {}
    for row, (j, ga) in enumerate(indices):
        for column, (k, gb) in enumerate(indices):
            m = abs(j - k)
            core = poly3_mul(absw2_pow[order - m], a_polys[m])
            core = poly3_mul(
                core,
                poly3_mul(
                    {g_basis[ga]: Fraction(1)},
                    {g_basis[gb]: Fraction(1)},
                ),
            )
            for key, coefficient in core.items():
                label, reduction = expectation_label(key)
                if label is None or reduction == 0:
                    continue
                matrix = A.setdefault(
                    label,
                    [[Fraction(0)] * size for _ in range(size)],
                )
                matrix[row][column] += coefficient * reduction
    if h2loc:
        A = h2_shift_matrices(A)
    return {
        "kind": "fiber_toeplitz",
        "sector": sector,
        "degree_cap": radial_cap,
        "order": order,
        "minor": False,
        "h2loc": h2loc,
        "basis": indices,
        "g_basis": g_basis,
        "size": size,
        "T": A,
        "G": {},
        "A": A,
    }


def pair_power_label(power: int) -> Label:
    return ("constant",) if power == 0 else ("pair", power)


def build_pair_hankel_localized(
    max_pair_degree: int = 6, h2loc: bool = False
) -> dict:
    """Localized Hankel of the pair variable: [p_{a+b} - p_{a+b+2}] >= 0.

    The moment matrix of the positive measure (1 - t^2) d nu with nu
    the law of t = X.Y -- valid for every measure, not implied by the
    g_l >= 0 rows or the Gram families (localization at t in [-1,1]).
    """
    degrees = list(range(0, max_pair_degree + 1, 2))
    size = len(degrees)
    A: dict[Label, list[list[Fraction]]] = {}
    for row, a in enumerate(degrees):
        for column, b in enumerate(degrees):
            for power, coefficient in (
                (a + b, Fraction(1)),
                (a + b + 2, Fraction(-1)),
            ):
                label = pair_power_label(power)
                matrix = A.setdefault(
                    label,
                    [[Fraction(0)] * size for _ in range(size)],
                )
                matrix[row][column] += coefficient
    if h2loc:
        A = h2_shift_matrices(A)
    return {
        "kind": "pair_hankel_loc",
        "sector": "empty_type",
        "degree_cap": max_pair_degree,
        "minor": True,
        "h2loc": h2loc,
        "basis": degrees,
        "size": size,
        "T": A,
        "G": {},
        "A": A,
    }


def build_pair_weighted_jensen(
    max_pair_degree: int = 4, h2loc: bool = False
) -> dict:
    """Covariance of the vector ((1 - t^2) t^a)_a over one pair sample.

    T_ab = p_{a+b} - 2 p_{a+b+2} + p_{a+b+4} (pair labels),
    G_ab = (p_a - p_{a+2})(p_b - p_{b+2})   (pair-product labels);
    T - G >= 0.  The h2-localized copy lands on p2 x pair and
    p2 x pair x pair -- the pair-product home of the v2 residual.
    """
    degrees = list(range(0, max_pair_degree + 1, 2))
    size = len(degrees)
    T: dict[Label, list[list[Fraction]]] = {}
    G: dict[Label, list[list[Fraction]]] = {}
    for row, a in enumerate(degrees):
        for column, b in enumerate(degrees):
            for power, coefficient in (
                (a + b, Fraction(1)),
                (a + b + 2, Fraction(-2)),
                (a + b + 4, Fraction(1)),
            ):
                matrix = T.setdefault(
                    pair_power_label(power),
                    [[Fraction(0)] * size for _ in range(size)],
                )
                matrix[row][column] += coefficient
            for pa, ca in ((a, Fraction(1)), (a + 2, Fraction(-1))):
                for pb, cb in (
                    (b, Fraction(1)),
                    (b + 2, Fraction(-1)),
                ):
                    label = multiply_labels(
                        pair_power_label(pa), pair_power_label(pb)
                    )
                    matrix = G.setdefault(
                        label,
                        [[Fraction(0)] * size for _ in range(size)],
                    )
                    matrix[row][column] += ca * cb
    if h2loc:
        T = h2_shift_matrices(T)
        G = h2_shift_matrices(G)
    A = subtract_matrices(T, G, size)
    return {
        "kind": "pair_jensen_minor",
        "sector": "empty_type",
        "degree_cap": max_pair_degree,
        "minor": True,
        "h2loc": h2loc,
        "basis": degrees,
        "size": size,
        "T": T,
        "G": G,
        "A": A,
    }


def family_name(family: dict) -> str:
    name = family["sector"]
    if family["kind"] == "two_root_jensen":
        name = f"jensen_{name}"
    elif family["kind"] == "h2_complement":
        name = (
            f"h2comp_{'cov' if family['which'] == 'A' else 'gram'}_"
            f"{name}"
        )
    elif family["kind"] == "h2_complement_pair":
        name = (
            f"h2comp_{'cov' if family['which'] == 'A' else 'gram'}_pair"
        )
    elif family["kind"] == "fiber_toeplitz":
        name = f"ftoep{family['order']}_{name}"
        prefix = "h2loc_" if family["h2loc"] else ""
        return prefix + name + f"_r{family['degree_cap']}"
    elif family["kind"] == "pair_hankel_loc":
        name = "pair_hankel_loc"
        prefix = "h2loc_" if family["h2loc"] else ""
        return prefix + name + f"_d{family['degree_cap']}"
    elif family["kind"] == "pair_jensen_minor":
        name = "pair_jensen_minor"
        prefix = "h2loc_" if family["h2loc"] else ""
        return prefix + name + f"_d{family['degree_cap']}"
    else:
        name = "jensen_pair"
    if family["minor"]:
        name += "_minor"
    if family.get("s2"):
        name += "_s2"
    if family["h2loc"]:
        name = "h2loc_" + name
    return name + f"_d{family['degree_cap']}"


def default_families(
    cap: int = 7,
    minor_cap: int = 6,
    include_minor: bool = True,
    include_complements: bool = True,
) -> list[dict]:
    """The h2-localized copies use the same caps as the plain families,
    matching --h2-localized-all (which h2-shifts every block at full
    size; confirmed against the block dump)."""
    families = []
    for h2loc in (False, True):
        families.append(build_pair_family(cap, h2loc=h2loc))
        for sector in ("even_00", "even_11", "odd_01", "odd_10"):
            families.append(
                build_two_root_family(sector, cap, h2loc=h2loc)
            )
        if include_minor:
            for sector in ("even_00", "even_11"):
                families.append(
                    build_two_root_family(
                        sector, minor_cap, minor=True, h2loc=h2loc
                    )
                )
    if include_complements:
        families.append(build_pair_complement_family(cap, "G"))
        for sector in ("even_00", "even_11", "odd_01", "odd_10"):
            families.append(
                build_h2_complement_family(sector, cap, which="G")
            )
        for sector in ("even_00", "even_11"):
            families.append(
                build_h2_complement_family(sector, cap, which="A")
            )
    return families


# ---------------------------------------------------------------------------
# Direct numerical evaluation (validation only)
# ---------------------------------------------------------------------------

def direct_family_matrix(family: dict, points, weights) -> np.ndarray:
    """The family matrix evaluated directly on a discrete measure."""
    count = len(points)
    gram = points @ points.T

    def h2_value() -> float:
        p2 = sum(
            weights[x] * weights[y] * gram[x, y] ** 2
            for x in range(count)
            for y in range(count)
        )
        return (3.0 * p2 - 1.0) / 2.0

    if family["kind"] == "fiber_toeplitz":
        order = family["order"]
        g_basis = family["g_basis"]
        indices = family["basis"]
        size = len(indices)
        matrix = np.zeros((size, size))
        for x1 in range(count):
            for x2 in range(count):
                s = gram[x1, x2]
                vectors = np.zeros((size, count), dtype=complex)
                for y in range(count):
                    t1, t2 = gram[x1, y], gram[x2, y]
                    det = np.linalg.det(
                        np.stack([points[x1], points[x2], points[y]])
                    )
                    w = (t2 - s * t1) + 1j * det
                    magnitude = abs(w)
                    if magnitude < 1e-14:
                        continue
                    phase2 = (w / magnitude) ** 2
                    for row, (j, ga) in enumerate(indices):
                        i_exp, j_exp, k_exp = g_basis[ga]
                        vectors[row, y] = (
                            magnitude**order
                            * phase2**j
                            * t1**i_exp
                            * t2**j_exp
                            * s**k_exp
                        )
                inner = vectors @ np.diag(weights) @ vectors.conj().T
                matrix += weights[x1] * weights[x2] * inner.real
        if family["h2loc"]:
            matrix = h2_value() * matrix
        return matrix

    if family["kind"] in ("pair_hankel_loc", "pair_jensen_minor"):
        degrees = family["basis"]
        omega_moment = {}
        for power in range(0, 2 * max(degrees) + 5):
            omega_moment[power] = sum(
                weights[x] * weights[y] * gram[x, y] ** power
                for x in range(count)
                for y in range(count)
            )
        size = len(degrees)
        matrix = np.zeros((size, size))
        for row, a in enumerate(degrees):
            for column, b in enumerate(degrees):
                if family["kind"] == "pair_hankel_loc":
                    matrix[row, column] = (
                        omega_moment[a + b] - omega_moment[a + b + 2]
                    )
                else:
                    second = (
                        omega_moment[a + b]
                        - 2 * omega_moment[a + b + 2]
                        + omega_moment[a + b + 4]
                    )
                    mean_a = omega_moment[a] - omega_moment[a + 2]
                    mean_b = omega_moment[b] - omega_moment[b + 2]
                    matrix[row, column] = second - mean_a * mean_b
        if family["h2loc"]:
            matrix = h2_value() * matrix
        return matrix

    pair_like = family["kind"] in ("pair_jensen", "h2_complement_pair")
    if pair_like:
        degrees = family["basis"]
        pair_values = np.array(
            [
                sum(
                    weights[x] * weights[y] * gram[x, y] ** d
                    for x in range(count)
                    for y in range(count)
                )
                for d in degrees
            ]
        )
        second = np.array(
            [
                [
                    sum(
                        weights[x] * weights[y] * gram[x, y] ** (a + b)
                        for x in range(count)
                        for y in range(count)
                    )
                    for b in degrees
                ]
                for a in degrees
            ]
        )
        outer = np.outer(pair_values, pair_values)
        cov_matrix, gram_matrix = second - outer, outer
    else:
        basis = family["basis"]
        odd = family["sector"].startswith("odd")
        size = len(basis)
        cov_matrix = np.zeros((size, size))
        gram_matrix = np.zeros((size, size))
        dets = None
        if odd:
            dets = np.zeros((count, count, count))
            for x1 in range(count):
                for x2 in range(count):
                    for y in range(count):
                        dets[x1, x2, y] = np.linalg.det(
                            np.stack(
                                [points[x1], points[x2], points[y]]
                            )
                        )
        for x1 in range(count):
            for x2 in range(count):
                s = gram[x1, x2]
                phi = np.zeros((size, count))
                for y in range(count):
                    t1, t2 = gram[x1, y], gram[x2, y]
                    base = 1.0 if not odd else dets[x1, x2, y]
                    for index, (i, j, k) in enumerate(basis):
                        phi[index, y] = base * t1**i * t2**j * s**k
                mean = phi @ weights
                second = phi @ np.diag(weights) @ phi.T
                rho = weights[x1] * weights[x2]
                if family["minor"]:
                    rho *= 1.0 - s * s
                if family.get("s2"):
                    rho *= s * s
                cov_matrix += rho * (second - np.outer(mean, mean))
                gram_matrix += rho * np.outer(mean, mean)
    p2 = sum(
        weights[x] * weights[y] * gram[x, y] ** 2
        for x in range(count)
        for y in range(count)
    )
    h2 = (3.0 * p2 - 1.0) / 2.0
    if family["kind"] in ("two_root_jensen", "pair_jensen"):
        matrix = cov_matrix
        if family["h2loc"]:
            matrix = h2 * matrix
        return matrix
    base = gram_matrix if family["which"] == "G" else cov_matrix
    return (1.0 - h2) * base


def assemble_from_labels(
    matrices: dict[Label, list[list[Fraction]]], points, weights
) -> np.ndarray:
    size = len(next(iter(matrices.values())))
    out = np.zeros((size, size))
    for label, matrix in matrices.items():
        value = label_moment(label, points, weights)
        out += value * np.array(
            [[float(entry) for entry in row] for row in matrix]
        )
    return out


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------

def self_test() -> None:
    rng = np.random.default_rng(29)
    failures: list[str] = []

    def check(name: str, condition: bool) -> None:
        print(("PASS " if condition else "FAIL ") + name)
        if not condition:
            failures.append(name)

    # 1. small families: exact label expansion == direct integration,
    #    and direct T - G is PSD (validity) at random atomic measures.
    small = [
        build_pair_family(6),
        build_pair_family(6, h2loc=True),
        build_two_root_family("even_00", 3),
        build_two_root_family("even_11", 3),
        build_two_root_family("odd_01", 3),
        build_two_root_family("odd_10", 3),
        build_two_root_family("even_00", 3, minor=True),
        build_two_root_family("even_00", 3, h2loc=True),
        build_two_root_family("odd_10", 3, minor=True, h2loc=True),
        build_two_root_family("even_00", 3, s2=True),
        build_two_root_family("even_11", 3, s2=True, h2loc=True),
        build_h2_complement_family("even_00", 3, which="G"),
        build_h2_complement_family("even_11", 3, which="A"),
        build_h2_complement_family("odd_01", 3, which="G"),
        build_pair_complement_family(6, "G"),
        build_fiber_toeplitz_family(2, 2, "even_00"),
        build_fiber_toeplitz_family(2, 1, "even_11", h2loc=True),
        build_fiber_toeplitz_family(3, 0, "even_00"),
        build_pair_hankel_localized(6),
        build_pair_hankel_localized(6, h2loc=True),
        build_pair_weighted_jensen(4),
        build_pair_weighted_jensen(4, h2loc=True),
    ]
    worst_gap = 0.0
    worst_eig = 0.0
    agree = True
    psd_ok = True
    for family in small:
        for _ in range(3):
            points, weights = random_antipodal_measure(rng, 4)
            direct = direct_family_matrix(family, points, weights)
            expanded = assemble_from_labels(
                family["A"], points, weights
            )
            gap = np.max(np.abs(direct - expanded)) / max(
                1.0, np.max(np.abs(direct))
            )
            worst_gap = max(worst_gap, gap)
            if gap > 1e-9:
                agree = False
            eigenvalues = np.linalg.eigvalsh(direct)
            worst_eig = min(worst_eig, float(eigenvalues[0]))
            if eigenvalues[0] < -1e-10:
                psd_ok = False
    check(
        f"exact A = T - G expansion == direct covariance "
        f"(worst rel gap {worst_gap:.2e})",
        agree,
    )
    check(
        f"T - G is PSD at random atomic measures "
        f"(worst eigenvalue {worst_eig:.2e})",
        psd_ok,
    )

    # 2. Toeplitz sharpness structure.  (i) Coherence: at a single-orbit
    #    measure {+-u} every leaf-even function is a.s. constant on the
    #    leaf, so every even-sector Jensen matrix vanishes -- the flag
    #    image of "|zeta_k| = 1 forces equality in the Toeplitz bound".
    #    (ii) At a generic circle measure the mode-2 direction
    #    C_1 = 2 t1 t2 - s has a strictly positive Jensen gap (the
    #    coupling has teeth).
    from math import cos, pi, sin

    u = np.array([[1.0, 0.0, 0.0]])
    orbit_points = np.vstack([u, -u])
    orbit_weights = np.array([0.5, 0.5])
    coherent_ok = True
    worst_coherent = 0.0
    for sector in ("even_00", "even_11"):
        family = build_two_root_family(sector, 4)
        direct = direct_family_matrix(
            family, orbit_points, orbit_weights
        )
        worst_coherent = max(worst_coherent, np.max(np.abs(direct)))
        if np.max(np.abs(direct)) > 1e-12:
            coherent_ok = False
    check(
        f"even-sector Jensen matrices vanish at the coherent one-orbit "
        f"measure (max |entry| {worst_coherent:.2e})",
        coherent_ok,
    )
    generic = np.array(
        [
            [cos(0.3), sin(0.3), 0.0],
            [cos(1.1), sin(1.1), 0.0],
            [cos(2.0), sin(2.0), 0.0],
        ]
    )
    generic_points = np.vstack([generic, -generic])
    generic_weights = np.array([0.5, 0.3, 0.2, 0.5, 0.3, 0.2]) / 2.0
    family = build_two_root_family("even_11", 4)
    basis_index = {e: i for i, e in enumerate(family["basis"])}
    coefficients = np.zeros(family["size"])
    for exponent, value in (((1, 1, 0), 2.0), ((0, 0, 1), -1.0)):
        coefficients[basis_index[exponent]] = value
    direct = direct_family_matrix(
        family, generic_points, generic_weights
    )
    generic_gap = float(coefficients @ direct @ coefficients)
    check(
        f"mode-2 (C_1) Jensen gap is strictly positive at a generic "
        f"circle measure (gap {generic_gap:.2e})",
        generic_gap > 1e-6,
    )

    # 3. h2-localization: direct h2 * (T - G) == label-shifted expansion
    #    is already covered in 1; here confirm h2 >= 0 on random measures.
    h2_ok = True
    for _ in range(5):
        points, weights = random_antipodal_measure(rng, 5)
        gram = points @ points.T
        p2 = sum(
            weights[x] * weights[y] * gram[x, y] ** 2
            for x in range(len(points))
            for y in range(len(points))
        )
        if 3.0 * p2 - 1.0 < -1e-12:
            h2_ok = False
    check("h2 = (3 p2 - 1)/2 >= 0 at random measures", h2_ok)

    # 4. the hat-generator coefficient vectors lie in the sector bases
    #    (T1(a): the modulated families are inside the existing spans).
    inside = True
    for family_index, orders in ((2, range(-5, 7)), (6, range(-1, 3))):
        for n in orders:
            polynomial = hat_generator(family_index, n)
            sector = (
                "even_00" if n % 2 == 0 else "even_11"
            )
            basis = set(sector_basis(sector, 7))
            if any(key not in basis for key in polynomial):
                inside = False
    check(
        "every within-degree hat C_n / hat S_n lies inside one "
        "two-root parity-sector monomial basis (degree cap 7)",
        inside,
    )

    print()
    if failures:
        raise SystemExit(f"{len(failures)} self-test failure(s): {failures}")
    print("all toeplitz-block self-tests passed")


# ---------------------------------------------------------------------------
# T2: pairing against a fingerprint direction
# ---------------------------------------------------------------------------

def family_pairing(family: dict, direction: dict[str, float]) -> dict:
    """Eigenvalues of M(D) = sum_L D_L A_L, with coverage statistics."""
    size = family["size"]
    matrix = np.zeros((size, size))
    hit_mass = 0.0
    miss_mass = 0.0
    hit = 0
    for label, entry in family["A"].items():
        weight = float(
            max(abs(value) for row in entry for value in row)
        )
        key = str(label)
        if key in direction:
            hit += 1
            hit_mass += weight
            matrix += direction[key] * np.array(
                [[float(value) for value in row] for row in entry]
            )
        else:
            miss_mass += weight
    eigenvalues = np.linalg.eigvalsh(matrix)
    negative = eigenvalues[eigenvalues < -1e-9]
    return {
        "family": family_name(family),
        "size": size,
        "labels": len(family["A"]),
        "labels_in_D": hit,
        "absent_label_matrix_mass": miss_mass,
        "present_label_matrix_mass": hit_mass,
        "min_eig": float(eigenvalues[0]),
        "max_eig": float(eigenvalues[-1]),
        "negative_count": int(negative.size),
        "negative_sum": float(negative.sum()),
        "frobenius": float(np.linalg.norm(matrix)),
        "eigenvalues_bottom5": [float(v) for v in eigenvalues[:5]],
    }


def corner_pairing(
    family: dict, direction: dict[str, float], corner: str
) -> dict:
    """Same eigen test on the T or G corner alone (baselines)."""
    size = family["size"]
    matrix = np.zeros((size, size))
    for label, entry in family[corner].items():
        key = str(label)
        if key in direction:
            matrix += direction[key] * np.array(
                [[float(value) for value in row] for row in entry]
            )
    eigenvalues = np.linalg.eigvalsh(matrix)
    return {
        "corner": corner,
        "min_eig": float(eigenvalues[0]),
        "max_eig": float(eigenvalues[-1]),
    }


def generator_jensen_rows(
    families: dict[str, dict], direction: dict[str, float]
) -> list[dict]:
    """Per-generator 1x1 Jensen pairings <t_hat - g_hat, D> for the
    modulated families (the sharp per-mode instruments)."""
    rows = []
    for family_index, orders in ((2, range(-5, 7)), (6, range(-1, 3))):
        for n in orders:
            polynomial = hat_generator(family_index, n)
            sector = "even_00" if n % 2 == 0 else "even_11"
            for h2loc in (False, True):
                name = ("h2loc_" if h2loc else "") + sector
                family = families.get(name)
                if family is None:
                    continue
                if any(
                    key not in set(family["basis"])
                    for key in polynomial
                ):
                    continue
                index = {e: i for i, e in enumerate(family["basis"])}
                coefficients = np.zeros(family["size"])
                for key, value in polynomial.items():
                    coefficients[index[key]] = float(value)
                total = 0.0
                for label, entry in family["A"].items():
                    key = str(label)
                    if key not in direction:
                        continue
                    entry_value = float(
                        coefficients
                        @ np.array(
                            [
                                [float(v) for v in row]
                                for row in entry
                            ]
                        )
                        @ coefficients
                    )
                    total += direction[key] * entry_value
                rows.append(
                    {
                        "generator": f"{'S' if family_index == 6 else 'C'}"
                        f"hat_{n}",
                        "h2loc": h2loc,
                        "jensen_pairing": total,
                    }
                )
    return rows


def run_pairing(direction_path: str, out_path: str | None) -> None:
    with open(direction_path) as handle:
        payload = json.load(handle)
    direction = payload["D"] if "D" in payload else payload
    print(
        f"direction: {direction_path} ({len(direction)} labels, "
        f"|D|_1 = {sum(abs(v) for v in direction.values()):.4e})"
    )

    print(
        "building exact families (two-root cap 7, minors cap 6, "
        "h2loc at full size matching --h2-localized-all) ..."
    )
    families = default_families()
    by_name: dict[str, dict] = {}
    for family in families:
        key = ("h2loc_" if family["h2loc"] else "") + family["sector"]
        if family["kind"] == "two_root_jensen" and not family["minor"]:
            by_name[key] = family

    print(
        f"\nM(D) = sum_L D_L (T_L - G_L) eigenvalue table "
        f"({len(families)} families):"
    )
    header = (
        f"{'family':>34} {'size':>5} {'min eig':>12} {'max eig':>12} "
        f"{'#neg':>5} {'neg sum':>12} {'miss-mass':>10}"
    )
    print(header)
    results = []
    for family in families:
        row = family_pairing(family, direction)
        row["corners"] = [
            corner_pairing(family, direction, "T"),
            corner_pairing(family, direction, "G"),
        ]
        results.append(row)
        print(
            f"{row['family']:>34} {row['size']:>5} "
            f"{row['min_eig']:>12.4e} {row['max_eig']:>12.4e} "
            f"{row['negative_count']:>5} {row['negative_sum']:>12.4e} "
            f"{row['absent_label_matrix_mass']:>10.3f}"
        )

    print("\ncorner baselines (T alone / G alone), min eigenvalues:")
    for row in results:
        t_corner, g_corner = row["corners"]
        print(
            f"{row['family']:>34}  T: {t_corner['min_eig']:>12.4e}  "
            f"G: {g_corner['min_eig']:>12.4e}"
        )

    print("\nper-generator Jensen pairings <t - g, D> (negative = "
          "the escape violates the conditional variance of that mode):")
    generator_rows = generator_jensen_rows(by_name, direction)
    for row in generator_rows:
        tag = "h2loc " if row["h2loc"] else "plain "
        print(
            f"  {tag}{row['generator']:>8}: "
            f"{row['jensen_pairing']:+.4e}"
        )

    if out_path:
        with open(out_path, "w") as handle:
            json.dump(
                {
                    "direction_file": direction_path,
                    "families": [
                        {
                            key: value
                            for key, value in row.items()
                        }
                        for row in results
                    ],
                    "generator_jensen_rows": generator_rows,
                },
                handle,
                indent=1,
            )
        print(f"\nwrote {out_path}")


# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------

def run_dump_check(blocks_path: str) -> None:
    """Cross-check the exact G corners against the solver's block dump.

    Confirms basis ordering and the h2-localization convention: my
    two-leaf Gram G must equal the dumped two_root_* matrices, and my
    h2-shifted G the dumped h2loc_two_root_* matrices, label by label.
    """
    with open(blocks_path) as handle:
        dump = json.load(handle)["blocks"]
    checks = [
        ("two_root_even_00", "even_00", 7, False, False),
        ("two_root_even_11", "even_11", 7, False, False),
        ("two_root_odd_01", "odd_01", 7, False, False),
        ("two_root_even_00_minor", "even_00", 6, True, False),
        ("h2loc_two_root_even_00", "even_00", 7, False, True),
        ("h2loc_two_root_even_11", "even_11", 7, False, True),
        ("h2loc_two_root_odd_10", "odd_10", 7, False, True),
        ("h2loc_two_root_even_11_minor", "even_11", 6, True, True),
    ]
    all_ok = True
    for name, sector, cap, minor, h2loc in checks:
        if name not in dump:
            print(f"  {name}: not in dump, skipped")
            continue
        family = build_two_root_family(
            sector, cap, minor=minor, h2loc=h2loc
        )
        mine = {str(label): matrix for label, matrix in family["G"].items()}
        theirs = dump[name]
        worst = 0.0
        size_ok = True
        for label, matrix in theirs.items():
            if len(matrix) != family["size"]:
                size_ok = False
                break
            reference = mine.get(label)
            for row in range(len(matrix)):
                for column in range(len(matrix)):
                    mine_value = (
                        float(reference[row][column]) if reference else 0.0
                    )
                    worst = max(
                        worst,
                        abs(mine_value - matrix[row][column]),
                    )
        extra = set(mine) - set(theirs)
        extra_mass = max(
            (
                abs(float(value))
                for label in extra
                for row in family["G"][eval(label)]
                for value in row
            ),
            default=0.0,
        )
        ok = size_ok and worst < 1e-9 and extra_mass < 1e-9
        all_ok = all_ok and ok
        print(
            f"  {name}: size {'ok' if size_ok else 'MISMATCH'}, "
            f"worst entry gap {worst:.2e}, extra-label mass "
            f"{extra_mass:.2e} -> {'PASS' if ok else 'FAIL'}"
        )
    print("dump cross-check " + ("PASSED" if all_ok else "FAILED"))


def run_coverage(blocks_path: str) -> None:
    with open(blocks_path) as handle:
        vocabulary = set(json.load(handle)["labels"])
    print(f"vocabulary: {len(vocabulary)} labels from {blocks_path}")
    for family in default_families():
        labels = {str(label) for label in family["A"]}
        missing = sorted(labels - vocabulary)
        print(
            f"{family_name(family):>34}: {len(labels):>4} labels, "
            f"{len(missing):>3} outside the vocabulary"
        )
        for label in missing[:6]:
            print(f"      missing: {label}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--pair-d",
        nargs="?",
        const="sdpa_runs/fingerprint_D_e3e4.json",
        metavar="FINGERPRINT_JSON",
    )
    parser.add_argument("--out")
    parser.add_argument(
        "--coverage",
        nargs="?",
        const="sdpa_runs/blocks_deg14_h2w_h2all.json",
        metavar="BLOCKS_JSON",
    )
    parser.add_argument(
        "--check-dump",
        nargs="?",
        const="sdpa_runs/blocks_deg14_h2w_h2all.json",
        metavar="BLOCKS_JSON",
    )
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return
    if args.check_dump:
        run_dump_check(args.check_dump)
        return
    if args.coverage:
        run_coverage(args.coverage)
        return
    if args.pair_d:
        run_pairing(args.pair_d, args.out)
        return
    parser.print_help()


if __name__ == "__main__":
    main()
