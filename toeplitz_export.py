#!/usr/bin/env python3
"""GMP-ready export: weighted degree-14 problem + Jensen/Toeplitz blocks.

Exports sdpa_runs/deg14_h2w_h2all_toep.dat-s: the KKT-inclusive weighted
degree-14 problem (same nine-toggle set + --h2-weighted-target
--h2-localized-all as sdpa_runs/deg14_h2w_h2all.dat-s) with the
averaging-contraction (conditional-Jensen / Toeplitz) PSD matrix
families of toeplitz_blocks.py adjoined as first-class blocks.

Soundness note (why this is NOT a post-hoc file append): the base
export eliminates equalities exactly and then drops 337 directions
whose images on the base blocks are linearly dependent.  A block
appended to the finished file would see only the kept directions --
silently pinning the dropped coordinates and potentially cutting
genuine measures.  Instead this driver patches
sos_search.export_sdpa_problem, hands the extra families to the
original exporter, and lets the exact elimination / image-selection
machinery redo the bookkeeping with the new blocks' images included
(m grows as previously-collapsed directions become distinguishable).
Everything stays exact end to end: the family entries are integers or
halves (asserted), so the float hand-off is reconstructed exactly by
rationalize_float.

Appended families (each M(y) = sum_L A_L y_L >= 0, valid for every
measure, all-measures cone; docs/ENRICHMENTS.md):

  jensen_pair / h2loc_jensen_pair / h2comp_gram_pair
  jensen_even_00 / jensen_even_11            (T - G, cap 7)
  h2loc_jensen_even_00 / h2loc_jensen_even_11
  h2comp_gram_even_00/11, h2comp_gram_odd_01/10   ((1-h2) x Gram)

Families are trimmed (basis elements removed greedily) until every
entry's label lies in the run's own ordered label set; each trimmed
family is verified before appending: exact equality against a direct
rational covariance computation and exact PSD-ness at the ONB measure
and at the orthonormal 4-point cross (h2 = 1/4 there, so localized and
complement families are non-trivially exercised).

Usage:
  .venv/bin/python toeplitz_export.py --out sdpa_runs/deg14_h2w_h2all_toep.dat-s
  .venv/bin/python toeplitz_export.py --selectors   # after the export
"""

from __future__ import annotations

import argparse
import json
import pickle
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

import sos_search
from toeplitz_blocks import (
    build_fiber_toeplitz_family,
    build_h2_complement_family,
    build_pair_complement_family,
    build_pair_family,
    build_pair_hankel_localized,
    build_pair_weighted_jensen,
    build_two_root_family,
    direct_family_matrix,
    family_name,
)

Label = tuple

def base_flags(degree: int, cone: str = "kkt") -> list[str]:
    """Base-problem flags per cone.

    kkt     : the KKT-inclusive nine-toggle weighted problem
              (deg14/16_h2w_h2all pattern);
    am      : the proof-carrying all-measures cone (no KKT toggles);
    am_we1  : all-measures + the weighted-(E1) two-layer projection
              (deg18_h2w_h2all_am_we1 pattern; docs/SHARP_STRUCTURE.md).
    """
    flags = [
        "--degree", str(degree), "--no-pointwise-sos",
        "--harmonics", "--three-point-flags", "--four-point-flags",
        "--two-root-flags",
    ]
    if cone == "kkt":
        flags += [
            "--gradient", "--potential", "--hessian",
            "--global-tangent-gaps",
        ]
    flags += [
        "--rank-relations",
        "--h2-weighted-target", "--h2-localized-all",
    ]
    if cone == "am_we1":
        # "both layers projected": pure layer by the weighted-(E1)
        # bases, h2loc layer by the unweighted (E1) bases
        # (docs/SHARP_STRUCTURE.md sections 6 and 8).
        flags += [
            "--e1-project",
            f"sdpa_runs/e1w_projection_deg{degree}.json",
            "--e1-project-localized",
            f"sdpa_runs/e1_projection_deg{degree}.json",
        ]
    flags += ["--summary-only"]
    return flags


# t0 = -eps - 2/3, 40 significant digits (REGEN_NOTES.md A4);
# 5em5 means eps = (16/3) 1e-5 = 1/18750 (t0 terminates exactly).
SELECTOR_BOUNDS = {
    "1em3": "-6.676666666666666666666666666666666666667E-1",
    "1em4": "-6.667666666666666666666666666666666666667E-1",
    "5em5": "-6.6672E-1",
    "2em5": "-6.666866666666666666666666666666666666667E-1",
    "1em5": "-6.666766666666666666666666666666666666667E-1",
    "5em6": "-6.666716666666666666666666666666666666667E-1",
    "1em6": "-6.666676666666666666666666666666666666667E-1",
    "5em7": "-6.666671666666666666666666666666666666667E-1",
}
SELECTOR_SETS = {
    ("deg14", "v1", "kkt"): ["1em3", "1em4"],
    ("deg14", "v2", "kkt"): ["1em3", "1em4", "5em5"],
    ("deg14", "v3", "kkt"): ["1em3", "1em4", "5em5", "2em5"],
    ("deg16", "v1", "kkt"): ["1em4", "1em5"],
    ("deg16", "v3", "kkt"): ["1em4", "1em5", "5em6"],
    ("deg16", "v3", "am"): ["1em4", "1em5", "5em6"],
    ("deg18", "v3", "am_we1"): ["1em5", "1em6", "5em7"],
}

SCRATCH = Path(
    "/private/tmp/claude-501/-Users-jcpaik-Documents-research-"
    "p2-kkt-flag-sos/2d5da291-4f4d-44a8-bb7b-d3ca80702a32/scratchpad"
)


def capture_path(degree: int, version: str) -> Path:
    return SCRATCH / f"toeplitz_capture_deg{degree}_h2w_{version}.pkl"


# Backward-compatible alias for the v1 deg-14 capture (toeplitz_ab.py).
CAPTURE_PATH = SCRATCH / "toeplitz_capture_deg14_h2w.pkl"


# ---------------------------------------------------------------------------
# Exact evaluation on rational atomic measures (verification)
# ---------------------------------------------------------------------------

CROSS_POINTS = [
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(-1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(-1), Fraction(0)),
]
CROSS_WEIGHTS = [Fraction(1, 4)] * 4

ONB_POINTS = [
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1)),
    (Fraction(-1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(-1), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(-1)),
]
ONB_WEIGHTS = [Fraction(1, 6)] * 6


def dot(u, v) -> Fraction:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def det3(u, v, w) -> Fraction:
    return (
        u[0] * (v[1] * w[2] - v[2] * w[1])
        - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0])
    )


def exact_label_moment(label: Label, points, weights) -> Fraction:
    import itertools

    if label == ("constant",):
        return Fraction(1)
    if label[0] == "product":
        value = Fraction(1)
        for factor in label[1:]:
            value *= exact_label_moment(factor, points, weights)
        return value
    if label[0] == "pair":
        vertex_count, exponents = 2, (label[1],)
    elif label[0] == "triangle":
        vertex_count, exponents = 3, tuple(label[1:])
    else:
        vertex_count = int(label[0].split("_")[1])
        exponents = tuple(label[1:])
    edges = [
        (left, right)
        for left in range(vertex_count)
        for right in range(left + 1, vertex_count)
    ]
    total = Fraction(0)
    for combo in itertools.product(
        range(len(points)), repeat=vertex_count
    ):
        term = Fraction(1)
        for vertex in combo:
            term *= weights[vertex]
        for (left, right), exponent in zip(edges, exponents):
            if exponent:
                term *= dot(
                    points[combo[left]], points[combo[right]]
                ) ** exponent
        total += term
    return total


def exact_direct_family_matrix(family: dict, points, weights):
    """Exact-rational direct evaluation of the family matrix."""
    from theta_atoms import poly3_mul as _poly3_mul

    count = len(points)
    gram = [
        [dot(points[a], points[b]) for b in range(count)]
        for a in range(count)
    ]
    p2 = sum(
        weights[x] * weights[y] * gram[x][y] ** 2
        for x in range(count)
        for y in range(count)
    )
    h2 = (3 * p2 - 1) / 2

    def poly3_eval_exact(polynomial, t1, t2, s) -> Fraction:
        return sum(
            (
                value * t1**i * t2**j * s**k
                for (i, j, k), value in polynomial.items()
            ),
            Fraction(0),
        )

    if family["kind"] == "fiber_toeplitz":
        from toeplitz_blocks import ABS_W2, w_even_power_re

        order = family["order"]
        a_polys = w_even_power_re(order)
        absw2_pow = [{(0, 0, 0): Fraction(1)}]
        for _ in range(order):
            absw2_pow.append(_poly3_mul(absw2_pow[-1], ABS_W2))
        indices = family["basis"]
        g_basis = family["g_basis"]
        size = len(indices)
        matrix = [[Fraction(0)] * size for _ in range(size)]
        for x1 in range(count):
            for x2 in range(count):
                s = gram[x1][x2]
                rho = weights[x1] * weights[x2]
                for y in range(count):
                    t1, t2 = gram[x1][y], gram[x2][y]
                    wy = rho * weights[y]
                    for row, (j, ga) in enumerate(indices):
                        for column, (k, gb) in enumerate(indices):
                            m = abs(j - k)
                            gi, gj, gk = g_basis[ga]
                            hi, hj, hk = g_basis[gb]
                            value = (
                                poly3_eval_exact(
                                    absw2_pow[order - m], t1, t2, s
                                )
                                * poly3_eval_exact(
                                    a_polys[m], t1, t2, s
                                )
                                * t1 ** (gi + hi)
                                * t2 ** (gj + hj)
                                * s ** (gk + hk)
                            )
                            matrix[row][column] += wy * value
        if family["h2loc"]:
            matrix = [[h2 * v for v in row] for row in matrix]
        return matrix

    if family["kind"] in ("pair_hankel_loc", "pair_jensen_minor"):
        degrees = family["basis"]
        moment = {}
        for power in range(0, 2 * max(degrees) + 5):
            moment[power] = sum(
                weights[x] * weights[y] * gram[x][y] ** power
                for x in range(count)
                for y in range(count)
            )
        size = len(degrees)
        matrix = [[Fraction(0)] * size for _ in range(size)]
        for row, a in enumerate(degrees):
            for column, b in enumerate(degrees):
                if family["kind"] == "pair_hankel_loc":
                    matrix[row][column] = (
                        moment[a + b] - moment[a + b + 2]
                    )
                else:
                    matrix[row][column] = (
                        moment[a + b]
                        - 2 * moment[a + b + 2]
                        + moment[a + b + 4]
                        - (moment[a] - moment[a + 2])
                        * (moment[b] - moment[b + 2])
                    )
        if family["h2loc"]:
            matrix = [[h2 * v for v in row] for row in matrix]
        return matrix
    if family["kind"] in ("pair_jensen", "h2_complement_pair"):
        degrees = family["basis"]
        pair_values = [
            sum(
                weights[x] * weights[y] * gram[x][y] ** d
                for x in range(count)
                for y in range(count)
            )
            for d in degrees
        ]
        second = [
            [
                sum(
                    weights[x] * weights[y] * gram[x][y] ** (a + b)
                    for x in range(count)
                    for y in range(count)
                )
                for b in degrees
            ]
            for a in degrees
        ]
        size = len(degrees)
        cov = [
            [
                second[r][c] - pair_values[r] * pair_values[c]
                for c in range(size)
            ]
            for r in range(size)
        ]
        outer = [
            [pair_values[r] * pair_values[c] for c in range(size)]
            for r in range(size)
        ]
    else:
        basis = family["basis"]
        odd = family["sector"].startswith("odd")
        size = len(basis)
        cov = [[Fraction(0)] * size for _ in range(size)]
        outer = [[Fraction(0)] * size for _ in range(size)]
        for x1 in range(count):
            for x2 in range(count):
                s = gram[x1][x2]
                phi = [[Fraction(0)] * count for _ in range(size)]
                for y in range(count):
                    t1, t2 = gram[x1][y], gram[x2][y]
                    base = (
                        det3(points[x1], points[x2], points[y])
                        if odd
                        else Fraction(1)
                    )
                    for index, (i, j, k) in enumerate(basis):
                        phi[index][y] = base * t1**i * t2**j * s**k
                mean = [
                    sum(phi[a][y] * weights[y] for y in range(count))
                    for a in range(size)
                ]
                rho = weights[x1] * weights[x2]
                if family["minor"]:
                    rho *= 1 - s * s
                if family.get("s2"):
                    rho *= s * s
                for a in range(size):
                    for b in range(size):
                        second_ab = sum(
                            phi[a][y] * phi[b][y] * weights[y]
                            for y in range(count)
                        )
                        cov[a][b] += rho * (
                            second_ab - mean[a] * mean[b]
                        )
                        outer[a][b] += rho * mean[a] * mean[b]
    if family["kind"] in ("two_root_jensen", "pair_jensen"):
        matrix = cov
        if family["h2loc"]:
            matrix = [[h2 * value for value in row] for row in matrix]
        return matrix
    base = outer if family["which"] == "G" else cov
    return [[(1 - h2) * value for value in row] for row in base]


def assemble_exact(family: dict, points, weights):
    size = family["size"]
    out = [[Fraction(0)] * size for _ in range(size)]
    for label, matrix in family["A"].items():
        value = exact_label_moment(label, points, weights)
        if value == 0:
            continue
        for row in range(size):
            for column in range(size):
                out[row][column] += value * matrix[row][column]
    return out


def exact_psd(matrix) -> bool:
    """Exact rational PSD test by pivoted Gaussian elimination."""
    size = len(matrix)
    work = [[Fraction(value) for value in row] for row in matrix]
    active = list(range(size))
    while active:
        if any(work[i][i] < 0 for i in active):
            return False
        zeros = [i for i in active if work[i][i] == 0]
        for i in zeros:
            if any(work[i][j] != 0 for j in active):
                return False
        active = [i for i in active if work[i][i] > 0]
        if not active:
            return True
        pivot = max(active, key=lambda i: work[i][i])
        d = work[pivot][pivot]
        active = [i for i in active if i != pivot]
        for i in active:
            for j in active:
                work[i][j] -= work[i][pivot] * work[pivot][j] / d
    return True


# ---------------------------------------------------------------------------
# Families and trimming
# ---------------------------------------------------------------------------

def restrict_family(family: dict, keep: list[int]) -> dict:
    out = dict(family)
    out["basis"] = [family["basis"][index] for index in keep]
    out["size"] = len(keep)
    for corner in ("T", "G", "A"):
        matrices = {}
        for label, matrix in family[corner].items():
            restricted = [[matrix[r][c] for c in keep] for r in keep]
            if any(any(v for v in row) for row in restricted):
                matrices[label] = restricted
        out[corner] = matrices
    return out


def trim_to_labels(family: dict, allowed: set) -> tuple[dict, int]:
    dropped = 0
    while True:
        missing = [
            label for label in family["A"] if label not in allowed
        ]
        if not missing:
            return family, dropped
        counts = [0] * family["size"]
        for label in missing:
            matrix = family["A"][label]
            for row in range(family["size"]):
                if any(
                    matrix[row][c] for c in range(family["size"])
                ):
                    counts[row] += 1
        worst = max(
            range(family["size"]), key=lambda i: (counts[i], i)
        )
        family = restrict_family(
            family,
            [i for i in range(family["size"]) if i != worst],
        )
        dropped += 1


def export_families(degree: int = 14, version: str = "v1") -> list[dict]:
    """The exported family list.

    v1 (the GMP-validated set: pole law 10.63x -> 9.59x/decade):
    pair triple + even Jensen (plain, h2loc) + h2comp_gram all sectors.
    v2 adds, per the T2 candidate tables (all measured with matrix
    kill signs; docs/ENRICHMENTS.md section 7): minor Jensen
    (plain + h2loc), h2comp_cov, s^2-localized Jensen (collision
    boundary; plain + h2loc), h2comp_gram minors (even + odd), and
    the cap-5 odd-sector Jensen blocks (plain + h2loc).
    Caps scale with the degree: cap = degree // 2, minors cap - 1.
    """
    cap = degree // 2
    families = [
        build_pair_family(cap),
        build_pair_family(cap, h2loc=True),
        build_pair_complement_family(cap, "G"),
        build_two_root_family("even_00", cap),
        build_two_root_family("even_11", cap),
        build_two_root_family("even_00", cap, h2loc=True),
        build_two_root_family("even_11", cap, h2loc=True),
        build_h2_complement_family("even_00", cap, which="G"),
        build_h2_complement_family("even_11", cap, which="G"),
        build_h2_complement_family("odd_01", cap, which="G"),
        build_h2_complement_family("odd_10", cap, which="G"),
    ]
    if version in ("v2", "v3"):
        minor_cap = cap - 1
        odd_cap = cap - 2
        for sector in ("even_00", "even_11"):
            families += [
                build_two_root_family(sector, minor_cap, minor=True),
                build_two_root_family(
                    sector, minor_cap, minor=True, h2loc=True
                ),
                build_h2_complement_family(
                    sector, minor_cap, which="A"
                ),
                build_two_root_family(sector, minor_cap, s2=True),
                build_two_root_family(
                    sector, minor_cap, s2=True, h2loc=True
                ),
                build_h2_complement_family(
                    sector, minor_cap, minor=True, which="G"
                ),
            ]
        for sector in ("odd_01", "odd_10"):
            families += [
                build_two_root_family(sector, odd_cap),
                build_two_root_family(sector, odd_cap, h2loc=True),
                build_h2_complement_family(
                    sector, minor_cap, minor=True, which="G"
                ),
            ]
    if version == "v3":
        # Fiber-Toeplitz blocks (trigonometric moment matrices of the
        # leaf's azimuthal fiber; docs/ENRICHMENTS.md section 8)
        # and the pair-sector moment families, targeting the v2
        # residual (p2 x triangle 50%, triangle 17%, pair products 28%).
        # Entry degree of FT(K, r) is 4K + 2r <= degree; at degree 16
        # the K = 4 tower opens up.
        if degree >= 18:
            ft_specs = [
                (4, 0, "even_00"),
                (3, 3, "even_00"),
                (3, 3, "even_11"),
                (2, 5, "even_00"),
                (2, 5, "even_11"),
            ]
        elif degree >= 16:
            ft_specs = [
                (4, 0, "even_00"),
                (3, 2, "even_00"),
                (3, 2, "even_11"),
                (2, 4, "even_00"),
                (2, 4, "even_11"),
            ]
        else:
            ft_specs = [
                (3, 0, "even_00"),
                (2, cap - 4, "even_00"),
                (2, cap - 4, "even_11"),
            ]
        for h2loc in (False, True):
            families += [
                build_fiber_toeplitz_family(
                    order, radial_cap, sector, h2loc=h2loc
                )
                for order, radial_cap, sector in ft_specs
            ]
            families += [
                build_pair_hankel_localized(cap - 1, h2loc=h2loc),
                build_pair_weighted_jensen(cap - 2, h2loc=h2loc),
            ]
    return families


def verify_family(family: dict) -> None:
    for points, weights, name in (
        (ONB_POINTS, ONB_WEIGHTS, "ONB"),
        (CROSS_POINTS, CROSS_WEIGHTS, "cross"),
    ):
        assembled = assemble_exact(family, points, weights)
        direct = exact_direct_family_matrix(family, points, weights)
        if assembled != direct:
            raise SystemExit(
                f"exact mismatch for {family_name(family)} at {name}"
            )
        if not exact_psd(assembled):
            raise SystemExit(
                f"{family_name(family)} not PSD at {name} (bug)"
            )
    float_points = np.array(
        [[float(v) for v in p] for p in CROSS_POINTS]
    )
    float_direct = direct_family_matrix(
        family, float_points, np.array([0.25] * 4)
    )
    exact_direct = exact_direct_family_matrix(
        family, CROSS_POINTS, CROSS_WEIGHTS
    )
    gap = max(
        abs(float(exact_direct[r][c]) - float_direct[r, c])
        for r in range(family["size"])
        for c in range(family["size"])
    )
    if gap > 1e-9:
        raise SystemExit(
            f"float/exact direct mismatch for {family_name(family)}"
        )


def family_float_matrices(family: dict) -> dict:
    """Exactly float-representable label matrices (asserted)."""
    matrices = {}
    for label, matrix in family["A"].items():
        array = np.zeros((family["size"], family["size"]))
        for row in range(family["size"]):
            for column in range(family["size"]):
                value = matrix[row][column]
                if value == 0:
                    continue
                recovered = Fraction(
                    float(value)
                ).limit_denominator(10**9)
                if recovered != value:
                    raise SystemExit(
                        f"entry {value} of {family_name(family)} not "
                        f"float-exact (denominator "
                        f"{value.denominator})"
                    )
                array[row, column] = float(value)
        matrices[label] = array
    return matrices


# ---------------------------------------------------------------------------
# Interception driver
# ---------------------------------------------------------------------------

def run_export(
    out_path: str,
    degree: int,
    version: str,
    extra_flags: list[str] | None = None,
    cone: str = "kkt",
    reference_map: str | None = None,
) -> None:
    original = sos_search.export_sdpa_problem
    report: dict = {}
    tag = version + ("" if cone == "kkt" else f"_{cone}") + (
        "_" + "".join(f.strip("-") for f in extra_flags)
        if extra_flags
        else ""
    )
    save_path = capture_path(degree, tag)
    reference_blocks = None
    if reference_map:
        with open(reference_map) as handle:
            reference_blocks = [
                (entry["name"], entry["size"])
                for entry in json.load(handle)["blocks"]
            ]

    def patched(
        path, digits, target, ordered_labels,
        psd_blocks, free_label_matrices, relations,
    ):
        allowed = set(ordered_labels)
        extra_blocks = []
        family_records = []
        print(
            f"[toeplitz] base problem: {len(psd_blocks)} blocks, "
            f"{len(ordered_labels)} labels",
            file=sys.stderr,
        )
        if reference_blocks is not None:
            built = [
                (
                    name,
                    next(iter(matrices.values())).shape[0],
                )
                for name, matrices in psd_blocks
            ]
            if sorted(built) != sorted(reference_blocks):
                raise SystemExit(
                    "[toeplitz] base-block mismatch against "
                    f"{reference_map}: rebuilt "
                    f"{len(built)} blocks vs reference "
                    f"{len(reference_blocks)} — wrong flags, aborting"
                )
            print(
                "[toeplitz] base blocks match the reference map "
                "(names and sizes)",
                file=sys.stderr,
            )
        for family in export_families(degree, version):
            before = family["size"]
            family, dropped = trim_to_labels(family, allowed)
            if family["size"] == 0:
                print(
                    f"[toeplitz] {family_name(family)} trimmed away "
                    f"entirely, skipped",
                    file=sys.stderr,
                )
                continue
            verify_family(family)
            extra_blocks.append(
                (family_name(family), family_float_matrices(family))
            )
            family_records.append(
                {
                    "name": family_name(family),
                    "size": family["size"],
                    "trimmed_from": before,
                    "labels": len(family["A"]),
                    "basis": [
                        list(b) if isinstance(b, tuple) else b
                        for b in family["basis"]
                    ],
                }
            )
            print(
                f"[toeplitz] verified {family_name(family)} "
                f"(size {family['size']}"
                + (f", trimmed from {before}" if dropped else "")
                + ")",
                file=sys.stderr,
            )
        with open(save_path, "wb") as handle:
            pickle.dump(
                {
                    "target": target,
                    "ordered_labels": ordered_labels,
                    "psd_blocks": psd_blocks,
                    "free_label_matrices": free_label_matrices,
                    "relations": relations,
                    "extra_blocks": extra_blocks,
                },
                handle,
            )
        print(
            f"[toeplitz] capture saved to {save_path}",
            file=sys.stderr,
        )
        report["families"] = family_records
        return original(
            path, digits, target, ordered_labels,
            list(psd_blocks) + extra_blocks,
            free_label_matrices, relations,
        )

    sos_search.export_sdpa_problem = patched
    try:
        argv_backup = sys.argv
        sys.argv = (
            ["sos_search.py", "--export-sdpa", out_path]
            + base_flags(degree, cone)
            + list(extra_flags or [])
        )
        try:
            result = sos_search.solve(sos_search.parse_args())
        finally:
            sys.argv = argv_backup
    finally:
        sos_search.export_sdpa_problem = original

    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    sidecar_path = Path(out_path + ".families.json")
    sidecar_path.write_text(
        json.dumps(
            {
                "note": (
                    "Jensen/Toeplitz averaging-contraction blocks "
                    "adjoined first-class; docs/ENRICHMENTS.md"
                    ".md.  bound = objValPrimal + 2/3."
                ),
                "version": version,
                "degree": degree,
                "cone": cone,
                "families": report.get("families", []),
                "base_flags": base_flags(degree, cone)
                + list(extra_flags or []),
            },
            indent=1,
        )
    )
    print(f"wrote {sidecar_path}")


def run_selectors(
    out_path: str, degree: int, version: str, cone: str = "kkt"
) -> None:
    prefix = "sel_toep" if degree == 14 else f"sel{degree}_toep"
    if cone == "am":
        prefix = prefix.replace("_toep", "_am_toep")
    elif cone == "am_we1":
        prefix = prefix.replace("_toep", "_we1_toep")
    if version == "v2":
        prefix = prefix.replace("toep", "toep2")
    elif version == "v3":
        prefix = prefix.replace("toep", "toep3")
    for tag in SELECTOR_SETS[(f"deg{degree}", version, cone)]:
        bound = SELECTOR_BOUNDS[tag]
        selector_path = str(
            Path(out_path).with_name(f"{prefix}_{tag}.dat-s")
        )
        completed = subprocess.run(
            [
                sys.executable, "sdpa_selector.py", out_path,
                selector_path, f"--bound={bound}",
            ],
            capture_output=True, text=True,
        )
        print(completed.stdout.strip())
        if completed.returncode != 0:
            print(completed.stderr, file=sys.stderr)
            raise SystemExit(f"selector build failed for {tag}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out")
    parser.add_argument("--degree", type=int, default=14)
    parser.add_argument(
        "--version", choices=["v1", "v2", "v3"], default="v1"
    )
    parser.add_argument("--selectors", action="store_true")
    parser.add_argument(
        "--cone", choices=["kkt", "am", "am_we1"], default="kkt"
    )
    parser.add_argument(
        "--reference-map",
        help=(
            "abort unless the rebuilt base blocks match this "
            "existing export's .map.json (names and sizes)"
        ),
    )
    parser.add_argument(
        "--extra-flag",
        action="append",
        default=[],
        help=(
            "additional sos_search flags for the base problem "
            "(e.g. --extra-flag=--gap-cut-e5 for the e5-cut variant)"
        ),
    )
    args = parser.parse_args()
    default_out = {
        (14, "v1", "kkt"): "sdpa_runs/deg14_h2w_h2all_toep.dat-s",
        (14, "v2", "kkt"): "sdpa_runs/deg14_h2w_h2all_toep2.dat-s",
        (14, "v3", "kkt"): "sdpa_runs/deg14_h2w_h2all_toep3.dat-s",
        (16, "v1", "kkt"): "sdpa_runs/deg16_h2w_h2all_toep.dat-s",
        (16, "v3", "kkt"): "sdpa_runs/deg16_h2w_h2all_toep3.dat-s",
        (16, "v3", "am"): "sdpa_runs/deg16_h2w_h2all_am_toep3.dat-s",
        (18, "v3", "am_we1"): "sdpa_runs/deg18_we1_toep3.dat-s",
    }[(args.degree, args.version, args.cone)]
    out_path = args.out or default_out
    if args.selectors:
        run_selectors(out_path, args.degree, args.version, args.cone)
        return
    run_export(
        out_path, args.degree, args.version, args.extra_flag,
        args.cone, args.reference_map,
    )


if __name__ == "__main__":
    main()
