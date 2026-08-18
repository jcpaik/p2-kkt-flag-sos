"""Exact isotropic audit of the weighted-(E1) projection.

Every pure-square family of a sharp all-measures certificate for the
weighted target W = h2*E must have vanishing expectation at every zero of
W, i.e. at every isotropic antipodal measure.  This script verifies, in
exact rational arithmetic, that the projected bases exported by

    solve_e1.py --export-projection e1w_projection_degD.json --weighted

vanish at isotropic test measures:

  (a) the ONB measure;
  (b) a two-parameter continuous family of rotated-ONB mixtures
      (rational rotations, rational mixture weights);
  (c) a one-parameter continuous pole-equator isotropic family
      (equal-mass angle pairs a quarter turn apart kill nu-hat(2) for
      every angle and every weight split);
  (d) cross mixtures of (b) and (c).

Checks per measure (all exact, `== 0` on Fractions):

  1. every two-root sector vector has zero leaf average at every sampled
     root pair, on-support and off-support (the weighted-(E1) condition
     holds for ALL root pairs);
  2. the deviatoric tensor Dev = M2 - I/3 vanishes -- this is the entire
     projected one-root layer (spin 0 {t^2-1/3}, spin 1 {t}, spin 2 {1})
     and the harmonic_flag weight (t^2 - 1/3);
  3. every projected spin2_flag combination V_c(mu) = 0;
  4. the pair-flag direction 3 p2 - 1 = 0 (also the h2 factor of every
     localized block, so the exempt h2-multiplied layer vanishes too);
  5. g_2(mu) = 0 (the harmonic_2 scalar block).

Negative controls (the audit has teeth):

  6. the unweighted-only leaf T6 + 1/3 (cut by the weighted projection)
     has nonzero average at some isotropic test measure and root;
  7. the unweighted-only two-root generator C^_2 = C_2 + T_2(s)/3
     likewise;
  8. a non-isotropic measure (pole weight 1/2) violates checks 1-5.

Run:  python3 audit_weighted_e1.py [--projection sdpa_runs/e1w_projection_deg14.json]
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction

Vec = tuple[Fraction, Fraction, Fraction]
Measure = list[tuple[Vec, Fraction]]  # projective atoms (antipodalized)


def dot(u: Vec, v: Vec) -> Fraction:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def det3(u: Vec, v: Vec, w: Vec) -> Fraction:
    return (
        u[0] * (v[1] * w[2] - v[2] * w[1])
        - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0])
    )


def rot_z(c: Fraction, s: Fraction):
    return lambda v: (c * v[0] - s * v[1], s * v[0] + c * v[1], v[2])


def rot_x(c: Fraction, s: Fraction):
    return lambda v: (v[0], c * v[1] - s * v[2], s * v[1] + c * v[2])


ONB: Measure = [
    ((Fraction(1), Fraction(0), Fraction(0)), Fraction(1, 3)),
    ((Fraction(0), Fraction(1), Fraction(0)), Fraction(1, 3)),
    ((Fraction(0), Fraction(0), Fraction(1)), Fraction(1, 3)),
]


def rotate_measure(measure: Measure, rotation) -> Measure:
    return [(rotation(v), w) for v, w in measure]


def mix(*parts: tuple[Measure, Fraction]) -> Measure:
    atoms: Measure = []
    for measure, weight in parts:
        atoms.extend((v, w * weight) for v, w in measure)
    return atoms


def pole_equator_isotropic(angle: tuple[Fraction, Fraction], u: Fraction) -> Measure:
    """(1/3) poles + (2/3) equator on {0, phi, pi/2, phi + pi/2}.

    The quarter-turn pairing kills every equator mode m == 2 (mod 4), in
    particular nu-hat(2), for every angle phi and split u, so the measure
    is isotropic for all parameters (nu-hat(4) != 0: it is not in the
    E = 0 family unless the parameters make nu-hat(6) vanish too).
    """

    c, s = angle
    return [
        ((Fraction(0), Fraction(0), Fraction(1)), Fraction(1, 3)),
        ((Fraction(1), Fraction(0), Fraction(0)), Fraction(2, 3) * u / 2),
        ((Fraction(0), Fraction(1), Fraction(0)), Fraction(2, 3) * u / 2),
        ((c, s, Fraction(0)), Fraction(2, 3) * (1 - u) / 2),
        ((-s, c, Fraction(0)), Fraction(2, 3) * (1 - u) / 2),
    ]


def second_moment_deviation(measure: Measure) -> list[list[Fraction]]:
    dev = [[Fraction(0)] * 3 for _ in range(3)]
    for v, w in measure:
        for r in range(3):
            for c in range(3):
                dev[r][c] += w * v[r] * v[c]
    for d in range(3):
        dev[d][d] -= Fraction(1, 3)
    return dev


def pair_power(measure: Measure, power: int) -> Fraction:
    return sum(
        (
            wi * wj * dot(vi, vj) ** power
            for vi, wi in measure
            for vj, wj in measure
        ),
        Fraction(0),
    )


def legendre_pair_energy(measure: Measure, order: int) -> Fraction:
    coefficients = {
        2: [(Fraction(3, 2), 2), (Fraction(-1, 2), 0)],
        4: [
            (Fraction(35, 8), 4),
            (Fraction(-30, 8), 2),
            (Fraction(3, 8), 0),
        ],
        6: [
            (Fraction(231, 16), 6),
            (Fraction(-315, 16), 4),
            (Fraction(105, 16), 2),
            (Fraction(-5, 16), 0),
        ],
    }[order]
    return sum(
        (weight * pair_power(measure, power) for weight, power in coefficients),
        Fraction(0),
    )


def energy(measure: Measure) -> Fraction:
    return (
        Fraction(-4, 3)
        + 20 * pair_power(measure, 2)
        - 48 * pair_power(measure, 4)
        + 32 * pair_power(measure, 6)
    )


def two_root_leaf_average(
    measure: Measure,
    x1: Vec,
    x2: Vec,
    basis: list[tuple[int, int, int]],
    vector: list[Fraction],
    det_sector: bool,
) -> Fraction:
    s_value = dot(x1, x2)
    total = Fraction(0)
    for y, w in measure:
        t1 = dot(x1, y)
        t2 = dot(x2, y)
        value = Fraction(0)
        for (i, j, k), coefficient in zip(basis, vector):
            if coefficient:
                value += coefficient * t1**i * t2**j * s_value**k
        if det_sector:
            value *= det3(x1, x2, y)
        total += w * value
    return total


def spin2_flag_vector(
    measure: Measure, leaves: tuple[int, ...], pair: int
) -> list[list[Fraction]]:
    """V_f(mu) = int (x x^T - I/3) prod (x.y_i)^{a_i} (y1.y2)^b dmu^{1+k}."""

    result = [[Fraction(0)] * 3 for _ in range(3)]
    for x, wx in measure:
        if len(leaves) == 0:
            factor = Fraction(1)
        elif len(leaves) == 1:
            factor = sum(
                (w * dot(x, y) ** leaves[0] for y, w in measure), Fraction(0)
            )
        else:
            factor = sum(
                (
                    w1
                    * w2
                    * dot(x, y1) ** leaves[0]
                    * dot(x, y2) ** leaves[1]
                    * dot(y1, y2) ** pair
                    for y1, w1 in measure
                    for y2, w2 in measure
                ),
                Fraction(0),
            )
        for r in range(3):
            for c in range(3):
                base = x[r] * x[c] - (Fraction(1, 3) if r == c else Fraction(0))
                result[r][c] += wx * base * factor
    return result


ROOT_POINTS: list[Vec] = [
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1)),
    (Fraction(3, 5), Fraction(4, 5), Fraction(0)),
    (Fraction(0), Fraction(3, 5), Fraction(4, 5)),
    (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3)),
    (Fraction(2, 7), Fraction(3, 7), Fraction(6, 7)),
    (Fraction(1, 9), Fraction(4, 9), Fraction(8, 9)),
    (Fraction(2, 11), Fraction(6, 11), Fraction(9, 11)),
]


def audit_measure(
    name: str,
    measure: Measure,
    sectors,
    spin2,
    require_zero: bool = True,
) -> dict[str, object]:
    mass = sum(w for _, w in measure)
    assert mass == 1, f"{name}: mass {mass}"
    root_pairs = [
        (ROOT_POINTS[i], ROOT_POINTS[j])
        for i in range(len(ROOT_POINTS))
        for j in range(i + 1, len(ROOT_POINTS))
    ]
    # include on-support pairs
    atoms = [v for v, _ in measure]
    root_pairs.extend(
        (atoms[i], atoms[j])
        for i in range(len(atoms))
        for j in range(i + 1, len(atoms))
    )

    failures: list[str] = []
    worst_two_root = Fraction(0)
    for sector_name, (basis, vectors) in sectors.items():
        det_sector = "odd" in sector_name
        for index, vector in enumerate(vectors):
            for x1, x2 in root_pairs:
                value = two_root_leaf_average(
                    measure, x1, x2, basis, vector, det_sector
                )
                if value:
                    worst_two_root = max(worst_two_root, abs(value))
                    failures.append(
                        f"{sector_name}[{index}] at pair "
                        f"({x1},{x2}): {value}"
                    )
    dev = second_moment_deviation(measure)
    dev_zero = all(value == 0 for row in dev for value in row)
    if not dev_zero:
        failures.append(f"Dev != 0: {dev}")

    spin2_zero = True
    if spin2 is not None:
        basis, vectors = spin2
        flag_values = [
            spin2_flag_vector(measure, tuple(leaves), pair)
            for leaves, pair in basis
        ]
        for index, vector in enumerate(vectors):
            combo = [[Fraction(0)] * 3 for _ in range(3)]
            for coefficient, value in zip(vector, flag_values):
                if coefficient:
                    for r in range(3):
                        for c in range(3):
                            combo[r][c] += coefficient * value[r][c]
            if any(entry != 0 for row in combo for entry in row):
                spin2_zero = False
                failures.append(f"spin2 combo {index} != 0: {combo}")

    h2_value = (3 * pair_power(measure, 2) - 1) / 2
    if h2_value != 0:
        failures.append(f"h2 != 0: {h2_value}")
    g2_value = legendre_pair_energy(measure, 2)
    if g2_value != 0:
        failures.append(f"g_2 != 0: {g2_value}")

    report = {
        "measure": name,
        "atoms": len(measure),
        "all_zero": not failures,
        "energy_E": str(energy(measure)),
        "g4": str(legendre_pair_energy(measure, 4)),
        "g6": str(legendre_pair_energy(measure, 6)),
    }
    if require_zero and failures:
        report["failures"] = failures[:5]
    if not require_zero:
        report["violations_found"] = len(failures)
    return report


def negative_controls(measure: Measure) -> dict[str, str]:
    """Unweighted-admissible directions cut by the weighted projection
    must NOT vanish on general isotropic measures.

    Note: the mode-2 tower C^_n = -2 U_{n-2}(s)(t1^2 - 1/3)
    + 2 U_{n-1}(s)(t1 t2 - s/3) is y-quadratic and therefore remains
    weighted-admissible for every n (checked in solve_e1.py Part F); the
    weighted cut falls on the mode-6 tower S^_n and on T6 + 1/3, which
    are the controls used here.
    """

    # T6 + 1/3 leaf average at root x: sum_y w (T6(x.y) + 1/3)
    def t6_average(x: Vec) -> Fraction:
        total = Fraction(0)
        for y, w in measure:
            t = dot(x, y)
            total += w * (32 * t**6 - 48 * t**4 + 18 * t**2 - 1 + Fraction(1, 3))
        return total

    # S^_1 = U_5(t1) t2 - U_4(t1) s + s/3 (docs/SHARP_STRUCTURE.md §5)
    def s1_hat_average(x1: Vec, x2: Vec) -> Fraction:
        s_value = dot(x1, x2)
        total = Fraction(0)
        for y, w in measure:
            t1 = dot(x1, y)
            t2 = dot(x2, y)
            u5 = 32 * t1**5 - 32 * t1**3 + 6 * t1
            u4 = 16 * t1**4 - 12 * t1**2 + 1
            total += w * (u5 * t2 - u4 * s_value + s_value / 3)
        return total

    worst_t6 = max(
        (abs(t6_average(x)) for x in ROOT_POINTS), default=Fraction(0)
    )
    worst_s1 = max(
        (
            abs(s1_hat_average(ROOT_POINTS[i], ROOT_POINTS[j]))
            for i in range(len(ROOT_POINTS))
            for j in range(i + 1, len(ROOT_POINTS))
        ),
        default=Fraction(0),
    )
    return {
        "max |<T6 + 1/3>| over roots": str(worst_t6),
        "max |<S^_1>| over root pairs": str(worst_s1),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--projection", default="sdpa_runs/e1w_projection_deg14.json"
    )
    args = parser.parse_args()

    with open(args.projection) as handle:
        data = json.load(handle)
    assert data.get("weighted"), "audit expects a weighted projection file"
    sectors = {
        name: (
            [tuple(exponent) for exponent in entry["basis"]],
            [[Fraction(v) for v in vec] for vec in entry["vectors"]],
        )
        for name, entry in data["sectors"].items()
    }
    spin2 = None
    if data.get("spin2_flag"):
        entry = data["spin2_flag"]
        spin2 = (
            [(tuple(leaves), pair) for leaves, pair in entry["basis"]],
            [[Fraction(v) for v in vec] for vec in entry["vectors"]],
        )

    r1 = rot_z(Fraction(3, 5), Fraction(4, 5))
    r2 = rot_x(Fraction(5, 13), Fraction(12, 13))
    onb_r1 = rotate_measure(ONB, r1)
    onb_r2 = rotate_measure(ONB, r2)

    tests: list[tuple[str, Measure]] = [("ONB", ONB)]
    for t, u in (
        (Fraction(1, 2), Fraction(0)),
        (Fraction(1, 3), Fraction(1, 4)),
        (Fraction(2, 7), Fraction(3, 7)),
        (Fraction(1, 5), Fraction(1, 5)),
        (Fraction(0), Fraction(3, 5)),
    ):
        tests.append(
            (
                f"rotated-ONB mixture t={t}, u={u}",
                mix((ONB, 1 - t - u), (onb_r1, t), (onb_r2, u)),
            )
        )
    for angle, u in (
        ((Fraction(3, 5), Fraction(4, 5)), Fraction(1, 2)),
        ((Fraction(3, 5), Fraction(4, 5)), Fraction(1, 3)),
        ((Fraction(5, 13), Fraction(12, 13)), Fraction(1, 4)),
        ((Fraction(5, 13), Fraction(12, 13)), Fraction(2, 5)),
    ):
        tests.append(
            (
                f"pole-equator isotropic cos={angle[0]}, u={u}",
                pole_equator_isotropic(angle, u),
            )
        )
    tests.append(
        (
            "cross mixture",
            mix(
                (pole_equator_isotropic((Fraction(3, 5), Fraction(4, 5)), Fraction(1, 3)), Fraction(1, 2)),
                (mix((ONB, Fraction(1, 2)), (onb_r1, Fraction(1, 2))), Fraction(1, 2)),
            ),
        )
    )

    all_pass = True
    for name, measure in tests:
        report = audit_measure(name, measure, sectors, spin2)
        marker = "ok " if report["all_zero"] else "FAIL"
        print(
            f"  [{marker}] {name}: E = {report['energy_E']}, "
            f"g4 = {report['g4']}, g6 = {report['g6']}"
        )
        if not report["all_zero"]:
            all_pass = False
            for failure in report.get("failures", []):
                print(f"        {failure}")

    print("  negative controls (must be nonzero on some isotropic measure):")
    control = negative_controls(
        mix((ONB, Fraction(1, 2)), (onb_r1, Fraction(1, 2)))
    )
    for key, value in control.items():
        nonzero = Fraction(value) != 0
        marker = "ok " if nonzero else "FAIL"
        all_pass = all_pass and nonzero
        print(f"  [{marker}] {key} = {value}")

    non_isotropic = [
        ((Fraction(0), Fraction(0), Fraction(1)), Fraction(1, 2)),
        ((Fraction(1), Fraction(0), Fraction(0)), Fraction(1, 4)),
        ((Fraction(0), Fraction(1), Fraction(0)), Fraction(1, 4)),
    ]
    report = audit_measure(
        "non-isotropic control (pole weight 1/2)",
        non_isotropic,
        sectors,
        spin2,
        require_zero=False,
    )
    violated = report["violations_found"] > 0
    marker = "ok " if violated else "FAIL"
    all_pass = all_pass and violated
    print(
        f"  [{marker}] non-isotropic control violates "
        f"{report['violations_found']} conditions (audit has teeth)"
    )

    print()
    if all_pass:
        print("Weighted-(E1) isotropic audit: all checks passed.")
    else:
        raise SystemExit("Weighted-(E1) isotropic audit FAILED")


if __name__ == "__main__":
    main()
