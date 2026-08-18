#!/usr/bin/env python3
"""Theta-atom machinery: exact generators, label expansions, tail majorants.

Implements docs/ENRICHMENTS.md.  The theta atom for family f in {2, 6}
and rational q in (0, 1) is the resummed diagonal series

    tau_{f,q}(mu) = sum_{n in Z} q^{n^2} Q[Ghat^f_n](mu),

where Q[A] is the 1x1 two-root Gram block
E[(int A(x1.y, x2.y, x1.x2) dmu(y))^2] and Ghat^2_n = Chat_n,
Ghat^6_n = Shat_n are the modulated (E1)-admissible generators of
docs/SHARP_STRUCTURE.md.  Everything here is exact rational except
the explicitly numeric self-tests and ray pairings.

Library + CLI:
  --self-test          exact and numeric validation (run before trusting)
  --majorants          print exact tail majorants T^f_q(N)
  --expand-json PATH   dump exact label expansions of Q[Ghat_n]
  --pair RAY_JSON      pairing table of atom-induced cuts against a ray
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction

Poly3 = dict[tuple[int, int, int], Fraction]  # (i,j,k) -> coeff of t1^i t2^j s^k

A_CONST = Fraction(4, 3)
SLOPE = {2: Fraction(4), 6: Fraction(12)}


# ---------------------------------------------------------------------------
# Chebyshev helpers (exact)
# ---------------------------------------------------------------------------

def chebyshev_t(n: int) -> list[Fraction]:
    """Coefficient list of T_n, index = power."""
    n = abs(n)
    previous, current = [Fraction(1)], [Fraction(0), Fraction(1)]
    if n == 0:
        return previous
    for _ in range(n - 1):
        nxt = [Fraction(0)] + [2 * value for value in current]
        for power, value in enumerate(previous):
            nxt[power] -= value
        previous, current = current, nxt
    return current


def chebyshev_u(n: int) -> list[Fraction]:
    """Coefficient list of U_n for n >= -1 (U_{-1} = 0)."""
    if n == -1:
        return [Fraction(0)]
    previous, current = [Fraction(1)], [Fraction(0), Fraction(2)]
    if n == 0:
        return previous
    for _ in range(n - 1):
        nxt = [Fraction(0)] + [2 * value for value in current]
        for power, value in enumerate(previous):
            nxt[power] -= value
        previous, current = current, nxt
    return current


# ---------------------------------------------------------------------------
# Poly3 arithmetic
# ---------------------------------------------------------------------------

def poly3_add(a: Poly3, b: Poly3, scale: Fraction = Fraction(1)) -> Poly3:
    out = dict(a)
    for key, value in b.items():
        out[key] = out.get(key, Fraction(0)) + scale * value
        if out[key] == 0:
            del out[key]
    return out


def poly3_mul(a: Poly3, b: Poly3) -> Poly3:
    out: Poly3 = {}
    for (i, j, k), va in a.items():
        for (p, r, t), vb in b.items():
            key = (i + p, j + r, k + t)
            out[key] = out.get(key, Fraction(0)) + va * vb
            if out[key] == 0:
                del out[key]
    return out


def poly3_degree(a: Poly3) -> int:
    return max((sum(key) for key in a), default=0)


def poly3_eval(a: Poly3, t1: float, t2: float, s: float) -> float:
    return sum(
        float(value) * t1**i * t2**j * s**k for (i, j, k), value in a.items()
    )


def s_poly(coefficients: list[Fraction]) -> Poly3:
    return {
        (0, 0, power): value
        for power, value in enumerate(coefficients)
        if value != 0
    }


def t1_poly(coefficients: list[Fraction]) -> Poly3:
    return {
        (power, 0, 0): value
        for power, value in enumerate(coefficients)
        if value != 0
    }


# ---------------------------------------------------------------------------
# Modulated generator families (docs/SHARP_STRUCTURE.md)
# ---------------------------------------------------------------------------

C0: Poly3 = {(0, 0, 0): Fraction(-1), (2, 0, 0): Fraction(2)}       # T2(t1)
C1: Poly3 = {(1, 1, 0): Fraction(2), (0, 0, 1): Fraction(-1)}       # 2 t1 t2 - s
S0: Poly3 = t1_poly(chebyshev_t(6))                                  # T6(t1)
# S1 = U5(t1) t2 - U4(t1) s
S1: Poly3 = poly3_add(
    poly3_mul(t1_poly(chebyshev_u(5)), {(0, 1, 0): Fraction(1)}),
    poly3_mul(t1_poly(chebyshev_u(4)), {(0, 0, 1): Fraction(1)}),
    Fraction(-1),
)
TWO_S: Poly3 = {(0, 0, 1): Fraction(2)}


def raw_generator(family: int, n: int) -> Poly3:
    """C_n (family 2) or S_n (family 6), any integer n, via the recurrence."""
    seed0, seed1 = (C0, C1) if family == 2 else (S0, S1)
    if n == 0:
        return dict(seed0)
    if n > 0:
        previous, current = seed0, seed1
        for _ in range(n - 1):
            previous, current = current, poly3_add(
                poly3_mul(TWO_S, current), previous, Fraction(-1)
            )
        return current
    # backwards: G_{m-1} = 2 s G_m - G_{m+1}
    following, current = seed1, seed0
    for _ in range(-n):
        following, current = current, poly3_add(
            poly3_mul(TWO_S, current), following, Fraction(-1)
        )
    return current


def hat_generator(family: int, n: int) -> Poly3:
    """Ghat_n = G_n + T_|n|(s)/3, the (E1)-admissible modulated generator."""
    correction = {
        key: value / 3 for key, value in s_poly(chebyshev_t(abs(n))).items()
    }
    return poly3_add(raw_generator(family, n), correction)


def closed_form_c(n: int) -> Poly3:
    """C_n = C_0 T_n(s) + 2 t1 (t2 - s t1) U_{n-1}(s)  (n in Z).

    For n < 0 use T_{-n} = T_n and U_{-n-1} = -U_{n-1}.
    """
    sign = 1 if n >= 0 else -1
    u_index = abs(n) - 1
    term1 = poly3_mul(C0, s_poly(chebyshev_t(abs(n))))
    factor: Poly3 = {(1, 1, 0): Fraction(2), (2, 0, 1): Fraction(-2)}
    term2 = poly3_mul(factor, s_poly(chebyshev_u(u_index)))
    return poly3_add(term1, term2, Fraction(sign))


# ---------------------------------------------------------------------------
# Label expansion of the 1x1 two-root block Q[A]
# ---------------------------------------------------------------------------
# Vertex order (X=0, Z=1, Y=2, W=3); edge order of sos_search.graph_edges(4):
# (XZ, XY, XW, ZY, ZW, YW).  A is a polynomial in t1 = X.Y, t2 = Z.Y,
# s = X.Z; the squared block pairs leaf Y against an independent leaf W.

def expand_q_block(polynomial: Poly3):
    """Exact label expansion of Q[A] as dict[Label, Fraction]."""
    from sos_search import graph_expectation_label

    expansion: dict[tuple, Fraction] = {}
    terms = list(polynomial.items())
    for (i, j, k), left in terms:
        for (p, r, t), right in terms:
            label, reduction = graph_expectation_label(
                4, (k + t, i, p, j, r, 0)
            )
            if label is None or reduction == 0:
                continue
            value = left * right * reduction
            expansion[label] = expansion.get(label, Fraction(0)) + value
            if expansion[label] == 0:
                del expansion[label]
    return expansion


def h2_localize(expansion):
    """Exact label expansion of h2 * F given the expansion of F.

    h2 = (3 p2 - 1)/2, and p2 x (k-point label) is the disconnected
    product label; h2 in [0, 1] on probability measures, so h2-localized
    atoms keep the same tail majorants (docs/ENRICHMENTS.md sec. 8).
    """
    from sos_search import multiply_labels

    p2 = ("pair", 2)
    out: dict = {}
    for label, value in expansion.items():
        for factor, shifted in (
            (Fraction(3, 2), multiply_labels(p2, label)),
            (Fraction(-1, 2), label),
        ):
            updated = out.get(shifted, Fraction(0)) + factor * value
            if updated:
                out[shifted] = updated
            else:
                out.pop(shifted, None)
    return out


# ---------------------------------------------------------------------------
# Tail majorants (exact rational; docs/ENRICHMENTS.md section 2.1)
# ---------------------------------------------------------------------------

def tail_first_moment(family: int, q: Fraction, cut: int) -> Fraction:
    """Exact rational bound on sum_{n>N} q^{n^2} (4/3 + c n)."""
    c = SLOPE[family]
    alpha = A_CONST + c * cut
    r = q ** (2 * cut + 1)
    geom1 = r / (1 - r)
    geom2 = r / (1 - r) ** 2
    return q ** (cut * cut) * (alpha * geom1 + c * geom2)


def tail_second_moment(family: int, q: Fraction, cut: int) -> Fraction:
    """Exact rational bound on sum_{n>N} q^{n^2} (4/3 + c n)^2."""
    c = SLOPE[family]
    alpha = A_CONST + c * cut
    r = q ** (2 * cut + 1)
    geom1 = r / (1 - r)
    geom2 = r / (1 - r) ** 2
    geom3 = r * (1 + r) / (1 - r) ** 3
    return q ** (cut * cut) * (
        alpha * alpha * geom1 + 2 * alpha * c * geom2 + c * c * geom3
    )


def tail_majorant(family: int, q: Fraction, cut: int) -> Fraction:
    """T^f_q(N): exact bound on sum_{|n|>N} q^{n^2} Q[Ghat_n] at any measure."""
    return 2 * tail_second_moment(family, q, cut)


# ---------------------------------------------------------------------------
# Numeric evaluation on discrete antipodal measures (validation only)
# ---------------------------------------------------------------------------

def random_antipodal_measure(rng, atoms: int):
    import numpy as np

    points = rng.normal(size=(atoms, 3))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    weights = rng.random(atoms)
    weights /= weights.sum()
    full_points = np.vstack([points, -points])
    full_weights = np.concatenate([weights, weights]) / 2.0
    return full_points, full_weights


def label_moment(label, points, weights) -> float:
    """Evaluate a multigraph moment label on a discrete measure."""
    import numpy as np

    if label == ("constant",):
        return 1.0
    if label[0] == "product":
        value = 1.0
        for factor in label[1:]:
            value *= label_moment(factor, points, weights)
        return value
    if label[0] == "pair":
        vertex_count, exponents = 2, (label[1],)
    elif label[0] == "triangle":
        vertex_count, exponents = 3, tuple(label[1:])
    elif isinstance(label[0], str) and label[0].startswith("graph_"):
        vertex_count = int(label[0].split("_")[1])
        exponents = tuple(label[1:])
    else:
        raise ValueError(f"unsupported label {label}")
    gram = points @ points.T
    edges = [
        (left, right)
        for left in range(vertex_count)
        for right in range(left + 1, vertex_count)
    ]
    total = 0.0
    for combo in itertools.product(range(len(points)), repeat=vertex_count):
        weight = 1.0
        for vertex in combo:
            weight *= weights[vertex]
        term = weight
        for (left, right), exponent in zip(edges, exponents):
            if exponent:
                term *= gram[combo[left], combo[right]] ** exponent
        total += term
    return total


def q_block_direct(polynomial: Poly3, points, weights) -> float:
    """Q[A] evaluated directly on a discrete measure."""
    import numpy as np

    gram = points @ points.T
    total = 0.0
    for x1 in range(len(points)):
        for x2 in range(len(points)):
            leaf = sum(
                weights[y]
                * poly3_eval(
                    polynomial, gram[x1, y], gram[x2, y], gram[x1, x2]
                )
                for y in range(len(points))
            )
            total += weights[x1] * weights[x2] * leaf * leaf
    return total


# ---------------------------------------------------------------------------
# Ray pairing (C2)
# ---------------------------------------------------------------------------

def monomial_multipliers(maximum_degree: int) -> list[Poly3]:
    """All monomials t1^i t2^j s^k of total degree <= maximum_degree."""
    return [
        {(i, j, total - i - j): Fraction(1)}
        for total in range(maximum_degree + 1)
        for i in range(total + 1)
        for j in range(total - i + 1)
    ]


def localized_leaves(
    family: int, n: int, degree: int
) -> list[tuple[tuple[int, int, int], Poly3]]:
    """(multiplier exponent, leaf Ghat_n * m) pairs within the degree cap.

    Each leaf is bounded by (4/3 + c|n|) on the cube (|m| <= 1 there), so
    the per-multiplier atoms tau_{f,q,m} = sum_n q^{n^2} Q[Ghat_n m]
    satisfy the same tail majorants as the plain atom.
    Q[leaf] has degree 2 deg(leaf) <= degree.
    """
    base = hat_generator(family, n)
    room = degree // 2 - poly3_degree(base)
    if room < 0:
        return []
    return [
        (next(iter(m)), poly3_mul(base, m))
        for m in monomial_multipliers(room)
    ]


def family_window(family: int, degree: int) -> tuple[int, int]:
    """Range [n_min, n_max] with Q[Ghat_n] inside the degree-d label algebra.

    deg Ghat_n = n + off for n >= 0 and |n| + off + 1 for n < 0, with
    off = 1 (family 2) or 5 (family 6); Q doubles the degree.
    """
    offset = 1 if family == 2 else 5
    n_max = degree // 2 - offset
    n_min = -(degree // 2 - offset - 1)
    return n_min, n_max


def ray_pairings(ray: dict[str, float], degree: int = 14):
    """l_n(r) for each family and each within-degree n, plus coverage info."""
    rows = []
    for family in (2, 6):
        n_min, n_max = family_window(family, degree)
        for n in range(n_min, n_max + 1):
            expansion = expand_q_block(hat_generator(family, n))
            pairing = 0.0
            missing_weight = 0.0
            hit = 0
            for label, coefficient in expansion.items():
                key = str(label)
                if key in ray:
                    pairing += float(coefficient) * ray[key]
                    hit += 1
                else:
                    missing_weight += abs(float(coefficient))
            rows.append(
                {
                    "family": family,
                    "n": n,
                    "pairing": pairing,
                    "labels": len(expansion),
                    "labels_in_ray_support": hit,
                    "absent_label_coefficient_mass": missing_weight,
                }
            )
    return rows


def window_table(rows, q_values, degree: int = 14):
    """Pairings of the atom window cuts W_{N',N} against the ray.

    Constraint (W_{N',N}): T^f_q(N') - sum_{N'<|n|<=N} q^{n^2} l_n(y) >= 0.
    Ray pairing = -sum_window q^{n^2} l_n(r); negative kills the ray.
    """
    by_family: dict[int, dict[int, float]] = {2: {}, 6: {}}
    for row in rows:
        by_family[row["family"]][row["n"]] = row["pairing"]
    table = []
    for family in (2, 6):
        n_min, n_max = family_window(family, degree)
        cap = max(n_max, -n_min)
        for q in q_values:
            for low in range(0, cap):
                for high in range(low + 1, cap + 1):
                    weighted = sum(
                        float(q) ** (n * n) * by_family[family][n]
                        for n in range(max(-high, n_min), min(high, n_max) + 1)
                        if abs(n) > low
                    )
                    table.append(
                        {
                            "family": family,
                            "q": str(q),
                            "window": [low, high],
                            "cut_pairing": -weighted,
                            "majorant": float(tail_majorant(family, q, low)),
                        }
                    )
    return table


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------

def self_test() -> None:
    import numpy as np

    failures = []

    def check(name: str, condition: bool) -> None:
        print(("PASS " if condition else "FAIL ") + name)
        if not condition:
            failures.append(name)

    # 1. closed form == recurrence for C_n, n in [-8, 8]
    check(
        "closed form C_n = C_0 T_n + 2 t1 (t2 - s t1) U_{n-1} for n in [-8,8]",
        all(
            raw_generator(2, n) == closed_form_c(n) for n in range(-8, 9)
        ),
    )

    # 2. pole-value lemma: G_n(0, 0, s) = -T_|n|(s), hence Ghat_n(0,0,s) = -(2/3) T_|n|(s)
    def pole_slice(polynomial: Poly3) -> Poly3:
        return {
            key: value
            for key, value in polynomial.items()
            if key[0] == 0 and key[1] == 0
        }

    pole_ok = True
    for family in (2, 6):
        for n in range(-6, 7):
            expected = {
                (0, 0, power): -value
                for power, value in enumerate(chebyshev_t(abs(n)))
                if value != 0
            }
            if pole_slice(raw_generator(family, n)) != expected:
                pole_ok = False
    check("pole-value lemma G_n(0,0,s) = -T_|n|(s), both families", pole_ok)

    # 3. explicit expansions from the docs
    check(
        "C_2 and C_-1 match the documented expansions",
        raw_generator(2, 2)
        == {
            (1, 1, 1): Fraction(4),
            (0, 0, 2): Fraction(-2),
            (2, 0, 0): Fraction(-2),
            (0, 0, 0): Fraction(1),
        }
        and raw_generator(2, -1)
        == {
            (2, 0, 1): Fraction(4),
            (1, 1, 0): Fraction(-2),
            (0, 0, 1): Fraction(-1),
        },
    )
    check(
        "degrees: deg Ghat_n = n + off (n >= 1), |n| + off + 1 (n <= -1), "
        "off = 1 / 5",
        all(
            poly3_degree(hat_generator(2, n)) == n + 1
            and poly3_degree(hat_generator(6, n)) == n + 5
            for n in range(1, 7)
        )
        and all(
            poly3_degree(hat_generator(2, n)) == abs(n) + 2
            and poly3_degree(hat_generator(6, n)) == abs(n) + 6
            for n in range(-6, 0)
        ),
    )

    # 4. sup bounds |Ghat_n| <= 4/3 + c |n| on the cube (random + corner grid)
    rng = np.random.default_rng(7)
    grid = [-1.0, -0.5, 0.0, 0.5, 1.0]
    points = [
        (t1, t2, s)
        for t1 in grid
        for t2 in grid
        for s in grid
    ] + [tuple(rng.uniform(-1, 1, size=3)) for _ in range(4000)]
    bound_ok = True
    worst = 0.0
    for family in (2, 6):
        for n in range(-8, 9):
            polynomial = hat_generator(family, n)
            bound = float(A_CONST + SLOPE[family] * abs(n))
            for t1, t2, s in points:
                ratio = abs(poly3_eval(polynomial, t1, t2, s)) / bound
                worst = max(worst, ratio)
                if ratio > 1 + 1e-12:
                    bound_ok = False
    check(
        f"|Ghat_n| <= 4/3 + c|n| on the cube (worst ratio {worst:.6f})",
        bound_ok,
    )

    # 5. exact tail majorants dominate high-precision partial tails
    from decimal import Decimal, getcontext

    getcontext().prec = 60
    tails_ok = True
    for family in (2, 6):
        c = SLOPE[family]
        for q in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(9, 10)):
            for cut in (0, 1, 2, 4, 6):
                exact_first = Decimal(0)
                exact_second = Decimal(0)
                qd = Decimal(q.numerator) / Decimal(q.denominator)
                for n in range(cut + 1, cut + 400):
                    weight = qd ** (n * n)
                    linear = Decimal(4) / 3 + int(c) * n
                    exact_first += weight * linear
                    exact_second += weight * linear * linear
                first_bound = tail_first_moment(family, q, cut)
                second_bound = tail_second_moment(family, q, cut)
                for exact, bound in (
                    (exact_first, first_bound),
                    (exact_second, second_bound),
                ):
                    bound_dec = Decimal(bound.numerator) / Decimal(
                        bound.denominator
                    )
                    if bound_dec < exact:
                        tails_ok = False
    check("rational tail majorants dominate 400-term partial tails", tails_ok)

    # 6. label expansions of Q[Ghat_n] agree with direct numeric integration
    rng = np.random.default_rng(11)
    expansion_ok = True
    worst_gap = 0.0
    for family, orders in ((2, (-3, -1, 0, 1, 2, 4)), (6, (-1, 0, 2))):
        for n in orders:
            polynomial = hat_generator(family, n)
            expansion = expand_q_block(polynomial)
            for _ in range(3):
                points, weights = random_antipodal_measure(rng, 4)
                direct = q_block_direct(polynomial, points, weights)
                expanded = sum(
                    float(value) * label_moment(label, points, weights)
                    for label, value in expansion.items()
                )
                gap = abs(direct - expanded) / max(1.0, abs(direct))
                worst_gap = max(worst_gap, gap)
                if gap > 1e-9:
                    expansion_ok = False
                if direct < -1e-12:
                    expansion_ok = False
    check(
        f"Q[Ghat_n] label expansion == direct integration (worst rel gap {worst_gap:.2e})",
        expansion_ok,
    )

    # 7. sandwich validity on random measures:
    #    sum_{|n|<=N} q^{n^2} Q_n <= tau <= sum_{|n|<=N} + T(N)
    rng = np.random.default_rng(23)
    sandwich_ok = True
    for family in (2, 6):
        for q in (Fraction(1, 2), Fraction(3, 4)):
            points, weights = random_antipodal_measure(rng, 4)
            values = {
                n: q_block_direct(hat_generator(family, n), points, weights)
                for n in range(-12, 13)
            }
            tau = sum(float(q) ** (n * n) * values[n] for n in values)
            for cut in (0, 1, 2, 3):
                partial = sum(
                    float(q) ** (n * n) * values[n]
                    for n in range(-cut, cut + 1)
                )
                upper = partial + float(tail_majorant(family, q, cut))
                if not (partial - 1e-9 <= tau <= upper + 1e-9):
                    sandwich_ok = False
    check("sandwich cuts hold at random discrete measures", sandwich_ok)

    print()
    if failures:
        raise SystemExit(f"{len(failures)} self-test failure(s): {failures}")
    print("all theta-atom self-tests passed")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--majorants", action="store_true")
    parser.add_argument("--expand-json")
    parser.add_argument("--pair", metavar="RAY_JSON")
    parser.add_argument("--degree", type=int, default=14)
    parser.add_argument(
        "--q-values",
        default="1/4,1/2,3/4,9/10",
        help="comma list of rationals",
    )
    args = parser.parse_args()

    q_values = [Fraction(text) for text in args.q_values.split(",")]

    if args.self_test:
        self_test()
        return

    if args.majorants:
        for family in (2, 6):
            for q in q_values:
                for cut in range(0, 8):
                    value = tail_majorant(family, q, cut)
                    print(
                        f"T^({family})_{q}(N={cut}) = "
                        f"{value.numerator}/{value.denominator} "
                        f"~ {float(value):.6e}"
                    )
        return

    if args.expand_json:
        payload = []
        for family in (2, 6):
            cap = family_window(family, args.degree)
            for n in range(-cap, cap + 1):
                expansion = expand_q_block(hat_generator(family, n))
                payload.append(
                    {
                        "family": family,
                        "n": n,
                        "degree": 2 * poly3_degree(hat_generator(family, n)),
                        "expansion": {
                            str(label): f"{value.numerator}/{value.denominator}"
                            for label, value in sorted(
                                expansion.items(), key=lambda item: str(item[0])
                            )
                        },
                    }
                )
        with open(args.expand_json, "w") as handle:
            json.dump(payload, handle, indent=1)
        print(f"wrote {len(payload)} exact expansions to {args.expand_json}")
        return

    if args.pair:
        with open(args.pair) as handle:
            data = json.load(handle)
        ray = data["ray"]
        print(
            f"ray: status={data.get('status')} "
            f"squared_norm={data.get('squared_norm')} "
            f"labels={data.get('labels')} support={len(ray)}"
        )
        rows = ray_pairings(ray, args.degree)
        print("\nper-generator diagonal pairings l_n(r) = <Q[Ghat_n], ray>:")
        print(f"{'family':>6} {'n':>3} {'l_n(r)':>14} {'labels':>7} "
              f"{'in-support':>10} {'absent-mass':>12}")
        for row in rows:
            print(
                f"{row['family']:>6} {row['n']:>3} {row['pairing']:>14.6f} "
                f"{row['labels']:>7} {row['labels_in_ray_support']:>10} "
                f"{row['absent_label_coefficient_mass']:>12.3f}"
            )
        table = window_table(rows, q_values, args.degree)
        print("\natom window cuts W_{N',N}: pairing = -sum_{N'<|n|<=N} q^{n^2} l_n(r)")
        print("(negative pairing = valid cut violated along the ray = ray killed)")
        print(f"{'family':>6} {'q':>5} {'window':>9} {'cut_pairing':>14} "
              f"{'majorant T(N-)':>14}")
        for row in table:
            print(
                f"{row['family']:>6} {row['q']:>5} "
                f"({row['window'][0]},{row['window'][1]}]"
                f"{'':>2} {row['cut_pairing']:>14.6f} "
                f"{row['majorant']:>14.6e}"
            )
        out = {
            "ray_file": args.pair,
            "ray_status": data.get("status"),
            "per_generator": rows,
            "window_cuts": table,
        }
        out_path = args.pair.replace(".json", "") + "_theta_pairings.json"
        with open(out_path, "w") as handle:
            json.dump(out, handle, indent=1)
        print(f"\nwrote {out_path}")
        return

    parser.print_help()


if __name__ == "__main__":
    main()
