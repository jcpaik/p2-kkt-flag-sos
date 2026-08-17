"""Exact solution of the (E1) circle-mode equations for certificate leaves.

Complementary slackness against the zero-energy family

    mu* = (1/3) delta_{+-e} + (2/3) nu,   nu on C = e^perp,
    nu antipodal, nu-hat(2) = nu-hat(6) = 0,

forces every flag-square leaf Phi(R; .) of a sharp certificate to satisfy

    int Phi(R; y) dmu*(y) = 0   for all roots R in (supp mu*)^r
                                and all admissible (e, nu).       (E1)

This script solves (E1) in closed form, in exact rational arithmetic:

  Part A: one-root leaves, every O(2) spin sector.  The admissible spaces
          are finite dimensional and degree independent:

            spin 0:  span{ T2 + 1/3,  T6 + 1/3 }          (dim 2)
            spin 1:  span{ t,  U5 }                       (dim 2)
            spin 2:  span{ 1,  (4t^2-1)^2 }  (radial)     (dim 2)
            spin 3:  span{ t^3 }                          (dim 1)
            spin >= 4:  {0}
          and K itself = (T2 + 1/3) + (T6 + 1/3).

  Part B: unrooted (arity-0) families -- harmonic energies, the pair-power
          flags, spin-2 deviatoric flags, and the harmonic-weighted pair
          flags of order 2/4/6.

  Part C: two-root leaves A(t1,t2,s) + det(x1,x2,y) B(t1,t2,s): exact
          admissible dimensions per degree, and closed-form generator
          families (the infinite-dimensional part of the certificate).

Method.  Fix the standard zero measure: poles +-e = (0,0,+-1), equator
C = {y(th) = (cos th, sin th, 0)}.  Every root configuration is, up to a
rotation about e and reflections, a tuple of points from {+-e} u C with
symbolic equator angles.  The unknown leaf integral against mu* is

    (1/6)[Phi(y=+e) + Phi(y=-e)] + (2/3) sum_m  c_m tau_m,

where c_m are the Fourier coefficients of Phi restricted to C and
tau_m = int e^{i m th} dnu.  Admissibility constraints on nu give
tau_0 = 1, tau_m = 0 for m odd and m in {+-2, +-6}; all remaining tau_m
are free complex symbols.  Since the admissible truncated moment sets
have nonempty interior around Haar measure, "= 0 for every admissible nu"
is equivalent to "= 0 identically as a polynomial in the free tau_m"
(and, for configurations with a symbolic root angle, identically in that
angle).  All computations below are Laurent-polynomial manipulations over
Gaussian rationals; no floating point is used anywhere.

Run with

    python3 solve_e1.py            # full check suite
    python3 solve_e1.py --table    # also print arity-2 dimension tables
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product as iter_product

# ----------------------------------------------------------------------------
# Gaussian-rational Laurent polynomials in several angle variables.
#
# A value is dict[key -> (re, im)] with key a tuple of integer Fourier modes,
# one slot per angle variable, and (re, im) a pair of Fractions.
# ----------------------------------------------------------------------------

CF = tuple[Fraction, Fraction]
Key = tuple[int, ...]
Laurent = dict[Key, CF]

CF_ZERO: CF = (Fraction(0), Fraction(0))


def cf(re=0, im=0) -> CF:
    return (Fraction(re), Fraction(im))


def cf_add(x: CF, y: CF) -> CF:
    return (x[0] + y[0], x[1] + y[1])


def cf_mul(x: CF, y: CF) -> CF:
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def cf_scale(x: CF, r: Fraction) -> CF:
    return (x[0] * r, x[1] * r)


def cf_is_zero(x: CF) -> bool:
    return x[0] == 0 and x[1] == 0


def l_const(slots: int, value: CF) -> Laurent:
    if cf_is_zero(value):
        return {}
    return {(0,) * slots: value}


def l_one(slots: int) -> Laurent:
    return l_const(slots, cf(1))


def l_add(p: Laurent, q: Laurent) -> Laurent:
    result = dict(p)
    for key, value in q.items():
        total = cf_add(result.get(key, CF_ZERO), value)
        if cf_is_zero(total):
            result.pop(key, None)
        else:
            result[key] = total
    return result


def l_scale(p: Laurent, r: Fraction) -> Laurent:
    if r == 0:
        return {}
    return {key: cf_scale(value, r) for key, value in p.items()}


def l_mul(p: Laurent, q: Laurent) -> Laurent:
    result: Laurent = {}
    for key_p, value_p in p.items():
        for key_q, value_q in q.items():
            key = tuple(a + b for a, b in zip(key_p, key_q))
            total = cf_add(result.get(key, CF_ZERO), cf_mul(value_p, value_q))
            if cf_is_zero(total):
                result.pop(key, None)
            else:
                result[key] = total
    return result


def l_pow(p: Laurent, exponent: int, slots: int) -> Laurent:
    result = l_one(slots)
    for _ in range(exponent):
        result = l_mul(result, p)
    return result


def l_cos(slots: int, slot: int, harmonic: int = 1) -> Laurent:
    key_plus = tuple(harmonic if i == slot else 0 for i in range(slots))
    key_minus = tuple(-harmonic if i == slot else 0 for i in range(slots))
    if harmonic == 0:
        return l_one(slots)
    return {key_plus: cf(Fraction(1, 2)), key_minus: cf(Fraction(1, 2))}


def l_sin(slots: int, slot: int, harmonic: int = 1) -> Laurent:
    if harmonic == 0:
        return {}
    key_plus = tuple(harmonic if i == slot else 0 for i in range(slots))
    key_minus = tuple(-harmonic if i == slot else 0 for i in range(slots))
    return {
        key_plus: cf(0, Fraction(-1, 2)),
        key_minus: cf(0, Fraction(1, 2)),
    }


Vector = list[Laurent]


def v_equator(slots: int, slot: int) -> Vector:
    return [l_cos(slots, slot), l_sin(slots, slot), {}]


def v_pole(slots: int, sign: int) -> Vector:
    return [{}, {}, l_const(slots, cf(sign))]


def l_i_times(p: Laurent) -> Laurent:
    """Multiply a Laurent polynomial by the imaginary unit."""

    return {key: (-value[1], value[0]) for key, value in p.items()}


def l_exp(slots: int, slot: int, harmonic: int = 1) -> Laurent:
    """e^{i harmonic theta_slot}."""

    return l_add(l_cos(slots, slot, harmonic), l_i_times(l_sin(slots, slot, harmonic)))


def v_dot(u: Vector, v: Vector) -> Laurent:
    result: Laurent = {}
    for cu, cv in zip(u, v):
        result = l_add(result, l_mul(cu, cv))
    return result


def v_det(u: Vector, v: Vector, w: Vector) -> Laurent:
    result: Laurent = {}
    for (i, j, k, sign) in (
        (0, 1, 2, 1), (1, 2, 0, 1), (2, 0, 1, 1),
        (0, 2, 1, -1), (2, 1, 0, -1), (1, 0, 2, -1),
    ):
        term = l_mul(l_mul(u[i], v[j]), w[k])
        result = l_add(result, l_scale(term, Fraction(sign)))
    return result


# ----------------------------------------------------------------------------
# Pairing an angle slot against nu.
#
# tau-monomials are canonical sorted tuples of the signed free modes.  tau_0
# contributes 1; constrained modes contribute 0.  `constrained` lists |m| that
# vanish (default {2, 6} plus all odd m).
# ----------------------------------------------------------------------------

TauMonomial = tuple[int, ...]


def pair_slots(
    poly: Laurent,
    circle_slots: list[int],
    constrained: frozenset[int] = frozenset({2, 6}),
) -> dict[tuple[Key, TauMonomial], CF]:
    """Integrate the given slots against independent copies of nu.

    Returns dict keyed by (remaining-slot key, tau monomial).  Slots not in
    `circle_slots` survive in the remaining key (symbolic root angles).
    """

    result: dict[tuple[Key, TauMonomial], CF] = {}
    for key, value in poly.items():
        taus: list[int] = []
        dead = False
        for slot in circle_slots:
            mode = key[slot]
            if mode == 0:
                continue
            if mode % 2 != 0 or abs(mode) in constrained:
                dead = True
                break
            taus.append(mode)
        if dead:
            continue
        remaining = tuple(
            m for i, m in enumerate(key) if i not in circle_slots
        )
        tau_key = tuple(sorted(taus))
        total = cf_add(result.get((remaining, tau_key), CF_ZERO), value)
        if cf_is_zero(total):
            result.pop((remaining, tau_key), None)
        else:
            result[(remaining, tau_key)] = total
    return result


def leaf_integral(
    leaf_on_circle: Laurent,
    leaf_at_pole_plus: Laurent,
    leaf_at_pole_minus: Laurent,
    theta_slot: int,
    constrained: frozenset[int] = frozenset({2, 6}),
) -> dict[tuple[Key, TauMonomial], CF]:
    """int Phi dmu* for one leaf slot: (1/6)(poles) + (2/3) circle-vs-nu.

    All inputs share the same slot layout; the pole values must not depend on
    the theta slot.  The output key drops the theta slot.
    """

    pole = l_scale(l_add(leaf_at_pole_plus, leaf_at_pole_minus), Fraction(1, 6))
    circle = l_scale(leaf_on_circle, Fraction(2, 3))
    paired = pair_slots(circle, [theta_slot], constrained)
    for key, value in pole.items():
        remaining = tuple(m for i, m in enumerate(key) if i != theta_slot)
        assert key[theta_slot] == 0
        entry = (remaining, ())
        total = cf_add(paired.get(entry, CF_ZERO), value)
        if cf_is_zero(total):
            paired.pop(entry, None)
        else:
            paired[entry] = total
    return paired


# ----------------------------------------------------------------------------
# Exact linear algebra: sparse rational RREF, nullspace, membership.
# ----------------------------------------------------------------------------

Row = dict[int, Fraction]


def add_condition_rows(
    rows: dict[object, Row],
    conditions: dict[tuple[Key, TauMonomial], CF],
    tag: object,
    column: int,
) -> None:
    """Turn a condition dictionary into linear rows on basis column `column`."""

    for (key, tau_key), value in conditions.items():
        for part, component in (("re", value[0]), ("im", value[1])):
            if component == 0:
                continue
            row = rows.setdefault((tag, key, tau_key, part), {})
            row[column] = row.get(column, Fraction(0)) + component


def nullspace(rows: list[Row], width: int) -> list[list[Fraction]]:
    """Nullspace basis of the sparse rational row system."""

    reduced: list[tuple[int, Row]] = []  # (pivot column, row)
    for row in rows:
        work = dict(row)
        for pivot, pivot_row in reduced:
            coefficient = work.get(pivot)
            if coefficient is None:
                continue
            for col, val in pivot_row.items():
                updated = work.get(col, Fraction(0)) - coefficient * val
                if updated == 0:
                    work.pop(col, None)
                else:
                    work[col] = updated
        if not work:
            continue
        pivot = min(work)
        inverse = 1 / work[pivot]
        work = {col: val * inverse for col, val in work.items()}
        reduced.append((pivot, work))
    # back-substitute to full RREF
    reduced.sort()
    for index in range(len(reduced) - 1, -1, -1):
        pivot, row = reduced[index]
        for j in range(index):
            other_pivot, other = reduced[j]
            coefficient = other.get(pivot)
            if coefficient is None:
                continue
            for col, val in row.items():
                updated = other.get(col, Fraction(0)) - coefficient * val
                if updated == 0:
                    other.pop(col, None)
                else:
                    other[col] = updated
    pivots = {pivot for pivot, _ in reduced}
    free_columns = [c for c in range(width) if c not in pivots]
    basis = []
    for free in free_columns:
        vector = [Fraction(0)] * width
        vector[free] = Fraction(1)
        for pivot, row in reduced:
            vector[pivot] = -row.get(free, Fraction(0))
        basis.append(vector)
    return basis


def in_span(vector: list[Fraction], basis: list[list[Fraction]]) -> bool:
    """Exact membership of `vector` in the span of `basis`."""

    rows = [dict((i, v) for i, v in enumerate(b) if v != 0) for b in basis]
    work = dict((i, v) for i, v in enumerate(vector) if v != 0)
    reduced: list[tuple[int, Row]] = []
    for row in rows:
        current = dict(row)
        for pivot, pivot_row in reduced:
            coefficient = current.get(pivot)
            if coefficient is None:
                continue
            for col, val in pivot_row.items():
                updated = current.get(col, Fraction(0)) - coefficient * val
                if updated == 0:
                    current.pop(col, None)
                else:
                    current[col] = updated
        if not current:
            continue
        pivot = min(current)
        inverse = 1 / current[pivot]
        reduced.append((pivot, {c: v * inverse for c, v in current.items()}))
    for pivot, pivot_row in reduced:
        coefficient = work.get(pivot)
        if coefficient is None:
            continue
        for col, val in pivot_row.items():
            updated = work.get(col, Fraction(0)) - coefficient * val
            if updated == 0:
                work.pop(col, None)
            else:
                work[col] = updated
    return not work


def same_span(a: list[list[Fraction]], b: list[list[Fraction]]) -> bool:
    return (
        len(a) == len(b)
        and all(in_span(v, b) for v in a)
        and all(in_span(v, a) for v in b)
    )


# ----------------------------------------------------------------------------
# Part A: one-root leaves.
#
# A spin-k leaf with radial monomial t^d is Phi(x; y) = t^d w^k, with
# t = x.y and w = y.u + i y.v in a right-handed frame (u, v, x).
# ----------------------------------------------------------------------------


def arity1_rows(
    spin: int,
    degrees: list[int],
    include_tilt_check: bool = False,
) -> list[Row]:
    rows: dict[object, Row] = {}
    for column, degree in enumerate(degrees):
        # Config A: root at the pole, frame u=(1,0,0), v=(0,1,0).
        # On C: t = 0 (so t^d = [d == 0]) and w = e^{i theta};
        # at y = +-x: t = +-1, w = 0.
        slots = 1
        w_circle = l_exp(slots, 0)
        if degree == 0:
            circle = l_pow(w_circle, spin, slots)
        else:
            circle = {}
        if spin == 0:
            pole_plus = l_const(slots, cf(1))          # (+1)^d * w^0
            pole_minus = l_const(slots, cf((-1) ** degree))
        else:
            pole_plus = {}
            pole_minus = {}
        conditions = leaf_integral(circle, pole_plus, pole_minus, 0)
        add_condition_rows(rows, conditions, "poleroot", column)

        # Config B: root x = (1,0,0) on the equator, frame u=(0,1,0), v=e.
        # On C: t = cos theta, w = sin theta; at y = +-e: t = 0, w = +-i.
        t_circle = l_cos(slots, 0)
        w_circle = l_sin(slots, 0)
        circle = l_mul(l_pow(t_circle, degree, slots), l_pow(w_circle, spin, slots))
        i_unit = l_const(slots, cf(0, 1))
        minus_i_unit = l_const(slots, cf(0, -1))
        zero_deg = l_one(slots) if degree == 0 else {}
        pole_plus = l_mul(zero_deg, l_pow(i_unit, spin, slots))
        pole_minus = l_mul(zero_deg, l_pow(minus_i_unit, spin, slots))
        conditions = leaf_integral(circle, pole_plus, pole_minus, 0)
        add_condition_rows(rows, conditions, "equatorroot", column)

        if include_tilt_check:
            # Config B tilted by a symbolic angle beta (slot 1): the circle
            # through x with poles +-(0, -sin beta, cos beta).  Completeness of
            # configs A and B is equivalent to these rows adding nothing.
            slots = 2
            t_circle = l_cos(slots, 0)
            w_circle = l_mul(l_sin(slots, 0), l_exp(slots, 1))  # sin theta e^{i beta}
            circle = l_mul(
                l_pow(t_circle, degree, slots), l_pow(w_circle, spin, slots)
            )
            e_beta = l_i_times(l_exp(slots, 1))  # w(+e_beta) = i e^{i beta}
            minus_e_beta = l_scale(e_beta, Fraction(-1))
            zero_deg = l_one(slots) if degree == 0 else {}
            pole_plus = l_mul(zero_deg, l_pow(e_beta, spin, slots))
            pole_minus = l_mul(zero_deg, l_pow(minus_e_beta, spin, slots))
            conditions = leaf_integral(circle, pole_plus, pole_minus, 0)
            add_condition_rows(rows, conditions, "tiltedroot", column)

    return list(rows.values())


def chebyshev_t(n: int) -> list[Fraction]:
    """Coefficient list (by power of t) of the Chebyshev polynomial T_n."""

    previous, current = [Fraction(1)], [Fraction(0), Fraction(1)]
    if n == 0:
        return previous
    for _ in range(n - 1):
        doubled = [Fraction(0)] + [2 * c for c in current]
        nxt = [
            doubled[i] - (previous[i] if i < len(previous) else Fraction(0))
            for i in range(len(doubled))
        ]
        previous, current = current, nxt
    return current


def chebyshev_u(n: int) -> list[Fraction]:
    previous, current = [Fraction(1)], [Fraction(0), Fraction(2)]
    if n == 0:
        return previous
    for _ in range(n - 1):
        doubled = [Fraction(0)] + [2 * c for c in current]
        nxt = [
            doubled[i] - (previous[i] if i < len(previous) else Fraction(0))
            for i in range(len(doubled))
        ]
        previous, current = current, nxt
    return current


def coefficients_to_vector(
    coefficients: list[Fraction], degrees: list[int]
) -> list[Fraction]:
    vector = [Fraction(0)] * len(degrees)
    for power, value in enumerate(coefficients):
        if value == 0:
            continue
        vector[degrees.index(power)] = value
    return vector


def poly_add(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    length = max(len(a), len(b))
    return [
        (a[i] if i < len(a) else Fraction(0))
        + (b[i] if i < len(b) else Fraction(0))
        for i in range(length)
    ]


def poly_scale(a: list[Fraction], r) -> list[Fraction]:
    return [Fraction(r) * c for c in a]


def poly_mul(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] += ai * bj
    return result


# ----------------------------------------------------------------------------
# Part C: two-root leaves.
# ----------------------------------------------------------------------------


def two_root_basis(max_degree: int, leaf_parity: int) -> list[tuple[int, int, int]]:
    """Monomials t1^i t2^j s^k with i+j+k <= max_degree, i+j = leaf_parity mod 2."""

    basis = []
    for i in range(max_degree + 1):
        for j in range(max_degree + 1 - i):
            for k in range(max_degree + 1 - i - j):
                if (i + j) % 2 == leaf_parity:
                    basis.append((i, j, k))
    return basis


def arity2_rows(
    basis: list[tuple[int, int, int]],
    det_sector: bool,
    configs: tuple[str, ...] = ("b2", "a2", "a2m", "c2p", "c2m"),
) -> list[Row]:
    rows: dict[object, Row] = {}
    for column, (i, j, k) in enumerate(basis):
        for name in configs:
            if name == "b2":
                slots = 2  # theta (slot 0), delta (slot 1)
                x1 = [l_one(slots), {}, {}]
                x2 = [l_cos(slots, 1), l_sin(slots, 1), {}]
            elif name == "a2":
                slots = 1
                x1 = v_pole(slots, 1)
                x2 = [l_one(slots), {}, {}]
            elif name == "a2m":
                slots = 1
                x1 = [l_one(slots), {}, {}]
                x2 = v_pole(slots, 1)
            elif name == "c2p":
                slots = 1
                x1 = v_pole(slots, 1)
                x2 = v_pole(slots, 1)
            elif name == "c2m":
                slots = 1
                x1 = v_pole(slots, 1)
                x2 = v_pole(slots, -1)
            else:
                raise ValueError(name)

            y_circle = v_equator(slots, 0)
            y_plus = v_pole(slots, 1)
            y_minus = v_pole(slots, -1)

            s_val = v_dot(x1, x2)
            s_power = l_pow(s_val, k, slots)

            def leaf_value(y: Vector) -> Laurent:
                t1 = v_dot(x1, y)
                t2 = v_dot(x2, y)
                value = l_mul(
                    l_mul(l_pow(t1, i, slots), l_pow(t2, j, slots)), s_power
                )
                if det_sector:
                    value = l_mul(value, v_det(x1, x2, y))
                return value

            conditions = leaf_integral(
                leaf_value(y_circle), leaf_value(y_plus), leaf_value(y_minus), 0
            )
            add_condition_rows(rows, conditions, name, column)
    return list(rows.values())


def evaluate_rows_on_vector(rows: list[Row], vector: list[Fraction]) -> bool:
    """True iff the vector satisfies every condition row."""

    for row in rows:
        total = Fraction(0)
        for col, val in row.items():
            total += val * vector[col]
        if total != 0:
            return False
    return True


def two_root_vector(
    poly: dict[tuple[int, int, int], Fraction],
    basis: list[tuple[int, int, int]],
) -> list[Fraction]:
    vector = [Fraction(0)] * len(basis)
    for monomial, value in poly.items():
        vector[basis.index(monomial)] = value
    return vector


def poly3_mul(
    a: dict[tuple[int, int, int], Fraction],
    b: dict[tuple[int, int, int], Fraction],
) -> dict[tuple[int, int, int], Fraction]:
    result: dict[tuple[int, int, int], Fraction] = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = (ka[0] + kb[0], ka[1] + kb[1], ka[2] + kb[2])
            result[key] = result.get(key, Fraction(0)) + va * vb
            if result[key] == 0:
                del result[key]
    return result


GRAM3 = {
    (0, 0, 0): Fraction(1),
    (1, 1, 1): Fraction(2),
    (2, 0, 0): Fraction(-1),
    (0, 2, 0): Fraction(-1),
    (0, 0, 2): Fraction(-1),
}


# ----------------------------------------------------------------------------
# Part B: unrooted families.
# ----------------------------------------------------------------------------


def multi_slot_integral(
    integrand_per_combo,
    slot_count: int,
    constrained: frozenset[int] = frozenset({2, 6}),
) -> dict[TauMonomial, CF]:
    """Integrate a fully-unrooted multilinear functional against mu*^slots.

    `integrand_per_combo(points, slots)` receives one vector per slot (pole or
    equator with its own angle slot) and returns a Laurent polynomial in the
    active angle slots; the function sums the 3^slot_count combinations with
    the correct weights and pairs every active angle against nu.
    """

    totals: dict[TauMonomial, CF] = {}
    choices = [("+", "-", "c")] * slot_count
    for combo in iter_product(*choices):
        weight = Fraction(1)
        points: Vector = []
        active = [index for index, c in enumerate(combo) if c == "c"]
        for index, choice in enumerate(combo):
            if choice == "+":
                weight *= Fraction(1, 6)
                points.append(v_pole(slot_count, 1))
            elif choice == "-":
                weight *= Fraction(1, 6)
                points.append(v_pole(slot_count, -1))
            else:
                weight *= Fraction(2, 3)
                points.append(v_equator(slot_count, index))
        integrand = integrand_per_combo(points, slot_count)
        integrand = l_scale(integrand, weight)
        paired = pair_slots(integrand, active, constrained)
        for (remaining, tau_key), value in paired.items():
            assert all(m == 0 for m in remaining)
            total = cf_add(totals.get(tau_key, CF_ZERO), value)
            if cf_is_zero(total):
                totals.pop(tau_key, None)
            else:
                totals[tau_key] = total
    return totals


def legendre_coefficients(order: int) -> list[Fraction]:
    """Coefficients of the Legendre polynomial P_order by power."""

    import math

    coefficients = [Fraction(0)] * (order + 1)
    for m in range(order // 2 + 1):
        power = order - 2 * m
        value = (
            Fraction((-1) ** m)
            * Fraction(math.comb(order, m))
            * Fraction(math.comb(2 * order - 2 * m, order))
            / Fraction(2**order)
        )
        coefficients[power] = value
    return coefficients


def pair_energy(order: int) -> dict[TauMonomial, CF]:
    """g_order(mu*) = int int P_order(x.y) dmu* dmu* as a tau polynomial."""

    coefficients = legendre_coefficients(order)

    def integrand(points, slots):
        inner = v_dot(points[0], points[1])
        total: Laurent = {}
        for power, value in enumerate(coefficients):
            if value == 0:
                continue
            total = l_add(total, l_scale(l_pow(inner, power, slots), value))
        return total

    return multi_slot_integral(integrand, 2)


def pair_power_moment(power: int) -> dict[TauMonomial, CF]:
    def integrand(points, slots):
        return l_pow(v_dot(points[0], points[1]), power, slots)

    return multi_slot_integral(integrand, 2)


def spin2_flag_value(flag: tuple[tuple[int, ...], int]) -> list[dict[TauMonomial, CF]]:
    """The 6 independent entries of V = int (x x^T - I/3) (leaf factors) dmu*^n."""

    leaves, pair_power = flag
    slot_count = 1 + len(leaves)
    entries = []
    for (r, c) in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)):
        def integrand(points, slots, r=r, c=c):
            x = points[0]
            value = l_mul(x[r], x[c])
            if r == c:
                value = l_add(value, l_const(slots, cf(Fraction(-1, 3))))
            for index, leaf_degree in enumerate(leaves):
                value = l_mul(
                    value,
                    l_pow(v_dot(x, points[1 + index]), leaf_degree, slots),
                )
            if pair_power:
                value = l_mul(
                    value,
                    l_pow(v_dot(points[1], points[2]), pair_power, slots),
                )
            return value

        entries.append(multi_slot_integral(integrand, slot_count))
    return entries


def harmonic_pair_flag_matrix(
    order: int, weights: list[int]
) -> dict[TauMonomial, list[list[Fraction]]]:
    """Coefficient matrices of G_ab = E[P_l(X.Z)(X.Y)^a(Z.W)^b] at mu*."""

    coefficients = legendre_coefficients(order)
    size = len(weights)
    matrices: dict[TauMonomial, list[list[Fraction]]] = {}
    for row, a in enumerate(weights):
        for col, b in enumerate(weights):
            def integrand(points, slots, a=a, b=b):
                x, z, y, w = points
                inner = v_dot(x, z)
                legendre: Laurent = {}
                for power, value in enumerate(coefficients):
                    if value == 0:
                        continue
                    legendre = l_add(
                        legendre, l_scale(l_pow(inner, power, slots), value)
                    )
                return l_mul(
                    legendre,
                    l_mul(
                        l_pow(v_dot(x, y), a, slots),
                        l_pow(v_dot(z, w), b, slots),
                    ),
                )

            totals = multi_slot_integral(integrand, 4)
            for tau_key, value in totals.items():
                assert value[1] == 0
                matrix = matrices.setdefault(
                    tau_key,
                    [[Fraction(0)] * size for _ in range(size)],
                )
                matrix[row][col] += value[0]
    return matrices


# ----------------------------------------------------------------------------
# Projection export for sos_search.py --e1-project.
# ----------------------------------------------------------------------------

SECTOR_PARITIES = {
    "two_root_even_00": (0, 0, 0),
    "two_root_even_11": (0, 1, 1),
    "two_root_odd_01": (1, 0, 1),
    "two_root_odd_10": (1, 1, 0),
}


def solver_monomials(max_degree: int) -> list[tuple[int, int, int]]:
    """Replicates sos_search.monomials ordering exactly."""

    if max_degree < 0:
        return []
    return [
        (i, j, total - i - j)
        for total in range(max_degree + 1)
        for i in range(total + 1)
        for j in range(total - i + 1)
    ]


def solver_spin2_basis(cap: int) -> list[tuple[tuple[int, ...], int]]:
    """Replicates sos_search.spin2_flag_basis exactly."""

    flags: list[tuple[tuple[int, ...], int]] = [((), 0)]
    for a in range(2, cap + 1, 2):
        flags.append(((a,), 0))
    for a1 in range(1, cap + 1):
        for a2 in range(a1, cap + 1):
            if (a1 + a2) % 2:
                continue
            for b in range(cap + 1):
                if (a1 + b) % 2 or (a2 + b) % 2:
                    continue
                if a1 + a2 + b > cap:
                    continue
                flags.append(((a1, a2), b))
    return flags


def export_projection(path: str, solver_degree: int) -> None:
    import json

    data: dict = {"solver_degree": solver_degree, "sectors": {}}
    flag_degree = solver_degree // 2
    localizing_degree = (solver_degree - 2) // 2
    for name, (p_ij, p_ik, p_jk) in SECTOR_PARITIES.items():
        det_sector = "odd" in name
        for minor in (False, True):
            degree = localizing_degree if minor else flag_degree
            basis = [
                (i, j, k)
                for (i, j, k) in solver_monomials(degree)
                if ((i + j) % 2, (i + k) % 2, (j + k) % 2)
                == (p_ij, p_ik, p_jk)
            ]
            # The minor sectors carry the root weight 1 - s^2, which vanishes
            # on the s = +-1 configurations, so those impose no condition.
            configs = (
                ("b2", "a2", "a2m") if minor else ("b2", "a2", "a2m", "c2p", "c2m")
            )
            rows = arity2_rows(basis, det_sector=det_sector, configs=configs)
            space = nullspace(rows, len(basis))
            key = name + ("_minor" if minor else "")
            data["sectors"][key] = {
                "basis": [list(exponent) for exponent in basis],
                "vectors": [[str(value) for value in vec] for vec in space],
            }
            print(
                f"  {key}: {len(space)} admissible of {len(basis)}"
            )
    spin2_cap = (solver_degree - 2) // 2
    spin2_basis = solver_spin2_basis(spin2_cap)
    rows_map: dict[object, Row] = {}
    for column, flag in enumerate(spin2_basis):
        entries = spin2_flag_value(flag)
        for entry_index, totals in enumerate(entries):
            add_condition_rows(
                rows_map,
                {((), key): value for key, value in totals.items()},
                ("spin2", entry_index),
                column,
            )
    spin2_space = nullspace(list(rows_map.values()), len(spin2_basis))
    data["spin2_flag"] = {
        "basis": [[list(leaves), pair] for (leaves, pair) in spin2_basis],
        "vectors": [[str(value) for value in vec] for vec in spin2_space],
    }
    print(f"  spin2_flag: {len(spin2_space)} admissible of {len(spin2_basis)}")
    with open(path, "w") as handle:
        json.dump(data, handle, indent=1)
    print(f"  wrote {path}")


# ----------------------------------------------------------------------------
# Check harness.
# ----------------------------------------------------------------------------

checks: list[tuple[str, bool]] = []


def check(name: str, condition: bool) -> None:
    checks.append((name, bool(condition)))
    marker = "ok " if condition else "FAIL"
    print(f"  [{marker}] {name}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", action="store_true")
    parser.add_argument("--max-table-degree", type=int, default=10)
    parser.add_argument(
        "--export-projection",
        metavar="PATH",
        help="write the (E1)-admissible bases for sos_search.py --e1-project",
    )
    parser.add_argument("--solver-degree", type=int, default=14)
    args = parser.parse_args()

    if args.export_projection:
        print(f"Exporting (E1) projection at solver degree {args.solver_degree}")
        export_projection(args.export_projection, args.solver_degree)
        return

    print("Part A: one-root (E1)-admissible leaves, by spin sector")

    predicted = {
        0: [
            poly_add(chebyshev_t(2), [Fraction(1, 3)]),
            poly_add(chebyshev_t(6), [Fraction(1, 3)]),
        ],
        1: [[Fraction(0), Fraction(1)], chebyshev_u(5)],
        2: [
            [Fraction(1)],
            poly_mul(
                [Fraction(-1), Fraction(0), Fraction(4)],
                [Fraction(-1), Fraction(0), Fraction(4)],
            ),
        ],
        3: [[Fraction(0), Fraction(0), Fraction(0), Fraction(1)]],
    }

    for spin in range(0, 9):
        cap = 24
        degrees = list(range(spin % 2, cap + 1, 2))
        rows = arity1_rows(spin, degrees)
        space = nullspace(rows, len(degrees))
        expected_dim = len(predicted.get(spin, []))
        check(f"spin {spin}: admissible dimension = {expected_dim}", len(space) == expected_dim)
        if spin in predicted:
            vectors = [
                coefficients_to_vector(p, degrees) for p in predicted[spin]
            ]
            check(
                f"spin {spin}: admissible space matches the closed form",
                same_span(space, vectors),
            )
        if spin <= 3:
            rows_tilt = arity1_rows(spin, degrees, include_tilt_check=True)
            space_tilt = nullspace(rows_tilt, len(degrees))
            check(
                f"spin {spin}: tilted circles through the root add no condition",
                same_span(space, space_tilt),
            )

    # Degree stability: the same spaces already at cap 8.
    for spin, expected_dim in ((0, 2), (1, 2), (2, 2), (3, 1), (4, 0)):
        degrees = list(range(spin % 2, 9, 2))
        space = nullspace(arity1_rows(spin, degrees), len(degrees))
        check(
            f"spin {spin}: dimension already {expected_dim} at degree cap 8",
            len(space) == expected_dim,
        )

    # K = (T2 + 1/3) + (T6 + 1/3): the kernel is the sum of the two
    # admissible spin-0 leaves.
    K = [Fraction(-4, 3), 0, Fraction(20), 0, Fraction(-48), 0, Fraction(32)]
    K = [Fraction(c) for c in K]
    total = poly_add(predicted[0][0], predicted[0][1])
    check(
        "K = (T2 + 1/3) + (T6 + 1/3)",
        poly_add(K, poly_scale(total, -1)) == [Fraction(0)] * 7,
    )

    # The spin-1 admissibles are exactly the tangential derivatives of the
    # spin-0 admissibles: {f' : f in spin0} spans {t, U5}.
    derivative_span = [
        [Fraction(i + 1) * f[i + 1] for i in range(len(f) - 1)]
        for f in predicted[0]
    ]
    deg_cap = 6
    degrees1 = list(range(1, deg_cap, 2))
    check(
        "spin-1 space = derivatives of the spin-0 space",
        same_span(
            [coefficients_to_vector(p, degrees1) for p in derivative_span],
            [coefficients_to_vector(p, degrees1) for p in predicted[1]],
        ),
    )

    print("Part B: unrooted families on the zero-energy family")

    # Harmonic energies g_l = int int P_l: only l = 2 vanishes identically.
    g2 = pair_energy(2)
    check("g_2(mu*) == 0 identically (h_2 face)", not g2)
    for order in (4, 6):
        values = pair_energy(order)
        check(
            f"g_{order}(mu*) != 0 (harmonic_{order} blocks are dead)",
            ((), ) and values.get((), CF_ZERO) != CF_ZERO,
        )

    # The target E = -4/3 + 20 p2 - 48 p4 + 32 p6 vanishes identically.
    target_weights = {
        0: Fraction(-4, 3),
        2: Fraction(20),
        4: Fraction(-48),
        6: Fraction(32),
    }
    target_value: dict[TauMonomial, CF] = {}
    moments = {p: pair_power_moment(p) for p in target_weights}
    for power, weight in target_weights.items():
        for tau_key, value in moments[power].items():
            total = cf_add(
                target_value.get(tau_key, CF_ZERO), cf_scale(value, weight)
            )
            if cf_is_zero(total):
                target_value.pop(tau_key, None)
            else:
                target_value[tau_key] = total
    check("T(mu*) == 0 identically in the free nu-modes", not target_value)

    # Pole-equator SOS identity with free weight w and free nu-hat(2), nu-hat(6):
    # E = 6 (w - 1/3)^2 + (1-w)^2 (|nu2|^2 + |nu6|^2).  Checked as an exact
    # polynomial identity in w and the now-unconstrained tau_{+-2}, tau_{+-6}
    # via the same machinery, at rational sample weights.
    K_weights = {
        0: Fraction(-4, 3),
        2: Fraction(20),
        4: Fraction(-48),
        6: Fraction(32),
    }
    unconstrained: frozenset[int] = frozenset()

    def energy_at_weight(w: Fraction) -> dict[TauMonomial, CF]:
        def integrand_factory(power):
            def integrand(points, slots):
                return l_pow(v_dot(points[0], points[1]), power, slots)

            return integrand

        total: dict[TauMonomial, CF] = {}
        for power, weight in K_weights.items():
            def integrand(points, slots, power=power):
                return l_pow(v_dot(points[0], points[1]), power, slots)

            totals: dict[TauMonomial, CF] = {}
            choices = [("+", "-", "c")] * 2
            for combo in iter_product(*choices):
                mass = Fraction(1)
                points = []
                active = []
                for index, choice in enumerate(combo):
                    if choice == "+":
                        mass *= w / 2
                        points.append(v_pole(2, 1))
                    elif choice == "-":
                        mass *= w / 2
                        points.append(v_pole(2, -1))
                    else:
                        mass *= 1 - w
                        points.append(v_equator(2, index))
                        active.append(index)
                integrand_value = l_scale(integrand(points, 2), mass * weight)
                for (remaining, tau_key), value in pair_slots(
                    integrand_value, active, unconstrained
                ).items():
                    total_value = cf_add(total.get(tau_key, CF_ZERO), value)
                    if cf_is_zero(total_value):
                        total.pop(tau_key, None)
                    else:
                        total[tau_key] = total_value
        return total

    identity_holds = True
    for w in (Fraction(0), Fraction(1, 3), Fraction(1, 2), Fraction(2, 7), Fraction(1)):
        energy = energy_at_weight(w)
        expected: dict[TauMonomial, CF] = {}
        constant = 6 * (w - Fraction(1, 3)) ** 2
        if constant:
            expected[()] = cf(constant)
        square = (1 - w) ** 2
        if square:
            # |nu2|^2 = tau_2 tau_{-2}, |nu6|^2 = tau_6 tau_{-6}
            expected[(-2, 2)] = cf(square)
            expected[(-6, 6)] = cf(square)
        if energy != expected:
            identity_holds = False
    check(
        "pole-equator identity E = 6(w-1/3)^2 + (1-w)^2(|nu2|^2+|nu6|^2)",
        identity_holds,
    )

    # Pair-power (empty type) flags: which combinations of p_a vanish on the
    # family.  Expected: exactly span{3 p_2 - 1, E}.
    pair_degrees = list(range(0, 13, 2))
    rows_map: dict[object, Row] = {}
    for column, power in enumerate(pair_degrees):
        add_condition_rows(
            rows_map,
            {((), key): value for key, value in pair_power_moment(power).items()},
            "pair",
            column,
        )
    space = nullspace(list(rows_map.values()), len(pair_degrees))
    h2_vector = [Fraction(0)] * len(pair_degrees)
    h2_vector[pair_degrees.index(0)] = Fraction(-1)
    h2_vector[pair_degrees.index(2)] = Fraction(3)
    t_vector = [Fraction(0)] * len(pair_degrees)
    for power, weight in target_weights.items():
        t_vector[pair_degrees.index(power)] = weight
    check(
        "pair flags: admissible = span{3 p2 - 1, E} (dim 2 at degree 12)",
        len(space) == 2 and same_span(space, [h2_vector, t_vector]),
    )

    # Spin-2 unrooted flags (deviatoric family).
    spin2_flags: list[tuple[tuple[int, ...], int]] = [((), 0)]
    for a in range(2, 7, 2):
        spin2_flags.append(((a,), 0))
    for a1 in range(1, 7):
        for a2 in range(a1, 7):
            if (a1 + a2) % 2:
                continue
            for b in range(0, 7):
                if (a1 + b) % 2 or (a2 + b) % 2:
                    continue
                if a1 + a2 + b > 6:
                    continue
                spin2_flags.append(((a1, a2), b))
    rows_map = {}
    for column, flag in enumerate(spin2_flags):
        entries = spin2_flag_value(flag)
        for entry_index, totals in enumerate(entries):
            add_condition_rows(
                rows_map,
                {((), key): value for key, value in totals.items()},
                ("spin2", entry_index),
                column,
            )
    spin2_space = nullspace(list(rows_map.values()), len(spin2_flags))
    empty_flag = [Fraction(0)] * len(spin2_flags)
    empty_flag[0] = Fraction(1)
    check(
        "spin-2 flags: the empty (deviatoric) flag is admissible",
        in_span(empty_flag, spin2_space),
    )
    print(
        f"    spin-2 unrooted flags: {len(spin2_space)} admissible directions "
        f"out of {len(spin2_flags)} (degree cap 6)"
    )

    # Harmonic-weighted pair flags of order 2, 4, 6.  A weight combination
    # c survives iff the inner average F_c(x) = int (sum_a c_a (x.y)^a) dmu*(y)
    # pairs to zero against the order-l harmonic moments of mu*.  Sufficient:
    # F_c vanishes on supp mu*, i.e. c is an (E1)-admissible spin-0 radial --
    # so span{T2 + 1/3, T6 + 1/3} always survives.  For order 2 the constant
    # weight also survives (mu* is isotropic).  The check below shows these
    # sufficient conditions are exactly the admissible spaces.
    f1_weights = coefficients_to_vector(
        poly_add(chebyshev_t(2), [Fraction(1, 3)]), list(range(0, 7, 2))
    )
    f2_weights = coefficients_to_vector(
        poly_add(chebyshev_t(6), [Fraction(1, 3)]), list(range(0, 7, 2))
    )
    const_weights = coefficients_to_vector([Fraction(1)], list(range(0, 7, 2)))
    for order, weight_cap in ((2, 6), (4, 6), (6, 6)):
        weights = list(range(0, weight_cap + 1, 2))
        matrices = harmonic_pair_flag_matrix(order, weights)
        rows_list: list[Row] = []
        for tau_key, matrix in matrices.items():
            for row_values in matrix:
                row = {
                    index: value
                    for index, value in enumerate(row_values)
                    if value != 0
                }
                if row:
                    rows_list.append(row)
        space = nullspace(rows_list, len(weights))
        if order == 2:
            check(
                "harmonic_flag order 2: admissible = span{1, T2+1/3, T6+1/3}",
                same_span(space, [const_weights, f1_weights, f2_weights]),
            )
        else:
            check(
                f"harmonic_flag order {order}: admissible = span{{T2+1/3, T6+1/3}}",
                same_span(space, [f1_weights, f2_weights]),
            )

    print("Part C: two-root leaves")

    # Even sector at a fixed working degree.
    working_degree = 8
    basis_even = two_root_basis(working_degree, 0)
    rows_even = arity2_rows(basis_even, det_sector=False)
    space_even = nullspace(rows_even, len(basis_even))
    print(
        f"    even sector, degree {working_degree}: admissible dimension "
        f"{len(space_even)} of {len(basis_even)}"
    )

    # Membership: the arity-1 embeddings f1, f2 in either root variable.
    f1 = poly_add(chebyshev_t(2), [Fraction(1, 3)])
    f2 = poly_add(chebyshev_t(6), [Fraction(1, 3)])
    member_checks = []
    for f, name in ((f1, "T2+1/3"), (f2, "T6+1/3")):
        for slot in (0, 1):
            poly = {}
            for power, value in enumerate(f):
                if value == 0:
                    continue
                key = (power, 0, 0) if slot == 0 else (0, power, 0)
                poly[key] = value
            member_checks.append(
                evaluate_rows_on_vector(
                    rows_even, two_root_vector(poly, basis_even)
                )
            )
    check(
        "even sector contains f(t1) and f(t2) for f = T2+1/3, T6+1/3",
        all(member_checks),
    )

    # Membership: the Gram-determinant ideal family D * t1 t2 * (anything).
    # deg D = 3 (the 2 t1 t2 s term), so multipliers go up to degree - 5.
    ideal_ok = True
    seed = poly3_mul(GRAM3, {(1, 1, 0): Fraction(1)})
    for (i, j, k) in two_root_basis(working_degree - 5, 0):
        candidate = poly3_mul(seed, {(i, j, k): Fraction(1)})
        vector = two_root_vector(candidate, basis_even)
        if not evaluate_rows_on_vector(rows_even, vector):
            ideal_ok = False
            break
    check(
        "even sector contains D * t1 t2 * m for every admissible-parity m",
        ideal_ok,
    )

    # Products f_a(t1) f_b(t2) are NOT admissible (the Delsarte obstruction
    # reappears as a mode-4 violation).
    product_poly: dict[tuple[int, int, int], Fraction] = {}
    for p1, v1 in enumerate(f1):
        for p2, v2 in enumerate(f1):
            if v1 and v2:
                key = (p1, p2, 0)
                product_poly[key] = product_poly.get(key, Fraction(0)) + v1 * v2
    check(
        "even sector does NOT contain (T2+1/3)(t1) (T2+1/3)(t2)",
        not evaluate_rows_on_vector(
            rows_even, two_root_vector(product_poly, basis_even)
        ),
    )

    # The Gram-determinant ideal part in closed form: D * Q is admissible iff
    #   Q(0, 0, s) == 0,
    #   Q(0, t, 0) in span{2 t^4 - t^2},   Q(t, 0, 0) in span{2 t^4 - t^2}
    # (the last from sin^2 th Q = beta cos 2th + gamma cos 6th with zero mean
    # and (T2 - T6)/(1 - t^2) = 8 t^2 (2 t^2 - 1)).
    q_degree = working_degree - 3
    q_basis = two_root_basis(q_degree, 0)
    # computed side: {Q : D*Q satisfies the configuration conditions}
    rows_on_q: list[Row] = []
    a_rows = rows_even
    for row in a_rows:
        new_row: Row = {}
        for q_col, q_mono in enumerate(q_basis):
            candidate = poly3_mul(GRAM3, {q_mono: Fraction(1)})
            total = Fraction(0)
            for mono, value in candidate.items():
                col = basis_even.index(mono)
                if col in row:
                    total += row[col] * value
            if total:
                new_row[q_col] = total
        if new_row:
            rows_on_q.append(new_row)
    q_space_computed = nullspace(rows_on_q, len(q_basis))
    # predicted side
    predicted_rows: list[Row] = []
    for q_col, (i, j, k) in enumerate(q_basis):
        if i == 0 and j == 0:
            predicted_rows.append({q_col: Fraction(1)})
    for slice_index in (0, 1):
        slice_coefficients: dict[int, dict[int, Fraction]] = {}
        for q_col, (i, j, k) in enumerate(q_basis):
            if k != 0:
                continue
            if slice_index == 0 and i == 0:
                slice_coefficients.setdefault(j, {})[q_col] = Fraction(1)
            if slice_index == 1 and j == 0:
                slice_coefficients.setdefault(i, {})[q_col] = Fraction(1)
        for power, cols in slice_coefficients.items():
            if power in (2, 4):
                continue
            predicted_rows.append(dict(cols))
        row: Row = {}
        for col, value in slice_coefficients.get(2, {}).items():
            row[col] = row.get(col, Fraction(0)) + 2 * value
        for col, value in slice_coefficients.get(4, {}).items():
            row[col] = row.get(col, Fraction(0)) + value
        if row:
            predicted_rows.append(row)
    q_space_predicted = nullspace(predicted_rows, len(q_basis))
    check(
        "even sector ideal part: D*Q admissible iff Q(0,0,s)=0 and both "
        "slices lie in span{2t^4 - t^2}",
        same_span(q_space_computed, q_space_predicted),
    )

    # Modulated generators.  C_n is the canonical polynomial lift of
    # cos(2 theta - n delta): C_0 = T2(t1), C_1 = 2 t1 t2 - s, and the
    # Chebyshev recurrence C_{n+1} = 2 s C_n - C_{n-1}.  S_n lifts
    # cos(6 theta - n delta): S_0 = T6(t1), S_1 = 2 T4(t1) C_1 - C_{-1}.
    # The corrected elements are  C^_n = C_n + T_n(s)/3  and likewise for S.
    def chebyshev_in_s(n: int) -> dict[tuple[int, int, int], Fraction]:
        return {
            (0, 0, power): value
            for power, value in enumerate(chebyshev_t(n))
            if value != 0
        }

    def lift_family(seed0, seed1, count):
        lifts = [seed0, seed1]
        two_s = {(0, 0, 1): Fraction(2)}
        while len(lifts) < count:
            nxt = poly3_mul(two_s, lifts[-1])
            for key, value in lifts[-2].items():
                nxt[key] = nxt.get(key, Fraction(0)) - value
                if nxt[key] == 0:
                    del nxt[key]
            lifts.append(nxt)
        return lifts

    c0 = {(0, 0, 0): Fraction(-1), (2, 0, 0): Fraction(2)}
    c1 = {(1, 1, 0): Fraction(2), (0, 0, 1): Fraction(-1)}
    c_family = lift_family(c0, c1, 9)
    # C_{-1} = 2 s C_0 - C_1
    c_minus1 = poly3_mul({(0, 0, 1): Fraction(2)}, c0)
    for key, value in c1.items():
        c_minus1[key] = c_minus1.get(key, Fraction(0)) - value
        if c_minus1[key] == 0:
            del c_minus1[key]
    t4_of_t1 = {
        (power, 0, 0): value
        for power, value in enumerate(chebyshev_t(4))
        if value != 0
    }
    s1 = poly3_mul(poly3_mul({(0, 0, 0): Fraction(2)}, t4_of_t1), c1)
    for key, value in c_minus1.items():
        s1[key] = s1.get(key, Fraction(0)) - value
        if s1[key] == 0:
            del s1[key]
    s0 = {
        (power, 0, 0): value
        for power, value in enumerate(chebyshev_t(6))
        if value != 0
    }
    s_family = lift_family(s0, s1, 9)

    modulated_degree = 12
    basis_mod = two_root_basis(modulated_degree, 0)
    rows_mod = arity2_rows(basis_mod, det_sector=False)
    c_pattern = []
    s_pattern = []
    for family, pattern, count in (
        (c_family, c_pattern, 9),
        (s_family, s_pattern, 8),  # deg S^_n = n + 5 for n >= 1
    ):
        for n in range(count):
            candidate = dict(family[n])
            for key, value in chebyshev_in_s(n).items():
                candidate[key] = (
                    candidate.get(key, Fraction(0)) + value / 3
                )
                if candidate[key] == 0:
                    del candidate[key]
            vector = two_root_vector(candidate, basis_mod)
            pattern.append(
                evaluate_rows_on_vector(rows_mod, vector)
            )
    check(
        "modulated families C^_n = C_n + T_n(s)/3 and S^_n = S_n + T_n(s)/3 "
        "are admissible for every n",
        all(c_pattern) and all(s_pattern),
    )

    # Completeness of the closed-form description at the working degree:
    # admissible = span{C^_n, S^_n} + {D Q - (1/3)(1-s^2) Q(0,0,s) : Q}_adm.
    def lift_backward(family, seed_next, count):
        # family[0] = C_0, family[1] = C_1; produce C_{-1}, ..., C_{-count}
        two_s = {(0, 0, 1): Fraction(2)}
        lifts = []
        previous, current = seed_next, family[0]
        for _ in range(count):
            nxt = poly3_mul(two_s, current)
            for key, value in previous.items():
                nxt[key] = nxt.get(key, Fraction(0)) - value
                if nxt[key] == 0:
                    del nxt[key]
            lifts.append(nxt)
            previous, current = current, nxt
        return lifts

    def poly3_degree(poly):
        return max((sum(key) for key in poly), default=0)

    generators: list[list[Fraction]] = []
    c_negative = lift_backward(c_family, c_family[1], 7)
    s_negative = lift_backward(s_family, s_family[1], 3)
    for family, negative in ((c_family, c_negative), (s_family, s_negative)):
        members = list(family) + list(negative)
        orders = list(range(len(family))) + [-(n + 1) for n in range(len(negative))]
        for n, member in zip(orders, members):
            candidate = dict(member)
            for key, value in chebyshev_in_s(abs(n)).items():
                candidate[key] = candidate.get(key, Fraction(0)) + value / 3
                if candidate[key] == 0:
                    del candidate[key]
            if poly3_degree(candidate) > working_degree:
                continue
            vector = two_root_vector(candidate, basis_even)
            assert evaluate_rows_on_vector(rows_even, vector)
            generators.append(vector)

    def mq_image(q_mono):
        i, j, k = q_mono
        image = poly3_mul(GRAM3, {q_mono: Fraction(1)})
        if i == 0 and j == 0:
            for key, value in (
                ((0, 0, k), Fraction(-1, 3)),
                ((0, 0, k + 2), Fraction(1, 3)),
            ):
                image[key] = image.get(key, Fraction(0)) + value
                if image[key] == 0:
                    del image[key]
        return image

    rows_on_q = []
    for row in rows_even:
        new_row: Row = {}
        for q_col, q_mono in enumerate(q_basis):
            total = Fraction(0)
            for mono, value in mq_image(q_mono).items():
                col = basis_even.index(mono)
                if col in row:
                    total += row[col] * value
            if total:
                new_row[q_col] = total
        if new_row:
            rows_on_q.append(new_row)
    for q_vector in nullspace(rows_on_q, len(q_basis)):
        image: dict[tuple[int, int, int], Fraction] = {}
        for q_col, weight in enumerate(q_vector):
            if weight == 0:
                continue
            for mono, value in mq_image(q_basis[q_col]).items():
                image[mono] = image.get(mono, Fraction(0)) + weight * value
                if image[mono] == 0:
                    del image[mono]
        vector = two_root_vector(image, basis_even)
        assert evaluate_rows_on_vector(rows_even, vector)
        generators.append(vector)

    def rank_of(vectors: list[list[Fraction]]) -> int:
        reduced: list[tuple[int, Row]] = []
        for vector in vectors:
            work = {i: v for i, v in enumerate(vector) if v != 0}
            for pivot, pivot_row in reduced:
                coefficient = work.get(pivot)
                if coefficient is None:
                    continue
                for col, val in pivot_row.items():
                    updated = work.get(col, Fraction(0)) - coefficient * val
                    if updated == 0:
                        work.pop(col, None)
                    else:
                        work[col] = updated
            if work:
                pivot = min(work)
                inverse = 1 / work[pivot]
                reduced.append(
                    (pivot, {c: v * inverse for c, v in work.items()})
                )
        return len(reduced)

    generator_rank = rank_of(generators)
    mq_dimension = len(nullspace(rows_on_q, len(q_basis)))
    print(
        f"    canonical closed-form generators span {generator_rank} of "
        f"{len(space_even)} admissible dimensions at degree {working_degree}"
    )
    print(
        f"    (M(s), D*Q) kernel family dimension: {mq_dimension}"
    )

    # Which delta-modulations of the circle modes 2 and 6 are realizable at
    # all inside the admissible space (profile realizability).
    def b2_profile(vector, basis):
        slots = 2
        x2 = [l_cos(slots, 1), l_sin(slots, 1), {}]
        y = v_equator(slots, 0)
        x1 = [l_one(slots), {}, {}]
        t1 = v_dot(x1, y)
        t2 = v_dot(x2, y)
        s_val = v_dot(x1, x2)
        total: Laurent = {}
        for value, (i, j, k) in zip(vector, basis):
            if value == 0:
                continue
            term = l_mul(
                l_mul(l_pow(t1, i, slots), l_pow(t2, j, slots)),
                l_pow(s_val, k, slots),
            )
            total = l_add(total, l_scale(term, value))
        return total

    realizable = {2: set(), 6: set()}
    for vector in nullspace(rows_mod, len(basis_mod)):
        profile = b2_profile(vector, basis_mod)
        for (m_theta, m_delta), value in profile.items():
            if m_theta in (2, 6) and not cf_is_zero(value):
                realizable[m_theta].add(-m_delta)
    for mode in (2, 6):
        mods = sorted(realizable[mode])
        print(f"    realizable delta-modulations of circle mode {mode}: {mods}")

    # Completeness theorem for the even sector: the kernel of the circle
    # profile map (A restricted to the coplanar circle, coefficients of the
    # theta-modes 2 and 6 as functions of delta) is exactly the
    # (M(s), D*Q) family, so
    #     dim admissible = dim{Q : D Q - (1/3)(1-s^2) Q(0,0,s) admissible}
    #                      + rank(profile map).
    profile_vectors = []
    profile_index: dict[tuple[int, int, str], int] = {}
    for vector in space_even:
        profile = b2_profile(vector, basis_even)
        coordinates: dict[int, Fraction] = {}
        for (m_theta, m_delta), value in profile.items():
            if m_theta not in (2, 6):
                continue
            for part, component in (("re", value[0]), ("im", value[1])):
                if component == 0:
                    continue
                index = profile_index.setdefault(
                    (m_theta, m_delta, part), len(profile_index)
                )
                coordinates[index] = component
        profile_vectors.append(coordinates)
    profile_rank = len(
        nullspace([], 0)
    )  # placeholder replaced by rank computation below
    reduced_rank: list[tuple[int, Row]] = []
    for coordinates in profile_vectors:
        work = dict(coordinates)
        for pivot, pivot_row in reduced_rank:
            coefficient = work.get(pivot)
            if coefficient is None:
                continue
            for col, val in pivot_row.items():
                updated = work.get(col, Fraction(0)) - coefficient * val
                if updated == 0:
                    work.pop(col, None)
                else:
                    work[col] = updated
        if work:
            pivot = min(work)
            inverse = 1 / work[pivot]
            reduced_rank.append(
                (pivot, {c: v * inverse for c, v in work.items()})
            )
    profile_rank = len(reduced_rank)
    print(
        f"    profile-map rank at degree {working_degree}: {profile_rank}"
    )
    check(
        "even sector: dim admissible = dim (M,Q)-family + profile rank",
        len(space_even) == mq_dimension + profile_rank,
    )

    # ------------------------------------------------------------------
    # Identities displayed in docs/TWO_ROOT_GENERATORS.md.
    # ------------------------------------------------------------------

    # Explicit expansions of the canonical mode-2 lifts.
    explicit_c = {
        2: {(1, 1, 1): Fraction(4), (0, 0, 2): Fraction(-2),
            (2, 0, 0): Fraction(-2), (0, 0, 0): Fraction(1)},
        3: {(1, 1, 2): Fraction(8), (0, 0, 3): Fraction(-4),
            (2, 0, 1): Fraction(-4), (1, 1, 0): Fraction(-2),
            (0, 0, 1): Fraction(3)},
        4: {(1, 1, 3): Fraction(16), (0, 0, 4): Fraction(-8),
            (2, 0, 2): Fraction(-8), (1, 1, 1): Fraction(-8),
            (0, 0, 2): Fraction(8), (2, 0, 0): Fraction(2),
            (0, 0, 0): Fraction(-1)},
        -1: {(2, 0, 1): Fraction(4), (1, 1, 0): Fraction(-2),
             (0, 0, 1): Fraction(-1)},
    }
    check(
        "explicit expansions of C_2, C_3, C_4, C_-1",
        all(
            (c_family[n] if n >= 0 else c_negative[-n - 1]) == poly
            for n, poly in explicit_c.items()
        ),
    )

    # Lift ambiguity: C_2 = T_2(t_2) + 2 D.
    t2_of_t2 = {(0, 2, 0): Fraction(2), (0, 0, 0): Fraction(-1)}
    difference = dict(c_family[2])
    for key, value in t2_of_t2.items():
        difference[key] = difference.get(key, Fraction(0)) - value
        if difference[key] == 0:
            del difference[key]
    check(
        "C_2 = T_2(t_2) + 2 D",
        difference == {key: 2 * value for key, value in GRAM3.items()},
    )

    # Closed forms of the mode-6 seeds in second-kind Chebyshev polynomials.
    def u_poly(n: int) -> dict[int, Fraction]:
        return {
            power: value
            for power, value in enumerate(chebyshev_u(n))
            if value != 0
        }

    s1_nice: dict[tuple[int, int, int], Fraction] = {}
    for power, value in u_poly(5).items():
        s1_nice[(power, 1, 0)] = value
    for power, value in u_poly(4).items():
        s1_nice[(power, 0, 1)] = s1_nice.get((power, 0, 1), Fraction(0)) - value
    s_minus1_nice: dict[tuple[int, int, int], Fraction] = {}
    for power, value in u_poly(6).items():
        s_minus1_nice[(power, 0, 1)] = value
    for power, value in u_poly(5).items():
        s_minus1_nice[(power, 1, 0)] = (
            s_minus1_nice.get((power, 1, 0), Fraction(0)) - value
        )
    check(
        "S_1 = U_5(t1) t2 - U_4(t1) s and S_-1 = U_6(t1) s - U_5(t1) t2",
        s_family[1] == s1_nice and s_negative[0] == s_minus1_nice,
    )

    # Pole-value lemma: C_n(0,0,s) = S_n(0,0,s) = -T_n(s), which fixes the
    # +T_n(s)/3 correction uniformly.
    def pole_slice(poly) -> dict[int, Fraction]:
        return {
            key[2]: value
            for key, value in poly.items()
            if key[0] == 0 and key[1] == 0
        }

    pole_lemma = True
    for n in range(-2, 6):
        member = c_family[n] if n >= 0 else c_negative[-n - 1]
        expected = {
            power: -value
            for power, value in enumerate(chebyshev_t(abs(n)))
            if value != 0
        }
        if pole_slice(member) != expected:
            pole_lemma = False
    for n in range(-1, 4):
        member = s_family[n] if n >= 0 else s_negative[-n - 1]
        expected = {
            power: -value
            for power, value in enumerate(chebyshev_t(abs(n)))
            if value != 0
        }
        if pole_slice(member) != expected:
            pole_lemma = False
    check("pole-value lemma: C_n(0,0,s) = S_n(0,0,s) = -T_n(s)", pole_lemma)

    # Negative-index corrected lifts are admissible too.
    negative_ok = True
    for member, order in (
        (c_negative[0], -1),
        (c_negative[1], -2),
        (s_negative[0], -1),
    ):
        candidate = dict(member)
        for key, value in chebyshev_in_s(abs(order)).items():
            candidate[key] = candidate.get(key, Fraction(0)) + value / 3
            if candidate[key] == 0:
                del candidate[key]
        vector = two_root_vector(candidate, basis_mod)
        if not evaluate_rows_on_vector(rows_mod, vector):
            negative_ok = False
    check("C^_-1, C^_-2, S^_-1 are admissible", negative_ok)

    # General one-step lifts, as exact identities on the coplanar circle:
    #   cos(m th - d) = U_{m-1}(t1) t2 - U_{m-2}(t1) s,
    #   cos(m th + d) = U_m(t1) s - U_{m-1}(t1) t2.
    def circle_restriction(poly) -> Laurent:
        slots = 2
        x1 = [l_one(slots), {}, {}]
        x2 = [l_cos(slots, 1), l_sin(slots, 1), {}]
        y = v_equator(slots, 0)
        t1v, t2v, sv = v_dot(x1, y), v_dot(x2, y), v_dot(x1, x2)
        total: Laurent = {}
        for (i, j, k), value in poly.items():
            term = l_mul(
                l_mul(l_pow(t1v, i, slots), l_pow(t2v, j, slots)),
                l_pow(sv, k, slots),
            )
            total = l_add(total, l_scale(term, value))
        return total

    lifts_ok = True
    for m in range(2, 8):
        minus_lift: dict[tuple[int, int, int], Fraction] = {}
        for power, value in u_poly(m - 1).items():
            minus_lift[(power, 1, 0)] = value
        for power, value in u_poly(m - 2).items():
            minus_lift[(power, 0, 1)] = (
                minus_lift.get((power, 0, 1), Fraction(0)) - value
            )
        plus_lift: dict[tuple[int, int, int], Fraction] = {}
        for power, value in u_poly(m).items():
            plus_lift[(power, 0, 1)] = value
        for power, value in u_poly(m - 1).items():
            plus_lift[(power, 1, 0)] = (
                plus_lift.get((power, 1, 0), Fraction(0)) - value
            )
        for lift, sign in ((minus_lift, -1), (plus_lift, 1)):
            profile = circle_restriction(lift)
            expected_profile = {
                (m, sign): cf(Fraction(1, 2)),
                (-m, -sign): cf(Fraction(1, 2)),
            }
            if profile != expected_profile:
                lifts_ok = False
    check(
        "one-step lifts cos(m th -+ d) = U_{m-+1}... for m = 2..7",
        lifts_ok,
    )

    # Closed-form solution of the recurrence (basis T_n(s), U_{n-1}(s)):
    #   C_n = C_0 T_n(s) + 2 t1 (t2 - s t1) U_{n-1}(s),
    #   S_n = S_0 T_n(s) + (S_1 - s S_0) U_{n-1}(s),
    # whence |C_n| <= 1 + 4n and |S_n| <= 1 + 12n on the Gram body: the
    # theta atoms sum_n q^{n^2} C^_n converge for every |q| < 1.
    def s_poly(coefficients: list[Fraction]) -> dict[tuple[int, int, int], Fraction]:
        return {
            (0, 0, power): value
            for power, value in enumerate(coefficients)
            if value != 0
        }

    closed_ok = True
    for family in (c_family, s_family):
        seed0, seed1 = family[0], family[1]
        s_seed0 = poly3_mul({(0, 0, 1): Fraction(1)}, seed0)
        residue = dict(seed1)
        for key, value in s_seed0.items():
            residue[key] = residue.get(key, Fraction(0)) - value
            if residue[key] == 0:
                del residue[key]
        for n in range(0, 6):
            expected = poly3_mul(seed0, s_poly(chebyshev_t(n)))
            if n >= 1:
                for key, value in poly3_mul(
                    residue, s_poly(chebyshev_u(n - 1))
                ).items():
                    expected[key] = expected.get(key, Fraction(0)) + value
                    if expected[key] == 0:
                        del expected[key]
            if expected != family[n]:
                closed_ok = False
    check(
        "closed form C_n = C_0 T_n(s) + 2 t1 (t2 - s t1) U_{n-1}(s) "
        "(and the S analogue)",
        closed_ok,
    )

    # Odd (orientation) sector: the coplanar configuration imposes nothing;
    # all conditions come from the pole-root slices.
    basis_odd = two_root_basis(working_degree, 1)
    rows_odd_b2 = arity2_rows(basis_odd, det_sector=True, configs=("b2",))
    check(
        "odd sector: the two-equator-root configuration imposes no condition",
        not any(row for row in rows_odd_b2),
    )
    rows_odd = arity2_rows(basis_odd, det_sector=True)
    space_odd = nullspace(rows_odd, len(basis_odd))
    print(
        f"    odd sector, degree {working_degree}: admissible dimension "
        f"{len(space_odd)} of {len(basis_odd)}"
    )

    # Odd-sector codimension in closed form: only the two pole-root slices
    # constrain B, forcing B(0,.,0) and B(.,0,0) into span{U1, U5}:
    # codim = 2 * #{even m not in {2,6}: 4 <= m <= degree + 1}.
    codim_ok = True
    for degree in range(2, 11):
        basis_o = two_root_basis(degree, 1)
        space_o = nullspace(
            arity2_rows(basis_o, det_sector=True), len(basis_o)
        )
        expected_codim = 2 * len(
            [m for m in range(4, degree + 2, 2) if m not in (2, 6)]
        )
        if len(basis_o) - len(space_o) != expected_codim:
            codim_ok = False
    check(
        "odd sector codim = 2 #{even m, 4 <= m <= deg+1, m != 6} for deg 2..10",
        codim_ok,
    )

    if args.table:
        print("Part C tables: admissible dimension by degree")
        print("    degree | even: dim/total | odd: dim/total")
        for degree in range(2, args.max_table_degree + 1):
            basis_e = two_root_basis(degree, 0)
            space_e = nullspace(
                arity2_rows(basis_e, det_sector=False), len(basis_e)
            )
            basis_o = two_root_basis(degree, 1)
            space_o = nullspace(
                arity2_rows(basis_o, det_sector=True), len(basis_o)
            )
            print(
                f"    {degree:6d} | {len(space_e):4d}/{len(basis_e):4d}"
                f"      | {len(space_o):4d}/{len(basis_o):4d}"
            )

        # The four parity sectors used by sos_search.py (flag degree = total
        # SDP degree // 2): even_00 = all exponents even, even_11 = i,j odd
        # and k even, odd_01 / odd_10 = the orientation-odd sectors.
        def sector_basis(max_degree, name):
            selected = []
            for i in range(max_degree + 1):
                for j in range(max_degree + 1 - i):
                    for k in range(max_degree + 1 - i - j):
                        ij, ik, jk = (i + j) % 2, (i + k) % 2, (j + k) % 2
                        if name == "even_00" and (ij, ik, jk) == (0, 0, 0):
                            selected.append((i, j, k))
                        if name == "even_11" and (ij, ik, jk) == (0, 1, 1):
                            selected.append((i, j, k))
                        if name == "odd_01" and (ij, ik, jk) == (1, 0, 1):
                            selected.append((i, j, k))
                        if name == "odd_10" and (ij, ik, jk) == (1, 1, 0):
                            selected.append((i, j, k))
            return selected

        print("    solver sectors (flag degree = SDP degree // 2):")
        print(
            "    flag degree | even_00 | even_11 | odd_01 | odd_10"
        )
        for degree in range(2, args.max_table_degree + 1):
            report = []
            for name in ("even_00", "even_11", "odd_01", "odd_10"):
                basis = sector_basis(degree, name)
                det_sector = name.startswith("odd")
                space = nullspace(
                    arity2_rows(basis, det_sector=det_sector), len(basis)
                )
                report.append(f"{len(space):3d}/{len(basis):3d}")
            print(
                f"    {degree:11d} | "
                + " | ".join(report)
            )

    print("Part D: the weighted target h_2 E (zero set adds all isotropic measures)")

    # A sharp all-measures certificate of h_2 E >= 0 must vanish on every
    # isotropic measure and on the pole-equator family.  For a pure square,
    # vanishing at Haar measure and its isotropic perturbations forces the
    # leaf, as a function of y, to be a pure degree-2 spherical harmonic
    # (see docs/E1_ADMISSIBLE.md §6): h_2-multiplied blocks are exempt (the
    # h_2 factor already vanishes on both families).  In the spin-k sector
    # the leaf radial f must satisfy
    #     f(t) (1-t^2)^{k/2}  in  span{ P_2^k(t) },
    # i.e. orthogonality to every associated Legendre P_l^k with l != 2.
    import math

    def poly_derivative(coefficients: list[Fraction], times: int) -> list[Fraction]:
        result = list(coefficients)
        for _ in range(times):
            result = [
                Fraction(power + 1) * value
                for power, value in enumerate(result[1:])
            ]
        return result

    def monomial_integral(power: int) -> Fraction:
        # int_{-1}^{1} t^power dt
        return Fraction(2, power + 1) if power % 2 == 0 else Fraction(0)

    def one_minus_t2_power(exponent: int) -> list[Fraction]:
        coefficients = [Fraction(0)] * (2 * exponent + 1)
        for m in range(exponent + 1):
            coefficients[2 * m] = Fraction((-1) ** m * math.comb(exponent, m))
        return coefficients

    def isotropic_rows(spin: int, degrees: list[int]) -> list[Row]:
        rows: list[Row] = []
        max_degree = max(degrees) + 4
        weight = one_minus_t2_power(spin)
        for order in range(spin + (spin % 2), max_degree + 1, 2):
            if order == 2:
                continue
            derivative = poly_derivative(legendre_coefficients(order), spin)
            kernel = poly_mul(weight, derivative) if derivative else []
            row: Row = {}
            for column, degree in enumerate(degrees):
                value = sum(
                    (
                        coefficient * monomial_integral(power + degree)
                        for power, coefficient in enumerate(kernel)
                    ),
                    Fraction(0),
                )
                if value != 0:
                    row[column] = value
            if row:
                rows.append(row)
        return rows

    weighted_predictions = {
        0: [[Fraction(-1, 3), Fraction(0), Fraction(1)]],  # t^2 - 1/3
        1: [[Fraction(0), Fraction(1)]],  # t
        2: [[Fraction(1)]],  # 1
        3: [],
    }
    for spin in range(0, 5):
        cap = 12
        degrees = list(range(spin % 2, cap + 1, 2))
        rows = arity1_rows(spin, degrees) + isotropic_rows(spin, degrees)
        space = nullspace(rows, len(degrees))
        predicted_vectors = [
            coefficients_to_vector(polynomial, degrees)
            for polynomial in weighted_predictions.get(spin, [])
        ]
        check(
            f"weighted spin {spin}: admissible = "
            + (
                "span{"
                + {0: "t^2 - 1/3", 1: "t", 2: "1"}.get(spin, "")
                + "}"
                if spin <= 2
                else "{0}"
            ),
            same_span(space, predicted_vectors),
        )

    # The three survivors are exactly the spin components of the deviatoric
    # tensor int (y y^T - I/3) dmu rooted at x, so the pure-square layer of
    # the weighted certificate is the deviatoric flag family and nothing
    # else; every other block must carry the h_2 factor explicitly (the
    # --h2-localized-all structure).

    # Pair-power layer: combinations of p_a vanishing on all isotropic
    # measures.  In Legendre coordinates p_a = sum_l lambda_{a,l} g_l with
    # g_0 = 1, g_2 free to be zero, g_{l>=4} free: conditions kill every
    # coefficient except the g_2 direction.
    pair_degrees = list(range(0, 13, 2))
    legendre_rows: list[Row] = []
    max_order = max(pair_degrees)
    for order in [0] + list(range(4, max_order + 1, 2)):
        legendre = legendre_coefficients(order)
        row: Row = {}
        for column, power in enumerate(pair_degrees):
            coefficient = sum(
                (
                    value * monomial_integral(power + lpower)
                    for lpower, value in enumerate(legendre)
                ),
                Fraction(0),
            )
            if coefficient != 0:
                row[column] = coefficient
        if row:
            legendre_rows.append(row)
    weighted_pair_space = nullspace(legendre_rows, len(pair_degrees))
    h2_direction = [Fraction(0)] * len(pair_degrees)
    h2_direction[pair_degrees.index(0)] = Fraction(-1)
    h2_direction[pair_degrees.index(2)] = Fraction(3)
    check(
        "weighted pair flags: admissible = span{3 p2 - 1} (E drops out)",
        same_span(weighted_pair_space, [h2_direction]),
    )

    print("Part E: the spin-2 operator gap I - A_2 >= 0")

    # A_2(mu) is the average of pi_2(rho_x) = (M -> rho_x M rho_x) over mu,
    # acting on the 5-dimensional space of symmetric traceless matrices.
    # Each pi_2(rho_x) is an orthogonal involution, so -I <= A_2 <= I is a
    # valid constraint for every measure.  Claims verified here:
    #   (i)  the axial quadrupole diag(1,1,-2) is an eigenvalue-1
    #        eigenvector of A_2(mu*) for EVERY member of the zero family
    #        (identically in the free nu-modes): I - A_2 is active on the
    #        whole family, so its multipliers survive sharpness;
    #   (ii) at the ONB, (A_2 - 1)(A_2 + 1/3) = 0 with trace 1:
    #        spectrum {1, 1, -1/3, -1/3, -1/3}.
    sym_basis = [
        [(0, 0, Fraction(1)), (1, 1, Fraction(-1))],
        [(0, 0, Fraction(1)), (1, 1, Fraction(1)), (2, 2, Fraction(-2))],
        [(0, 1, Fraction(1)), (1, 0, Fraction(1))],
        [(0, 2, Fraction(1)), (2, 0, Fraction(1))],
        [(1, 2, Fraction(1)), (2, 1, Fraction(1))],
    ]

    def basis_matrix(index, slots):
        matrix = [[{} for _ in range(3)] for _ in range(3)]
        for (r, c, value) in sym_basis[index]:
            matrix[r][c] = l_const(slots, cf(value))
        return matrix

    def rho_of(point, slots):
        # 2 x x^T - I with Laurent entries
        matrix = [[l_scale(l_mul(point[r], point[c]), Fraction(2)) for c in range(3)] for r in range(3)]
        for d in range(3):
            matrix[d][d] = l_add(matrix[d][d], l_const(slots, cf(-1)))
        return matrix

    def mat_mul(a, b):
        return [
            [
                l_add(l_add(l_mul(a[r][0], b[0][c]), l_mul(a[r][1], b[1][c])), l_mul(a[r][2], b[2][c]))
                for c in range(3)
            ]
            for r in range(3)
        ]

    def sym_coordinates(matrix):
        # coordinates in the (non-orthogonal) sym_basis: solve by structure
        # basis duals: b0 = (m00 - m11)/2, b1 = (m00 + m11 - 2 m22)/6,
        # b2 = m01, b3 = m02, b4 = m12  (valid on traceless symmetric input)
        m = matrix
        half = Fraction(1, 2)
        sixth = Fraction(1, 6)
        return [
            l_scale(l_add(m[0][0], l_scale(m[1][1], Fraction(-1))), half),
            l_scale(
                l_add(l_add(m[0][0], m[1][1]), l_scale(m[2][2], Fraction(-2))),
                sixth,
            ),
            m[0][1],
            m[0][2],
            m[1][2],
        ]

    # (i) family-wide activity, identically in the free nu-modes
    slots = 1
    pole = v_pole(slots, 1)
    equator_point = v_equator(slots, 0)
    axial = basis_matrix(1, slots)  # diag(1,1,-2)
    rho_pole = rho_of(pole, slots)
    rho_eq = rho_of(equator_point, slots)
    image_pole = mat_mul(mat_mul(rho_pole, axial), rho_pole)
    image_eq = mat_mul(mat_mul(rho_eq, axial), rho_eq)
    activity_ok = True
    for r in range(3):
        for c in range(3):
            pole_part = l_scale(image_pole[r][c], Fraction(1, 3))
            paired = pair_slots(l_scale(image_eq[r][c], Fraction(2, 3)), [0])
            total: dict = {}
            for (key, tau_key), value in paired.items():
                total[tau_key] = cf_add(total.get(tau_key, CF_ZERO), value)
            for key, value in pole_part.items():
                total[()] = cf_add(total.get((), CF_ZERO), value)
            expected = axial[r][c].get((0,), CF_ZERO)
            total[()] = cf_add(total.get((), CF_ZERO), cf_scale(expected, Fraction(-1)))
            if any(not cf_is_zero(v) for v in total.values()):
                activity_ok = False
    check(
        "A_2(mu*)[diag(1,1,-2)] = diag(1,1,-2) identically on the zero family",
        activity_ok,
    )

    # (ii) ONB spectrum via the minimal polynomial
    onb_matrix = [[Fraction(0)] * 5 for _ in range(5)]
    for axis in range(3):
        point = [l_const(0, cf(1)) if d == axis else {} for d in range(3)]
        rho_axis = rho_of(point, 0)
        for a in range(5):
            image = mat_mul(mat_mul(rho_axis, basis_matrix(a, 0)), rho_axis)
            for b, coordinate in enumerate(sym_coordinates(image)):
                value = coordinate.get((), CF_ZERO)
                assert value[1] == 0
                onb_matrix[b][a] += Fraction(1, 3) * value[0]
    def mat5_mul(a, b):
        return [
            [sum(a[r][k] * b[k][c] for k in range(5)) for c in range(5)]
            for r in range(5)
        ]
    identity5 = [[Fraction(1) if r == c else Fraction(0) for c in range(5)] for r in range(5)]
    shifted_down = [
        [onb_matrix[r][c] - identity5[r][c] for c in range(5)] for r in range(5)
    ]
    shifted_up = [
        [onb_matrix[r][c] + Fraction(1, 3) * identity5[r][c] for c in range(5)]
        for r in range(5)
    ]
    product = mat5_mul(shifted_down, shifted_up)
    trace = sum(onb_matrix[d][d] for d in range(5))
    check(
        "A_2(ONB): (A_2 - 1)(A_2 + 1/3) = 0 and tr = 1 "
        "(spectrum {1,1,-1/3,-1/3,-1/3})",
        all(value == 0 for row in product for value in row) and trace == 1,
    )

    failed = [name for name, ok in checks if not ok]
    print()
    if failed:
        print(f"{len(checks) - len(failed)}/{len(checks)} checks passed; FAILED:")
        for name in failed:
            print(f"  - {name}")
        raise SystemExit(1)
    print(f"All {len(checks)} checks passed.")


if __name__ == "__main__":
    main()
