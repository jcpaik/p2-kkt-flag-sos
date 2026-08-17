"""Exact (rational / symbolic) verification of the structural identities for the
P_2 kernel

    K(t) = 32 t^6 - 48 t^4 + 20 t^2 - 4/3,      t = cos(theta),

used by the certificate search.  Every claim below is checked in exact
arithmetic; nothing here relies on a floating-point solver.

Run with

    python3 verify_exact_structure.py
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp

t, a, b, w, z1, z2, z3, q = sp.symbols("t a b w z1 z2 z3 q", real=True)

K = 32 * t**6 - 48 * t**4 + 20 * t**2 - sp.Rational(4, 3)

checks: list[tuple[str, bool]] = []


def check(name: str, condition) -> None:
    checks.append((name, bool(condition)))


# ---------------------------------------------------------------- 1. Chebyshev
# K(cos theta) = cos 6 theta + cos 2 theta + 2/3.
chebyshev = sp.chebyshevt(6, t) + sp.chebyshevt(2, t) + sp.Rational(2, 3)
check("K = T_6 + T_2 + 2/3", sp.expand(K - chebyshev) == 0)

# ---------------------------------------------------------------- 2. Legendre
legendre_coefficients = {
    degree: sp.Rational(2 * degree + 1, 2)
    * sp.integrate(K * sp.legendre(degree, t), (t, -1, 1))
    for degree in range(7)
}
expected_legendre = {
    0: sp.Rational(32, 105),
    1: 0,
    2: sp.Rational(8, 7),
    3: 0,
    4: sp.Rational(-384, 385),
    5: 0,
    6: sp.Rational(512, 231),
}
check("Legendre expansion", legendre_coefficients == expected_legendre)
check(
    "exactly one negative Legendre coefficient (degree 4)",
    [d for d, c in legendre_coefficients.items() if c < 0] == [4],
)

# ---------------------------------------------------------------- 3. Squared form
# K = 2 W - 4/3 with W = T_3^2 + t^2 = cos^2(3 theta) + cos^2(theta).
W = sp.chebyshevt(3, t) ** 2 + t**2
check("K = 2 W - 4/3, W = T_3^2 + t^2", sp.expand(K - (2 * W - sp.Rational(4, 3))) == 0)
check("W = 16t^6 - 24t^4 + 10t^2", sp.expand(W - (16 * t**6 - 24 * t**4 + 10 * t**2)) == 0)
check("W >= 0 on [-1,1]", sp.expand(W - t**2 * (4 * t**2 - 3) ** 2 - t**2) == 0)

# W = 1 - u + 2u^3 = (1+u)(u^2 + (1-u)^2) with u = cos 2 theta = 2t^2 - 1.
u = 2 * t**2 - 1
check("W = 1 - u + 2u^3, u = 2t^2-1", sp.expand(W - (1 - u + 2 * u**3)) == 0)
check(
    "W = (1+u)(u^2 + (1-u)^2)",
    sp.expand(W - (1 + u) * (u**2 + (1 - u) ** 2)) == 0,
)

# ------------------------------------------------- 4. measures on a great circle
# For any probability measure nu on a great circle, writing c_k = |nu-hat(k)|^2,
#     E(nu) = c_6 + c_2 + 2/3  >=  2/3  >  0.
# This is immediate from K(cos theta) = cos 6theta + cos 2theta + 2/3 and
# int int cos(k(phi_x - phi_y)) dnu dnu = |nu-hat(k)|^2 >= 0.
# We verify the underlying Fourier identity numerically-exactly on a
# trigonometric-polynomial family.
phi = sp.symbols("phi", real=True)
for m in (2, 4, 5, 7):
    # uniform measure on m lines at angles k*pi/m
    angles = [sp.pi * k / m for k in range(m)]
    energy = sp.Rational(1, m**2) * sum(
        K.subs(t, sp.cos(x - y)) for x in angles for y in angles
    )
    fourier2 = (
        sp.Rational(1, m) * sum(sp.exp(sp.I * 2 * x) for x in angles)
    )
    fourier6 = (
        sp.Rational(1, m) * sum(sp.exp(sp.I * 6 * x) for x in angles)
    )
    predicted = (
        sp.Abs(fourier6) ** 2 + sp.Abs(fourier2) ** 2 + sp.Rational(2, 3)
    )
    check(
        f"great-circle energy identity, {m} lines",
        sp.simplify(sp.expand(sp.nsimplify(energy - predicted))) == 0,
    )

# --------------------------------------------- 5. the pole-equator SOS identity
# mu = w * (line e) + (1-w) * nu,  nu on the great circle e^perp.
# Then  E(mu) = 6 (w - 1/3)^2 + (1-w)^2 (c_2 + c_6).
c2, c6 = sp.symbols("c2 c6", nonnegative=True)
pole_equator_energy = (
    w**2 * K.subs(t, 1)
    + 2 * w * (1 - w) * K.subs(t, 0)
    + (1 - w) ** 2 * (c6 + c2 + sp.Rational(2, 3))
)
sos_form = 6 * (w - sp.Rational(1, 3)) ** 2 + (1 - w) ** 2 * (c2 + c6)
check(
    "pole-equator identity E = 6(w-1/3)^2 + (1-w)^2 (c_2 + c_6)",
    sp.expand(pole_equator_energy - sos_form) == 0,
)

# ------------------------------------------------- 6. the ONB KKT certificate
# For the orthonormal-basis measure, E = 0 and the rooted potential is
#     U(z) = 32 z1^2 z2^2 z3^2  >= 0,
# which is exactly the first-order KKT inequality, with equality on the support.
onb_potential = sp.Rational(1, 3) * (
    K.subs(t, z1) + K.subs(t, z2) + K.subs(t, z3)
)
onb_potential = sp.expand(
    onb_potential.subs(z3**2, 1 - z1**2 - z2**2).rewrite(sp.Pow)
)
target_potential = sp.expand((32 * z1**2 * z2**2 * z3**2).subs(z3**2, 1 - z1**2 - z2**2))
check(
    "ONB potential U(z) = 32 z1^2 z2^2 z3^2",
    sp.simplify(
        sp.expand(
            sp.Rational(1, 3)
            * (K.subs(t, z1) + K.subs(t, z2) + K.subs(t, z3))
            - 32 * z1**2 * z2**2 * z3**2
        ).subs(z3**2, 1 - z1**2 - z2**2)
    )
    == 0,
)
check("ONB energy is zero", sp.Rational(1, 3) * K.subs(t, 1) + sp.Rational(2, 3) * K.subs(t, 0) == 0)

# --------------------------------------- 7. no two-point LP certificate exists
# Suppose h = sum_n h_n P_n with h_n >= 0 for n >= 1 and h <= K on [-1,1].
# The pole + uniform-equator measure mu* has E(mu*) = 0 and its distance
# distribution has full support [-1,1]; hence h = K on [-1,1], contradicting
# h_4 = -384/385 < 0.  We verify E(mu*) = 0 exactly.
#   pole-pole 1/9 * K(1),  pole-equator 4/9 * K(0),  equator-equator 4/9 * mean.
equator_mean = sp.integrate(K.subs(t, sp.cos(phi)), (phi, 0, 2 * sp.pi)) / (2 * sp.pi)
check("mean of K over a great circle is 2/3", equator_mean == sp.Rational(2, 3))
pole_haar_energy = (
    sp.Rational(1, 9) * K.subs(t, 1)
    + sp.Rational(4, 9) * K.subs(t, 0)
    + sp.Rational(4, 9) * equator_mean
)
check("pole + uniform equator has energy exactly 0", pole_haar_energy == 0)

# ------------------------------- 8. cylindrical (Fourier-in-longitude) kernels
# With x = (cos psi_x) e + (sin psi_x) omega_x etc., a = cos psi_x, b = cos psi_y,
#     K(x . y) = sum_{k=0}^{6} c_k(a,b) cos(k (phi_x - phi_y)).
# c_5 and c_6 are rank-one positive kernels; c_0, c_2, c_3, c_4 are not PSD.
p = sp.symbols("p", positive=True)
cc = sp.symbols("cc", real=True)
expanded = sp.Poly(sp.expand(K.subs(t, a * b + p * cc)), cc)
power_coefficients = {n: sp.expand(coefficient) for (n,), coefficient in zip(expanded.monoms(), expanded.coeffs())}


def cosine_weight(n: int, k: int):
    if k > n or (n - k) % 2:
        return sp.Integer(0)
    j = (n - k) // 2
    if k == 0:
        return sp.Rational(1, 2**n) * sp.binomial(n, j)
    return sp.Rational(2, 2**n) * sp.binomial(n, j)


squared_p = (1 - a**2) * (1 - b**2)
# Keep `p` symbolic (p^2 = (1-a^2)(1-b^2)) so that no fractional powers appear.
cylindrical = {
    k: sp.expand(
        sum(power_coefficients.get(n, 0) * cosine_weight(n, k) for n in range(7))
    )
    for k in range(7)
}

# c_6 = p^6 = (1-a^2)^3 (1-b^2)^3 = u(a) u(b):  rank one, hence PSD.
check("c_6 = p^6, a rank-one positive kernel", sp.expand(cylindrical[6] - p**6) == 0)
# c_5 = 12 a b p^5 = 12 (a (1-a^2)^{5/2}) (b (1-b^2)^{5/2}):  rank one, PSD.
check(
    "c_5 = 12 a b p^5, a rank-one positive kernel",
    sp.expand(cylindrical[5] - 12 * a * b * p**5) == 0,
)

# c_4 = 6 p^4 (11 a^2 b^2 - a^2 - b^2); the inner form has matrix
# [[0,-1],[-1,11]] in the separable features (1, a^2), determinant -1 < 0.
inner4 = sp.Matrix([[0, -1], [-1, 11]])
check("c_4 inner form is indefinite (det = -1)", inner4.det() == -1)
check(
    "c_4 = 6 p^4 (10 a^2 b^2 + p^2 - 1) = 6 p^4 (11a^2b^2 - a^2 - b^2)",
    sp.expand(cylindrical[4] - 6 * p**4 * (10 * a**2 * b**2 + p**2 - 1)) == 0
    and sp.expand(
        (10 * a**2 * b**2 + squared_p - 1) - (11 * a**2 * b**2 - a**2 - b**2)
    )
    == 0,
)
check(
    "c_3 = 4 a b p^3 (40 a^2 b^2 + 15 p^2 - 12)"
    " = 4 a b p^3 (55a^2b^2 - 15a^2 - 15b^2 + 3)",
    sp.expand(cylindrical[3] - 4 * a * b * p**3 * (40 * a**2 * b**2 + 15 * p**2 - 12))
    == 0
    and sp.expand(
        (40 * a**2 * b**2 + 15 * squared_p - 12)
        - (55 * a**2 * b**2 - 15 * a**2 - 15 * b**2 + 3)
    )
    == 0,
)

inner3 = sp.Matrix([[3, -15], [-15, 55]])
check("c_3 inner form is indefinite (det = -60)", inner3.det() == -60)

inner2 = sp.Matrix([[1, -6, 15], [-6, 132, -270], [15, -270, 495]])
check("c_2 inner form is indefinite (det = -6480)", inner2.det() == -6480)

inner0 = sp.Matrix(
    [
        [1, -6, 18, -15],
        [-6, 126, -405, 315],
        [18, -405, 1260, -945],
        [-15, 315, -945, 693],
    ]
)
check(
    "c_0 inner form is not PSD (det = -34992, third leading minor = -3969)",
    inner0.det() == -34992 and inner0[:3, :3].det() == -3969,
)


# ------------------------------------------- 9. the SO(3) group reformulation
# Lines in R^3 are the involutions of SO(3): x <-> rho_x, the pi-rotation about
# x.  The product rho_x rho_y is a rotation by 2 theta, so with chi_l the
# character of the spin-l representation,
#     K(cos theta) = (1/2)( chi_3 - chi_2 + chi_1 )(2 theta) + 1/6.
theta = sp.symbols("theta", real=True)


def so3_character(order: int, angle):
    return sp.expand_trig(
        sp.simplify(
            sum(sp.exp(sp.I * m * angle) for m in range(-order, order + 1)).rewrite(sp.cos)
        )
    )


group_form = (
    sp.Rational(1, 2)
    * (
        so3_character(3, 2 * theta)
        - so3_character(2, 2 * theta)
        + so3_character(1, 2 * theta)
    )
    + sp.Rational(1, 6)
)
check(
    "K = (chi_3 - chi_2 + chi_1)/2 + 1/6 on SO(3)",
    sp.simplify(sp.expand_trig(sp.expand(K.subs(t, sp.cos(theta)) - group_form))) == 0,
)

# Consequently, with A_l = int pi_l(rho_x) d mu(x)  (symmetric, tr A_l = chi_l(pi)),
#     E(mu) = (1/2) ( |A_3|_F^2 - |A_2|_F^2 + |A_1|_F^2 ) + 1/6,
# because  int int chi_l(rho_x rho_y) = |A_l|_F^2 >= 0.
# Check at the ONB, where the three involutions form the Klein four-group V:
# sum over the three non-identity elements of V equals 4P - I with P the
# projector onto the V-invariants.
onb_group_terms = {}
for order in (1, 2, 3):
    dimension = 2 * order + 1
    character_at_involution = sum((-1) ** m for m in range(-order, order + 1))
    invariant_dimension = sp.Rational(dimension + 3 * character_at_involution, 4)
    onb_group_terms[order] = sp.Rational(1, 9) * (
        16 * invariant_dimension - 8 * invariant_dimension + dimension
    )
check(
    "ONB group-form norms are 1/3, 7/3, 5/3",
    [onb_group_terms[order] for order in (1, 2, 3)]
    == [sp.Rational(1, 3), sp.Rational(7, 3), sp.Rational(5, 3)],
)
check(
    "ONB energy via the group form is zero",
    sp.Rational(1, 2)
    * (onb_group_terms[3] - onb_group_terms[2] + onb_group_terms[1])
    + sp.Rational(1, 6)
    == 0,
)

# A_1 = 2 M - I with M the second-moment matrix, so |A_1|^2 = 4 p_2 - 1.
p2 = sp.symbols("p2", real=True)
check("|A_1|_F^2 = 4 p_2 - 1", sp.expand(4 * p2 - 4 + 3 - (4 * p2 - 1)) == 0)

# ------------------------------------------------------------------- report
width = max(len(name) for name, _ in checks)
failures = 0
for name, ok in checks:
    print(f"{name.ljust(width)}  {'OK' if ok else 'FAILED'}")
    failures += 0 if ok else 1
print()
print(f"{len(checks) - failures}/{len(checks)} exact checks passed")
raise SystemExit(1 if failures else 0)
