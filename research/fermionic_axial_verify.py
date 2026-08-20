"""Exact symbolic verifier for fermionic_axial_proof.md."""

from __future__ import annotations

import sys

sys.path.insert(0, "research")

import sympy as sp

import fermionic_axial_scalar_exact as axial


(m1, m2, m3), averaged, f, delta, squares = axial.invariant_moment_form()
_, _, ublock, cblock, _, _, _ = axial.split_value(averaged)

a = sp.expand(m1 - m2)
ell = sp.expand(3 * m1 - 13 * m2 + 12 * m3)
vee = sp.expand(m1 - 4 * m2 + 3 * m3)

expected_u = ell * sp.diag(1, -1, 1, -1) / 4
assert sp.simplify(ublock - expected_u) == sp.zeros(4)

nonzero_c = [
    sp.factor(cblock[i, j])
    for i in range(4)
    for j in range(6)
    if cblock[i, j] != 0
]
assert len(nonzero_c) == 4
assert all(sp.factor(value**2 - 3 * vee**2 / 16) == 0 for value in nonzero_c)

expected_squares = sp.Rational(1, 2) * ell**2 + 6 * vee**2
assert sp.factor(squares - expected_squares) == 0
assert sp.factor(f[0, 0] - 3 * a) == 0
assert sp.factor(delta - (9 * a - 2)) == 0

t = sp.symbols("t")
l_point = 3 * t - 13 * t**2 + 12 * t**3
v_point = t - 4 * t**2 + 3 * t**3
pointwise_slack = sp.factor(-l_point - 2 * v_point - 9 * t * (1 - t) + 2)
assert sp.expand(pointwise_slack - 2 * (1 - t) * (3 * t - 1) ** 2) == 0

# Exact weighted Cauchy identity:
# q - 3/8(L+2V)^2 = (L-6V)^2/8.
assert sp.factor(
    expected_squares
    - sp.Rational(3, 8) * (ell + 2 * vee) ** 2
    - (ell - 6 * vee) ** 2 / 8
) == 0

print("Ubar =", expected_u)
print("nonzero Cbar entries =", nonzero_c)
print("square term =", sp.factor(expected_squares))
print("delta =", 9 * a - 2)
print("pointwise slack =", pointwise_slack)
print("Cauchy slack =", sp.factor((ell - 6 * vee) ** 2 / 8))
