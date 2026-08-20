#!/usr/bin/env python3
"""Exact symbolic checks for research/tensor_zero_chord.md."""

import sympy as sp


q, u, r = sp.symbols("q u r", real=True)
s = 1 - q

A = 2 * q**2 - 4 * q + 3
C = 5 * q**2 - 2 * q + 3
B = 7 * q**2 + 2 * q + 3

# Upper bound for L after Re(n*mbar)<=u and Re(m^2)<=u^2.
L_upper = 2 * C - 2 * A * r - B * u**2 + 2 * u + s**2 * u**2
remainder = sp.expand(8 * (3 - 2 * q * s) * (1 - r) - L_upper)

assert sp.simplify(
    sp.diff(remainder, r) + 2 * (6 * q**2 - 4 * q + 9)
) == 0
assert sp.simplify(
    remainder.subs(r, u)
    - 2
    * (u - 1)
    * ((3 * q**2 + 2 * q + 1) * u - 3 * q**2 + 6 * q - 9)
) == 0
assert sp.simplify(
    ((3 * q**2 + 2 * q + 1) * u - 3 * q**2 + 6 * q - 9).subs(u, 1)
    - 8 * (q - 1)
) == 0

# Potential expansion from the pole and equatorial parts.
K = lambda x: 32 * x**3 - 48 * x**2 + 20 * x - sp.Rational(4, 3)
equator = (
    -sp.Rational(4, 3)
    + 10 * s
    - 6 * s**2 * (3 + r)
    + s**3 * (10 + 6 * r)
)
potential = sp.factor(sp.Rational(1, 3) * K(q) + sp.Rational(2, 3) * equator)
assert sp.simplify(potential - 4 * q * s**2 * (1 - r)) == 0

print("zero-face chord identities verified exactly")
