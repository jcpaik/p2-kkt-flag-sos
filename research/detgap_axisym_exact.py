"""Exact tests of the normalized spin-2 determinant-gap conjecture.

For an axisymmetric projective measure, write z=x_3 and let
    m_k = E[z^k].
The normalized operator is Bhat=(5/4)(I-A2), equivalently
    <S,Bhat S> = 5 E[|Sx|^2-(x^T S x)^2]
on Frobenius-normalized traceless symmetric matrices.  This script
derives its three SO(2)-isotypic eigenvalues and tests harmonic-density
families exactly with SymPy.
"""

import sympy as sp


z, phi, eps = sp.symbols("z phi eps", real=True)
r = sp.sqrt(1 - z**2)
x = sp.Matrix([r * sp.cos(phi), r * sp.sin(phi), z])

Q0 = sp.diag(-1, -1, 2) / sp.sqrt(6)
Q1 = sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / sp.sqrt(2)
Q2 = sp.diag(1, -1, 0) / sp.sqrt(2)


def sphere_average_axisym(expr, density):
    """Average against density(z) relative to Haar probability."""
    return sp.factor(
        sp.integrate(
            sp.integrate(sp.expand_trig(expr), (phi, 0, 2 * sp.pi))
            * density
            / (4 * sp.pi),
            (z, -1, 1),
        )
    )


def bhat_eigenvalue(Q, density):
    qx = (x.T * Q * x)[0]
    sx2 = (x.T * Q**2 * x)[0]
    return sp.factor(5 * sphere_average_axisym(sx2 - qx**2, density))


def report(density, energy):
    vals = [bhat_eigenvalue(Q, density) for Q in (Q0, Q1, Q2)]
    det = sp.factor(vals[0] * vals[1] ** 2 * vals[2] ** 2)
    gap = sp.factor(energy - sp.Rational(32, 105) * det**2)
    print("density =", density)
    print("lambda_0, lambda_1, lambda_2 =", vals)
    print("det =", det)
    print("E =", energy)
    print("gap factor =", gap)
    print("gap series =", sp.series(gap, eps, 0, 6))


P2 = sp.legendre(2, z)
P4 = sp.legendre(4, z)

# If density=1+eps P_l, the squared Legendre moment is eps^2/(2l+1)^2.
E_P4 = sp.Rational(32, 105) - sp.Rational(384, 385) * eps**2 / 81
E_P2 = sp.Rational(32, 105) + sp.Rational(8, 7) * eps**2 / 25

report(1 + eps * P4, E_P4)
print()
report(1 + eps * P2, E_P2)
