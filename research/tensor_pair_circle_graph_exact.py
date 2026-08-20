"""Exact four-sample graph expansion of J - 108 F - 288 A4.

Here A4 is the denominator-cleared projected-circle square averaged with
root weight (x.y)^4 (1-(x.y)^2)^6.  The factor (1-a^2)^6 is already
absorbed when the orthogonal projector is cleared, leaving the polynomial
``a^4 * Phi`` below.
"""

from collections import defaultdict
from fractions import Fraction
import sys

import sympy as sp

sys.path.insert(0, ".")
import sos_search as ss


# Edge order used by sos_search: 01, 02, 03, 12, 13, 23.
a, u, r, v, t, c = sp.symbols("a u r v t c")
VARIABLES = (a, u, r, v, t, c)


def projected_circle_integrand(root_power=2):
    U = u**2 + v**2 - 2*a*u*v
    V = r**2 + t**2 - 2*a*r*t
    L = u*r + v*t - a*(u*t + v*r)
    UV = U*V
    return sp.expand(
        a ** (2*root_power)
        * (32*L**6 - 48*L**4*UV + 20*L**2*UV**2 - 2*UV**3)
    )


def determinant_integrand():
    # The unsymmetrized rooted form has the same iid expectation as F.
    D = 1 + 2*a*u*v - a*a - u*u - v*v
    return sp.expand(8*a*a*u*u*(v-a*u)**2*D)


def j_integrand():
    return 144*a**6 - 216*a**4 + 87*a**2 - 5


def frame_self_integrand():
    """Polynomial four-sample representative of ``J(T_mu,T_mu)``.

    The two ordered root pairs are ``(x,y)`` and ``(z,w)``.  On the
    rank-three Gram variety this equals

        576 w_xy w_zw I_J(nu_Rxy, nu_Rzw).

    The apparently rational frame formula has a removable factor
    ``1-c**2``.  We remove it modulo the 4-by-4 Gram determinant; adding a
    multiple of that determinant does not change any spherical evaluation.
    """
    s = 1-a*a
    d = 1-c*c
    Dxyz = s-u*u-v*v+2*a*u*v
    Dxyw = s-r*r-t*t+2*a*r*t
    Dxzw = d-u*u-r*r+2*c*u*r
    Dyzw = d-v*v-t*t+2*c*v*t
    p12 = sp.det(sp.Matrix([[1, a, r], [a, 1, t], [u, v, c]]))
    p34 = sp.det(sp.Matrix([[a, u, r], [v, 1, c], [t, c, 1]]))

    n1 = u*u*(v-a*u)**2*Dxyz
    n2 = (r-c*u)**2*(t-c*v-a*r+a*c*u)**2 * (
        Dxyw+c*c*Dxyz-2*c*p12
    )
    n3 = Dxzw*(Dyzw+a*a*Dxzw-2*a*p34)*(u*t-r*v)**2
    numerator = sp.expand(n2+n3)

    gram4 = sp.det(sp.Matrix([
        [1, a, u, r],
        [a, 1, v, t],
        [u, v, 1, c],
        [r, t, c, 1],
    ]))
    qplus = (r-u)**2*(a*(r-u)-t+v)**2
    qminus = (r+u)**2*(a*(r+u)-t-v)**2
    multiplier = sp.expand(((1+c)*qplus+(1-c)*qminus)/2)
    quotient, remainder = sp.div(
        sp.Poly(sp.expand(numerator-multiplier*gram4), c),
        sp.Poly(d, c),
    )
    assert remainder.is_zero
    # 27648*w_xy*w_zw times the relative-frame column product; the 1/16
    # in w_xy*w_zw leaves the factor 1728 below.
    return sp.expand(1728*a*a*c*c*(d*d*n1+quotient.as_expr()))


def polynomial_terms(expression):
    polynomial = sp.Poly(sp.expand(expression), *VARIABLES)
    terms = []
    for exponent, coefficient in polynomial.terms():
        coefficient = sp.Rational(coefficient)
        terms.append(
            (
                Fraction(int(coefficient.p), int(coefficient.q)),
                tuple(map(int, exponent)),
            )
        )
    return terms


def exact_expectation_vector(expression):
    out = defaultdict(Fraction)
    for coefficient, exponent in polynomial_terms(expression):
        label, reduction = ss.graph_expectation_label(4, exponent)
        if label is not None and reduction:
            out[label] += coefficient*reduction
    return {label: value for label, value in out.items() if value}


def expansion(constant=288, root_power=2, frame_constant=0):
    expression = (
        j_integrand()
        - 108*determinant_integrand()
        - constant*projected_circle_integrand(root_power)
        - frame_constant*frame_self_integrand()
    )
    return expression, exact_expectation_vector(expression)


def audit():
    A = projected_circle_integrand(2)
    F = determinant_integrand()
    G, vector = expansion()
    print("expanded monomials A,F,G", len(sp.Poly(A,*VARIABLES).terms()), len(sp.Poly(F,*VARIABLES).terms()), len(sp.Poly(G,*VARIABLES).terms()))
    print("canonical expectation labels", len(vector))
    by_type = defaultdict(int)
    for label in vector:
        by_type[label[0]] += 1
    print("label types", dict(by_type))
    for label, coefficient in sorted(vector.items(), key=lambda item: str(item[0])):
        print(repr(label), coefficient)

    T = frame_self_integrand()
    H, hvector = expansion(frame_constant=3)
    print("expanded monomials T,H", len(sp.Poly(T,*VARIABLES).terms()), len(sp.Poly(H,*VARIABLES).terms()))
    print("H canonical expectation labels", len(hvector))


if __name__ == "__main__":
    audit()
