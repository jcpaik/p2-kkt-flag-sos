"""Exact rational counterexample to ``J-108F-288A4 >= 0``.

The escape is a very small insertion into a non-Haar pole--equator zero
measure, so ordinary atom optimizers readily miss it.  Every coordinate,
weight, and calculation below is rational.
"""

from fractions import Fraction as Q

from tensor_pair_circle_graph_exact import (
    determinant_integrand,
    j_integrand,
    polynomial_terms,
    projected_circle_integrand,
)


POINTS = (
    (Q(0), Q(0), Q(1)),
    (Q(1), Q(0), Q(0)),
    (Q(0), Q(1), Q(0)),
    (Q(3, 5), Q(4, 5), Q(0)),
    (Q(-4, 5), Q(3, 5), Q(0)),
    (Q(36, 901), Q(48, 901), Q(899, 901)),
)
BASE_WEIGHTS = (
    Q(1, 3), Q(3, 100), Q(3, 100), Q(91, 300), Q(91, 300), Q(0),
)
EPSILON = Q(1, 100_000_000)
WEIGHTS = tuple((1-EPSILON)*weight for weight in BASE_WEIGHTS[:-1]) + (EPSILON,)


def expectation(expression, weights=WEIGHTS):
    gram = tuple(
        tuple(sum(x*y for x, y in zip(left, right, strict=True))
              for right in POINTS)
        for left in POINTS
    )
    terms = polynomial_terms(expression)
    result = Q(0)
    for i in range(6):
        for j in range(6):
            for k in range(6):
                for l in range(6):
                    edges = (
                        gram[i][j], gram[i][k], gram[i][l],
                        gram[j][k], gram[j][l], gram[k][l],
                    )
                    value = Q(0)
                    for coefficient, exponents in terms:
                        monomial = coefficient
                        for edge, exponent in zip(edges, exponents, strict=True):
                            monomial *= edge**exponent
                        value += monomial
                    result += weights[i]*weights[j]*weights[k]*weights[l]*value
    return result


def gateaux_derivative(expression):
    """Exact derivative from the five-atom zero face toward atom 6.

    The expectation is represented with four iid labels even when the
    expression ignores some labels.  Differentiating the four weight factors
    is therefore uniform for J, F, and A4.
    """
    base = expectation(expression, BASE_WEIGHTS)
    terms = polynomial_terms(expression)
    gram = tuple(
        tuple(sum(x*y for x, y in zip(left, right, strict=True))
              for right in POINTS)
        for left in POINTS
    )
    derivative = -4*base
    for insertion_slot in range(4):
        subtotal = Q(0)
        for i in range(5):
            for j in range(5):
                for k in range(5):
                    old_indices = (i, j, k)
                    indices = []
                    cursor = 0
                    for slot in range(4):
                        if slot == insertion_slot:
                            indices.append(5)
                        else:
                            indices.append(old_indices[cursor])
                            cursor += 1
                    edges = (
                        gram[indices[0]][indices[1]],
                        gram[indices[0]][indices[2]],
                        gram[indices[0]][indices[3]],
                        gram[indices[1]][indices[2]],
                        gram[indices[1]][indices[3]],
                        gram[indices[2]][indices[3]],
                    )
                    value = Q(0)
                    for coefficient, exponents in terms:
                        monomial = coefficient
                        for edge, exponent in zip(edges, exponents, strict=True):
                            monomial *= edge**exponent
                        value += monomial
                    subtotal += (BASE_WEIGHTS[i]*BASE_WEIGHTS[j]
                                 *BASE_WEIGHTS[k]*value)
        derivative += subtotal
    return derivative


def audit():
    assert sum(BASE_WEIGHTS) == 1
    assert sum(WEIGHTS) == 1
    assert all(sum(entry*entry for entry in point) == 1 for point in POINTS)

    # The first five atoms are a convex mixture of the ONBs
    # (e1,e2,e3) and (f1,f2,e3), so the unperturbed target vanishes.
    assert expectation(j_integrand(), BASE_WEIGHTS) == 0
    assert expectation(determinant_integrand(), BASE_WEIGHTS) == 0
    assert expectation(projected_circle_integrand(2), BASE_WEIGHTS) == 0

    J = expectation(j_integrand())
    F = expectation(determinant_integrand())
    A4 = expectation(projected_circle_integrand(2))
    gap = J-108*F-288*A4
    assert J > 0 and F > 0 and A4 > 0
    assert gap < Q(-4, 100_000_000_000_000)

    print("J,F,A4", float(J), float(F), float(A4))
    print("gap", float(gap))
    print("exact gap", gap)

    Jp = gateaux_derivative(j_integrand())
    Fp = gateaux_derivative(determinant_integrand())
    A4p = gateaux_derivative(projected_circle_integrand(2))
    numerator = Jp-108*Fp
    critical_coefficient = numerator/A4p
    assert numerator > 0 and A4p > 0 and critical_coefficient < 288
    print("Gateaux J',F',A4'", Jp, Fp, A4p)
    print("Gateaux (J-108F)'", numerator, float(numerator))
    print("critical A4 coefficient", critical_coefficient,
          float(critical_coefficient))


if __name__ == "__main__":
    audit()
