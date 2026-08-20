"""Exact obstruction to bounding the A4 residual by ``J(T_mu,T_mu)``.

This does *not* refute ``J-108F-288A4 >= 0``.  It refutes the tempting
Schur strengthening with coefficient 27/10 (and hence every larger
coefficient).  All points and weights below are rational.
"""

from fractions import Fraction as Q

from tensor_pair_circle_graph_exact import (
    determinant_integrand,
    frame_self_integrand,
    j_integrand,
    polynomial_terms,
    projected_circle_integrand,
)


STEREOGRAPHIC = (
    (Q(-3181, 5000), Q(33, 1000)),
    (Q(-2521, 5000), Q(-749, 625)),
    (Q(107, 2500), Q(4597, 5000)),
    (Q(-6339, 1000), Q(17873, 5000)),
    (Q(502, 125), Q(93, 200)),
    (Q(-1019, 625), Q(1133, 5000)),
)
WEIGHTS = (
    Q(8409, 100000),
    Q(8623, 100000),
    Q(29231, 100000),
    Q(1871, 25000),
    Q(2741, 10000),
    Q(18843, 100000),
)


def sphere_point(p, q):
    denominator = 1+p*p+q*q
    return (2*p/denominator, 2*q/denominator, (p*p+q*q-1)/denominator)


def expectation(terms, gram):
    result = Q(0)
    size = len(WEIGHTS)
    for i in range(size):
        for j in range(size):
            for k in range(size):
                for l in range(size):
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
                    result += WEIGHTS[i]*WEIGHTS[j]*WEIGHTS[k]*WEIGHTS[l]*value
    return result


def audit():
    points = tuple(sphere_point(*pair) for pair in STEREOGRAPHIC)
    gram = tuple(tuple(sum(x*y for x, y in zip(left, right, strict=True))
                       for right in points) for left in points)
    values = {
        "J": expectation(polynomial_terms(j_integrand()), gram),
        "F": expectation(polynomial_terms(determinant_integrand()), gram),
        "A4": expectation(polynomial_terms(projected_circle_integrand(2)), gram),
        "JTT": expectation(polynomial_terms(frame_self_integrand()), gram),
    }
    gap = values["J"]-108*values["F"]-288*values["A4"]
    strengthened = gap-Q(27, 10)*values["JTT"]
    assert sum(WEIGHTS) == 1
    assert gap > Q(1, 200)                 # 0.0054527...
    assert values["JTT"] > Q(1, 500)      # 0.0022969...
    assert strengthened < Q(-7, 10000)    # -0.0007490...
    print({name: float(value) for name, value in values.items()})
    print("gap", float(gap), "strengthened", float(strengthened))


if __name__ == "__main__":
    audit()
