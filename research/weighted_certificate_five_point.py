"""Explore the fully symmetrized five-point kernel for F(mu) E(mu).

For five projective points, average h on a triple times K on its
complementary pair.  Its iid expectation is exactly F(mu) E(mu).
"""

import itertools
import numpy as np


def K(t):
    return 32*t**6 - 48*t**4 + 20*t**2 - 4/3


def h_from_gram(g, tri):
    i, j, k = tri
    a, b, c = g[i, j], g[i, k], g[j, k]
    D = 1 + 2*a*b*c - a*a - b*b - c*c
    cyc = (
        a*a*b*b*(c-a*b)**2
        + a*a*c*c*(b-a*c)**2
        + b*b*c*c*(a-b*c)**2
    )
    return (8/3)*D*cyc


def sym_kernel(x):
    g = x @ x.T
    value = 0.0
    terms = []
    all_indices = set(range(5))
    for tri in itertools.combinations(range(5), 3):
        pair = tuple(sorted(all_indices.difference(tri)))
        term = h_from_gram(g, tri) * K(g[pair])
        value += term
        terms.append(term)
    return value / 10, np.array(terms)


def stationary_four_kernel(x):
    """Four-point symmetrization of h(X,Y,Z)K(X,W).

    At a stationary measure, its expectation equals F(mu)E(mu), because
    U_mu(X)=E(mu) on the support.  The expression is averaged over the
    choice of the three-point h-kernel and over the anchored vertex in that
    triple.
    """
    g = x @ x.T
    total = 0.0
    terms = []
    all_indices = set(range(4))
    for tri in itertools.combinations(range(4), 3):
        (outside,) = all_indices.difference(tri)
        hv = h_from_gram(g, tri)
        for anchor in tri:
            term = hv * K(g[anchor, outside])
            total += term
            terms.append(term)
    return total / 12, np.array(terms)


def main():
    rng = np.random.default_rng(20260820)
    best = (np.inf, None, None)
    count_negative = 0
    for _ in range(200_000):
        x = rng.normal(size=(5, 3))
        x /= np.linalg.norm(x, axis=1, keepdims=True)
        value, terms = sym_kernel(x)
        if value < best[0]:
            best = (value, x.copy(), terms.copy())
        count_negative += value < -1e-12
    print("negative", count_negative, "best", best[0])
    if best[1] is not None:
        print(best[1])
        print(best[2])

    best4 = (np.inf, None, None)
    count_negative4 = 0
    for _ in range(200_000):
        x = rng.normal(size=(4, 3))
        x /= np.linalg.norm(x, axis=1, keepdims=True)
        value, terms = stationary_four_kernel(x)
        if value < best4[0]:
            best4 = (value, x.copy(), terms.copy())
        count_negative4 += value < -1e-12
    print("stationary4 negative", count_negative4, "best", best4[0])
    print(best4[1])
    print(best4[2])


if __name__ == "__main__":
    main()
