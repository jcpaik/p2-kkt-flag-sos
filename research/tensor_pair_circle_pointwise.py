"""Audit the symmetrized four-sample integrand for the A4 target."""

import itertools

import numpy as np
import sympy as sp

from tensor_pair_circle_graph_exact import VARIABLES, expansion


EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EDGE_INDEX = {edge: index for index, edge in enumerate(EDGES)}


def permuted(expression, permutation):
    substitution = {}
    for index, edge in enumerate(EDGES):
        image = tuple(sorted((permutation[edge[0]], permutation[edge[1]])))
        substitution[VARIABLES[index]] = VARIABLES[EDGE_INDEX[image]]
    return expression.xreplace(substitution)


def symmetrized_target():
    expression, _ = expansion()
    return sp.expand(
        sum(permuted(expression, permutation)
            for permutation in itertools.permutations(range(4))) / 24
    )


def scan(samples=200_000, seed=20260822):
    polynomial = symmetrized_target()
    function = sp.lambdify(VARIABLES, polynomial, "numpy")
    rng = np.random.default_rng(seed)
    minimum = None
    witness = None
    for _ in range((samples + 9999) // 10000):
        points = rng.normal(size=(min(10000, samples), 4, 3))
        points /= np.linalg.norm(points, axis=2)[:, :, None]
        gram = np.einsum("bik,bjk->bij", points, points)
        values = function(*(gram[:, edge[0], edge[1]] for edge in EDGES))
        index = int(np.argmin(values))
        if minimum is None or values[index] < minimum:
            minimum = float(values[index])
            witness = gram[index]
        samples -= len(points)
        if samples <= 0:
            break
    print("monomials", len(sp.Poly(polynomial, *VARIABLES).terms()))
    print("minimum", minimum)
    print("gram", witness)
    return polynomial, minimum, witness


if __name__ == "__main__":
    scan()
