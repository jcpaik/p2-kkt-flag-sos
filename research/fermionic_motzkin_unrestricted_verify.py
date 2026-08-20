"""Exact checks for the unrestricted fermionic--Motzkin decomposition."""

from fractions import Fraction as Q
from itertools import combinations
from random import Random

import sympy as sp


def check(weights):
    """`weights[(i,j)]` are rational diagonal wedge weights on K5."""
    edges = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    assert set(weights) == set(edges)
    assert sum(weights.values(), Q(0)) == 1

    degrees = [
        sum((weights[e] for e in edges if i in e), Q(0)) for i in range(5)
    ]
    # Relabeling is immaterial; the caller puts a maximal vertex at zero.
    assert degrees[0] == max(degrees)
    c = degrees[0]
    a = [weights[0, i] for i in range(1, 5)]
    outer = {(i, j): weights[i, j] for i in range(1, 5) for j in range(i + 1, 5)}
    r = [
        sum((value for (i, j), value in outer.items() if k in (i, j)), Q(0))
        for k in range(1, 5)
    ]
    delta = 3 * c - 2
    u0 = [a[i] - r[i] - delta / 4 for i in range(4)]
    disjoint = outer[1, 2] * outer[3, 4]
    disjoint += outer[1, 3] * outer[2, 4]
    disjoint += outer[1, 4] * outer[2, 3]

    adjacency = sum(
        weights[e] * weights[f]
        for n, e in enumerate(edges)
        for f in edges[n + 1 :]
        if set(e) & set(f)
    )
    claimed = Q(1, 3) - adjacency
    decomposed = 2 * disjoint + sum(v * v for v in u0) / 2 - delta * delta / 24
    assert claimed == decomposed


def main():
    rng = Random(20260820)
    edges = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    for _ in range(100):
        raw = [rng.randrange(1, 100) for _ in edges]
        total = sum(raw)
        weights = {edge: Q(value, total) for edge, value in zip(edges, raw)}
        degrees = {
            i: sum((weights[e] for e in edges if i in e), Q(0)) for i in range(5)
        }
        top = max(degrees, key=degrees.get)
        permutation = [top] + [i for i in range(5) if i != top]
        inverse = {old: new for new, old in enumerate(permutation)}
        relabeled = {
            tuple(sorted((inverse[i], inverse[j]))): value
            for (i, j), value in weights.items()
        }
        check(relabeled)

    # The generic-Grassmann four-star is the exact missing obstruction.
    star = {edge: Q(0) for edge in edges}
    for i in range(1, 5):
        star[0, i] = Q(1, 4)
    check(star)

    verify_hodge_contraction()
    print("unrestricted fermionic--Motzkin identities verified exactly")


def verify_hodge_contraction():
    """Verify the four-dimensional identity used in the invariant form."""
    edges = list(combinations(range(4), 2))
    edge_index = {edge: index for index, edge in enumerate(edges)}

    def oriented_edge(i, j):
        if i < j:
            return edge_index[i, j], 1
        if j < i:
            return edge_index[j, i], -1
        return None, 0

    hodge = sp.zeros(6)
    for column, (i, j) in enumerate(edges):
        for row, (k, ell) in enumerate(edges):
            if len({i, j, k, ell}) == 4:
                permutation = [i, j, k, ell]
                inversions = sum(
                    permutation[a] > permutation[b]
                    for a in range(4)
                    for b in range(a + 1, 4)
                )
                hodge[row, column] = (-1) ** inversions
    assert hodge * hodge == sp.eye(6)

    variables = sp.symbols("b0:21")
    matrix = sp.zeros(6)
    position = 0
    for i in range(6):
        for j in range(i, 6):
            matrix[i, j] = matrix[j, i] = variables[position]
            position += 1

    contraction = sp.zeros(4)
    for i in range(4):
        for k in range(4):
            for j in range(4):
                left, left_sign = oriented_edge(i, j)
                right, right_sign = oriented_edge(k, j)
                if left_sign and right_sign:
                    contraction[i, k] += left_sign * right_sign * matrix[left, right]

    claimed = sp.trace(contraction * contraction)
    decomposed = (
        sp.trace(matrix * matrix)
        + sp.trace(matrix) ** 2
        - sp.trace(matrix * hodge * matrix * hodge)
    )
    assert sp.expand(claimed - decomposed) == 0


if __name__ == "__main__":
    main()
