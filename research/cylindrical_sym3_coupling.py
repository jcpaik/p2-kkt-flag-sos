"""Exact singular values of the Jordan contraction Sym^3(Sym^2_0 R^3)->Sym^2.

This is an auxiliary calculation for the cylindrical/Toeplitz proof search.
"""

import itertools
import sympy as sp


sqrt = sp.sqrt
z = sp.Integer(0)
E = [
    sp.diag(1, -1, 0) / sqrt(2),
    sp.diag(1, 1, -2) / sqrt(6),
    sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / sqrt(2),
    sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / sqrt(2),
    sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / sqrt(2),
]

# Fully symmetric Jordan structure constants d[a,b,c] = <e_a,e_b*e_c>.
d = sp.MutableDenseNDimArray.zeros(5, 5, 5)
for a, b, c in itertools.product(range(5), repeat=3):
    d[a, b, c] = sp.simplify(sp.trace(E[a] * (E[b] * E[c] + E[c] * E[b]) / 2))

triples = list(itertools.combinations_with_replacement(range(5), 3))
pairs = list(itertools.product(range(5), repeat=2))
A = sp.zeros(len(pairs), len(triples))

for col, triple in enumerate(triples):
    perms = sorted(set(itertools.permutations(triple)))
    coeff = 1 / sqrt(len(perms))
    full = {}
    for i, j, k in perms:
        full[i, j, k] = coeff
    for row, (i, a) in enumerate(pairs):
        value = sum(d[a, j, k] * full.get((i, j, k), z) for j in range(5) for k in range(5))
        A[row, col] = sp.simplify(value)

G = sp.simplify(A * A.T)
print("shape", A.shape)
print("rank", A.rank())
print("AA* eigenvalues")
print(G.eigenvals())
print("charpoly")
print(sp.factor(G.charpoly().as_expr()))
