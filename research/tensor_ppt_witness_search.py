"""Search partial-transpose averages of the purity witness on Sym^3(R^3)^2."""

import itertools
import math

import cvxpy as cp
import numpy as np


triples = list(itertools.product(range(3), repeat=3))
triple_index = {a: i for i, a in enumerate(triples)}
comps = [(a, b, 3-a-b) for a in range(4) for b in range(4-a)]
S = np.zeros((27, 10))
for col, alpha in enumerate(comps):
    words = [w for w in triples if tuple(w.count(i) for i in range(3)) == alpha]
    for w in words:
        S[triple_index[w], col] = 1 / math.sqrt(len(words))

# Work in axis order A1,A2,A3,B1,B2,B3.
basis6 = list(itertools.product(range(3), repeat=6))
index6 = {a: i for i, a in enumerate(basis6)}


def local_operator(kind, particle):
    out = np.zeros((729, 729))
    for ket in basis6:
        if kind == "F":
            bra = list(ket)
            bra[particle], bra[particle+3] = bra[particle+3], bra[particle]
            out[index6[tuple(bra)], index6[ket]] = 1
        else:
            # |Omega><Omega| on Ai Bi and identity elsewhere.
            if ket[particle] != ket[particle+3]:
                continue
            for value in range(3):
                bra = list(ket)
                bra[particle] = value
                bra[particle+3] = value
                out[index6[tuple(bra)], index6[ket]] = 1
    return out


V = np.zeros((729, 100))
for a in range(10):
    for b in range(10):
        # np.kron has order A triple, then B triple.
        V[:, 10*a+b] = np.kron(S[:, a], S[:, b])

F = [local_operator("F", i) for i in range(3)]
P = [local_operator("P", i) for i in range(3)]


def witness(mask):
    X = [P[i] if mask[i] else F[i] for i in range(3)]
    raw = 2 * (X[0] @ X[1] @ X[2])
    raw -= X[0] @ X[1] + X[0] @ X[2] + X[1] @ X[2]
    raw += (X[0] + X[1] + X[2]) / 3
    return V.T @ raw @ V


Ws = []
for weight in range(4):
    mats = [witness(mask) for mask in itertools.product([0, 1], repeat=3) if sum(mask) == weight]
    Ws.append(sum(mats) / len(mats))
    print("weight", weight, "eig", np.linalg.eigvalsh(Ws[-1])[:10])

lam = cp.Variable(4, nonneg=True)
combo = sum(lam[i] * Ws[i] for i in range(4))
problem = cp.Problem(cp.Maximize(cp.lambda_min(combo)), [cp.sum(lam) == 1])
problem.solve(solver=cp.CLARABEL)
print(problem.status, problem.value, lam.value)
print(np.linalg.eigvalsh(sum(lam.value[i]*Ws[i] for i in range(4))))
