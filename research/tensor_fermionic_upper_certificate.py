"""Search simple quadratic certificates using the exact upper marginal LMI.

On the block-trace plus spin-4 Clebsch slice, genuine tangent moments obey

    G >= 0,   L=P(G11)-F(G) >= 0.

We test whether the target is an ordinary square plus nonnegative trace
pairings tr(F L), tr(P L), and block-derived positive pairings.
"""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

import tensor_fermionic_general_relaxation_opt as rel
import tensor_fermionic_relaxation_opt as base


def symmetric_basis(n):
    out = []
    for i in range(n):
        E = np.zeros((n, n)); E[i, i] = 1; out.append(E)
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n)); E[i, j] = E[j, i] = 1 / np.sqrt(2); out.append(E)
    return out


SYM = symmetric_basis(10)
EXPANSION = np.stack([x.reshape(-1) for x in SYM], axis=1)
trace_ratio = np.zeros((1, 100))
for i in range(3): trace_ratio[0, 10 * i + i] = 5
for i in range(10): trace_ratio[0, 10 * i + i] -= 1
constraints = np.vstack([trace_ratio, rel.L4_REDUCED]) @ EXPANSION
N = null_space(constraints)
BASIS = [sum(N[a, k] * SYM[a] for a in range(55)) for k in range(N.shape[1])]
DIM = len(BASIS)


def maps(G):
    physical = rel.T @ G @ rel.T.T
    F = (base.GAMMA.numpy() @ physical.reshape(-1)).reshape(5, 5)
    P = (rel.UPPER_NP @ G[:3, :3].reshape(-1)).reshape(5, 5)
    return F, P, P - F


def bilinear(form):
    out = np.zeros((DIM, DIM))
    diagonals = [form(A, A) for A in BASIS]
    for i, A in enumerate(BASIS):
        out[i, i] = diagonals[i]
        for j in range(i + 1, DIM):
            value = (form(A + BASIS[j], A + BASIS[j]) - diagonals[i] - diagonals[j]) / 2
            out[i, j] = out[j, i] = value
    return out


TARGET = bilinear(
    lambda A, _: np.sum(A * A)
    - 0.5 * np.sum(maps(A)[0] ** 2)
    + np.trace(A) ** 2 / 3
)
I_FL = bilinear(lambda A, _: np.sum(maps(A)[0] * maps(A)[2]))
I_PL = bilinear(lambda A, _: np.sum(maps(A)[1] * maps(A)[2]))


def scalar_search():
    lam = cp.Variable(2, nonneg=True)
    slack = TARGET - lam[0] * I_FL - lam[1] * I_PL
    problem = cp.Problem(cp.Maximize(cp.lambda_min(slack)))
    value = problem.solve(solver="CLARABEL")
    print("dim", DIM, "target eig", np.linalg.eigvalsh(TARGET))
    print("scalar", problem.status, value, lam.value)
    if lam.value is not None:
        print("slack eig", np.linalg.eigvalsh(slack.value))


if __name__ == "__main__":
    scalar_search()
