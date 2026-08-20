"""Search a small exact-looking scalar certificate for the 40D J lemma.

This is diagnostic: all multiplier forms below are nonnegative for a PSD
two-fermion matrix G.  A positive semidefinite residual would give a very
small certificate that can then be rationalized by hand.
"""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np

import tensor_fermionic_l4_compound_sdp as model
import tensor_wedge_relation as wedge


BASIS = model.BASIS
N = len(BASIS)
T = model.relations.T


def hole(G):
    F = model.contraction(G)
    physical = T @ G @ T.T
    Hp = np.trace(G) * np.eye(10) - wedge.additive_wedge(F) + physical
    return (T.T @ Hp @ T + T.T @ Hp.T @ T) / 2


def polarized(form):
    out = np.zeros((N, N))
    diagonal = [form(A, A) for A in BASIS]
    for i, A in enumerate(BASIS):
        out[i, i] = diagonal[i]
        for j in range(i + 1, N):
            B = BASIS[j]
            out[i, j] = out[j, i] = (
                form(A + B, A + B) - diagonal[i] - diagonal[j]
            ) / 2
    return (out + out.T) / 2


def e2(X):
    return (np.trace(X) ** 2 - np.sum(X * X)) / 2


QJ = model.Q.copy()
for i, A in enumerate(BASIS):
    A0 = A[:3, :3] - np.eye(3) * np.trace(A) / 15
    for j, B in enumerate(BASIS):
        B0 = B[:3, :3] - np.eye(3) * np.trace(B) / 15
        QJ[i, j] -= (25 / 12) * np.sum(A0 * B0)


FORMS = {
    "e2G": polarized(lambda G, _: e2(G)),
    "e2A": polarized(lambda G, _: e2(G[:3, :3])),
    "e2D": polarized(lambda G, _: e2(G[3:, 3:])),
    "crossG": polarized(
        lambda G, _: np.trace(G[:3, :3]) * np.trace(G[3:, 3:])
        - np.sum(G[:3, 3:] ** 2)
    ),
    "GH": polarized(lambda G, _: np.sum(G * hole(G))),
    "e2H": polarized(lambda G, _: e2(hole(G))),
    "e2HA": polarized(lambda G, _: e2(hole(G)[:3, :3])),
    "e2HD": polarized(lambda G, _: e2(hole(G)[3:, 3:])),
    "crossH": polarized(
        lambda G, _: np.trace(hole(G)[:3, :3])
        * np.trace(hole(G)[3:, 3:])
        - np.sum(hole(G)[:3, 3:] ** 2)
    ),
    "AAh": polarized(lambda G, _: np.sum(G[:3, :3] * hole(G)[:3, :3])),
    "DDh": polarized(lambda G, _: np.sum(G[3:, 3:] * hole(G)[3:, 3:])),
    "FL": polarized(
        lambda G, _: np.trace(G) * np.trace(model.contraction(G))
        - np.sum(model.contraction(G) ** 2)
    ),
}


def solve():
    names = list(FORMS)
    lam = cp.Variable(len(names), nonneg=True)
    residual = QJ - sum(lam[i] * FORMS[name] for i, name in enumerate(names))
    problem = cp.Problem(cp.Maximize(cp.lambda_min(residual)))
    value = problem.solve(solver="CLARABEL")
    print("status/value", problem.status, value)
    if lam.value is not None:
        print("multipliers", dict(zip(names, lam.value)))
        print("residual eig", np.linalg.eigvalsh(residual.value))


if __name__ == "__main__":
    solve()
