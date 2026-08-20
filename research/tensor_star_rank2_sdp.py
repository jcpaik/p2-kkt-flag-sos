"""Quadratic certificate search for the rank-two top-orbital chart.

Fix S=diag(1,-1,0)/sqrt(2), impose only F S = c S on the exact
28-dimensional tangent span, and retain the top slack
K=c I-F|S^perp.  Search

    P = square + <Wg,C2(G)> + <Wk,C2(K)> + <Wgk,G tensor K>,

where every cone multiplier is PSD.  Feasibility is a diagnostic until a
solution is rationalized exactly.
"""

import itertools
import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
import sympy as sp

import tensor_star_canonical_exact as exact


def data():
    coordinates, span = exact.tangent_coefficient_span()
    Gmap = sp.zeros(100, 28)
    for row, (a, b) in enumerate(coordinates):
        for k in range(28):
            value = span[row, k] if a == b else span[row, k] / sp.sqrt(2)
            Gmap[10 * a + b, k] = Gmap[10 * b + a, k] = value

    pairs = exact.PAIRS
    edge = {pair: i for i, pair in enumerate(pairs)}

    def signed_pair(i, j):
        return ((i, j), 1) if i < j else ((j, i), -1)

    Fmap = sp.zeros(25, 28)
    for i in range(5):
        for k in range(5):
            for j in range(5):
                if j == i or j == k:
                    continue
                pi, si = signed_pair(i, j)
                pk, sk = signed_pair(k, j)
                Fmap[5 * i + k, :] += si * sk * Gmap[
                    10 * edge[pi] + edge[pk], :
                ]

    eigen = sp.Matrix([Fmap[i, :] for i in range(1, 5)])
    null = sp.Matrix.hstack(*eigen.nullspace())
    Gmap = Gmap * null
    Fmap = Fmap * null
    dimension = null.shape[1]
    Gs = [
        np.array(sp.Matrix(10, 10, Gmap[:, i]), dtype=float)
        for i in range(dimension)
    ]
    Fs = [
        np.array(sp.Matrix(5, 5, Fmap[:, i]), dtype=float)
        for i in range(dimension)
    ]
    Ks = [F[0, 0] * np.eye(4) - F[1:, 1:] for F in Fs]
    return Gs, Fs, Ks


GS, FS, KS = data()
N = len(GS)


def hodge():
    J = np.zeros((6, 6))
    J[0, 5] = J[5, 0] = 1
    J[1, 4] = J[4, 1] = -1
    J[2, 3] = J[3, 2] = 1
    return J


J = hodge()
OUTER_PAIRS = list(itertools.combinations(range(4), 2))


def contraction_outer(B):
    out = np.zeros((4, 4))
    index = {pair: i for i, pair in enumerate(OUTER_PAIRS)}

    def signed(i, j):
        return ((i, j), 1) if i < j else ((j, i), -1)

    for i in range(4):
        for k in range(4):
            for j in range(4):
                if j == i or j == k:
                    continue
                pi, si = signed(i, j)
                pk, sk = signed(k, j)
                out[i, k] += si * sk * B[index[pi], index[pk]]
    return out


def target(G):
    A, C, B = G[:4, :4], G[:4, 4:], G[4:, 4:]
    R = contraction_outer(B)
    mass = np.trace(G)
    c = np.trace(A)
    delta = 3 * c - 2 * mass
    U = A - R - delta * np.eye(4) / 4
    return (
        np.sum(U * U) / 2
        + 2 * np.sum(C * C)
        + np.trace(B @ J @ B @ J)
        - delta * delta / 24
    )


def bilinear_matrix(form):
    out = np.zeros((N, N))
    diagonal = [form(A, A) for A in range(N)]
    for i in range(N):
        out[i, i] = diagonal[i]
        for j in range(i + 1, N):
            out[i, j] = out[j, i] = form(i, j)
    return out


Q = np.zeros((N, N))
for i in range(N):
    Q[i, i] = target(GS[i])
    for j in range(i + 1, N):
        Q[i, j] = Q[j, i] = (
            target(GS[i] + GS[j]) - target(GS[i]) - target(GS[j])
        ) / 2


def polarized_compound(A, B):
    pairs = list(itertools.combinations(range(A.shape[0]), 2))
    out = np.zeros((len(pairs), len(pairs)))
    for p, (i, j) in enumerate(pairs):
        for q, (k, ell) in enumerate(pairs):
            out[p, q] = (
                A[i, k] * B[j, ell]
                + B[i, k] * A[j, ell]
                - A[i, ell] * B[j, k]
                - B[i, ell] * A[j, k]
            ) / 2
    return (out + out.T) / 2


def solve(use_g=True, use_k=True, use_gk=True):
    S = cp.Variable((N, N), symmetric=True)
    Wg = cp.Variable((45, 45), symmetric=True) if use_g else None
    Wk = cp.Variable((6, 6), symmetric=True) if use_k else None
    Wgk = cp.Variable((40, 40), symmetric=True) if use_gk else None
    constraints = [S >> 0]
    objective = cp.trace(S)
    for variable in (Wg, Wk, Wgk):
        if variable is not None:
            constraints.append(variable >> 0)
            objective += cp.trace(variable)

    for i in range(N):
        for j in range(i, N):
            expression = S[i, j]
            if Wg is not None:
                expression += cp.sum(
                    cp.multiply(Wg, polarized_compound(GS[i], GS[j]))
                )
            if Wk is not None:
                expression += cp.sum(
                    cp.multiply(Wk, polarized_compound(KS[i], KS[j]))
                )
            if Wgk is not None:
                tensor = (
                    np.kron(GS[i], KS[j]) + np.kron(GS[j], KS[i])
                ) / 2
                expression += cp.sum(cp.multiply(Wgk, tensor))
            constraints.append(expression == Q[i, j])

    problem = cp.Problem(cp.Minimize(objective), constraints)
    try:
        value = problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-9,
            tol_feas=2e-9,
            tol_gap_rel=2e-9,
            max_iter=3000,
        )
    except cp.error.SolverError:
        value = problem.solve(solver="SCS", eps=3e-7, max_iters=500000)
    print("dimension", N, "target eig", np.linalg.eigvalsh(Q))
    print("status/value", problem.status, value)
    if S.value is not None:
        print("S eig", np.linalg.eigvalsh(S.value))
        for name, variable in (("Wg", Wg), ("Wk", Wk), ("Wgk", Wgk)):
            if variable is not None:
                print(name, "eig", np.linalg.eigvalsh(variable.value))


if __name__ == "__main__":
    solve()
