"""Least-norm off-diagonal bound from the true tangent linear span.

Diagnostic for the Motzkin-star route.  In an arbitrary orbital basis of
H2, write a fermionic state as diagonal edge weights plus an off-diagonal
matrix.  A genuine tangent moment lies in the 28-dimensional SO(3) orbit
span and its marginal is diagonal in an F-eigenbasis.  These linear facts
force a least possible off-diagonal Hilbert--Schmidt norm, computed here as
an explicit quadratic form in the ten edge weights.
"""

import sys

sys.path.insert(0, "research")

import numpy as np
from scipy.linalg import null_space

import tensor_fermionic_general_relaxation_opt as rel
import tensor_fermionic_relaxation_opt as base


PAIRS = rel.wedge.pairs


def wedge_matrix(R):
    """Matrix of wedge^2 R, with columns in the transformed wedge basis."""
    U = np.zeros((10, 10))
    for a, (i, j) in enumerate(PAIRS):
        for b, (k, ell) in enumerate(PAIRS):
            U[a, b] = R[i, k] * R[j, ell] - R[i, ell] * R[j, k]
    return U


def tangent_span(samples=120, seed=71):
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(samples):
        x = rng.normal(size=3)
        x /= np.linalg.norm(x)
        block = rel.tangent_data(x)[1]
        physical = rel.T @ block @ rel.T.T
        rows.append(physical.reshape(-1))
    rows = np.stack(rows)
    _, singular, vh = np.linalg.svd(rows, full_matrices=True)
    rank = int(np.sum(singular > 1e-9))
    # Constraint matrices are orthonormal in the full (not symmetric-
    # coordinate) Frobenius metric.
    return vh[rank:].reshape(-1, 10, 10), singular, rank


SPAN_NULL, SPAN_SINGULAR, SPAN_RANK = tangent_span()


def offdiag_basis():
    out = []
    for i in range(10):
        for j in range(i + 1, 10):
            E = np.zeros((10, 10))
            E[i, j] = E[j, i] = 1 / np.sqrt(2)
            out.append(E)
    return out


OFF = offdiag_basis()
DIAG = [np.diag(np.eye(10)[i]) for i in range(10)]


def least_norm_matrix(R):
    """Return the forced least-norm form and any compatibility equations.

    If ``compatibility @ p == 0``, the least possible off-diagonal norm is
    exactly ``p.T @ M @ p``.  Special orbital bases can make compatibility
    nonempty; dropping those equations would produce a spurious bound.
    """
    U = wedge_matrix(R)

    def physical(X):
        return U @ X @ U.T

    # True-span equations.
    A_span = np.array(
        [[np.sum(N * physical(E)) for E in OFF] for N in SPAN_NULL]
    )
    B_span = np.array(
        [[np.sum(N * physical(E)) for E in DIAG] for N in SPAN_NULL]
    )

    # The chosen orbital basis diagonalizes F.  D already has diagonal
    # contraction, so the off-diagonal contraction must vanish.
    marginal_rows = []
    for i in range(5):
        for j in range(i + 1, 5):
            marginal_rows.append((i, j))
    A_marginal = []
    for i, j in marginal_rows:
        row = []
        for E in OFF:
            Fp = (base.GAMMA.numpy() @ physical(E).reshape(-1)).reshape(5, 5)
            Fn = R.T @ Fp @ R
            row.append(np.sqrt(2) * Fn[i, j])
        A_marginal.append(row)
    A = np.vstack([A_span, np.array(A_marginal)])
    B = np.vstack([B_span, np.zeros((len(marginal_rows), 10))])

    # Compress dependent equations, then use the minimum-norm right inverse.
    u, singular, _ = np.linalg.svd(A, full_matrices=True)
    rank = int(np.sum(singular > 1e-9))
    Ar = u[:, :rank].T @ A
    Br = u[:, :rank].T @ B
    right = Ar.T @ np.linalg.inv(Ar @ Ar.T)
    K = -right @ Br
    compatibility_raw = u[:, rank:].T @ B
    _, csingular, cvh = np.linalg.svd(compatibility_raw, full_matrices=False)
    crank = int(np.sum(csingular > 1e-9))
    compatibility = cvh[:crank]
    residual = np.linalg.norm((A @ K + B) - u[:, rank:] @ compatibility_raw)
    return (K.T @ K + K.T @ K) / 2, compatibility, residual, rank, K


def graph_remainder_matrix():
    """Quadratic matrix for 2 Db + ||u0||^2/2 - Delta^2/24."""
    def value(p):
        a = np.array([p[PAIRS.index((0, i))] for i in range(1, 5)])
        b = {(i, j): p[PAIRS.index((i, j))] for i in range(1, 5) for j in range(i + 1, 5)}
        r = np.array([sum(b[tuple(sorted((i, j)))] for j in range(1, 5) if j != i) for i in range(1, 5)])
        Db = b[(1, 2)] * b[(3, 4)] + b[(1, 3)] * b[(2, 4)] + b[(1, 4)] * b[(2, 3)]
        c = np.sum(a)
        delta = 3 * c - 2 * np.sum(p)
        u = a - r
        u0 = u - delta * np.ones(4) / 4
        return 2 * Db + np.sum(u0 * u0) / 2 - delta * delta / 24

    Q = np.zeros((10, 10))
    eye = np.eye(10)
    diagonal = [value(eye[i]) for i in range(10)]
    for i in range(10):
        Q[i, i] = diagonal[i]
        for j in range(i + 1, 10):
            Q[i, j] = Q[j, i] = (value(eye[i] + eye[j]) - diagonal[i] - diagonal[j]) / 2
    return Q


GRAPH = graph_remainder_matrix()


def orbital_basis(theta=0.0, outer=None):
    """An H2 orbital basis with diagonal canonical first vector."""
    S = (
        np.cos(theta) * np.diag([1, -1, 0]) / np.sqrt(2)
        + np.sin(theta) * np.diag([1, 1, -2]) / np.sqrt(6)
    )
    e0 = np.array([np.sum(S * E) for E in rel.H2_MATRICES])
    complement = null_space(e0[None, :])
    if outer is not None:
        complement = complement @ outer
    return np.column_stack([e0, complement])


def report(theta=0.0, outer=None):
    R = orbital_basis(theta, outer)
    M, compatibility, residual, rank, _ = least_norm_matrix(R)
    Q = (M + GRAPH + M.T + GRAPH.T) / 2
    print("span rank/null", SPAN_RANK, len(SPAN_NULL))
    print("forced system rank/residual", rank, residual)
    print("compatibility rank/rows", len(compatibility), compatibility)
    print("M eig", np.linalg.eigvalsh(M))
    print("full quadratic eig", np.linalg.eigvalsh(Q))
    return M, Q


if __name__ == "__main__":
    report()
    rng = np.random.default_rng(8)
    raw = rng.normal(size=(4, 4))
    outer, _ = np.linalg.qr(raw)
    report(theta=0.31, outer=outer)
