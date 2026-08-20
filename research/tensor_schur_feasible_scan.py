"""Numerically audit the BB* pairing bound on the *feasible* isotropic Schur cone.

For a unit spin-4 direction h (normalized so ||B(h)||_HS=1), determine by
bisection the largest b for which there is a spin-6 completion A6 with

    (2/35) I + sqrt(b) A4(h) + A6 - 5 b B(h)B(h)^T >= 0.

Then compare the exact scalar lower bound obtained by pairing this Schur
complement with BB* against the desired energy inequality.
"""

import itertools
import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

from tensor_hankel_fw import maps


words = list(itertools.product(range(3), repeat=3))
word_index = {w: i for i, w in enumerate(words)}
comps = [(a, b, 3 - a - b) for a in range(4) for b in range(4 - a)]
S = np.zeros((27, 10))
for j, alpha in enumerate(comps):
    matches = [w for w in words if tuple(w.count(i) for i in range(3)) == alpha]
    for w in matches:
        S[word_index[w], j] = 1 / np.sqrt(len(matches))

eps = np.zeros((3, 3, 3))
eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1
gens = []
for axis in range(3):
    L = np.array([[eps[i, axis, j] for j in range(3)] for i in range(3)])
    G = (
        np.kron(np.kron(L, np.eye(3)), np.eye(3))
        + np.kron(np.kron(np.eye(3), L), np.eye(3))
        + np.kron(np.kron(np.eye(3), np.eye(3)), L)
    )
    gens.append(S.T @ G @ S)

J = np.zeros((10, 3))
for a in range(3):
    T = np.zeros((3, 3, 3))
    for i, j, k in words:
        T[i, j, k] = (
            (i == j and k == a) + (i == k and j == a) + (j == k and i == a)
        ) / np.sqrt(15)
    J[:, a] = S.T @ T.reshape(-1)
U = null_space(J.T)
g3 = [U.T @ g @ U for g in gens]


def cas_full(X):
    return -sum(g @ (g @ X - X @ g) - (g @ X - X @ g) @ g for g in gens)


def proj_full(X, ell):
    out = X.copy()
    for j in [0, 2, 4, 6]:
        if j != ell:
            out = (cas_full(out) - j * (j + 1) * out) / (
                ell * (ell + 1) - j * (j + 1)
            )
    return out


def cas3(X):
    return -sum(g @ (g @ X - X @ g) - (g @ X - X @ g) @ g for g in g3)


def proj3(X, ell):
    out = X.copy()
    for j in [0, 2, 4, 6]:
        if j != ell:
            out = (cas3(out) - j * (j + 1) * out) / (
                ell * (ell + 1) - j * (j + 1)
            )
    return out


raw4 = [proj_full(maps[3][:, :, k], 4) for k in range(28)]
V4 = np.stack([x.reshape(-1) for x in raw4], axis=1)
u4, s4, vh4 = np.linalg.svd(V4, full_matrices=False)
basis4 = [u4[:, i].reshape(10, 10) for i in range(9)]
# Pure moment-coordinate representatives for the same orthonormal H4 matrices.
# Inverting the *full* Hankel map avoids contaminating a spin-4 matrix by
# domain components annihilated by the spin-4 projection.
Vfull = np.stack([maps[3][:, :, k].reshape(-1) for k in range(28)], axis=1)
Vfull_pinv = np.linalg.pinv(Vfull)
ybasis4 = [Vfull_pinv @ basis4[i].reshape(-1) for i in range(9)]

comps2 = [(a, b, 2 - a - b) for a in range(3) for b in range(3 - a)]
j2 = np.array(
    [1 / np.sqrt(3) if max(alpha) == 2 else 0.0 for alpha in comps2]
)[:, None]
U2 = null_space(j2.T)
C4basis = [
    U2.T @ np.einsum("abk,k->ab", maps[2], y) @ U2 for y in ybasis4
]

raw6 = [proj_full(maps[3][:, :, k], 6) for k in range(28)]
V6 = np.stack([x.reshape(-1) for x in raw6], axis=1)
u6, s6, _ = np.linalg.svd(V6, full_matrices=False)
basis6 = [U.T @ u6[:, i].reshape(10, 10) @ U for i in range(13)]


def direction(coeff):
    H4 = sum(c * M for c, M in zip(coeff, basis4))
    A4 = U.T @ H4 @ U
    B = U.T @ H4 @ J
    scale = np.sqrt(np.sum(B * B))
    A4 /= scale
    B /= scale
    C4 = sum(c * M for c, M in zip(coeff, C4basis)) / scale
    Q = B @ B.T
    return A4, B, C4, (
        np.sum(Q * Q),
        np.sum(A4 * proj3(Q, 4)),
        np.sum(proj3(Q, 6) ** 2),
    )


def pairing_gap(b, inv):
    q, c, d = inv
    rhs = 5 * b * b * q - (2 / 35) * b - b ** 1.5 * c
    lower_q = max(0.0, rhs) ** 2 / (b * b * d) if d > 1e-14 else 0.0

    # A second necessary consequence of R >= 0 is ||R||^2 <= tr(R)^2.
    # With y=||A6|| and |<A6,Q6>| <= y ||Q6|| this gives
    # y^2 - 10 y ||Q6|| <= D.  If D<0 its smaller positive root is
    # another rigorous lower bound for y.
    q_actual = b * b * q
    c_actual = b ** 1.5 * c
    d_actual = b * b * d
    trace_r = 2 / 5 - 5 * b
    D = (
        trace_r * trace_r
        - (4 / 175 + 3 * b / 11)
        + 10 * ((2 / 35) * b + c_actual)
        - 25 * q_actual
    )
    lower_trace = 0.0
    disc = 25 * d_actual + D
    if D < 0 and disc >= 0:
        lower_trace = (5 * np.sqrt(d_actual) - np.sqrt(disc)) ** 2
    lower = max(lower_q, lower_trace)
    return 2 / 105 - 5 * b / 11 + 2 * lower


def audit(coeff):
    A4, B, C4, inv = direction(coeff)
    Q = B @ B.T
    z = cp.Variable(13)
    A6 = sum(z[i] * basis6[i] for i in range(13))
    t = cp.Parameter(nonneg=True)
    R = (2 / 35) * np.eye(7) + t * A4 + A6 - 5 * t * t * Q
    problem = cp.Problem(cp.Minimize(cp.sum_squares(A6)), [R >> 0])

    def solve(b):
        t.value = np.sqrt(b)
        try:
            value = problem.solve(solver="CLARABEL", warm_start=True)
        except cp.error.SolverError:
            value = problem.solve(
                solver="SCS", eps=2e-7, max_iters=100000, warm_start=True
            )
        return problem.status in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE), value

    lo, hi = 0.0, 0.08000001
    for _ in range(18):
        mid = (lo + hi) / 2
        feasible, _ = solve(mid)
        if feasible:
            lo = mid
        else:
            hi = mid
    bmax = lo
    grid = np.linspace(1e-6, max(1e-6, bmax), 300)
    gaps = np.array([pairing_gap(b, inv) for b in grid])
    j = int(np.argmin(gaps))
    feasible, true_min = solve(grid[j])
    true_gap = 2 / 105 - 5 * grid[j] / 11 + 2 * true_min
    determinant = np.linalg.det(np.eye(5) - 5 * np.sqrt(grid[j]) * C4)
    det_true_gap = true_gap - (2 / 105) * determinant**2
    return bmax, grid[j], gaps[j], true_gap, det_true_gap, determinant, inv


def main():
    rng = np.random.default_rng(20260819)
    worst = None
    for i in range(160):
        result = audit(rng.normal(size=9))
        if worst is None or result[2] < worst[0]:
            worst = (result[2], i, result)
            print("worst", worst, flush=True)


if __name__ == "__main__":
    main()
