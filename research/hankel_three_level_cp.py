"""Test the natural CP extension of the three-level purity inequality.

The Pluecker cubic feature z_x and the symmetric cubic feature h_x are
related by an invertible linear filter A.  Thus an arbitrary G >= 0 defines
H=A^{-1} G A^{-T} >= 0, and the ordinary bosonic partial traces of H give
positive lower-level matrices.  On the Hankel linear span these coincide
with the tangent quantities L and M.  This script tests whether

    ||G||^2 - ||L||^2/2 + ||M||^2/6 + mass^2/18 >= 0

already holds on the larger filtered-bosonic PSD cone (with the homogeneous
H1/H3 trace ratio).
"""

import math
import sys

sys.path.insert(0, "research")

import numpy as np
import torch
import cvxpy as cp

import hankel_affine_fermion_full as fermion
import hankel_reduction_identities as reduction
import tensor_hankel_fw as hankel
import tensor_wedge_relation as wedge


DEG3 = hankel.compositions(3, 3)


def cubic_feature(x):
    return np.array(
        [math.sqrt(hankel.multinomial(alpha)) * np.prod(x ** np.array(alpha)) for alpha in DEG3]
    )


def fit_filter(seed=17):
    rng = np.random.default_rng(seed)
    points = rng.normal(size=(100, 3))
    points /= np.linalg.norm(points, axis=1)[:, None]
    source = np.stack([cubic_feature(x) for x in points])
    target = np.stack([fermion.tangent_z(x) for x in points])
    transpose, *_ = np.linalg.lstsq(source, target, rcond=None)
    A = transpose.T
    residual = np.linalg.norm(source @ transpose - target) / np.linalg.norm(target)
    return A, residual


A, FILTER_RESIDUAL = fit_filter()
AINV = np.linalg.inv(A)
P1 = fermion.W1 @ fermion.W1.T
P3 = np.eye(10) - P1


def levels(G):
    H = AINV @ G @ AINV.T
    R2 = reduction.reduce32(H)
    M = reduction.reduce21(R2)
    L = 2 * wedge.U2.T @ R2 @ wedge.U2
    mass = np.trace(H)
    value = (
        np.sum(G * G)
        - 0.5 * np.sum(L * L)
        + np.sum(M * M) / 6
        + mass * mass / 18
    )
    return value, H, R2, L, M, mass


def random_fixed_ratio(rng, rank=10):
    # A block-diagonal congruence can impose tr(G11):tr(G33)=1:4 while
    # preserving PSD and allowing arbitrary cross blocks.
    raw = rng.normal(size=(10, rank))
    G = raw @ raw.T
    top = np.trace(P1 @ G)
    bottom = np.trace(P3 @ G)
    scale = np.sqrt(4 * top / bottom)
    D = P1 + scale * P3
    G = D @ G @ D
    return G / np.trace(G)


def scan(samples=100000, seed=91):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(samples):
        G = random_fixed_ratio(rng, rank=int(rng.integers(1, 11)))
        record = (levels(G)[0], G)
        if best is None or record[0] < best[0]:
            best = record
    value, G = best
    result = levels(G)
    print("filter residual/determinant", FILTER_RESIDUAL, np.linalg.det(A))
    print("best", value, "eig G", np.linalg.eigvalsh(G))
    print("mass/block traces", result[-1], np.trace(P1 @ G), np.trace(P3 @ G))
    print("norms", np.sum(G * G), np.sum(result[3] ** 2), np.sum(result[4] ** 2))
    return best


torch.set_default_dtype(torch.float64)


def factor_opt(rank=10, restarts=30, steps=10000, penalty=1e4, seed=117):
    inverse = torch.tensor(AINV)
    p1 = torch.tensor(P1)
    # Partial-trace maps as explicit matrices for autograd.
    map32 = torch.zeros((36, 100))
    for i in range(10):
        for j in range(10):
            E = np.zeros((10, 10))
            E[i, j] = 1
            map32[:, 10 * i + j] = torch.tensor(reduction.reduce32(E).reshape(-1))
    map21 = torch.zeros((9, 36))
    for i in range(6):
        for j in range(6):
            E = np.zeros((6, 6))
            E[i, j] = 1
            map21[:, 6 * i + j] = torch.tensor(reduction.reduce21(E).reshape(-1))
    U2 = torch.tensor(wedge.U2)
    generator = torch.Generator().manual_seed(seed + rank)
    best = None
    for restart in range(restarts):
        raw = torch.randn((10, rank), generator=generator, requires_grad=True)
        optimizer = torch.optim.Adam([raw], lr=0.01)
        for _ in range(steps):
            optimizer.zero_grad()
            G0 = raw @ raw.T
            G = G0 / torch.trace(G0)
            H = inverse @ G @ inverse.T
            R2 = (map32 @ H.reshape(-1)).reshape(6, 6)
            M = (map21 @ R2.reshape(-1)).reshape(3, 3)
            L = 2 * U2.T @ R2 @ U2
            mass = torch.trace(H)
            value = (
                torch.sum(G * G)
                - torch.sum(L * L) / 2
                + torch.sum(M * M) / 6
                + mass * mass / 18
            )
            constraint = (5 * torch.trace(p1 @ G) - 1) ** 2
            loss = value + penalty * constraint
            loss.backward()
            optimizer.step()
        with torch.no_grad():
            G0 = raw @ raw.T
            G = G0 / torch.trace(G0)
            record = (
                float(levels(G.numpy())[0]),
                float((5 * torch.trace(p1 @ G) - 1) ** 2),
                np.linalg.eigvalsh(G.numpy()),
                G.numpy(),
            )
            if best is None or record[0] < best[0]:
                best = record
                print("factor best", restart, best[:-1], flush=True)
    return best


def linear_partial_trace_maps():
    map32 = np.zeros((36, 100))
    for i in range(10):
        for j in range(10):
            E = np.zeros((10, 10))
            E[i, j] = 1
            map32[:, 10 * i + j] = reduction.reduce32(E).reshape(-1)
    map21 = np.zeros((9, 36))
    for i in range(6):
        for j in range(6):
            E = np.zeros((6, 6))
            E[i, j] = 1
            map21[:, 6 * i + j] = reduction.reduce21(E).reshape(-1)
    return map32, map21


MAP32, MAP21 = linear_partial_trace_maps()


def level_maps_from_g():
    """Explicit C-order maps G -> (L, M, mass).

    Keeping these as single affine maps also avoids CVXPY's faulty trace and
    matrix-product canonicalizations in version 1.8.1.
    """
    map_l = np.zeros((25, 100))
    map_m = np.zeros((9, 100))
    map_mass = np.zeros(100)
    for i in range(10):
        for j in range(10):
            E = np.zeros((10, 10))
            E[i, j] = 1
            _, H, _, L, M, mass = levels(E)
            column = 10 * i + j
            map_l[:, column] = L.reshape(-1)
            map_m[:, column] = M.reshape(-1)
            map_mass[column] = mass
    return map_l, map_m, map_mass


MAP_GL, MAP_GM, MAP_GMASS = level_maps_from_g()


def cccp(initial, iterations=50, verbose=False):
    G = cp.Variable((10, 10), symmetric=True)
    gvec = cp.reshape(G, (100,), order="C")
    lvec = MAP_GL @ gvec
    mvec = MAP_GM @ gvec
    mass = MAP_GMASS @ gvec
    constraints = [
        G >> 0,
        np.eye(10).reshape(-1) @ gvec == 1,
        5 * P1.reshape(-1) @ gvec == 1,
    ]
    current = initial.copy()
    for iteration in range(iterations):
        l0 = levels(current)[3].reshape(-1)
        # Linearization of -||L||^2/2 has gradient -L0.
        objective = cp.Minimize(
            cp.sum_squares(gvec)
            + cp.sum_squares(mvec) / 6
            + cp.square(mass) / 18
            - l0 @ lvec
        )
        problem = cp.Problem(objective, constraints)
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=1e-10,
            tol_feas=1e-10,
            tol_gap_rel=1e-10,
            max_iter=1000,
        )
        current = G.value
        if verbose:
            print("cccp", iteration, levels(current)[0], np.linalg.eigvalsh(current))
    return levels(current)[0], current


def cccp_restarts(restarts=30, iterations=50, seed=331):
    rng = np.random.default_rng(seed)
    best = None
    for restart in range(restarts):
        initial = random_fixed_ratio(rng, rank=int(rng.integers(1, 11)))
        result = cccp(initial, iterations=iterations)
        if best is None or result[0] < best[0]:
            best = result
            print(
                "cccp best",
                restart,
                best[0],
                np.linalg.eigvalsh(best[1]),
                flush=True,
            )
    return best


if __name__ == "__main__":
    scan()
