"""Numerical diagnostics for the quantitative tangent-star lemma.

For an actual finite tangent-Veronese mixture this computes every term in
Eq. (3)--(4) of ``fermionic_motzkin_unrestricted.md`` after diagonalizing
the one-particle marginal.  The search is diagnostic only; the identities
it checks are exact consequences of the fermionic construction.
"""

from __future__ import annotations

import numpy as np
import torch

import tensor_fermionic_general_relaxation_opt as rel
import tensor_fermionic_relaxation_opt as base
import tensor_wedge_relation as wedge


T = np.hstack([wedge.W1, wedge.W3])
PAIRS = wedge.pairs


def wedge_matrix(orthogonal: np.ndarray) -> np.ndarray:
    """Matrix of wedge^2(orthogonal), with columns in the new basis."""
    out = np.empty((10, 10))
    for a, (i, j) in enumerate(PAIRS):
        for b, (k, ell) in enumerate(PAIRS):
            out[a, b] = (
                orthogonal[i, k] * orthogonal[j, ell]
                - orthogonal[i, ell] * orthogonal[j, k]
            )
    return out


def mixture(xs: np.ndarray, weights: np.ndarray):
    block = sum(weight * rel.tangent_data(x)[1] for x, weight in zip(xs, weights))
    G = T @ block @ T.T
    F = (base.GAMMA.numpy() @ G.reshape(-1)).reshape(5, 5)
    return G, F, block


def star_terms(G: np.ndarray, F: np.ndarray):
    values, vectors = np.linalg.eigh(F)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    compound = wedge_matrix(vectors)
    diagonalized = compound.T @ G @ compound
    p = np.diag(diagonalized)

    c = values[0]
    delta = 3 * c - 2
    a = np.array([p[PAIRS.index((0, i))] for i in range(1, 5)])
    b = {
        (i, j): p[PAIRS.index((i, j))]
        for i in range(1, 5)
        for j in range(i + 1, 5)
    }
    r = np.array(
        [sum(value for edge, value in b.items() if i in edge) for i in range(1, 5)]
    )
    u0 = a - r - delta / 4
    disjoint = b[1, 2] * b[3, 4] + b[1, 3] * b[2, 4] + b[1, 4] * b[2, 3]
    off = np.sum(diagonalized**2) - np.sum(p**2)
    rhs = 24 * off + 48 * disjoint + 12 * np.sum(u0**2)
    q = off + 2 * disjoint + np.sum(u0**2) / 2
    purity = np.sum(G**2) - np.sum(F**2) / 2 + 1 / 3
    identity_error = purity - (q - delta**2 / 24)
    return {
        "c": c,
        "delta": delta,
        "rhs": rhs,
        "ratio": delta**2 / rhs if delta > 0 and rhs > 1e-15 else 0.0,
        "off": off,
        "disjoint": disjoint,
        "u0_squared": np.sum(u0**2),
        "purity": purity,
        "identity_error": identity_error,
        "occupations": values,
    }


def random_search(samples=10_000, seed=20260820):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(samples):
        count = int(rng.integers(1, 9))
        xs = rng.normal(size=(count, 3))
        xs /= np.linalg.norm(xs, axis=1, keepdims=True)
        weights = rng.dirichlet(np.ones(count) * 0.4)
        G, F, block = mixture(xs, weights)
        record = star_terms(G, F)
        if record["delta"] > 0 and (best is None or record["ratio"] > best[0]):
            second = block[:3, :3]
            best = (
                record["ratio"],
                record,
                xs,
                weights,
                np.sum((second - np.eye(3) / 15) ** 2),
            )
    print("best", best[:2] if best else None)
    if best:
        print("A deviation", best[4])
        print("weights", best[3])
        print("points", best[2])
    return best


if __name__ == "__main__":
    random_search()


def torch_orbit(raw_x: torch.Tensor):
    """Return Q_x and wedge^2(Q_x) for a batch of normalized roots."""
    x = raw_x / torch.linalg.norm(raw_x, dim=1, keepdim=True)
    h2 = torch.tensor(np.stack(rel.H2_MATRICES), dtype=raw_x.dtype)
    sx = torch.einsum("aij,nj->nai", h2, x)
    quadratic = torch.einsum("ni,nai->na", x, sx)
    q = 2 * (
        torch.einsum("nai,nbi->nab", sx, sx)
        - torch.einsum("na,nb->nab", quadratic, quadratic)
    )
    pair_array = PAIRS
    r = torch.empty((len(x), 10, 10), dtype=raw_x.dtype)
    for a, (i, j) in enumerate(pair_array):
        for b, (k, ell) in enumerate(pair_array):
            r[:, a, b] = q[:, i, k] * q[:, j, ell] - q[:, i, ell] * q[:, j, k]
    return q, r


def torch_star_ratio(raw_x: torch.Tensor, raw_weights: torch.Tensor):
    q, r = torch_orbit(raw_x)
    weights = torch.softmax(raw_weights, dim=0)
    F = torch.einsum("n,nij->ij", weights, q)
    G = torch.einsum("n,nij->ij", weights, r)
    values, vectors = torch.linalg.eigh(F)
    order = torch.arange(4, -1, -1)
    values = values[order]
    vectors = vectors[:, order]
    compound = torch.empty((10, 10), dtype=raw_x.dtype)
    for a, (i, j) in enumerate(PAIRS):
        for b, (k, ell) in enumerate(PAIRS):
            compound[a, b] = (
                vectors[i, k] * vectors[j, ell] - vectors[i, ell] * vectors[j, k]
            )
    ge = compound.T @ G @ compound
    p = torch.diag(ge)
    c = values[0]
    delta = 3 * c - 2
    a = torch.stack([p[PAIRS.index((0, i))] for i in range(1, 5)])
    b = {(i, j): p[PAIRS.index((i, j))] for i in range(1, 5) for j in range(i + 1, 5)}
    degrees = torch.stack(
        [sum(value for edge, value in b.items() if i in edge) for i in range(1, 5)]
    )
    u0 = a - degrees - delta / 4
    disjoint = b[1, 2] * b[3, 4] + b[1, 3] * b[2, 4] + b[1, 4] * b[2, 3]
    off = torch.sum(ge * ge) - torch.sum(p * p)
    rhs = 24 * off + 48 * disjoint + 12 * torch.sum(u0 * u0)
    return delta * delta / (rhs + 1e-14), delta, rhs, F, G


def optimize_ratio(count=4, restarts=20, steps=4000, seed=771):
    best = None
    for restart in range(restarts):
        torch.manual_seed(seed + 100 * count + restart)
        roots = torch.randn((count, 3), dtype=torch.float64, requires_grad=True)
        weights = torch.randn(count, dtype=torch.float64, requires_grad=True)
        optimizer = torch.optim.Adam([roots, weights], lr=0.015)
        for _ in range(steps):
            optimizer.zero_grad()
            ratio, delta, *_ = torch_star_ratio(roots, weights)
            # Stay in the only relevant (uncapped) region.
            loss = -ratio + 100 * torch.relu(torch.tensor(1e-4) - delta) ** 2
            loss.backward()
            optimizer.step()
        with torch.no_grad():
            ratio, delta, rhs, F, _ = torch_star_ratio(roots, weights)
            record = (
                float(ratio),
                float(delta),
                float(rhs),
                torch.softmax(weights, 0).numpy(),
                (roots / torch.linalg.norm(roots, dim=1, keepdim=True)).numpy(),
                torch.linalg.eigvalsh(F).numpy(),
            )
            if best is None or record[0] > best[0]:
                best = record
                print("ratio best", count, restart, best, flush=True)
    return best
