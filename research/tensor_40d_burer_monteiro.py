"""Low-rank boundary search for the 40D fermionic spectrahedron.

This uses an augmented Lagrangian in a factor X=Y Y^T.  It is intended as
an adversarial counterexample search: every reported state is checked against
the affine constraints after the optimization, and promising states can be
polished/rationalized separately.
"""

import sys

sys.path.insert(0, "research")

import numpy as np
import torch

import tensor_fermionic_l4_compound_sdp as model


torch.set_default_dtype(torch.float64)


def affine_system():
    # model.EXPANSION has orthonormal symmetric coordinates.
    trace = np.array([np.trace(E) for E in model.SYM])[None, :]
    homogeneous = model.constraint_coordinates
    raw = np.vstack([trace, homogeneous])
    rhs = np.r_[1.0, np.zeros(homogeneous.shape[0])]
    # Whiten the independent constraint rows.  This makes multiplier and
    # penalty updates insensitive to the numerical CG normalizations.
    u, singular, _ = np.linalg.svd(raw, full_matrices=False)
    rank = np.sum(singular > 1e-9)
    transform = (u[:, :rank] / singular[:rank]).T
    return transform @ raw, transform @ rhs


AFFINE, RHS = affine_system()
SYM = torch.tensor(np.stack(model.SYM))
AFFINE_T = torch.tensor(AFFINE)
RHS_T = torch.tensor(RHS)
F_LINEAR = torch.tensor(
    np.stack([model.contraction(E).reshape(-1) for E in model.SYM], axis=1)
)


def coordinates(X):
    return torch.einsum("aij,ij->a", SYM, X)


def objective(X, weighted=False):
    z = coordinates(X)
    F = (F_LINEAR @ z).reshape(5, 5)
    value = torch.sum(X * X) - torch.sum(F * F) / 2 + torch.trace(X) ** 2 / 3
    if weighted:
        Adev = X[:3, :3] - torch.eye(3) * torch.trace(X) / 15
        value = value - torch.tensor(25 / 12) * torch.sum(Adev * Adev)
    return value


def constraints(X):
    return AFFINE_T @ coordinates(X) - RHS_T


def optimize(
    rank,
    restarts=40,
    outer_steps=12,
    inner_steps=800,
    seed=20260820,
    weighted=False,
):
    generator = torch.Generator().manual_seed(seed + 101 * rank + int(weighted))
    best = None
    for restart in range(restarts):
        Y = torch.randn((10, rank), generator=generator)
        Y /= torch.linalg.norm(Y)
        Y.requires_grad_()
        multiplier = torch.zeros(len(RHS))
        penalty = 10.0
        for outer in range(outer_steps):
            optimizer = torch.optim.LBFGS(
                [Y],
                lr=0.8,
                max_iter=inner_steps,
                tolerance_grad=1e-12,
                tolerance_change=1e-14,
                line_search_fn="strong_wolfe",
            )

            def closure():
                optimizer.zero_grad()
                X = Y @ Y.T
                residual = constraints(X)
                loss = (
                    objective(X, weighted=weighted)
                    + torch.dot(multiplier, residual)
                    + penalty * torch.sum(residual * residual) / 2
                )
                loss.backward()
                return loss

            optimizer.step(closure)
            with torch.no_grad():
                X = Y @ Y.T
                residual = constraints(X)
                multiplier += penalty * residual
                violation = float(torch.linalg.norm(residual))
                if violation > 0.2 / penalty:
                    penalty *= 5
        with torch.no_grad():
            X = (Y @ Y.T).numpy()
            residual = np.linalg.norm(AFFINE @ np.array([
                np.sum(E * X) for E in model.SYM
            ]) - RHS)
            value = float(objective(torch.tensor(X), weighted=weighted))
            record = value, residual, X
            if best is None or (residual < 1e-7 and value < best[0]) or best[1] >= 1e-7:
                best = record
                print(
                    "weighted" if weighted else "plain",
                    "rank",
                    rank,
                    "restart",
                    restart,
                    "value",
                    value,
                    "residual",
                    residual,
                    "spectrum",
                    np.linalg.eigvalsh(X),
                    flush=True,
                )
    if best is not None:
        np.savez(
            f"research/tensor_40d_bm_r{rank}_{'weighted' if weighted else 'plain'}.npz",
            value=best[0],
            residual=best[1],
            X=best[2],
        )
    return best


if __name__ == "__main__":
    for weighted in (False, True):
        for rank in (2, 3, 4, 5):
            optimize(rank, weighted=weighted)
