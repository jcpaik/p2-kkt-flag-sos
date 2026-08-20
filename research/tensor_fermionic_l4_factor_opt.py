"""Nonconvex factor audit of the spin-2/spin-4 PSD relaxation."""

import sys

sys.path.insert(0, "research")

import numpy as np
import torch

import tensor_fermionic_general_relaxation_opt as rel


torch.set_default_dtype(torch.float64)

F_LINEAR = torch.tensor(rel.F_LINEAR)
SPIN2 = torch.tensor(rel.LINEAR_CONSTRAINT_REDUCED)
L4 = torch.tensor(rel.L4_REDUCED)


def values(raw, penalty):
    X0 = raw @ raw.T
    X = X0 / torch.trace(X0)
    x = X.reshape(-1)
    F = (F_LINEAR @ x).reshape(5, 5)
    gap = torch.sum(X * X) - torch.sum(F * F) / 2 + torch.tensor(1 / 3)
    error = (5 * torch.trace(X[:3, :3]) - 1) ** 2
    error = error + torch.sum((SPIN2 @ x) ** 2) + torch.sum((L4 @ x) ** 2)
    return gap + penalty * error, gap, error, X, F


def optimize(rank=7, restarts=30, steps=10000, penalty=2e5, seed=123):
    rng = torch.Generator().manual_seed(seed + rank)
    best = None
    for restart in range(restarts):
        raw = torch.randn((10, rank), generator=rng, requires_grad=True)
        opt = torch.optim.Adam([raw], lr=0.01)
        for _ in range(steps):
            opt.zero_grad()
            loss, *_ = values(raw, penalty)
            loss.backward()
            opt.step()
        with torch.no_grad():
            _, gap, error, X, F = values(raw, penalty)
            record = (
                float(gap),
                float(error),
                restart,
                np.linalg.eigvalsh(X.numpy()),
                np.linalg.eigvalsh(F.numpy()),
            )
            if best is None or record[0] < best[0]:
                best = record
                print("rank", rank, "best", best, flush=True)
    return best


if __name__ == "__main__":
    for rank in (3, 5, 7, 10):
        optimize(rank=rank)
