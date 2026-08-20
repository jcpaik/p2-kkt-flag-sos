"""Adversarial search for the normalized determinant-gap conjecture.

The group operator is A2[S] = E rho_x S rho_x on traceless symmetric
3-by-3 matrices and Bhat=(5/4)(I-A2), so Bhat=I for Haar measure.
We minimize

    E(mu) - (32/105) det(Bhat)^2

over finite atomic probability measures using PyTorch/Adam.
"""

import argparse

import torch


torch.set_default_dtype(torch.float64)


def sym_basis():
    mats = []
    mats.append(torch.diag(torch.tensor([1.0, -1.0, 0.0])) / 2**0.5)
    mats.append(torch.diag(torch.tensor([1.0, 1.0, -2.0])) / 6**0.5)
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        M = torch.zeros(3, 3)
        M[i, j] = M[j, i] = 1 / 2**0.5
        mats.append(M)
    return torch.stack(mats)


Q = sym_basis()


def values(raw_x, logits):
    x = raw_x / torch.linalg.vector_norm(raw_x, dim=1, keepdim=True)
    w = torch.softmax(logits, dim=0)
    t = x @ x.T
    p2 = torch.sum((w[:, None] * w[None, :]) * t**2)
    p4 = torch.sum((w[:, None] * w[None, :]) * t**4)
    p6 = torch.sum((w[:, None] * w[None, :]) * t**6)
    energy = -4 / 3 + 20 * p2 - 48 * p4 + 32 * p6

    rho = 2 * x[:, :, None] * x[:, None, :] - torch.eye(3)
    # action[n,a,b] = <Q_a, rho_n Q_b rho_n>
    action = torch.einsum("aij,nip,bpq,nqj->nab", Q, rho, Q, rho)
    A2 = torch.einsum("n,nab->ab", w, action)
    Bhat = (5 / 4) * (torch.eye(5) - A2)
    det = torch.linalg.det(Bhat)
    gap = energy - (32 / 105) * det**2
    return gap, energy, det, x, w, torch.linalg.eigvalsh(Bhat)


def run(n, restarts, steps, seed):
    generator = torch.Generator().manual_seed(seed)
    best = None
    for restart in range(restarts):
        raw_x = torch.randn(n, 3, generator=generator, requires_grad=True)
        logits = torch.randn(n, generator=generator, requires_grad=True)
        optimizer = torch.optim.Adam([raw_x, logits], lr=0.025)
        for step in range(steps):
            optimizer.zero_grad()
            gap, *_ = values(raw_x, logits)
            gap.backward()
            optimizer.step()
        result = values(raw_x, logits)
        scalar = float(result[0].detach())
        if best is None or scalar < best[0]:
            best = (scalar, restart, [v.detach().cpu() for v in result[1:]])
            energy, det, x, w, eig = best[2]
            print(
                "best",
                n,
                restart,
                scalar,
                "E",
                float(energy),
                "det",
                float(det),
                "eig",
                eig.numpy(),
                flush=True,
            )
    energy, det, x, w, eig = best[2]
    print("weights", w.numpy())
    print("points", x.numpy())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--restarts", type=int, default=100)
    parser.add_argument("--steps", type=int, default=2500)
    parser.add_argument("--seed", type=int, default=20260820)
    args = parser.parse_args()
    run(args.n, args.restarts, args.steps, args.seed)
