"""Optimize the fermionic purity gap over the fixed-H1-block PSD relaxation.

We work in Lambda^2(H2)=H1+H3 and impose only

    G = [[I3/15,C],[C^T,D]] >= 0,  tr(D)=4/5.

The one-particle contraction F=2 Tr_2 G is automatic.  The target is

    Q(G)=||G||^2 - ||F||^2/2 + 1/3 >= 0.

Schur-parameterize D=15 C^T C+R with R>=0 and normalize the common lower
scale so tr(D)=4/5.  Adam then searches the entire interior of the relaxation.
"""

import sys

sys.path.insert(0, "research")

import numpy as np
import torch

import tensor_wedge_relation as wedge


torch.set_default_dtype(torch.float64)

T = torch.tensor(np.hstack([wedge.W1, wedge.W3]))
pairs = wedge.pairs
GAMMA = np.zeros((25, 100))
for a, (i, j) in enumerate(pairs):
    for b, (k, ell) in enumerate(pairs):
        col = 10 * a + b
        GAMMA[5 * i + k, col] += j == ell
        GAMMA[5 * i + ell, col] -= j == k
        GAMMA[5 * j + k, col] -= i == ell
        GAMMA[5 * j + ell, col] += i == k
GAMMA = torch.tensor(GAMMA)

# Orthogonal projection onto the spin-2 operator sector in End(H2).
def casimir_h2(X):
    return -sum(
        torch.tensor(g) @ (torch.tensor(g) @ X - X @ torch.tensor(g))
        - (torch.tensor(g) @ X - X @ torch.tensor(g)) @ torch.tensor(g)
        for g in wedge.gens2
    )


def spin2_part(X):
    out = X
    for ell in (0, 4):
        out = (casimir_h2(out) - ell * (ell + 1) * out) / (6 - ell * (ell + 1))
    return out


def feasible(raw_c, raw_l):
    total = 15 * torch.sum(raw_c**2) + torch.sum(raw_l**2)
    scale = torch.sqrt(torch.tensor(4 / 5) / total)
    C = scale * raw_c
    Rfactor = scale * raw_l
    D = 15 * C.T @ C + Rfactor @ Rfactor.T
    top = torch.eye(3) / 15
    block = torch.cat(
        [torch.cat([top, C], dim=1), torch.cat([C.T, D], dim=1)], dim=0
    )
    G = T @ block @ T.T
    return G, block


def gap(raw_c, raw_l):
    G, block = feasible(raw_c, raw_l)
    F = (GAMMA @ G.reshape(-1)).reshape(5, 5)
    return torch.sum(G**2) - torch.sum(F**2) / 2 + torch.tensor(1 / 3), G, F, block


def optimize(restarts=100, steps=5000, learning_rate=0.015, seed=20260820, spin2_penalty=0.0):
    generator = torch.Generator().manual_seed(seed)
    best = None
    for restart in range(restarts):
        C = torch.randn((3, 7), generator=generator, requires_grad=True)
        L = torch.randn((7, 7), generator=generator, requires_grad=True)
        optimizer = torch.optim.Adam([C, L], lr=learning_rate)
        for step in range(steps):
            optimizer.zero_grad()
            value, _, F, _ = gap(C, L)
            loss = value + spin2_penalty * torch.sum(spin2_part(F) ** 2)
            loss.backward()
            optimizer.step()
        with torch.no_grad():
            value, G, F, block = gap(C, L)
            record = (
                float(value),
                float(torch.linalg.norm(spin2_part(F))),
                restart,
                np.linalg.eigvalsh(G.numpy()),
                np.linalg.eigvalsh(F.numpy()),
                np.linalg.eigvalsh(block[3:, 3:].numpy()),
                C.detach().numpy(),
                L.detach().numpy(),
            )
            if best is None or record[0] < best[0]:
                best = record
                print("best", best[:-2], flush=True)
    return best


if __name__ == "__main__":
    optimize()
