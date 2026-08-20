"""Search for an isotropic tight frame below the axial det-gap envelope."""

import argparse

import torch


torch.set_default_dtype(torch.float64)
C = 2 / 105


def axial(v):
    lower = 3 * (528 * v**2 - 126 * v - 77) / 700
    upper = -(198 * v**2 - 21 * v - 77) / 175
    u6 = torch.where(lower > 0, lower, torch.where(upper < 0, upper, 0 * v))
    determinant = (1 - 12 * v / 7) * (1 + 8 * v / 7) ** 2 * (1 - 2 * v / 7) ** 2
    return C - 24 * v**2 / 385 + 32 * u6**2 / 231 - C * determinant**2


def values(raw):
    # The columns of q are orthonormal.  Its transpose is a Parseval frame;
    # column squared norms give weights 3w_i.
    q, _ = torch.linalg.qr(raw)
    y = q.T
    norms = torch.linalg.vector_norm(y, dim=0)
    x = (y / norms).T
    w = norms**2 / 3
    t = x @ x.T
    ww = w[:, None] * w[None, :]
    p4leg = (35 * t**4 - 30 * t**2 + 3) / 8
    p6leg = (231 * t**6 - 315 * t**4 + 105 * t**2 - 5) / 16
    a4 = torch.sum(ww * p4leg)
    a6 = torch.sum(ww * p6leg)

    # Bhat quadratic form on a Frobenius-orthonormal H2 basis.
    basis = [
        torch.diag(torch.tensor([1.0, -1.0, 0.0])) / 2**0.5,
        torch.diag(torch.tensor([1.0, 1.0, -2.0])) / 6**0.5,
    ]
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        matrix = torch.zeros(3, 3)
        matrix[i, j] = matrix[j, i] = 1 / 2**0.5
        basis.append(matrix)
    basis = torch.stack(basis)
    rho = 2 * x[:, :, None] * x[:, None, :] - torch.eye(3)
    action = torch.einsum("aij,nip,bpq,nqj->nab", basis, rho, basis, rho)
    a2 = torch.einsum("n,nab->ab", w, action)
    determinant = torch.linalg.det(1.25 * (torch.eye(5) - a2))
    gap = C - 24 * a4 / 385 + 32 * a6 / 231 - C * determinant**2

    root = torch.sqrt(torch.clamp(a4, min=0))
    positive = axial(root)
    negative = axial(-root)
    # The negative axial moment is feasible only down to -7/18.
    envelope = torch.where(root <= 7 / 18, torch.minimum(positive, negative), positive)
    return gap - envelope, gap, envelope, a4, a6, determinant, x, w


def run(n, restarts, steps, seed):
    generator = torch.Generator().manual_seed(seed)
    best = None
    for restart in range(restarts):
        raw = torch.randn(n, 3, generator=generator, requires_grad=True)
        optimizer = torch.optim.Adam([raw], lr=0.02)
        for _ in range(steps):
            optimizer.zero_grad()
            result = values(raw)
            # Keep the search in the axial norm range A4 <= 49/144.
            penalty = 100 * torch.relu(result[3] - 49 / 144) ** 2
            (result[0] + penalty).backward()
            optimizer.step()
        result = values(raw)
        scalar = float(result[0].detach())
        if best is None or scalar < best[0]:
            best = (scalar, restart, [item.detach() for item in result[1:]])
            print("best", best[0], restart, [float(v) for v in best[2][:5]], flush=True)
    print("best restart", best[1])
    print("points", best[2][5].cpu().numpy())
    print("weights", best[2][6].cpu().numpy())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--restarts", type=int, default=50)
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260820)
    args = parser.parse_args()
    run(args.n, args.restarts, args.steps, args.seed)
