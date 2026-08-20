"""Adversarial atom search for J-108F-c times the projected-circle square."""

import argparse

import numpy as np
import torch


torch.set_default_dtype(torch.float64)


def objective(raw_x, logits, constant=1152.0, root_power=2):
    x = raw_x / torch.linalg.vector_norm(raw_x, dim=1, keepdim=True)
    w = torch.softmax(logits, dim=0)
    g = x @ x.T
    w2 = w[:, None] * w[None, :]
    p2 = torch.sum(w2 * g**2)
    p4 = torch.sum(w2 * g**4)
    p6 = torch.sum(w2 * g**6)
    J = 144 * p6 - 216 * p4 + 87 * p2 - 5

    # Symmetric determinant certificate F.
    a = g[:, :, None]
    b = g[:, None, :]
    c = g[None, :, :]
    D = 1 + 2 * a * b * c - a * a - b * b - c * c
    cyc = (
        a * a * b * b * (c - a * b) ** 2
        + a * a * c * c * (b - a * c) ** 2
        + b * b * c * c * (a - b * c) ** 2
    )
    w3 = w[:, None, None] * w[None, :, None] * w[None, None, :]
    F = torch.sum(w3 * (8 / 3) * D * cyc)

    aa = g[:, :, None, None]
    u = g[:, None, :, None]
    v = g[None, :, :, None]
    r = g[:, None, None, :]
    t = g[None, :, None, :]
    U = u * u + v * v - 2 * aa * u * v
    V = r * r + t * t - 2 * aa * r * t
    L = u * r + v * t - aa * (u * t + v * r)
    UV = U * V
    phi = aa ** (2 * root_power) * (
        32 * L**6 - 48 * L**4 * UV + 20 * L**2 * UV**2 - 2 * UV**3
    )
    w4 = (
        w[:, None, None, None]
        * w[None, :, None, None]
        * w[None, None, :, None]
        * w[None, None, None, :]
    )
    A = torch.sum(w4 * phi)
    return J - 108 * F - constant * A, (J, F, A), x, w


def run(n=8, starts=32, steps=5000, seed=123, constant=1152, root_power=2):
    rng = np.random.default_rng(seed)
    best = None
    for start in range(starts):
        raw = torch.tensor(rng.normal(size=(n, 3)), requires_grad=True)
        logits = torch.tensor(rng.normal(size=n), requires_grad=True)
        opt = torch.optim.Adam([raw, logits], lr=0.025)
        for step in range(steps):
            opt.zero_grad()
            value, _, _, _ = objective(raw, logits, constant, root_power)
            value.backward()
            opt.step()
            if step in (2500, 4000):
                for group in opt.param_groups:
                    group["lr"] *= 0.25
        value, data, x, w = objective(raw, logits, constant, root_power)
        record = (
            float(value.detach()),
            tuple(float(v.detach()) for v in data),
            x.detach().numpy(),
            w.detach().numpy(),
        )
        if best is None or record[0] < best[0]:
            best = record
            print("best", start, record[:2], "weights", record[3])
    return best


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--starts", type=int, default=16)
    parser.add_argument("--steps", type=int, default=4000)
    parser.add_argument("--constant", type=float, default=1152)
    parser.add_argument("--root-power", type=int, default=2)
    args = parser.parse_args()
    run(**vars(args))
