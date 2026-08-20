"""Test the tempting rooted strengthening U_J(x) >= 108 D_x."""

import torch


torch.set_default_dtype(torch.float64)


def optimize(batch=384, atoms=8, steps=5000, seed=812):
    generator = torch.Generator().manual_seed(seed)
    raw_y = torch.randn(batch, atoms, 3, generator=generator, requires_grad=True)
    raw_x = torch.randn(batch, 3, generator=generator, requires_grad=True)
    raw_w = torch.randn(batch, atoms, generator=generator, requires_grad=True)
    optimizer = torch.optim.Adam([raw_y, raw_x, raw_w], lr=0.025)
    for step in range(steps + 1):
        optimizer.zero_grad()
        y = raw_y / torch.linalg.vector_norm(raw_y, dim=2, keepdim=True)
        x = raw_x / torch.linalg.vector_norm(raw_x, dim=1, keepdim=True)
        w = torch.softmax(raw_w, dim=1)
        t = torch.einsum("bik,bk->bi", y, x)
        q = t.square()
        U = torch.sum(w * (144 * q**3 - 216 * q**2 + 87 * q - 5), dim=1)
        seed1 = torch.tensor([1.0, 0.0, 0.0]).expand_as(x)
        seed2 = torch.tensor([0.0, 1.0, 0.0]).expand_as(x)
        seedv = torch.where((x[:, :1].abs() < 0.8), seed1, seed2)
        u = torch.linalg.cross(x, seedv)
        u = u / torch.linalg.vector_norm(u, dim=1, keepdim=True)
        v = torch.linalg.cross(x, u)
        z = torch.einsum("bik,bk->bi", y, u) + 1j * torch.einsum(
            "bik,bk->bi", y, v
        )
        A = torch.sum(w * q * z.abs() ** 4, dim=1)
        B = torch.sum(w.to(torch.complex128) * q * z**4, dim=1)
        D = A.square() - B.abs().square()
        gap = U - 108 * D
        gap.sum().backward()
        optimizer.step()
        if step % 500 == 0:
            i = torch.argmin(gap)
            print(step, float(gap[i]), float(U[i]), float(D[i]), flush=True)


if __name__ == "__main__":
    optimize()
