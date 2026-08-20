"""Batched adversarial atom search for W=E-(4/9)h2."""

import torch


torch.set_default_dtype(torch.float64)


def optimize(atoms=10, batch=384, steps=8000, seed=20260820):
    generator = torch.Generator().manual_seed(seed + atoms)
    raw_x = torch.randn(batch, atoms, 3, generator=generator, requires_grad=True)
    raw_w = torch.randn(batch, atoms, generator=generator, requires_grad=True)
    optimizer = torch.optim.Adam([raw_x, raw_w], lr=0.02)
    best = None
    for step in range(steps + 1):
        optimizer.zero_grad()
        x = raw_x / torch.linalg.vector_norm(raw_x, dim=2, keepdim=True)
        w = torch.softmax(raw_w, dim=1)
        gram = torch.einsum("bik,bjk->bij", x, x)
        q = gram.square()
        H = 144 * q**3 - 216 * q**2 + 87 * q - 5
        value = torch.einsum("bi,bij,bj->b", w, H, w)
        value.sum().backward()
        optimizer.step()
        if step % 500 == 0:
            i = torch.argmin(value)
            record = float(value[i]), x[i].detach().numpy(), w[i].detach().numpy()
            if best is None or record[0] < best[0]:
                best = record
            print(atoms, step, record[0], flush=True)
    return best


if __name__ == "__main__":
    for atoms in (4, 5, 6, 8, 10, 14):
        optimize(atoms=atoms)
