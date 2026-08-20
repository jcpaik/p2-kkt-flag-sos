import torch


torch.set_default_dtype(torch.float64)
device = "cpu"
batch, atoms = 256, 7
raw_x = torch.randn(batch, atoms, 3, device=device, requires_grad=True)
raw_w = torch.randn(batch, atoms, device=device, requires_grad=True)
optimizer = torch.optim.Adam([raw_x, raw_w], lr=0.025)


def values():
    x = raw_x / torch.linalg.vector_norm(raw_x, dim=2, keepdim=True)
    w = torch.softmax(raw_w, dim=1)
    gram = torch.einsum("bik,bjk->bij", x, x)
    q = gram.square()
    kernel = 32 * q**3 - 48 * q**2 + 20 * q - 4 / 3
    energy = torch.einsum("bi,bij,bj->b", w, kernel, w)

    a = gram[:, :, :, None]
    b = gram[:, :, None, :]
    c = gram[:, None, :, :]
    determinant = 1 + 2 * a * b * c - a.square() - b.square() - c.square()
    integrand = a.square() * b.square() * (c - a * b).square() * determinant
    f = 8 * torch.einsum("bi,bj,bk,bijk->b", w, w, w, integrand)
    return energy, f, x, w


for step in range(5001):
    optimizer.zero_grad()
    energy, f, x, w = values()
    gap = energy - 24 * f
    gap.sum().backward()
    optimizer.step()
    if step % 500 == 0:
        i = torch.argmin(gap)
        print(step, gap[i].item(), energy[i].item(), f[i].item())

energy, f, x, w = values()
i = torch.argmin(energy - 24 * f)
print("best", (energy[i] - 24 * f[i]).item(), energy[i].item(), f[i].item())
print(w[i].detach().numpy())
print(x[i].detach().numpy())
