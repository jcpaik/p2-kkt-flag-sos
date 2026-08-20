import torch

torch.set_default_dtype(torch.float64)
batch, atoms = 256, 8
raw_x = torch.randn(batch, atoms, 3, requires_grad=True)
raw_w = torch.randn(batch, atoms, requires_grad=True)
opt = torch.optim.Adam([raw_x, raw_w], lr=0.02)


def calc():
    x = raw_x / torch.linalg.vector_norm(raw_x, dim=2, keepdim=True)
    w = torch.softmax(raw_w, dim=1)
    gram = torch.einsum("bik,bjk->bij", x, x)
    q = gram.square()
    E = torch.einsum(
        "bi,bij,bj->b", w, 32 * q**3 - 48 * q**2 + 20 * q - 4 / 3, w
    )
    seed1 = torch.tensor([1.0, 0.0, 0.0]).expand_as(x)
    seed2 = torch.tensor([0.0, 1.0, 0.0]).expand_as(x)
    seed = torch.where((x[:, :, :1].abs() < 0.8), seed1, seed2)
    tangent1 = torch.linalg.cross(x, seed)
    tangent1 = tangent1 / torch.linalg.vector_norm(tangent1, dim=2, keepdim=True)
    tangent2 = torch.linalg.cross(x, tangent1)
    t = gram
    z = torch.einsum("bic,bjc->bij", tangent1, x) + 1j * torch.einsum(
        "bic,bjc->bij", tangent2, x
    )
    A = torch.einsum("bj,bij->bi", w, t.square() * z.abs() ** 4)
    B = torch.einsum("bj,bij->bi", w.to(torch.complex128), t.square() * z**4)
    D = (A.square() - B.abs().square()).clamp_min(0)
    C = torch.einsum("bi,bi->b", w, A - B.abs())
    F = torch.einsum("bi,bi->b", w, D)
    return E, C, F, x, w


for step in range(3001):
    opt.zero_grad()
    E, C, F, *_ = calc()
    gap = E - 4 * C
    gap.sum().backward()
    opt.step()
    if step % 500 == 0:
        i = torch.argmin(gap)
        print(step, gap[i].item(), E[i].item(), C[i].item(), E[i].item() / F[i].item() if F[i] > 1e-14 else 0)

E, C, F, x, w = calc()
i = torch.argmin(E - 4 * C)
print("best", (E[i] - 4 * C[i]).item(), E[i].item(), C[i].item(), F[i].item())
print(w[i].detach().numpy())
print(x[i].detach().numpy())
