import torch

torch.set_default_dtype(torch.float64)
batch, atoms = 128, 8
raw_x = torch.randn(batch, atoms, 3, requires_grad=True)
raw_w = torch.randn(batch, atoms, requires_grad=True)
opt = torch.optim.Adam([raw_x, raw_w], lr=0.02)
eye = torch.eye(atoms, dtype=torch.bool)[None, :, :]


def calc():
    x = raw_x / torch.linalg.vector_norm(raw_x, dim=2, keepdim=True)
    w = torch.softmax(raw_w, dim=1)
    gram = torch.einsum("bik,bjk->bij", x, x)
    q = gram.square()
    E = torch.einsum(
        "bi,bij,bj->b", w, 32 * q**3 - 48 * q**2 + 20 * q - 4 / 3, w
    )
    cross = torch.linalg.cross(x[:, :, None, :], x[:, None, :, :])
    norm2 = cross.square().sum(dim=3)
    # Fixed distinct indices; near-parallel pairs use their limiting normal.
    normal = cross / torch.sqrt(norm2.clamp_min(1e-14))[:, :, :, None]
    qn = torch.einsum("bijc,bkc->bijk", normal, x).square()
    val = 4 * qn * (1 - qn).square()
    val = torch.where(eye[:, :, :, None], torch.zeros_like(val), val)
    R = torch.einsum("bi,bj,bk,bijk->b", w, w, w, val)
    return E, R, x, w


for step in range(601):
    opt.zero_grad()
    E, R, *_ = calc()
    gap = E - R
    gap.sum().backward()
    opt.step()
    if step % 500 == 0:
        i = torch.argmin(gap)
        print(step, gap[i].item(), E[i].item(), R[i].item())

E, R, x, w = calc()
i = torch.argmin(E - R)
print("best", (E[i] - R[i]).item(), E[i].item(), R[i].item())
print(w[i].detach().numpy())
print(x[i].detach().numpy())
