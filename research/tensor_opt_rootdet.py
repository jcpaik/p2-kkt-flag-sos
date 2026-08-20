import torch

torch.set_default_dtype(torch.float64)
batch, atoms = 512, 8
raw_y = torch.randn(batch, atoms, 3, requires_grad=True)
raw_x = torch.randn(batch, 3, requires_grad=True)
raw_w = torch.randn(batch, atoms, requires_grad=True)
opt = torch.optim.Adam([raw_y, raw_x, raw_w], lr=0.025)


def calc():
    y = raw_y / torch.linalg.vector_norm(raw_y, dim=2, keepdim=True)
    x = raw_x / torch.linalg.vector_norm(raw_x, dim=1, keepdim=True)
    w = torch.softmax(raw_w, dim=1)
    gram = torch.einsum("bik,bjk->bij", y, y)
    q = gram.square()
    energy = torch.einsum(
        "bi,bij,bj->b", w, 32 * q**3 - 48 * q**2 + 20 * q - 4 / 3, w
    )
    # Coordinate-free determinant: for each x choose an arbitrary tangent
    # frame. Householder-stable projected global vectors make u,v.
    seed1 = torch.tensor([1.0, 0.0, 0.0]).expand_as(x)
    seed2 = torch.tensor([0.0, 1.0, 0.0]).expand_as(x)
    seed = torch.where((x[:, :1].abs() < 0.8), seed1, seed2)
    u = torch.linalg.cross(x, seed)
    u = u / torch.linalg.vector_norm(u, dim=1, keepdim=True)
    v = torch.linalg.cross(x, u)
    t = torch.einsum("bik,bk->bi", y, x)
    z = torch.einsum("bik,bk->bi", y, u) + 1j * torch.einsum("bik,bk->bi", y, v)
    A = torch.sum(w * t.square() * z.abs() ** 4, dim=1)
    B = torch.sum(w * t.square() * z**4, dim=1)
    det = A.square() - B.abs().square()
    return energy, det, y, x, w


for step in range(5001):
    opt.zero_grad()
    E, D, *_ = calc()
    loss = E / 16 - D
    loss.sum().backward()
    opt.step()
    if step % 500 == 0:
        i = torch.argmin(loss)
        print(step, loss[i].item(), E[i].item(), D[i].item())

E, D, y, x, w = calc()
i = torch.argmin(E / 16 - D)
print("best", (E[i] / 16 - D[i]).item(), E[i].item(), D[i].item())
print(x[i].detach().numpy())
print(w[i].detach().numpy())
print(y[i].detach().numpy())
