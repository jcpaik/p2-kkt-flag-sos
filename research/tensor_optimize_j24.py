"""Batched atom search for a counterexample to E >= 24 F."""

import torch

torch.set_default_dtype(torch.float64)


def objective(raw, logits):
    x = raw / raw.norm(dim=-1, keepdim=True)
    w = logits.softmax(dim=-1)
    g = torch.einsum("bik,bjk->bij", x, x)
    p2 = torch.einsum("bi,bj,bij->b", w, w, g**2)
    p4 = torch.einsum("bi,bj,bij->b", w, w, g**4)
    p6 = torch.einsum("bi,bj,bij->b", w, w, g**6)
    E = -4/3 + 20*p2 - 48*p4 + 32*p6
    a = g[:, :, :, None]
    b = g[:, :, None, :]
    c = g[:, None, :, :]
    D = 1 + 2*a*b*c - a*a - b*b - c*c
    cyc = a*a*b*b*(c-a*b)**2 + a*a*c*c*(b-a*c)**2 + b*b*c*c*(a-b*c)**2
    h = (8/3)*D*cyc
    F = torch.einsum("bi,bj,bk,bijk->b", w, w, w, h)
    return E - 24*F, E, F, x, w


def main(B=128, N=8, steps=5000):
    gen = torch.Generator().manual_seed(101)
    raw = torch.randn(B, N, 3, generator=gen, requires_grad=True)
    logits = torch.randn(B, N, generator=gen, requires_grad=True)
    opt = torch.optim.Adam([raw, logits], lr=0.02)
    for step in range(steps):
        opt.zero_grad()
        J, E, F, _, _ = objective(raw, logits)
        J.sum().backward()
        opt.step()
        if step % 500 == 0 or step == steps-1:
            k = J.argmin().item()
            print(step, J[k].item(), E[k].item(), F[k].item())
    J,E,F,x,w=objective(raw,logits)
    for k in torch.argsort(J)[:10]:
        print("RESULT", k.item(), J[k].item(), E[k].item(), F[k].item())
        print(w[k].detach().numpy())
        print(x[k].detach().numpy())


if __name__ == "__main__":
    main()
