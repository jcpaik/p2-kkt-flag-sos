"""Minimax atom search for a counterexample to E >= min_ONB C_R."""

import argparse

import numpy as np
import torch


torch.set_default_dtype(torch.float64)


def quaternion_matrix(raw):
    q = raw / torch.linalg.vector_norm(raw, dim=-1, keepdim=True)
    w, x, y, z = q.unbind(-1)
    return torch.stack(
        [
            1 - 2 * (y*y + z*z), 2 * (x*y-z*w), 2 * (x*z+y*w),
            2 * (x*y+z*w), 1 - 2 * (x*x+z*z), 2 * (y*z-x*w),
            2 * (x*z-y*w), 2 * (y*z+x*w), 1 - 2 * (x*x+y*y),
        ], dim=-1
    ).reshape(q.shape[:-1] + (3,3))


def data(raw_x, logits, raw_q):
    x = raw_x / torch.linalg.vector_norm(raw_x, dim=1, keepdim=True)
    w = torch.softmax(logits, dim=0)
    g = x @ x.T
    w2 = w[:,None]*w[None,:]
    E = torch.sum(w2*(32*g**6-48*g**4+20*g**2-4/3))
    R = quaternion_matrix(raw_q)
    coordinates = torch.einsum("ni,kij->knj", x, R)
    C = 32*torch.sum(w[None,:]*torch.prod(coordinates**2,dim=-1),dim=-1)
    return E,C,x,w,R


def run(n=12, frames=32, starts=8, steps=6000, seed=89, temperature=300):
    rng=np.random.default_rng(seed);best=None
    for start in range(starts):
        raw_x=torch.tensor(rng.normal(size=(n,3)),requires_grad=True)
        logits=torch.tensor(rng.normal(size=n),requires_grad=True)
        raw_q=torch.tensor(rng.normal(size=(frames,4)),requires_grad=True)
        outer=torch.optim.Adam([raw_x,logits],lr=.015)
        inner=torch.optim.Adam([raw_q],lr=.025)
        for step in range(steps):
            # Several descent steps find low-C frames (maximizers of E-C).
            for _ in range(2):
                inner.zero_grad();E,C,*_=data(raw_x.detach(),logits.detach(),raw_q)
                # Keep all frame particles active, with a tiny repulsion omitted.
                C.sum().backward();inner.step()
            outer.zero_grad();E,C,*_=data(raw_x,logits,raw_q.detach())
            gaps=E-C
            loss=torch.logsumexp(temperature*gaps,dim=0)/temperature
            loss.backward();outer.step()
            if step in (3000,5000):
                for group in outer.param_groups:group['lr']*=.25
                for group in inner.param_groups:group['lr']*=.25
        E,C,x,w,R=data(raw_x,logits,raw_q)
        record=(float((E-C.min()).detach()),float(E.detach()),float(C.min().detach()),x.detach().numpy(),w.detach().numpy(),R.detach().numpy())
        if best is None or record[0]<best[0]:
            best=record;print('best',start,record[:3],'weights',record[4])
    return best


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=12);p.add_argument('--frames',type=int,default=32);p.add_argument('--starts',type=int,default=8);p.add_argument('--steps',type=int,default=6000);p.add_argument('--seed',type=int,default=89);p.add_argument('--temperature',type=float,default=300);a=p.parse_args();run(**vars(a))
