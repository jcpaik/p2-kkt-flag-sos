"""Batched Adam stress test for G_c=J-108F-c A4."""

import argparse
import numpy as np
import torch

torch.set_default_dtype(torch.float64)


def objective(raw, logits, constant=288.0, frame_constant=0.0):
    x=raw/torch.linalg.vector_norm(raw,dim=-1,keepdim=True)
    n=x.shape[1]
    w=torch.softmax(logits,dim=-1)
    g=torch.einsum('bni,bmi->bnm',x,x)
    w2=w[:,:,None]*w[:,None,:]
    p2=torch.sum(w2*g**2,dim=(1,2));p4=torch.sum(w2*g**4,dim=(1,2));p6=torch.sum(w2*g**6,dim=(1,2))
    J=144*p6-216*p4+87*p2-5
    a=g[:,:, :,None];b=g[:,:,None,:];c=g[:,None,:,:]
    D=1+2*a*b*c-a*a-b*b-c*c
    cyc=a*a*b*b*(c-a*b)**2+a*a*c*c*(b-a*c)**2+b*b*c*c*(a-b*c)**2
    w3=w[:,:,None,None]*w[:,None,:,None]*w[:,None,None,:]
    F=torch.sum(w3*(8/3)*D*cyc,dim=(1,2,3))
    aa=g[:,:, :,None,None];u=g[:,:,None,:,None];v=g[:,None,:,:,None];r=g[:,:,None,None,:];t=g[:,None,:,None,:]
    U=u*u+v*v-2*aa*u*v;V=r*r+t*t-2*aa*r*t;L=u*r+v*t-aa*(u*t+v*r);UV=U*V
    phi=aa**4*(32*L**6-48*L**4*UV+20*L**2*UV**2-2*UV**3)
    w4=w[:,:,None,None,None]*w[:,None,:,None,None]*w[:,None,None,:,None]*w[:,None,None,None,:]
    A=torch.sum(w4*phi,dim=(1,2,3,4))
    # T_mu has 3*n*n atoms of weights 8*w_i*w_j*[a^2(1-a^2)^2/4].
    # Frames at coincident pairs have zero transform weight; clamp only the
    # normalization denominator, not the polynomial weight.
    s=torch.clamp(1-g*g,min=1e-14)
    roots=x[:,:,None,:].expand(-1,-1,n,-1)
    partners=x[:,None,:,:].expand(-1,n,-1,-1)
    uu=(partners-g[:,:,:,None]*roots)/torch.sqrt(s[:,:,:,None])
    nn=torch.linalg.cross(roots,partners)/torch.sqrt(s[:,:,:,None])
    tx=torch.stack((roots,uu,nn),dim=3).reshape(raw.shape[0],3*n*n,3)
    pair_weights=2*w[:,:,None]*w[:,None,:]*g*g*s*s
    tw=pair_weights[:,:,:,None].expand(-1,-1,-1,3).reshape(raw.shape[0],3*n*n)
    tg=torch.einsum('bni,bmi->bnm',tx,tx)
    tk=144*tg**6-216*tg**4+87*tg**2-5
    JTT=torch.einsum('bi,bij,bj->b',tw,tk,tw)
    return J-108*F-constant*A-frame_constant*JTT,J,F,A,JTT,x,w


def run(n=10,batch=32,steps=8000,seed=919,constant=288.0,frame_constant=0.0):
    rng=np.random.default_rng(seed+n)
    raw=torch.tensor(rng.normal(size=(batch,n,3)),requires_grad=True)
    logits=torch.tensor(rng.normal(size=(batch,n)),requires_grad=True)
    opt=torch.optim.Adam([raw,logits],lr=.025)
    for step in range(steps):
        opt.zero_grad();gap,*_=objective(raw,logits,constant,frame_constant);gap.sum().backward();opt.step()
        if step in (steps//2,3*steps//4):
            for group in opt.param_groups:group['lr']*=.2
    values=objective(raw,logits,constant,frame_constant)
    index=int(torch.argmin(values[0]))
    print('n,batch,steps,index',n,batch,steps,index)
    print('gap,J,F,A,JTT',*[float(value[index].detach()) for value in values[:5]])
    print('weights',values[6][index].detach().numpy())
    print('points',values[5][index].detach().numpy())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=10);p.add_argument('--batch',type=int,default=32);p.add_argument('--steps',type=int,default=8000);p.add_argument('--seed',type=int,default=919);p.add_argument('--constant',type=float,default=288);p.add_argument('--frame-constant',type=float,default=0);a=p.parse_args();run(**vars(a))
