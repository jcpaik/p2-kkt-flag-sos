"""Adversarial search for the canonical ONB-transform inequality."""

import torch

torch.set_default_dtype(torch.float64)


def kval(t):
    return 32*t**6-48*t**4+20*t**2-4/3


def objective(raw,logits,constant=1):
    x=raw/raw.norm(dim=-1,keepdim=True)
    w=logits.softmax(dim=-1)
    B,N,_=x.shape
    g=torch.einsum('bic,bjc->bij',x,x)
    pairw=w[:,:,None]*w[:,None,:]
    s2=(1-g*g).clamp_min(1e-12)
    omega=g*g*s2*s2
    m=.25*torch.sum(pairw*omega,dim=(1,2))
    E=torch.einsum('bi,bij,bj->b',w,kval(g),w)
    # Frames indexed by ordered pair and axis.
    xx=x[:,:,None,:].expand(B,N,N,3)
    yy=x[:,None,:,:].expand(B,N,N,3)
    u=(yy-g[:,:,:,None]*xx)/torch.sqrt(s2)[:,:,:,None]
    v=torch.linalg.cross(xx,u)
    frames=torch.stack((xx,u,v),dim=3).reshape(B,3*N*N,3)
    masses=(pairw*omega/12)[:,:,:,None].expand(B,N,N,3).reshape(B,3*N*N)
    gm=torch.einsum('bic,bjc->bij',frames,frames)
    Eeta=torch.einsum('bi,bij,bj->b',masses,kval(gm),masses)
    gxeta=torch.einsum('bic,bjc->bij',x,frames)
    F=torch.einsum('bi,bij,bj->b',w,kval(gxeta),masses)
    violation=Eeta-2*m*F-constant*m*m*E
    return violation,E,m,F,Eeta,x,w


def main(B=64,N=5,steps=3000,constant=4):
    gen=torch.Generator().manual_seed(55)
    raw=torch.randn(B,N,3,generator=gen,requires_grad=True)
    logits=torch.randn(B,N,generator=gen,requires_grad=True)
    opt=torch.optim.Adam([raw,logits],lr=.02)
    for step in range(steps):
        opt.zero_grad();out=objective(raw,logits,constant);loss=-out[0].sum();loss.backward();opt.step()
        if step%250==0:
            i=out[0].argmax();print(step,*[z[i].item() for z in out[:5]])
    out=objective(raw,logits,constant);i=out[0].argmax();print('RESULT',*[z[i].item() for z in out[:5]])
    print(out[-1][i].detach().numpy());print(out[-2][i].detach().numpy())


if __name__=='__main__':main()
