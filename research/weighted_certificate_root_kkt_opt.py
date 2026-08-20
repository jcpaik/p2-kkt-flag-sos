"""Adversarial search for D_x <= E/16 at a first-order KKT measure."""

import torch

torch.set_default_dtype(torch.float64)


def kval(t): return 32*t**6-48*t**4+20*t**2-4/3
def kp(t): return 192*t**5-192*t**3+40*t


def calc(raw,logits,probes):
    x=raw/raw.norm(dim=-1,keepdim=True);w=logits.softmax(dim=-1)
    g=torch.einsum('bik,bjk->bij',x,x);K=kval(g);U=torch.einsum('bij,bj->bi',K,w);E=torch.einsum('bi,bi->b',w,U)
    # Tangential gradient of U at every support atom.
    coeff=kp(g)*w[:,None,:]
    grad=torch.einsum('bij,bjk->bik',coeff,x)
    grad=grad-torch.einsum('bik,bik->bi',grad,x)[:,:,None]*x
    # D at atom 0.
    root=x[:,0];seed1=torch.tensor([1.,0,0]).expand_as(root);seed2=torch.tensor([0.,1.,0]).expand_as(root)
    seed=torch.where((root[:,:1].abs()<.8),seed1,seed2);u=torch.linalg.cross(root,seed);u=u/u.norm(dim=-1,keepdim=True);v=torch.linalg.cross(root,u)
    t=torch.einsum('bik,bk->bi',x,root);z=torch.einsum('bik,bk->bi',x,u)+1j*torch.einsum('bik,bk->bi',x,v)
    A=torch.sum(w*t*t*z.abs()**4,dim=1);B=torch.sum(w.to(torch.complex128)*t*t*z**4,dim=1);D=A*A-B.abs()**2
    probe=probes/probes.norm(dim=-1,keepdim=True)
    gp=torch.einsum('bik,bpk->bip',x,probe);Up=torch.einsum('bi,bip->bp',w,kval(gp))
    station=torch.mean((U-E[:,None])**2,dim=1)+torch.mean(grad*grad,dim=(1,2))+torch.mean(torch.relu(E[:,None]-Up)**2,dim=1)
    return E/16-D,station,E,D,x,w,U,grad,Up


def main(B=128,N=8,P=64,steps=8000,penalty=1e5):
    gen=torch.Generator().manual_seed(314)
    raw=torch.randn(B,N,3,generator=gen,requires_grad=True);logits=torch.randn(B,N,generator=gen,requires_grad=True)
    probes=torch.randn(B,P,3,generator=gen)
    opt=torch.optim.Adam([raw,logits],lr=.015)
    for step in range(steps):
        opt.zero_grad();out=calc(raw,logits,probes);loss=(out[0]+penalty*out[1]).sum();loss.backward();opt.step()
        if step%500==0:
            score=out[0]+penalty*out[1];i=score.argmin();print(step,out[0][i].item(),out[1][i].item(),out[2][i].item(),out[3][i].item())
    out=calc(raw,logits,probes);score=out[0]+penalty*out[1];i=score.argmin();print('RESULT',out[0][i].item(),out[1][i].item(),out[2][i].item(),out[3][i].item());print(out[5][i].detach().numpy());print(out[4][i].detach().numpy())


if __name__=='__main__':main()
