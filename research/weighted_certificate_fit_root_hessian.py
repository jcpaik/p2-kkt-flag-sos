"""Fit E/16-D_x as a linear functional of the radial spin-4 transport block."""

import numpy as np

from weighted_certificate_fermion_cap import rooted_spin4_radial_transport


def kernel(t): return 32*t**6-48*t**4+20*t**2-4/3


def gap(x,w,root=0):
    g=x@x.T;E=w@kernel(g)@w;p=x[root]
    seed=np.eye(3)[np.argmin(np.abs(p))];u=np.cross(p,seed);u/=np.linalg.norm(u);v=np.cross(p,u)
    t=x@p;z=x@u+1j*(x@v);A=np.sum(w*t*t*np.abs(z)**4);B=np.sum(w*t*t*z**4);D=A*A-abs(B)**2
    return E/16-D


def main():
    rng=np.random.default_rng(121)
    rows=[];targets=[]
    tri=np.triu_indices(6)
    for _ in range(500):
        n=5;x=rng.normal(size=(n,3));x/=np.linalg.norm(x,axis=1)[:,None];w=rng.dirichlet(np.ones(n))
        _,_,H=rooted_spin4_radial_transport(x,w,0)
        p=x[0];seed=np.eye(3)[np.argmin(np.abs(p))];u=np.cross(p,seed);u/=np.linalg.norm(u);v=np.cross(p,u)
        t=x@p;yu=x@u;yv=x@v;g=x@x.T;K=kernel(g);U=K@w;E=w@U
        q=U[0]-E
        kp=lambda s:192*s**5-192*s**3+40*s
        grad=np.array([np.sum(w*kp(t)*yu),np.sum(w*kp(t)*yv)])
        moments=[]
        for total in range(0,9):
            for i in range(total+1):
                for j in range(total-i+1):
                    k=total-i-j
                    if total%2==0:
                        moments.append(np.sum(w*t**i*yu**j*yv**k))
        moments=np.array(moments)
        rows.append(np.r_[H[tri],q*moments,grad[0]*moments,grad[1]*moments]);targets.append(gap(x,w))
    A=np.array(rows);b=np.array(targets)
    c,res,rank,s=np.linalg.lstsq(A,b,rcond=None)
    print('rank',rank,'maxerr',np.max(abs(A@c-b)),'rms',np.sqrt(np.mean((A@c-b)**2)))
    ch=c[:len(tri[0])]
    C=np.zeros((6,6));C[np.diag_indices(6)]=ch[np.arange(6)] if False else 0
    for coeff,i,j in zip(ch,*tri):
        if i==j:C[i,j]=coeff
        else:C[i,j]=C[j,i]=coeff/2
    print('coeff H upper',ch)
    print('functional matrix eig',np.linalg.eigvalsh(C))
    print(C)


if __name__=='__main__':main()
