"""Adversarial test of the root-conditional pair-circle inequality.

For fixed leaves z,w, symmetrize the four-sample target only in the two
root variables x,y.  If the resulting kernel is copositive for every z,w,
integrating z,w against the original measure proves the global target.
"""

import argparse
import itertools
import numpy as np
import torch

torch.set_default_dtype(torch.float64)


def canonical_expression(edge):
    """Canonical four-vertex expression, with edge[b,...,i,j]."""
    a=edge[...,0,1];u=edge[...,0,2];r=edge[...,0,3]
    v=edge[...,1,2];t=edge[...,1,3]
    U=u*u+v*v-2*a*u*v
    V=r*r+t*t-2*a*r*t
    L=u*r+v*t-a*(u*t+v*r)
    UV=U*V
    circle=a**4*(32*L**6-48*L**4*UV+20*L**2*UV**2-2*UV**3)
    determinant=1+2*a*u*v-a*a-u*u-v*v
    rooted=8*a*a*u*u*(v-a*u)**2*determinant
    base=144*a**6-216*a**4+87*a**2-5
    return base-108*rooted-288*circle


def values(raw_roots, logits, raw_leaves, fully_symmetric=False):
    roots=raw_roots/torch.linalg.vector_norm(raw_roots,dim=-1,keepdim=True)
    leaves=raw_leaves/torch.linalg.vector_norm(raw_leaves,dim=-1,keepdim=True)
    z,w=leaves[:,0],leaves[:,1]
    weights=torch.softmax(logits,dim=-1)
    a=torch.einsum("bni,bmi->bnm",roots,roots)
    uz=torch.einsum("bni,bi->bn",roots,z)
    uw=torch.einsum("bni,bi->bn",roots,w)
    left_z,right_z=uz[:,:,None],uz[:,None,:]
    left_w,right_w=uw[:,:,None],uw[:,None,:]
    U=left_z**2+right_z**2-2*a*left_z*right_z
    V=left_w**2+right_w**2-2*a*left_w*right_w
    L=left_z*left_w+right_z*right_w-a*(left_z*right_w+right_z*left_w)
    UV=U*V
    circle=a**4*(32*L**6-48*L**4*UV+20*L**2*UV**2-2*UV**3)
    determinant=1+2*a*left_z*right_z-a*a-left_z**2-right_z**2
    rooted=8*a*a*left_z**2*(right_z-a*left_z)**2*determinant
    rooted=(rooted+rooted.transpose(1,2))/2
    base=144*a**6-216*a**4+87*a**2-5
    kernel=base-108*rooted-288*circle
    if fully_symmetric:
        # Build the 4 by 4 Gram matrix for (root_i,root_j,z,w), retaining
        # the two root indices as batch axes, and average the canonical
        # expression over all label permutations.
        batch,size,_=roots.shape
        edge=torch.ones((batch,size,size,4,4),dtype=roots.dtype)
        edge[...,0,1]=edge[...,1,0]=a
        edge[...,0,2]=edge[...,2,0]=left_z
        edge[...,1,2]=edge[...,2,1]=right_z
        edge[...,0,3]=edge[...,3,0]=left_w
        edge[...,1,3]=edge[...,3,1]=right_w
        zw=torch.einsum("bi,bi->b",z,w)[:,None,None]
        edge[...,2,3]=edge[...,3,2]=zw
        kernel=sum(
            canonical_expression(edge[...,permutation,:][...,:,permutation])
            for permutation in itertools.permutations(range(4))
        )/24
    value=torch.einsum("bi,bij,bj->b",weights,kernel,weights)
    return value,roots,weights,leaves,kernel


def run(n=8,batch=96,steps=8000,seed=912,fully_symmetric=False):
    rng=np.random.default_rng(seed+n)
    raw_roots=torch.tensor(rng.normal(size=(batch,n,3)),requires_grad=True)
    logits=torch.tensor(rng.normal(size=(batch,n)),requires_grad=True)
    raw_leaves=torch.tensor(rng.normal(size=(batch,2,3)),requires_grad=True)
    optimizer=torch.optim.Adam([raw_roots,logits,raw_leaves],lr=.025)
    for step in range(steps):
        optimizer.zero_grad();value,*_=values(raw_roots,logits,raw_leaves,fully_symmetric)
        value.sum().backward();optimizer.step()
        if step in (steps//2,3*steps//4):
            for group in optimizer.param_groups:group["lr"]*=.2
    value,roots,weights,leaves,kernel=values(raw_roots,logits,raw_leaves,fully_symmetric)
    index=int(torch.argmin(value))
    print("n,batch,steps,index,value",n,batch,steps,index,float(value[index]))
    print("kernel min eig",float(torch.linalg.eigvalsh(kernel[index])[0]))
    print("weights",weights[index].detach().numpy())
    print("roots",roots[index].detach().numpy())
    print("leaves",leaves[index].detach().numpy())


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--n",type=int,default=8)
    parser.add_argument("--batch",type=int,default=96)
    parser.add_argument("--steps",type=int,default=8000)
    parser.add_argument("--seed",type=int,default=912)
    parser.add_argument("--fully-symmetric",action="store_true")
    run(**vars(parser.parse_args()))
