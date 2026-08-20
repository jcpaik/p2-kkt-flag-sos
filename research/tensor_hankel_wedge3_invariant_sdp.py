"""Exact-form search mass*E=Tr(W wedge^3 H), W SO(3)-invariant PSD."""

import itertools
import sys
sys.path.insert(0,"research")
import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

from tensor_hankel_fw import maps, compositions

def sym_generator(degree,L):
    basis=compositions(degree,3);index={a:i for i,a in enumerate(basis)}
    G=np.zeros((len(basis),len(basis)))
    for col,beta in enumerate(basis):
        for source in range(3):
            for target in range(3):
                if beta[source]==0 or target==source or L[target,source]==0:continue
                alpha=list(beta);alpha[source]-=1;alpha[target]+=1
                G[index[tuple(alpha)],col]+=L[target,source]*np.sqrt(beta[source]*(beta[target]+1))
    return G
eps=np.zeros((3,3,3));eps[0,1,2]=eps[1,2,0]=eps[2,0,1]=1;eps[0,2,1]=eps[2,1,0]=eps[1,0,2]=-1
Ls=[np.array([[eps[i,a,j] for j in range(3)] for i in range(3)]) for a in range(3)]
GV=[sym_generator(3,L) for L in Ls]

inds=list(itertools.combinations(range(10),3));ix={a:i for i,a in enumerate(inds)}

def wedgegen(G):
    W=np.zeros((120,120))
    for col,I in enumerate(inds):
        for pos,src in enumerate(I):
            for tgt in range(10):
                if G[tgt,src]==0 or (tgt in I and tgt!=src):continue
                J=list(I);J[pos]=tgt
                inv=sum(J[a]>J[b] for a in range(3) for b in range(a+1,3))
                W[ix[tuple(sorted(J))],col]+=(-1)**inv*G[tgt,src]
    return W

WG=[wedgegen(g) for g in GV];Cas=-sum(g@g for g in WG)
ce,cv=np.linalg.eigh(Cas)

def symmetric_basis(d):
    out=[]
    for i in range(d):
        E=np.zeros((d,d));E[i,i]=1;out.append(E)
        for j in range(i+1,d):
            E=np.zeros((d,d));E[i,j]=E[j,i]=1/np.sqrt(2);out.append(E)
    return out

blocks={}
for j in range(7):
    U=cv[:,np.abs(ce-j*(j+1))<1e-7]
    if not U.shape[1]:continue
    GG=[U.T@g@U for g in WG];sb=symmetric_basis(U.shape[1])
    equations=np.empty((3*U.shape[1]**2,len(sb)))
    for k,B in enumerate(sb):equations[:,k]=np.concatenate([(g@B-B@g).reshape(-1) for g in GG])
    if np.linalg.norm(equations)<1e-9:
        ns=np.eye(len(sb))
    else:
        ns=null_space(equations,rcond=1e-9)
    BB=[]
    for k in range(ns.shape[1]):BB.append(sum(ns[p,k]*sb[p] for p in range(len(sb))))
    blocks[j]=(U,BB)
    print('block',j,U.shape[1],'comm',len(BB))

Iarr=np.array(inds)
perms=list(itertools.permutations(range(3)))
sign=[]
for p in perms:sign.append((-1)**sum(p[a]>p[b] for a in range(3) for b in range(a+1,3)))
def compound(H):
    out=np.zeros((120,120))
    for p,sgn in zip(perms,sign):
        out += sgn*(H[Iarr[:,None,0],Iarr[None,:,p[0]]]
                    *H[Iarr[:,None,1],Iarr[None,:,p[1]]]
                    *H[Iarr[:,None,2],Iarr[None,:,p[2]]])
    return out

gram={d:np.einsum('abk,abl->kl',M,M) for d,M in maps.items()}
mass=np.trace(maps[3],axis1=0,axis2=1)
Q=32*gram[3]-48*gram[2]+20*gram[1]-(4/3)*np.outer(mass,mass)

theta={j:cp.Variable(len(BB)) for j,(U,BB) in blocks.items()}
X={j:sum(theta[j][k]*B for k,B in enumerate(BB)) for j,(U,BB) in blocks.items()}
cons=[X[j]>>0 for j in blocks]
rng=np.random.default_rng(1881)
Z=rng.normal(size=(100,28))
def coeffs(z):
    H=np.einsum('abk,k->ab',maps[3],z);C=compound(H);out={}
    for j,(U,BB) in blocks.items():
        D=U.T@C@U;out[j]=np.array([np.sum(B*D) for B in BB])
    return out
rows=[coeffs(z) for z in Z]
b=np.array([(mass@z)*(z@Q@z) for z in Z])
for row,val in zip(rows,b):cons.append(sum(theta[j]@row[j] for j in blocks)==val)
prob=cp.Problem(cp.Minimize(sum(cp.trace(X[j]) for j in blocks)),cons)
for solver in (cp.CLARABEL,cp.SCS):
 try:
  val=prob.solve(solver=solver,verbose=False);print(solver,prob.status,val)
  if prob.status in ('optimal','optimal_inaccurate'):
   for j in blocks:print(j,np.linalg.eigvalsh(X[j].value),theta[j].value)
   err=[]
   for z in rng.normal(size=(100,28)):
    row=coeffs(z);err.append(sum(theta[j].value@row[j] for j in blocks)-(mass@z)*(z@Q@z))
   print('extra',max(abs(np.array(err))))
 except cp.error.SolverError as e:print(solver,e)
