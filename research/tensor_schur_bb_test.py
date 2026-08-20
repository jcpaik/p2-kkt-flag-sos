"""Test the sharp Schur-complement pairing against BB* in spin coordinates."""

import itertools, math, sys
sys.path.insert(0,'research')
import numpy as np
import cvxpy as cp
from scipy.linalg import null_space

from tensor_hankel_fw import maps

words=list(itertools.product(range(3),repeat=3));wi={w:i for i,w in enumerate(words)}
comps=[(a,b,3-a-b) for a in range(4) for b in range(4-a)]
S=np.zeros((27,10))
for j,a in enumerate(comps):
 W=[w for w in words if tuple(w.count(i) for i in range(3))==a]
 for w in W:S[wi[w],j]=1/np.sqrt(len(W))
eps=np.zeros((3,3,3));eps[0,1,2]=eps[1,2,0]=eps[2,0,1]=1;eps[0,2,1]=eps[2,1,0]=eps[1,0,2]=-1
gens=[]
for axis in range(3):
 L=np.array([[eps[i,axis,j] for j in range(3)] for i in range(3)])
 G=np.kron(np.kron(L,np.eye(3)),np.eye(3))+np.kron(np.kron(np.eye(3),L),np.eye(3))+np.kron(np.kron(np.eye(3),np.eye(3)),L)
 gens.append(S.T@G@S)
J=np.zeros((10,3))
for a in range(3):
 T=np.zeros((3,3,3))
 for i,j,k in words:T[i,j,k]=((i==j and k==a)+(i==k and j==a)+(j==k and i==a))/np.sqrt(15)
 J[:,a]=S.T@T.reshape(-1)
U=null_space(J.T);g3=[U.T@g@U for g in gens]
def cas_full(X):return -sum(g@(g@X-X@g)-(g@X-X@g)@g for g in gens)
def proj_full(X,l):
 o=X.copy()
 for j in [0,2,4,6]:
  if j!=l:o=(cas_full(o)-j*(j+1)*o)/(l*(l+1)-j*(j+1))
 return o
def cas3(X):return -sum(g@(g@X-X@g)-(g@X-X@g)@g for g in g3)
def proj3(X,l):
 o=X.copy()
 for j in [0,2,4,6]:
  if j!=l:o=(cas3(o)-j*(j+1)*o)/(l*(l+1)-j*(j+1))
 return o

# An orthonormal basis for the unique spin-4 Hankel component.
spin4=[]
for k in range(28):spin4.append(proj_full(maps[3][:,:,k],4))
V=np.stack([x.reshape(-1) for x in spin4],axis=1);uu,ss,_=np.linalg.svd(V,full_matrices=False);basis=[uu[:,i].reshape(10,10) for i in range(9)]

def invariants(coeff):
 H=sum(c*M for c,M in zip(coeff,basis));A4=U.T@H@U;B=U.T@H@J
 b=np.sum(B*B);A4=A4/np.sqrt(b);B=B/np.sqrt(b) # normalize b=1
 Q=B@B.T;q=np.sum(Q*Q);c=np.sum(A4*proj3(Q,4));d=np.sum(proj3(Q,6)**2)
 return q,c,d

rng=np.random.default_rng(772);worst=(1e9,None)
for trial in range(10000):
 coeff=rng.normal(size=9);inv=invariants(coeff)
 for b in np.linspace(.0001,.08,400):
  # Scale B by sqrt(b), A4 likewise. RHS of Schur pairing.
  q,c,d=inv
  rhs=5*b*b*q-(2/35)*b-(b**1.5)*c
  lower=max(0,rhs)**2/(b*b*d) if d>1e-16 else 0
  gap=2/105-(5/11)*b+2*lower
  if gap<worst[0]:worst=(gap,(trial,b,inv,rhs,lower,coeff.copy()))
print('worst',worst)

# Solve the exact spin-6 completion for the worst raw direction.
coeff=worst[1][-1];H4=sum(c*M for c,M in zip(coeff,basis));A4=U.T@H4@U;B0=U.T@H4@J
sc=np.sqrt(np.sum(B0*B0));A4/=sc;B0/=sc
spin6=[proj_full(maps[3][:,:,k],6) for k in range(28)]
VV=np.stack([x.reshape(-1) for x in spin6],axis=1);u6,s6,_=np.linalg.svd(VV,full_matrices=False)
A6basis=[U.T@u6[:,i].reshape(10,10)@U for i in range(13)]
z=cp.Variable(13);A6=sum(z[i]*A6basis[i] for i in range(13))
bp=cp.Parameter(nonneg=True);R=(2/35)*np.eye(7)+cp.sqrt(bp)*A4+A6-5*bp*(B0@B0.T)
prob=cp.Problem(cp.Minimize(cp.sum_squares(A6)),[R>>0])
for btest in [worst[1][1],.04,.02,.01]:
 bp.value=btest
 try:
  val=prob.solve(solver='CLARABEL');print('completion',btest,prob.status,val,'gap',2/105-5*btest/11+2*val)
 except Exception as e:print(e)
