"""Test the induced spin-3 principal-block inequality."""
import itertools,math
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize

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
U=null_space(J.T);gs=[U.T@g@U for g in gens]
def cas(X):return -sum(g@(g@X-X@g)-(g@X-X@g)@g for g in gs)
def proj(X,l):
 o=X.copy()
 for j in [0,2,4,6]:
  if j!=l:o=(cas(o)-j*(j+1)*o)/(l*(l+1)-j*(j+1))
 return o
def unpack(z,r):
 Z=z.reshape(7,r);A=Z@Z.T;return .4*A/np.trace(A)
def q(A):return 2/105+75/32*np.sum(proj(A,2)**2)-5/3*np.sum(proj(A,4)**2)+2*np.sum(proj(A,6)**2)
rng=np.random.default_rng(33)
for r in range(1,8):
 best=(1e9,None)
 for st in range(1):
  z=rng.normal(size=7*r);res=minimize(lambda zz:q(unpack(zz,r)),z,method='L-BFGS-B',options={'maxiter':300,'gtol':1e-8})
  if res.fun<best[0]:best=(res.fun,np.linalg.eigvalsh(unpack(res.x,r)))
 print(r,best)
