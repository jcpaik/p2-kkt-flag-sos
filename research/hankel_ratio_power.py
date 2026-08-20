"""Alternating conic power method for the Hankel Q positive/negative split."""
import math
import cvxpy as cp
import numpy as np

def comps(n,k=3):
 if k==1:return [(n,)]
 return [(a,)+t for a in range(n+1) for t in comps(n-a,k-1)]
def multinomial(a):
 o=math.factorial(sum(a))
 for q in a:o//=math.factorial(q)
 return o
d6=comps(6);ix={a:i for i,a in enumerate(d6)}
def mmap(d):
 b=comps(d);M=np.zeros((len(b),len(b),28));rem=comps(3-d)
 for i,a in enumerate(b):
  for j,c in enumerate(b):
   sc=np.sqrt(multinomial(a)*multinomial(c))
   for g in rem:
    ex=tuple(a[q]+c[q]+2*g[q] for q in range(3));M[i,j,ix[ex]]+=sc*multinomial(g)
 return M
maps={d:mmap(d) for d in [1,2,3]};grams={d:np.einsum('abk,abl->kl',M,M) for d,M in maps.items()};mass=np.trace(maps[3],axis1=0,axis2=1)
Q=32*grams[3]-48*grams[2]+20*grams[1]-4/3*np.outer(mass,mass)
ev,V=np.linalg.eigh(Q);P=(np.sqrt(ev[ev>1e-8])[:,None]*V[:,ev>1e-8].T);N=(np.sqrt(-ev[ev<-1e-8])[:,None]*V[:,ev<-1e-8].T)
if __name__ == "__main__":
 print('split',P.shape,N.shape)
 y=cp.Variable(28);direction=cp.Parameter(9);H=sum(y[k]*maps[3][:,:,k] for k in range(28))
 prob=cp.Problem(cp.Maximize(direction@(N@y)),[H>>0,cp.norm(P@y)<=1])
 rng=np.random.default_rng(122)
 best=(0,None)
 for st in range(50):
  d=rng.normal(size=9);d/=np.linalg.norm(d)
  for it in range(50):
   direction.value=d;val=prob.solve(solver='CLARABEL');yy=y.value;nn=N@yy;d2=nn/np.linalg.norm(nn)
   if np.linalg.norm(d2-d)<1e-10:break
   d=d2
  ratio=np.linalg.norm(N@yy)/np.linalg.norm(P@yy)
  if ratio>best[0]:best=(ratio,(yy,d,it));print('best',st,it,ratio,'Q',yy@Q@yy,'mass',mass@yy,'eig',np.linalg.eigvalsh(sum(yy[k]*maps[3][:,:,k] for k in range(28))))
