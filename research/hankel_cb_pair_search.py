"""Search mixtures of rank-seven Cayley--Bacharach Hankel rays."""
import math
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
maps={d:mmap(d) for d in [1,2,3]}
grams={d:np.einsum('abk,abl->kl',M,M) for d,M in maps.items()}
mass=np.trace(maps[3],axis1=0,axis2=1)
Q=32*grams[3]-48*grams[2]+20*grams[1]-4/3*np.outer(mass,mass)

def grid_data(xs,ys):
 pts=[];u=[]
 for i,x in enumerate(xs):
  px=np.prod([x-xs[k] for k in range(3) if k!=i])
  for j,y in enumerate(ys):
   py=np.prod([y-ys[k] for k in range(3) if k!=j])
   raw=np.array([x,y,1.]);r=np.linalg.norm(raw);pts.append(raw/r);u.append(r**3/(px*py))
 return np.array(pts),np.array(u)
def ray(rng):
 while True:
  xs=np.sort(rng.normal(size=3));ys=np.sort(rng.normal(size=3));pts,u=grid_data(xs,ys);neg=int(rng.integers(9))
  # A generic projective image substantially enlarges the reducible-grid family.
  U,_,Vt=np.linalg.svd(rng.normal(size=(3,3)));A=U@np.diag(np.exp(rng.normal(scale=1.2,size=3)))@Vt
  raw=pts@A.T;sc=np.linalg.norm(raw,axis=1);pts=raw/sc[:,None];u=u*sc**3
  pos=np.exp(rng.normal(scale=1.5,size=8));inds=[i for i in range(9) if i!=neg]
  w=np.zeros(9);w[inds]=pos;w[neg]=-u[neg]**2/np.sum(u[inds]**2/pos);sm=w.sum()
  if sm>0 and -w[neg]/sm>.001:
   w/=sm;y=np.array([sum(wi*np.prod(p**np.array(a)) for wi,p in zip(w,pts)) for a in d6])
   return y,(xs,ys,neg,w,pts)
rng=np.random.default_rng(4433);N=12000
Y=[];meta=[]
for i in range(N):
 y,m=ray(rng);Y.append(y);meta.append(m)
Y=np.array(Y);q=np.einsum('bi,ij,bj->b',Y,Q,Y)
print('self min',q.min(),q.max())
best=(np.inf,None,None,None)
for lo in range(0,N,400):
 B=Y[lo:lo+400]@Q@Y.T;score=B+np.sqrt(np.maximum(q[lo:lo+400,None]*q[None,:],0))
 for k in range(score.shape[0]):score[k,lo+k]=np.inf
 aa,bb=np.unravel_index(np.argmin(score),score.shape)
 if score[aa,bb]<best[0]:best=(score[aa,bb],lo+aa,bb,B[aa,bb])
_,ii,jj,bij=best
print('pair',ii,jj,'score',best[0],'q',q[ii],q[jj],'b',bij)
# exact line minimizer on normalized convex segment
d=Y[jj]-Y[ii];aa=d@Q@d;bb=2*Y[ii]@Q@d;t=np.clip(-bb/(2*aa),0,1) if aa>0 else (1 if q[jj]<q[ii] else 0)
print('mix',t,(Y[ii]+t*d)@Q@(Y[ii]+t*d));print('meta1',meta[ii]);print('meta2',meta[jj])
