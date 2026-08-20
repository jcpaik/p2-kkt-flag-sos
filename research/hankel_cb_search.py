"""Search rank-7 extreme PSD Hankel rays from 3x3 Cayley--Bacharach grids."""
import numpy as np
from scipy.optimize import minimize

def K(t):return 32*t**6-48*t**4+20*t**2-4/3

def grid_data(xs,ys):
    pts=[];u=[]
    for i,x in enumerate(xs):
      px=np.prod([x-xs[k] for k in range(3) if k!=i])
      for j,y in enumerate(ys):
       py=np.prod([y-ys[k] for k in range(3) if k!=j])
       raw=np.array([x,y,1.]);r=np.linalg.norm(raw)
       pts.append(raw/r);u.append(r**3/(px*py))
    return np.array(pts),np.array(u)

def weights_from_z(u,z,negative=8):
    inds=[i for i in range(9) if i!=negative]
    pos=np.exp(z)
    w=np.zeros(9);w[inds]=pos
    w[negative]=-u[negative]**2/np.sum(u[inds]**2/pos)
    return w

def objective(z,pts,u,negative):
    w=weights_from_z(u,z,negative)
    mass=w.sum()
    if mass<=1e-12:return 1e3+(mass-1)**2
    G=pts@pts.T
    return (w@K(G)@w)/(mass*mass)

rng=np.random.default_rng(912)
best=(1e9,None)
# Broad Monte Carlo first; require a genuinely negative pseudo-weight.
for trial in range(20000):
 xs=np.sort(rng.normal(size=3));ys=np.sort(rng.normal(size=3))
 pts,u=grid_data(xs,ys);neg=int(rng.integers(9));z=rng.normal(scale=1.5,size=8)
 w=weights_from_z(u,z,neg);mass=w.sum()
 if mass>0 and -w[neg]/mass>.001:
  val=objective(z,pts,u,neg)
  if val<best[0]:best=(val,(xs,ys,neg,w/mass,None));print('random best',trial,val,-w[neg]/mass)
for trial in range(12):
 xs=np.sort(rng.normal(size=3));ys=np.sort(rng.normal(size=3))
 pts,u=grid_data(xs,ys)
 for neg in range(9):
  z0=rng.normal(scale=2,size=8)
  res=minimize(objective,z0,args=(pts,u,neg),method='L-BFGS-B',bounds=[(-8,8)]*8,options={'maxiter':500,'ftol':1e-13,'gtol':1e-10})
  if res.fun<best[0]:
   w=weights_from_z(u,res.x,neg);best=(res.fun,(xs,ys,neg,w/w.sum(),res));print('best',trial,best[0],xs,ys,neg,w/w.sum(),res.success)
