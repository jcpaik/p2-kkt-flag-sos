import itertools
import numpy as np
from scipy.optimize import lsq_linear


def exponents(total):
    return [(i,j,total-i-j) for i in range(total+1)
            for j in range(total-i+1)]


exps=exponents(4)+exponents(6)


def feature(points, weights=None):
    if weights is None: weights=np.ones(len(points))/len(points)
    vals=[]
    for e in exps:
        vals.append(np.dot(weights,np.prod(points**np.array(e),axis=1)))
    return np.array(vals)


rng=np.random.default_rng(84)
frames=[]
for _ in range(1000):
    q,_=np.linalg.qr(rng.normal(size=(3,3)))
    frames.append(feature(q.T))
A=np.asarray(frames).T
Aaug=np.vstack([A,10*np.ones(A.shape[1])])

tetra=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3)
phi=(1+np.sqrt(5))/2
ico=[]
for a in [-1,1]:
    for b in [-phi,phi]: ico += [[0,a,b],[a,b,0],[b,0,a]]
ico=np.unique(np.asarray(ico),axis=0);ico/=np.linalg.norm(ico,axis=1)[:,None]

for name,p in [('tetra',tetra),('ico',ico)]:
    target=feature(p)
    sol=lsq_linear(Aaug,np.r_[target,10],bounds=(0,np.inf),
                   tol=1e-10,lsmr_tol=1e-10,max_iter=500,verbose=0)
    err=np.linalg.norm(A@sol.x-target)
    print(name,sol.cost,err,sol.x.sum(),np.count_nonzero(sol.x>1e-7),sol.optimality)
