"""Sample-screen mass*E as a PSD-weighted third compound of H."""

import itertools
import cvxpy as cp
import numpy as np

from tensor_hankel_fw import maps

inds=list(itertools.combinations(range(10),3))
gram={d:np.einsum('abk,abl->kl',M,M) for d,M in maps.items()}
mass=np.trace(maps[3],axis1=0,axis2=1)
Q=32*gram[3]-48*gram[2]+20*gram[1]-(4/3)*np.outer(mass,mass)

def compound(H):
    C=np.empty((120,120))
    for i,I in enumerate(inds):
        for j,J in enumerate(inds):
            C[i,j]=np.linalg.det(H[np.ix_(I,J)])
    return (C+C.T)/2

rng=np.random.default_rng(665)
Z=rng.normal(size=(450,28))
Cs=[];b=[]
for k,z in enumerate(Z):
    H=np.einsum('abk,k->ab',maps[3],z)
    Cs.append(compound(H));b.append((mass@z)*(z@Q@z))
print('built')
W=cp.Variable((120,120),symmetric=True)
cons=[W>>0]+[cp.sum(cp.multiply(W,C))==v for C,v in zip(Cs,b)]
prob=cp.Problem(cp.Minimize(cp.trace(W)),cons)
for solver in (cp.CLARABEL,cp.SCS):
 try:
  val=prob.solve(solver=solver,verbose=False)
  print(solver,prob.status,val)
  if W.value is not None:
   print('eig',np.linalg.eigvalsh(W.value)[:20]);np.save('/tmp/tensor_wedge3_W.npy',W.value)
   err=[]
   for z in rng.normal(size=(100,28)):
    H=np.einsum('abk,k->ab',maps[3],z);err.append(np.sum(W.value*compound(H))-(mass@z)*(z@Q@z))
   print('extra',max(abs(np.array(err))))
 except cp.error.SolverError as e:print(e)
