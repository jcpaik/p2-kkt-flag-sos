import numpy as np
from scipy.optimize import linprog


def data(p, w, degree=6):
    g = p @ p.T
    e = np.einsum('i,j,ij->', w, w, 32*g**6-48*g**4+20*g**2-4/3)
    a = g[:, :, None]
    b = g[:, None, :]
    c = g[None, :, :]
    d = 1 + 2*a*b*c-a*a-b*b-c*c
    base = 32*b*b*(c-a*b)**2*d
    vals = [np.einsum('i,j,k,ijk->',w,w,w,(a*a)**k*base)
            for k in range(degree+1)]
    return e, vals


rng=np.random.default_rng(904)
rows=[]; target=[]
for trial in range(200):
    n=18
    p=rng.normal(size=(n,3));p/=np.linalg.norm(p,axis=1)[:,None]
    feats=np.vstack([np.ones(n),p[:,0]**2,p[:,1]**2,
                     p[:,0]*p[:,1],p[:,0]*p[:,2],p[:,1]*p[:,2]])
    sol=linprog(rng.normal(size=n),A_eq=feats,
                b_eq=[1,1/3,1/3,0,0,0],bounds=(0,None),method='highs')
    if sol.success:
        keep=sol.x>1e-9
        e,v=data(p[keep],sol.x[keep])
        rows.append(v);target.append(e)

rows=np.asarray(rows);target=np.asarray(target)
for m in range(1,8):
    coef,res,rank,s=np.linalg.lstsq(rows[:,:m],target,rcond=None)
    err=np.max(np.abs(rows[:,:m]@coef-target))
    print(m,rank,err,coef)

