"""Canonical ONB transform eta_mu with I(mu,eta_mu)=F(mu).

For an ordered pair (x,y), let R(x,y) be the orthonormal frame consisting
of x, the normalized tangent projection of y, and x cross y.  With
w=t^2(1-t^2)^2 and nu_R the uniform measure on the three frame axes,

    eta_mu = (1/4) E_{x,y}[w nu_R].

Then eta has mass m=(1/4)Ew and I(mu,eta)=F.
"""

import numpy as np


def kernel(t):
    return 32*t**6 - 48*t**4 + 20*t**2 - 4/3


def eta_atoms(x, weights, tol=1e-14):
    atoms = []
    masses = []
    n = len(weights)
    for i in range(n):
        for j in range(n):
            t = float(x[i] @ x[j])
            s2 = 1-t*t
            pair_weight = weights[i]*weights[j]*t*t*s2*s2/12
            if pair_weight <= tol:
                continue
            u = x[j]-t*x[i]
            u /= np.linalg.norm(u)
            v = np.cross(x[i],u)
            atoms.extend((x[i],u,v))
            masses.extend((pair_weight,)*3)
    if not atoms:
        return np.empty((0,3)),np.empty(0)
    return np.array(atoms),np.array(masses)


def interaction(x,w,y,v):
    return w @ kernel(x@y.T) @ v


def quantities(x,w):
    eta,v=eta_atoms(x,w)
    E=interaction(x,w,x,w)
    if len(v)==0:
        return E,0.,0.,0.,0.
    m=v.sum()
    F=interaction(x,w,eta,v)
    Eeta=interaction(eta,v,eta,v)
    return E,m,F,Eeta,Eeta-2*m*F


def main():
    rng=np.random.default_rng(20260820)
    worst=(-np.inf,None)
    for n in (3,4,5,8,12):
        vals=[]
        for _ in range(1000):
            x=rng.normal(size=(n,3));x/=np.linalg.norm(x,axis=1)[:,None]
            w=rng.dirichlet(np.ones(n))
            q=quantities(x,w)
            vals.append(q[-1])
            if q[-1]>worst[0]:worst=(q[-1],(x,w,q))
        print(n,min(vals),max(vals),np.mean(vals))
    print('worst',worst[0],worst[1][2])
    print(worst[1][1]);print(worst[1][0])


if __name__=='__main__':main()
