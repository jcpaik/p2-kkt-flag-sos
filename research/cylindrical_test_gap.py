import numpy as np


def K(t): return 32*t**6-48*t**4+20*t**2-4/3


def evaluate(p,w):
    g=p@p.T
    energy=np.einsum('i,j,ij->',w,w,K(g))
    gap=0
    for i,x in enumerate(p):
        u=np.cross(x,[1.,0,0]) if abs(x[0])<.9 else np.cross(x,[0.,1,0])
        u/=np.linalg.norm(u);v=np.cross(x,u)
        z=g[i]**2
        A=np.dot(w,z*(1-z)**2)
        B=np.dot(w,z*((p@u)+1j*(p@v))**4)
        gap+=w[i]*(A-abs(B))
    return energy,gap


rng=np.random.default_rng(22);best=(1e9,None)
for n in range(3,21):
    for trial in range(1000):
        p=rng.normal(size=(n,3));p/=np.linalg.norm(p,axis=1)[:,None]
        w=rng.dirichlet(np.ones(n))
        e,g=evaluate(p,w)
        d=e-4*g
        if d<best[0]:best=(d,(n,e,g))
    print(n,best)

