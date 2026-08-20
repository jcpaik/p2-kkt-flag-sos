import numpy as np


def K(t): return 32*t**6-48*t**4+20*t**2-4/3


def residual_on_support(p,w):
    n=len(p);g=p@p.T
    U=g*0
    # U_mu at support
    Um=K(g)@w
    Un=np.zeros(n)
    for i,x in enumerate(p):
        for j,y in enumerate(p):
            a=g[i,j]; q=a*a*(1-a*a)**2
            if q<1e-20: continue
            e=(y-a*x)/np.sqrt(1-a*a)
            nn=np.cross(x,y)/np.sqrt(1-a*a)
            # potential of the ONB probability measure at every support root
            dots=np.c_[p@x,p@e,p@nn]
            pot=np.mean(K(dots),axis=1)
            Un += w[i]*w[j]*q*pot
    return Um-4*Un,Um,Un


rng=np.random.default_rng(676)
best=(1e9,None)
for n in range(2,16):
    for trial in range(500):
        p=rng.normal(size=(n,3));p/=np.linalg.norm(p,axis=1)[:,None]
        w=rng.dirichlet(np.ones(n))
        r,um,un=residual_on_support(p,w)
        if r.min()<best[0]:best=(r.min(),(n,trial,r,um,un,p,w))
    print(n,best[0])

