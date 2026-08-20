import itertools
import math
import numpy as np


triples=list(itertools.combinations_with_replacement(range(3),3))


def sym3(x):
    vals=[]
    for t in triples:
        counts=[t.count(i) for i in range(3)]
        mult=math.factorial(3)
        for c in counts: mult//=math.factorial(c)
        vals.append(math.sqrt(mult)*np.prod(x[list(t)]))
    return np.array(vals)


def moment(points,w):
    X=np.array([sym3(x) for x in points])
    return np.einsum('i,ia,ib->ab',w,X,X)


def samples(name):
    if name=='onb': return np.eye(3),np.ones(3)/3
    if name=='polehaar':
        m=1000; th=np.arange(m)*2*np.pi/m
        p=np.vstack([[0,0,1],np.c_[np.cos(th),np.sin(th),np.zeros(m)]])
        return p,np.r_[1/3,np.ones(m)*2/(3*m)]
    if name=='tetra':
        p=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3)
        return p,np.ones(4)/4


for name in ['onb','polehaar','tetra']:
    p,w=samples(name);M=moment(p,w)
    ev=np.linalg.eigvalsh(M)
    print(name,'rank',sum(ev>1e-9),'ev',np.round(ev,10),'purity',np.sum(M*M))

