import numpy as np


def kernel(t):
    return 32*t**6 - 48*t**4 + 20*t**2 - 4/3


def evals(points, weights):
    g = points @ points.T
    e = np.sum(weights[:, None] * weights[None, :] * kernel(g))
    f = 0.0
    for i in range(len(points)):
        for j in range(len(points)):
            a = g[i, j]
            for k in range(len(points)):
                b, c = g[i, k], g[j, k]
                d = 1 + 2*a*b*c - a*a - b*b - c*c
                f += weights[i]*weights[j]*weights[k] * 8*a*a*b*b*(c-a*b)**2*d
    return e, f


rng = np.random.default_rng(1024)
for nbase in range(1, 7):
    for trial in range(3):
        pts=[]
        base_weights=rng.dirichlet(np.ones(nbase))
        ws=[]
        for wb in base_weights:
            q,_=np.linalg.qr(rng.normal(size=(3,3)))
            pts.extend(q.T)
            ws.extend([wb/3]*3)
        e,f=evals(np.asarray(pts),np.asarray(ws))
        print(nbase, trial, e, f, e/f if f>1e-13 else np.nan)

