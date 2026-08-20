import numpy as np
from scipy.optimize import linprog


def evaluate(p, w):
    g = p @ p.T
    e = np.einsum('i,j,ij->', w, w,
                  32*g**6 - 48*g**4 + 20*g**2 - 4/3)
    a = g[:, :, None]
    b = g[:, None, :]
    c = g[None, :, :]
    d = 1 + 2*a*b*c - a*a - b*b - c*c
    term = 8*a*a*b*b*(c-a*b)**2*d
    f = np.einsum('i,j,k,ijk->', w, w, w, term)
    return e, f


rng = np.random.default_rng(20260819)
best = (1e100, None)
for trial in range(300):
    n = 16
    p = rng.normal(size=(n, 3))
    p /= np.linalg.norm(p, axis=1)[:, None]
    feats = np.vstack([
        np.ones(n),
        p[:, 0]**2, p[:, 1]**2,
        p[:, 0]*p[:, 1], p[:, 0]*p[:, 2], p[:, 1]*p[:, 2],
    ])
    target = np.array([1, 1/3, 1/3, 0, 0, 0])
    sol = linprog(rng.normal(size=n), A_eq=feats, b_eq=target,
                  bounds=(0, None), method='highs')
    if not sol.success:
        continue
    keep = sol.x > 1e-9
    e, f = evaluate(p[keep], sol.x[keep])
    ratio = e/f if f > 1e-13 else 1e100
    if ratio < best[0]:
        best = ratio, (trial, keep.sum(), e, f, p[keep], sol.x[keep])
        print('best', ratio, best[1][:4])
print('final', best[0], best[1][:4] if best[1] else None)

