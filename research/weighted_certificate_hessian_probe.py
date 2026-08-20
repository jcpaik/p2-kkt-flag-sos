"""Probe products of the rooted spin-4 determinant with support Hessians."""

import numpy as np


def kernel(t):
    return 32*t**6 - 48*t**4 + 20*t**2 - 4/3


def kp(t):
    return 192*t**5 - 192*t**3 + 40*t


def kpp(t):
    return 960*t**4 - 576*t**2 + 40


def values(x, w):
    g = x @ x.T
    n = len(w)
    potential = kernel(g) @ w
    Droot = np.zeros(n)
    trace_hessian = np.zeros(n)
    det_hessian = np.zeros(n)
    adj_pair = np.zeros(n)
    for i in range(n):
        # deterministic tangent frame
        seed = np.eye(3)[np.argmin(np.abs(x[i]))]
        u = np.cross(x[i], seed)
        u /= np.linalg.norm(u)
        v = np.cross(x[i], u)
        t = g[i]
        zu = x @ u
        zv = x @ v
        fc = t*(zu**2-zv**2)
        fs = 2*t*zu*zv
        M = np.array([
            [np.sum(w*fc*fc), np.sum(w*fc*fs)],
            [np.sum(w*fc*fs), np.sum(w*fs*fs)],
        ])
        # In this real normalization det M = D_x/4; record D_x itself.
        Droot[i] = 4*np.linalg.det(M)
        H = np.zeros((2, 2))
        for a, za in enumerate((zu, zv)):
            for b, zb in enumerate((zu, zv)):
                H[a, b] = np.sum(w*(kpp(t)*za*zb - (kp(t)*t if a == b else 0)))
        trace_hessian[i] = np.trace(H)
        det_hessian[i] = np.linalg.det(H)
        adjM = np.array([[M[1,1], -M[0,1]],[-M[1,0],M[0,0]]])
        adj_pair[i] = np.sum(adjM*H)
    F = np.sum(w*Droot)
    E = w @ potential
    return {
        "E": E,
        "F": F,
        "EF_support": np.sum(w*Droot*potential),
        "D_traceH": np.sum(w*Droot*trace_hessian),
        "D_detH": np.sum(w*Droot*det_hessian),
        "adjM_H": np.sum(w*adj_pair),
        "min_traceH": trace_hessian.min(),
    }


def commutator_transport_mc(x, w, samples=200_000, seed=0):
    """MC value of the projectively valid field phi_x P_x r.

    The sum over an orthonormal basis of ambient vectors r has local kernel
    L and cross kernel C.  Shared leaves Y,Z define
    phi_x=(x.y)(x.z)((P_x y).(P_x z)) det(x,y,z).
    """
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(w), size=(samples, 4), p=w)
    X, W, Y, Z = (x[idx[:,i]] for i in range(4))
    s = np.einsum("ij,ij->i", X, W)
    a = np.einsum("ij,ij->i", X, Y)
    b = np.einsum("ij,ij->i", X, Z)
    c = np.einsum("ij,ij->i", Y, Z)
    dx = np.einsum("ij,ij->i", X, np.cross(Y,Z))
    phix = a*b*(c-a*b)*dx
    aw = np.einsum("ij,ij->i", W, Y)
    bw = np.einsum("ij,ij->i", W, Z)
    dw = np.einsum("ij,ij->i", W, np.cross(Y,Z))
    phiw = aw*bw*(c-aw*bw)*dw
    L = (1-s*s)*kpp(s)-2*s*kp(s)
    C = -s*(1-s*s)*kpp(s)+(1+s*s)*kp(s)
    q = np.mean(phix*phix*L+phix*phiw*C)
    target = np.mean(8*phix*phix*kernel(s))
    f = np.mean(8*phix*phix)
    return q, target, f


def main():
    rng = np.random.default_rng(90210)
    for n in (3, 5, 8, 12):
        rows = []
        for _ in range(1000):
            x = rng.normal(size=(n,3)); x /= np.linalg.norm(x,axis=1)[:,None]
            w = rng.dirichlet(np.ones(n))
            rows.append(values(x,w))
        for key in ("D_traceH", "adjM_H"):
            ratio = [r[key]/r["EF_support"] for r in rows if abs(r["EF_support"]) > 1e-8]
            print(n,key,min(ratio),max(ratio),np.median(ratio))


if __name__ == "__main__":
    main()
