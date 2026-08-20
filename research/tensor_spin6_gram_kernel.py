"""Recover the exact contraction formula for <Pi6(vv*),Pi6(ww*)>.

Here v,w are spin-3 vectors, realized as symmetric trace-free rank-three
tensors on R^3.  The bidegree-(2,2) invariants are squared norms of their
r-fold contractions, r=0,1,2,3.  Numerical fitting followed by rational
reconstruction identifies the exact coefficients; a separate exact symbolic
verifier can then be built from rational STF bases.
"""

import sys

sys.path.insert(0, "research")

import numpy as np
import sympy as sp

import tensor_schur_feasible_scan as data


def tensor(v):
    """Spin-3 coordinate vector -> normalized symmetric rank-three tensor."""
    return (data.S @ data.U @ v).reshape(3, 3, 3)


def invariants(v, w):
    V = tensor(v)
    W = tensor(w)
    c0 = np.sum(V * V) * np.sum(W * W)
    c1_tensor = np.einsum("aij,akl->ijkl", V, W)
    c1 = np.sum(c1_tensor * c1_tensor)
    c2_tensor = np.einsum("abi,abj->ij", V, W)
    c2 = np.sum(c2_tensor * c2_tensor)
    c3 = np.sum(V * W) ** 2
    return np.array([c0, c1, c2, c3])


def gram(v, w):
    gv = data.proj3(np.outer(v, v), 6)
    gw = data.proj3(np.outer(w, w), 6)
    return np.sum(gv * gw)


def fit(samples=200):
    rng = np.random.default_rng(20260820)
    X = []
    y = []
    for _ in range(samples):
        v = rng.normal(size=7)
        w = rng.normal(size=7)
        X.append(invariants(v, w))
        y.append(gram(v, w))
    X = np.array(X)
    y = np.array(y)
    coeff, *_ = np.linalg.lstsq(X, y, rcond=None)
    print("coeff", coeff)
    print("rational", [sp.Rational(float(x)).limit_denominator(10**6) for x in coeff])
    print("residual", np.max(np.abs(X @ coeff - y)))


if __name__ == "__main__":
    fit()
