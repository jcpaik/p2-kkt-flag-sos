"""Numerical audit of E(mu) >= min_R 32 int prod_i (r_i.x)^2 dmu."""

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation


# Robust counterexample found by tensor_onb_minimax_opt.py.  The displayed
# precision already leaves a 3.28e-2 gap after 500 independent frame starts.
COUNTEREXAMPLE_X = np.array([
    [-0.07442228, -0.18480423,  0.97995343],
    [ 0.36216506,  0.93199640,  0.01480492],
    [ 0.95835542,  0.24795516,  0.14167971],
    [-0.65988406,  0.32568766, -0.67711194],
    [-0.38194692,  0.36212165, -0.85028493],
    [ 0.15039476, -0.94817296, -0.27990973],
    [-0.43797216,  0.62442433,  0.64674156],
    [-0.78085683, -0.62463822,  0.00947162],
    [ 0.89120584, -0.13962169,  0.43157610],
    [ 0.38432886, -0.18532284, -0.90440410],
])
COUNTEREXAMPLE_X /= np.linalg.norm(COUNTEREXAMPLE_X, axis=1)[:, None]
COUNTEREXAMPLE_W = np.array([
    0.10490798, 0.11713372, 0.09407234, 0.05145806, 0.09700014,
    0.12153769, 0.11364985, 0.09665728, 0.10319469, 0.10038825,
])
COUNTEREXAMPLE_W /= COUNTEREXAMPLE_W.sum()


def energy(xs, weights):
    g = xs @ xs.T
    return float(np.sum(weights[:, None] * weights[None, :] * (32*g**6-48*g**4+20*g**2-4/3)))


def frame_value(rotvec, xs, weights):
    R = Rotation.from_rotvec(rotvec).as_matrix()
    coordinates = xs @ R
    return float(32 * np.sum(weights * np.prod(coordinates**2, axis=1)))


def minimize_frame(xs, weights, starts=30, seed=0):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(starts):
        v = Rotation.random(random_state=rng).as_rotvec()
        result = minimize(frame_value, v, args=(xs, weights), method="BFGS", options={"gtol":1e-11,"maxiter":1000})
        if best is None or result.fun < best.fun:
            best = result
    return best.fun, Rotation.from_rotvec(best.x).as_matrix()


def audit(samples=100, seed=33):
    rng = np.random.default_rng(seed)
    worst = None
    for sample in range(samples):
        n = rng.integers(3, 15)
        xs = rng.normal(size=(n, 3))
        xs /= np.linalg.norm(xs, axis=1)[:, None]
        weights = rng.dirichlet(np.ones(n))
        E = energy(xs, weights)
        C, R = minimize_frame(xs, weights, starts=12, seed=sample)
        row = E-C,E,C,xs,weights,R
        if worst is None or row[0]<worst[0]: worst=row
    print("worst gap,E,C",worst[:3])
    return worst


if __name__ == "__main__":
    audit()
