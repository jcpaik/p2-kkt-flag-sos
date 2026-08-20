"""Exact/numerical audit of the alpha=6 projective-plus-weight variation.

At a KKT measure, for S in H2 use the projective velocity
v_S(x)=Sx-(x^T S x)x and the centered mass velocity
h_S(x)=alpha(x^T S x-tr(SM)).  The second variation is a quadratic form
Q_alpha(S).  At alpha=6 its trace kernel drops from degree eight to degree
six.  This script constructs the full 5 by 5 quadratic form.
"""

import sys

sys.path.insert(0, "research")

import numpy as np

from tensor_fermionic_general_relaxation_opt import (
    F_LINEAR,
    H2_MATRICES,
    T,
    tangent_data,
)
from tensor_fermionic_relaxation_opt import GAMMA


def kernel(t):
    return 32 * t**6 - 48 * t**4 + 20 * t**2 - 4 / 3


def kernel_prime(t):
    return 192 * t**5 - 192 * t**3 + 40 * t


def kernel_second(t):
    return 960 * t**4 - 576 * t**2 + 40


def polarized_pair(x, y, S, T, alpha=6):
    """Polarized pair integrand for the KKT second-variation form."""
    a = x @ y

    def data(A):
        sx = x @ A @ x
        sy = y @ A @ y
        vx = A @ x - sx * x
        vy = A @ y - sy * y
        da = vx @ y + x @ vy
        return sx, sy, vx, vy, da

    sx, sy, vx, vy, ds = data(S)
    tx, ty, wx, wy, dt = data(T)
    positional = (
        0.5 * kernel_second(a) * ds * dt
        + 0.5 * kernel_prime(a) * (vx @ wy + wx @ vy)
        - 0.5
        * a
        * kernel_prime(a)
        * (vx @ wx + vy @ wy)
    )
    mixed = 0.5 * alpha * kernel_prime(a) * (
        (sx + sy) * dt + (tx + ty) * ds
    )
    mass = 0.5 * alpha**2 * kernel(a) * (sx * ty + tx * sy)
    return positional + mixed + mass


def variation_matrix(points, weights, alpha=6):
    points = np.asarray(points)
    weights = np.asarray(weights)
    M = np.einsum("n,ni,nj->ij", weights, points, points)
    E = sum(
        weights[i] * weights[j] * kernel(points[i] @ points[j])
        for i in range(len(points))
        for j in range(len(points))
    )
    out = np.zeros((5, 5))
    means = np.array([np.trace(S @ M) for S in H2_MATRICES])
    for a, S in enumerate(H2_MATRICES):
        for b, T in enumerate(H2_MATRICES):
            out[a, b] = sum(
                weights[i]
                * weights[j]
                * polarized_pair(points[i], points[j], S, T, alpha)
                for i in range(len(points))
                for j in range(len(points))
            )
            out[a, b] -= alpha**2 * E * means[a] * means[b]
    return out, E, M


def uncentered_variation_matrix(points, weights, alpha=6):
    """The pair part before subtracting alpha^2 E m m^T."""
    out, E, M = variation_matrix(points, weights, alpha)
    means = np.array([np.trace(S @ M) for S in H2_MATRICES])
    return out + alpha**2 * E * np.outer(means, means), E, M


def fermionic_moments(points, weights):
    block = sum(w * tangent_data(x)[1] for x, w in zip(points, weights))
    G = T @ block @ T.T
    F = (GAMMA.numpy() @ G.reshape(-1)).reshape(5, 5)
    contraction_g2 = (GAMMA.numpy() @ (G @ G).reshape(-1)).reshape(5, 5)
    return F, G, contraction_g2


def fit_simple_fermionic_formula(samples=100, seed=10):
    """Diagnostic fit to low-complexity matrices built from F and G."""
    rng = np.random.default_rng(seed)
    design = []
    response = []
    for _ in range(samples):
        n = 12
        points = rng.normal(size=(n, 3))
        points /= np.linalg.norm(points, axis=1)[:, None]
        weights = rng.dirichlet(np.ones(n))
        Q, E, M = uncentered_variation_matrix(points, weights)
        F, G, contraction_g2 = fermionic_moments(points, weights)
        candidates = [
            np.eye(5),
            F,
            F @ F,
            contraction_g2,
            np.trace(F @ F) * np.eye(5),
            np.trace(G @ G) * np.eye(5),
            E * np.eye(5),
        ]
        design.append(np.stack([x.reshape(-1) for x in candidates], axis=1))
        response.append(Q.reshape(-1))
    design = np.concatenate(design)
    response = np.concatenate(response)
    coeff, *_ = np.linalg.lstsq(design, response, rcond=None)
    residual = np.linalg.norm(design @ coeff - response) / np.linalg.norm(response)
    return coeff, residual


def fusion_cap(points, weights):
    F = sum(w * tangent_data(x)[0] for x, w in zip(points, weights))
    return (2 / 3) * np.eye(5) - F, F


def examples():
    onb = np.eye(3)
    rows = [("onb", onb, np.ones(3) / 3)]
    # Pole plus an eight-line equator, an exact zero measure.
    theta = np.arange(8) * np.pi / 8
    pole_eq = np.vstack([np.array([[0.0, 0.0, 1.0]]), np.c_[np.cos(theta), np.sin(theta), np.zeros(8)]])
    rows.append(("pole-equator", pole_eq, np.r_[1 / 3, np.ones(8) / 12]))
    for name, points, weights in rows:
        Q, E, M = variation_matrix(points, weights)
        C, F = fusion_cap(points, weights)
        print(name, "E", E)
        print("Q eig", np.linalg.eigvalsh(Q))
        print("cap eig", np.linalg.eigvalsh(C))
        print("commutator", np.linalg.norm(Q @ C - C @ Q))
        print("Q", Q)
        print("cap", C)


if __name__ == "__main__":
    examples()
    print("simple fermionic fit", fit_simple_fermionic_formula())
