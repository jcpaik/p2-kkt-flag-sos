"""Search positive group-correlation inequalities for the Hankel target.

For H>=0 and a cubic p, f_p(R)=<R p,H R p> is nonnegative.  Hence its
right-translation autocorrelation, averaged over a rotation conjugacy class,
is nonnegative.  On the unique Hankel spin-l sectors this gives

  sum_l ||H_l||^2 ||Proj_l(pp^T)||^2 chi_l(theta)/(2l+1)^2 >= 0.

We test whether positive combinations of these exact inequalities, plus
ordinary harmonic squares, reproduce the target.
"""

import sys

sys.path.insert(0, "research")

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog

import tensor_schur_feasible_scan as rep
from tensor_hankel_compound_sdp import TARGET
from tensor_hankel_fw import maps


SPINS = (0, 2, 4, 6)
DIMENSIONS = np.array([2 * ell + 1 for ell in SPINS])
VFULL = np.stack([maps[3][:, :, k].reshape(-1) for k in range(28)], axis=1)
VFULL_PINV = np.linalg.pinv(VFULL)


def hankel_spin_basis(ell):
    raw = [rep.proj_full(maps[3][:, :, k], ell).reshape(-1) for k in range(28)]
    u, singular, _ = np.linalg.svd(np.stack(raw, axis=1), full_matrices=False)
    return [u[:, i].reshape(10, 10) for i in range(2 * ell + 1)]


BASES = {ell: hankel_spin_basis(ell) for ell in SPINS}


def target_coefficients():
    out = []
    for ell in SPINS:
        matrix = BASES[ell][0]
        y = VFULL_PINV @ matrix.reshape(-1)
        out.append(y @ TARGET @ y)
    return np.array(out)


TARGET_COEFFICIENTS = target_coefficients()


def character(ell, theta):
    if abs(np.sin(theta / 2)) < 1e-10:
        return 2 * ell + 1
    return np.sin((ell + 0.5) * theta) / np.sin(theta / 2)


def generator(p, theta):
    p = p / np.linalg.norm(p)
    A = np.outer(p, p)
    projection_norms = np.array(
        [sum(np.sum(A * B) ** 2 for B in BASES[ell]) for ell in SPINS]
    )
    chars = np.array([character(ell, theta) for ell in SPINS])
    return projection_norms * chars / DIMENSIONS**2


def search(samples=20000, angles=300, seed=4):
    rng = np.random.default_rng(seed)
    generators = []
    metadata = []
    theta_grid = np.linspace(1e-4, 2 * np.pi - 1e-4, angles)
    for sample in range(samples):
        p = rng.normal(size=10)
        p /= np.linalg.norm(p)
        A = np.outer(p, p)
        norms = np.array(
            [sum(np.sum(A * B) ** 2 for B in BASES[ell]) for ell in SPINS]
        )
        for theta in theta_grid:
            chars = np.array([character(ell, theta) for ell in SPINS])
            row = norms * chars / DIMENSIONS**2
            if row[2] < -1e-12:
                generators.append(row)
                metadata.append((sample, theta, p.copy(), norms.copy()))
    G = np.array(generators).T
    positive = [0, 1, 3]
    result = linprog(
        G[2],
        A_ub=G[positive],
        b_ub=TARGET_COEFFICIENTS[positive],
        bounds=(0, None),
        method="highs",
    )
    print("target coefficients", TARGET_COEFFICIENTS)
    print("LP", result.success, result.fun, "needed", TARGET_COEFFICIENTS[2])
    if result.success:
        active = np.where(result.x > 1e-8)[0]
        for i in active:
            sample, theta, p, norms = metadata[i]
            print("active", result.x[i], "theta", theta, "p", p, "norms", norms, "gen", G[:, i])
    return result, metadata, G


if __name__ == "__main__":
    print("basis dimensions", {ell: len(BASES[ell]) for ell in SPINS})
    search(samples=3000, angles=160)
