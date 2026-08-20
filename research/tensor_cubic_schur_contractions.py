"""Identify tensor-contraction normalizations in the isotropic cubic Schur core."""

import sys

sys.path.insert(0, "research")

import numpy as np
import sympy as sp

import tensor_schur_feasible_scan as data
import tensor_spin6_gram_kernel as gram
from tensor_hankel_fw import maps


def spin4_tensor(coeff):
    A4, B, C4, _ = data.direction(coeff)
    # Recreate the same normalized pure-spin-4 sextic moment coordinate.
    H4 = sum(c * M for c, M in zip(coeff, data.basis4))
    B0 = data.U.T @ H4 @ data.J
    scale = np.sqrt(np.sum(B0 * B0))
    y = sum(c * yy for c, yy in zip(coeff, data.ybasis4)) / scale
    M2 = np.einsum("abk,k->ab", maps[2], y)
    # data.S for degree 3; build the analogous isometry for degree 2.
    words2 = [(i, j) for i in range(3) for j in range(3)]
    S2 = np.zeros((9, 6))
    for col, alpha in enumerate(data.comps2):
        matches = [w for w in words2 if tuple(w.count(i) for i in range(3)) == alpha]
        for w in matches:
            S2[3 * w[0] + w[1], col] = 1 / np.sqrt(len(matches))
    H = (S2 @ M2 @ S2.T).reshape(3, 3, 3, 3)
    return H, A4, B, C4


def report(seed=11, samples=20):
    rng = np.random.default_rng(seed)
    coeff = rng.normal(size=9)
    H, A4, B, C4 = spin4_tensor(coeff)
    print("trace H", np.max(np.abs(np.einsum("iikl->kl", H))))
    ratios_a = []
    ratios_b = []
    for _ in range(samples):
        v = rng.normal(size=7)
        V = gram.tensor(v)
        raw_a = np.einsum("ijab,ijk,abk", H, V, V)
        ratios_a.append((v @ A4 @ v) / raw_a)
        raw_b = np.einsum("ijkl,ijk->l", H, V)
        ratios_b.extend((B.T @ v / raw_b)[np.abs(raw_b) > 1e-8])
    for name, values in (("A", ratios_a), ("B", ratios_b)):
        value = float(np.median(values))
        print(name, value, sp.nsimplify(value), "spread", np.ptp(values))
    # C4 is direct contraction with the normalized H2 tensor basis.
    print("H norm", np.sum(H * H), "C4 norm", np.sum(C4 * C4))


if __name__ == "__main__":
    report()
