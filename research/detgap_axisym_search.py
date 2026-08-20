"""Numerical stress test of the determinant gap on axisymmetric measures."""

import argparse

import numpy as np


def legendre_even(y):
    p2 = (3 * y - 1) / 2
    p4 = (35 * y**2 - 30 * y + 3) / 8
    p6 = (231 * y**3 - 315 * y**2 + 105 * y - 5) / 16
    return p2, p4, p6


def values(y, w):
    p2, p4, p6 = legendre_even(y)
    u2, u4, u6 = [np.dot(w, p) for p in (p2, p4, p6)]
    lam0 = 1 + 5 * u2 / 7 - 12 * u4 / 7
    lam1 = 1 + 5 * u2 / 14 + 8 * u4 / 7
    lam2 = 1 - 5 * u2 / 7 - 2 * u4 / 7
    det = lam0 * lam1**2 * lam2**2
    energy = 32 / 105 + 8 * u2**2 / 7 - 384 * u4**2 / 385 + 512 * u6**2 / 231
    return energy - (32 / 105) * det**2, energy, det, (u2, u4, u6), (lam0, lam1, lam2)


def run(n, batches, batch_size, seed):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(batches):
        y = rng.random((batch_size, n))
        w = rng.dirichlet(np.ones(n) / 2, batch_size)
        p2, p4, p6 = legendre_even(y)
        u2 = np.sum(w * p2, axis=1)
        u4 = np.sum(w * p4, axis=1)
        u6 = np.sum(w * p6, axis=1)
        l0 = 1 + 5 * u2 / 7 - 12 * u4 / 7
        l1 = 1 + 5 * u2 / 14 + 8 * u4 / 7
        l2 = 1 - 5 * u2 / 7 - 2 * u4 / 7
        det = l0 * l1**2 * l2**2
        energy = 32 / 105 + 8 * u2**2 / 7 - 384 * u4**2 / 385 + 512 * u6**2 / 231
        gap = energy - (32 / 105) * det**2
        i = int(np.argmin(gap))
        if best is None or gap[i] < best[0]:
            best = (gap[i], y[i], w[i])
    print("best", best[0])
    print("y", best[1])
    print("w", best[2])
    print("values", values(best[1], best[2]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=4)
    parser.add_argument("--batches", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20260820)
    args = parser.parse_args()
    run(args.n, args.batches, args.batch_size, args.seed)
