"""Audit the one-ray Schur separator Y=t(-M)_+.

For M+A6 >= 0 and A6 in spin 6, put N=(-M)_+.  Pairing with N gives

    <Proj6 N, A6> >= ||N||^2,

and hence ||A6||^2 >= ||N||^4 / ||Proj6 N||^2.  This file compares that
closed lower bound with E and with the determinant-strengthened conjecture.
"""

import itertools
import sys

sys.path.insert(0, "research")

import numpy as np

import tensor_schur_feasible_scan as data


def negative_part_bound(M):
    eigenvalues, eigenvectors = np.linalg.eigh(M)
    negative = np.maximum(-eigenvalues, 0.0)
    if np.max(negative) < 1e-13:
        return 0.0, eigenvalues, 0.0, 0.0
    N = (eigenvectors * negative) @ eigenvectors.T
    N6 = data.proj3(N, 6)
    n2 = np.sum(N * N)
    n62 = np.sum(N6 * N6)
    lower = n2 * n2 / n62 if n62 > 1e-28 else np.inf
    return lower, eigenvalues, n2, n62


def frame_qp_bound(M):
    """Exact dual bound from each negative eigenvector constraint separately."""
    eigenvalues, eigenvectors = np.linalg.eigh(M)
    indices = np.flatnonzero(eigenvalues < -1e-12)
    if len(indices) == 0:
        return 0.0, eigenvalues, np.empty((0, 0)), np.empty(0)
    a = -eigenvalues[indices]
    Z = [
        data.proj3(np.outer(eigenvectors[:, i], eigenvectors[:, i]), 6)
        for i in indices
    ]
    K = np.array([[np.sum(x * y) for y in Z] for x in Z])
    best = 0.0
    best_eta = np.zeros(len(indices))
    for size in range(1, len(indices) + 1):
        for subset_tuple in itertools.combinations(range(len(indices)), size):
            subset = np.array(subset_tuple)
            Ks = K[np.ix_(subset, subset)]
            try:
                eta = np.linalg.solve(Ks, a[subset])
            except np.linalg.LinAlgError:
                continue
            if np.min(eta) < -1e-10:
                continue
            value = float(a[subset] @ eta)
            if value > best:
                best = value
                best_eta = np.zeros(len(indices))
                best_eta[subset] = eta
    return best, eigenvalues, K, best_eta


def scan(seed=20260820, trials=20000, scales=120, verbose=True):
    rng = np.random.default_rng(seed)
    worst_e = None
    worst_det = None
    worst_ratio = None
    worst_kappa = None
    for trial in range(trials):
        coeff = rng.normal(size=9)
        A4, B, C4, _ = data.direction(coeff)
        Q = B @ B.T
        for b in np.linspace(0.0, 0.08, scales + 1)[1:]:
            t = np.sqrt(b)
            bhat_eigenvalues = np.linalg.eigvalsh(np.eye(5) - 5 * t * C4)
            if bhat_eigenvalues[0] < -1e-10 or bhat_eigenvalues[-1] > 5 / 3 + 1e-10:
                continue
            M = (2 / 35) * np.eye(7) + t * A4 - 5 * b * Q
            lower, meig, n2, n62 = negative_part_bound(M)
            e_lower = 2 / 105 - 5 * b / 11 + 2 * lower
            determinant = np.prod(bhat_eigenvalues)
            det_gap = e_lower - (2 / 105) * determinant**2
            rec = (e_lower, trial, b, lower, determinant, bhat_eigenvalues, meig, n2, n62, coeff)
            rec_det = (det_gap,) + rec[1:]
            if worst_e is None or rec[0] < worst_e[0]:
                worst_e = rec
                if verbose:
                    print("worst E", worst_e[:-1], flush=True)
            if worst_det is None or rec_det[0] < worst_det[0]:
                worst_det = rec_det
                if verbose:
                    print("worst det", worst_det[:-1], flush=True)
            deficit = 5 * b / 11 - 2 / 105
            if deficit > 0:
                ratio = 2 * lower / deficit
                rec_ratio = (ratio,) + rec[1:]
                if worst_ratio is None or ratio < worst_ratio[0]:
                    worst_ratio = rec_ratio
                    if verbose:
                        print("worst ratio", worst_ratio[:-1], flush=True)
            required = (5 * b / 11 - 2 / 105 + (2 / 105) * determinant**2) / 2
            if required > 1e-14 and n2 > 0:
                # Any universal ||P6 N||^2 <= kappa ||N||^2 closes this
                # sample provided kappa <= this threshold.
                kappa_threshold = n2 / required
                rec_kappa = (kappa_threshold, n62 / n2) + rec[1:]
                if worst_kappa is None or kappa_threshold < worst_kappa[0]:
                    worst_kappa = rec_kappa
                    if verbose:
                        print("worst kappa", worst_kappa[:-1], flush=True)
    if not verbose:
        print("FINAL E", None if worst_e is None else worst_e[:-1])
        print("FINAL det", None if worst_det is None else worst_det[:-1])
        print("FINAL ratio", None if worst_ratio is None else worst_ratio[:-1])
        print("FINAL kappa", None if worst_kappa is None else worst_kappa[:-1])
    return worst_e, worst_det, worst_ratio, worst_kappa


if __name__ == "__main__":
    scan()
