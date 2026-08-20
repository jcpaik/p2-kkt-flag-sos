"""Adversarial audit of the independent negative-eigenvector Schur bound."""

import sys

sys.path.insert(0, "research")

import numpy as np
import cvxpy as cp
from scipy.optimize import minimize

import tensor_negative_part_scan as sep
import tensor_schur_feasible_scan as data


def scan(seed=20260821, trials=20000, scales=80):
    rng = np.random.default_rng(seed)
    worst_e = None
    worst_det = None
    max_negative = 0
    for trial in range(trials):
        coeff = rng.normal(size=9)
        A4, B, C4, _ = data.direction(coeff)
        Q = B @ B.T
        # Only the high-spin-purity range can require a degree-six lower bound.
        for b in np.linspace(22 / 525, 0.08, scales + 1)[1:]:
            t = np.sqrt(b)
            bhat = np.linalg.eigvalsh(np.eye(5) - 5 * t * C4)
            if bhat[0] < -1e-10 or bhat[-1] > 5 / 3 + 1e-10:
                continue
            M = (2 / 35) * np.eye(7) + t * A4 - 5 * b * Q
            lower, meig, K, eta = sep.frame_qp_bound(M)
            max_negative = max(max_negative, K.shape[0])
            base = 2 / 105 - 5 * b / 11
            e_gap = base + 2 * lower
            determinant = np.prod(bhat)
            det_gap = e_gap - (2 / 105) * determinant**2
            rec = (e_gap, trial, b, lower, determinant, bhat, meig, K, eta, coeff)
            rec_det = (det_gap,) + rec[1:]
            if worst_e is None or e_gap < worst_e[0]:
                worst_e = rec
            if worst_det is None or det_gap < worst_det[0]:
                worst_det = rec_det
    print("max negative", max_negative)
    print("worst E", None if worst_e is None else worst_e[:-1])
    print("worst det", None if worst_det is None else worst_det[:-1])
    if worst_e is not None:
        print("worst E coeff", repr(worst_e[-1]))
    if worst_det is not None:
        print("worst det coeff", repr(worst_det[-1]))
    return worst_e, worst_det


def boundary_value(coeff, determinant=False, return_record=False):
    A4, B, C4, _ = data.direction(coeff)
    ceig = np.linalg.eigvalsh(C4)
    limits = [np.sqrt(0.08)]
    if ceig[-1] > 1e-12:
        limits.append(1 / (5 * ceig[-1]))
    if ceig[0] < -1e-12:
        limits.append(-2 / (15 * ceig[0]))
    t = min(limits) * (1 - 1e-10)
    b = t * t
    bhat = np.linalg.eigvalsh(np.eye(5) - 5 * t * C4)
    M = (2 / 35) * np.eye(7) + t * A4 - 5 * b * (B @ B.T)
    lower, meig, K, eta = sep.frame_qp_bound(M)
    value = 2 / 105 - 5 * b / 11 + 2 * lower
    determinant_value = np.prod(bhat)
    if determinant:
        value -= (2 / 105) * determinant_value**2
    if return_record:
        return value, b, lower, determinant_value, bhat, meig, K, eta
    return value


def polish(coeff, determinant=False):
    result = minimize(
        lambda x: boundary_value(x, determinant=determinant),
        coeff,
        method="Powell",
        options={"maxiter": 3000, "xtol": 1e-11, "ftol": 1e-13},
    )
    print(result)
    print(boundary_value(result.x, determinant=determinant, return_record=True))
    print("coeff", repr(result.x))
    return result


def compression_bound(M):
    """Min ||A6||^2 subject only to PSD on the negative eigenspace of M."""
    eigenvalues, eigenvectors = np.linalg.eigh(M)
    W = eigenvectors[:, eigenvalues < -1e-11]
    if W.shape[1] == 0:
        return 0.0, eigenvalues
    z = cp.Variable(13)
    A6 = sum(z[i] * data.basis6[i] for i in range(13))
    problem = cp.Problem(cp.Minimize(cp.sum_squares(A6)), [W.T @ (M + A6) @ W >> 0])
    value = problem.solve(solver="CLARABEL")
    return value, eigenvalues


def scan_compression(seed=20260822, trials=300, scales=30):
    rng = np.random.default_rng(seed)
    worst_e = None
    worst_det = None
    for trial in range(trials):
        coeff = rng.normal(size=9)
        A4, B, C4, _ = data.direction(coeff)
        Q = B @ B.T
        for b in np.linspace(22 / 525, 0.08, scales + 1)[1:]:
            t = np.sqrt(b)
            bhat = np.linalg.eigvalsh(np.eye(5) - 5 * t * C4)
            if bhat[0] < -1e-10 or bhat[-1] > 5 / 3 + 1e-10:
                continue
            M = (2 / 35) * np.eye(7) + t * A4 - 5 * b * Q
            lower, meig = compression_bound(M)
            base = 2 / 105 - 5 * b / 11
            e_gap = base + 2 * lower
            determinant = np.prod(bhat)
            det_gap = e_gap - (2 / 105) * determinant**2
            rec = (e_gap, trial, b, lower, determinant, bhat, meig, coeff)
            rec_det = (det_gap,) + rec[1:]
            if worst_e is None or e_gap < worst_e[0]:
                worst_e = rec
            if worst_det is None or det_gap < worst_det[0]:
                worst_det = rec_det
    print("compression E", worst_e[:-1], repr(worst_e[-1]))
    print("compression det", worst_det[:-1], repr(worst_det[-1]))
    return worst_e, worst_det


if __name__ == "__main__":
    scan()
