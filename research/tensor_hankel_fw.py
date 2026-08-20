"""Frank-Wolfe search over normalized PSD ternary sextic catalecticants."""

import itertools
import math

import cvxpy as cp
import numpy as np


def compositions(total, length):
    if length == 1:
        return [(total,)]
    out = []
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            out.append((first,) + tail)
    return out


def multinomial(alpha):
    result = math.factorial(sum(alpha))
    for value in alpha:
        result //= math.factorial(value)
    return result


deg6 = compositions(6, 3)
index6 = {alpha: i for i, alpha in enumerate(deg6)}


def add(*alphas):
    return tuple(sum(values) for values in zip(*alphas))


def matrix_map(degree):
    basis = compositions(degree, 3)
    result = np.zeros((len(basis), len(basis), len(deg6)))
    remainder = 3 - degree
    remainder_basis = compositions(remainder, 3)
    for row, alpha in enumerate(basis):
        for col, beta in enumerate(basis):
            scale = math.sqrt(multinomial(alpha) * multinomial(beta))
            for gamma in remainder_basis:
                # Expansion of (x_1^2+x_2^2+x_3^2)^remainder.
                coefficient = multinomial(gamma)
                exponent = add(alpha, beta, tuple(2 * value for value in gamma))
                result[row, col, index6[exponent]] += scale * coefficient
    return result


maps = {degree: matrix_map(degree) for degree in (1, 2, 3)}


def matrices(y):
    return {degree: np.einsum("abk,k->ab", tensor, y) for degree, tensor in maps.items()}


def energy(y):
    rho = matrices(y)
    return (
        32 * np.sum(rho[3] ** 2)
        - 48 * np.sum(rho[2] ** 2)
        + 20 * np.sum(rho[1] ** 2)
        - 4 / 3
    )


def gradient(y):
    rho = matrices(y)
    return (
        64 * np.einsum("ab,abk->k", rho[3], maps[3])
        - 96 * np.einsum("ab,abk->k", rho[2], maps[2])
        + 40 * np.einsum("ab,abk->k", rho[1], maps[1])
    )


def uniform_moment(alpha):
    if any(value % 2 for value in alpha):
        return 0.0
    # S^2 even moment: product (alpha_i-1)!! / (7!!), total degree 6.
    numerator = 1
    for value in alpha:
        for odd in range(1, value, 2):
            numerator *= odd
    return numerator / 105


y_variable = cp.Variable(len(deg6))
hankel = sum(y_variable[k] * maps[3][:, :, k] for k in range(len(deg6)))
normalizer = np.trace(maps[3], axis1=0, axis2=1)
direction = cp.Parameter(len(deg6))
problem = cp.Problem(
    cp.Minimize(direction @ y_variable),
    [hankel >> 0, normalizer @ y_variable == 1],
)


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    starts = [np.array([uniform_moment(alpha) for alpha in deg6])]
    for _ in range(8):
        direction.value = rng.normal(size=len(deg6))
        problem.solve(solver=cp.CLARABEL)
        starts.append(y_variable.value.copy())

    for start_number, y in enumerate(starts):
        best = energy(y)
        for iteration in range(300):
            direction.value = gradient(y)
            problem.solve(solver=cp.CLARABEL)
            vertex = y_variable.value.copy()
            delta = vertex - y
            # Exact line minimizer for the quadratic energy.
            e0 = energy(y)
            e1 = energy(y + delta)
            em = energy(y + 0.5 * delta)
            curvature = 2 * (e1 + e0 - 2 * em)
            slope = 4 * em - e1 - 3 * e0
            step = np.clip(-slope / (2 * curvature), 0, 1) if curvature > 1e-14 else (1 if e1 < e0 else 0)
            y += step * delta
            best = min(best, energy(y))
            if step < 1e-10:
                break
        eigenvalues = np.linalg.eigvalsh(matrices(y)[3])
        print(start_number, iteration, best, eigenvalues)
