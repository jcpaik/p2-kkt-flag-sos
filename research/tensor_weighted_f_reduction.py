"""Test the canonical S5 <= Sym(S3 tensor I tensor I) reduction for J>=108F."""

import itertools
import math
import sys

sys.path.insert(0, "research")

import numpy as np

import tensor_schur_feasible_scan as rep
import tensor_weighted_f_wedge3_sdp as rooted


def symmetric_embedding(degree):
    words = list(itertools.product(range(3), repeat=degree))
    index = {word: i for i, word in enumerate(words)}
    comps = [
        (a, b, degree - a - b)
        for a in range(degree + 1)
        for b in range(degree + 1 - a)
    ]
    S = np.zeros((3**degree, len(comps)))
    for column, alpha in enumerate(comps):
        matches = [
            word for word in words if tuple(word.count(i) for i in range(3)) == alpha
        ]
        for word in matches:
            S[index[word], column] = 1 / np.sqrt(len(matches))
    return comps, S


COMPS5, SYM5 = symmetric_embedding(5)


def feature5(x):
    return np.array(
        [
            np.sqrt(
                math.factorial(5)
                / np.prod([math.factorial(entry) for entry in alpha])
            )
            * np.prod(x ** np.array(alpha))
            for alpha in COMPS5
        ]
    )


def overlap5(t):
    return t * (2 * t * t - 1) * (3 * t * t - 1) / 2


def fit_filter(samples=500, seed=410):
    rng = np.random.default_rng(seed)
    anchor = rng.normal(size=3)
    anchor /= np.linalg.norm(anchor)
    wa = rooted.tangent_wedge(anchor)
    domain, target = [], []
    for _ in range(samples):
        x = rng.normal(size=3)
        x /= np.linalg.norm(x)
        w = rooted.tangent_wedge(x)
        expected = overlap5(anchor @ x)
        actual = wa @ w
        if abs(expected) < 1e-8:
            continue
        w *= np.sign(expected * actual)
        domain.append(feature5(x))
        target.append(w)
    domain, target = np.stack(domain), np.stack(target)
    linear, *_ = np.linalg.lstsq(domain, target, rcond=None)
    residual = np.linalg.norm(domain @ linear - target) / np.linalg.norm(target)
    return linear.T, residual


FILTER, FILTER_RESIDUAL = fit_filter()


def upper_point(x):
    P = np.outer(x, x)
    full = np.kron(np.kron(np.kron(np.kron(P, P), P), np.eye(3)), np.eye(3))
    B = SYM5.T @ full @ SYM5
    return FILTER @ B @ FILTER.T


def state_data(xs, weights):
    rho3 = np.zeros((10, 10))
    W = np.zeros((21, 21))
    upper = np.zeros((21, 21))
    y6 = np.zeros(28)
    for x, weight in zip(xs, weights):
        v3 = np.array(
            [
                np.sqrt(
                    math.factorial(3)
                    / np.prod([math.factorial(entry) for entry in alpha])
                )
                * np.prod(x ** np.array(alpha))
                for alpha in rep.comps
            ]
        )
        w = rooted.tangent_wedge(x)
        rho3 += weight * np.outer(v3, v3)
        W += weight * np.outer(w, w)
        upper += weight * upper_point(x)
        y6 += weight * rooted.degree_six_evaluation(x)
    A = rep.U.T @ rho3 @ rep.U
    C2 = rooted.compound(A, np.array(rooted.PAIRS7))
    F = (16 / 9) * np.sum(C2 * W)
    bound = (16 / 9) * np.sum(C2 * upper)
    J = rooted.j_value(y6)
    return rho3, W, upper, C2, F, bound, J


def audit(samples=100, seed=618):
    print("filter residual", FILTER_RESIDUAL)
    print("filter singular values", np.linalg.svd(FILTER, compute_uv=False))
    rng = np.random.default_rng(seed)
    worst = None
    for _ in range(samples):
        xs = rng.normal(size=(8, 3))
        xs /= np.linalg.norm(xs, axis=1)[:, None]
        weights = rng.dirichlet(np.ones(len(xs)))
        rho3, W, upper, C2, F, bound, J = state_data(xs, weights)
        record = J / 108 - bound, J / 108 - F, F, bound, J
        if worst is None or record[0] < worst[0]:
            worst = record
        if np.linalg.eigvalsh(upper - W)[0] < -2e-8:
            raise RuntimeError("reduction failed")
    print("worst J/108-upper, J/108-F, F, upper, J", worst)


if __name__ == "__main__":
    audit()
