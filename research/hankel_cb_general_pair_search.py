"""Search mixtures of general real 3x3-line Cayley--Bacharach rays.

Two products of three real projective lines meet in nine real points.  Their
cubic evaluation functionals have a unique Cayley--Bacharach relation.  A
diagonal form with one negative coefficient, chosen from that relation, is
PSD on the eight-dimensional evaluation image and gives a rank-seven Hankel
extreme ray.
"""
import math
import numpy as np


def comps(n, k=3):
    if k == 1:
        return [(n,)]
    return [(a,) + t for a in range(n + 1) for t in comps(n - a, k - 1)]


def multinomial(a):
    out = math.factorial(sum(a))
    for q in a:
        out //= math.factorial(q)
    return out


d3 = comps(3)
d6 = comps(6)
ix = {a: i for i, a in enumerate(d6)}


def moment_map(d):
    basis = comps(d)
    out = np.zeros((len(basis), len(basis), 28))
    remainder = comps(3 - d)
    for i, a in enumerate(basis):
        for j, b in enumerate(basis):
            scale = np.sqrt(multinomial(a) * multinomial(b))
            for g in remainder:
                exponent = tuple(a[q] + b[q] + 2 * g[q] for q in range(3))
                out[i, j, ix[exponent]] += scale * multinomial(g)
    return out


maps = {d: moment_map(d) for d in (1, 2, 3)}
grams = {d: np.einsum("abk,abl->kl", M, M) for d, M in maps.items()}
mass = np.trace(maps[3], axis1=0, axis2=1)
Q = 32 * grams[3] - 48 * grams[2] + 20 * grams[1] - 4 / 3 * np.outer(mass, mass)


def eval_rows(points, degree):
    basis = comps(degree)
    return np.array([[np.prod(p ** np.array(a)) for a in basis] for p in points])


def sample_ray(rng):
    while True:
        # Row normalization removes irrelevant line scalings.  Independent
        # random triples are in general position almost surely.
        left = rng.normal(size=(3, 3))
        right = rng.normal(size=(3, 3))
        left /= np.linalg.norm(left, axis=1)[:, None]
        right /= np.linalg.norm(right, axis=1)[:, None]
        raw = np.array([np.cross(l, r) for l in left for r in right])
        norms = np.linalg.norm(raw, axis=1)
        if norms.min() < 2e-3:
            continue
        points = raw / norms[:, None]
        V = eval_rows(points, 3)
        _, singular, U_t = np.linalg.svd(V.T, full_matrices=True)
        # V has row rank eight: its right/left nullspace is the CB relation.
        u = U_t[-1]
        if singular[-1] > 2e-8 or np.min(np.abs(u)) < 1e-6:
            continue
        neg = int(rng.integers(9))
        inds = [i for i in range(9) if i != neg]
        positive = np.exp(rng.normal(scale=1.8, size=8))
        weights = np.zeros(9)
        weights[inds] = positive
        weights[neg] = -u[neg] ** 2 / np.sum(u[inds] ** 2 / positive)
        total = weights.sum()
        if total <= 0 or -weights[neg] / total < 5e-4:
            continue
        weights /= total
        y = weights @ eval_rows(points, 6)
        return y, (left, right, neg, weights, points, u)


rng = np.random.default_rng(20260819)
N = 20000
Y = []
meta = []
for _ in range(N):
    y, data = sample_ray(rng)
    Y.append(y)
    meta.append(data)
Y = np.asarray(Y)
q = np.einsum("bi,ij,bj->b", Y, Q, Y)
print("self min/max", q.min(), q.max(), "index", int(np.argmin(q)))

best = (np.inf, None, None, None)
for lo in range(0, N, 250):
    B = Y[lo : lo + 250] @ Q @ Y.T
    score = B + np.sqrt(np.maximum(q[lo : lo + 250, None] * q[None, :], 0))
    for k in range(score.shape[0]):
        score[k, lo + k] = np.inf
    a, b = np.unravel_index(np.argmin(score), score.shape)
    if score[a, b] < best[0]:
        best = (score[a, b], lo + a, b, B[a, b])

_, i, j, cross = best
print("pair", i, j, "Lorentz score", best[0], "q", q[i], q[j], "cross", cross)
d = Y[j] - Y[i]
aa = d @ Q @ d
bb = 2 * Y[i] @ Q @ d
t = np.clip(-bb / (2 * aa), 0, 1) if aa > 0 else (1 if q[j] < q[i] else 0)
mixture = Y[i] + t * d
print("mix", t, mixture @ Q @ mixture)
print("meta1", meta[i])
print("meta2", meta[j])
