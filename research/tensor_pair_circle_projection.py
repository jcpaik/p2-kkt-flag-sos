"""Audit the pair-projected great-circle square for the weighted J target.

For a nonparallel ordered root pair (x,y), put P=span{x,y}.  Push mu to
the unit circle of P with weight |Pz|^6.  If xi_z is the oriented complex
coordinate of Pz, the exact circle deficit is

    S_P = |int |Pz|^4 xi_z^2 dmu|^2 + |int xi_z^6 dmu|^2 >= 0.

The root average below uses a^2(1-a^2)^m.  The first polynomial value is
m=6; all calculations here are diagnostic, with exact interpolation done in
a separate script.
"""

import sys

import numpy as np

sys.path.insert(0, "research")
import tensor_weighted_f_wedge3_sdp as rooted


def circle_square(xs, weights, x, y):
    """Return the projected circle deficit for one nonparallel root pair."""
    a = float(x @ y)
    s = 1.0 - a * a
    if s < 1e-14:
        return 0.0
    e1 = x
    e2 = (y - a * x) / np.sqrt(s)
    xi = xs @ e1 + 1j * (xs @ e2)
    r2 = np.abs(xi) ** 2
    m2 = np.sum(weights * r2**2 * xi**2)
    m6 = np.sum(weights * xi**6)
    return float(abs(m2) ** 2 + abs(m6) ** 2)


def root_average(xs, weights, m=6):
    value = 0.0
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            a = float(x @ y)
            s = max(0.0, 1.0 - a * a)
            if s < 1e-14:
                continue
            value += (
                weights[i]
                * weights[j]
                * a**2
                * s**m
                * circle_square(xs, weights, x, y)
            )
    return value


def polynomial_root_average(xs, weights, root_power=1):
    """Denominator-cleared m=6 average, evaluated as a four-point polynomial.

    ``root_power=1`` is the a^2 weight proposed originally; ``root_power=2``
    is the a^4 variant whose pole-splitting scaling matches J-108F.
    """
    g = np.asarray(xs) @ np.asarray(xs).T
    a = g[:, :, None, None]
    u = g[:, None, :, None]
    v = g[None, :, :, None]
    r = g[:, None, None, :]
    t = g[None, :, None, :]
    U = u * u + v * v - 2 * a * u * v
    V = r * r + t * t - 2 * a * r * t
    L = u * r + v * t - a * (u * t + v * r)
    UV = U * V
    phi = a ** (2 * root_power) * (
        32 * L**6 - 48 * L**4 * UV + 20 * L**2 * UV**2 - 2 * UV**3
    )
    w4 = np.einsum("i,j,k,l->ijkl", weights, weights, weights, weights)
    return float(np.sum(w4 * phi))


def energy_data(xs, weights):
    gram = xs @ xs.T
    p2 = np.sum(weights[:, None] * weights[None, :] * gram**2)
    p4 = np.sum(weights[:, None] * weights[None, :] * gram**4)
    p6 = np.sum(weights[:, None] * weights[None, :] * gram**6)
    J = 144 * p6 - 216 * p4 + 87 * p2 - 5
    F = rooted.direct_f(xs, weights)
    return J, F, J - 108 * F


def audit(seed=812):
    rng = np.random.default_rng(seed)
    examples = []
    examples.append(("ONB", np.eye(3), np.ones(3) / 3))
    n = 12
    pe = np.vstack(
        [
            np.array([[0.0, 0.0, 1.0]]),
            np.array(
                [
                    [np.cos(2 * np.pi * k / n), np.sin(2 * np.pi * k / n), 0]
                    for k in range(n)
                ]
            ),
        ]
    )
    examples.append(("pole-circle", pe, np.r_[1 / 3, np.ones(n) * 2 / (3 * n)]))
    for k in range(5):
        xs = rng.normal(size=(7, 3))
        xs /= np.linalg.norm(xs, axis=1)[:, None]
        examples.append((f"random-{k}", xs, rng.dirichlet(np.ones(7))))
    for name, xs, weights in examples:
        J, F, gap = energy_data(xs, weights)
        av = root_average(xs, weights)
        print(name, "A6", av, "J", J, "F", F, "J-108F", gap, "ratio", gap / av if av else np.nan)


if __name__ == "__main__":
    audit()
