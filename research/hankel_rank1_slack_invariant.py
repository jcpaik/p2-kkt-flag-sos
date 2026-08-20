"""Invariant quartic for the rank-one KKT slack branch.

If the normalized KKT slack is q=U-E=s f^2, with f a ternary cubic,
harmonic coefficient matching and complementarity give

    s = -N/D,       E = 32/105 + N^2/D,

where N=int_S2 f^2 and

    D=sum_{l=2,4,6} ||Pi_l(f^2)||^2/c_l,
    (c_2,c_4,c_6)=(8/35,-128/1155,512/3003).

Thus the negative-D branch is harmless exactly when D/N^2 <= -105/32.
This file constructs the quartic D exactly (rational matrices) and provides
an optimizer used to identify all candidate critical orbits.
"""

import math
from fractions import Fraction

import numpy as np
import scipy.optimize as so
import sympy as sp


def compositions(n, k=3):
    if k == 1:
        return [(n,)]
    return [(a,) + tail for a in range(n + 1) for tail in compositions(n - a, k - 1)]


C3 = compositions(3)
C6 = compositions(6)
I6 = {a: i for i, a in enumerate(C6)}


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sphere_moment(alpha):
    if any(a % 2 for a in alpha):
        return Fraction(0)
    numerator = 1
    for a in alpha:
        for odd in range(1, a, 2):
            numerator *= odd
    denominator = 1
    for odd in range(1, sum(alpha) + 2, 2):
        denominator *= odd
    return Fraction(numerator, denominator)


def legendre_coefficients(ell):
    t = sp.symbols("t")
    poly = sp.Poly(sp.legendre(ell, t), t)
    return {power[0]: Fraction(int(value.p), int(value.q)) for power, value in poly.terms()}


def projection_gram(ell):
    """B_l with g^T B_l g=||Pi_l g||^2 for a degree-six g."""
    B = [[Fraction(0) for _ in C6] for _ in C6]
    for power, coefficient in legendre_coefficients(ell).items():
        for delta in compositions(power):
            multinomial = math.factorial(power)
            for d in delta:
                multinomial //= math.factorial(d)
            factor = (2 * ell + 1) * coefficient * multinomial
            for i, alpha in enumerate(C6):
                mx = sphere_moment(add(alpha, delta))
                if not mx:
                    continue
                for j, beta in enumerate(C6):
                    my = sphere_moment(add(beta, delta))
                    if my:
                        B[i][j] += factor * mx * my
    return B


B_EXACT = {ell: projection_gram(ell) for ell in (0, 2, 4, 6)}
C = {2: Fraction(8, 35), 4: -Fraction(128, 1155), 6: Fraction(512, 3003)}
D_EXACT = [
    [sum(B_EXACT[ell][i][j] / C[ell] for ell in (2, 4, 6)) for j in range(28)]
    for i in range(28)
]
B = {ell: np.array(B_EXACT[ell], dtype=float) for ell in B_EXACT}
D_MATRIX = np.array(D_EXACT, dtype=float)
MOMENT6 = np.array([float(sphere_moment(a)) for a in C6])

# g=SQUARE_MAP @ vec(f f^T), using row-major vectorization.  The derivative
# is assembled explicitly below because it is only 28 by 10.
PRODUCT_INDEX = np.array([[I6[add(a, b)] for b in C3] for a in C3])
GRAM3 = np.array(
    [[float(sphere_moment(add(a, b))) for b in C3] for a in C3]
)


def square_coefficients(f):
    g = np.zeros(28)
    for i, alpha in enumerate(C3):
        for j, beta in enumerate(C3):
            g[I6[add(alpha, beta)]] += f[i] * f[j]
    return g


def square_jacobian(f):
    out = np.zeros((28, 10))
    for i in range(10):
        for j in range(10):
            out[PRODUCT_INDEX[i, j], i] += f[j]
            out[PRODUCT_INDEX[i, j], j] += f[i]
    return out


def value(f, details=False):
    g = square_coefficients(f)
    n = MOMENT6 @ g
    norms = {ell: g @ B[ell] @ g for ell in (0, 2, 4, 6)}
    d = g @ D_MATRIX @ g
    r = d / n**2
    if details:
        return r, n, d, norms
    return r


def normalize(f):
    n = value(f, details=True)[1]
    return f / np.sqrt(n)


def critical_search(restarts=300, seed=7321):
    rng = np.random.default_rng(seed)
    records = []
    # Optimize R on the Euclidean unit sphere.  Scale invariance makes this
    # equivalent to N=1 and avoids a nonlinear equality constraint.
    for _ in range(restarts):
        x0 = rng.normal(size=10)
        x0 /= np.linalg.norm(x0)

        def objective(x):
            return value(x / np.linalg.norm(x))

        result = so.minimize(objective, x0, method="BFGS", options={"maxiter": 3000, "gtol": 1e-11})
        x = normalize(result.x)
        r, n, d, norms = value(x, details=True)
        records.append((r, np.linalg.norm(result.jac), x, norms, result.success))
    records.sort(key=lambda row: row[0])
    clusters = []
    for row in records:
        if not any(abs(row[0] - old[0]) < 1e-7 for old in clusters):
            clusters.append(row)
    return clusters


def critical_root_search(restarts=500, seed=8128):
    """Find stationary orbits from the exact Lagrange equations.

    We solve grad D=lambda grad N and N=1.  SO(3) makes each solution a
    three-dimensional orbit, so least_squares (rather than a square Newton
    solve) is deliberately used with its rank-deficient Jacobian support.
    """
    rng = np.random.default_rng(seed)
    records = []

    def equations(z):
        f, lam = z[:10], z[10]
        g = square_coefficients(f)
        jac = square_jacobian(f)
        n = MOMENT6 @ g
        d = g @ D_MATRIX @ g
        grad_d = 2 * jac.T @ D_MATRIX @ g
        grad_n = 2 * GRAM3 @ f
        return np.r_[grad_d - lam * grad_n, n - 1]

    for _ in range(restarts):
        f = normalize(rng.normal(size=10))
        z0 = np.r_[f, 2 * value(f)]
        result = so.least_squares(
            equations,
            z0,
            xtol=1e-13,
            ftol=1e-13,
            gtol=1e-13,
            max_nfev=5000,
        )
        residual = np.linalg.norm(equations(result.x))
        if residual > 2e-7:
            continue
        f = normalize(result.x[:10])
        r, _, _, norms = value(f, details=True)
        # Reconstruct the pseudo-moment middle catalecticant.  Raw monomial
        # coordinates are sufficient for its inertia.
        H = GRAM3.copy()
        if r < -1e-10:
            g = square_coefficients(f)
            s = -1 / r
            for i in range(10):
                for j in range(10):
                    H[i, j] += s * D_MATRIX[PRODUCT_INDEX[i, j]] @ g
        eig = np.linalg.eigvalsh(H)
        row = (r, residual, eig, f, norms)
        if not any(abs(r - old[0]) < 2e-7 and np.sum(eig < -1e-7) == np.sum(old[2] < -1e-7) for old in records):
            records.append(row)
    records.sort(key=lambda row: row[0])
    return records


if __name__ == "__main__":
    print("projection ranks", {ell: np.linalg.matrix_rank(B[ell], tol=1e-10) for ell in B})
    rng = np.random.default_rng(4)
    for _ in range(3):
        f = normalize(rng.normal(size=10))
        print("sample", value(f, details=True))
    for row in critical_search():
        print("critical", row[0], "grad", row[1], "norms", row[3], "success", row[4])
        print(row[2])
    for row in critical_root_search():
        print("root", row[0], "res", row[1], "inertia", np.linalg.eigvalsh(np.diag(row[2])))
        print("H eig", row[2], "norms", row[4])
        print(row[3])
