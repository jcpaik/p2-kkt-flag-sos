import itertools

import numpy as np
from scipy.linalg import null_space


# Orthonormal monomial basis of Sym^3(R^3): x^alpha/sqrt(alpha!/3!).
alphas = [(a, b, 3 - a - b) for a in range(4) for b in range(4 - a)]
index = {a: i for i, a in enumerate(alphas)}


def generator(axis):
    """Infinitesimal rotation acting on homogeneous cubics."""
    # Vector-field generators: (e_axis cross x).grad.
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
    eps[1, 0, 2] = eps[2, 1, 0] = eps[0, 2, 1] = -1
    raw = np.zeros((10, 10))
    facts = [1, 1, 2, 6]
    for col, alpha in enumerate(alphas):
        for j in range(3):
            if alpha[j] == 0:
                continue
            for k in range(3):
                coefficient = eps[axis, k, j] * alpha[j]
                if coefficient == 0:
                    continue
                beta = list(alpha)
                beta[j] -= 1
                beta[k] += 1
                beta = tuple(beta)
                # Convert between normalized monomial bases.
                norm_ratio = np.sqrt(
                    np.prod([facts[beta[q]] for q in range(3)])
                    / np.prod([facts[alpha[q]] for q in range(3)])
                )
                raw[index[beta], col] += coefficient * norm_ratio
    return raw


# Laplacian/trace map Sym^3 -> Sym^1. Its nullspace is H_3.
trace = np.zeros((3, 10))
facts = [1, 1, 2, 6]
for col, alpha in enumerate(alphas):
    for j in range(3):
        if alpha[j] >= 2:
            beta = list(alpha)
            beta[j] -= 2
            out = beta.index(1)
            trace[out, col] += alpha[j] * (alpha[j] - 1) / np.sqrt(
                6 / np.prod([facts[q] for q in alpha])
            )

H = null_space(trace)
gens = [H.T @ generator(i) @ H for i in range(3)]
assert all(np.linalg.norm(g + g.T) < 1e-10 for g in gens)


def casimir(X):
    return -sum(g @ (g @ X - X @ g) - (g @ X - X @ g) @ g for g in gens)


def spin_projection(X, ell):
    out = X.copy()
    for j in (0, 2, 4, 6):
        if j != ell:
            out = (casimir(out) - j * (j + 1) * out) / (
                ell * (ell + 1) - j * (j + 1)
            )
    return out


weights = {0: 32 / 15, 2: 6, 4: -64 / 15, 6: 128 / 25}


def qform(X):
    return sum(weights[l] * np.sum(spin_projection(X, l) ** 2) for l in weights)


rng = np.random.default_rng(3)
best = (1e9, None)
for rank in range(1, 1):
    for _ in range(500):
        A = rng.normal(size=(7, rank))
        X = A @ A.T
        X /= np.trace(X)
        q = qform(X)
        if q < best[0]:
            best = (q, np.linalg.eigvalsh(X))
    print(rank, best)


if __name__ == "__main__":
    from scipy.optimize import minimize

    def unpack(z, rank):
        A = z.reshape(7, rank)
        X = A @ A.T
        return X / np.trace(X)

    for rank in range(2, 8):
        opt_best = (1e9, None)
        for _ in range(4):
            z0 = rng.normal(size=7 * rank)
            result = minimize(
                lambda z: qform(unpack(z, rank)),
                z0,
                method="BFGS",
                options={"maxiter": 1500, "gtol": 1e-10},
            )
            if result.fun < opt_best[0]:
                X = unpack(result.x, rank)
                opt_best = (result.fun, np.linalg.eigvalsh(X))
        print("opt", rank, opt_best)
