"""Exact/numerical identities for bosonic reduction residuals."""
import itertools
import math
import numpy as np

from tensor_hankel_fw import deg6, maps, matrices


def comps(n, k=3):
    if k == 1:
        return [(n,)]
    return [(a,) + rest for a in range(n + 1) for rest in comps(n - a, k - 1)]


def sym_embedding(n):
    basis = comps(n)
    U = np.zeros((3 ** n, len(basis)))
    for col, alpha in enumerate(basis):
        words = set(itertools.permutations(sum(([i] * alpha[i] for i in range(3)), [])))
        for word in words:
            row = sum(word[k] * 3 ** (n - 1 - k) for k in range(n))
            U[row, col] = 1 / math.sqrt(len(words))
    return U


U1 = sym_embedding(1)
U2 = sym_embedding(2)
U3 = sym_embedding(3)


def lift32(R2):
    full = np.kron(np.eye(3), U2 @ R2 @ U2.T)
    return U3.T @ full @ U3


def lift21(R1):
    return U2.T @ np.kron(np.eye(3), U1 @ R1 @ U1.T) @ U2


def reduce32(X3):
    full = (U3 @ X3 @ U3.T).reshape(3, 9, 3, 9)
    return U2.T @ np.einsum("iaib->ab", full) @ U2


def reduce21(X2):
    full = (U2 @ X2 @ U2.T).reshape(3, 3, 3, 3)
    return U1.T @ np.einsum("iaib->ab", full) @ U1


def invariants(y):
    rho = matrices(y)
    H, R2, R1 = rho[3], rho[2], rho[1]
    L32, L21 = lift32(R2), lift21(R1)
    D3, D2 = L32 - H, L21 - R2
    mass = np.trace(H)
    return np.array([
        np.sum(H * H), np.sum(R2 * R2), np.sum(R1 * R1),
        mass * mass,
        np.sum(L32 * L32), np.sum(L21 * L21),
        np.sum(H * L32), np.sum(R2 * L21),
        np.trace(D3), np.trace(D2), np.sum(D3 * D3), np.sum(D2 * D2),
    ])


rng = np.random.default_rng(8)
data = np.array([invariants(rng.normal(size=28)) for _ in range(100)])
# Regress each derived quadratic against a,b,c,mass^2; traces separately.
for j, name in zip((4, 5, 6, 7, 10, 11),
                   ("||L32||2", "||L21||2", "<H,L32>", "<R2,L21>", "||D3||2", "||D2||2")):
    coef, *_ = np.linalg.lstsq(data[:, :4], data[:, j], rcond=None)
    print(name, coef, "err", np.max(np.abs(data[:, j] - data[:, :4] @ coef)))
print("sample traces", data[0, 8:10])

# Linear equivariant identities for partial traces of the residuals.
lhs3 = []
rhs3 = []
lhs2 = []
rhs2 = []
for _ in range(40):
    y = rng.normal(size=28)
    H, R2, R1 = matrices(y)[3], matrices(y)[2], matrices(y)[1]
    D3, D2 = lift32(R2) - H, lift21(R1) - R2
    lhs3.append(reduce32(D3).ravel())
    rhs3.append(np.stack((R2, lift21(R1), np.trace(H) * np.eye(6)), axis=-1).reshape(-1, 3))
    lhs2.append(reduce21(D2).ravel())
    rhs2.append(np.stack((R1, np.trace(H) * np.eye(3)), axis=-1).reshape(-1, 2))
coef3, *_ = np.linalg.lstsq(np.concatenate(rhs3), np.concatenate(lhs3), rcond=None)
coef2, *_ = np.linalg.lstsq(np.concatenate(rhs2), np.concatenate(lhs2), rcond=None)
print("Tr D3 coefficients R2,L21,I", coef3)
print("Tr D2 coefficients R1,I", coef2)
