"""Compare the isotropic cubic Schur block with wedge^2(Bhat) on H2.

Representation theory gives wedge^2(V_2)=V_1+V_3.  Since the normalized
group gap Bhat acts on V_2 while the cubic catalecticant splits as V_3+V_1,
their dimensions strongly suggest a compound-matrix identity.  This script
constructs the exact numerical intertwiners and audits each homogeneous term.
"""

import itertools
import sys

sys.path.insert(0, "research")

import numpy as np
from scipy.linalg import null_space

import tensor_schur_feasible_scan as data


def sym_rep(degree):
    words = list(itertools.product(range(3), repeat=degree))
    word_index = {w: i for i, w in enumerate(words)}
    comps = [
        (a, b, degree - a - b)
        for a in range(degree + 1)
        for b in range(degree + 1 - a)
    ]
    S = np.zeros((3**degree, len(comps)))
    for col, alpha in enumerate(comps):
        matches = [
            w for w in words if tuple(w.count(i) for i in range(3)) == alpha
        ]
        for w in matches:
            S[word_index[w], col] = 1 / np.sqrt(len(matches))
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1
    generators = []
    for axis in range(3):
        L = np.array([[eps[i, axis, j] for j in range(3)] for i in range(3)])
        G = np.zeros((3**degree, 3**degree))
        for slot in range(degree):
            factors = [np.eye(3)] * degree
            factors[slot] = L
            term = factors[0]
            for factor in factors[1:]:
                term = np.kron(term, factor)
            G += term
        generators.append(S.T @ G @ S)
    return comps, generators


comps2, gens_sym2 = sym_rep(2)
j2 = np.array(
    [1 / np.sqrt(3) if max(alpha) == 2 else 0.0 for alpha in comps2]
)[:, None]
U2 = null_space(j2.T)
gens2 = [U2.T @ g @ U2 for g in gens_sym2]

pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]


def additive_wedge(G):
    out = np.zeros((10, 10))
    for a, (i, j) in enumerate(pairs):
        for b, (k, ell) in enumerate(pairs):
            out[a, b] = (
                G[i, k] * (j == ell)
                + G[j, ell] * (i == k)
                - G[i, ell] * (j == k)
                - G[j, k] * (i == ell)
            )
    return out


gens_wedge = [additive_wedge(g) for g in gens2]
casimir = -sum(g @ g for g in gens_wedge)
evals, evecs = np.linalg.eigh(casimir)
W1 = evecs[:, np.abs(evals - 2) < 1e-7]
W3 = evecs[:, np.abs(evals - 12) < 1e-7]
gw3 = [W3.T @ g @ W3 for g in gens_wedge]
gw1 = [W1.T @ g @ W1 for g in gens_wedge]


def intertwiner(left, right):
    # Solve left[a] T = T right[a].
    equations = []
    for L, R in zip(left, right):
        equations.append(np.kron(np.eye(R.shape[0]), L) - np.kron(R.T, np.eye(L.shape[0])))
    kernel = null_space(np.vstack(equations))
    T = kernel[:, 0].reshape(left[0].shape[0], right[0].shape[0], order="F")
    scale = np.sqrt(np.trace(T.T @ T) / T.shape[1])
    return T / scale


T33 = intertwiner(gw3, data.g3)
# data.gens are on Sym^3; the trace embedding J carries the standard V1
# action.  Extract it in the orthonormal J basis.
g1_cubic = [data.J.T @ g @ data.J for g in data.gens]
T11 = intertwiner(gw1, g1_cubic)


def compound(B):
    out = np.zeros((10, 10))
    for a, (i, j) in enumerate(pairs):
        for b, (k, ell) in enumerate(pairs):
            out[a, b] = B[i, k] * B[j, ell] - B[i, ell] * B[j, k]
    return out


def h3_block(B):
    block = W3.T @ compound(B) @ W3
    return T33.T @ block @ T33


def all_blocks(B):
    C = compound(B)
    C33 = T33.T @ (W3.T @ C @ W3) @ T33
    C11 = T11.T @ (W1.T @ C @ W1) @ T11
    C31 = T33.T @ (W3.T @ C @ W1) @ T11
    return C33, C31, C11


def report(coeff):
    A4, B, C4, _ = data.direction(coeff)
    Q = B @ B.T
    I5 = np.eye(5)
    W0 = h3_block(I5)
    Wplus = h3_block(I5 - 5 * C4)
    Wminus = h3_block(I5 + 5 * C4)
    W1term = (Wplus - Wminus) / 2
    W2term = (Wplus + Wminus) / 2 - W0
    C0 = all_blocks(I5)
    Cp = all_blocks(I5 - 5 * C4)
    Cm = all_blocks(I5 + 5 * C4)
    linear_blocks = tuple((p - m) / 2 for p, m in zip(Cp, Cm))
    quadratic_blocks = tuple((p + m) / 2 - z for p, m, z in zip(Cp, Cm, C0))
    print("intertwiner errors", [np.linalg.norm(L @ T33 - T33 @ R) for L, R in zip(gw3, data.g3)])
    print("W0 scalar", np.trace(W0) / 7, np.linalg.norm(W0 - np.eye(7)))
    print("linear fit A4", np.sum(W1term * A4) / np.sum(W1term * W1term), "residual", np.linalg.norm(A4 - (np.sum(W1term*A4)/np.sum(W1term*W1term))*W1term))
    print("quadratic fit Q", np.sum(W2term * Q) / np.sum(W2term * W2term), "residual", np.linalg.norm(Q - (np.sum(W2term*Q)/np.sum(W2term*W2term))*W2term))
    print("spin norms W2", [np.linalg.norm(data.proj3(W2term, ell)) for ell in (0, 2, 4, 6)])
    print("spin norms Q", [np.linalg.norm(data.proj3(Q, ell)) for ell in (0, 2, 4, 6)])
    print(
        "quadratic spin ratios W2/Q",
        [
            np.sum(data.proj3(W2term, ell) * data.proj3(Q, ell))
            / np.sum(data.proj3(Q, ell) ** 2)
            for ell in (0, 2, 4, 6)
        ],
    )
    print("cross linear fit B", np.sum(linear_blocks[1]*B)/np.sum(linear_blocks[1]**2), "residual", np.linalg.norm(B-(np.sum(linear_blocks[1]*B)/np.sum(linear_blocks[1]**2))*linear_blocks[1]))
    print("cross quadratic norm", np.linalg.norm(quadratic_blocks[1]))
    print("11 linear norm", np.linalg.norm(linear_blocks[2]), "11 quadratic", quadratic_blocks[2])
    BtB = B.T @ B
    print("11 quadratic fit BtB", np.sum(quadratic_blocks[2]*BtB)/np.sum(BtB*BtB), "residual", np.linalg.norm(quadratic_blocks[2]-(np.sum(quadratic_blocks[2]*BtB)/np.sum(BtB*BtB))*BtB))


if __name__ == "__main__":
    report(np.random.default_rng(3).normal(size=9))
