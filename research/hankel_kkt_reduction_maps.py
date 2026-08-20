"""Structural reduction-map diagnostics for the all-rank J KKT system.

This uses an L2(S^2)-orthonormal cubic basis for the self-consistent
operator C_p and the normalized symmetric-tensor basis for bosonic partial
traces.  Every constructed map is exact at the representation-theoretic
level; floating point is used here only to identify candidate identities.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import fractional_matrix_power
from scipy.optimize import nnls

from hankel_rank1_slack_invariant import (
    B_EXACT,
    C3,
    C6,
    I6,
    add,
    sphere_moment,
)
from hankel_reduction_identities import U3, lift21, lift32, reduce21, reduce32
from tensor_hankel_fw import maps


def multinomial(alpha: tuple[int, ...]) -> int:
    import math

    out = math.factorial(sum(alpha))
    for entry in alpha:
        out //= math.factorial(entry)
    return out


# tensor_hankel_fw and hankel_rank1_slack_invariant use the same graded-lex
# ordering (first exponent increasing, then second exponent increasing).
assert C3 == [(a, b, 3 - a - b) for a in range(4) for b in range(4 - a)]

M = np.array(
    [[float(sphere_moment(add(a, b))) for b in C3] for a in C3], dtype=float
)
M_HALF = fractional_matrix_power(M, 0.5)
M_INV_HALF = fractional_matrix_power(M, -0.5)
J_COEFFICIENTS = {2: 22 / 35, 4: -192 / 385, 6: 768 / 1001}
D6 = np.array(
    [
        [
            sum(float(B_EXACT[ell][i][j]) / J_COEFFICIENTS[ell] for ell in (2, 4, 6))
            for j in range(28)
        ]
        for i in range(28)
    ]
)


def h1_projector() -> np.ndarray:
    raw = np.zeros((10, 3))
    for coordinate in range(3):
        for radial in range(3):
            alpha = [0, 0, 0]
            alpha[coordinate] += 1
            alpha[radial] += 2
            raw[C3.index(tuple(alpha)), coordinate] += 1
    columns = M_HALF @ raw
    q, _ = np.linalg.qr(columns)
    return q @ q.T


PI1 = h1_projector()
PI3 = np.eye(10) - PI1


def block_coordinates(S: np.ndarray):
    # Work in an adapted orthonormal basis only locally.
    values, vectors = np.linalg.eigh(PI1)
    U1 = vectors[:, values > 0.5]
    U3 = vectors[:, values < 0.5]
    A = U1.T @ S @ U1
    X = U1.T @ S @ U3
    B = U3.T @ S @ U3
    return U1, U3, A, X, B


def elementary_positive_maps(S: np.ndarray) -> list[np.ndarray]:
    U1, U3, A, X, B = block_coordinates(S)
    a, b = np.trace(A), np.trace(B)
    out = []
    for block1, cross, block3 in (
        (a * np.eye(3) - A, np.zeros_like(X), np.zeros((7, 7))),
        (np.zeros((3, 3)), np.zeros_like(X), b * np.eye(7) - B),
        (b * np.eye(3), X, a * np.eye(7)),
        (b * np.eye(3), -X, a * np.eye(7)),
    ):
        out.append(
            U1 @ block1 @ U1.T
            + U1 @ cross @ U3.T
            + U3 @ cross.T @ U1.T
            + U3 @ block3 @ U3.T
        )
    return out


def search_elementary_identity(samples: int = 300, seed: int = 91) -> None:
    """Test whether the four elementary positive maps close the sharp bound."""
    rng = np.random.default_rng(seed)
    rows = []
    target = []
    alpha = 35 / 48
    for _ in range(samples):
        Z = rng.normal(size=(10, 10))
        S = Z @ Z.T
        N = np.trace(S)
        C = cp_operator(S)
        D = np.sum(S * C)
        P = N * C - D * np.eye(10)
        rows.append([np.sum(P * value) for value in elementary_positive_maps(S)])
        target.append(N * (-D - alpha * N * N))
    rows = np.asarray(rows)
    target = np.asarray(target)
    coefficients, residual = nnls(rows, target)
    least_squares, *_ = np.linalg.lstsq(rows, target, rcond=None)
    print("elementary nnls", coefficients, "relative residual", residual / np.linalg.norm(target))
    print(
        "elementary ls",
        least_squares,
        "relative residual",
        np.linalg.norm(rows @ least_squares - target) / np.linalg.norm(target),
    )


def product_coefficients(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    out = np.zeros(len(C6))
    for i, alpha in enumerate(C3):
        for j, beta in enumerate(C3):
            out[I6[add(alpha, beta)]] += left[i] * right[j]
    return out


# e = M^{-1/2} m is the orthonormal vector of cubic functions, where m is
# the raw cubic monomial vector.  products[i,j] is the sextic coefficient
# vector of e_i e_j.
products = np.empty((10, 10, 28))
for i in range(10):
    for j in range(10):
        products[i, j] = product_coefficients(M_INV_HALF[i], M_INV_HALF[j])


def cp_operator(S: np.ndarray) -> np.ndarray:
    """Return C_p for J, with p=e^T S e and N=tr(S)."""
    p = np.einsum("ij,ijk->k", S, products)
    return np.einsum("ijk,kl,l->ij", products, D6, p)


def to_tensor_gram(S: np.ndarray) -> np.ndarray:
    """Gram matrix of the same SOS in normalized Sym^3 tensor coordinates."""
    # normalized tensor monomials v_alpha=sqrt(multinomial(alpha))*m_alpha
    scale = np.diag([np.sqrt(multinomial(alpha)) for alpha in C3])
    # e=M^{-1/2}m=M^{-1/2}scale^{-1}v.
    transform = M_INV_HALF @ np.linalg.inv(scale)
    return transform.T @ S @ transform


def functional_to_tensor(P: np.ndarray) -> np.ndarray:
    """Moment matrix of a functional from its orthonormal cubic matrix."""
    scale = np.diag([np.sqrt(multinomial(alpha)) for alpha in C3])
    transform = M_INV_HALF @ np.linalg.inv(scale)
    return np.linalg.inv(transform) @ P @ np.linalg.inv(transform).T


def normalized_rank_one(raw_coefficients: np.ndarray) -> np.ndarray:
    """L2-normalized rank-one Gram from raw monomial cubic coefficients."""
    norm = raw_coefficients @ M @ raw_coefficients
    # f=c^T m=(M^(1/2)c)^T e.
    e_coefficients = M_HALF @ raw_coefficients
    return np.outer(e_coefficients, e_coefficients) / norm


def kkt_data(S: np.ndarray):
    N = np.trace(S)
    C = cp_operator(S)
    D = np.sum(S * C)
    P = N * C - D * np.eye(10)
    St = to_tensor_gram(S)
    Pt = functional_to_tensor(P)
    R2P = reduce32(Pt)
    R1P = reduce21(R2P)
    q3 = np.sum(St * lift32(R2P))
    R2S = reduce32(St)
    q2 = np.sum(R2S * lift21(R1P))
    # Residual pairings are equal because <S,P>=0 at exact complementarity.
    residual3 = np.sum(St * (lift32(R2P) - Pt))
    residual2 = np.sum(R2S * (lift21(R1P) - R2P))
    Sfull = U3 @ St @ U3.T
    Pfull = U3 @ Pt @ U3.T
    shadows = []
    current = Pfull
    for subsystem in range(3):
        current = reduction_on_subsystem(current, subsystem)
        shadows.append(np.sum(Sfull * current))
    return {
        "N": N,
        "D": D,
        "C_eigs": np.linalg.eigvalsh(C),
        "P_eigs": np.linalg.eigvalsh(P),
        "SP": np.linalg.norm(S @ P),
        "q3": q3,
        "q2": q2,
        "residual3": residual3,
        "residual2": residual2,
        "shadow_pairings": shadows,
        "traces_tensor": (np.trace(St), np.trace(Pt), np.trace(R2P), np.trace(R1P)),
    }


def reduction_on_subsystem(X: np.ndarray, subsystem: int) -> np.ndarray:
    """Apply R=Tr(.)I-(.) to one qutrit factor of a 3-qutrit operator."""
    tensor = X.reshape(3, 3, 3, 3, 3, 3)
    # axes 0,1,2 are ket indices and 3,4,5 their bra partners.
    reduced = np.trace(tensor, axis1=subsystem, axis2=subsystem + 3)
    lifted = np.zeros_like(tensor)
    other = [axis for axis in range(3) if axis != subsystem]
    for value in range(3):
        index = [slice(None)] * 6
        index[subsystem] = value
        index[subsystem + 3] = value
        # np.trace preserves the other ket axes followed by other bra axes.
        lifted[tuple(index)] = reduced
    return lifted.reshape(27, 27) - X


def equality_states():
    # Raw basis is z^3,yz^2,y^2z,y^3,xz^2,xyz,xy^2,x^2z,x^2y,x^3.
    pole = np.zeros(10)
    pole[C3.index((2, 0, 1))] = 2
    pole[C3.index((0, 2, 1))] = 2
    xyz = np.zeros(10)
    xyz[C3.index((1, 1, 1))] = 4 * np.sqrt(2)
    return normalized_rank_one(pole), normalized_rank_one(xyz)


if __name__ == "__main__":
    for name, state in zip(("pole", "xyz"), equality_states(), strict=True):
        print(name)
        for key, value in kkt_data(state).items():
            print(" ", key, value)
    search_elementary_identity()
