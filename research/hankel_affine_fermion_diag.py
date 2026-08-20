"""Search the exact affine tangent-plane relaxation on diagonal 2-fermion states.

Use the concrete orthonormal basis of H_2 consisting of two diagonal and
three off-diagonal real quadrics.  The spin-1 summand in wedge^2 H_2 is the
normalized image of the three infinitesimal rotation generators.  For a
diagonal density G=sum p_ij |Ei wedge Ej><...|, all constraints are linear
in the ten edge weights p_ij.

The actual tangent orbit obeys

    M = 5 G_11,                    tr M = 1,
    Pi_2(F) = (3/7) T_{M-I/3},

where F is the one-particle contraction and

    T_A(S) = AS + SA - (2/3) tr(AS) I.

This script searches that affine slice and prints a candidate that can then
be rationalized and checked symbolically.
"""

import itertools

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize


rt2 = np.sqrt(2.0)
rt6 = np.sqrt(6.0)

E = [
    np.diag([1.0, -1.0, 0.0]) / rt2,
    np.diag([1.0, 1.0, -2.0]) / rt6,
    np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]) / rt2,
    np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]) / rt2,
    np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]) / rt2,
]

L = [
    np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]),
    np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]),
    np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
]

pairs = list(itertools.combinations(range(5), 2))


def generator(axis):
    """Matrix of S -> [L_axis,S] in the orthonormal H2 basis."""
    return np.array(
        [[np.trace(E[i] @ (L[axis] @ E[j] - E[j] @ L[axis])) for j in range(5)] for i in range(5)]
    )


J = [generator(a) for a in range(3)]
W1 = np.array([[J[a][i, j] / np.sqrt(5.0) for a in range(3)] for i, j in pairs])


def T_matrix(A):
    """Matrix of T_A on H2."""
    return np.array(
        [
            [
                np.trace(E[i] @ (A @ E[j] + E[j] @ A - (2.0 / 3.0) * np.trace(A @ E[j]) * np.eye(3)))
                for j in range(5)
            ]
            for i in range(5)
        ]
    )


# Orthogonal projector from symmetric End(H2) onto its spin-2 sector.  The
# five matrices T_Ea are mutually orthogonal with a common norm.
TB = np.stack([T_matrix(A) for A in E])
TGRAM = np.einsum("aij,bij->ab", TB, TB)


def pi2(X):
    coeff = np.linalg.solve(TGRAM, np.einsum("aij,ij->a", TB, X))
    return np.einsum("a,aij->ij", coeff, TB)


def data(p):
    G11 = W1.T @ np.diag(p) @ W1
    M = 5.0 * G11
    degrees = np.zeros(5)
    for value, (i, j) in zip(p, pairs):
        degrees[i] += value
        degrees[j] += value
    F = np.diag(degrees)
    relation = pi2(F) - (3.0 / 7.0) * T_matrix(M - np.eye(3) / 3.0)
    gap = np.sum(p**2) - 0.5 * np.sum(degrees**2) + 1.0 / 3.0
    return M, F, relation, gap


def upper_slack(p):
    """The pointwise-valid affine marginal slack P_M-F."""
    M, F, _, _ = data(p)
    P = np.array(
        [
            [np.trace(M @ (E[i] @ E[j] + E[j] @ E[i])) for j in range(5)]
            for i in range(5)
        ]
    )
    return P - F


def linear_constraints():
    # Trace p=1, trace G11=1/5, and all entries of the spin-2 relation.
    rows = [np.ones(10), np.sum(W1**2, axis=1)]
    rhs = [1.0, 0.2]
    _, _, rel0, _ = data(np.zeros(10))
    for i in range(5):
        for j in range(i, 5):
            row = np.array([data(np.eye(10)[k])[2][i, j] - rel0[i, j] for k in range(10)])
            rows.append(row)
            rhs.append(-rel0[i, j])
    A = np.stack(rows)
    b = np.array(rhs)
    # Select independent rows of [A|b] using an SVD of A^T.
    q, r = np.linalg.qr(A.T)
    del q
    keep = np.where(np.abs(np.diag(r)) > 1e-10)[0]
    return A[keep], b[keep]


def objective(p):
    return data(p)[3]


def main():
    print("W1 orthogonality", W1.T @ W1)
    print("T Gram", TGRAM)
    A, b = linear_constraints()
    print("constraint shape/rank", A.shape, np.linalg.matrix_rank(A), "residual Haar", np.linalg.norm(A @ (np.ones(10) / 10) - b))
    # Get an equality-feasible center and null-space coordinates.
    center = np.linalg.lstsq(A, b, rcond=None)[0]
    N = null_space(A)
    print("center", center, "nullity", N.shape[1], "residual", np.linalg.norm(A @ center - b))

    rng = np.random.default_rng(20260820)
    for impose_upper in (False, True):
        best = None
        for restart in range(100):
            z0 = rng.normal(scale=0.05, size=N.shape[1])
            cons = [{"type": "ineq", "fun": lambda z, k=k: (center + N @ z)[k]} for k in range(10)]
            if impose_upper:
                # For this diagonal ansatz the upper slack is diagonal.
                cons.extend(
                    {
                        "type": "ineq",
                        "fun": lambda z, k=k: upper_slack(center + N @ z)[k, k],
                    }
                    for k in range(5)
                )
            result = minimize(
                lambda z: objective(center + N @ z),
                z0,
                method="SLSQP",
                constraints=cons,
                options={"ftol": 1e-14, "maxiter": 2000},
            )
            p = center + N @ result.x
            record = (objective(p), p, result.success)
            if best is None or record[0] < best[0]:
                best = record
        gap, p, success = best
        M, F, relation, _ = data(p)
        print("upper", impose_upper, "best gap/success", gap, success)
        print("p", repr(p))
        print("M", M)
        print("F", F)
        print(
            "relation residual",
            np.linalg.norm(relation),
            "min p/slack",
            p.min(),
            np.linalg.eigvalsh(upper_slack(p)).min(),
        )


if __name__ == "__main__":
    main()
