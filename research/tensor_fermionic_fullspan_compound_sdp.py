"""Minor/SOS audits after restoring the omitted Pluecker relations.

The minimal 40-dimensional relaxation uses trace ratio, one spin-2 relation,
and the spin-4 matching relation.  A genuine tangent moment also obeys

* the first Bianchi identity (five further spin-2 equations), and
* absence of the spin-3 cross block (seven equations).

Together these cut the symmetric 10-by-10 space to the true 28-dimensional
linear span H0+H2+H4+H6.  This script tests whether ordinary squares plus
2-by-2 minors certify the target after either or both missing relations are
restored.
"""

import itertools
import sys
import os

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

import tensor_fermionic_general_relaxation_opt as rel
import tensor_fermionic_l4_compound_sdp as old


def bianchi_rows():
    pair_index = {pair: a for a, pair in enumerate(rel.wedge.pairs)}
    rows = []
    for i, j, k, ell in itertools.combinations(range(5), 4):
        row = np.zeros((10, 10))
        row[pair_index[i, j], pair_index[k, ell]] += 1
        row[pair_index[i, k], pair_index[j, ell]] -= 1
        row[pair_index[i, ell], pair_index[j, k]] += 1
        rows.append(row.reshape(-1) @ rel.TBASIS)
    return np.stack(rows)


BIANCHI = bianchi_rows()


def cross_spin3_rows():
    rows = []
    for constraint in rel.CROSS_NULL:
        row = np.zeros((10, 10))
        row[:3, 3:] = constraint.reshape(3, 7)
        rows.append(row.reshape(-1))
    return np.stack(rows)


CROSS3 = cross_spin3_rows()


def build_basis(use_bianchi=True, use_cross3=True):
    trace_ratio = np.zeros((1, 100))
    for i in range(3):
        trace_ratio[0, 10 * i + i] = 5
    for i in range(10):
        trace_ratio[0, 10 * i + i] -= 1
    blocks = [trace_ratio, rel.LINEAR_CONSTRAINT_REDUCED, rel.L4_REDUCED]
    if use_bianchi:
        blocks.append(BIANCHI)
    if use_cross3:
        blocks.append(CROSS3)
    coordinates = np.vstack(blocks) @ old.EXPANSION
    null = null_space(coordinates)
    basis = [
        sum(null[a, k] * old.SYM[a] for a in range(len(old.SYM)))
        for k in range(null.shape[1])
    ]
    return coordinates, basis


def target_matrix(basis):
    n = len(basis)
    out = np.zeros((n, n))
    weighted = os.environ.get("FERMIONIC_WEIGHTED", "0") == "1"
    for a, A in enumerate(basis):
        FA = old.contraction(A)
        Adev = A[:3, :3] - np.eye(3) * np.trace(A) / 15
        for b, B in enumerate(basis):
            FB = old.contraction(B)
            Bdev = B[:3, :3] - np.eye(3) * np.trace(B) / 15
            out[a, b] = (
                np.sum(A * B)
                - 0.5 * np.sum(FA * FB)
                + np.trace(A) * np.trace(B) / 3
                - (25 / 12) * np.sum(Adev * Bdev) * weighted
            )
    return (out + out.T) / 2


def solve(use_bianchi=True, use_cross3=True):
    coordinates, basis = build_basis(use_bianchi, use_cross3)
    target = target_matrix(basis)
    n = len(basis)
    compounds = {}
    for a in range(n):
        for b in range(a, n):
            compounds[a, b] = old.polarized_compound(basis[a], basis[b])
    W = cp.Variable((45, 45), symmetric=True)
    S = cp.Variable((n, n), symmetric=True)
    constraints = [W >> 0, S >> 0]
    for a in range(n):
        for b in range(a, n):
            constraints.append(
                S[a, b] + cp.sum(cp.multiply(W, compounds[a, b]))
                == target[a, b]
            )
    problem = cp.Problem(cp.Minimize(cp.trace(W) + cp.trace(S)), constraints)
    try:
        value = problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-9,
            tol_feas=2e-9,
            tol_gap_rel=2e-9,
            max_iter=3000,
        )
    except cp.error.SolverError:
        value = problem.solve(solver="SCS", eps=3e-7, max_iters=500000)
    print(
        "variant", use_bianchi, use_cross3,
        "dimension/rank", n, np.linalg.matrix_rank(coordinates),
        "target eig", np.linalg.eigvalsh(target),
        "status/value", problem.status, value,
    )
    if W.value is not None:
        print("W eig", np.linalg.eigvalsh(W.value))
        print("S eig", np.linalg.eigvalsh(S.value))


if __name__ == "__main__":
    solve(use_bianchi=True, use_cross3=False)
    solve(use_bianchi=True, use_cross3=True)
