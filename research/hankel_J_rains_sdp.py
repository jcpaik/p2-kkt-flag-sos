"""Decomposable/PPT-shadow search for the strengthened two-copy witness.

On Sym^3(R^3) tensor Sym^3(R^3), the witness

    W = 144 F123 - 72(F12+F13+F23) + 29(F1+F2+F3) - 5 I

has expectation 144 p6-216 p4+87 p2-5 on H tensor H.  A fully Hankel H
is invariant under every local partial transpose.  We therefore search the
strictly larger Rains/decomposable ansatz

    W = sum_{w=0}^3 average_{|S|=w} Phi_S(Z_w),     Z_w >= 0,

where Phi_S embeds in the six qutrits, partially transposes the A particles
in S, and projects back.  By SO(3) averaging, each Z_w may be restricted to
the 30-dimensional symmetric commutant.
"""

import itertools
import math

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space


TRIPLES = list(itertools.product(range(3), repeat=3))
TRIPLE_INDEX = {word: i for i, word in enumerate(TRIPLES)}
COMPS = [(a, b, 3 - a - b) for a in range(4) for b in range(4 - a)]
S = np.zeros((27, 10))
for column, alpha in enumerate(COMPS):
    words = [word for word in TRIPLES if tuple(word.count(i) for i in range(3)) == alpha]
    for word in words:
        S[TRIPLE_INDEX[word], column] = 1 / math.sqrt(len(words))


def triple_generators():
    generators = []
    for i, j in ((0, 1), (0, 2), (1, 2)):
        L = np.zeros((3, 3))
        L[i, j], L[j, i] = 1, -1
        full = np.kron(np.kron(L, np.eye(3)), np.eye(3))
        full += np.kron(np.kron(np.eye(3), L), np.eye(3))
        full += np.kron(np.kron(np.eye(3), np.eye(3)), L)
        generators.append(S.T @ full @ S)
    return generators


GENERATORS_V = triple_generators()
GENERATORS_2 = [np.kron(g, np.eye(10)) + np.kron(np.eye(10), g) for g in GENERATORS_V]


def symmetric_basis(n):
    out = []
    for i in range(n):
        E = np.zeros((n, n))
        E[i, i] = 1
        out.append(E)
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n))
            E[i, j] = E[j, i] = 1 / np.sqrt(2)
            out.append(E)
    return out


def invariant_basis():
    casimir = -sum(g @ g for g in GENERATORS_2)
    values, vectors = np.linalg.eigh(casimir)
    full = []
    labels = []
    for ell in range(7):
        E = vectors[:, np.abs(values - ell * (ell + 1)) < 2e-7]
        restricted = [E.T @ g @ E for g in GENERATORS_2]
        symmetric = symmetric_basis(E.shape[1])
        equations = np.stack(
            [
                np.concatenate([(g @ B - B @ g).reshape(-1) for g in restricted])
                for B in symmetric
            ],
            axis=1,
        )
        kernel = (
            np.eye(len(symmetric))
            if np.linalg.norm(equations) < 1e-7
            else null_space(equations, rcond=1e-8)
        )
        labels.append((ell, E.shape[1], kernel.shape[1]))
        for column in range(kernel.shape[1]):
            B = sum(kernel[row, column] * symmetric[row] for row in range(len(symmetric)))
            full.append(E @ B @ E.T)
    print("commutant", labels, "total", len(full))
    return full


BASIS = invariant_basis()

# Orbit averaging also permits every PSD block to commute with interchange
# of the two Sym^3 copies.  Restricting now makes every member of a
# partial-transpose orbit equivalent to one representative below.
COPY_SWAP = np.zeros((100, 100))
for a in range(10):
    for b in range(10):
        COPY_SWAP[10 * b + a, 10 * a + b] = 1
copy_equations = np.stack([(COPY_SWAP @ B @ COPY_SWAP - B).reshape(-1) for B in BASIS], axis=1)
copy_kernel = null_space(copy_equations, rcond=1e-8)
BASIS = [
    sum(copy_kernel[row, column] * BASIS[row] for row in range(len(BASIS)))
    for column in range(copy_kernel.shape[1])
]
print("copy-symmetric commutant", len(BASIS))

# Each six-qutrit word belongs to exactly one symmetric occupation pair, so
# embedding/projection can be represented by a column label and one weight.
WORDS6 = list(itertools.product(range(3), repeat=6))
WORD_INDEX = {word: i for i, word in enumerate(WORDS6)}
LABEL = np.zeros(729, dtype=int)
WEIGHT = np.zeros(729)
for row, word in enumerate(WORDS6):
    left, right = word[:3], word[3:]
    alpha = tuple(left.count(i) for i in range(3))
    beta = tuple(right.count(i) for i in range(3))
    a, b = COMPS.index(alpha), COMPS.index(beta)
    LABEL[row] = 10 * a + b
    WEIGHT[row] = S[TRIPLE_INDEX[left], a] * S[TRIPLE_INDEX[right], b]


def transpose_map(mask):
    """Sparse-index data for B -> V^T (V B V^T)^T_mask V."""
    out_index = np.empty(729 * 729, dtype=int)
    in_index = np.empty(729 * 729, dtype=int)
    coefficient = np.empty(729 * 729)
    cursor = 0
    for bra_index, bra0 in enumerate(WORDS6):
        for ket_index, ket0 in enumerate(WORDS6):
            bra, ket = list(bra0), list(ket0)
            for particle in mask:
                bra[particle], ket[particle] = ket[particle], bra[particle]
            bt, kt = WORD_INDEX[tuple(bra)], WORD_INDEX[tuple(ket)]
            out_index[cursor] = 100 * LABEL[bt] + LABEL[kt]
            in_index[cursor] = 100 * LABEL[bra_index] + LABEL[ket_index]
            coefficient[cursor] = WEIGHT[bt] * WEIGHT[kt] * WEIGHT[bra_index] * WEIGHT[ket_index]
            cursor += 1
    return out_index, in_index, coefficient


# The 64 subsets have eight orbits under simultaneous particle permutations,
# copy interchange, and complementation (full transpose).  Binary masks are
# shown here as tuples of the local axes to transpose.
ORBIT_REPRESENTATIVES = [0, 1, 3, 7, 9, 10, 11, 14]


def mask_tuple(mask):
    return tuple(i for i in range(6) if (mask >> i) & 1)


MAP_DATA = {mask_tuple(mask): transpose_map(mask_tuple(mask)) for mask in ORBIT_REPRESENTATIVES}


def phi(B, mask):
    out_index, in_index, coefficient = MAP_DATA[mask]
    flat = np.bincount(
        out_index,
        weights=coefficient * B.reshape(-1)[in_index],
        minlength=10000,
    )
    return flat.reshape(100, 100)


IMAGES = [[phi(B, mask_tuple(mask)) for B in BASIS] for mask in ORBIT_REPRESENTATIVES]


def local_swap(particle):
    out = np.zeros((729, 729))
    for ket_index, ket0 in enumerate(WORDS6):
        bra = list(ket0)
        bra[particle], bra[particle + 3] = bra[particle + 3], bra[particle]
        out[WORD_INDEX[tuple(bra)], ket_index] = 1
    return out


F = [local_swap(i) for i in range(3)]
V = np.zeros((729, 100))
for row in range(729):
    V[row, LABEL[row]] = WEIGHT[row]
RAW = 144 * (F[0] @ F[1] @ F[2])
RAW -= 72 * (F[0] @ F[1] + F[0] @ F[2] + F[1] @ F[2])
RAW += 29 * (F[0] + F[1] + F[2])
RAW -= 5 * np.eye(729)
TARGET = V.T @ RAW @ V


def solve():
    variables = [cp.Variable(len(BASIS)) for _ in ORBIT_REPRESENTATIVES]
    blocks = [sum(variable[i] * BASIS[i] for i in range(len(BASIS))) for variable in variables]
    represented = sum(
        sum(variables[w][i] * IMAGES[w][i] for i in range(len(BASIS)))
        for w in range(len(ORBIT_REPRESENTATIVES))
    )
    constraints = [block >> 0 for block in blocks]
    constraints.append(represented == TARGET)
    problem = cp.Problem(cp.Minimize(sum(cp.trace(block) for block in blocks)), constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                result = problem.solve(
                    solver=solver,
                    tol_gap_abs=2e-8,
                    tol_feas=2e-8,
                    tol_gap_rel=2e-8,
                    max_iter=3000,
                )
            else:
                result = problem.solve(solver=solver, eps=2e-6, max_iters=500000)
            print(solver, problem.status, result)
            if all(variable.value is not None for variable in variables):
                residual = sum(
                    sum(variables[w].value[i] * IMAGES[w][i] for i in range(len(BASIS)))
                    for w in range(len(ORBIT_REPRESENTATIVES))
                ) - TARGET
                print("residual", np.linalg.norm(residual), np.max(np.abs(residual)))
                for w, block in enumerate(blocks):
                    value = sum(variables[w].value[i] * BASIS[i] for i in range(len(BASIS)))
                    print("block", w, "trace/min", np.trace(value), np.linalg.eigvalsh(value)[0])
            if problem.status in ("optimal", "optimal_inaccurate"):
                break
        except cp.error.SolverError as error:
            print(solver, error)


if __name__ == "__main__":
    print("target eig", np.linalg.eigvalsh(TARGET)[:20])
    solve()
