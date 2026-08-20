"""Degree-eight Positivstellensatz search for the rank-one slack branch.

For a cubic ``f`` in the spherical-L2 orthonormal coordinates, let

    P_f = N C_f - D I,   N = ||f||^2,   D = <f,C_f f>.

The reconstructed Hankel matrix is positive semidefinite exactly when
``P_f >= 0`` on the negative-D branch.  Since ``f.T P_f f`` vanishes
identically, positivity already implies the rank-one KKT equation.

The degree-six certificate with linear vector multipliers is infeasible.  This
script tests the next complete invariant level

    -N^2 (N^2 + c0 D)
      = sum_a z_a(f)^T P_f z_a(f) + SOS_8(f),

where the ``z_a`` range over arbitrary quadratic vector-valued forms.  Their
Gram covariance is block diagonalized under SO(3) on
``Hom(Sym^2(V),V)``.  The octic SOS Gram is likewise block diagonalized on
``Sym^4(V)``.  A feasible identity would prove the entire rank-one slack
branch (not merely a numerical lower bound).
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from itertools import combinations_with_replacement

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np

import hankel_rank1_covariant_sdp as cov
import hankel_rank1_slack_invariant as data


def symmetric_power_data(dimension: int, degree: int, generators):
    monomials = list(combinations_with_replacement(range(dimension), degree))
    counts = np.array([np.bincount(m, minlength=dimension) for m in monomials])
    index = {tuple(row): i for i, row in enumerate(counts)}

    def feature(z):
        out = np.empty(len(monomials))
        for row, count in enumerate(counts):
            multiplicity = math.factorial(degree)
            for value in count:
                multiplicity //= math.factorial(int(value))
            out[row] = math.sqrt(multiplicity) * np.prod(z ** count)
        return out

    lifted = []
    for generator in generators:
        out = np.zeros((len(monomials), len(monomials)))
        for column, count in enumerate(counts):
            for source in range(dimension):
                if count[source] == 0:
                    continue
                for target in range(dimension):
                    value = generator[target, source]
                    if abs(value) < 1e-14:
                        continue
                    changed = count.copy()
                    changed[source] -= 1
                    changed[target] += 1
                    row = index[tuple(changed)]
                    if source == target:
                        factor = count[source]
                    else:
                        factor = math.sqrt(count[source] * (count[target] + 1))
                    out[row, column] += value * factor
        lifted.append(out)
    return feature, lifted


def spin_frames(generators, max_spin=15):
    """Return aligned weight frames for every SO(3) isotypic component."""
    casimir = -sum(generator @ generator for generator in generators)
    values, vectors = np.linalg.eigh(casimir)
    jz, jx, jy = generators[0], generators[2], -generators[1]
    result = {}
    for ell in range(max_spin + 1):
        E = vectors[:, np.abs(values - ell * (ell + 1)) < 3e-6]
        if E.shape[1] == 0:
            continue
        multiplicity = E.shape[1] // (2 * ell + 1)
        hz = 1j * (E.T @ jz @ E)
        weights, weight_vectors = np.linalg.eigh(hz)
        local = weight_vectors[:, np.abs(weights + ell) < 3e-6]
        if local.shape[1] != multiplicity:
            raise RuntimeError((ell, E.shape, local.shape, weights))
        candidates = [E.T @ (jx + 1j * jy) @ E, E.T @ (jx - 1j * jy) @ E]
        errors = [np.linalg.norm(hz @ L - L @ hz - L) for L in candidates]
        raising = candidates[int(np.argmin(errors))]
        frames = [E @ local]
        for weight in range(-ell, ell):
            local = raising @ local / math.sqrt((ell - weight) * (ell + weight + 1))
            frames.append(E @ local)
        combined = np.hstack(frames)
        error = np.linalg.norm(combined.conj().T @ combined - np.eye(combined.shape[1]))
        if error > 2e-5:
            raise RuntimeError((ell, "frame orthogonality", error))
        result[ell] = frames
        print("spin", ell, "dimension", E.shape[1], "multiplicity", multiplicity)
    return result


SQUARE_FEATURE, SQUARE_GENERATORS = symmetric_power_data(10, 2, cov.GENERATORS_V)
QUARTIC_FEATURE, QUARTIC_GENERATORS = symmetric_power_data(10, 4, cov.GENERATORS_V)

# Row-major vectorization of a 10 by 55 map L: Sym^2(V) -> V.
MAP_GENERATORS = [
    np.kron(generator, np.eye(55)) - np.kron(np.eye(10), lifted.T)
    for generator, lifted in zip(cov.GENERATORS_V, SQUARE_GENERATORS)
]


def kernel_data(name):
    if name == "E":
        coefficients = data.C
        c0 = Fraction(32, 105)
    elif name == "J":
        coefficients = {
            2: Fraction(22, 35),
            4: -Fraction(192, 385),
            6: Fraction(768, 1001),
        }
        c0 = Fraction(48, 35)
    else:
        raise ValueError(name)
    exact = [
        [
            sum(data.B_EXACT[ell][i][j] / coefficients[ell] for ell in (2, 4, 6))
            for j in range(28)
        ]
        for i in range(28)
    ]
    return np.asarray(exact, dtype=float), float(c0)


def c_matrix(z, d_matrix):
    raw = cov.INV_SQRT_G @ z
    square = data.square_coefficients(raw)
    raw_matrix = np.array(
        [[d_matrix[data.PRODUCT_INDEX[i, j]] @ square for j in range(10)] for i in range(10)]
    )
    return cov.INV_SQRT_G.T @ raw_matrix @ cov.INV_SQRT_G, square


def multiplier_row(frames, q, P):
    """Multiplicity-space coefficient of an invariant quadratic-map Gram."""
    multiplicity = frames[0].shape[1]
    out = np.zeros((multiplicity, multiplicity), dtype=np.complex128)
    for frame in frames:
        # Each column is a row-major 10 by 55 map.
        values = np.stack([frame[:, j].reshape(10, 55) @ q for j in range(multiplicity)], axis=1)
        out += values.conj().T @ P @ values
    return np.real((out + out.conj().T) / 2)


def sos_row(frames, feature):
    multiplicity = frames[0].shape[1]
    out = np.zeros((multiplicity, multiplicity), dtype=np.complex128)
    for frame in frames:
        value = frame.conj().T @ feature
        out += np.outer(value.conj(), value)
    return np.real((out + out.conj().T) / 2)


def solve(name="J", samples=500, seed=20260820):
    d_matrix, c0 = kernel_data(name)
    print("decomposing quadratic vector multipliers")
    multiplier_frames = spin_frames(MAP_GENERATORS, 9)
    print("decomposing octic SOS features")
    sos_frames = spin_frames(QUARTIC_GENERATORS, 12)

    rng = np.random.default_rng(seed + (name == "J"))
    multiplier_rows = {ell: [] for ell in multiplier_frames}
    sos_rows = {ell: [] for ell in sos_frames}
    targets = []
    for sample in range(samples):
        z = rng.normal(size=10)
        z /= np.linalg.norm(z)
        C, square = c_matrix(z, d_matrix)
        D = square @ d_matrix @ square
        N = z @ z
        P = N * C - D * np.eye(10)
        q = SQUARE_FEATURE(z)
        quartic = QUARTIC_FEATURE(z)
        for ell, frames in multiplier_frames.items():
            multiplier_rows[ell].append(multiplier_row(frames, q, P))
        for ell, frames in sos_frames.items():
            sos_rows[ell].append(sos_row(frames, quartic))
        targets.append(-N * N * (N * N + c0 * D))
        if sample and sample % 50 == 0:
            print("sample", sample)

    multiplier_rows = {ell: np.asarray(rows) for ell, rows in multiplier_rows.items()}
    sos_rows = {ell: np.asarray(rows) for ell, rows in sos_rows.items()}
    targets = np.asarray(targets)

    blocks_m = {
        ell: cp.Variable((frames[0].shape[1], frames[0].shape[1]), symmetric=True)
        for ell, frames in multiplier_frames.items()
    }
    blocks_s = {
        ell: cp.Variable((frames[0].shape[1], frames[0].shape[1]), symmetric=True)
        for ell, frames in sos_frames.items()
    }

    scalar_variables = []
    design_parts = []
    for blocks, rows_by_spin in ((blocks_m, multiplier_rows), (blocks_s, sos_rows)):
        for ell, block in blocks.items():
            variables = []
            columns = []
            rows = rows_by_spin[ell]
            for i in range(block.shape[0]):
                for j in range(i, block.shape[1]):
                    variables.append(block[i, j])
                    columns.append(rows[:, i, j] * (2 if i != j else 1))
            scalar_variables.extend(variables)
            design_parts.extend(columns)
    variables = cp.hstack(scalar_variables)
    design = np.stack(design_parts, axis=1)
    U, singular, Vh = np.linalg.svd(design, full_matrices=False)
    rank = int(np.sum(singular > 2e-9 * singular[0]))
    residual = np.linalg.norm(targets - U[:, :rank] @ (U[:, :rank].T @ targets))
    print(name, "design", design.shape, "rank", rank, "target residual", residual)
    equations = Vh[:rank] @ variables == (U[:, :rank].T @ targets) / singular[:rank]
    constraints = [equations]
    constraints += [block >> 0 for block in blocks_m.values()]
    constraints += [block >> 0 for block in blocks_s.values()]
    objective = cp.Minimize(
        sum((2 * ell + 1) * cp.trace(block) for ell, block in blocks_m.items())
        + sum((2 * ell + 1) * cp.trace(block) for ell, block in blocks_s.items())
    )
    problem = cp.Problem(objective, constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                value = problem.solve(
                    solver=solver,
                    tol_gap_abs=2e-8,
                    tol_feas=2e-8,
                    tol_gap_rel=2e-8,
                    max_iter=4000,
                )
            else:
                value = problem.solve(solver=solver, eps=2e-6, max_iters=500000)
            print(name, solver, problem.status, value)
            if any(block.value is not None for block in blocks_m.values()):
                for family, blocks in (("M", blocks_m), ("S", blocks_s)):
                    for ell, block in blocks.items():
                        print(family, ell, np.linalg.eigvalsh(block.value)[:8])
            if problem.status in ("optimal", "optimal_inaccurate"):
                break
        except cp.error.SolverError as error:
            print(name, solver, error)


if __name__ == "__main__":
    solve(sys.argv[1] if len(sys.argv) > 1 else "J")
