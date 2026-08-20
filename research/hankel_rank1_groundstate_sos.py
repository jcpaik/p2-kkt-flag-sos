"""Search an SOS certificate for the rank-one KKT ground-state condition.

For a cubic ``f`` (in the spherical-L2 orthonormal basis), put

    N = ||f||^2,
    D = <f^2,T f^2>,
    C_f(g,h) = <gh,T f^2>,
    P_f = N C_f - D I.

At a rank-one KKT point, positivity of the reconstructed Hankel matrix is
equivalent to ``P_f >= 0``.  The sharp desired conclusion is

    -N (N^2 + c0 D) >= 0,

where c0 is the constant harmonic coefficient of the kernel.  This file
searches the first nontrivial invariant Positivstellensatz

    -N (N^2+c0 D)
       = sum_a (B_a f)^T P_f (B_a f) + SOS_6(f).

The covariance of the linear maps B_a is an invariant PSD operator on
End(Sym^3 R^3).  The sextic SOS Gram is block-diagonalized exactly by spin;
its largest PSD multiplicity block is only 7 by 7 (rather than 220 by 220).
Coefficient equality is imposed by generic samples and then independently
validated.  This is a certificate search, not by itself an exact proof.
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


MONOMIALS3 = list(combinations_with_replacement(range(10), 3))
COUNTS3 = np.array([np.bincount(monomial, minlength=10) for monomial in MONOMIALS3])


def symmetric_cube_feature(z: np.ndarray) -> np.ndarray:
    """ON coordinates of z^(symmetric tensor 3); norm squared is ||z||^6."""
    out = []
    for monomial, counts in zip(MONOMIALS3, COUNTS3):
        multiplicity = math.factorial(3)
        for count in counts:
            multiplicity //= math.factorial(int(count))
        out.append(math.sqrt(multiplicity) * np.prod(z[list(monomial)]))
    return np.asarray(out)


def symmetric_cube_generator(generator: np.ndarray) -> np.ndarray:
    """Derived action on the ON symmetric-cube (bosonic Fock) basis."""
    index = {monomial: i for i, monomial in enumerate(MONOMIALS3)}
    out = np.zeros((len(MONOMIALS3), len(MONOMIALS3)))
    for column, counts in enumerate(COUNTS3):
        for source in range(10):
            if counts[source] == 0:
                continue
            for target in range(10):
                value = generator[target, source]
                if abs(value) < 1e-14:
                    continue
                changed = counts.copy()
                changed[source] -= 1
                changed[target] += 1
                monomial = tuple(
                    coordinate
                    for coordinate, count in enumerate(changed)
                    for _ in range(int(count))
                )
                row = index[monomial]
                if source == target:
                    factor = counts[source]
                else:
                    factor = math.sqrt(counts[source] * (counts[target] + 1))
                out[row, column] += value * factor
    return out


GENERATORS = [symmetric_cube_generator(generator) for generator in cov.GENERATORS_V]


def spin_frames() -> dict[int, list[np.ndarray]]:
    """Aligned weight frames U_m for every spin isotypic component.

    For spin ell with multiplicity r, each returned U_m is 220 by r and the
    invariant Gram associated with A>=0 is sum_m U_m A U_m^*.
    """
    casimir = -sum(generator @ generator for generator in GENERATORS)
    values, vectors = np.linalg.eigh(casimir)

    # The raw rotation ordering is xy, xz, yz.  Thus Jz=Jxy, Jx=Jyz,
    # and Jy=-Jxz gives the standard oriented commutation relations.
    jz, jx, jy = GENERATORS[0], GENERATORS[2], -GENERATORS[1]
    result = {}
    for ell in range(10):
        E = vectors[:, np.abs(values - ell * (ell + 1)) < 2e-7]
        if E.shape[1] == 0:
            continue
        multiplicity = E.shape[1] // (2 * ell + 1)
        hz = 1j * (E.T @ jz @ E)
        weights, weight_vectors = np.linalg.eigh(hz)
        start = weight_vectors[:, np.abs(weights + ell) < 2e-7]
        if start.shape[1] != multiplicity:
            raise RuntimeError((ell, E.shape, weights, start.shape))

        plus_candidates = [
            E.T @ (jx + 1j * jy) @ E,
            E.T @ (jx - 1j * jy) @ E,
        ]
        errors = [np.linalg.norm(hz @ L - L @ hz - L) for L in plus_candidates]
        raising = plus_candidates[int(np.argmin(errors))]
        if min(errors) > 1e-6:
            raise RuntimeError((ell, errors))

        local = start
        frames = [E @ local]
        for weight in range(-ell, ell):
            factor = math.sqrt((ell - weight) * (ell + weight + 1))
            local = raising @ local / factor
            frames.append(E @ local)
        combined = np.hstack(frames)
        if np.linalg.norm(combined.conj().T @ combined - np.eye(combined.shape[1])) > 2e-6:
            raise RuntimeError((ell, "nonorthogonal frames"))
        result[ell] = frames
        print("spin", ell, "multiplicity", multiplicity)
    return result


SPIN_FRAMES = spin_frames()


def sos_rows(z: np.ndarray) -> dict[int, np.ndarray]:
    feature = symmetric_cube_feature(z)
    out = {}
    for ell, frames in SPIN_FRAMES.items():
        matrix = sum(
            np.outer(frame.conj().T @ feature, frame.T @ feature)
            for frame in frames
        )
        out[ell] = np.real((matrix + matrix.conj().T) / 2)
    return out


def kernel_data(name: str):
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


def solve(name="E", samples=900, seed=7781):
    D_matrix, c0 = kernel_data(name)
    data.D_MATRIX = D_matrix
    rng = np.random.default_rng(seed + (name == "J"))

    covariance_rows = []
    targets = []
    block_rows = {ell: [] for ell in SPIN_FRAMES}
    saved = []
    for _ in range(samples):
        z = rng.normal(size=10)
        raw = cov.INV_SQRT_G @ z
        square = data.square_coefficients(raw)
        N = z @ z
        D = square @ D_matrix @ square
        row0, row1, _, _ = cov.sample_rows(z)
        covariance_rows.append(N * row1 - D * row0)
        targets.append(-N * (N * N + c0 * D))
        rows = sos_rows(z)
        for ell in SPIN_FRAMES:
            block_rows[ell].append(rows[ell])
        saved.append((z, N, D))

    covariance_rows = np.asarray(covariance_rows)
    targets = np.asarray(targets)
    block_rows = {ell: np.asarray(rows) for ell, rows in block_rows.items()}

    covariance_coefficients = cp.Variable(len(cov.BASIS))
    covariance = sum(
        covariance_coefficients[i] * cov.BASIS[i] for i in range(len(cov.BASIS))
    )
    blocks = {
        ell: cp.Variable((frames[0].shape[1], frames[0].shape[1]), symmetric=True)
        for ell, frames in SPIN_FRAMES.items()
    }
    # Only 17 invariant sextic coefficient equations are independent.  The
    # raw sample equations are extremely redundant and make conic solvers
    # incorrectly report failure, so compress them by an SVD first.
    scalar_variables = [covariance_coefficients]
    design_parts = [covariance_rows]
    for ell, block in blocks.items():
        local_variables = []
        local_rows = []
        for i in range(block.shape[0]):
            for j in range(i, block.shape[1]):
                local_variables.append(block[i, j])
                local_rows.append(block_rows[ell][:, i, j] * (2 if i != j else 1))
        scalar_variables.append(cp.hstack(local_variables))
        design_parts.append(np.stack(local_rows, axis=1))
    scalar_variables = cp.hstack(scalar_variables)
    design = np.hstack(design_parts)
    U, singular, Vh = np.linalg.svd(design, full_matrices=False)
    rank = int(np.sum(singular > 1e-8 * singular[0]))
    compressed_matrix = Vh[:rank]
    compressed_target = (U[:, :rank].T @ targets) / singular[:rank]
    print(name, "coefficient rank", rank, "span residual", np.linalg.norm(
        targets - U[:, :rank] @ (U[:, :rank].T @ targets)
    ))
    constraints = [covariance >> 0]
    constraints += [block >> 0 for block in blocks.values()]
    constraints.append(compressed_matrix @ scalar_variables == compressed_target)
    objective = cp.trace(covariance) + sum(
        (2 * ell + 1) * cp.trace(block) for ell, block in blocks.items()
    )
    problem = cp.Problem(cp.Minimize(objective), constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                value = problem.solve(
                    solver=solver,
                    tol_gap_abs=2e-9,
                    tol_feas=2e-9,
                    tol_gap_rel=2e-9,
                    max_iter=4000,
                )
            else:
                value = problem.solve(solver=solver, eps=5e-7, max_iters=500000)
            print(name, solver, problem.status, value)
            if covariance_coefficients.value is not None:
                Z = sum(
                    covariance_coefficients.value[i] * cov.BASIS[i]
                    for i in range(len(cov.BASIS))
                )
                print("covariance eig", np.linalg.eigvalsh(Z)[:20])
                for ell, block in blocks.items():
                    print("SOS", ell, np.linalg.eigvalsh(block.value))

                errors = []
                for _ in range(200):
                    z = rng.normal(size=10)
                    raw = cov.INV_SQRT_G @ z
                    square = data.square_coefficients(raw)
                    N = z @ z
                    D = square @ D_matrix @ square
                    row0, row1, _, _ = cov.sample_rows(z)
                    actual = (N * row1 - D * row0) @ covariance_coefficients.value
                    rows = sos_rows(z)
                    actual += sum(
                        np.sum(rows[ell] * block.value) for ell, block in blocks.items()
                    )
                    errors.append(actual + N * (N * N + c0 * D))
                print("validation", np.max(np.abs(errors)))
            if problem.status in ("optimal", "optimal_inaccurate"):
                break
        except cp.error.SolverError as error:
            print(name, solver, error)


if __name__ == "__main__":
    solve("E")
    solve("J")
