"""Search a finite contact-plus-Hessian certificate for the fusion cap.

For a hypothetical global minimizer the potential gap Q(z)=U(z)-E is
nonnegative and the indexed mass/position Hessian is PSD.  This script asks
whether the cap operator can be written as a positive localized contact Gram
plus an equivariant CP image of that Hessian.  The contact feature basis is

    (z.S.z) (z.x)^k,  (x.S.x) (z.x)^k,  (z.S.x) (z.x)^k.

Sphere integration is deterministic Gauss--Legendre x Fourier quadrature,
exact for the polynomial degrees used (up to floating evaluation error).
"""

from __future__ import annotations

import argparse

import cvxpy as cp
import numpy as np
from numpy.polynomial.legendre import leggauss

from cap_hessian_cp_fit import CHANNELS, contact_matrices, cp_image
from weighted_certificate_fermion_cap import (
    BASIS,
    cap_matrix,
    coupled_mass_position_hessian,
    kernel_prime,
)


def sphere_quadrature(maximum_degree: int):
    radial_count = maximum_degree // 2 + 2
    azimuth_count = 2 * maximum_degree + 3
    z_coordinate, radial_weight = leggauss(radial_count)
    azimuth = 2 * np.pi * np.arange(azimuth_count) / azimuth_count
    points = []
    weights = []
    for value, weight in zip(z_coordinate, radial_weight, strict=True):
        radius = np.sqrt(max(0.0, 1 - value * value))
        for angle in azimuth:
            points.append([radius * np.cos(angle), radius * np.sin(angle), value])
            weights.append(weight / (2 * azimuth_count))
    return np.asarray(points), np.asarray(weights)


def contact_block_coefficients(
    points: np.ndarray,
    weights: np.ndarray,
    maximum_power: int,
):
    maximum_degree = 10 + 2 * maximum_power
    roots, root_weights = sphere_quadrature(maximum_degree)
    inner = roots @ points.T
    kernel = 32 * inner**6 - 48 * inner**4 + 20 * inner**2 - 4 / 3
    atom_gram = points @ points.T
    atom_kernel = 32 * atom_gram**6 - 48 * atom_gram**4 + 20 * atom_gram**2 - 4 / 3
    energy = weights @ atom_kernel @ weights
    gap = kernel @ weights - energy

    q_root = np.array(
        [np.einsum("zi,ij,zj->z", roots, matrix, roots) for matrix in BASIS]
    )
    q_atom = np.array(
        [np.einsum("ni,ij,nj->n", points, matrix, points) for matrix in BASIS]
    )
    mixed = np.array(
        [np.einsum("zi,ij,nj->zn", roots, matrix, points) for matrix in BASIS]
    )

    linear_forms = []
    names = []
    for power in range(maximum_power + 1):
        factor = inner**power
        # Each array has shape (root, atom, spin-coordinate).
        linear_forms.extend(
            [
                factor[:, :, None] * q_root.T[:, None, :],
                factor[:, :, None] * q_atom.T[None, :, :],
                factor[:, :, None] * mixed.transpose(1, 2, 0),
            ]
        )
        names.extend((f"qz_t{power}", f"qx_t{power}", f"b_t{power}"))

    size = len(linear_forms)
    coefficients = [[None for _ in range(size)] for _ in range(size)]
    measure_weights = root_weights[:, None] * gap[:, None] * weights[None, :]
    for left in range(size):
        for right in range(left, size):
            matrix = np.einsum(
                "zn,zna,znb->ab",
                measure_weights,
                linear_forms[left],
                linear_forms[right],
                optimize=True,
            )
            if left != right:
                matrix = (matrix + matrix.T) / 2
            coefficients[left][right] = matrix
            coefficients[right][left] = matrix.T
    return names, coefficients


def support_contact_equalities(
    points: np.ndarray,
    weights: np.ndarray,
    maximum_power: int,
):
    """Flag-matrix consequences of Q(z)=0 on the support.

    This is the free-sign analogue of :func:`contact_block_coefficients`:
    the root is now sampled from the measure itself.  At a KKT measure its
    potential gap vanishes root by root, so every returned matrix is zero.
    Keeping an independent leaf makes this strictly richer than the three
    root-only spin-two contractions in :func:`root_contact_equalities`.
    """
    inner = points @ points.T
    kernel = 32 * inner**6 - 48 * inner**4 + 20 * inner**2 - 4 / 3
    potential = kernel @ weights
    energy = weights @ potential
    gap = potential - energy
    q = np.array(
        [np.einsum("ni,ij,nj->n", points, matrix, points) for matrix in BASIS]
    )
    mixed = np.array(
        [np.einsum("zi,ij,nj->zn", points, matrix, points) for matrix in BASIS]
    )
    linear_forms = []
    for power in range(maximum_power + 1):
        factor = inner**power
        linear_forms.extend(
            [
                factor[:, :, None] * q.T[:, None, :],
                factor[:, :, None] * q.T[None, :, :],
                factor[:, :, None] * mixed.transpose(1, 2, 0),
            ]
        )
    measure_weights = weights[:, None] * gap[:, None] * weights[None, :]
    coefficients = []
    for left in range(len(linear_forms)):
        for right in range(left, len(linear_forms)):
            matrix = np.einsum(
                "zn,zna,znb->ab",
                measure_weights,
                linear_forms[left],
                linear_forms[right],
                optimize=True,
            )
            if left != right:
                matrix = matrix + matrix.T
            else:
                matrix = (matrix + matrix.T) / 2
            coefficients.append(matrix)
    return coefficients


def random_measure(rng: np.random.Generator, atoms: int):
    points = rng.normal(size=(atoms, 3))
    points /= np.linalg.norm(points, axis=1)[:, None]
    weights = rng.dirichlet(np.ones(atoms))
    return points, weights


def root_contact_equalities(points: np.ndarray, weights: np.ndarray):
    gram = points @ points.T
    kernel = 32 * gram**6 - 48 * gram**4 + 20 * gram**2 - 4 / 3
    potential = kernel @ weights
    energy = weights @ potential
    gap = potential - energy
    q = np.array(
        [np.einsum("ni,ij,nj->n", points, matrix, points) for matrix in BASIS]
    )
    scalar = np.einsum("an,bn->nab", q, q)
    tangent = np.zeros((len(points), 5, 5))
    for index, point in enumerate(points):
        fields = [matrix @ point - q[basis_index, index] * point
                  for basis_index, matrix in enumerate(BASIS)]
        tangent[index] = 2 * np.array(
            [[left @ right for right in fields] for left in fields]
        )
    remainder = np.eye(5)[None, :, :] - 1.5 * scalar - tangent
    return [
        np.einsum("n,nab->ab", weights * gap, family)
        for family in (scalar, tangent, remainder)
    ]


def gradient_stationarity_matrix(points: np.ndarray, weights: np.ndarray):
    """The indexed support-gradient equality D+D^T=0.

    At a stationary measure the tangent gradient

        g(x) = int K'(x.y) P_{x perp} y dmu(y)

    vanishes at every support point.  Contract it with the quadratic density
    q_S(x)=x.S.x and projective field V_T(x)=P_{x perp}Tx.  Thus the full
    matrix D_{S,T}=int q_S(x) g(x).V_T(x) dmu(x) is zero, not merely its
    trace.  Its symmetric part is the portion relevant to a symmetric cap
    operator; equivariant CP images below retain all five spin channels.
    """
    gram = points @ points.T
    projected = points[None, :, :] - gram[:, :, None] * points[:, None, :]
    gradient = np.einsum(
        "ij,j,ij,ijk->ik",
        np.ones_like(gram),
        weights,
        kernel_prime(gram),
        projected,
        optimize=True,
    )
    q = np.array(
        [np.einsum("ni,ij,nj->n", points, matrix, points) for matrix in BASIS]
    )
    fields = np.array(
        [points @ matrix - q[index, :, None] * points
         for index, matrix in enumerate(BASIS)]
    )
    contractions = np.einsum("ni,bni->bn", gradient, fields, optimize=True)
    matrix = np.einsum("n,an,bn->ab", weights, q, contractions, optimize=True)
    return matrix + matrix.T


def localized_gradient_equalities(
    points: np.ndarray,
    weights: np.ndarray,
    maximum_power: int,
):
    """Low-degree two-root covariants generated by support gradients.

    For a support root ``z`` and an independent leaf ``x``, every scalar
    below is linear in ``S``.  The right families contain one contraction
    with the tangent gradient ``g(z)`` and hence vanish pointwise at a KKT
    measure.  Pairing them against the ordinary left families retains the
    indexed (rather than merely traced) gradient equations.
    """
    inner = points @ points.T
    projected = points[None, :, :] - inner[:, :, None] * points[:, None, :]
    gradient = np.einsum(
        "j,ij,ijk->ik",
        weights,
        kernel_prime(inner),
        projected,
        optimize=True,
    )
    q = np.array(
        [np.einsum("ni,ij,nj->n", points, matrix, points) for matrix in BASIS]
    )
    mixed = np.array(
        [np.einsum("zi,ij,nj->zn", points, matrix, points) for matrix in BASIS]
    )
    gradient_dot_leaf = gradient @ points.T
    gradient_s_root = np.array(
        [np.einsum("ni,ij,nj->n", gradient, matrix, points)
         for matrix in BASIS]
    )
    gradient_s_leaf = np.array(
        [np.einsum("zi,ij,nj->zn", gradient, matrix, points)
         for matrix in BASIS]
    )

    left_forms = []
    right_forms = []
    for power in range(maximum_power + 1):
        factor = inner**power
        left_forms.extend(
            [
                factor[:, :, None] * q.T[:, None, :],
                factor[:, :, None] * q.T[None, :, :],
                factor[:, :, None] * mixed.transpose(1, 2, 0),
            ]
        )
        right_forms.extend(
            [
                factor[:, :, None] * gradient_s_root.T[:, None, :],
                factor[:, :, None] * gradient_s_leaf.transpose(1, 2, 0),
                factor[:, :, None]
                * gradient_dot_leaf[:, :, None]
                * q.T[:, None, :],
                factor[:, :, None]
                * gradient_dot_leaf[:, :, None]
                * q.T[None, :, :],
                factor[:, :, None]
                * gradient_dot_leaf[:, :, None]
                * mixed.transpose(1, 2, 0),
            ]
        )
    measure_weights = weights[:, None] * weights[None, :]
    coefficients = []
    for left in left_forms:
        for right in right_forms:
            matrix = np.einsum(
                "zn,zna,znb->ab", measure_weights, left, right, optimize=True
            )
            coefficients.append((matrix + matrix.T) / 2)
    return coefficients


def solve(maximum_power: int, sample_count: int, atoms: int, seed: int):
    rng = np.random.default_rng(seed)
    contact_variable = cp.Variable((3 * (maximum_power + 1),) * 2, symmetric=True)
    hessian_variables = [cp.Variable((2, 2), symmetric=True) for _ in range(5)]
    equality_variable = cp.Variable(3)
    gradient_variables = cp.Variable(5)
    support_variable = cp.Variable(
        (3 * (maximum_power + 1)) * (3 * (maximum_power + 1) + 1) // 2
    )
    localized_gradient_variable = cp.Variable(
        15 * (maximum_power + 1) ** 2
    )
    scalar_variable = cp.Variable(2)
    root_contact_variable = cp.Variable(3, nonneg=True)
    constraints = [contact_variable >> 0]
    constraints.extend(variable >> 0 for variable in hessian_variables)
    # A scalar remainder a+bE is automatically nonnegative on a
    # hypothetical negative-energy branch when a>=0 and b<=0.
    constraints.extend((scalar_variable[0] >= 0, scalar_variable[1] <= 0))
    residuals = []
    numerical_design = []
    numerical_target = []
    names = None
    upper = np.triu_indices(5)
    for _ in range(sample_count):
        points, weights = random_measure(rng, atoms)
        names, contact_coefficients = contact_block_coefficients(
            points, weights, maximum_power
        )
        root_contacts, _ = contact_matrices(points, weights)
        _, mass, cross, position = coupled_mass_position_hessian(points, weights)
        gram = points @ points.T
        atom_kernel = 32 * gram**6 - 48 * gram**4 + 20 * gram**2 - 4 / 3
        energy = weights @ atom_kernel @ weights
        equalities = root_contact_equalities(points, weights)
        gradient_equality = gradient_stationarity_matrix(points, weights)
        support_equalities = support_contact_equalities(
            points, weights, maximum_power
        )
        localized_gradient = localized_gradient_equalities(
            points, weights, maximum_power
        )
        represented = 0
        columns = []
        for left in range(len(names)):
            for right in range(left, len(names)):
                coefficient = contact_coefficients[left][right]
                scale = 1 if left == right else 2
                represented += scale * contact_variable[left, right] * coefficient
                columns.append((scale * coefficient)[upper])
        for index, matrix in enumerate(root_contacts):
            represented += root_contact_variable[index] * matrix
            columns.append(matrix[upper])
        for channel, variable in zip(CHANNELS, hessian_variables, strict=True):
            channel_matrices = [
                cp_image(channel, mass),
                cp_image(channel, cross + cross.T),
                cp_image(channel, position),
            ]
            represented += variable[0, 0] * channel_matrices[0]
            represented += variable[0, 1] * channel_matrices[1]
            represented += variable[1, 1] * channel_matrices[2]
            columns.extend(matrix[upper] for matrix in channel_matrices)
        for index, matrix in enumerate(equalities):
            represented += equality_variable[index] * matrix
            columns.append(matrix[upper])
        for index, channel in enumerate(CHANNELS):
            matrix = cp_image(channel, gradient_equality)
            represented += gradient_variables[index] * matrix
            columns.append(matrix[upper])
        for index, matrix in enumerate(support_equalities):
            represented += support_variable[index] * matrix
            columns.append(matrix[upper])
        for index, matrix in enumerate(localized_gradient):
            represented += localized_gradient_variable[index] * matrix
            columns.append(matrix[upper])
        for index, scalar_matrix in enumerate((np.eye(5), energy * np.eye(5))):
            represented += scalar_variable[index] * scalar_matrix
            columns.append(scalar_matrix[upper])
        target = cap_matrix(points, weights) / 2
        residuals.append(cp.vec(represented - target, order="C"))
        numerical_design.append(np.stack(columns, axis=1))
        numerical_target.append(target[upper])

    design = np.concatenate(numerical_design)
    target_vector = np.concatenate(numerical_target)
    _, _, rank, singular = np.linalg.lstsq(design, target_vector, rcond=None)
    augmented_rank = np.linalg.matrix_rank(
        np.c_[design, target_vector], tol=1e-8
    )
    print(
        "features", names,
        "design", design.shape,
        "rank", rank,
        "augmented", augmented_rank,
        "singular tail", singular[-3:],
    )

    residual = cp.hstack(residuals)
    problem = cp.Problem(cp.Minimize(cp.sum_squares(residual)), constraints)
    value = problem.solve(
        solver="MOSEK",
        mosek_params={"MSK_DPAR_INTPNT_CO_TOL_REL_GAP": 1e-10},
    )
    print("status", problem.status, "squared residual", value)
    if contact_variable.value is not None:
        print("contact eig", np.linalg.eigvalsh(contact_variable.value))
        for ell, variable in enumerate(hessian_variables):
            print("hessian", ell, variable.value, np.linalg.eigvalsh(variable.value))
        print("root equalities", equality_variable.value)
        print("gradient equalities", gradient_variables.value)
        print("support equalities", support_variable.value)
        print("localized gradient", localized_gradient_variable.value)
        print("scalar I, E I", scalar_variable.value)
        print("root contact", root_contact_variable.value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-power", type=int, default=0)
    parser.add_argument("--samples", type=int, default=20)
    parser.add_argument("--atoms", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260820)
    args = parser.parse_args()
    solve(args.maximum_power, args.samples, args.atoms, args.seed)


if __name__ == "__main__":
    main()
