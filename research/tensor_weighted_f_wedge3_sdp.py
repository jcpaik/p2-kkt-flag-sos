"""Third-compound audit for the sharp conjecture J >= 108 F.

Let H be the normalized ternary sextic middle catalecticant on Sym^3(R^3).
Both mass*J and F are homogeneous cubics in its 28 Hankel coordinates.  This
script searches the natural exact-form certificate

    mass(H) J(H) - 108 F(H)
      = mass(H) <W2, C2(H)> + <W3, C3(H)>,

with SO(3)-invariant W2,W3 positive semidefinite.  Feasibility would reduce
the proof to positivity of second and third compound matrices of H.
"""

import itertools
import math
import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

import tensor_hankel_fw as moments
import tensor_schur_feasible_scan as rep


PAIRS = np.array(list(itertools.combinations(range(10), 2)))
TRIPLES = np.array(list(itertools.combinations(range(10), 3)))
PAIRS7 = list(itertools.combinations(range(7), 2))


def compound(H, indices):
    submatrices = H[indices[:, None, :, None], indices[None, :, None, :]]
    return np.linalg.det(submatrices)


def polynomial_product(*linear_forms):
    coefficients = {(0, 0, 0): 1.0}
    for form in linear_forms:
        updated = {}
        for alpha, value in coefficients.items():
            for coordinate in range(3):
                beta = list(alpha)
                beta[coordinate] += 1
                beta = tuple(beta)
                updated[beta] = updated.get(beta, 0.0) + value * form[coordinate]
        coefficients = updated
    return coefficients


def cubic_vector(coefficients):
    vector = np.zeros(10)
    for index, alpha in enumerate(rep.comps):
        multinomial = math.factorial(3)
        for entry in alpha:
            multinomial //= math.factorial(entry)
        vector[index] = coefficients.get(alpha, 0.0) / np.sqrt(multinomial)
    return vector


def add_polynomials(*terms):
    out = {}
    for scale, polynomial in terms:
        for alpha, value in polynomial.items():
            out[alpha] = out.get(alpha, 0.0) + scale * value
    return out


def tangent_wedge(x):
    seed = np.eye(3)[np.argmin(np.abs(x))]
    u = np.cross(x, seed)
    u /= np.linalg.norm(u)
    v = np.cross(x, u)
    fc = add_polynomials(
        (1.0, polynomial_product(x, u, u)),
        (-1.0, polynomial_product(x, v, v)),
    )
    fs = polynomial_product(x, u, v)
    ec = np.sqrt(3 / 2) * (rep.U.T @ cubic_vector(fc))
    es = 2 * np.sqrt(3 / 2) * (rep.U.T @ cubic_vector(fs))
    wedge = np.array(
        [ec[i] * es[j] - ec[j] * es[i] for i, j in PAIRS7]
    )
    # Orientation can flip, but the projector below is invariant.
    return wedge / np.linalg.norm(wedge)


def degree_six_evaluation(x):
    return np.array([np.prod(x ** np.array(alpha)) for alpha in moments.deg6])


def fit_w_map(samples=500, seed=919):
    rng = np.random.default_rng(seed)
    domain = []
    target = []
    for _ in range(samples):
        x = rng.normal(size=3)
        x /= np.linalg.norm(x)
        w = tangent_wedge(x)
        domain.append(degree_six_evaluation(x))
        target.append(np.outer(w, w).reshape(-1))
    domain = np.stack(domain)
    target = np.stack(target)
    linear, *_ = np.linalg.lstsq(domain, target, rcond=None)
    residual = np.linalg.norm(domain @ linear - target) / np.linalg.norm(target)
    return linear.T, residual


W_LINEAR, W_RESIDUAL = fit_w_map()


def f_value(y):
    H = np.einsum("abk,k->ab", moments.maps[3], y)
    A = rep.U.T @ H @ rep.U
    C2 = compound(A, np.array(PAIRS7))
    W = (W_LINEAR @ y).reshape(21, 21)
    return (16 / 9) * np.sum(C2 * W)


def j_value(y):
    rho = moments.matrices(y)
    mass = np.trace(rho[3])
    return (
        144 * np.sum(rho[3] ** 2)
        - 216 * np.sum(rho[2] ** 2)
        + 87 * np.sum(rho[1] ** 2)
        - 5 * mass * mass
    )


def direct_f(xs, weights):
    gram = xs @ xs.T
    a = gram[:, :, None, None]
    b = gram[:, None, :, None]
    c = gram[None, :, :, None]
    # Direct loops are clearer for this one-time normalization audit.
    value = 0.0
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            for k, z in enumerate(xs):
                aa, bb, cc = x @ y, x @ z, y @ z
                determinant = 1 + 2 * aa * bb * cc - aa * aa - bb * bb - cc * cc
                value += (
                    8
                    * weights[i]
                    * weights[j]
                    * weights[k]
                    * aa**2
                    * bb**2
                    * (cc - aa * bb) ** 2
                    * determinant
                )
    return value


def exterior_generator(generator, degree):
    indices = list(itertools.combinations(range(10), degree))
    index = {entry: i for i, entry in enumerate(indices)}
    out = np.zeros((len(indices), len(indices)))
    for column, source_tuple in enumerate(indices):
        for position, source in enumerate(source_tuple):
            for target in range(10):
                value = generator[target, source]
                if abs(value) < 1e-15:
                    continue
                target_tuple = list(source_tuple)
                target_tuple[position] = target
                if len(set(target_tuple)) < degree:
                    continue
                inversions = sum(
                    target_tuple[a] > target_tuple[b]
                    for a in range(degree)
                    for b in range(a + 1, degree)
                )
                out[index[tuple(sorted(target_tuple))], column] += (-1) ** inversions * value
    return out


def symmetric_basis(dimension):
    out = []
    for i in range(dimension):
        E = np.zeros((dimension, dimension))
        E[i, i] = 1
        out.append(E)
    for i in range(dimension):
        for j in range(i + 1, dimension):
            E = np.zeros((dimension, dimension))
            E[i, j] = E[j, i] = 1 / np.sqrt(2)
            out.append(E)
    return out


def invariant_exterior_data(degree):
    generators = [exterior_generator(g, degree) for g in rep.gens]
    casimir = -sum(g @ g for g in generators)
    eigenvalues, eigenvectors = np.linalg.eigh(casimir)
    full_basis = []
    local_blocks = []
    labels = []
    for ell in range(16):
        E = eigenvectors[:, np.abs(eigenvalues - ell * (ell + 1)) < 2e-7]
        if E.shape[1] == 0:
            continue
        restricted = [E.T @ g @ E for g in generators]
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
            if np.linalg.norm(equations) < 1e-9
            else null_space(equations, rcond=1e-9)
        )
        local = []
        for column in range(kernel.shape[1]):
            B = sum(kernel[row, column] * symmetric[row] for row in range(len(symmetric)))
            B = (B + B.T) / 2
            local.append(B)
            full_basis.append(E @ B @ E.T)
        local_blocks.append(local)
        labels.append((ell, E.shape[1], len(local)))
    print("exterior", degree, labels, flush=True)
    return full_basis, local_blocks


def solve(samples=260, seed=20260820):
    basis2, blocks2 = invariant_exterior_data(2)
    basis3, blocks3 = invariant_exterior_data(3)
    rng = np.random.default_rng(seed)
    rows2, rows3, targets = [], [], []
    for _ in range(samples):
        y = rng.normal(size=28)
        H = np.einsum("abk,k->ab", moments.maps[3], y)
        mass = np.trace(H)
        C2 = compound(H, PAIRS)
        C3 = compound(H, TRIPLES)
        rows2.append(mass * np.array([np.sum(B * C2) for B in basis2]))
        rows3.append(np.array([np.sum(B * C3) for B in basis3]))
        targets.append(mass * j_value(y) - 108 * f_value(y))
    rows2, rows3, targets = np.stack(rows2), np.stack(rows3), np.array(targets)
    raw = np.concatenate([rows2, rows3], axis=1)
    least, *_ = np.linalg.lstsq(raw, targets, rcond=None)
    print(
        "span",
        raw.shape,
        np.linalg.matrix_rank(raw),
        "relative residual",
        np.linalg.norm(raw @ least - targets) / np.linalg.norm(targets),
        flush=True,
    )

    variables2 = [cp.Variable(len(block)) for block in blocks2]
    variables3 = [cp.Variable(len(block)) for block in blocks3]
    constraints = []
    for variable, block in zip(variables2, blocks2):
        constraints.append(sum(variable[k] * block[k] for k in range(len(block))) >> 0)
    for variable, block in zip(variables3, blocks3):
        constraints.append(sum(variable[k] * block[k] for k in range(len(block))) >> 0)
    vector2 = cp.hstack(variables2)
    vector3 = cp.hstack(variables3)
    constraints.append(rows2 @ vector2 + rows3 @ vector3 == targets)
    objective = 0
    for variable, block in zip(variables2, blocks2):
        objective += cp.trace(sum(variable[k] * block[k] for k in range(len(block))))
    for variable, block in zip(variables3, blocks3):
        objective += cp.trace(sum(variable[k] * block[k] for k in range(len(block))))
    problem = cp.Problem(cp.Minimize(objective), constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                result = problem.solve(
                    solver=solver,
                    tol_gap_abs=2e-9,
                    tol_feas=2e-9,
                    tol_gap_rel=2e-9,
                    max_iter=3000,
                )
            else:
                result = problem.solve(solver=solver, eps=5e-7, max_iters=500000)
            print(solver, problem.status, result, flush=True)
            if problem.status in ("optimal", "optimal_inaccurate"):
                coefficients2 = np.concatenate([v.value for v in variables2])
                coefficients3 = np.concatenate([v.value for v in variables3])
                errors = []
                for _ in range(200):
                    y = rng.normal(size=28)
                    H = np.einsum("abk,k->ab", moments.maps[3], y)
                    mass = np.trace(H)
                    C2, C3 = compound(H, PAIRS), compound(H, TRIPLES)
                    r2 = mass * np.array([np.sum(B * C2) for B in basis2])
                    r3 = np.array([np.sum(B * C3) for B in basis3])
                    target = mass * j_value(y) - 108 * f_value(y)
                    errors.append(r2 @ coefficients2 + r3 @ coefficients3 - target)
                print("validation", np.max(np.abs(errors)), flush=True)
                np.savez(
                    "research/tensor_weighted_f_wedge3_solution.npz",
                    coefficients2=coefficients2,
                    coefficients3=coefficients3,
                )
                break
        except cp.error.SolverError as error:
            print(solver, error, flush=True)


def audit_normalization(seed=77):
    rng = np.random.default_rng(seed)
    xs = rng.normal(size=(7, 3))
    xs /= np.linalg.norm(xs, axis=1)[:, None]
    weights = rng.dirichlet(np.ones(len(xs)))
    y = sum(w * degree_six_evaluation(x) for w, x in zip(weights, xs))
    print("W map residual", W_RESIDUAL)
    print("F compound/direct", f_value(y), direct_f(xs, weights))
    print("mass/J", np.trace(np.einsum("abk,k->ab", moments.maps[3], y)), j_value(y))


if __name__ == "__main__":
    audit_normalization()
    solve()
