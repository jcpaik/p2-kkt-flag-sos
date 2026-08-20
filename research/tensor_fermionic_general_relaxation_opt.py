"""Test the full linear fermionic relaxation for tangent Veronese planes.

For a genuine tangent-plane mixture, the spin-2 part of the one-particle
marginal F is an equivariant linear image of the H1 block of the Pluecker
density G.  This script fits that exact (up to floating point basis choices)
linear map from orbit samples, then minimizes

    ||G||^2 - ||F||^2/2 + 1/3

over arbitrary PSD G of trace one subject only to

    tr G_11 = 1/5,
    Pi_2(F) = L(G_11).

If this relaxation is nonnegative it is a promising abstract fermionic
lemma that would prove the unrestricted problem.
"""

import sys

sys.path.insert(0, "research")

import numpy as np
import torch
import cvxpy as cp

import tensor_fermionic_relaxation_opt as base
import tensor_wedge_relation as wedge


torch.set_default_dtype(torch.float64)

T = np.hstack([wedge.W1, wedge.W3])


def sym2_coordinates(S):
    """Coordinates in wedge.comps2's orthonormal symmetric-tensor basis."""
    out = []
    for alpha in wedge.comps2:
        inds = [i for i, count in enumerate(alpha) for _ in range(count)]
        if inds[0] == inds[1]:
            out.append(S[inds[0], inds[1]])
        else:
            out.append(np.sqrt(2) * S[inds[0], inds[1]])
    return np.array(out)


def sym2_matrix(coordinates):
    out = np.zeros((3, 3))
    for value, alpha in zip(coordinates, wedge.comps2):
        inds = [i for i, count in enumerate(alpha) for _ in range(count)]
        if inds[0] == inds[1]:
            out[inds[0], inds[1]] = value
        else:
            out[inds[0], inds[1]] = value / np.sqrt(2)
            out[inds[1], inds[0]] = value / np.sqrt(2)
    return out


H2_MATRICES = [sym2_matrix(wedge.U2[:, a]) for a in range(5)]


def tangent_data(x):
    seed = np.eye(3)[np.argmin(np.abs(x))]
    u = np.cross(x, seed)
    u /= np.linalg.norm(u)
    v = np.cross(x, u)
    cols = []
    for y in (u, v):
        S = (np.outer(x, y) + np.outer(y, x)) / np.sqrt(2)
        cols.append(wedge.U2.T @ sym2_coordinates(S))
    W = np.stack(cols, axis=1)
    Q = W @ W.T
    z = np.array(
        [W[i, 0] * W[j, 1] - W[i, 1] * W[j, 0] for i, j in wedge.pairs]
    )
    z /= np.linalg.norm(z)
    zb = T.T @ z
    Rb = np.outer(zb, zb)
    return Q, Rb


def marginal_upper(x):
    """P_x with <S,P_x S>=2|Sx|^2; pointwise P_x-Q_x is PSD."""
    vectors = np.stack([S @ x for S in H2_MATRICES])
    return 2 * vectors @ vectors.T


def fit_relation(samples=300, seed=5):
    rng = np.random.default_rng(seed)
    xs = rng.normal(size=(samples, 3))
    xs /= np.linalg.norm(xs, axis=1)[:, None]
    domain = []
    target = []
    for x in xs:
        Q, Rb = tangent_data(x)
        domain.append(Rb[:3, :3].reshape(-1))
        target.append(base.spin2_part(torch.tensor(Q)).numpy().reshape(-1))
    domain = np.stack(domain)
    target = np.stack(target)
    linear, *_ = np.linalg.lstsq(domain, target, rcond=None)
    residual = np.linalg.norm(domain @ linear - target) / np.linalg.norm(target)
    return linear.T, residual


RELATION_NP, RELATION_RESIDUAL = fit_relation()
RELATION = torch.tensor(RELATION_NP)


def fit_upper_relation(samples=300, seed=6):
    rng = np.random.default_rng(seed)
    xs = rng.normal(size=(samples, 3))
    xs /= np.linalg.norm(xs, axis=1)[:, None]
    domain = []
    target = []
    for x in xs:
        _, Rb = tangent_data(x)
        domain.append(Rb[:3, :3].reshape(-1))
        target.append(marginal_upper(x).reshape(-1))
    domain = np.stack(domain)
    target = np.stack(target)
    linear, *_ = np.linalg.lstsq(domain, target, rcond=None)
    residual = np.linalg.norm(domain @ linear - target) / np.linalg.norm(target)
    return linear.T, residual


UPPER_NP, UPPER_RESIDUAL = fit_upper_relation()


def fit_cross_span(samples=300, seed=7):
    rng = np.random.default_rng(seed)
    xs = rng.normal(size=(samples, 3))
    xs /= np.linalg.norm(xs, axis=1)[:, None]
    rows = np.stack([tangent_data(x)[1][:3, 3:].reshape(-1) for x in xs])
    _, singular, vh = np.linalg.svd(rows, full_matrices=True)
    rank = int(np.sum(singular > 1e-10))
    return vh[:rank], vh[rank:], singular


CROSS_SPAN, CROSS_NULL, CROSS_SINGULAR = fit_cross_span()


def casimir_operator_block(X, generators):
    return -sum(g @ (g @ X - X @ g) - (g @ X - X @ g) @ g for g in generators)


def project_operator_spin(X, generators, ell, spins):
    out = X.copy()
    for other in spins:
        if other != ell:
            out = (casimir_operator_block(out, generators) - other * (other + 1) * out) / (
                ell * (ell + 1) - other * (other + 1)
            )
    return out


def fit_diagonal_from_cross(ell, samples=500, seed=11):
    rng = np.random.default_rng(seed + ell)
    xs = rng.normal(size=(samples, 3))
    xs /= np.linalg.norm(xs, axis=1)[:, None]
    domain = []
    target = []
    for x in xs:
        _, Rb = tangent_data(x)
        domain.append(Rb[:3, 3:].reshape(-1))
        target.append(
            project_operator_spin(Rb[3:, 3:], wedge.gw3, ell, (0, 2, 4, 6)).reshape(-1)
        )
    domain = np.stack(domain)
    target = np.stack(target)
    linear, *_ = np.linalg.lstsq(domain, target, rcond=None)
    residual = np.linalg.norm(domain @ linear - target) / np.linalg.norm(target)
    return linear.T, residual


DIAGONAL4_FROM_CROSS, DIAGONAL4_RESIDUAL = fit_diagonal_from_cross(4)

# Linear projector on the spin-4 sector of symmetric End(H3).
PROJECT4_33 = np.zeros((49, 49))
for a in range(7):
    for b in range(7):
        E = np.zeros((7, 7))
        E[a, b] = 1
        PROJECT4_33[:, 7 * a + b] = project_operator_spin(
            E, wedge.gw3, 4, (0, 2, 4, 6)
        ).reshape(-1)


def independent_symmetric_rows(linear_map, block_shape=(10, 10), tolerance=1e-9):
    expansion = []
    rows, cols = block_shape
    for a in range(rows):
        for b in range(a, cols):
            E = np.zeros(block_shape)
            E[a, b] = 1
            E[b, a] = 1
            if a == b:
                E[a, b] = 1
            expansion.append(E.reshape(-1))
    expansion = np.stack(expansion, axis=1)
    u, singular, _ = np.linalg.svd(linear_map @ expansion, full_matrices=False)
    rank = int(np.sum(singular > tolerance))
    return u[:, :rank].T @ linear_map, singular


# Build the nine independent l=4 matching constraints on the full block.
L4_FULL = np.zeros((49, 100))
for a in range(7):
    for b in range(7):
        L4_FULL[:, 10 * (a + 3) + (b + 3)] = PROJECT4_33[:, 7 * a + b]
for a in range(3):
    for b in range(7):
        L4_FULL[:, 10 * a + (b + 3)] -= DIAGONAL4_FROM_CROSS[:, 7 * a + b]
L4_REDUCED, L4_SINGULAR = independent_symmetric_rows(L4_FULL)

# Linear coordinate maps in the block basis, used by the constrained CCCP
# audit below.
TBASIS = np.zeros((100, 100))
for a in range(10):
    for b in range(10):
        E = np.zeros((10, 10))
        E[a, b] = 1
        TBASIS[:, 10 * a + b] = (T @ E @ T.T).reshape(-1)
F_LINEAR = base.GAMMA.numpy() @ TBASIS


def spin2_numpy(X):
    return base.spin2_part(torch.tensor(X)).numpy()


SPIN2_LINEAR = np.zeros((25, 25))
for a in range(5):
    for b in range(5):
        E = np.zeros((5, 5))
        E[a, b] = 1
        SPIN2_LINEAR[:, 5 * a + b] = spin2_numpy(E).reshape(-1)

RELATION_FULL = np.zeros((25, 100))
for a in range(3):
    for b in range(3):
        RELATION_FULL[:, 10 * a + b] = RELATION_NP[:, 3 * a + b]
LINEAR_CONSTRAINT = SPIN2_LINEAR @ F_LINEAR - RELATION_FULL

# Remove redundant rows after restricting the block variable to symmetric
# matrices.  Keeping all 25 floating-point rows makes conic solvers interpret
# numerical null rows as inconsistent equalities.
SYM_EXPANSION = []
for a in range(10):
    for b in range(a, 10):
        E = np.zeros((10, 10))
        E[a, b] = 1
        E[b, a] = 1
        if a == b:
            E[a, b] = 1
        SYM_EXPANSION.append(E.reshape(-1))
SYM_EXPANSION = np.stack(SYM_EXPANSION, axis=1)
u_constraint, s_constraint, _ = np.linalg.svd(
    LINEAR_CONSTRAINT @ SYM_EXPANSION, full_matrices=False
)
CONSTRAINT_RANK = int(np.sum(s_constraint > 1e-9))
LINEAR_CONSTRAINT_REDUCED = u_constraint[:, :CONSTRAINT_RANK].T @ LINEAR_CONSTRAINT


def contraction(block):
    G = torch.tensor(T) @ block @ torch.tensor(T).T
    F = (base.GAMMA @ G.reshape(-1)).reshape(5, 5)
    return G, F


def objective(raw, penalty=1e4):
    block0 = raw @ raw.T
    block = block0 / torch.trace(block0)
    G, F = contraction(block)
    relation_error = base.spin2_part(F).reshape(-1) - RELATION @ block[:3, :3].reshape(-1)
    constraints = (torch.trace(block[:3, :3]) - torch.tensor(1 / 5)) ** 2
    constraints = constraints + torch.sum(relation_error**2)
    gap = torch.sum(G**2) - torch.sum(F**2) / 2 + torch.tensor(1 / 3)
    return gap + penalty * constraints, gap, constraints, block, F


def optimize(restarts=60, steps=12000, learning_rate=0.01, penalty=2e5, seed=8):
    rng = torch.Generator().manual_seed(seed)
    best = None
    for restart in range(restarts):
        raw = torch.randn((10, 10), generator=rng, requires_grad=True)
        opt = torch.optim.Adam([raw], lr=learning_rate)
        for _ in range(steps):
            opt.zero_grad()
            loss, *_ = objective(raw, penalty)
            loss.backward()
            opt.step()
        with torch.no_grad():
            _, gap, violation, block, F = objective(raw, penalty)
            record = (
                float(gap),
                float(violation),
                restart,
                np.linalg.eigvalsh(block.numpy()),
                np.linalg.eigvalsh(F.numpy()),
                float(torch.trace(block[:3, :3])),
            )
            if best is None or record[0] < best[0]:
                best = record
                print("best", best, flush=True)
    return best


def orbit_mixture(xs, weights=None):
    if weights is None:
        weights = np.ones(len(xs)) / len(xs)
    return sum(w * tangent_data(x)[1] for x, w in zip(xs, weights))


def cccp(
    initial,
    iterations=40,
    verbose=False,
    use_upper=True,
    use_spin2=True,
    use_cross_parity=True,
    use_l4=True,
):
    """Convex-concave iterations with the linear constraints imposed exactly."""
    X = cp.Variable((10, 10), symmetric=True)
    xvec = cp.reshape(X, (100,), order="C")
    fvec = F_LINEAR @ xvec
    upper_vec = UPPER_NP @ cp.reshape(X[:3, :3], (9,), order="C")
    upper = cp.reshape(upper_vec, (5, 5), order="C")
    fmatrix = cp.reshape(fvec, (5, 5), order="C")
    constraints = [
        X >> 0,
        cp.trace(X) == 1,
        cp.trace(X[:3, :3]) == cp.Constant(1 / 5),
    ]
    if use_upper:
        constraints.append(upper - fmatrix >> 0)
    if use_spin2:
        constraints.append(LINEAR_CONSTRAINT_REDUCED @ xvec == 0)
    if use_cross_parity:
        constraints.append(
            CROSS_NULL @ cp.reshape(X[:3, 3:], (21,), order="C") == 0
        )
    if use_l4:
        constraints.append(L4_REDUCED @ xvec == 0)
    current = initial.copy()
    for it in range(iterations):
        f0 = (F_LINEAR @ current.reshape(-1)).reshape(5, 5)
        # Write the linearized concave term as a matrix Frobenius pairing.
        # This avoids a CVXPY canonicalization bug seen for C-order flattened
        # symmetric variables under Python 3.14.
        linear_matrix = (F_LINEAR.T @ f0.reshape(-1)).reshape(10, 10)
        center = (linear_matrix + linear_matrix.T) / 4
        # Equivalent to ||X||^2-<linear_matrix,X>, up to a constant.
        objective = cp.Minimize(cp.sum_squares(X - center))
        problem = cp.Problem(objective, constraints)
        try:
            value = problem.solve(
                solver="CLARABEL",
                tol_gap_abs=1e-9,
                tol_feas=1e-9,
                tol_gap_rel=1e-9,
                max_iter=1000,
            )
        except cp.error.SolverError:
            value = None
        if X.value is None:
            value = problem.solve(solver="SCS", eps=2e-8, max_iters=200000)
        current = X.value
        F = (F_LINEAR @ current.reshape(-1)).reshape(5, 5)
        true_gap = np.sum(current**2) - np.sum(F**2) / 2 + 1 / 3
        if verbose:
            print("cccp", it, value, true_gap, np.linalg.eigvalsh(current))
    return true_gap, current, F


def cccp_restarts(restarts=20, iterations=30, seed=100):
    rng = np.random.default_rng(seed)
    starts = []
    # Haar-like feasible point from a large orbit sample.
    xs = rng.normal(size=(1000, 3))
    xs /= np.linalg.norm(xs, axis=1)[:, None]
    starts.append(orbit_mixture(xs))
    # Genuine random orbit mixtures are automatically feasible.
    for _ in range(restarts - 1):
        xs = rng.normal(size=(20, 3))
        xs /= np.linalg.norm(xs, axis=1)[:, None]
        starts.append(orbit_mixture(xs, rng.dirichlet(np.ones(len(xs)))))
    best = None
    for i, start in enumerate(starts):
        result = cccp(start, iterations=iterations)
        if best is None or result[0] < best[0]:
            best = result
            print(
                "cccp best",
                i,
                best[0],
                np.linalg.eigvalsh(best[1]),
                np.linalg.eigvalsh(best[2]),
                flush=True,
            )
    return best


def extreme_restarts(restarts=40, iterations=120, seed=101):
    """Start CCCP at random extreme points of the minimal 40D relaxation."""
    rng = np.random.default_rng(seed)
    X = cp.Variable((10, 10), symmetric=True)
    xvec = cp.reshape(X, (100,), order="C")
    direction = cp.Parameter((10, 10), symmetric=True)
    constraints = [
        X >> 0,
        cp.trace(X) == 1,
        cp.trace(X[:3, :3]) == cp.Constant(1 / 5),
        LINEAR_CONSTRAINT_REDUCED @ xvec == 0,
        L4_REDUCED @ xvec == 0,
    ]
    problem = cp.Problem(cp.Minimize(cp.sum(cp.multiply(direction, X))), constraints)
    best = None
    for restart in range(restarts):
        raw = rng.normal(size=(10, 10))
        direction.value = (raw + raw.T) / 2
        try:
            problem.solve(solver="CLARABEL")
        except cp.error.SolverError:
            problem.solve(solver="SCS", eps=1e-7, max_iters=200000)
        result = cccp(
            X.value,
            iterations=iterations,
            use_upper=False,
            use_spin2=True,
            use_cross_parity=False,
            use_l4=True,
        )
        if best is None or result[0] < best[0]:
            best = result
            print(
                "extreme best",
                restart,
                result[0],
                np.linalg.eigvalsh(result[1]),
                np.linalg.eigvalsh(result[2]),
                flush=True,
            )
    return best


if __name__ == "__main__":
    print("relation residual", RELATION_RESIDUAL, "upper residual", UPPER_RESIDUAL)
    # Audit pointwise orbit normalization.
    q, rb = tangent_data(np.array([0.2, -0.3, np.sqrt(0.87)]))
    print("orbit block traces", np.trace(rb[:3, :3]), np.trace(rb[3:, 3:]))
    print(
        "orbit relation error",
        np.linalg.norm(
            base.spin2_part(torch.tensor(q)).numpy().reshape(-1)
            - RELATION_NP @ rb[:3, :3].reshape(-1)
        ),
    )
    cccp_restarts()
