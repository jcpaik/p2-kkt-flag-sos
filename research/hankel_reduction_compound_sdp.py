"""Search an exact degree-two certificate using all reduction-map PSD images.

For a fully symmetric/partial-transpose-invariant three-qutrit state H, both
  D3 = P_sym(I tensor R2)P_sym - H,
  D2 = P_sym(I tensor R1)P_sym - R2
are PSD by the reduction criterion.  We test whether Q is a sum of positive
linear functionals on the second compounds of H,D3,R2,D2,R1, plus squares of
linear moment forms.
"""
import cvxpy as cp
import numpy as np

from tensor_hankel_fw import maps
from hankel_reduction_identities import lift32, lift21


def polarized_wedge(A, B):
    n = A.shape[0]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    out = np.empty((len(pairs), len(pairs)))
    for row, (i, j) in enumerate(pairs):
        for col, (k, ell) in enumerate(pairs):
            out[row, col] = 0.5 * (
                A[i, k] * B[j, ell]
                + B[i, k] * A[j, ell]
                - A[i, ell] * B[j, k]
                - B[i, ell] * A[j, k]
            )
    return out


Hmap = np.moveaxis(maps[3], 2, 0)
R2map = np.moveaxis(maps[2], 2, 0)
R1map = np.moveaxis(maps[1], 2, 0)
D3map = np.array([lift32(R2map[k]) - Hmap[k] for k in range(28)])
D2map = np.array([lift21(R1map[k]) - R2map[k] for k in range(28)])
linear_maps = {
    "H": Hmap,
    "D3": D3map,
    "R2": R2map,
    "D2": D2map,
    "R1": R1map,
}

grams = {d: np.einsum("abk,abl->kl", M, M) for d, M in maps.items()}
mass = np.trace(maps[3], axis1=0, axis2=1)
target = 32 * grams[3] - 48 * grams[2] + 20 * grams[1] - 4 / 3 * np.outer(mass, mass)

weights = {}
constraints = []
for name, M in linear_maps.items():
    size = M.shape[1] * (M.shape[1] - 1) // 2
    weights[name] = cp.Variable((size, size), symmetric=True, name="W_" + name)
    constraints.append(weights[name] >> 0)

S = cp.Variable((28, 28), symmetric=True, name="S")
constraints.append(S >> 0)

for left in range(28):
    for right in range(left, 28):
        expression = S[left, right]
        for name, M in linear_maps.items():
            coefficient = polarized_wedge(M[left], M[right])
            expression += cp.sum(cp.multiply(weights[name], coefficient))
        constraints.append(expression == target[left, right])

objective = cp.Minimize(cp.trace(S) + sum(cp.trace(W) for W in weights.values()))
problem = cp.Problem(objective, constraints)
for solver in (cp.CLARABEL, cp.SCS):
    try:
        value = problem.solve(solver=solver, verbose=False, max_iter=200000)
        print(solver, problem.status, value)
        if S.value is not None:
            print("S min eig", np.linalg.eigvalsh(S.value).min())
            for name, W in weights.items():
                print(name, "min eig", np.linalg.eigvalsh(W.value).min(), "trace", np.trace(W.value))
            print("max eq residual", max(float(np.max(np.abs(c.violation()))) for c in constraints[6:]))
    except (cp.error.SolverError, TypeError) as error:
        print(solver, error)
