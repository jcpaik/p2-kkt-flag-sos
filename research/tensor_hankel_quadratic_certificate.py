"""Search E = SOS in sextic moments + PSD-weighted 2x2 catalecticant minors."""

import cvxpy as cp
import numpy as np

from tensor_hankel_fw import maps


pairs = [(i, j) for i in range(10) for j in range(i + 1, 10)]


def wedge_bilinear(A, B):
    out = np.zeros((45, 45))
    for row, (i, j) in enumerate(pairs):
        for col, (k, ell) in enumerate(pairs):
            out[row, col] = 0.5 * (
                A[i, k] * B[j, ell]
                + B[i, k] * A[j, ell]
                - A[i, ell] * B[j, k]
                - B[i, ell] * A[j, k]
            )
    return out


gram = {d: np.einsum("abk,abl->kl", M, M) for d, M in maps.items()}
mass = np.trace(maps[3], axis1=0, axis2=1)
target = 32*gram[3] - 48*gram[2] + 20*gram[1] - (4/3)*np.outer(mass, mass)

W = cp.Variable((45, 45), symmetric=True)
R = cp.Variable((28, 28), symmetric=True)
constraints = [W >> 0, R >> 0]
for i in range(28):
    for j in range(i, 28):
        C = wedge_bilinear(maps[3][:, :, i], maps[3][:, :, j])
        constraints.append(R[i, j] + cp.sum(cp.multiply(W, C)) == target[i, j])

problem = cp.Problem(cp.Minimize(cp.trace(W)+cp.trace(R)), constraints)
for solver in (cp.CLARABEL, cp.SCS):
    try:
        value = problem.solve(solver=solver, verbose=False)
        print(solver, problem.status, value)
        if W.value is not None:
            print("W eig", np.linalg.eigvalsh(W.value)[:10])
            print("R eig", np.linalg.eigvalsh(R.value)[:10])
            print("res", max(np.max(np.abs(c.violation())) for c in constraints[2:]))
            np.savez("/tmp/tensor_hankel_qcert.npz", W=W.value, R=R.value)
    except cp.error.SolverError as error:
        print(solver, error)
