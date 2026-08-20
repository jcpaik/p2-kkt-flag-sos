"""Search E(H)=sum_r tr(H A_r H A_r^T) on the Hankel subspace."""

import cvxpy as cp
import numpy as np

from tensor_hankel_fw import maps


M = maps[3]
gram = {d: np.einsum("abk,abl->kl", X, X) for d, X in maps.items()}
mass = np.trace(M, axis1=0, axis2=1)
target = 32*gram[3] - 48*gram[2] + 20*gram[1] - (4/3)*np.outer(mass, mass)

C = cp.Variable((100, 100), symmetric=True)
constraints = [C >> 0]
for alpha in range(28):
    A = M[:, :, alpha]
    for beta in range(alpha, 28):
        B = M[:, :, beta]
        # D[(j,k),(i,l)] contracts with C[(j,k),(i,l)].
        D4 = 0.5 * (
            np.einsum("ij,kl->jkil", A, B)
            + np.einsum("ij,kl->jkil", B, A)
        )
        D = D4.reshape(100, 100)
        D = 0.5*(D+D.T)
        constraints.append(cp.sum(cp.multiply(C, D)) == target[alpha, beta])

problem = cp.Problem(cp.Minimize(cp.trace(C)), constraints)
for solver in (cp.CLARABEL, cp.SCS):
    try:
        value = problem.solve(solver=solver, verbose=False)
        print(solver, problem.status, value)
        if C.value is not None:
            eig = np.linalg.eigvalsh(C.value)
            print("eig", eig[:20], eig[-20:])
            print("res", max(np.max(np.abs(c.violation())) for c in constraints[1:]))
            np.save("/tmp/tensor_hankel_cpmap_C.npy", C.value)
    except cp.error.SolverError as error:
        print(solver, error)
