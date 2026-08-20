"""Test E >= (32/105) det((5/4)(I-A2))^2 on Hankel pseudo-moments."""
import cvxpy as cp
import numpy as np

from tensor_hankel_fw import maps, matrices, energy, uniform_moment, deg6


def traceless_basis():
    out = []
    out.append(np.diag([1.0, -1.0, 0.0]) / np.sqrt(2))
    out.append(np.diag([1.0, 1.0, -2.0]) / np.sqrt(6))
    for i, j in ((0, 1), (0, 2), (1, 2)):
        S = np.zeros((3, 3))
        S[i, j] = S[j, i] = 1 / np.sqrt(2)
        out.append(S)
    return out


TB = traceless_basis()

# The degree-one monomial basis used by maps is (z,y,x), while matrices
# below use Cartesian order (x,y,z).
P1 = np.eye(3)[::-1]

# The degree-two normalized monomial basis is
# (z^2,sqrt2 yz,y^2,sqrt2 xz,sqrt2 xy,x^2).
def sym2_coordinates(S):
    return np.array([S[2, 2], np.sqrt(2) * S[1, 2], S[1, 1],
                     np.sqrt(2) * S[0, 2], np.sqrt(2) * S[0, 1], S[0, 0]])


TB2 = [sym2_coordinates(S) for S in TB]


def spin2_gap(y):
    rho = matrices(y)
    R2 = rho[2]
    R1 = P1 @ rho[1] @ P1
    B = np.empty((5, 5))
    for a, Sa in enumerate(TB):
        for b, Sb in enumerate(TB):
            B[a, b] = (
                2 * np.trace(R1 @ (Sa @ Sb + Sb @ Sa))
                - 4 * TB2[a] @ R2 @ TB2[b]
            )
    determinant = np.linalg.det(1.25 * B)
    return energy(y) - (32 / 105) * determinant * determinant, B, determinant


def gradient_fd(y):
    # det(B)^2 derivative is computed stably by finite differences only for
    # exploratory Frank-Wolfe; all reported identities use direct evaluation.
    out = np.empty(28)
    h = 2e-6
    for k in range(28):
        step = h * max(1.0, abs(y[k]))
        e = np.zeros(28)
        e[k] = step
        out[k] = (spin2_gap(y + e)[0] - spin2_gap(y - e)[0]) / (2 * step)
    return out


normalizer = np.trace(maps[3], axis1=0, axis2=1)
yvar = cp.Variable(28)
hankel = sum(yvar[k] * maps[3][:, :, k] for k in range(28))
direction = cp.Parameter(28)
oracle = cp.Problem(cp.Minimize(direction @ yvar), [hankel >> 0, normalizer @ yvar == 1])


def main():
    uniform = np.array([uniform_moment(a) for a in deg6])
    print("uniform", spin2_gap(uniform), "energy", energy(uniform))
    rng = np.random.default_rng(9173)
    feasible = [uniform]
    # Collect exposed Hankel rays and random convex mixtures.
    for _ in range(250):
        direction.value = rng.normal(size=28)
        oracle.solve(solver="CLARABEL")
        feasible.append(yvar.value.copy())
    for _ in range(2000):
        inds = rng.choice(len(feasible), size=rng.integers(2, 12), replace=True)
        weights = rng.dirichlet(np.ones(len(inds)) * 0.35)
        feasible.append(sum(w * feasible[i] for w, i in zip(weights, inds)))
    values = np.array([spin2_gap(y)[0] for y in feasible])
    j = int(np.argmin(values))
    print("random min", values[j], "E/B/det", energy(feasible[j]), spin2_gap(feasible[j])[1:])

    # Difference-of-convex-style Frank-Wolfe starts with exact line scans.
    for start in range(20):
        y = feasible[int(rng.integers(len(feasible)))].copy()
        best = spin2_gap(y)[0]
        for iteration in range(150):
            direction.value = gradient_fd(y)
            oracle.solve(solver="CLARABEL")
            vertex = yvar.value.copy()
            # Global one-dimensional grid, then local golden search.
            grid = np.linspace(0, 1, 101)
            vals = np.array([spin2_gap((1 - t) * y + t * vertex)[0] for t in grid])
            k = int(np.argmin(vals))
            if vals[k] >= spin2_gap(y)[0] - 1e-10:
                break
            lo = grid[max(0, k - 1)]
            hi = grid[min(100, k + 1)]
            for _ in range(30):
                t1 = lo + (hi - lo) / 3
                t2 = hi - (hi - lo) / 3
                if spin2_gap((1 - t1) * y + t1 * vertex)[0] < spin2_gap((1 - t2) * y + t2 * vertex)[0]:
                    hi = t2
                else:
                    lo = t1
            t = (lo + hi) / 2
            y = (1 - t) * y + t * vertex
            best = min(best, spin2_gap(y)[0])
        print("FW", start, iteration, best, "E/det", energy(y), spin2_gap(y)[2],
              "rank", np.linalg.eigvalsh(matrices(y)[3]))


if __name__ == "__main__":
    main()
