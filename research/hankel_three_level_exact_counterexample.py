"""Exact counterexample to the trace-ratio filtered-bosonic relaxation.

This does *not* counterexample the fully Hankel conjecture.  It shows that
positivity of the cubic state and its ordinary bosonic reductions, together
with the H1:H3 trace ratio of the Pluecker image, is insufficient.

All calculations are over Q(sqrt(2), sqrt(3), sqrt(5), ...), through SymPy.
"""

import math

import sympy as sp


Q = sp.Rational
sqrt = sp.sqrt


def compositions(n, k=3):
    if k == 1:
        return [(n,)]
    return [(a,) + rest for a in range(n + 1) for rest in compositions(n - a, k - 1)]


C3 = compositions(3)
C2 = compositions(2)
C1 = compositions(1)


def pluecker_filter():
    r"""The exact matrix A with z_x=A x^{\otimes 3} in normalized bases."""
    A = sp.zeros(10)
    r3 = 1 / sqrt(3)
    A[0, 5] = sqrt(2)
    A[1, 2] = A[1, 7] = r3
    A[2, 1], A[2, 8] = r3, -2 * r3
    A[3, 4], A[3, 6] = r3, -2 * r3
    A[4, 2], A[4, 7] = -1, 1
    A[5, 1] = -1
    A[6, 4] = 1
    A[7, 4], A[7, 6], A[7, 9] = -r3, -r3, 1
    A[8, 1], A[8, 3], A[8, 8] = r3, -1, r3
    A[9, 0], A[9, 2], A[9, 7] = 1, -r3, -r3
    return A


A = pluecker_filter()


def reduce_one(X, upper, lower, n):
    """Normalized bosonic one-particle partial trace Sym^n -> Sym^(n-1)."""
    out = sp.zeros(len(lower))
    for a, beta in enumerate(lower):
        for b, gamma in enumerate(lower):
            for i in range(3):
                alpha = list(beta)
                delta = list(gamma)
                alpha[i] += 1
                delta[i] += 1
                out[a, b] += (
                    sqrt((beta[i] + 1) * (gamma[i] + 1))
                    * X[upper.index(tuple(alpha)), upper.index(tuple(delta))]
                    / n
                )
    return sp.simplify(out)


trace_vector = sp.zeros(6, 1)
for i, alpha in enumerate(C2):
    if max(alpha) == 2:
        trace_vector[i] = 1 / sqrt(3)
P2 = sp.eye(6) - trace_vector * trace_vector.T


def levels(H):
    G = sp.simplify(A * H * A.T)
    R2 = reduce_one(H, C3, C2, 3)
    M = reduce_one(R2, C2, C1, 2)
    # This is L represented on the full Sym^2 space.  Its Frobenius norm is
    # the same as that of the 5-by-5 traceless-quadratic representation.
    L = sp.simplify(2 * P2 * R2 * P2)
    return G, R2, M, L, sp.trace(H)


def frobenius(X, Y):
    return sp.trace(X.T * Y)


e = [sp.eye(10)[:, i] for i in range(10)]

# SO(2)-adapted unit vectors about the first coordinate axis.
u0 = (e[4] + e[6]) / sqrt(2)
v0 = (u0 + Q(7, 2) * sqrt(6) * e[9]) / sqrt(Q(149, 2))

u1c = (sqrt(3) * e[0] + e[2]) / 2
u1s = (sqrt(3) * e[3] + e[1]) / 2
v1c = (3 * u1c - e[7]) / sqrt(10)
v1s = (3 * u1s - e[8]) / sqrt(10)

v3c = (e[0] - sqrt(3) * e[2]) / 2
v3s = (e[3] - sqrt(3) * e[1]) / 2

H0 = v0 * v0.T
H1 = v1c * v1c.T + v1s * v1s.T
H3 = v3c * v3c.T + v3s * v3s.T

a = Q(3129, 6260)
b = Q(7671, 31300)
c = Q(1, 200)
H = sp.simplify(a * H0 + b * H1 + c * H3)
G, R2, M, L, mass = levels(H)

# Exact projector onto the spin-1 block of wedge^2(H2), in the edge basis.
P1 = sp.zeros(10)
P1[1, 1], P1[1, 9], P1[9, 1], P1[9, 9] = Q(4, 5), Q(2, 5), Q(2, 5), Q(1, 5)
for indices, signs in [((2, 5, 8), (1, 1)), ((3, 6, 7), (-1, -1))]:
    i, j, k = indices
    s1, s2 = signs
    P1[i, i], P1[j, j], P1[k, k] = Q(1, 5), Q(3, 5), Q(1, 5)
    P1[i, j] = P1[j, i] = s1 * sqrt(3) / 5
    P1[i, k] = P1[k, i] = s2 * Q(1, 5)
    P1[j, k] = P1[k, j] = sqrt(3) / 5

gap = sp.factor(
    frobenius(G, G)
    - frobenius(L, L) / 2
    + frobenius(M, M) / 6
    + mass**2 / 18
)


if __name__ == "__main__":
    print("weights", a, b, c)
    print("PSD decomposition ranks", 1, 2, 2)
    print("tr H, tr G, tr(P1 G)", mass, sp.trace(G), sp.trace(P1 * G))
    print("||G||^2", sp.factor(frobenius(G, G)))
    print("||L||^2", sp.factor(frobenius(L, L)))
    print("||M||^2", sp.factor(frobenius(M, M)))
    print("gap", gap, float(gap))
    # Fully Hankel would require H_55/6 = H_46/3, since both entries
    # represent the degree-six moment with multi-index (2,2,2).
    print("H[5,5], H[4,6], 2H[4,6]-H[5,5]", H[5, 5], H[4, 6], 2 * H[4, 6] - H[5, 5])
