"""Exact D4-symmetric audit of the negative-eigenvector frame separator.

The spin-4 tensor

    H = H0/60 + H4/36,

with H0=35 z^4-30 z^2 r^2+3 r^4 and H4=Re(x+iy)^4,
is the exact orbit found by the adversarial frame-QP search.  This script
builds orthonormal STF tensor bases symbolically and computes Bhat, the cubic
Schur core M, and the spin-6 Gram data without floating-point fitting.
"""

import itertools

import sympy as sp


x, y, z = sp.symbols("x y z")
vars3 = (x, y, z)


def tensor_from_poly(poly, degree):
    P = sp.Poly(sp.expand(poly), *vars3)
    shape = (3,) * degree
    out = sp.MutableDenseNDimArray.zeros(*shape)
    for indices in itertools.product(range(3), repeat=degree):
        alpha = tuple(indices.count(i) for i in range(3))
        coefficient = P.coeff_monomial(x ** alpha[0] * y ** alpha[1] * z ** alpha[2])
        multinomial = sp.factorial(degree)
        for value in alpha:
            multinomial /= sp.factorial(value)
        out[indices] = coefficient / multinomial
    return out


def inner(A, B):
    return sp.simplify(
        sum(A[index] * B[index] for index in itertools.product(*[range(n) for n in A.shape]))
    )


def normalize_polys(polys, degree):
    tensors = [tensor_from_poly(p, degree) for p in polys]
    # These magnetic bases are already pairwise orthogonal.
    assert all(inner(tensors[i], tensors[j]) == 0 for i in range(len(tensors)) for j in range(i))
    return [sp.ImmutableDenseNDimArray(T / sp.sqrt(inner(T, T))) for T in tensors]


H2 = normalize_polys(
    [2 * z**2 - x**2 - y**2, x * z, y * z, x**2 - y**2, x * y], 2
)
H3 = normalize_polys(
    [
        2 * z**3 - 3 * z * (x**2 + y**2),
        x * (4 * z**2 - x**2 - y**2),
        y * (4 * z**2 - x**2 - y**2),
        z * (x**2 - y**2),
        2 * x * y * z,
        x**3 - 3 * x * y**2,
        3 * x**2 * y - y**3,
    ],
    3,
)

H0 = 35 * z**4 - 30 * z**2 * (x**2 + y**2 + z**2) + 3 * (x**2 + y**2 + z**2) ** 2
H4 = x**4 - 6 * x**2 * y**2 + y**4
H = sp.ImmutableDenseNDimArray(tensor_from_poly(H0 / 60 + H4 / 36, 4))


def c4_matrix(Htensor=H):
    return sp.Matrix(
        5,
        5,
        lambda a, b: sum(
            Htensor[i, j, k, ell] * H2[a][i, j] * H2[b][k, ell]
            for i, j, k, ell in itertools.product(range(3), repeat=4)
        ),
    )


def a4_matrix(Htensor=H):
    return sp.Matrix(
        7,
        7,
        lambda a, b: sp.Rational(9, 11)
        * sum(
            Htensor[i, j, p, q] * H3[a][i, j, k] * H3[b][p, q, k]
            for i, j, k, p, q in itertools.product(range(3), repeat=5)
        ),
    )


def b_matrix(Htensor=H):
    # B*:H3 -> H1, with H1 in the standard Cartesian basis.
    return sp.Matrix(
        7,
        3,
        lambda a, ell: sp.sqrt(sp.Rational(3, 5))
        * sum(
            Htensor[i, j, k, ell] * H3[a][i, j,k]
            for i, j, k in itertools.product(range(3), repeat=3)
        ),
    )


def contraction_A(V):
    return sp.Matrix(
        3,
        3,
        lambda a, b: sum(
            V[a, i, j] * V[b, i, j] for i, j in itertools.product(range(3), repeat=2)
        ),
    )


def contraction_C(V, W):
    return sp.Matrix(
        3,
        3,
        lambda i, j: sum(
            V[a, b, i] * W[a, b, j] for a, b in itertools.product(range(3), repeat=2)
        ),
    )


def spin6_gram(V, W):
    av = contraction_A(V)
    aw = contraction_A(W)
    c = contraction_C(V, W)
    vv = inner(V, V)
    ww = inner(W, W)
    vw = inner(V, W)
    return sp.factor(
        -sp.Rational(493, 462) * vv * ww
        + sp.Rational(115, 22) * sp.trace(av * aw)
        - sp.Rational(45, 11) * sp.trace(c.T * c)
        + vw**2
    )


def symmetric_product(V, W):
    out = sp.MutableDenseNDimArray.zeros(3, 3, 3, 3, 3, 3)
    for subset in itertools.combinations(range(6), 3):
        other = [i for i in range(6) if i not in subset]
        for index in itertools.product(range(3), repeat=6):
            out[index] += V[tuple(index[i] for i in subset)] * W[
                tuple(index[i] for i in other)
            ] / 20
    return sp.ImmutableDenseNDimArray(out)


def trace_tensor(A):
    degree = len(A.shape)
    out = sp.MutableDenseNDimArray.zeros(*((3,) * (degree - 2)))
    for index in itertools.product(range(3), repeat=degree - 2):
        out[index] = sum(A[(a, a) + index] for a in range(3))
    return sp.ImmutableDenseNDimArray(out)


def stf6_inner(A, B):
    traces_a = [A]
    traces_b = [B]
    for _ in range(3):
        traces_a.append(trace_tensor(traces_a[-1]))
        traces_b.append(trace_tensor(traces_b[-1]))
    coefficients = [1, -sp.Rational(15, 11), sp.Rational(5, 11), -sp.Rational(5, 231)]
    return sp.factor(sum(c * inner(a, b) for c, a, b in zip(coefficients, traces_a, traces_b)))


def report():
    C = c4_matrix()
    Bhat = sp.eye(5) - 5 * C
    A = a4_matrix()
    B = b_matrix()
    M = sp.Rational(2, 35) * sp.eye(7) + A - 5 * B * B.T
    print("H norm", inner(H, H))
    print("B norm2", sp.trace(B * B.T))
    print("Bhat", Bhat)
    print("Bhat eigen", Bhat.eigenvals())
    print("M", M)
    print("M eigen", M.eigenvals())
    for eigenvalue, multiplicity, vectors in M.eigenvects():
        if eigenvalue < 0:
            normalized = []
            for vector in vectors:
                vector = vector / sp.sqrt((vector.T * vector)[0])
                V = sum((vector[i] * H3[i] for i in range(7)), sp.ImmutableDenseNDimArray.zeros(3, 3, 3))
                normalized.append(V)
            print("negative", eigenvalue, multiplicity)
            print("Gram", sp.Matrix([[spin6_gram(V, W) for W in normalized] for V in normalized]))
            if multiplicity == 2:
                V, W = normalized
                products = [
                    symmetric_product(V, V),
                    symmetric_product(W, W),
                    sp.sqrt(2) * symmetric_product(V, W),
                ]
                compression_gram = sp.Matrix(
                    [[stf6_inner(P, Q) for Q in products] for P in products]
                )
                print("compression Gram", compression_gram)
                print("compression Gram inverse", compression_gram.inv())


if __name__ == "__main__":
    report()
