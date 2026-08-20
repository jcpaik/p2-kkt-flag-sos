"""Small flag-SOS experiment for the conjectural inequality E/16 >= F."""

from fractions import Fraction

import cvxpy as cp
import numpy as np
import sympy as sp

from sos_search import (
    GRAM_DETERMINANT,
    ONE,
    PRINCIPAL_MINORS,
    expectation_label,
    flag_expectation_matrix,
    make_psd_block,
    monomials,
    tangent_harmonic_polynomials,
)


def add_triangle(target, exponents, coefficient):
    # Code convention is (XY,YZ,ZX); mathematical convention below is
    # (XY,XZ,YZ), hence swap the last two entries.
    r, s, t = exponents
    label, factor = expectation_label((r, t, s))
    assert label is not None
    target[label] = target.get(label, Fraction(0)) + coefficient * factor


target = {
    ("constant",): Fraction(-1, 12),
    ("pair", 2): Fraction(5, 4),
    ("pair", 4): Fraction(-3),
    ("pair", 6): Fraction(2),
}

# Subtract F from E/16.
f_terms = {
    (0, 4, 6): -16,
    (1, 5, 5): 16,
    (1, 3, 5): 32,
    (2, 4, 4): -40,
    (0, 4, 4): 8,
    (2, 2, 4): -24,
    (3, 3, 3): 32,
    (1, 3, 3): -16,
    (2, 2, 2): 8,
}
for exponents, coefficient in f_terms.items():
    add_triangle(target, exponents, Fraction(-coefficient))
target = {label: value for label, value in target.items() if value}


degree = 14
blocks = []
pointwise_terms = []
for name, multiplier, multiplier_degree in pointwise_terms or [
    # Keep the expression syntactically iterable while disabling the much
    # larger generic pointwise module in this higher-degree experiment.
]:
    basis = monomials((degree - multiplier_degree) // 2)
    Q, matrices = make_psd_block(name, basis, multiplier)
    blocks.append((name, Q, matrices))
for order, harmonic in enumerate(tangent_harmonic_polynomials(degree // 2)):
    max_leaf = (degree - 2 * order) // 2
    leaves = list(range(order % 2, max_leaf + 1, 2))
    if leaves:
        blocks.append(
            (
                f"flag_{order}",
                cp.Variable((len(leaves), len(leaves)), symmetric=True),
                flag_expectation_matrix(leaves, harmonic),
            )
        )

# Ordinary harmonic intensities I_l = E P_l(X.Y).
harmonic_blocks = []
t = sp.symbols("t")
for ell in range(2, degree + 1, 2):
    vector = {}
    for (power,), coefficient in sp.Poly(sp.legendre(ell, t), t).terms():
        label = ("constant",) if power == 0 else ("pair", int(power))
        vector[label] = Fraction(int(coefficient.p), int(coefficient.q))
    harmonic_blocks.append((ell, cp.Variable(nonneg=True), vector))

labels = set(target)
for _, _, matrices in blocks:
    labels.update(matrices)
for _, _, vector in harmonic_blocks:
    labels.update(vector)

constraints = [Q >> 0 for _, Q, _ in blocks]
for label in labels:
    rhs = 0
    for _, Q, matrices in blocks:
        if label in matrices:
            rhs += cp.sum(cp.multiply(Q, matrices[label]))
    for _, scalar, vector in harmonic_blocks:
        if label in vector:
            rhs += float(vector[label]) * scalar
    constraints.append(rhs == float(target.get(label, 0)))

problem = cp.Problem(
    cp.Minimize(
        sum(cp.trace(Q) for _, Q, _ in blocks)
        + sum(scalar for _, scalar, _ in harmonic_blocks)
    ),
    constraints,
)
for solver in (cp.CLARABEL, cp.SCS):
    try:
        value = problem.solve(solver=solver, verbose=False)
        print(solver, problem.status, value)
    except cp.error.SolverError as error:
        print(solver, error)

for name, Q, _ in blocks:
    if Q.value is not None:
        print(name, np.linalg.eigvalsh(Q.value))
for ell, scalar, _ in harmonic_blocks:
    if scalar.value is not None:
        print("harmonic", ell, scalar.value)
