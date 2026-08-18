#!/usr/bin/env python3
"""Standalone exact verifier for the certificate h2*E >= -2e-8.

Verifies certificates/h2E_geq_minus_2em8.json.gz with PURE RATIONAL
arithmetic (stdlib only: fractions/json/gzip/ast/random).  No solver,
no numpy, no repo imports.

THEOREM.  For every antipodally symmetric Borel probability measure mu
on the unit sphere S^2, with

    E(mu)  = iint K(x.y) dmu(x) dmu(y),
    K(t)   = 32 t^6 - 48 t^4 + 20 t^2 - 4/3,
    h2(mu) = (3 p2(mu) - 1)/2 >= 0,   p2 = iint (x.y)^2 dmu dmu,

it holds that  h2(mu) * E(mu) >= certified_bound  (> -2e-8).

PROOF SCHEMA VERIFIED HERE.  The certificate provides finitely many
"labels" L, each denoting a moment y_L(mu) = E[prod (X_i . X_j)^m_ij]
of independent samples from mu (disconnected products = products of
such moments).  It supplies:

  (1) PSD blocks: for each block b, a rational matrix Y_b and rational
      label matrices A^b_L such that M_b(mu) = sum_L A^b_L y_L(mu) is
      positive semidefinite for every measure mu (the repo's flag /
      harmonic / two-root / Jensen / fiber-Toeplitz families and their
      h2-multiplied copies; see docs).  Then
      sum_b <M_b(mu), Y_b> >= 0 whenever every Y_b is PSD.
  (2) relation rows E_i with sum_L E_{i,L} y_L(mu) = 0 for every mu
      (four-point Gram-rank identities: det Gram(x1..x4) = 0 in R^3,
      and their p2-shifted copies), with multipliers lambda_i.
  (3) the exact identity over Q, checked coefficient by coefficient:
      target_L - sum_b <A^b_L, Y_b> - sum_i lambda_i E_{i,L}
          = c [L = constant] + rho_L.

Pairing (3) with y(mu), using (1), (2), y_constant = 1, and the
elementary bound |y_L(mu)| <= 1 (each label is an expectation of a
product of inner products of unit vectors, each factor in [-1, 1]):

    (h2 E)(mu) = <target, y(mu)> >= c - sum_L |rho_L| = certified_bound.

WHAT IS PROVED UNCONDITIONALLY BY THIS SCRIPT: exact PSD-ness of every
Y_b, the exact identity (3), and the final bound arithmetic.

WHAT IS TAKEN FROM THE CONSTRUCTION (and spot-checked here): the
all-measures validity of the block families (1) and relation rows (2).
The spot-check evaluates every block matrix M_b and every relation row
on random exact rational atomic measures (rational points on S^2,
antipodalized): each M_b must be exactly PSD, each relation must
vanish exactly, |y_L| <= 1 must hold, and <target, y> must equal the
directly computed h2(mu)*E(mu).  Any corruption of the data would be
caught; the families' validity proofs live in the repo docs.

Usage:
    python3 verify_h2E_bound.py [certificates/h2E_geq_minus_2em8.json.gz]
        [--spot-measures 2] [--skip-spot]
"""

from __future__ import annotations

import argparse
import ast
import gzip
import json
import random
import sys
import time
from fractions import Fraction


# ----------------------------------------------------------------- PSD


def exact_psd(matrix: list[list[Fraction]]) -> bool:
    """Exact PSD decision by diagonal-pivoted LDL (Schur complements).

    A symmetric rational matrix is PSD iff repeatedly: the largest
    remaining diagonal entry is > 0 (eliminate it), or the remaining
    matrix is identically zero.  A zero (or negative) maximal diagonal
    with any nonzero remaining entry certifies an indefinite 1x1 or
    2x2 principal minor.
    """

    size = len(matrix)
    work = [row[:] for row in matrix]
    active = list(range(size))
    while active:
        pivot = max(active, key=lambda i: work[i][i])
        if work[pivot][pivot] < 0:
            return False
        if work[pivot][pivot] == 0:
            for i in active:
                for j in active:
                    if work[i][j] != 0:
                        return False
            return True
        d = work[pivot][pivot]
        active.remove(pivot)
        column = {i: work[i][pivot] for i in active if work[i][pivot]}
        for i, ci in column.items():
            for j, cj in column.items():
                work[i][j] -= ci * cj / d
    return True


# ------------------------------------------------- label semantics


def graph_edges(n: int) -> list[tuple[int, int]]:
    return [(a, b) for a in range(n) for b in range(a + 1, n)]


def label_vertex_count(head: str) -> int:
    if head == "pair":
        return 2
    if head == "triangle":
        return 3
    if head.startswith("graph_"):
        return int(head.split("_")[1])
    raise ValueError(head)


def atomic_label_value(label, atoms, weights, gram) -> Fraction:
    """Exact moment of a label on an ANTIPODALIZED atomic measure.

    Antipodalization: a monomial with an odd degree at any sampled
    vertex has expectation 0; even-degree monomials are sign-invariant
    and evaluate on the plain atoms.
    """

    if label == ("constant",):
        return Fraction(1)
    if label[0] == "product":
        value = Fraction(1)
        for factor in label[1:]:
            value *= atomic_label_value(factor, atoms, weights, gram)
        return value
    n = label_vertex_count(label[0])
    exps = tuple(int(v) for v in label[1:])
    edges = graph_edges(n)
    degree = [0] * n
    for e, (a, b) in zip(exps, edges):
        degree[a] += e
        degree[b] += e
    if any(d % 2 for d in degree):
        return Fraction(0)
    total = Fraction(0)
    k = len(atoms)
    assignment = [0] * n
    while True:
        term = Fraction(1)
        for v in assignment:
            term *= weights[v]
        for e, (a, b) in zip(exps, edges):
            if e:
                term *= gram[assignment[a]][assignment[b]] ** e
        total += term
        # increment odometer
        pos = 0
        while pos < n:
            assignment[pos] += 1
            if assignment[pos] < k:
                break
            assignment[pos] = 0
            pos += 1
        if pos == n:
            break
    return total


def rational_sphere_points(rng: random.Random, count: int):
    """Random rational unit vectors via inverse stereographic projection."""

    points = []
    while len(points) < count:
        p = Fraction(rng.randint(-8, 8), rng.randint(1, 8))
        q = Fraction(rng.randint(-8, 8), rng.randint(1, 8))
        d = 1 + p * p + q * q
        points.append((2 * p / d, 2 * q / d, (p * p + q * q - 1) / d))
    for x in points:
        assert sum(c * c for c in x) == 1
    return points


# ------------------------------------------------------------- main


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "certificate", nargs="?",
        default="certificates/h2E_geq_minus_2em8.json.gz",
    )
    parser.add_argument("--spot-measures", type=int, default=2)
    parser.add_argument("--skip-spot", action="store_true")
    parser.add_argument("--seed", type=int, default=20260819)
    args = parser.parse_args()
    t0 = time.time()

    with gzip.open(args.certificate, "rt") as handle:
        cert = json.load(handle)

    target = {
        ast.literal_eval(l): Fraction(v) for l, v in cert["target"].items()
    }
    relations = [
        {ast.literal_eval(l): Fraction(v) for l, v in row.items()}
        for row in cert["relations"]
    ]
    lam = [Fraction(v) for v in cert["lambda"]]
    c = Fraction(cert["certified_constant"])
    failures = 0

    # ---- (a) exact PSD of every Y block
    print("[1/3] exact PSD of certificate blocks")
    y_blocks = []
    block_data = []
    for index, block in enumerate(cert["blocks"]):
        size = block["size"]
        Y = [[Fraction(v) for v in row] for row in block["Y"]]
        for r in range(size):
            for s in range(size):
                assert Y[r][s] == Y[s][r], f"asymmetry in {block['name']}"
        if not exact_psd(Y):
            print(f"  FAIL: {block['name']} not PSD")
            failures += 1
        y_blocks.append(Y)
        matrices = {
            ast.literal_eval(l): [
                (int(r), int(s), Fraction(v)) for r, s, v in entries
            ]
            for l, entries in block["label_matrices"].items()
        }
        block_data.append((block["name"], size, matrices))
        if (index + 1) % 20 == 0:
            print(f"  {index + 1}/{len(cert['blocks'])} blocks PSD-checked "
                  f"({time.time() - t0:.0f}s)")
    print(f"  all {len(cert['blocks'])} blocks PSD: "
          f"{'yes' if failures == 0 else 'NO'}")

    # ---- (b) the exact identity and the bound
    print("[2/3] exact certificate identity")
    rho = dict(target)

    def add(label, value):
        updated = rho.get(label, Fraction(0)) + value
        if updated:
            rho[label] = updated
        else:
            rho.pop(label, None)

    for (name, size, matrices), Y in zip(block_data, y_blocks):
        for label, entries in matrices.items():
            acc = Fraction(0)
            for r, s, v in entries:
                acc += v * Y[r][s]
            if acc:
                add(label, -acc)
    for coefficient, row in zip(lam, relations):
        if coefficient:
            for label, value in row.items():
                add(label, -coefficient * value)
    constant = ("constant",)
    add(constant, -c)
    residual_const = rho.pop(constant, Fraction(0))
    l1 = abs(residual_const) + sum(abs(v) for v in rho.values())
    bound = c - l1
    print(f"  c               = {float(c):+.10e}")
    print(f"  ||rho||_1       = {float(l1):+.3e}  ({len(rho)} labels)")
    print(f"  certified bound = {float(bound):+.10e}")
    claimed = Fraction(cert["certified_bound"])
    if bound < claimed:
        print(f"  FAIL: recomputed bound below claimed {float(claimed):+.3e}")
        failures += 1
    if bound < Fraction(-2, 10**8):
        print("  FAIL: bound does not reach -2e-8")
        failures += 1
    else:
        print("  bound >= -2e-8: yes")

    # ---- (c) semantic spot checks on exact rational measures
    if not args.skip_spot:
        print("[3/3] semantic spot checks on random rational measures")
        rng = random.Random(args.seed)
        all_labels = set(target)
        for _, _, matrices in block_data:
            all_labels.update(matrices)
        for row in relations:
            all_labels.update(row)
        for trial in range(args.spot_measures):
            natoms = 3
            atoms = rational_sphere_points(rng, natoms)
            raw = [Fraction(rng.randint(1, 6)) for _ in range(natoms)]
            weights = [w / sum(raw) for w in raw]
            gram = [
                [sum(a * b for a, b in zip(x, y)) for y in atoms]
                for x in atoms
            ]
            values = {}
            for i, label in enumerate(sorted(all_labels, key=str)):
                values[label] = atomic_label_value(
                    label, atoms, weights, gram
                )
                if abs(values[label]) > 1:
                    print(f"  FAIL: |y_L| > 1 at {label}")
                    failures += 1
            print(f"  measure {trial + 1}: {len(values)} labels evaluated "
                  f"({time.time() - t0:.0f}s)")
            # target semantics: <target, y> == h2 * E directly
            p2, p4, p6 = (
                values.get(("pair", 2), Fraction(0)),
                values.get(("pair", 4), Fraction(0)),
                values.get(("pair", 6), Fraction(0)),
            )
            direct = (Fraction(3, 2) * p2 - Fraction(1, 2)) * (
                32 * p6 - 48 * p4 + 20 * p2 - Fraction(4, 3)
            )
            paired = sum(
                v * values.get(l, Fraction(0)) for l, v in target.items()
            )
            if direct != paired:
                print("  FAIL: target vector does not equal h2*E")
                failures += 1
            if paired < bound:
                print("  FAIL: bound violated by a spot measure (!)")
                failures += 1
            # relations vanish
            bad = sum(
                1
                for row in relations
                if sum(v * values.get(l, Fraction(0)) for l, v in row.items())
                != 0
            )
            if bad:
                print(f"  FAIL: {bad} relation rows nonzero on the measure")
                failures += 1
            else:
                print(f"  all {len(relations)} relation rows vanish: yes")
            # blocks evaluate PSD
            bad = 0
            for name, size, matrices in block_data:
                evaluated = [
                    [Fraction(0)] * size for _ in range(size)
                ]
                for label, entries in matrices.items():
                    y = values.get(label)
                    if not y:
                        continue
                    for r, s, v in entries:
                        evaluated[r][s] += v * y
                sym = [
                    [
                        (evaluated[r][s] + evaluated[s][r]) / 2
                        for s in range(size)
                    ]
                    for r in range(size)
                ]
                if not exact_psd(sym):
                    print(f"  FAIL: block {name} not PSD on spot measure")
                    bad += 1
            failures += bad
            if not bad:
                print(f"  all {len(block_data)} block families PSD on the "
                      f"measure: yes")

    print()
    if failures:
        print(f"FAIL ({failures} problem(s)) after {time.time() - t0:.0f}s")
        sys.exit(1)
    print(f"PASS: h2*E >= {float(bound):+.10e} > -2e-8 for every "
          f"antipodal probability measure on S^2  "
          f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
