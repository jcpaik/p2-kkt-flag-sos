# Implementation guide

## Overview

`sos_search.py` constructs either side of a finite semidefinite relaxation:

- the **primal** searches for SOS Gram matrices and KKT multipliers whose
  expansion equals the target \(T\);
- the **dual** searches for a formal moment functional satisfying every
  encoded positivity/equality constraint while minimizing \(T\).

The dual is numerically much more stable in the current singular problem.

## Moment labels

A Gram monomial is stored as a multigraph. Edge multiplicity is the exponent
of the corresponding inner product. `reduce_graph_matrix` applies:

- antipodal parity;
- isotropic contraction of degree-two vertices;
- factorization over disconnected components;
- canonical relabeling under vertex permutations.

Unreduced connected moments become labels such as `triangle`, `graph_4`,
`graph_5`, or `graph_6`.

## Positive blocks

The main block families are:

- `flag_*`: one-root \(O(2)\) harmonic flag squares;
- `two_root_*`: two-root even and orientation-odd sectors;
- `star_flag_*`: higher-arity flags with root-to-leaf factors;
- `weighted_flag_*`: higher-arity flags with additional root-root factors;
- `harmonic_*`: ordinary nonnegative spherical-harmonic moments;
- `hessian_*`: scalar and matrix-valued KKT Hessian multipliers;
- `global_*_gap`: nonnegative multipliers of the global first-variation gap.

Every primal variable in these families is constrained to be PSD.

## Free blocks and relations

The following are equality multipliers and therefore need not be PSD:

- `potential_flag_*`: \(U_\mu(x)-E(\mu)=0\) on the support;
- gradient relations: \(\nabla U_\mu(x)=0\);
- `rank_flag_*` and scalar rank relations:
  \(\det\operatorname{Gram}_4=0\);
- exact isotropic and exchangeability reductions built into the labels.

## Degree and arity

`--degree D` bounds the total Gram-polynomial degree.

`--max-flag-arity A` adds systematic conditional squares through \(A\)
sampled vertices after gluing. Five-point constraints begin at arity 5;
six-point constraints begin at arity 6.

`--max-root-factor-degree R` permits total root-root degree at most \(R\) in
the higher-arity flags. The current strongest run uses \(A=5\) and \(R=2\).

`--max-hessian-arity A` adds matrix-valued Hessian multipliers with additional
shared conditioning roots.

## Equality-face audit

`--check-onb` evaluates every PSD block and equality relation on the uniform
orthonormal-basis measure. The target must be zero, PSD blocks must have no
negative eigenvalues, and equality blocks must vanish.

`--facial-reduce-onb` parametrizes each primal PSD matrix directly on the
nullspace forced by equality at the ONB. This is valid only for the sharp
target with `--target-epsilon 0`.

## Constraint scaling

`--scale-constraints` rescales mathematically equivalent constraints before
calling MOSEK. This is essential near the singular equality face. Since
different scalings produced substantially different floating-point answers
in early runs, no objective is accepted without residual and nesting checks.

## Reproducing the best run

```sh
python sos_search.py \
  --dual --summary-only --scale-constraints --rank-relations \
  --degree 14 --no-pointwise-sos \
  --harmonics --three-point-flags --four-point-flags --two-root-flags \
  --max-flag-arity 5 --max-root-factor-degree 2 \
  --gradient --potential --potential-matrices \
  --hessian --four-point-hessian \
  --global-gap --global-tangent-gaps \
  --tolerance 1e-10
```

The JSON result includes:

- `objective`;
- `minimum_block_eigenvalue`;
- `maximum_free_residual`;
- `maximum_relation_residual`;
- the number of formal moment labels.

## Verification policy

A solver status of `optimal` means only that MOSEK satisfied its numerical
termination criteria. It is not a proof.

An exact certificate should be distributed as:

1. rational Gram matrices and free multipliers;
2. a machine-readable basis/label manifest;
3. an independent exact-arithmetic verifier;
4. a short human-readable derivation of every positive and zero term.

No such certificate file is currently present.
