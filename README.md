# KKT-infused flag algebra search for the \(P_2\) kernel

This repository contains a computer-assisted search for a sum-of-squares
(SOS) certificate for the copositivity of the even polynomial kernel

\[
K(t)=32t^6-48t^4+20t^2-\frac43,\qquad -1\le t\le 1.
\]

For a probability measure \(\mu\) on the real projective plane
\(\mathbb{RP}^2=S^2/\{x\sim -x\}\), define its \(K\)-energy by

\[
E(\mu)=
\iint K(x\cdot y)\,d\mu(x)\,d\mu(y).
\]

The desired copositivity statement is

\[
E(\mu)\ge 0
\qquad
\text{for every probability measure }\mu\text{ on }\mathbb{RP}^2.
\]

The program searches for a proof by combining:

1. **flag-algebra squares:** positive semidefinite moment matrices obtained
   from conditional expectations of polynomial functions of sampled points;
2. **KKT conditions:** necessary first- and second-order conditions satisfied
   by a hypothetical energy-minimizing measure;
3. **geometric identities:** antipodal symmetry, isotropy, and the fact that
   every Gram matrix of four vectors in \(\mathbb R^3\) has determinant zero;
4. **semidefinite programming:** MOSEK chooses PSD Gram matrices and free
   equality multipliers whose coefficient expansion should equal the target
   energy.

This mixture is called **KKT-infused flag algebra** in this repository. It is
not a standard named theorem. It means that ordinary flag-square positivity is
augmented with positive multipliers of the KKT Hessian and global
first-variation inequalities.

## Current status

No exact certificate is claimed, and copositivity is **not** proved.

The implementation no longer assumes isotropy. The moment reducer keeps
every antipodally even expectation as an independent label — in particular
\(p_2=\iint(x\cdot y)^2\,d\mu\,d\mu\) is a genuine variable — and the target
is the energy itself,

\[
E(\mu)=-\frac43+20p_2-48p_4+32p_6,
\qquad
p_j=\iint(x\cdot y)^j\,d\mu(x)d\mu(y).
\]

The isotropy deficit is subsumed into the flag decomposition through
harmonic squares: the scalar block \(h_2=\tfrac{3p_2-1}{2}
=\sum_m|\hat\mu_{2m}|^2\ge0\) and a spin-2 Gram block containing the
deviatoric second moment \(D=\int xx^T d\mu-\tfrac13I\) paired against
multi-leaf spin-2 flags, so that \(h_2\to0\) forces the classical
contraction identities row-by-row.

With this general hierarchy, the strongest stable configuration (degree
\(14\), five-point root-weighted flags, matrix rank identities, spin-2
blocks) yields

\[
E\ \ge\ -4.77\times10^{-3},
\]

and every degree-14 variant lands in the band
\(-(4.5\text{–}5.0)\times10^{-3}\). Degree-16/18 runs exceed MOSEK's
numerical range for this formulation. The isotropic-branch near-certificate
\(E\gtrsim-2.6\times10^{-7}\) therefore does **not** survive the removal of
the isotropy assumption at the same degree and arity.

The moment diagnostics identify the mechanism: the optimal pseudo-moment
pays a small positive harmonic energy \(h_2\approx3.3\times10^{-3}\) and
buys contraction violations of order \(\sqrt{h_2}\), which any convex
spin-2 Gram relaxation permits by Cauchy–Schwarz. Real measures forbid this
because all spin-2 correlation vectors lie in a fixed five-dimensional
space (a non-convex rank constraint). Closing the gap at finite degree
would require certificate multipliers that annihilate every spin-2 residual
direction at the minimizing face — an open structural question. See
[Numerical results](RESULTS.md) for the full record.

## Start here

- [Mathematical background](docs/MATHEMATICAL_BACKGROUND.md) defines the
  energy, copositivity, KKT conditions, flags, SOS blocks, degree, and arity
  from the ground up.
- [Certificate theorem](docs/CERTIFICATE_THEOREM.md) gives the exact
  certificate formula and proves what a rationally verified certificate would
  imply.
- [Implementation guide](docs/IMPLEMENTATION.md) connects the mathematics to
  the Python code and command-line options.
- [Degree-14 command reference](docs/DEGREE14_COMMAND_REFERENCE.md) derives
  every flag, KKT kernel, rank identity, block size, and dual constraint used
  by the strongest five-point command.
- [Numerical results](RESULTS.md) records the hierarchy and solver outcomes.

## Installation

Python 3.11 or later is recommended.

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

MOSEK requires a separate license. See the
[MOSEK licensing documentation](https://www.mosek.com/products/academic-licenses/)
for academic-license options.

Run the tests:

```sh
python -m pytest -q
```

Reproduce the strongest (non-isotropic) dual run:

```sh
python sos_search.py \
  --dual --summary-only --scale-constraints --rank-relations \
  --higher-rank-matrices \
  --degree 14 --no-pointwise-sos \
  --harmonics --three-point-flags --four-point-flags --two-root-flags \
  --max-flag-arity 5 --max-root-factor-degree 2 \
  --gradient --potential --potential-matrices \
  --hessian --four-point-hessian \
  --global-gap --global-tangent-gaps \
  --tolerance 1e-9
```

Expected objective (≈16/3 × the legacy value; last digits solver-dependent):

```text
-4.77479077359651e-03
```

Small changes in solver version, scaling, or tolerances can change the last
digits. A numerical value near zero must not be reported as an exact SOS
certificate.

## Repository contents

```text
sos_search.py                 SDP construction and MOSEK interface
test_sos_search.py            exact-identity and equality-face tests
RESULTS.md                    numerical search record
docs/MATHEMATICAL_BACKGROUND.md
docs/CERTIFICATE_THEOREM.md
docs/IMPLEMENTATION.md
```

## Related work

- D. Bilyk and R. W. Matzke,
  [On the Fejes Tóth Problem about the Sum of Angles Between Lines](https://arxiv.org/abs/1801.07837).
- A. A. Razborov,
  [Flag Algebras](https://doi.org/10.2178/jsl/1203350785).
- H. Cohn and J. Woo,
  [Three-point bounds for energy minimization](https://arxiv.org/abs/1103.0485).

## License

No open-source license has been selected yet. Public visibility alone does not
grant permission to copy, modify, or redistribute the code. A license should
be added by the repository owner before inviting external reuse.
