# KKT-infused flag algebra search for the \(P_2\) kernel

This repository contains a computer-assisted search for a sum-of-squares
(SOS) certificate for the copositivity of the even polynomial kernel

\[
K(t)=32t^6-48t^4+20t^2-\frac43,\qquad -1\le t\le 1.
\]

For a probability measure \(\mu\) on the real projective plane
\(\mathbb{RP}^2=S^2/\{x\sim -x\}\), define its \(K\)-energy by

\[
E_K(\mu)=
\iint K(x\cdot y)\,d\mu(x)\,d\mu(y).
\]

The desired copositivity statement is

\[
E_K(\mu)\ge 0
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

No exact certificate is claimed.

The strongest run currently uses degree \(14\), five-point root-weighted
flags, and root-factor degree \(2\). In the isotropic normalization

\[
T(\mu)=1-9p_4+6p_6=\frac{3}{16}E_K(\mu),
\qquad
p_j=\iint(x\cdot y)^j\,d\mu(x)d\mu(y),
\]

MOSEK reports

\[
T\ge -4.778814766126516\times10^{-8}.
\]

The maximum equality residual is \(2.66\times10^{-9}\), and the worst
moment-matrix eigenvalue is \(-9.89\times10^{-10}\). These numbers are strong
evidence for an exact zero bound, but floating-point output is not a proof.
The primal PSD matrices have not been recovered in a form that can be checked
over the rational numbers.

An exactification audit now reconstructs the ONB face over
\(\mathbb Q\), eliminates every unrestricted KKT/rank multiplier by rational
row reduction, and intersects with the exact pole–equator equality face. The
result has 503 independent coefficient equations. Its MOSEK primal iterates
drive feasibility errors down only by sending the trace and coefficient norms
to infinity. Thus the reported degree-14 solution cannot be rounded into a
finite rational certificate. See [Numerical results](RESULTS.md) for the
diagnostic data.

There is a second important limitation: the present implementation uses the
isotropy relation

\[
\int xx^T\,d\mu(x)=\frac13I.
\]

Therefore, exactifying the current certificate would prove the residual
**isotropic KKT case**. To deduce unrestricted copositivity, one must
additionally prove that a hypothetical negative global minimizer can be taken
isotropic, or extend the certificate algebra so that it does not assume
isotropy. This distinction is made precise in
[Certificate theorem](docs/CERTIFICATE_THEOREM.md).

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

Reproduce the strongest dual run:

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

Expected objective:

```text
-4.778814766126516e-08
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
