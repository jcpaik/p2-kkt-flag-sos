#!/usr/bin/env python
"""Cylindrical mode-domination route for the P2 kernel.

Mission: E(mu) = int int K(x.y) dmu dmu >= 0 for all (antipodal) probability
measures mu on S^2, with

    K(t) = 32 t^6 - 48 t^4 + 20 t^2 - 4/3 = cos(6 theta) + cos(2 theta) + 2/3.

This module implements a route independent of the flag-algebra hierarchy:
condition on the measure's own principal axis e, write x.y = ab + p cos(dphi)
with a = x.e, p^2 = (1-a^2)(1-b^2), expand

    K(x.y) = sum_{k=0}^{6} c_k(a,b) cos(k (phi_x - phi_y)),

and dominate the azimuthal Fourier modes k = 1..4 by explicit penalties in the
height marginal sigma = mu_0 (modes 5, 6 are PSD and dropped).  The master
claim is B(sigma) >= 0 for every even probability sigma on [-1,1] with
int a^2 dsigma >= 1/3, where

    B(sigma) = int int c_0 dsigma dsigma - sum_{k=1}^{4} pen_k(sigma).

Every inequality used is a theorem; see docs/SUBCASES_AND_RECORD.md for
statements, proofs, and the validity ledger.  All mode data below is exact
(rational); `verify_mode_decomposition()` re-proves the decomposition as a
polynomial identity in pure sympy arithmetic, and `derive_modes_by_integration`
re-derives it from scratch by symbolic Fourier integrals (slow, ~2 min).

Sub-programs
------------
* M1 (axisymmetric sub-case): E restricted to axisymmetric mu equals
  int int C(s,t) dtau dtau with s = a^2 and tau an arbitrary probability
  measure on [0,1]; `m1_build_and_solve` searches for an exact rational
  certificate C = PSD-Gram + Handelman with facial (sharpness) constraints,
  `m1_rationalize` rounds it, and `verify_m1_certificate` re-verifies it in
  pure rational arithmetic.
* M2 (falsification): `master_bound` evaluates B(sigma) for discrete even
  sigma; `m2_hunt` searches for violations B < 0; `true_E_min_given_profile`
  computes the actual minimum of E over atomic measures with the given
  |a|-profile, to separate "scheme leak" from "conjecture threat".

CLI:  python cylinder_cert.py verify | derive | m1 [--solve] [--verify PATH]
      | m2 [--quick] | family
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
from fractions import Fraction

import sympy as sp

# ----------------------------------------------------------------------------
# 0.  Symbols and the kernel
# ----------------------------------------------------------------------------

a, b, z, s, t, w = sp.symbols("a b z s t w", real=True)

#: K(t) = 32 t^6 - 48 t^4 + 20 t^2 - 4/3  (variable named z here to avoid
#: clashing with the latitude symbols).
K_POLY = 32 * z**6 - 48 * z**4 + 20 * z**2 - sp.Rational(4, 3)

HALF = sp.Rational(1, 2)

# ----------------------------------------------------------------------------
# 1.  Exact mode decomposition (hard-coded, re-verified by
#     verify_mode_decomposition below).
#
# For k = 0..6 the mode kernel is
#     c_k(a,b) = (ab)^eps * ((1-a^2)(1-b^2))^{k/2} * sum_{ij} M_k[i][j] s^i t^j
# with eps = k mod 2, s = a^2, t = b^2.  Equivalently
#     c_k(a,b) = sum_{ij} M_k[i][j] f^k_i(a) f^k_j(b),
#     f^k_i(a) = a^{eps + 2 i} (1 - a^2)^{k/2},   i = 0..d_k - 1.
# ----------------------------------------------------------------------------

MODE_DIMS = [4, 3, 3, 2, 2, 1, 1]

MODE_MATRICES = {
    0: sp.Matrix(
        [
            [sp.Rational(2, 3), -4, 12, -10],
            [-4, 84, -270, 210],
            [12, -270, 840, -630],
            [-10, 210, -630, 462],
        ]
    ),
    1: sp.Matrix([[16, -96, 120], [-96, 624, -720], [120, -720, 792]]),
    2: sp.Matrix([[1, -6, 15], [-6, 132, -270], [15, -270, 495]]),
    3: sp.Matrix([[12, -60], [-60, 220]]),
    4: sp.Matrix([[0, -6], [-6, 66]]),
    5: sp.Matrix([[12]]),
    6: sp.Matrix([[1]]),
}


def mode_basis_function(k: int, i: int, var=a):
    """f^k_i(a) = a^{(k mod 2) + 2 i} (1-a^2)^{k/2} (exact sympy expr)."""
    eps = k % 2
    return var ** (eps + 2 * i) * (1 - var**2) ** sp.Rational(k, 2)


def mode_kernel(k: int, va=a, vb=b):
    """c_k(a,b) as an exact sympy expression."""
    M = MODE_MATRICES[k]
    d = MODE_DIMS[k]
    return sp.expand(
        sum(
            M[i, j] * mode_basis_function(k, i, va) * mode_basis_function(k, j, vb)
            for i in range(d)
            for j in range(d)
        )
    )


def c0_st_poly():
    """C(s,t): the mode-0 kernel in the squared variables s = a^2, t = b^2."""
    M = MODE_MATRICES[0]
    return sp.expand(sum(M[i, j] * s**i * t**j for i in range(4) for j in range(4)))


def verify_mode_decomposition() -> list[tuple[str, bool]]:
    """Fast exact proof of  K(ab + p z) = sum_k c_k(a,b) T_k(z)  in
    Q[a,b,z][p]/(p^2 - (1-a^2)(1-b^2)),  plus structural facts.

    Since T_k(cos psi) = cos(k psi), this polynomial identity is equivalent to
    the Fourier-mode decomposition, by uniqueness of Fourier coefficients.
    Returns a list of (check-name, bool).
    """
    checks: list[tuple[str, bool]] = []
    P = sp.expand((1 - a**2) * (1 - b**2))  # p^2

    # LHS: K(ab + p z) reduced to even + odd*p parts.
    even = sp.Integer(0)
    odd = sp.Integer(0)
    Kc = sp.Poly(K_POLY, z).all_coeffs()[::-1]  # coeff of z^0..z^6
    for n, coeff in enumerate(Kc):
        if coeff == 0:
            continue
        for j in range(n + 1):
            term = coeff * sp.binomial(n, j) * (a * b) ** (n - j) * z**j
            if j % 2 == 0:
                even += term * P ** (j // 2)
            else:
                odd += term * P ** ((j - 1) // 2)
    even, odd = sp.expand(even), sp.expand(odd)

    # RHS: sum_k c_k T_k(z), with p^k split into even/odd parts the same way.
    reven = sp.Integer(0)
    rodd = sp.Integer(0)
    for k in range(7):
        M = MODE_MATRICES[k]
        d = MODE_DIMS[k]
        poly_st = sum(
            M[i, j] * a ** (2 * i) * b ** (2 * j) for i in range(d) for j in range(d)
        )
        eps = k % 2
        core = (a * b) ** eps * poly_st * sp.chebyshevt(k, z)
        if k % 2 == 0:
            reven += core * P ** (k // 2)
        else:
            rodd += core * P ** ((k - 1) // 2)
    reven, rodd = sp.expand(reven), sp.expand(rodd)

    checks.append(("mode decomposition, even part", sp.expand(even - reven) == 0))
    checks.append(("mode decomposition, odd part", sp.expand(odd - rodd) == 0))

    # Parity: f^k_i(-a) = (-1)^k f^k_i(a), hence c_k(-a,b) = (-1)^k c_k(a,b).
    parity_ok = True
    for k in range(7):
        for i in range(MODE_DIMS[k]):
            f = mode_basis_function(k, i)
            if sp.simplify(f.subs(a, -a) - (-1) ** k * f) != 0:
                parity_ok = False
    checks.append(("mode parity f^k_i(-a) = (-1)^k f^k_i(a)", parity_ok))

    # Determinant fingerprints from PLAN.md section 2 (basis-dependent scalings
    # recorded there: c0 monomial basis, c2 as here, c3/c4 in the (a^2,1) form).
    # PLAN.md records dets in integer-normalized bases: (3/2) c0, c2, c3/4, c4/6.
    checks.append(
        ("det (3/2) c0 basis matrix = -34992", (sp.Rational(3, 2) * MODE_MATRICES[0]).det() == -34992)
    )
    checks.append(("det c2 basis matrix = -6480", MODE_MATRICES[2].det() == -6480))
    checks.append(("det c3/4 basis matrix = -60", (MODE_MATRICES[3] / 4).det() == -60))
    checks.append(("det c4/6 basis matrix = -1", (MODE_MATRICES[4] / 6).det() == -1))

    # PSD modes: c5 = 12 u (x) u, c6 = g (x) g -- rank one, positive coefficient.
    checks.append(("c5 rank-one positive", MODE_MATRICES[5] == sp.Matrix([[12]])))
    checks.append(("c6 rank-one positive", MODE_MATRICES[6] == sp.Matrix([[1]])))

    # Anchors: C(0,0) = 2/3, C(1,t) = K(sqrt(t)), C(1,1) = 8/3.
    C = c0_st_poly()
    checks.append(("c0(0,0) = 2/3", C.subs({s: 0, t: 0}) == sp.Rational(2, 3)))
    checks.append(
        (
            "c0(1,t) = K(sqrt(t))",
            sp.expand(C.subs(s, 1) - K_POLY.subs(z, sp.sqrt(t))) == 0,
        )
    )
    checks.append(("c0(1,1) = 8/3", C.subs({s: 1, t: 1}) == sp.Rational(8, 3)))

    # Pole-equator family: int int c0 d sigma_w^2 = 6 (w - 1/3)^2 with
    # sigma_w = (1-w) delta_0 + w (delta_1 + delta_-1)/2  (s-coords: (1-w) at 0, w at 1).
    fam = (
        w**2 * C.subs({s: 1, t: 1})
        + 2 * w * (1 - w) * C.subs({s: 1, t: 0})
        + (1 - w) ** 2 * C.subs({s: 0, t: 0})
    )
    checks.append(
        ("family energy = 6 (w - 1/3)^2", sp.expand(fam - 6 * (w - sp.Rational(1, 3)) ** 2) == 0)
    )

    # First-order stationarity at sigma*: F(t) = (2/3) C(0,t) + (1/3) C(1,t)
    # = 4 t (1-t)^2 >= 0 on [0,1], vanishing exactly on supp sigma*.
    F = sp.expand(sp.Rational(2, 3) * C.subs(s, 0) + sp.Rational(1, 3) * C.subs(s, 1))
    checks.append(("stationarity F(t) = 4 t (1-t)^2", sp.expand(F - 4 * t * (1 - t) ** 2) == 0))

    # Mode-4 split identity: c4 + 6 (a^2 + b^2) p^4 = 66 a^2 b^2 p^4.
    p4 = ((1 - a**2) * (1 - b**2)) ** 2
    checks.append(
        (
            "mode-4 split identity",
            sp.expand(mode_kernel(4) + 6 * (a**2 + b**2) * p4 - 66 * a**2 * b**2 * p4) == 0,
        )
    )

    # Mode-3 split: c3 = 220 u3(x)u3 + 12 u1(x)u1 - 60 (u3(x)u1 + u1(x)u3),
    # u3 = a^3 (1-a^2)^{3/2}, u1 = a (1-a^2)^{3/2}.  (Same content as M_3.)
    u3a, u1a = mode_basis_function(3, 1), mode_basis_function(3, 0)
    u3b, u1b = mode_basis_function(3, 1, b), mode_basis_function(3, 0, b)
    split3 = 220 * u3a * u3b + 12 * u1a * u1b - 60 * (u3a * u1b + u1a * u3b)
    checks.append(("mode-3 split identity", sp.expand(mode_kernel(3) - split3) == 0))

    return checks


def derive_modes_by_integration():
    """Slow independent derivation of all c_k by symbolic Fourier integration.

    Returns {k: c_k expression}.  Used by `python cylinder_cert.py derive`;
    takes ~2 minutes.
    """
    psi = sp.symbols("psi", real=True)
    p = sp.sqrt((1 - a**2) * (1 - b**2))
    expr = sp.expand_trig(sp.expand(K_POLY.subs(z, a * b + p * sp.cos(psi))))
    out = {}
    for k in range(8):
        weight = HALF if k == 0 else 1
        ck = weight / sp.pi * sp.integrate(expr * sp.cos(k * psi), (psi, 0, 2 * sp.pi))
        out[k] = sp.simplify(sp.expand(ck))
    return out


# ----------------------------------------------------------------------------
# 2.  Exact box-constrained QP:  min u^T Q u  over  0 <= u_i <= m_i.
#
# Homogeneous quadratic on a box; the minimum is attained at a KKT point.  We
# enumerate the 3^n activity patterns (coordinate at 0, at its upper bound, or
# free with stationarity), solve the free block exactly, and keep the best
# feasible candidate.  Singular free blocks are skipped: on such a face the
# minimum is also attained on the face's boundary, which other patterns cover.
# ----------------------------------------------------------------------------


def box_qp_min(Q, m, exact: bool = True):
    """Return (min value, argmin list) of u^T Q u over the box prod [0, m_i].

    Q: sympy Matrix (symmetric).  m: list of nonnegative sympy numbers (exact)
    or floats (exact=False).  Always <= 0 since u = 0 is feasible.
    """
    n = len(m)
    if exact:
        Qm = sp.Matrix(Q)
        zero = sp.Integer(0)
        mvec = [sp.nsimplify(mi) if not isinstance(mi, sp.Expr) else mi for mi in m]
    else:
        import numpy as np

        Qm = np.array([[float(Q[i, j]) for j in range(n)] for i in range(n)])
        zero = 0.0
        mvec = [float(mi) for mi in m]

    best_val, best_u = None, None
    for pattern in itertools.product((0, 1, 2), repeat=n):
        free = [i for i in range(n) if pattern[i] == 2]
        u = [zero if pattern[i] == 0 else mvec[i] for i in range(n)]
        if free:
            fixed = [i for i in range(n) if pattern[i] != 2]
            if exact:
                A = Qm[free, free]
                if A.det() == 0:
                    continue
                rhs = -sp.Matrix([
                    sum(Qm[i, jj] * u[jj] for jj in fixed) for i in free
                ])
                sol = A.LUsolve(rhs)
                ok = True
                for idx, i in enumerate(free):
                    val = sp.nsimplify(sol[idx])
                    lo = sp.simplify(val)
                    if (lo < 0) == True or (lo > mvec[i]) == True:  # noqa: E712
                        ok = False
                        break
                    u[i] = lo
                if not ok:
                    continue
            else:
                import numpy as np

                A = Qm[np.ix_(free, free)]
                if abs(np.linalg.det(A)) < 1e-14:
                    continue
                rhs = -np.array([sum(Qm[i, jj] * u[jj] for jj in fixed) for i in free])
                sol = np.linalg.solve(A, rhs)
                if any(v < -1e-12 or v > mvec[i] + 1e-12 for v, i in zip(sol, free)):
                    continue
                for idx, i in enumerate(free):
                    u[i] = float(sol[idx])
        # value
        if exact:
            uv = sp.Matrix(u)
            val = sp.expand((uv.T * Qm * uv)[0, 0])
        else:
            val = sum(u[i] * Qm[i, j] * u[j] for i in range(n) for j in range(n))
        if best_val is None:
            best_val, best_u = val, list(u)
        else:
            if exact:
                if sp.simplify(val - best_val) < 0:
                    best_val, best_u = val, list(u)
            elif val < best_val:
                best_val, best_u = val, list(u)
    return best_val, best_u


def signed_abs(Q):
    """M-tilde: keep the diagonal, replace off-diagonal entries by -|entry|."""
    n = Q.shape[0]
    out = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            out[i, j] = Q[i, j] if i == j else -abs(Q[i, j])
    return out


# ----------------------------------------------------------------------------
# 3.  Penalties and the master bound B(sigma) for discrete even sigma.
#
# sigma is given as atoms [(alpha_j, w_j)] with alpha_j in [0,1]: each atom
# stands for the even pair (delta_{alpha} + delta_{-alpha})/2 (just delta_0
# when alpha = 0).  All moment functions below are even, so only alpha >= 0
# matters.  Validity of each penalty is proved in
# docs/SUBCASES_AND_RECORD.md (validity ledger).
# ----------------------------------------------------------------------------

#: Reduced mode-1 matrix: principal-axis constraint kills the first coordinate
#: (int a sqrt(1-a^2) dmu_1 = 0), leaving the (a^3, a^5) block of M_1.
MODE1_REDUCED = sp.Matrix([[624, -720], [-720, 792]])


def _moments(atoms, funcs, exact: bool):
    """sum_j w_j f(alpha_j) for each f in funcs (f given as sympy expr in a)."""
    out = []
    for f in funcs:
        if exact:
            out.append(sp.simplify(sum(wj * f.subs(a, al) for al, wj in atoms)))
        else:
            fl = sp.lambdify(a, f, "math")
            out.append(sum(float(wj) * fl(float(al)) for al, wj in atoms))
    return out


def penalty_mode(k: int, atoms, exact: bool = True, use_axis_constraint: bool = True):
    """pen_k(sigma) >= 0 with Q_k >= -pen_k(sigma): signed-abs box-QP penalty.

    Modes 5, 6 are PSD (penalty 0).  Mode 1 uses the principal-axis constraint
    z_1 = 0 when use_axis_constraint (the definition of B assumes e is the
    principal axis).
    """
    if k in (0, 5, 6):
        return sp.Integer(0) if exact else 0.0
    eps = k % 2
    d = MODE_DIMS[k]
    if k == 1 and use_axis_constraint:
        Q = MODE1_REDUCED
        idxs = [1, 2]
    else:
        Q = MODE_MATRICES[k]
        idxs = list(range(d))
    # dominated moments of |f^k_i|; all basis functions have |f| even in a and
    # equal to alpha^{eps+2i} (1-alpha^2)^{k/2} at |a| = alpha.
    funcs = [a ** (eps + 2 * i) * (1 - a**2) ** sp.Rational(k, 2) for i in idxs]
    mom = _moments(atoms, funcs, exact)
    val, _ = box_qp_min(signed_abs(Q), mom, exact=exact)
    return -val


def c0_energy(atoms, exact: bool = True):
    """int int c0 dsigma dsigma for the discrete even sigma."""
    C = c0_st_poly()
    if exact:
        vals = [(sp.nsimplify(al) ** 2, wj) for al, wj in atoms]
        return sp.simplify(
            sum(
                wi * wj * C.subs({s: si, t: sj})
                for si, wi in vals
                for sj, wj in vals
            )
        )
    Cf = sp.lambdify((s, t), C, "math")
    vals = [(float(al) ** 2, float(wj)) for al, wj in atoms]
    return sum(wi * wj * Cf(si, sj) for si, wi in vals for sj, wj in vals)


def master_bound(atoms, exact: bool = True, use_axis_constraint: bool = True):
    """B(sigma) and its parts for discrete even sigma = sum w_j pair(alpha_j).

    Returns dict with keys: c0, pen1..pen4, B, lambda1 (= int a^2 dsigma),
    feasible (lambda1 >= 1/3).
    """
    total = sum(wj for _, wj in atoms)
    if exact:
        assert sp.simplify(sp.nsimplify(total) - 1) == 0, "weights must sum to 1"
    else:
        assert abs(float(total) - 1) < 1e-9, "weights must sum to 1"
    out = {}
    out["c0"] = c0_energy(atoms, exact)
    pens = {}
    for k in (1, 2, 3, 4):
        pens[k] = penalty_mode(k, atoms, exact, use_axis_constraint)
        out[f"pen{k}"] = pens[k]
    if exact:
        out["B"] = sp.simplify(out["c0"] - sum(pens.values()))
        out["lambda1"] = sp.simplify(sum(wj * sp.nsimplify(al) ** 2 for al, wj in atoms))
        out["feasible"] = bool(sp.simplify(out["lambda1"] - sp.Rational(1, 3)) >= 0)
    else:
        out["B"] = out["c0"] - sum(pens.values())
        out["lambda1"] = sum(float(wj) * float(al) ** 2 for al, wj in atoms)
        out["feasible"] = out["lambda1"] >= 1 / 3 - 1e-12
    return out


# ----------------------------------------------------------------------------
# 4.  M1: the axisymmetric sub-case as 1-D copositivity of C(s,t) on [0,1].
#
# Certificate ansatz (all data rational):
#     C(s,t) = Phi(s)^T G Phi(t) + sum_beta lambda_beta B_beta(s,t),
# Phi = (1, s, s^2, s^3), G PSD with G v* = 0 for the zero-measure moment
# vector v* = (1, 1/3, 1/3, 1/3); B_beta are symmetrized Handelman products
# s^i (1-s)^j t^k (1-t)^l that vanish on the corner set {0,1}^2; lambda >= 0.
# Then int int C dtau dtau >= 0 for every probability tau on [0,1], with
# equality structure matching tau* = (2/3) delta_0 + (1/3) delta_1.
# ----------------------------------------------------------------------------

#: rational basis of the orthogonal complement of v* = (1, 1/3, 1/3, 1/3),
#: i.e. of (3,1,1,1)^perp, as columns.
FACE_P = sp.Matrix([[1, 0, 0], [-3, 1, 0], [0, -1, 1], [0, 0, -1]])


def handelman_pairs(deg: int):
    """Symmetrized corner-vanishing Handelman index pairs up to degree `deg`.

    Elements: frozenset-like ordered pair (u, v) with u = (i,j), v = (k,l),
    u <= v lexicographically, such that s^i(1-s)^j t^k(1-t)^l vanishes at all
    four corners of [0,1]^2, which (for a product of nonnegative factors)
    happens iff u or v has both exponents >= 1.
    """
    anyset = [(i, j) for i in range(deg + 1) for j in range(deg + 1) if 0 < i + j <= deg]
    anyset.append((0, 0))
    vset = [(i, j) for (i, j) in anyset if i >= 1 and j >= 1]
    pairs = set()
    for u in vset:
        for v_ in anyset:
            pairs.add(tuple(sorted((u, v_))))
    return sorted(pairs)


def handelman_poly(pair):
    """Symmetrized basis polynomial for an index pair ((i,j),(k,l))."""
    (i, j), (k, l) = pair
    f = s**i * (1 - s) ** j * t**k * (1 - t) ** l
    g = s**k * (1 - s) ** l * t**i * (1 - t) ** j
    return sp.expand(HALF * (f + g))


def m1_system(deg: int):
    """Exact linear system for the certificate: unknowns x = (vech W, lambda).

    Column order: the 6 entries of the 3x3 symmetric W (00,01,02,11,12,22),
    then one lambda per Handelman pair.  Row per monomial s^i t^j, i <= j <=
    D = max(3, deg).  Returns (A, c, pairs, monoms) with A, c exact rational
    sympy Matrices such that  A x = c  <=>  C = Phi^T (P W P^T) Phi + sum
    lambda_beta B_beta  as polynomials.
    """
    pairs = handelman_pairs(deg)
    D = max(3, deg)
    monoms = [(i, j) for i in range(D + 1) for j in range(i, D + 1)]
    midx = {m: k for k, m in enumerate(monoms)}

    def sym_coeffs(expr):
        v = [sp.Integer(0)] * len(monoms)
        p_ = sp.Poly(sp.expand(expr), s, t)
        for (i, j), cc in p_.terms():
            if i <= j:
                v[midx[(i, j)]] = cc
        return v

    cols = []
    Phi_s = sp.Matrix([1, s, s**2, s**3])
    Phi_t = sp.Matrix([1, t, t**2, t**3])
    for i in range(3):
        for j in range(i, 3):
            E = sp.zeros(3, 3)
            E[i, j] = 1
            E[j, i] = 1
            G = FACE_P * E * FACE_P.T
            cols.append(sym_coeffs((Phi_s.T * G * Phi_t)[0, 0]))
    for pr_ in pairs:
        cols.append(sym_coeffs(handelman_poly(pr_)))
    A = sp.Matrix(cols).T
    c_vec = sp.Matrix(sym_coeffs(c0_st_poly()))
    return A, c_vec, pairs, monoms


def m1_build_and_solve(deg: int = 8, verbose: bool = False, support_tol: float = 1e-6):
    """Numeric (double) solve of the M1 certificate SDP with cvxpy.

    Stage 1: maximize lambda_min(W) subject to the exact linear identity and
    lambda >= 0.  Stage 2: restrict to the numeric support of lambda and
    maximize a joint interior margin (W >= gamma I, lambda_S >= gamma), which
    is what rational rounding needs.  Returns dict with W, lam, support,
    margin, status.
    """
    import cvxpy as cp
    import numpy as np

    A, c_vec, pairs, monoms = m1_system(deg)
    Anp = np.array([[float(A[i, j]) for j in range(A.cols)] for i in range(A.rows)])
    cnp = np.array([float(x) for x in c_vec])

    def base_problem(support=None, extra=None, objective="margin"):
        W = cp.Variable((3, 3), symmetric=True)
        gamma = cp.Variable()
        cols = list(range(len(pairs))) if support is None else list(support)
        lam = cp.Variable(len(cols), nonneg=(objective != "joint"))
        Asub = np.hstack([Anp[:, :6], Anp[:, [6 + i for i in cols]]])
        x = cp.hstack([W[0, 0], W[0, 1], W[0, 2], W[1, 1], W[1, 2], W[2, 2], lam])
        cons = [Asub @ x == cnp, W - gamma * np.eye(3) >> 0]
        if objective == "margin":
            prob = cp.Problem(cp.Maximize(gamma), cons)
        elif objective == "sparsify":
            cons.append(gamma >= extra)
            prob = cp.Problem(cp.Minimize(cp.sum(lam)), cons)
        else:  # joint margin
            cons.append(lam >= gamma)
            prob = cp.Problem(cp.Maximize(gamma), cons)
        for solver in ("MOSEK", "CLARABEL", "SCS"):
            try:
                prob.solve(solver=solver, verbose=verbose)
            except Exception:
                continue
            if prob.status in ("optimal", "optimal_inaccurate"):
                break
        gv = None if gamma.value is None else float(gamma.value)
        return prob.status, W.value, lam.value, gv, cols

    # Stage 1: pure feasibility margin.
    status, Wv, lamv, margin, _ = base_problem(objective="margin")
    result = {"status": status, "deg": deg, "pairs": pairs, "stage1_margin": margin}
    if Wv is None or margin is None or margin <= 0:
        return result
    # Stage 1b: sparsify at half the margin.
    st_s, Wv_s, lam_s, _, _ = base_problem(extra=0.5 * margin, objective="sparsify")
    if lam_s is not None:
        support = [i for i, v in enumerate(lam_s) if v > support_tol]
    else:
        support = [i for i, v in enumerate(lamv) if v > support_tol]
    # Stage 2: joint interior margin on the support.
    st2, Wv2, lam2, margin2, cols = base_problem(support=support, objective="joint")
    if Wv2 is not None and margin2 is not None and margin2 > 0:
        result.update(
            {
                "status": st2,
                "W": Wv2.tolist(),
                "lam_support": {int(i): float(v) for i, v in zip(cols, lam2)},
                "margin": margin2,
            }
        )
    else:
        result.update(
            {
                "W": Wv.tolist(),
                "lam_support": {
                    int(i): float(v) for i, v in enumerate(lamv) if v > support_tol
                },
                "margin": margin,
            }
        )
    return result


def m1_rationalize(numeric: dict, denom_bound: int = 3600):
    """Round the numeric certificate onto the exact affine set {A x = c}.

    RREF scheme (keeps denominators small): compute the exact reduced row
    echelon form of the support-restricted system, round only the FREE
    variables to small rationals, and solve the pivot variables exactly.  The
    strict interior margin of the numeric point absorbs the perturbation,
    keeping W PSD and lambda >= 0; `verify_m1_certificate` re-checks all of it
    in pure rational arithmetic.
    """
    deg = numeric["deg"]
    A, c_vec, pairs, _ = m1_system(deg)
    support = sorted(int(i) for i in numeric["lam_support"])
    lam_val = dict(numeric["lam_support"])

    def build(support_):
        keep_ = list(range(6)) + [6 + i for i in support_]
        Ak_ = A[:, keep_]
        aug_ = Ak_.row_join(c_vec)
        R_, piv_ = aug_.rref()
        consistent = (Ak_.cols) not in piv_
        return Ak_, R_, [p_ for p_ in piv_ if p_ < Ak_.cols], consistent

    Ak, R, pivots, ok = build(support)
    if not ok:
        # restricted system exactly inconsistent: grow the support with the
        # largest remaining numeric lambdas until consistent.
        rest = sorted(
            (i for i in range(len(pairs)) if i not in support),
            key=lambda i: -lam_val.get(i, 0.0),
        )
        for i in rest:
            support = sorted(support + [i])
            lam_val.setdefault(i, 0.0)
            Ak, R, pivots, ok = build(support)
            if ok:
                break
        if not ok:
            raise RuntimeError("could not reach an exactly consistent support")
    Wn = numeric["W"]
    xnum = [Wn[0][0], Wn[0][1], Wn[0][2], Wn[1][1], Wn[1][2], Wn[2][2]]
    xnum += [lam_val[i] for i in support]
    free = [j for j in range(Ak.cols) if j not in pivots]
    x_rat = [None] * Ak.cols
    for j in free:
        x_rat[j] = sp.Rational(Fraction(float(xnum[j])).limit_denominator(denom_bound))
    # back-substitute pivots: row r reads  x_{piv r} + sum_f R[r,f] x_f = R[r,-1]
    for r, pj in enumerate(pivots):
        val = R[r, Ak.cols]
        for f in free:
            if R[r, f] != 0:
                val -= R[r, f] * x_rat[f]
        x_rat[pj] = sp.nsimplify(val)
    x_rat = sp.Matrix(x_rat)
    assert all((Ak * x_rat - c_vec)[i] == 0 for i in range(c_vec.rows))
    Wr = sp.Matrix(
        [
            [x_rat[0], x_rat[1], x_rat[2]],
            [x_rat[1], x_rat[3], x_rat[4]],
            [x_rat[2], x_rat[4], x_rat[5]],
        ]
    )
    lam_exact = {pairs[i]: x_rat[6 + jj] for jj, i in enumerate(support)}
    lam_exact = {pr_: v for pr_, v in lam_exact.items() if v != 0}
    return {"deg": deg, "W": Wr, "lam": lam_exact}


def verify_m1_certificate(cert: dict) -> list[tuple[str, bool]]:
    """Pure-rational verification of an exact M1 certificate.

    cert: {"deg": int, "W": sympy 3x3 rational Matrix, "lam": {pair: Rational}}
    Checks: (i) exact polynomial identity, (ii) W PSD, (iii) lam >= 0 and
    corner-vanishing basis, (iv) G v* = 0.  No floats anywhere.
    """
    checks = []
    W = sp.Matrix(cert["W"])
    lam = cert["lam"]
    G = FACE_P * W * FACE_P.T
    Phi_s = sp.Matrix([1, s, s**2, s**3])
    Phi_t = sp.Matrix([1, t, t**2, t**3])
    gram = sp.expand((Phi_s.T * G * Phi_t)[0, 0])
    hand = sp.Integer(0)
    corner_ok = True
    nonneg_ok = True
    for pair, coef in lam.items():
        (i, j), (k, l) = pair
        if not (i >= 1 and j >= 1) and not (k >= 1 and l >= 1):
            corner_ok = False
        if sp.Rational(coef) < 0:
            nonneg_ok = False
        hand += sp.Rational(coef) * handelman_poly(pair)
    resid = sp.expand(c0_st_poly() - gram - hand)
    checks.append(("polynomial identity C = Gram + Handelman", resid == 0))
    checks.append(("W is PSD (exact)", W.is_positive_semidefinite is True))
    checks.append(("all lambda >= 0", nonneg_ok))
    checks.append(("Handelman basis corner-vanishing", corner_ok))
    vstar = sp.Matrix([1, sp.Rational(1, 3), sp.Rational(1, 3), sp.Rational(1, 3)])
    checks.append(("G v* = 0 (sharp at sigma*)", (G * vstar).norm() == 0))
    # G PSD follows from W PSD; assert anyway.
    checks.append(("G is PSD (exact)", G.is_positive_semidefinite is True))
    return checks


def m1_cert_to_json(cert: dict) -> dict:
    return {
        "deg": cert["deg"],
        "W": [[str(sp.Rational(cert["W"][i, j])) for j in range(3)] for i in range(3)],
        "lam": {str(k): str(sp.Rational(v)) for k, v in cert["lam"].items()},
        "face_P": [[str(FACE_P[i, j]) for j in range(3)] for i in range(4)],
        "vstar": ["1", "1/3", "1/3", "1/3"],
        "meaning": "C(s,t) = Phi^T (P W P^T) Phi + sum lam * sym Handelman; "
        "Phi=(1,s,s^2,s^3); proves int int C dtau dtau >= 0 on prob(0,1)",
    }


def m1_cert_from_json(d: dict) -> dict:
    W = sp.Matrix(3, 3, lambda i, j: sp.Rational(d["W"][i][j]))
    lam = {}
    for kstr, vstr in d["lam"].items():
        pair = tuple(tuple(int(x) for x in half) for half in sp.sympify(kstr))
        lam[pair] = sp.Rational(vstr)
    return {"deg": d["deg"], "W": W, "lam": lam}


# ----------------------------------------------------------------------------
# 4b.  Exact single-orbit analysis (the M2 verdict's backbone).
#
# For the single-pair profile sigma = (delta_alpha + delta_-alpha)/2 and a
# measure mu with that height profile and upper-fiber distribution rho,
#     E(mu) = C(s,s) + sum_{k=1}^{6} |rho-hat(k)|^2 v_k(s),      s = alpha^2,
# where v_k(s) = f^k(alpha)^T M_k f^k(alpha) is an exact polynomial in s
# (all half-powers square away).  This identity is exact and makes the
# per-mode scalar-domination limit computable in closed form.
# ----------------------------------------------------------------------------


def single_orbit_mode_value(k: int, sv=s):
    """v_k(s) = s^eps (1-s)^k * sum_ij M_k[i,j] s^{i+j}, exact polynomial."""
    eps = k % 2
    M = MODE_MATRICES[k]
    d = MODE_DIMS[k]
    return sp.expand(
        sv**eps * (1 - sv) ** k * sum(M[i, j] * sv ** (i + j) for i in range(d) for j in range(d))
    )


def b_tight_single_orbit(sv=s):
    """The best bound any per-mode scalar-domination scheme can give on the
    single-pair profile (exact expression in s = alpha^2):

        B_tight(s) = C(s,s) + sum_{k=2}^{4} min(0, v_k(s)).

    Modes 5, 6 have v_k >= 0; mode 1 is killed by the principal-axis
    constraint z_1 = 0 (single orbit => zeta = 0).  A scheme knowing only
    |mu_k| <= sigma (+ parity + axis constraints) cannot do better, because
    the coherent measure mu_k = zeta * (fiber) with |zeta| = 1 is dominated.
    """
    C = c0_st_poly()
    val = C.subs({s: sv, t: sv})
    parts = {0: sp.expand(val)}
    for k in (2, 3, 4):
        parts[k] = single_orbit_mode_value(k, sv)
    return parts


def b_tight_single_orbit_value(s_val):
    """Exact value of B_tight at rational s_val (min(0, .) resolved exactly)."""
    s_val = sp.nsimplify(s_val)
    C = c0_st_poly().subs({s: s_val, t: s_val})
    tot = C
    for k in (2, 3, 4):
        vk = single_orbit_mode_value(k).subs(s, s_val)
        tot += sp.Min(0, vk)
    return sp.nsimplify(sp.simplify(tot))


def verify_circle_pair_theorem() -> list[tuple[str, bool]]:
    """Exact proof of the circle-pair theorem:

        E(mu) >= 0 for every antipodal measure mu supported on a pair of
        antipodal circles {a = +-alpha}, with arbitrary fiber measures;
        equality iff s = alpha^2 = 1/3 with |rho-hat(3)| = |rho-hat(6)| = 1
        (the ONB orbit) -- or the trivial zero-energy witnesses on the family.

    Proof scheme (docs/SUBCASES_AND_RECORD.md section 6): with
    E = C2(s) + sum_k y_k v_k(s), y_k = |rho-hat(k)|^2 in [0,1], and the
    Toeplitz couplings y_{2k} >= ((2 y_k - 1)^+)^2 for k = 2, 3, the minimum
    over admissible y is bounded below by explicit piecewise-rational
    formulas whose nonnegativity reduces to seven polynomial facts on
    rational intervals, checked here by exact Sturm root counts.
    """
    checks: list[tuple[str, bool]] = []
    C2 = sp.expand(c0_st_poly().subs(t, s))
    v = {k: single_orbit_mode_value(k) for k in range(1, 7)}
    g = 55 * s**2 - 30 * s + 3
    q1 = 99 * s**4 - 180 * s**3 + 108 * s**2 - 24 * s + 2
    q2 = 495 * s**4 - 540 * s**3 + 162 * s**2 - 12 * s + 1

    # structural identities
    checks.append(("v1 = 8 s (1-s) q1", sp.expand(v[1] - 8 * s * (1 - s) * q1) == 0))
    checks.append(("v2 = (1-s)^2 q2", sp.expand(v[2] - (1 - s) ** 2 * q2) == 0))
    checks.append(("v3 = 4 s (1-s)^3 g", sp.expand(v[3] - 4 * s * (1 - s) ** 3 * g) == 0))
    checks.append(
        ("v4 = 6 s (1-s)^4 (11 s - 2)", sp.expand(v[4] - 6 * s * (1 - s) ** 4 * (11 * s - 2)) == 0)
    )
    checks.append(("v5 = 12 s (1-s)^5 >= 0", sp.expand(v[5] - 12 * s * (1 - s) ** 5) == 0))
    checks.append(("v6 = (1-s)^6 >= 0", sp.expand(v[6] - (1 - s) ** 6) == 0))
    checks.append(
        (
            "(3,6) branch identity 4 v6 + v3 = 4 (1-s)^3 (3s-1)^2 (6s+1)",
            sp.expand(4 * v[6] + v[3] - 4 * (1 - s) ** 3 * (3 * s - 1) ** 2 * (6 * s + 1)) == 0,
        )
    )
    # I36 = C2 + v3/2 - s^2 g^2 with the sharp (3s-1)^2 factor
    h_poly = 891 * s**4 - 216 * s**3 - 81 * s**2 - 6 * s - 2
    checks.append(
        (
            "C2 + v3/2 - s^2 g^2 = -(3s-1)^2 h / 3",
            sp.expand(C2 + v[3] / 2 - s**2 * g**2 + (3 * s - 1) ** 2 * h_poly / 3) == 0,
        )
    )

    def nonneg_on(expr, lo, hi, name, strict_negative=False):
        """expr >= 0 on [lo, hi] (or <= 0 if strict_negative) by exact Sturm
        count of interior roots plus a sign sample; endpoint zeros allowed
        only when the polynomial factors them out explicitly (we avoid that
        case by choosing endpoints that are not roots)."""
        p_ = sp.Poly(sp.expand(expr), s)
        nroots = p_.count_roots(lo, hi)
        mid = (sp.Rational(lo) + sp.Rational(hi)) / 2
        sample = p_.eval(mid)
        if strict_negative:
            ok = nroots == 0 and sample < 0
        else:
            ok = nroots == 0 and sample > 0
        checks.append((name, bool(ok)))

    # window brackets (rational certificates)
    lo3, hi3 = sp.Rational(263, 2000), sp.Rational(42, 100)  # covers {g <= 0}
    checks.append(("g > 0 at 263/2000", g.subs(s, lo3) > 0))
    checks.append(("g > 0 at 42/100", g.subs(s, hi3) > 0))
    checks.append(("g < 0 at 33/250", g.subs(s, sp.Rational(33, 250)) < 0))
    checks.append(
        ("g-window inside (263/2000, 42/100)",
         sp.Poly(g, s).count_roots(0, lo3) == 0 and sp.Poly(g, s).count_roots(hi3, 1) == 0)
    )
    lo2, hi2 = sp.Rational(2, 5), sp.Rational(13, 20)  # covers {q2 <= 0}
    checks.append(("q2 > 0 at 2/5", q2.subs(s, lo2) > 0))
    checks.append(("q2 > 0 at 13/20", q2.subs(s, hi2) > 0))
    checks.append(
        ("q2-window inside (2/5, 13/20)",
         sp.Poly(q2, s).count_roots(0, lo2) == 0 and sp.Poly(q2, s).count_roots(hi2, 1) == 0)
    )
    checks.append(("q2 < 0 at 1/2", q2.subs(s, sp.Rational(1, 2)) < 0))

    # (P0) C2 > 0 on [0,1]
    nonneg_on(C2, 0, 1, "(P0) C2 > 0 on [0,1]")
    # (P1) q1 > 0 on [0,1]  => v1 >= 0
    nonneg_on(q1, 0, 1, "(P1) q1 > 0 on [0,1] (v1 >= 0)")
    # (P2b) beta2 = q2 + 24 s (1-s)^2 (11s-2) > 0 on the q2-window cover
    beta2 = sp.expand(q2 + 24 * s * (1 - s) ** 2 * (11 * s - 2))
    nonneg_on(beta2, lo2, hi2, "(P2b) (2,4) coupling in smooth branch on q2-window")
    # (PA1) C2 + v4 > 0 on [0, 33/250]
    nonneg_on(C2 + v[4], 0, sp.Rational(33, 250), "(PA1) C2 + v4 > 0 on [0, 33/250]")
    # (PA2) C2 + v4 + v3/2 - s^2 g^2 > 0 on [263/2000, 2/11]
    nonneg_on(
        C2 + v[4] + v[3] / 2 - s**2 * g**2,
        lo3,
        sp.Rational(2, 11),
        "(PA2) C2 + v4 + P36 > 0 on [263/2000, 2/11]",
    )
    # (PB) h < 0 on [263/2000, 42/100]  => C2 + P36 = -(3s-1)^2 h/3 >= 0 there,
    # sharp exactly at s = 1/3.
    nonneg_on(h_poly, lo3, hi3, "(PB) h < 0 on [263/2000, 42/100]", strict_negative=True)
    # (PC) J >= 0 on [2/5, 42/100] (both couplings active)
    J = sp.expand(
        96 * s * (11 * s - 2) * (C2 + v[3] / 2 - s**2 * g**2)
        + 48 * s * (11 * s - 2) * (1 - s) ** 2 * q2
        - q2**2
    )
    nonneg_on(J, lo2, hi3, "(PC) J > 0 on [2/5, 42/100]")
    # (PD) I24 >= 0 on [2/5, 13/20] (only (2,4) active)
    I24 = sp.expand(
        96 * s * (11 * s - 2) * C2 + 48 * s * (11 * s - 2) * (1 - s) ** 2 * q2 - q2**2
    )
    nonneg_on(I24, lo2, hi2, "(PD) I24 > 0 on [2/5, 13/20]")
    # sharpness bookkeeping: at s = 1/3 the bound and the ONB energy vanish
    third = sp.Rational(1, 3)
    onb_energy = C2.subs(s, third) + v[3].subs(s, third) + v[6].subs(s, third)
    checks.append(("ONB energy C2 + v3 + v6 = 0 at s = 1/3", sp.nsimplify(onb_energy) == 0))
    checks.append(
        ("B_tight(1/3) = -64/729 (scalar-domination undershoot)",
         b_tight_single_orbit_value(third) == sp.Rational(-64, 729))
    )
    checks.append(
        ("B_tight(1/2) = -1/192 (mode-2 undershoot at the cross profile)",
         b_tight_single_orbit_value(sp.Rational(1, 2)) == sp.Rational(-1, 192))
    )
    return checks


# ----------------------------------------------------------------------------
# 5.  M2: violation hunt for the master bound B(sigma).
# ----------------------------------------------------------------------------


def _atoms_from_vector(x, n_atoms):
    """Map an unconstrained parameter vector to atoms [(alpha, weight)]."""
    import numpy as np

    alphas = 1 / (1 + np.exp(-np.asarray(x[:n_atoms])))  # (0,1)
    raw = np.abs(np.asarray(x[n_atoms:])) + 1e-12
    weights = raw / raw.sum()
    return list(zip(alphas.tolist(), weights.tolist()))


def m2_objective(x, n_atoms, boundary: bool = False):
    """Float B(sigma) with an infeasibility penalty for lambda1 < 1/3."""
    atoms = _atoms_from_vector(x, n_atoms)
    lam1 = sum(wj * al**2 for al, wj in atoms)
    res = master_bound(atoms, exact=False)
    pen = 0.0
    if lam1 < 1 / 3:
        pen = 1e3 * (1 / 3 - lam1) ** 2 + 10.0 * (1 / 3 - lam1)
    if boundary:
        pen += 1e3 * (lam1 - 1 / 3) ** 2
    return res["B"] + pen


def m2_structured_scan(verbose: bool = True):
    """Scan structured configurations; return list of records."""
    import numpy as np

    records = []

    def rec(name, atoms):
        r = master_bound(atoms, exact=False)
        r["name"] = name
        r["atoms"] = atoms
        records.append(r)
        return r

    # (a) single even pair at alpha, alpha^2 >= 1/3
    for al in np.linspace(math.sqrt(1 / 3), 1.0, 41):
        rec("single-pair", [(float(al), 1.0)])
    # (b) pole-equator + one latitude atom
    for wpole in np.linspace(0.0, 0.6, 13):
        for v in np.linspace(0.0, 0.5, 11):
            weq = 1 - wpole - v
            if weq < -1e-9:
                continue
            for al in np.linspace(0.05, 0.95, 19):
                atoms = [(1.0, wpole), (0.0, max(weq, 0.0)), (float(al), v)]
                lam1 = wpole + v * al**2
                if lam1 < 1 / 3 - 1e-9:
                    continue
                rec("family+latitude", atoms)
    # (c) two free pairs + optional pole/equator masses
    rng = np.random.default_rng(0)
    for _ in range(200):
        al1, al2 = rng.uniform(0, 1, 2)
        w1 = rng.uniform(0, 1)
        w2 = 1 - w1
        atoms = [(float(al1), float(w1)), (float(al2), float(w2))]
        if w1 * al1**2 + w2 * al2**2 >= 1 / 3:
            rec("two-pair", atoms)
    records.sort(key=lambda r: r["B"])
    if verbose:
        for r in records[:10]:
            print(
                f"{r['name']:>16}  B={r['B']:+.6f}  c0={r['c0']:+.6f} "
                f"pens=({r['pen1']:.4f},{r['pen2']:.4f},{r['pen3']:.4f},{r['pen4']:.4f}) "
                f"lam1={r['lambda1']:.4f} atoms={[(round(al,4), round(wj,4)) for al, wj in r['atoms']]}"
            )
    return records


def m2_optimize(n_atoms: int = 3, boundary: bool = False, seed: int = 0, maxiter: int = 200):
    """differential_evolution + SLSQP polish on B(sigma); returns best record."""
    import numpy as np
    from scipy.optimize import differential_evolution, minimize

    bounds = [(-6, 6)] * n_atoms + [(-4, 4)] * n_atoms
    de = differential_evolution(
        m2_objective,
        bounds,
        args=(n_atoms, boundary),
        seed=seed,
        maxiter=maxiter,
        tol=1e-12,
        polish=True,
    )
    x = de.x
    loc = minimize(
        m2_objective,
        x,
        args=(n_atoms, boundary),
        method="Nelder-Mead",
        options={"maxiter": 4000, "xatol": 1e-12, "fatol": 1e-14},
    )
    if loc.fun < de.fun:
        x = loc.x
    atoms = _atoms_from_vector(x, n_atoms)
    res = master_bound(atoms, exact=False)
    res["atoms"] = atoms
    res["opt_value"] = min(de.fun, loc.fun)
    return res


def true_E_min_given_profile(atoms, n_phi: int = 6, seed: int = 0, restarts: int = 8):
    """Numerically minimize the true E over antipodal atomic measures with the
    given |a|-profile sigma (atoms = [(alpha_j, w_j)]).

    Each orbit j becomes n_phi points on the circle a = +alpha_j at angles
    phi_{jl} with weights w_j q_{jl} (q a distribution), plus the antipodal
    copies at -alpha_j, phi + pi, each at half weight.  E is computed directly
    from K.  Returns (best E, details).  This is a diagnostic (float) tool:
    if B(sigma) < 0 <= true min E, the leak is in the domination scheme.
    """
    import numpy as np
    from scipy.optimize import minimize

    Kf = sp.lambdify(z, K_POLY, "numpy")
    alphas = np.array([al for al, _ in atoms])
    wts = np.array([wj for _, wj in atoms])
    J = len(atoms)

    def unpack(x):
        phis = x[: J * n_phi].reshape(J, n_phi)
        raw = np.abs(x[J * n_phi :].reshape(J, n_phi)) + 1e-9
        q = raw / raw.sum(axis=1, keepdims=True)
        return phis, q

    def energy(x):
        phis, q = unpack(x)
        # points: for orbit j, both circles; antipodal map phi -> phi + pi
        pts = []
        pw = []
        for j in range(J):
            aj = alphas[j]
            rj = math.sqrt(max(0.0, 1 - aj**2))
            for l in range(n_phi):
                for sgn, dphi in ((1, 0.0), (-1, math.pi)):
                    pts.append(
                        (
                            rj * math.cos(phis[j, l] + dphi),
                            rj * math.sin(phis[j, l] + dphi),
                            sgn * aj,
                        )
                    )
                    pw.append(0.5 * wts[j] * q[j, l])
        P = np.array(pts)
        Wv = np.array(pw)
        Gram = P @ P.T
        np.clip(Gram, -1, 1, out=Gram)
        return float(Wv @ Kf(Gram) @ Wv)

    rng = np.random.default_rng(seed)
    best = None
    for _ in range(restarts):
        x0 = np.concatenate(
            [rng.uniform(0, 2 * math.pi, J * n_phi), rng.uniform(0.5, 1.5, J * n_phi)]
        )
        r = minimize(energy, x0, method="Nelder-Mead", options={"maxiter": 20000, "fatol": 1e-14})
        if best is None or r.fun < best.fun:
            best = r
    return best.fun, best


def realizable_mode_min(atoms, k: int):
    """Tight minimum of Q_k over complex measures nu with |nu| <= sigma
    supported on the atom set (per-atom phase relaxation): for discrete sigma,
    min over zeta_j in C, |zeta_j| <= 1 of  (sum_j w_j zeta_j f(alpha_j))^H
    M_k (same).  Because the form factors through the single vector
    v_j = f(alpha_j), this is a small phase optimization; we solve it by
    real parametrization zeta_j = r_j e^{i th_j} numerically (diagnostic).
    """
    import numpy as np
    from scipy.optimize import minimize

    eps = k % 2
    d = MODE_DIMS[k]
    M = np.array([[float(MODE_MATRICES[k][i, j]) for j in range(d)] for i in range(d)])
    F = np.array(
        [
            [float(al) ** (eps + 2 * i) * (1 - float(al) ** 2) ** (k / 2) for i in range(d)]
            for al, _ in atoms
        ]
    )
    wts = np.array([float(wj) for _, wj in atoms])

    def val(x):
        J = len(atoms)
        r = 1 / (1 + np.exp(-x[:J]))
        th = x[J:]
        zeta = r * np.exp(1j * th)
        zvec = (wts * zeta) @ F
        return float(np.real(np.conj(zvec) @ M @ zvec))

    best = None
    rng = np.random.default_rng(0)
    for _ in range(12):
        x0 = np.concatenate([rng.normal(0, 2, len(atoms)), rng.uniform(0, 2 * math.pi, len(atoms))])
        r = minimize(val, x0, method="Nelder-Mead", options={"maxiter": 8000})
        if best is None or r.fun < best.fun:
            best = r
    return best.fun


# ----------------------------------------------------------------------------
# 6.  CLI
# ----------------------------------------------------------------------------


def _print_checks(checks):
    ok = True
    for name, good in checks:
        print(("PASS " if good else "FAIL ") + name)
        ok = ok and bool(good)
    print(f"{sum(bool(g) for _, g in checks)}/{len(checks)} checks passed")
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("verify", help="fast exact verification of the mode decomposition")
    sub.add_parser("derive", help="slow independent re-derivation by Fourier integrals")
    sub.add_parser("family", help="print B(sigma) on the pole-equator family (exact)")
    m1p = sub.add_parser("m1", help="axisymmetric certificate")
    m1p.add_argument("--solve", action="store_true")
    m1p.add_argument("--deg", type=int, default=8)
    m1p.add_argument("--out", default="sdpa_runs/cylinder/m1_certificate.json")
    m1p.add_argument("--verify", metavar="PATH", default=None)
    m2p = sub.add_parser("m2", help="violation hunt for the master bound")
    m2p.add_argument("--quick", action="store_true")
    m2p.add_argument("--out", default="sdpa_runs/cylinder/m2_scan.json")
    sub.add_parser(
        "pair-theorem",
        help="exact verification of the circle-pair theorem (Toeplitz-coupled bound)",
    )
    args = ap.parse_args(argv)

    if args.cmd == "verify":
        ok = _print_checks(verify_mode_decomposition())
        return 0 if ok else 1
    if args.cmd == "pair-theorem":
        ok = _print_checks(verify_circle_pair_theorem())
        return 0 if ok else 1
    if args.cmd == "derive":
        modes = derive_modes_by_integration()
        for k in range(8):
            print(f"c_{k} =", sp.factor(modes[k]))
        for k in range(7):
            diff = sp.simplify(sp.expand(modes[k] - mode_kernel(k)))
            print(f"c_{k} matches hard-coded matrix:", diff == 0)
        return 0
    if args.cmd == "family":
        for wv in (0, sp.Rational(1, 6), sp.Rational(1, 3), sp.Rational(1, 2), 1):
            atoms = [(1, wv), (0, 1 - wv)]
            res = master_bound(atoms, exact=True)
            print(
                f"w = {wv}: B = {res['B']}  (expected 6(w-1/3)^2 = {sp.expand(6*(wv-sp.Rational(1,3))**2)}), "
                f"pens = {[res[f'pen{k}'] for k in (1,2,3,4)]}"
            )
        return 0
    if args.cmd == "m1":
        if args.verify:
            with open(args.verify) as fh:
                cert = m1_cert_from_json(json.load(fh))
            ok = _print_checks(verify_m1_certificate(cert))
            return 0 if ok else 1
        if args.solve:
            num = m1_build_and_solve(deg=args.deg)
            print("numeric status:", num["status"], "margin:", num["margin"])
            if num["W"] is None:
                print("infeasible at this degree")
                return 1
            cert = m1_rationalize(num)
            checks = verify_m1_certificate(cert)
            ok = _print_checks(checks)
            if ok:
                os.makedirs(os.path.dirname(args.out), exist_ok=True)
                with open(args.out, "w") as fh:
                    json.dump(m1_cert_to_json(cert), fh, indent=1)
                print("exact certificate written to", args.out)
            return 0 if ok else 1
        ap.error("m1 needs --solve or --verify")
    if args.cmd == "m2":
        records = m2_structured_scan()
        out = []
        for r in records:
            rr = dict(r)
            rr.pop("feasible", None)
            out.append(rr)
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(out[:200], fh, indent=1, default=float)
        print("scan written to", args.out)
        if not args.quick:
            for n_atoms in (2, 3, 4):
                r = m2_optimize(n_atoms=n_atoms)
                print(f"n_atoms={n_atoms}: min B = {r['B']:+.6f} atoms={r['atoms']}")
        return 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
