"""Tests for the cylindrical mode-domination track (cylinder_cert.py).

Everything checked here is either an exact (rational / symbolic) identity or
an exact certificate verification; float paths are compared against exact
ones.  Run: python -m pytest test_cylinder_cert.py -q
"""

from __future__ import annotations

import json
import math
import os

import pytest
import sympy as sp

import cylinder_cert as cc
from cylinder_cert import (
    MODE_DIMS,
    MODE_MATRICES,
    b_tight_single_orbit_value,
    box_qp_min,
    c0_st_poly,
    master_bound,
    mode_basis_function,
    mode_kernel,
    s,
    signed_abs,
    single_orbit_mode_value,
    t,
    verify_circle_pair_theorem,
    verify_m1_certificate,
    verify_mode_decomposition,
)

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- M0: modes


def test_mode_decomposition_exact():
    checks = verify_mode_decomposition()
    failed = [name for name, ok in checks if not ok]
    assert not failed, f"failed: {failed}"


def test_c1_is_not_zero():
    assert sp.expand(mode_kernel(1)) != 0


def test_mode_dimensions_consistent():
    for k in range(7):
        assert MODE_MATRICES[k].shape == (MODE_DIMS[k], MODE_DIMS[k])
        assert MODE_MATRICES[k] == MODE_MATRICES[k].T


def test_family_sharpness_exact():
    """Pole-equator family: all penalties vanish, B = 6(w - 1/3)^2 exactly."""
    for wv in (sp.Integer(0), sp.Rational(1, 3), sp.Rational(1, 2)):
        res = master_bound([(1, wv), (0, 1 - wv)], exact=True)
        assert all(res[f"pen{k}"] == 0 for k in (1, 2, 3, 4))
        assert sp.simplify(res["B"] - 6 * (wv - sp.Rational(1, 3)) ** 2) == 0
    zero = master_bound([(1, sp.Rational(1, 3)), (0, sp.Rational(2, 3))], exact=True)
    assert zero["B"] == 0 and zero["feasible"]


def test_cross_ledger_exact():
    """The orthonormal 4-point cross: mode ledger sums to E = 2/3 exactly,
    matching the direct kernel sum."""
    half = sp.Rational(1, 2)
    C = c0_st_poly().subs({s: half, t: half})
    v2 = single_orbit_mode_value(2).subs(s, half)
    v4 = single_orbit_mode_value(4).subs(s, half)
    v6 = single_orbit_mode_value(6).subs(s, half)
    assert (C, v2, v4, v6) == (
        sp.Rational(25, 96),
        sp.Rational(-17, 64),
        sp.Rational(21, 32),
        sp.Rational(1, 64),
    )
    assert C + v2 + v4 + v6 == sp.Rational(2, 3)
    # direct: mu uniform on {±u, ±v}, u·v = 0: E = (1/16)(8 K(1) + 8 K(0))
    z = sp.Symbol("z")
    K = 32 * z**6 - 48 * z**4 + 20 * z**2 - sp.Rational(4, 3)
    direct = sp.Rational(1, 16) * (8 * K.subs(z, 1) + 8 * K.subs(z, 0))
    assert direct == sp.Rational(2, 3)


def test_onb_ledger_exact():
    """ONB seen along the cube diagonal: E = C2 + v3 + v6 = 0 at s = 1/3."""
    third = sp.Rational(1, 3)
    C = c0_st_poly().subs({s: third, t: third})
    v3 = single_orbit_mode_value(3).subs(s, third)
    v6 = single_orbit_mode_value(6).subs(s, third)
    assert (C, v3, v6) == (
        sp.Rational(64, 243),
        sp.Rational(-256, 729),
        sp.Rational(64, 729),
    )
    assert C + v3 + v6 == 0


# ---------------------------------------------------------------- box QP


def test_box_qp_exact_simple():
    # mode-4 shape: [[0,-6],[-6,66]] on box [0, m1] x [0, 0] -> 0
    Q = sp.Matrix([[0, -6], [-6, 66]])
    val, _ = box_qp_min(Q, [sp.Rational(1, 2), sp.Integer(0)], exact=True)
    assert val == 0
    # indefinite 2x2 with interior edge minimum:
    # min 66 x^2 - 12 x y over [0, 1] x [0, 1] at (x, y) = (1/11, 1)
    Q2 = sp.Matrix([[66, -6], [-6, 0]])
    val2, arg2 = box_qp_min(Q2, [sp.Integer(1), sp.Integer(1)], exact=True)
    assert val2 == sp.Rational(-6, 11)
    assert arg2 == [sp.Rational(1, 11), 1]


def test_box_qp_float_matches_exact():
    Q = signed_abs(MODE_MATRICES[2])
    m = [sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(1, 8)]
    val_e, _ = box_qp_min(Q, m, exact=True)
    val_f, _ = box_qp_min(Q, [float(x) for x in m], exact=False)
    assert abs(float(val_e) - val_f) < 1e-10


def test_master_bound_float_matches_exact():
    atoms_e = [(1 / sp.sqrt(2), sp.Rational(1, 2)), (0, sp.Rational(1, 2))]
    atoms_f = [(1 / math.sqrt(2), 0.5), (0.0, 0.5)]
    re_ = master_bound(atoms_e, exact=True)
    rf = master_bound(atoms_f, exact=False)
    assert abs(float(re_["B"]) - rf["B"]) < 1e-9


# ---------------------------------------------------------------- M2 verdict


def test_designed_bound_violated_exactly():
    """B_box at the single pair alpha = 1/sqrt2 equals -128159/27456 < 0."""
    res = master_bound([(1 / sp.sqrt(2), 1)], exact=True)
    assert sp.nsimplify(res["B"]) == sp.Rational(-128159, 27456)


def test_scalar_domination_impossibility_exact():
    assert b_tight_single_orbit_value(sp.Rational(1, 2)) == sp.Rational(-1, 192)
    assert b_tight_single_orbit_value(sp.Rational(1, 3)) == sp.Rational(-64, 729)


def test_coupled_bound_value_at_half():
    """(2,4)-coupled single-orbit bound at s = 1/2 equals 1733/14336 > 0."""
    half = sp.Rational(1, 2)
    C2 = c0_st_poly().subs(t, s)
    v2 = single_orbit_mode_value(2)
    q2 = 495 * s**4 - 540 * s**3 + 162 * s**2 - 12 * s + 1
    P24 = v2 / 2 - q2**2 / (96 * s * (11 * s - 2))
    val = sp.nsimplify(sp.simplify((C2 + P24).subs(s, half)))
    assert val == sp.Rational(1733, 14336)


# ---------------------------------------------------------------- theorems


def test_circle_pair_theorem_exact():
    checks = verify_circle_pair_theorem()
    failed = [name for name, ok in checks if not ok]
    assert not failed, f"failed: {failed}"


def test_m1_certificate_exact():
    path = os.path.join(HERE, "sdpa_runs", "cylinder", "m1_certificate.json")
    if not os.path.exists(path):
        pytest.skip("m1 certificate artifact not present (run: cylinder_cert.py m1 --solve)")
    with open(path) as fh:
        cert = cc.m1_cert_from_json(json.load(fh))
    checks = verify_m1_certificate(cert)
    failed = [name for name, ok in checks if not ok]
    assert not failed, f"failed: {failed}"
