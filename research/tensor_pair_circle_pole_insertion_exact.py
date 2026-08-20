"""Exact pole--Haar-equator insertion polynomial for the c=288 gap.

The inserted line is q=(sqrt(z),0,sqrt(1-z)); thus z is its squared
equatorial component.  Every circle integral is evaluated by an exact
Fourier-flow count, not numerical quadrature.
"""

from collections import defaultdict
from fractions import Fraction
import itertools
import math
import sys

import sympy as sp

sys.path.insert(0, "research")
import tensor_pair_circle_graph_exact as graph


EPSILON, Z = sp.symbols("epsilon z", real=True)
EDGES = ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))


def angular_flow_average(status, exponents):
    """Average one Gram monomial for statuses Q(inserted), P(pole), E(equator)."""
    equator_vertices = [i for i,value in enumerate(status) if value == "E"]
    equator_index = {vertex:index for index,vertex in enumerate(equator_vertices)}
    flows = {(0,)*len(equator_vertices): 1}
    denominator_power = 0
    z_half_power = 0
    one_minus_z_half_power = 0

    for exponent,(left,right) in zip(exponents,EDGES,strict=True):
        if not exponent:
            continue
        sl,sr = status[left],status[right]
        if "P" in (sl,sr):
            other = sr if sl == "P" else sl
            if other == "E":
                return sp.Integer(0)
            if other == "Q":
                one_minus_z_half_power += exponent
            # P-P contributes one.
            continue
        if sl == "Q" and sr == "Q":
            continue
        if "Q" in (sl,sr):
            z_half_power += exponent

        # Every remaining edge is a cosine: E-E is cos(theta_l-theta_r),
        # Q-E is cos(theta_E) since theta_Q=0.
        denominator_power += exponent
        updated = defaultdict(int)
        for flow,weight in flows.items():
            for choice in range(exponent+1):
                frequency = exponent - 2*choice
                shifted = list(flow)
                if sl == "E":
                    shifted[equator_index[left]] += frequency
                if sr == "E":
                    shifted[equator_index[right]] -= frequency
                updated[tuple(shifted)] += weight*math.comb(exponent,choice)
        flows = updated

    if z_half_power % 2 or one_minus_z_half_power % 2:
        # Antipodal parity should force these terms to vanish/cancel before
        # reaching here; a nonzero half power would signal a normalization bug.
        numerator = flows.get((0,)*len(equator_vertices),0)
        if numerator:
            raise AssertionError((status,exponents,z_half_power,one_minus_z_half_power))
        return sp.Integer(0)
    coefficient = sp.Rational(
        flows.get((0,)*len(equator_vertices),0), 2**denominator_power
    )
    return coefficient * Z**(z_half_power//2) * (1-Z)**(one_minus_z_half_power//2)


def status_weight(status):
    out = sp.Integer(1)
    for value in status:
        if value == "Q":
            out *= EPSILON
        elif value == "P":
            out *= (1-EPSILON)/3
        else:
            out *= 2*(1-EPSILON)/3
    return out


def insertion_polynomial(constant=288, frame_constant=0):
    expression,_ = graph.expansion(
        constant=constant, root_power=2, frame_constant=frame_constant
    )
    terms = graph.polynomial_terms(expression)
    result = sp.Integer(0)
    for status in itertools.product(("Q","P","E"), repeat=4):
        local = sp.Integer(0)
        for coefficient,exponents in terms:
            local += sp.Rational(coefficient.numerator,coefficient.denominator) * angular_flow_average(status,exponents)
        result += status_weight(status)*local
    return sp.factor(sp.expand(result))


def bernstein_coefficients(polynomial, variables, degrees):
    """Exact tensor-product Bernstein coefficients on the unit cube."""
    ordinary = sp.Poly(sp.expand(polynomial), *variables)
    output = {}
    for index in itertools.product(*(range(degree+1) for degree in degrees)):
        value = sp.Integer(0)
        for power in itertools.product(*(range(entry+1) for entry in index)):
            coefficient = ordinary.coeff_monomial(
                sp.prod(variable**entry for variable,entry in zip(variables,power,strict=True))
            )
            multiplier = sp.prod(
                sp.binomial(index_entry,power_entry)/sp.binomial(degree,power_entry)
                for index_entry,power_entry,degree in zip(index,power,degrees,strict=True)
            )
            value += coefficient*multiplier
        output[index] = sp.factor(value)
    return output


def audit():
    polynomial = insertion_polynomial()
    print("G_288(epsilon,z) =")
    print(polynomial)
    print("coefficients in epsilon:")
    p = sp.Poly(sp.expand(polynomial),EPSILON)
    for (degree,),coefficient in reversed(p.terms()):
        print(degree,sp.factor(coefficient))
    print("Gateaux coefficient",sp.factor(sp.diff(polynomial,EPSILON).subs(EPSILON,0)))
    quotient = sp.cancel(polynomial/EPSILON)
    bernstein = bernstein_coefficients(quotient,(EPSILON,Z),(3,11))
    assert all(value >= 0 for value in bernstein.values())
    print("Bernstein rows for G/epsilon (epsilon degree 3, z degree 11):")
    for row in range(4):
        print(row,[bernstein[row,column] for column in range(12)])


if __name__ == "__main__":
    audit()
