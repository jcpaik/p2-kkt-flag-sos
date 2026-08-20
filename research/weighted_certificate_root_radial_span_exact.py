"""Exact span obstruction for a *linear* rooted radial Hessian certificate.

The calculation uses the projectively correct parities:

* scalar mass spin-2 flags: radial degrees 0,2,4,6;
* tangent positional spin-2 flags: radial degrees 1,3,5.

It asks whether the rooted target ``E/16-D_x`` is a constant linear
compression of the mass/position trace kernels, modulo the root contact
``q(x)=0`` and every polynomial radial contraction of ``grad q(x)=0`` up to
the degree visible here.  The augmented rank is one larger, proving exact
infeasibility of this linear ansatz.  A nonlinear Schur/adjugate step is
therefore essential.
"""

import sympy as sp


a, b, c = sp.symbols("a b c")  # a=Y.Z, b=X.Z, c=X.Y


def kernel(value):
    return 32*value**6-48*value**4+20*value**2-sp.Rational(4, 3)


def position_entry(left, right):
    kp = sp.diff(kernel(a), a)
    kpp = sp.diff(kernel(a), a, 2)
    spin = 2*(a-b*c)**2-(1-b*b)*(1-c*c)
    local = (1-c*c)**2*(
        kpp*(b-a*c)**2-a*kp*(1-c*c)
    )
    cross = spin*(
        kpp*(b-a*c)*(c-a*b)
        + kp*(1-b*b-c*c+a*b*c)
    )
    return sp.expand(
        c**(left+right)*local
        + sp.Rational(1, 2)*(c**left*b**right+c**right*b**left)*cross
    )


def mass_entry(left, right):
    spin = 2*(a-b*c)**2-(1-b*b)*(1-c*c)
    return sp.expand(
        kernel(a)*spin*sp.Rational(1, 2)*(
            c**left*b**right+c**right*b**left
        )
    )


def cross_entry(mass_degree, position_degree):
    kp = sp.diff(kernel(a), a)
    spin = 2*(a-b*c)**2-(1-b*b)*(1-c*c)
    raw = kp*(
        c**(mass_degree+position_degree)*(1-c*c)**2*(b-a*c)
        + c**mass_degree*b**position_degree*spin*(c-a*b)
    )
    return sp.expand((raw+raw.xreplace({b: c, c: b}))/2)


def target():
    determinant = 1+2*a*b*c-a*a-b*b-c*c
    return sp.expand(
        kernel(a)/16-8*b*b*c*c*(a-b*c)**2*determinant
    )


def audit():
    positional = (1, 3, 5)
    mass = (0, 2, 4, 6)
    features = []
    for index, left in enumerate(positional):
        for right in positional[index:]:
            features.append(position_entry(left, right))
    for index, left in enumerate(mass):
        for right in mass[index:]:
            features.append(mass_entry(left, right))
    for left in mass:
        for right in positional:
            features.append(cross_entry(left, right))

    # q(X)=U(X)-E, represented as a two-leaf expectation.
    features.append(sp.expand((kernel(b)+kernel(c))/2-kernel(a)))
    kp = sp.diff(kernel(c), c)
    kpb = sp.diff(kernel(b), b)
    # Dot grad q(X) with every projectively valid radial vector moment.
    for degree in range(1, 18, 2):
        features.append(sp.expand(
            sp.Rational(1, 2)*(kp*b**degree+kpb*c**degree)*(a-b*c)
        ))

    polynomials = [sp.Poly(value, a, b, c) for value in features+[target()]]
    monomials = sorted(set().union(*(poly.monoms() for poly in polynomials)))
    matrix = sp.Matrix([
        [poly.coeff_monomial(monomial) for poly in polynomials[:-1]]
        for monomial in monomials
    ])
    rhs = sp.Matrix([
        polynomials[-1].coeff_monomial(monomial) for monomial in monomials
    ])
    rank = matrix.rank()
    augmented_rank = matrix.row_join(rhs).rank()
    assert rank == 38
    assert augmented_rank == 39
    print("shape", matrix.shape, "rank", rank, "augmented", augmented_rank)


if __name__ == "__main__":
    audit()
