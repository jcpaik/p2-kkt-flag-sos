"""SO(3)-trace transport-Hessian inequalities and an LP proof search."""

import sympy as sp
import numpy as np
from scipy.optimize import linprog

s = sp.symbols("s")
K = 32*s**6 - 48*s**4 + 20*s**2 - sp.Rational(4, 3)
kp = sp.diff(K, s)
kpp = sp.diff(kp, s)


def transport(l):
    P = sp.legendre(l, s)
    Pp = sp.diff(P, s)
    Ppp = sp.diff(Pp, s)
    p1 = sp.Rational(l*(l+1), 2)
    A = p1*(1-s**2)
    B = l*(l+1)
    C = (1-s**2)**2*Ppp - s*(1-s**2)*Pp
    D = -s*(1-s**2)*Ppp + (1+s**2)*Pp
    return sp.expand(kpp*(A+C) + kp*(-s*B+D))


def transport_toroidal(l):
    P = sp.legendre(l, s)
    Pp = sp.diff(P, s)
    Ppp = sp.diff(Pp, s)
    p1 = sp.Rational(l*(l+1), 2)
    A = p1*(1-s**2)
    B = l*(l+1)
    C = -(1-s**2)*Pp
    D = -(1-s**2)*Ppp + 2*s*Pp
    return sp.expand(kpp*(A+C) + kp*(-s*B+D))


def legendre_coefficients(poly, max_degree):
    # Integral projection, exact.
    return [sp.factor(sp.Rational(2*n+1, 2)*sp.integrate(poly*sp.legendre(n,s),(s,-1,1)))
            for n in range(max_degree+1)]


if __name__ == "__main__":
    ls = list(range(2, 22, 2))
    maxdeg = max(ls)+6
    tls = list(range(1, 22, 2))
    local_trace = sp.expand(kpp*(1-s**2)/2 - s*kp)
    labels = [("L",0)] + [("G",l) for l in ls] + [("T",l) for l in tls]
    polys = [local_trace] + [transport(l) for l in ls] + [transport_toroidal(l) for l in tls]
    coeffs = [legendre_coefficients(p,maxdeg) for p in polys]
    target = legendre_coefficients(K,maxdeg)
    for label,c in zip(labels,coeffs):
        print(label,[(n,c[n]) for n in range(maxdeg+1) if c[n]])

    evens = list(range(0,maxdeg+1,2))
    # Variables alpha_l and beta_n; target = sum alpha T_l + sum beta P_n.
    Aeq=[];beq=[]
    for n in evens:
        Aeq.append([float(c[n]) for c in coeffs] + [1.0 if m==n else 0.0 for m in evens])
        beq.append(float(target[n]))
    result=linprog(np.zeros(len(labels)+len(evens)),A_eq=np.array(Aeq),b_eq=np.array(beq),
                   bounds=[(0,None)]*(len(labels)+len(evens)),method="highs")
    print(result.success,result.message)
    if result.success:
        for label,v in zip(labels,result.x[:len(labels)]):
            if v>1e-9: print("alpha",label,v)
        for n,v in zip(evens,result.x[len(labels):]):
            if v>1e-9: print("beta",n,v)
