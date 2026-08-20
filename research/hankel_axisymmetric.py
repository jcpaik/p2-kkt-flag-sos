"""Exact O(2)-invariant degree-six Hankel/catalectic calculations."""

import math
import sympy as s


def comps(n, k=3):
    if k == 1:
        return [(n,)]
    return [(a,) + tail for a in range(n + 1) for tail in comps(n-a, k-1)]


def multinomial(a):
    out = math.factorial(sum(a))
    for q in a:
        out //= math.factorial(q)
    return out


alpha3 = comps(3)
s0, s1, s2, s3 = s.symbols("s0 s1 s2 s3", real=True)
ss = [s0, s1, s2, s3]


def odd_double_factorial(n):
    if n <= 0:
        return 1
    return math.prod(range(1, n+1, 2))


def moment(alpha):
    ax, ay, az = alpha
    if ax % 2 or ay % 2 or az % 2:
        return s.S.Zero
    a, b, k = ax//2, ay//2, az//2
    n = a+b
    angular = s.Rational(odd_double_factorial(2*a-1)*odd_double_factorial(2*b-1), 2**n*math.factorial(n))
    return angular*ss[k]


H = s.zeros(10)
for i,a in enumerate(alpha3):
    for j,b in enumerate(alpha3):
        scale = s.sqrt(multinomial(a)*multinomial(b))
        H[i,j] = scale*moment(tuple(a[q]+b[q] for q in range(3)))


def partial_trace_matrix(degree):
    aa = comps(degree)
    M = s.zeros(len(aa))
    rem = comps(3-degree)
    for i,a in enumerate(aa):
        for j,b in enumerate(aa):
            scale = s.sqrt(multinomial(a)*multinomial(b))
            for g in rem:
                exponent=tuple(a[q]+b[q]+2*g[q] for q in range(3))
                M[i,j] += scale*multinomial(g)*moment(exponent)
    return M


R2=partial_trace_matrix(2)
R1=partial_trace_matrix(1)
norm=s.expand(s.trace(H))
Q=s.expand(32*s.trace(H*H)-48*s.trace(R2*R2)+20*s.trace(R1*R1)-s.Rational(4,3)*norm**2)

print("basis",alpha3)
print("H="); print(H)
print("norm",norm)
print("Q",s.factor(Q))
print("eigen blocks by support")
# print nonzero graph connected components
seen=set()
for i in range(10):
    if i in seen: continue
    todo=[i]; comp=[]; seen.add(i)
    while todo:
        u=todo.pop(); comp.append(u)
        for v in range(10):
            if v not in seen and H[u,v] != 0:
                seen.add(v); todo.append(v)
    print(comp,[alpha3[k] for k in comp]); print(H.extract(comp,comp))

