"""Exact signed-permutation invariant Hankel/catalectic calculations."""
import math
import sympy as s

def comps(n,k=3):
    if k==1:return [(n,)]
    return [(a,)+t for a in range(n+1) for t in comps(n-a,k-1)]
def multinomial(a):
    o=math.factorial(sum(a))
    for q in a:o//=math.factorial(q)
    return o
a,b,c=s.symbols('a b c',real=True)
def moment(alpha):
    if any(q%2 for q in alpha):return 0
    pat=sorted(alpha,reverse=True)
    return {(6,0,0):a,(4,2,0):b,(2,2,2):c}[tuple(pat)]
def mat(deg):
    aa=comps(deg); M=s.zeros(len(aa)); rem=comps(3-deg)
    for i,u in enumerate(aa):
      for j,v in enumerate(aa):
       scale=s.sqrt(multinomial(u)*multinomial(v))
       for g in rem:
        exp=tuple(u[k]+v[k]+2*g[k] for k in range(3))
        M[i,j]+=scale*multinomial(g)*moment(exp)
    return M
H,R2,R1=mat(3),mat(2),mat(1)
n=s.trace(H)
Q=s.factor(32*s.trace(H*H)-48*s.trace(R2*R2)+20*s.trace(R1*R1)-s.Rational(4,3)*n*n)
print('basis',comps(3)); print(H); print('n',n,'Q',Q)
print('eig',H.eigenvals())
