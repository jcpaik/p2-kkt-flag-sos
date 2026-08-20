"""Search the D4/sign-invariant subcone of ternary sextic Hankel PSD matrices."""
import math
import cvxpy as cp
import numpy as np
import sympy as sp

def comps(n,k=3):
    if k==1:return [(n,)]
    return [(a,)+t for a in range(n+1) for t in comps(n-a,k-1)]
def multinomial(a):
    o=math.factorial(sum(a))
    for q in a:o//=math.factorial(q)
    return o
b3=comps(3)

# [x6=y6, z6, x4y2=x2y4, x4z2=y4z2, x2z4=y2z4, x2y2z2]
def idx(alpha):
    if any(q%2 for q in alpha):return None
    x,y,z=alpha
    if (x,y,z) in [(6,0,0),(0,6,0)]:return 0
    if (x,y,z)==(0,0,6):return 1
    if (x,y,z) in [(4,2,0),(2,4,0)]:return 2
    if (x,y,z) in [(4,0,2),(0,4,2)]:return 3
    if (x,y,z) in [(2,0,4),(0,2,4)]:return 4
    if (x,y,z)==(2,2,2):return 5
    raise ValueError(alpha)

def mapmat(deg):
    bb=comps(deg); M=np.zeros((len(bb),len(bb),6)); rem=comps(3-deg)
    for i,a in enumerate(bb):
      for j,b in enumerate(bb):
       sc=math.sqrt(multinomial(a)*multinomial(b))
       for g in rem:
        ex=tuple(a[k]+b[k]+2*g[k] for k in range(3)); ii=idx(ex)
        if ii is not None:M[i,j,ii]+=sc*multinomial(g)
    return M
maps={d:mapmat(d) for d in [1,2,3]}
grams={d:np.einsum('abk,abl->kl',M,M) for d,M in maps.items()}
n=np.trace(maps[3],axis1=0,axis2=1)
Q=32*grams[3]-48*grams[2]+20*grams[1]-4/3*np.outer(n,n)

def energy(y):return y@Q@y
y=cp.Variable(6); direction=cp.Parameter(6)
Hm=sum(y[k]*maps[3][:,:,k] for k in range(6))
prob=cp.Problem(cp.Minimize(direction@y),[Hm>>0,n@y==1])
rng=np.random.default_rng(435)
starts=[]
# uniform sphere
starts.append(np.array([1/7,1/7,1/35,1/35,1/35,1/105]))
for _ in range(100):
 direction.value=rng.normal(size=6);prob.solve(solver='CLARABEL');starts.append(y.value.copy())
best=(1e9,None)
for no,cur in enumerate(starts):
 for it in range(500):
  direction.value=2*Q@cur;prob.solve(solver='CLARABEL');v=y.value.copy();d=v-cur
  aa=d@Q@d;bb=2*cur@Q@d
  step=np.clip(-bb/(2*aa),0,1) if aa>1e-14 else (1 if energy(v)<energy(cur) else 0)
  cur=cur+step*d
  if step<1e-11:break
 val=energy(cur)
 if val<best[0]:best=(val,cur.copy());print('best',no,it,val,cur,np.linalg.eigvalsh(sum(cur[k]*maps[3][:,:,k] for k in range(6))))
print('Qmat');print(sp.Matrix(Q).applyfunc(lambda z:sp.Rational(str(z)).limit_denominator()))
