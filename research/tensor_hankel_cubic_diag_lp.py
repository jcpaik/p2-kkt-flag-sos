"""Search diagonal SO(3)-equivariant cubic certificates mass*E=sum v(y)^T H(y)v(y)."""

import itertools
import math

import numpy as np
from scipy.optimize import linprog
from scipy.linalg import null_space
import cvxpy as cp

from tensor_hankel_fw import maps, multinomial, deg6


def sym_generator(degree, L):
    basis = [a for a in itertools.product(range(degree + 1), repeat=3) if sum(a) == degree]
    # Match tensor_hankel_fw.compositions ordering.
    from tensor_hankel_fw import compositions
    basis = compositions(degree, 3)
    index = {a:i for i,a in enumerate(basis)}
    G = np.zeros((len(basis), len(basis)))
    for col, beta in enumerate(basis):
        for source in range(3):
            if beta[source] == 0:
                continue
            for target in range(3):
                if target == source or L[target, source] == 0:
                    continue
                alpha = list(beta)
                alpha[source] -= 1
                alpha[target] += 1
                G[index[tuple(alpha)], col] += L[target, source] * math.sqrt(
                    beta[source] * (beta[target] + 1)
                )
    return G


eps = np.zeros((3,3,3))
eps[0,1,2]=eps[1,2,0]=eps[2,0,1]=1
eps[0,2,1]=eps[2,1,0]=eps[1,0,2]=-1
Ls = [np.array([[eps[i,a,j] for j in range(3)] for i in range(3)]) for a in range(3)]
GV = [sym_generator(3,L) for L in Ls]
GY = [sym_generator(6,L) for L in Ls]


def casimir_projectors(gens, spins):
    C = -sum(G@G for G in gens)
    values, vectors = np.linalg.eigh(C)
    return {l: vectors[:,np.abs(values-l*(l+1))<1e-7] @ vectors[:,np.abs(values-l*(l+1))<1e-7].T
            for l in spins}


PV = casimir_projectors(GV, [1,3])
PY = casimir_projectors(GY, [0,2,4,6])
GH = [np.kron(GY[a], np.eye(10)) + np.kron(np.eye(28), GV[a]) for a in range(3)]
CH = -sum(G@G for G in GH)

projectors=[]; labels=[]; channel_isometries=[]
for av,Pa in PV.items():
    for ly,Pl in PY.items():
        channel=np.kron(Pl,Pa)
        # Orthonormal basis of the channel.
        val,vec=np.linalg.eigh(channel)
        U=vec[:,val>.5]
        Csmall=U.T@CH@U
        ce,cv=np.linalg.eigh(Csmall)
        for j in range(abs(av-ly),av+ly+1):
            take=np.abs(ce-j*(j+1))<1e-6
            if np.count_nonzero(take) != 2*j+1:
                continue
            W=U@cv[:,take]
            projectors.append(W@W.T)
            labels.append((av,ly,j))
            channel_isometries.append(W)
print("channels",len(labels),labels)

# Convert normalized Sym^6 coordinate z to raw moments y.
scale=np.sqrt(np.array([multinomial(a) for a in deg6],float))
M3=maps[3]/scale[None,None,:]
M2=maps[2]/scale[None,None,:]
M1=maps[1]/scale[None,None,:]
mass=np.trace(M3,axis1=0,axis2=1)
Q=(32*np.einsum('abk,abl->kl',M3,M3)-48*np.einsum('abk,abl->kl',M2,M2)
   +20*np.einsum('abk,abl->kl',M1,M1)-(4/3)*np.outer(mass,mass))


def row(z):
    H=np.einsum('abk,k->ab',M3,z)
    A=np.kron(z.reshape(1,-1),np.eye(10))
    return np.array([np.trace(H@A@P@A.T) for P in projectors])

rng=np.random.default_rng(923)
Z=rng.normal(size=(300,28))
A=np.array([row(z) for z in Z])
b=np.array([(mass@z)*(z@Q@z) for z in Z])
res=linprog(np.ones(len(labels)),A_eq=A,b_eq=b,bounds=[(0,None)]*len(labels),method='highs',
            options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
print(res.success,res.message)
if res.success:
    for lab,x in zip(labels,res.x):
        if x>1e-8:print(lab,x)
    print("residual",np.max(np.abs(A@res.x-b)))

# Full invariant Gram: allow cross terms between every pair of copies of spin j.
groups={}
for lab,W in zip(labels,channel_isometries):groups.setdefault(lab[2],[]).append((lab,W))
aligned={}
for j,channels in groups.items():
    ref=channels[0][1]
    Gref=[ref.T@G@ref for G in GH]
    arr=[]
    for lab,W in channels:
        Gcur=[W.T@G@W for G in GH]
        d=2*j+1
        equations=np.vstack([np.kron(np.eye(d),Gc)-np.kron(Gr.T,np.eye(d))
                             for Gc,Gr in zip(Gcur,Gref)])
        ns=null_space(equations)
        if ns.shape[1] != 1:print("bad intertwiner",j,lab,ns.shape)
        X=ns[:,0].reshape((d,d),order='F')
        uu,_,vv=np.linalg.svd(X)
        O=uu@vv
        err=max(np.linalg.norm(Gc@O-O@Gr) for Gc,Gr in zip(Gcur,Gref))
        if err>1e-6:print("align error",j,lab,err)
        arr.append((lab,W@O))
    aligned[j]=arr

Cvars={j:cp.Variable((len(arr),len(arr)),symmetric=True) for j,arr in aligned.items()}
cons=[C >> 0 for C in Cvars.values()]
def cross_rows(z):
    H=np.einsum('abk,k->ab',M3,z); AA=np.kron(z.reshape(1,-1),np.eye(10))
    out={}
    for j,arr in aligned.items():
        R=np.empty((len(arr),len(arr)))
        for c,(_,Ic) in enumerate(arr):
            for d,(_,Id) in enumerate(arr):
                R[c,d]=np.trace(H@AA@Ic@Id.T@AA.T)
        out[j]=(R+R.T)/2
    return out

RR=[cross_rows(z) for z in Z]
for k,z in enumerate(Z):
    cons.append(sum(cp.sum(cp.multiply(Cvars[j],RR[k][j])) for j in aligned)==b[k])
prob=cp.Problem(cp.Minimize(sum((2*j+1)*cp.trace(Cvars[j]) for j in aligned)),cons)
for solver in (cp.CLARABEL,cp.SCS):
    try:
        val=prob.solve(solver=solver,verbose=False)
        print("SDP",solver,prob.status,val)
        if prob.status in ('optimal','optimal_inaccurate'):
            for j,C in Cvars.items():
                ev=np.linalg.eigvalsh(C.value)
                print("j",j,"eig",ev,"matrix",C.value)
            ZZ=rng.normal(size=(100,28));err=[]
            for z in ZZ:
                rr=cross_rows(z);lhs=sum(np.sum(Cvars[j].value*rr[j]) for j in aligned)
                err.append(lhs-(mass@z)*(z@Q@z))
            print("extra residual",max(abs(np.array(err))))
    except cp.error.SolverError as e:print("SDP",solver,e)
