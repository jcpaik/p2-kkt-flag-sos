"""Degree-10 root-state completion relaxation for J >= 108 F.

The tangent Pluecker feature w_x lies in

    wedge^2 H_3 = H_1 + H_3 + H_5

and has squared component weights 2/7, 1/3, 8/21.  Its moment W_mu is a
degree-10 coherent state.  The sextic moment fixes the full 11, 13, 15, and
33 blocks, the spin <=6 part of the 35 block, and the spin <=6 part of the
55 block.  Only spin 8 (cross and 55) and spin 10 (55) require moments above
degree six.  This script maximizes F over the resulting PSD completion.

It is a diagnostic for an exact Schur/positive-map lemma; no numerical SDP
output is itself used as a proof.
"""

import itertools
import sys

sys.path.insert(0,"research")

import cvxpy as cp
import numpy as np

import tensor_hankel_fw as moments
import tensor_schur_feasible_scan as rep
import tensor_weighted_f_wedge3_sdp as rooted


PAIRS=list(itertools.combinations(range(7),2))


def additive_wedge(generator):
    out=np.zeros((21,21))
    for a,(i,j) in enumerate(PAIRS):
        for b,(k,l) in enumerate(PAIRS):
            out[a,b]=(
                generator[i,k]*(j==l)+generator[j,l]*(i==k)
                -generator[i,l]*(j==k)-generator[j,k]*(i==l)
            )
    return out


WEDGE_GENERATORS=[additive_wedge(g) for g in rep.g3]
CASIMIR=-sum(g@g for g in WEDGE_GENERATORS)
EIGENVALUES,EIGENVECTORS=np.linalg.eigh(CASIMIR)
IRREP={ell:EIGENVECTORS[:,np.abs(EIGENVALUES-ell*(ell+1))<1e-7]
       for ell in (1,3,5)}
CHANGE=np.hstack([IRREP[ell] for ell in (1,3,5)])
SLICES={1:slice(0,3),3:slice(3,10),5:slice(10,21)}
GENERATORS={
    ell:[IRREP[ell].T@g@IRREP[ell] for g in WEDGE_GENERATORS]
    for ell in (1,3,5)
}


def rectangular_spin_basis(left,right):
    """Orthonormal vec bases for Hom(H_right,H_left) spin sectors."""
    dl=2*left+1;dr=2*right+1
    actions=[np.kron(np.eye(dr),gl)-np.kron(gr.T,np.eye(dl))
             for gl,gr in zip(GENERATORS[left],GENERATORS[right])]
    casimir=-sum(action@action for action in actions)
    eigenvalues,eigenvectors=np.linalg.eigh(casimir)
    return {
        ell:eigenvectors[:,np.abs(eigenvalues-ell*(ell+1))<1e-7]
        for ell in range(abs(left-right),left+right+1)
    }


def symmetric_spin_basis(ell):
    """Orthonormal vec bases in symmetric operators on H_ell."""
    dimension=2*ell+1
    sym=[]
    for i in range(dimension):
        E=np.zeros((dimension,dimension));E[i,i]=1
        sym.append(E.reshape(-1,order="F"))
    for i in range(dimension):
        for j in range(i+1,dimension):
            E=np.zeros((dimension,dimension));E[i,j]=E[j,i]=1/np.sqrt(2)
            sym.append(E.reshape(-1,order="F"))
    embedding=np.stack(sym,axis=1)
    actions=[np.kron(np.eye(dimension),g)-np.kron(g.T,np.eye(dimension))
             for g in GENERATORS[ell]]
    restricted=[embedding.T@action@embedding for action in actions]
    casimir=-sum(action@action for action in restricted)
    eigenvalues,eigenvectors=np.linalg.eigh(casimir)
    return {
        spin:embedding@eigenvectors[:,np.abs(eigenvalues-spin*(spin+1))<1e-7]
        for spin in range(0,2*ell+1,2)
    }


CROSS35=rectangular_spin_basis(3,5)
SYMMETRIC55=symmetric_spin_basis(5)


def actual_data(points,weights):
    y=sum(weight*rooted.degree_six_evaluation(x)
          for x,weight in zip(points,weights))
    raw=np.einsum("abk,k->ab",moments.maps[3],y)
    H=rep.U.T@raw@rep.U
    C2=rooted.compound(H,np.array(rooted.PAIRS7))
    W=sum(weight*np.outer(rooted.tangent_wedge(x),rooted.tangent_wedge(x))
          for x,weight in zip(points,weights))
    return y,H,CHANGE.T@C2@CHANGE,CHANGE.T@W@CHANGE


def completion_bound(points,weights,solver="CLARABEL"):
    y,H,C,W=actual_data(points,weights)
    Z=cp.Variable((21,21),symmetric=True)
    constraints=[Z>>0]
    # These blocks use moments of degree at most six and are fixed fully.
    for left,right in ((1,1),(1,3),(1,5),(3,1),(3,3),(5,1)):
        constraints.append(Z[SLICES[left],SLICES[right]]
                           ==W[SLICES[left],SLICES[right]])

    # H3-H5 has degrees 2,4,6,8.  Fix precisely its <=6 components.
    difference35=cp.vec(
        Z[SLICES[3],SLICES[5]]-W[SLICES[3],SLICES[5]],order="F"
    )
    for spin in (2,3,4,5,6):
        # Odd sectors vanish for an antipodally symmetric measure, but fixing
        # them explicitly makes the relaxation basis-independent.
        basis=CROSS35.get(spin)
        if basis is not None and basis.shape[1]:
            constraints.append(basis.T@difference35==0)

    # H5-H5 contains even spins 0,...,10.  Fix the degree-six-visible part.
    difference55=cp.vec(
        Z[SLICES[5],SLICES[5]]-W[SLICES[5],SLICES[5]],order="F"
    )
    for spin in (0,2,4,6):
        constraints.append(SYMMETRIC55[spin].T@difference55==0)

    problem=cp.Problem(cp.Maximize(cp.sum(cp.multiply(C,Z))),constraints)
    problem.solve(solver=solver)
    Fmax=(16/9)*problem.value
    Factual=(16/9)*np.sum(C*W)
    J=rooted.j_value(y)
    return J,Factual,Fmax,problem.status


def examples(seed=71):
    rng=np.random.default_rng(seed)
    rows=[]
    rows.append(("ONB",np.eye(3),np.ones(3)/3))
    theta=np.arange(16)*np.pi/16
    pe=np.vstack([[0,0,1],np.c_[np.cos(theta),np.sin(theta),np.zeros(16)]])
    rows.append(("pole-equator",pe,np.r_[1/3,np.ones(16)/24]))
    for index in range(5):
        points=rng.normal(size=(8,3));points/=np.linalg.norm(points,axis=1)[:,None]
        rows.append((f"random-{index}",points,rng.dirichlet(np.ones(8))))
    for name,points,weights in rows:
        J,F,Fmax,status=completion_bound(points,weights)
        print(name,status,"J/108",J/108,"F",F,"completion",Fmax,
              "gap",J/108-Fmax)


if __name__=="__main__":
    examples()
