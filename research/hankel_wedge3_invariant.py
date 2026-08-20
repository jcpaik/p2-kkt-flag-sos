"""Test invariant third-compound certificates for trace(H) Q(H).

Set ``HANKEL_TARGET=J`` to test the strengthened weighted target
``144 tr(H^2)-216 tr(R2^2)+87 tr(R1^2)-5 tr(H)^2``.  The default remains
the original target E.
"""
import itertools, math, os
import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

from hankel_ratio_power import maps, mass, Q

if os.environ.get("HANKEL_TARGET", "E").upper() == "J":
    grams = {degree: np.einsum("abk,abl->kl", matrix, matrix)
             for degree, matrix in maps.items()}
    Q = (144 * grams[3] - 216 * grams[2] + 87 * grams[1]
         - 5 * np.outer(mass, mass))
    print("target J")
else:
    print("target E")

# Robust SO(3) generators on Sym^3 R^3 from the tensor-cube model.
words=list(itertools.product(range(3),repeat=3));wi={w:i for i,w in enumerate(words)}
comps=[(a,b,3-a-b) for a in range(4) for b in range(4-a)]
S=np.zeros((27,10))
for j,a in enumerate(comps):
 W=[w for w in words if tuple(w.count(i) for i in range(3))==a]
 for w in W:S[wi[w],j]=1/np.sqrt(len(W))
eps=np.zeros((3,3,3));eps[0,1,2]=eps[1,2,0]=eps[2,0,1]=1;eps[0,2,1]=eps[2,1,0]=eps[1,0,2]=-1
gens=[]
for axis in range(3):
 L=np.array([[eps[i,axis,j] for j in range(3)] for i in range(3)])
 G=np.kron(np.kron(L,np.eye(3)),np.eye(3))+np.kron(np.kron(np.eye(3),L),np.eye(3))+np.kron(np.kron(np.eye(3),np.eye(3)),L)
 gens.append(S.T@G@S)

trip=list(itertools.combinations(range(10),3));ti={a:i for i,a in enumerate(trip)}
trip_array=np.array(trip)
def wedge_gen(g):
 out=np.zeros((120,120))
 for col,I in enumerate(trip):
  for pos,i in enumerate(I):
   for j in range(10):
    if g[j,i]==0:continue
    J=list(I);J[pos]=j
    if len(set(J))<3:continue
    inv=sum(J[a]>J[b] for a in range(3) for b in range(a+1,3));Js=tuple(sorted(J))
    out[ti[Js],col]+=(-1)**inv*g[j,i]
 return out
wg=[wedge_gen(g) for g in gens];Cas=-sum(g@g for g in wg)
vals=np.linalg.eigvalsh(Cas);ells=[]
for l in range(13):
 if np.sum(abs(vals-l*(l+1))<1e-6):ells.append(l)
print('decomp',[(l,int(np.sum(abs(vals-l*(l+1))<1e-6))) for l in ells])
Ps=[]
I=np.eye(120)
for l in ells:
 P=I.copy()
 for j in ells:
  if j!=l:P=P@(Cas-j*(j+1)*I)/(l*(l+1)-j*(j+1))
 Ps.append((P+P.T)/2)

def compound3(H):
 return np.linalg.det(H[trip_array[:,None,:,None],trip_array[None,:,None,:]])
rng=np.random.default_rng(99);rows=[];rhs=[]
for s in range(100):
 y=rng.normal(size=28);H=np.einsum('abk,k->ab',maps[3],y);C=compound3(H)
 rows.append([np.sum(P*C) for P in Ps]);rhs.append((mass@y)*(y@Q@y))
rows=np.array(rows);rhs=np.array(rhs)
w,res,rank,sv=np.linalg.lstsq(rows,rhs,rcond=None);print('ls',w,'res',np.linalg.norm(rows@w-rhs),rank)
ww=cp.Variable(len(Ps),nonneg=True);pr=cp.Problem(cp.Minimize(cp.norm(rows@ww-rhs)));pr.solve(solver='CLARABEL');print('nonneg',pr.value,ww.value)

# Full symmetric SO(3) commutant, including multiplicity couplings.
comm_basis=[]; block_basis=[]
ce,CU=np.linalg.eigh(Cas)
for l in ells:
 E=CU[:,abs(ce-l*(l+1))<1e-6];d=E.shape[1];rg=[E.T@g@E for g in wg]
 sym=[]
 for i in range(d):
  for j in range(i,d):
   X=np.zeros((d,d));X[i,j]=X[j,i]=1 if i==j else 1/np.sqrt(2)
   sym.append(X)
 L=np.stack([np.concatenate([(g@X-X@g).reshape(-1) for g in rg]) for X in sym],axis=1)
 Z=np.eye(len(sym)) if np.linalg.norm(L)<1e-7 else null_space(L,rcond=1e-9)
 Bs=[]
 for k in range(Z.shape[1]):
  B=sum(Z[a,k]*sym[a] for a in range(len(sym)));Bs.append((B+B.T)/2);comm_basis.append(E@B@E.T)
 block_basis.append(Bs)
 print('comm',l,d,len(Bs))

rng=np.random.default_rng(991);RR=[];bb=[]
for ss in range(80):
 yy=rng.normal(size=28);HH=np.einsum('abk,k->ab',maps[3],yy);CC=compound3(HH)
 RR.append([np.sum(PP*CC) for PP in comm_basis]);bb.append((mass@yy)*(yy@Q@yy))
RR=np.array(RR);bb=np.array(bb);coef,res,rank,sv=np.linalg.lstsq(RR,bb,rcond=None);print('full ls rank/res',rank,np.linalg.norm(RR@coef-bb),coef)
vars=[];constraints=[];off=0
for Bs in block_basis:
 vv=cp.Variable(len(Bs));vars.append(vv);constraints.append(sum(vv[k]*Bs[k] for k in range(len(Bs))) >> 0);off+=len(Bs)
allv=cp.hstack(vars)
constraints.append(RR@allv==bb)
pp=cp.Problem(cp.Minimize(cp.sum([cp.trace(sum(vars[j][k]*block_basis[j][k] for k in range(len(block_basis[j])))) for j in range(len(block_basis))])),constraints)
for sol in ('CLARABEL','SCS'):
 try:
  value=pp.solve(solver=sol,verbose=False);print('full psd',sol,pp.status,value,'res',np.linalg.norm(RR@allv.value-bb) if allv.value is not None else None)
 except Exception as ex:print(sol,ex)

# Joint degree-three certificate: tr(H) * <W2,wedge^2 H> + <W3,wedge^3 H>.
def exterior_data(k):
 inds=list(itertools.combinations(range(10),k));ii={a:i for i,a in enumerate(inds)};arr=np.array(inds);dim=len(inds)
 def eg(g):
  out=np.zeros((dim,dim))
  for col,Ix in enumerate(inds):
   for pos,i in enumerate(Ix):
    for j in range(10):
     if g[j,i]==0:continue
     J=list(Ix);J[pos]=j
     if len(set(J))<k:continue
     inv=sum(J[a]>J[b] for a in range(k) for b in range(a+1,k));out[ii[tuple(sorted(J))],col]+=(-1)**inv*g[j,i]
  return out
 egs=[eg(g) for g in gens];cas=-sum(g@g for g in egs);ce,CU=np.linalg.eigh(cas);ls=[l for l in range(13) if np.sum(abs(ce-l*(l+1))<1e-6)]
 allB=[];blocks=[]
 for l in ls:
  E=CU[:,abs(ce-l*(l+1))<1e-6];d=E.shape[1];rg=[E.T@g@E for g in egs];sym=[]
  for i in range(d):
   for j in range(i,d):
    X=np.zeros((d,d));X[i,j]=X[j,i]=1 if i==j else 1/np.sqrt(2);sym.append(X)
  L=np.stack([np.concatenate([(g@X-X@g).reshape(-1) for g in rg]) for X in sym],axis=1)
  Z=np.eye(len(sym)) if np.linalg.norm(L)<1e-7 else null_space(L,rcond=1e-9);Bs=[]
  for z in range(Z.shape[1]):
   B=sum(Z[a,z]*sym[a] for a in range(len(sym)));B=(B+B.T)/2;Bs.append(B);allB.append(E@B@E.T)
  blocks.append(Bs)
 print('ext',k,[(l,int(np.sum(abs(ce-l*(l+1))<1e-6))) for l in ls],'comm',list(map(len,blocks)))
 return arr,allB,blocks

pair_array,B2s,blocks2=exterior_data(2)
def compound2(H):return np.linalg.det(H[pair_array[:,None,:,None],pair_array[None,:,None,:]])
rng=np.random.default_rng(177);R2=[];R3=[];bt=[]
for ss in range(100):
 yy=rng.normal(size=28);HH=np.einsum('abk,k->ab',maps[3],yy);C2=compound2(HH);C3=compound3(HH);mm=mass@yy
 R2.append([mm*np.sum(B*C2) for B in B2s]);R3.append([np.sum(B*C3) for B in comm_basis]);bt.append(mm*(yy@Q@yy))
R2=np.array(R2);R3=np.array(R3);bt=np.array(bt);RJ=np.c_[R2,R3]
cc,rr,rk,sv=np.linalg.lstsq(RJ,bt,rcond=None);print('joint ls',rk,np.linalg.norm(RJ@cc-bt),len(cc))
v2=[];v3=[];con=[]
for Bs in blocks2:
 vv=cp.Variable(len(Bs));v2.append(vv);con.append(sum(vv[k]*Bs[k] for k in range(len(Bs)))>>0)
for Bs in block_basis:
 vv=cp.Variable(len(Bs));v3.append(vv);con.append(sum(vv[k]*Bs[k] for k in range(len(Bs)))>>0)
av2=cp.hstack(v2);av3=cp.hstack(v3);con.append(R2@av2+R3@av3==bt)
pj=cp.Problem(cp.Minimize(sum(cp.trace(sum(v2[j][k]*blocks2[j][k] for k in range(len(blocks2[j])))) for j in range(len(blocks2)))+sum(cp.trace(sum(v3[j][k]*block_basis[j][k] for k in range(len(block_basis[j])))) for j in range(len(block_basis)))),con)
for sol in ('CLARABEL','SCS'):
 try:
  value=pj.solve(solver=sol,verbose=False);print('joint psd',sol,pj.status,value,np.linalg.norm(R2@av2.value+R3@av3.value-bt) if av2.value is not None else None)
 except Exception as ex:print(sol,ex)

# Add invariant SOS-in-moments quadratic remainders, one per harmonic sector.
def op_cas(X):return -sum(g@(g@X-X@g)-(g@X-X@g)@g for g in gens)
def op_proj(X,l):
 out=X.copy()
 for j in (0,2,4,6):
  if j!=l:out=(op_cas(out)-j*(j+1)*out)/(l*(l+1)-j*(j+1))
 return out
GL=[]
for l in (0,2,4,6):
 PM=np.stack([op_proj(maps[3][:,:,k],l) for k in range(28)],axis=2);GL.append(np.einsum('abk,abl->kl',PM,PM))
RS=[]
rng=np.random.default_rng(177)
for ss in range(100):
 yy=rng.normal(size=28);RS.append([(mass@yy)*(yy@G@yy) for G in GL])
RS=np.array(RS)
ssv=cp.Variable(4,nonneg=True);con2=list(con[:-1])+[R2@av2+R3@av3+RS@ssv==bt]
pj2=cp.Problem(cp.Minimize(sum(cp.trace(sum(v2[j][k]*blocks2[j][k] for k in range(len(blocks2[j])))) for j in range(len(blocks2)))+sum(cp.trace(sum(v3[j][k]*block_basis[j][k] for k in range(len(block_basis[j])))) for j in range(len(block_basis)))+cp.sum(ssv)),con2)
for sol in ('CLARABEL','SCS'):
 try:
  value=pj2.solve(solver=sol,verbose=False);print('joint+sos',sol,pj2.status,value,'s',ssv.value,'res',np.linalg.norm(R2@av2.value+R3@av3.value+RS@ssv.value-bt) if av2.value is not None else None)
 except Exception as ex:print(sol,ex)
