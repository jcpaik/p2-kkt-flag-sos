"""Locally optimize the best pair of rank-7 Cayley--Bacharach Hankel rays."""

import itertools
import numpy as np
import torch

torch.set_default_dtype(torch.float64)

P1=np.array([[ .91837175,-.13528845,.37187412],[.92263712,-.13567894,.3610152 ],[.9277229,-.13613493,.34755934],
 [.91617559,-.15787466,.36837193],[.91991059,-.15694875,.35934885],[.92449988,-.15575815,.34790714],
 [.91253725,-.18755912,.36345198],[.91554424,-.18561064,.35682969],[.91938797,-.18305361,.3481625 ]])
W1=np.array([.21655908,.05038459,.0101303,-.00619285,.05316188,.06300365,.02473642,.49575779,.09245914])
P2=np.array([[-.39154487,-.70127311,.59574209],[-.37318725,-.76879297,.51931556],[.17006857,-.66508611,-.72714314],
 [-.36866015,-.71965867,.58837156],[-.34785172,-.7869931,.50954983],[.18362207,-.66383827,-.72498392],
 [-.34435103,-.73819217,.58008162],[-.32105309,-.80506535,.49879324],[.19717192,-.66246103,-.72268155]])
W2=np.array([.10462045,.01410919,.25780987,.3648756,.01278476,.01303482,.02605989,-.00190409,.2086095])

exps=[a for a in itertools.product(range(4),repeat=3) if sum(a)==3]
def relation(P):
 V=np.array([[np.prod(p**np.array(a)) for p in P] for a in exps]);_,_,vh=np.linalg.svd(V)
 u=vh[-1];return u/np.linalg.norm(u)

def init(P,W,neg):
 u=relation(P);inds=[i for i in range(9) if i!=neg]
 pred=-u[neg]**2/sum(u[i]**2/W[i] for i in inds)
 print('relation residual/pred',np.linalg.norm(np.array([[np.prod(p**np.array(a)) for p in P] for a in exps])@u),pred,W[neg])
 return torch.tensor(P),torch.tensor(u),torch.tensor(np.log(W[inds])),inds

p1,u1,z10,i1=init(P1,W1,3);p2,u2,z20,i2=init(P2,W2,7)

def ray(p,u,z,inds,B,neg):
 A=torch.matrix_exp(B)
 raw=p@A.T;r=raw.norm(dim=1);pp=raw/r[:,None];uu=u*r**3
 pos=torch.exp(z);w=torch.zeros(9);w[inds]=pos
 w[neg]=-uu[neg]**2/torch.sum(uu[inds]**2/pos)
 mass=w.sum()
 return pp,w/mass,mass

def kval(g):return 32*g**6-48*g**4+20*g**2-4/3
def objective(par,report=False):
 B1=par[:9].reshape(3,3);B2=par[9:18].reshape(3,3);z1=par[18:26]+z10;z2=par[26:34]+z20;t=torch.sigmoid(par[34])
 q1,w1,m1=ray(p1,u1,z1,i1,B1,3);q2,w2,m2=ray(p2,u2,z2,i2,B2,7)
 pts=torch.cat([q1,q2]);w=torch.cat([t*w1,(1-t)*w2]);val=w@kval(pts@pts.T)@w
 penalty=1e-8*(torch.sum(B1**2)+torch.sum(B2**2))+100*torch.relu(1e-4-m1)**2+100*torch.relu(1e-4-m2)**2
 if report:return val,pts,w,m1,m2,t
 return val+penalty

par=torch.zeros(35,requires_grad=True);par.data[34]=torch.logit(torch.tensor(1-.6658291819728679))
opt=torch.optim.Adam([par],lr=.001)
for it in range(30000):
 opt.zero_grad();v=objective(par);v.backward();torch.nn.utils.clip_grad_norm_([par],10);opt.step()
 if it%1000==0:
  out=objective(par,True);print(it,*[x.item() if x.numel()==1 else None for x in out[::3]])
val,pts,w,m1,m2,t=objective(par,True)
print('FINAL',val.item(),'masses',m1.item(),m2.item(),'t',t.item(),'weights',w.detach().numpy(),'points',pts.detach().numpy(),'par',par.detach().numpy())
