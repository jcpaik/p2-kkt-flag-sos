"""Search a rank-seven Hankel pseudo-state mixed with a rotated copy."""
import numpy as np
from scipy.optimize import differential_evolution, minimize

def K(t):return 32*t**6-48*t**4+20*t**2-4/3
def grid_data(xs,ys):
    pts=[];u=[]
    for i,x in enumerate(xs):
      px=np.prod([x-xs[k] for k in range(3) if k!=i])
      for j,y in enumerate(ys):
       py=np.prod([y-ys[k] for k in range(3) if k!=j])
       raw=np.array([x,y,1.]);r=np.linalg.norm(raw)
       pts.append(raw/r);u.append(r**3/(px*py))
    return np.array(pts),np.array(u)
def weights_from_z(u,z,negative):
    inds=[i for i in range(9) if i!=negative];pos=np.exp(z);w=np.zeros(9);w[inds]=pos
    w[negative]=-u[negative]**2/np.sum(u[inds]**2/pos)
    return w

def rot(v):
 a,b,c=v
 ca,sa=np.cos(a),np.sin(a);cb,sb=np.cos(b),np.sin(b);cc,sc=np.cos(c),np.sin(c)
 Rz1=np.array([[ca,-sa,0],[sa,ca,0],[0,0,1.]])
 Ry=np.array([[cb,0,sb],[0,1,0],[-sb,0,cb]])
 Rz2=np.array([[cc,-sc,0],[sc,cc,0],[0,0,1.]])
 return Rz1@Ry@Rz2
def selfq(w,p):return w@K(p@p.T)@w
def cross(ang,w,p):return w@K(p@rot(ang)@p.T)@w
rng=np.random.default_rng(195)
best=(1e9,None)
for tr in range(200):
 xs=np.sort(rng.normal(size=3));ys=np.sort(rng.normal(size=3));p,u=grid_data(xs,ys);neg=int(rng.integers(9));z=rng.normal(scale=.6,size=8);w=weights_from_z(u,z,neg);sm=w.sum()
 if sm<=0 or -w[neg]/sm<.01:continue
 w/=sm;q=selfq(w,p)
 res=differential_evolution(lambda x:cross(x,w,p),[(0,2*np.pi),(0,np.pi),(0,2*np.pi)],popsize=8,maxiter=80,polish=True,tol=1e-10)
 mix=(q+res.fun)/2
 if mix<best[0]:best=(mix,(xs,ys,neg,w,q,res.fun,res.x,z));print('best',tr,mix,'q',q,'cross',res.fun,'angle',res.x,'negwt',w[neg]);print(xs,ys,w)

# Jointly optimize positive CB weights and relative rotation for the best grid.
xs,ys,neg,w,q,cr,ang,z=best[1];p,u=grid_data(xs,ys)
def joint(v):
 ww=weights_from_z(u,v[:8],neg);sm=ww.sum()
 if sm<=1e-5:return 100+(sm-1)**2
 ww/=sm
 return (selfq(ww,p)+cross(v[8:],ww,p))/2
rr=minimize(joint,np.r_[z,ang],method='L-BFGS-B',bounds=[(-8,8)]*8+[(0,2*np.pi),(0,np.pi),(0,2*np.pi)],options={'maxiter':10000,'ftol':1e-15,'gtol':1e-12,'maxls':100})
ww=weights_from_z(u,rr.x[:8],neg);ww/=ww.sum()
print('JOINT',rr.fun,rr.success,rr.message);print('grid',xs,ys,'neg',neg);print('z',rr.x[:8]);print('w',ww);print('ang',rr.x[8:],'q',selfq(ww,p),'cross',cross(rr.x[8:],ww,p))
