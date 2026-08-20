"""Test the fermionic occupation cap and S-gradient transport Hessian."""

import numpy as np


def kernel_prime(t): return 192*t**5-192*t**3+40*t
def kernel_second(t): return 960*t**4-576*t**2+40


def h2_basis():
    out=[]
    out.append(np.diag([1.,-1.,0.])/np.sqrt(2))
    out.append(np.diag([1.,1.,-2.])/np.sqrt(6))
    for i,j in ((0,1),(0,2),(1,2)):
        S=np.zeros((3,3));S[i,j]=S[j,i]=1/np.sqrt(2);out.append(S)
    return out


BASIS=h2_basis()


def cap_matrix(x,w):
    """Matrix of (2/3)I-F in an ON basis of H2."""
    F=np.zeros((5,5))
    for z,weight in zip(x,w):
        vals=[]
        for S in BASIS:
            vals.append(S@z-(z@S@z)*z)
        for i in range(5):
            for j in range(5):
                F[i,j]+=2*weight*(vals[i]@vals[j])
    return 2*np.eye(5)/3-F


def transport_matrix(x,w):
    """Half collective second variation for V_S=P_x Sx."""
    H=np.zeros((5,5))
    fields=[]
    for S in BASIS:
        alpha=np.einsum('ni,ij,nj->n',x,S,x)
        fields.append(x@S-alpha[:,None]*x)
    for a in range(5):
        for b in range(5):
            value=0.
            for i,X in enumerate(x):
                for j,Y in enumerate(x):
                    t=X@Y;va=fields[a][i];vb=fields[b][i]
                    wa=fields[a][j];wb=fields[b][j]
                    # Symmetric bilinearization of local and cross terms.
                    local=kernel_second(t)*(va@Y)*(vb@Y)-t*kernel_prime(t)*(va@vb)
                    cross=.5*kernel_second(t)*((va@Y)*(X@wb)+(vb@Y)*(X@wa))
                    cross+=.5*kernel_prime(t)*(va@wb+vb@wa)
                    value+=w[i]*w[j]*(local+cross)
            H[a,b]=value
    return (H+H.T)/2


def coupled_mass_position_hessian(x,w):
    """Full KKT Hessian for quadratic reweightings and S-gradient transports."""
    g=x@x.T;K=32*g**6-48*g**4+20*g**2-4/3;U=K@w;E=w@U
    vals=np.array([np.einsum('ni,ij,nj->n',x,S,x) for S in BASIS]);means=vals@w;dens=vals-means[:,None]
    fields=[]
    for S in BASIS:
        alpha=np.einsum('ni,ij,nj->n',x,S,x);fields.append(x@S-alpha[:,None]*x)
    fields=np.array(fields)
    M=np.array([[(w*dens[a])@K@(w*dens[b]) for b in range(5)] for a in range(5)])
    P=transport_matrix(x,w);C=np.zeros((5,5))
    for a in range(5):
      for b in range(5):
       total=0.
       for i,X in enumerate(x):
        for j,Y in enumerate(x):
         s=X@Y
         total+=w[i]*w[j]*dens[a,i]*kernel_prime(s)*((fields[b,i]@Y)+(X@fields[b,j]))
       C[a,b]=total
    H=np.block([[M,C],[C.T,P]])
    return (H+H.T)/2,M,C,P


def quadratic_mass_matrix(x,w,stationary_simplification=True):
    """Mass Hessian compressed to densities x^T S x.

    With a stationary measure this is <f_S,K f_T>-E <f_S><f_T>.
    The general centered-density formula is also available.
    """
    g=x@x.T
    K=32*g**6-48*g**4+20*g**2-4/3
    U=K@w;E=w@U
    vals=np.array([np.einsum('ni,ij,nj->n',x,S,x) for S in BASIS])
    means=vals@w
    H=np.zeros((5,5))
    for i in range(5):
        for j in range(5):
            if stationary_simplification:
                H[i,j]=(w*vals[i])@K@(w*vals[j])-E*means[i]*means[j]
            else:
                fi=vals[i]-means[i];fj=vals[j]-means[j]
                H[i,j]=(w*fi)@K@(w*fj)
    return (H+H.T)/2


def rooted_spin4_transport(x,w,root_index=0):
    """Collective Hessian for the two spin-4 meridional fields at a root.

    For root p and tangent frame (u,v), put
      f_c=t((u.y)^2-(v.y)^2), f_s=2t(u.y)(v.y),
      e_theta=(p-t y)/sqrt(1-t^2), V_{c,s}=f_{c,s}e_theta.
    Their L2 Gram is exactly the real spin-4 moment matrix whose determinant
    is D_p/4.
    """
    p=x[root_index]
    seed=np.eye(3)[np.argmin(np.abs(p))]
    u=np.cross(p,seed);u/=np.linalg.norm(u);v=np.cross(p,u)
    t=x@p;zu=x@u;zv=x@v;r=np.sqrt(np.maximum(0,1-t*t))
    et=np.zeros_like(x);mask=r>1e-10;et[mask]=(p-t[mask,None]*x[mask])/r[mask,None]
    f=np.array([t*(zu*zu-zv*zv),2*t*zu*zv])
    fields=f[:,:,None]*et[None,:,:]
    M=np.einsum('n,ain,b in->ab',w,fields,fields,optimize=True) if False else None
    M=np.array([[np.sum(w*np.einsum('ni,ni->n',fields[a],fields[b])) for b in range(2)] for a in range(2)])
    H=np.zeros((2,2))
    for a in range(2):
      for b in range(2):
       total=0.
       for i,X in enumerate(x):
        for j,Y in enumerate(x):
         s=X@Y;va=fields[a,i];vb=fields[b,i];wa=fields[a,j];wb=fields[b,j]
         local=kernel_second(s)*(va@Y)*(vb@Y)-s*kernel_prime(s)*(va@vb)
         cross=.5*kernel_second(s)*((va@Y)*(X@wb)+(vb@Y)*(X@wa))+.5*kernel_prime(s)*(va@wb+vb@wa)
         total+=w[i]*w[j]*(local+cross)
       H[a,b]=total
    return M,(H+H.T)/2


def rooted_spin4_radial_transport(x,w,root_index=0,degrees=(1,3,5)):
    """Radial hierarchy t^r zeta^2 e_theta, split into real/imag parts."""
    p=x[root_index];seed=np.eye(3)[np.argmin(np.abs(p))];u=np.cross(p,seed);u/=np.linalg.norm(u);v=np.cross(p,u)
    t=x@p;zu=x@u;zv=x@v;r=np.sqrt(np.maximum(0,1-t*t));et=np.zeros_like(x);mask=r>1e-10;et[mask]=(p-t[mask,None]*x[mask])/r[mask,None]
    funcs=[];names=[]
    for d in degrees:
        funcs.extend((t**d*(zu*zu-zv*zv),2*t**d*zu*zv));names.extend((f'c{d}',f's{d}'))
    fields=np.array(funcs)[:,:,None]*et[None,:,:];n=len(funcs);H=np.zeros((n,n));M=np.zeros((n,n))
    for a in range(n):
      for b in range(n):
       M[a,b]=np.sum(w*np.einsum('ni,ni->n',fields[a],fields[b]));total=0.
       for i,X in enumerate(x):
        for j,Y in enumerate(x):
         s=X@Y;va=fields[a,i];vb=fields[b,i];wa=fields[a,j];wb=fields[b,j]
         total+=w[i]*w[j]*(kernel_second(s)*(va@Y)*(vb@Y)-s*kernel_prime(s)*(va@vb)+.5*kernel_second(s)*((va@Y)*(X@wb)+(vb@Y)*(X@wa))+.5*kernel_prime(s)*(va@wb+vb@wa))
       H[a,b]=total
    return names,M,(H+H.T)/2


def rooted_spin4_radial_polynomial_transport(
    x,w,root_index=0,degrees=(1,3,5)
):
    """Polynomial radial hierarchy ``t**d*zeta**2 P_{y perp}root``.

    Unlike :func:`rooted_spin4_radial_transport`, this uses the
    *unnormalized* meridional vector ``root-t*y``.  Its traced entries agree
    with :func:`sos_search.rooted_radial_trace_entry_vector` and therefore
    admit exact Gram-polynomial expansion.
    """
    p=x[root_index]
    seed=np.eye(3)[np.argmin(np.abs(p))]
    u=np.cross(p,seed);u/=np.linalg.norm(u);v=np.cross(p,u)
    t=x@p;zu=x@u;zv=x@v
    meridional=p[None,:]-t[:,None]*x
    funcs=[];names=[]
    for d in degrees:
        funcs.extend((t**d*(zu*zu-zv*zv),2*t**d*zu*zv))
        names.extend((f'c{d}',f's{d}'))
    fields=np.array(funcs)[:,:,None]*meridional[None,:,:]
    n=len(funcs);H=np.zeros((n,n));M=np.zeros((n,n))
    for a in range(n):
      for b in range(n):
       M[a,b]=np.sum(w*np.einsum('ni,ni->n',fields[a],fields[b]));total=0.
       for i,X in enumerate(x):
        for j,Y in enumerate(x):
         s=X@Y;va=fields[a,i];vb=fields[b,i];wa=fields[a,j];wb=fields[b,j]
         total+=w[i]*w[j]*(
             kernel_second(s)*(va@Y)*(vb@Y)-s*kernel_prime(s)*(va@vb)
             +.5*kernel_second(s)*((va@Y)*(X@wb)+(vb@Y)*(X@wa))
             +.5*kernel_prime(s)*(va@wb+vb@wa)
         )
       H[a,b]=total
    return names,M,(H+H.T)/2


def rooted_spin4_polynomial_coupled(
    x,w,root_index=0,mass_degrees=(0,2,4),position_degrees=(1,3,5)
):
    """Projectively well-defined rooted mass/position spin-2 block.

    Scalar mass flags have even radial parity, while tangent positional
    flags have odd radial parity.  The positional fields use the polynomial
    meridional vector ``root-t*y``.  Densities are centered so the block is
    a valid zero-mass/positional second variation for any probability
    measure; at a stationary minimizer the whole block is PSD.
    """
    p=x[root_index]
    seed=np.eye(3)[np.argmin(np.abs(p))]
    u=np.cross(p,seed);u/=np.linalg.norm(u);v=np.cross(p,u)
    t=x@p;zu=x@u;zv=x@v;meridional=p[None,:]-t[:,None]*x
    mass_funcs=[];mass_names=[]
    for d in mass_degrees:
        mass_funcs.extend((t**d*(zu*zu-zv*zv),2*t**d*zu*zv))
        mass_names.extend((f'mc{d}',f'ms{d}'))
    position_funcs=[];position_names=[]
    for d in position_degrees:
        position_funcs.extend((t**d*(zu*zu-zv*zv),2*t**d*zu*zv))
        position_names.extend((f'pc{d}',f'ps{d}'))
    mass_funcs=np.array(mass_funcs);means=mass_funcs@w
    densities=mass_funcs-means[:,None]
    position_funcs=np.array(position_funcs)
    fields=position_funcs[:,:,None]*meridional[None,:,:]
    gram=x@x.T;K=32*gram**6-48*gram**4+20*gram**2-4/3
    M=np.einsum('ai,i,ij,j,bj->ab',densities,w,K,w,densities,optimize=True)
    C=np.zeros((len(mass_funcs),len(position_funcs)))
    for a in range(len(mass_funcs)):
      for b in range(len(position_funcs)):
       total=0.
       for i,X in enumerate(x):
        for j,Y in enumerate(x):
         s=X@Y
         total+=w[i]*w[j]*densities[a,i]*kernel_prime(s)*(
             fields[b,i]@Y+X@fields[b,j]
         )
       C[a,b]=total
    _,L,P=rooted_spin4_radial_polynomial_transport(
        x,w,root_index,position_degrees
    )
    H=np.block([[M,C],[C.T,P]])
    return mass_names+position_names,L,(H+H.T)/2,M,C,P


def rooted_spin4_radial_coupled(x,w,root_index=0,degrees=(1,3,5)):
    """Centered mass/position Hessian on the rooted radial spin-4 flags.

    The mass densities are the same real/imaginary functions
    ``t**d*zeta**2`` that multiply the meridional positional field in
    :func:`rooted_spin4_radial_transport`.  Centering makes every density a
    legitimate zero-mass variation for an arbitrary probability measure.
    At a stationary minimizer the returned 12-by-12 block is PSD.
    """
    p=x[root_index]
    seed=np.eye(3)[np.argmin(np.abs(p))]
    u=np.cross(p,seed);u/=np.linalg.norm(u);v=np.cross(p,u)
    t=x@p;zu=x@u;zv=x@v;r=np.sqrt(np.maximum(0,1-t*t))
    et=np.zeros_like(x);mask=r>1e-10
    et[mask]=(p-t[mask,None]*x[mask])/r[mask,None]
    funcs=[];names=[]
    for d in degrees:
        funcs.extend((t**d*(zu*zu-zv*zv),2*t**d*zu*zv))
        names.extend((f'c{d}',f's{d}'))
    funcs=np.array(funcs);means=funcs@w;dens=funcs-means[:,None]
    fields=funcs[:,:,None]*et[None,:,:]
    gram=x@x.T;K=32*gram**6-48*gram**4+20*gram**2-4/3
    M=np.einsum('ai,i,ij,j,bj->ab',dens,w,K,w,dens,optimize=True)
    C=np.zeros((len(funcs),len(funcs)))
    for a in range(len(funcs)):
      for b in range(len(funcs)):
       total=0.
       for i,X in enumerate(x):
        for j,Y in enumerate(x):
         s=X@Y
         total+=w[i]*w[j]*dens[a,i]*kernel_prime(s)*(
             fields[b,i]@Y+X@fields[b,j]
         )
       C[a,b]=total
    _,l2,P=rooted_spin4_radial_transport(x,w,root_index,degrees)
    H=np.block([[M,C],[C.T,P]])
    return names,l2,(H+H.T)/2,M,C,P


def robinson():
    pts=[]
    for zero in range(3):
        other=[i for i in range(3) if i!=zero]
        for sign in (1,-1):
            v=np.zeros(3);v[other]=[1,sign];pts.append(v/np.sqrt(2))
    for s in (1,-1):
        for t in (1,-1):pts.append(np.array([1,s,t])/np.sqrt(3))
    w=np.array([184/2265]*6+[387/3020]*4)
    return np.array(pts),w


def latitude(a,n=360):
    theta=2*np.pi*np.arange(n)/n
    x=np.c_[np.sqrt(1-a)*np.cos(theta),np.sqrt(1-a)*np.sin(theta),np.sqrt(a)*np.ones(n)]
    return x,np.ones(n)/n


def main():
    import scipy.optimize as so
    poly=lambda a:693*a**5-1575*a**4+1260*a**3-420*a**2+54*a-2
    a=so.brentq(poly,.421,8/19)
    for name,(x,w) in [('Robinson',robinson()),('latitude',latitude(a,72))]:
        C=cap_matrix(x,w);H=transport_matrix(x,w);M=quadratic_mass_matrix(x,w)
        print(name,'a' if name=='latitude' else '',a if name=='latitude' else '')
        print('cap eig',np.linalg.eigvalsh(C))
        print('F eig',2/3-np.linalg.eigvalsh(C)[::-1])
        print('transport eig',np.linalg.eigvalsh(H))
        print('quadratic mass eig',np.linalg.eigvalsh(M))
        print('joint cap in most-neg transport',end=' ')
        _,V=np.linalg.eigh(H);v=V[:,0];print(v@C@v)
        RM,RH=rooted_spin4_transport(x,w,0)
        # For latitude quadrature the first point is a valid support root.
        print('root M',RM,'det D',4*np.linalg.det(RM))
        print('root spin4 transport',RH,'eig',np.linalg.eigvalsh(RH))


if __name__=='__main__':main()
