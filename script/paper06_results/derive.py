"""Derive the curvature identities, action, constraints and orbit data."""
from .geometry import *

def require(checks,name,residual):
    ok=(residual=={} if isinstance(residual,dict) else simp(residual)==0)
    if not ok:raise ArithmeticError(f'{name}: {residual}')
    checks.append(name)

def derive_geometry():
    g=class_a();R=g['R'];checks=[]
    for i in range(4):require(checks,f'd_squared_e{i}',exterior(g['de'][i],g['de']))
    for i,j,k in CYCLIC:
        expected=add(scale(H-b,wedge(E[0],E[i])),scale(2*eta/q,wedge(E[j],E[k])))
        require(checks,f'geometric_torsion_{i}',add(g['T'][i],scale(-1,expected)))
    require(checks,'metric_compatibility',sum((SIGN[a]*v+SIGN[bv]*g['omega'][bv][a].get(idx,0))**2 for a,bv in product(range(4),repeat=2) for idx,v in g['omega'][a][bv].items()))
    spatial=R.applyfunc(lambda val:simp(s.sympify(val).subs({H:0,b:0,bD:0,eta:0,etaD:0})))
    spatial_ricci=ricci(spatial)[1:4,1:4]
    chi=simp(scalar(spatial)*q**2/6)
    ell=g['ell'];k=[simp(ns[i]*ell[i]-ell[(i+1)%3]*ell[(i+2)%3]) for i in range(3)]
    require(checks,'mean_sectional_curvature',sum(k)/3-chi)
    require(checks,'structure_matrix_scalar',6*chi-(sum(ns)**2/2-sum(v*v for v in ns)))
    x=bD+H*b;z=-etaD/q;w=2*b*eta/q;y=[b*b+(ki-eta**2)/q**2 for ki in k];ym=b*b+(chi-eta**2)/q**2
    for i,jj,kk in CYCLIC:
        require(checks,f'boost_block_{i}',add(g['Rforms'][0][i],scale(-1,add(scale(x,wedge(E[0],E[i])),scale(w,wedge(E[jj],E[kk]))))))
        require(checks,f'rotation_block_{i}',add(g['Rforms'][jj][kk],scale(-1,add(scale(z,wedge(E[0],E[i])),scale(y[i-1],wedge(E[jj],E[kk]))))))
    require(checks,'complete_NY',add(g['TT'],scale(-1,g['EER']),scale(-1,g['NY'])))
    require(checks,'scalar',scalar(R)-6*(x+ym))
    P=g['P'].get((0,1,2,3),0);BI=internal(R);BF=internal(R,False)
    require(checks,'Pontryagin_blocks',P-12*(ym*z-x*w))
    require(checks,'internal_blocks',BI-3*(w*ym-x*z))
    require(checks,'frame_sum_blocks',BF-3*(w*ym+x*z))
    phi=simp(q**3*g['CS'].get((1,2,3),0));phi0=simp(phi.subs({eta:0}))
    rel=simp(q**3*g['transgression'].get((1,2,3),0))
    require(checks,'CS_derivative',D(phi)-q**3*P)
    require(checks,'relative_polynomial',rel-(-12*chi*eta+4*eta**3-12*eta*q*q*b*b))
    require(checks,'CS_relative_difference',phi-phi0-rel)
    require(checks,'relative_form_derivative',add(exterior(g['transgression'],g['de']),scale(-1,g['P']),g['reference_P']))
    require(checks,'LC_Pontryagin',g['reference_P'])
    require(checks,'NY_scalar',g['NY'][(0,1,2,3)]-6*(etaD+2*H*eta)/q)
    V=s.Symbol('V',real=True)
    static=BI.subs({H:0,bD:0,etaD:0,b:-V/3})
    coefficient=simp(s.Poly(s.expand(static),V,eta).coeff_monomial(V*eta)*9*q**3/2)
    require(checks,'static_C_bridge',coefficient+9*chi)
    require(checks,'static_template',static-2*V*eta*(-V*V*q*q+9*eta*eta+coefficient)/(9*q**3))
    tensor={','.join(map(str,idx)):R[idx] for idx in product(range(4),repeat=4) if R[idx]!=0}
    return g,dict(convention='geometric torsion; positive e0123; signature (-,+,+,+); constant class-A structure matrix',chi=chi,ell=ell,sectional_curvatures=k,spatial_Ricci=spatial_ricci,spatial_tracefree=s.simplify(spatial_ricci-2*chi/q**2*s.eye(3)),x=x,y=y,y_mean=ym,z=z,w=w,R_EC=scalar(R),NY=g['NY'][(0,1,2,3)],TT=g['TT'][(0,1,2,3)],EER=g['EER'][(0,1,2,3)],Q_slice=6*q*q*eta,Pontryagin=P,internal_covariant=BI,legacy_frame_sum=BF,CS_slice=phi,CS_reference=phi0,relative_slice=rel,static_internal=static,C_static=coefficient,tensor=tensor,identities=checks)

def derive_reduction(geo):
    N,kap,vol=s.symbols('N kappa volume',positive=True)
    qd,qdd,Nd,ed,V,Vd,theta,p=s.symbols('qd qdd Nd ed V Vd theta p',real=True)
    pairs=((q,qd),(qd,qdd),(N,Nd),(eta,ed),(V,Vd))
    dt=lambda f:sum(s.diff(f,v)*vd for v,vd in pairs)
    bv=qd/(N*q)-V/3
    sub={H:qd/(N*q),b:bv,bD:dt(bv)/N,etaD:ed/N}
    sc=simp(geo['R_EC'].subs(sub));ny=simp(geo['NY'].subs(sub))
    raw=vol*N*q**3*(sc/(2*kap**2)+theta*ny)
    beh=3*vol*q*q*qd/(kap*kap*N);bec=-vol*q**3*V/kap**2;bny=6*theta*vol*q*q*eta
    L=simp(raw-dt(beh+bec+bny));A=3*vol/kap**2
    F=simp((L+A*q*qd**2/N)/(A*q*N));chi=geo['chi']
    Eqs={name:simp(s.diff(L,v)-dt(s.diff(L,vd))) for name,v,vd in [('q',q,qd),('N',N,Nd),('eta',eta,ed),('V',V,Vd)]}
    moment=s.diff(L,qd);ham=simp((p*qd-L).subs(qd,-N*p/(2*A*q)));constraint=simp(ham/N)
    checks=[]
    require(checks,'bulk_F',F-(chi-eta*eta+q*q*V*V/9))
    require(checks,'Noether_identity',qd*Eqs['q']+ed*Eqs['eta']+Vd*Eqs['V']-N*dt(Eqs['N']))
    hess=s.hessian(F,[eta,V]);require(checks,'auxiliary_determinant',hess.det()+4*q*q/9)
    shell={eta:0,V:0};require(checks,'shell_F',F.subs(shell)-chi)
    require(checks,'Hamiltonian_constraint',constraint-(-p*p/(4*A*q)-A*q*F))
    require(checks,'radius_Euler_shell',Eqs['q'].subs(shell)-A*(N*chi+qd*qd/N+2*q*qdd/N-2*q*qd*Nd/N**2))
    # EH endpoint cancellation computed from its raw acceleration density.
    ehraw=simp(vol*N*q**3*sc.subs({eta:0,ed:0,V:0,Vd:0})/(2*kap**2))
    require(checks,'GHY_delta_qdot',s.diff(ehraw,qdd)-s.diff(beh,qd))
    require(checks,'GHY_delta_N',s.diff(ehraw,Nd)-s.diff(beh,N))
    require(checks,'GHY_delta_q',s.diff(ehraw,qd)-dt(s.diff(ehraw,qdd))-s.diff(beh,q)-moment)
    # Canonical brackets for (p_eta,p_V,eta,V).
    bracket=s.Matrix([[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]])
    require(checks,'auxiliary_second_class_det',bracket.det()-1)
    velocity=s.hessian(L,[qd,Nd,ed,Vd])
    require(checks,'velocity_rank_one_coefficient',velocity[0,0]+2*A*q/N)
    return dict(R_EC=sc,NY=ny,raw=simp(raw),bulk=L,A=A,F=F,BEH=beh,BEC=bec,BNY=bny,
       auxiliary_Hessian=hess,auxiliary_determinant=simp(hess.det()),auxiliary_solution=s.solve([s.diff(F,eta),s.diff(F,V)],[eta,V],dict=True),
       p_q=moment,Hamiltonian=ham,constraint=constraint,Euler=Eqs,
       auxiliary_bracket=bracket,canonical_shift={str(v):s.diff(bec+bny,v) for v in (q,eta,V,N)},
       velocity_Hessian=velocity,torus_constraint_note='After eliminating torsion, C=-p^2/(4 A q) at chi=0 has vanishing gradient on p=0; no regular-constraint dimension count is inferred.',identities=checks)

def derive_classification(geo):
    values={'T3':(0,0,0),'Nil':(0,0,-1),'Sol':(0,1,-1),'S3':(2,2,2)}
    vols={'T3':(2*s.pi)**3,'Nil':(2*s.pi)**4,'Sol':s.log((3+s.sqrt(5))/2),'S3':2*s.pi**2}
    rows={};checks=[]
    for name,n in values.items():
        sub=dict(zip(ns,n));chi=simp(geo['chi'].subs(sub));ric=geo['spatial_Ricci'].subs(sub).applyfunc(simp)
        tf=geo['spatial_tracefree'].subs(sub).applyfunc(simp)
        rows[name]=dict(n=n,chi=chi,volume=vols[name],rank=s.diag(*n).rank(),homogeneous_closed_dimension=3-s.diag(*n).rank(),Ricci=ric,tracefree=tf,C_static=simp(geo['C_static'].subs(sub)),CS_reference=simp(geo['CS_reference'].subs(sub)),full_vacuum_in_scale_family=(name=='T3'),orbit='no_real' if chi>0 else 'static' if chi==0 else 'linear_two_sheets')
        for mode in ('EH','AX','VT','MX'):
            require(checks,f'{name}_{mode}_shell_P',geo['Pontryagin'].subs(sub).subs({eta:0,etaD:0}))
    # Nil and a scaled Sol have exactly the same chi=-1/12.
    equal={name:{key:(val.subs(dict(zip(ns,n))) if hasattr(val,'subs') else val) for key,val in [('chi',geo['chi']),('tracefree',geo['spatial_tracefree']),('CS_reference',geo['CS_reference'])]} for name,n in [('Nil',(0,0,-1)),('scaled_Sol',(0,s.Rational(1,2),-s.Rational(1,2)))]}
    return dict(geometries=rows,equal_chi_comparison=equal,identities=checks)
