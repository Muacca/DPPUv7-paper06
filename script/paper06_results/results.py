"""Exact recomputation and portable data for the article."""
import csv,hashlib,json,math,platform
from importlib.metadata import version
from pathlib import Path
from .derive import *

def serial(x):
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    if isinstance(x,s.MatrixBase):return [[serial(v) for v in row] for row in x.tolist()]
    if isinstance(x,(bool,int,float)) or x is None:return x
    return str(x)

def dump(path,data):path.write_text(json.dumps(serial(data),ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')

def recompute(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    g,geo=derive_geometry();red=derive_reduction(geo);cls=derive_classification(geo)
    dump(output/'geometry.json',geo);dump(output/'reduction.json',red);dump(output/'classification.json',cls)
    # All static, dynamic, and smooth-shell restrictions are explicit.
    modes={'EH':{eta:0,etaD:0,b:H},'AX':{b:H},'VT':{eta:0,etaD:0},'MX':{}}
    HD=s.Symbol('HD',real=True)
    rows=[]
    for name,info in cls['geometries'].items():
        nsub=dict(zip(ns,info['n']))
        for mode,ms in modes.items():
            ms=dict(ms)
            if mode in ('EH','AX'):ms[bD]=HD
            for regime in ('dynamic','static','auxiliary_shell'):
                sub={}
                if regime=='static':sub={H:0,HD:0,bD:0,etaD:0}
                elif regime=='auxiliary_shell':sub={eta:0,etaD:0,b:H,bD:HD}
                vals={key:simp(geo[key].subs(nsub).subs(ms,simultaneous=True).subs(sub,simultaneous=True)) for key in ('NY','Pontryagin','internal_covariant','Q_slice')}
                rows.append(dict(geometry=name,mode=mode,regime=regime,**vals))
    dump(output/'regimes.json',rows)
    # Curve data are generated from the computed expressions, not figure constants.
    t=s.Symbol('t',real=True);profile=s.Rational(2,5)*s.exp(-t*t)
    pulse=[]
    for name,info in cls['geometries'].items():
        sub={q:1,H:0,b:0,bD:0,eta:profile,etaD:s.diff(profile,t),**dict(zip(ns,info['n']))}
        funcs={key:s.lambdify(t,geo[key].subs(sub),'math') for key in ('NY','Pontryagin','relative_slice','Q_slice')}
        for i in range(161):
            tv=-2+i/40
            pulse.append(dict(geometry=name,t=format(tv,'.16g'),eta=format(float(profile.subs(t,tv)),'.16g'),**{key:format(fn(tv),'.16g') for key,fn in funcs.items()}))
    with (output/'pulse.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(pulse[0]),lineterminator='\n');writer.writeheader();writer.writerows(pulse)
    orbits=[]
    for name,info in cls['geometries'].items():
        chi=info['chi']
        if chi>0:continue
        for direction in ([0] if chi==0 else [-1,1]):
            slope=direction*math.sqrt(float(-chi))
            for i in range(201):
                tv=i/50;qq=1+slope*tv
                if qq<=0:continue
                orbits.append(dict(geometry=name,sheet=direction,tau=format(tv,'.16g'),q=format(qq,'.16g'),qprime=format(slope,'.16g')))
    with (output/'orbits.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(orbits[0]),lineterminator='\n');writer.writeheader();writer.writerows(orbits)
    # Source-free coframe equations, with an independent proper-time acceleration HD.
    lc=g['R'].applyfunc(lambda val:simp(s.sympify(val).subs({eta:0,etaD:0,b:H,bD:HD})))
    ri=ricci(lc);sc=scalar(lc);ein=ri-sc*s.diag(*SIGN)/2
    chi=geo['chi'];tracefree=geo['spatial_tracefree'];checks=[]
    require(checks,'Einstein_00',ein[0,0]-3*(H*H+chi/q**2))
    require(checks,'Einstein_0i',sum(ein[0,i]**2 for i in range(1,4)))
    for i in range(3):require(checks,f'Einstein_spatial_{i}',ein[i+1,i+1]-(-2*HD-3*H*H-chi/q**2+tracefree[i,i]))
    shell_ein=ein.subs(HD,-H**2).applyfunc(lambda val:simp(s.expand(val).subs(H**2,-chi/q**2)))
    dump(output/'field_equations.json',dict(Einstein=ein,shell_Einstein=shell_ein,identities=checks))
    # Shape, invariant, global-frame and canonical results are derived from
    # the current coframe and action, without reading previously saved results.
    from .shape import derive_first_variation, derive_shape_action
    from .torus import derive_fixed_shape, derive_point_coordinates, derive_endpoints
    from .invariants import derive_quadratic, derive_cubic
    from .evolution import derive_short_time
    from .frames import derive_global_frames
    from .canonical import derive_representatives
    first=derive_first_variation(); shaped=derive_shape_action()
    fixed=derive_fixed_shape(); point=derive_point_coordinates(); ends=derive_endpoints()
    quad=derive_quadratic(geo); cubic=derive_cubic(geo)
    evolution=derive_short_time(); frames=derive_global_frames(); canonical=derive_representatives(red)
    extra=[first,shaped,fixed,point,ends,quad,cubic,evolution,frames,canonical]
    shape_keys=('shape','Gram','U2','shape_action')
    dump(output/'shape_variation.json',dict(first_order=first,second_order={k:shaped[k] for k in shape_keys}))
    mixed={k:v for k,v in shaped.items() if k not in shape_keys}
    mixed['figure_slice']=dict(basis=['diag(1,-1,0)/sqrt(2)','diag(1,1,-2)/sqrt(6)'],
        constraint='sigma_1^2+sigma_2^2=6*h^2',h_range=[-1,1],h_samples=61,angle_samples=101,
        normalization=6,description='Two orthonormal shear directions; a three-dimensional slice of the five-shear constraint cone.')
    dump(output/'torus_deformations.json',dict(fixed_shape=fixed,point_coordinates=point,endpoints=ends,mixed_shape=mixed))
    dump(output/'curvature_invariants.json',dict(quadratic=quad,cubic=cubic))
    numerical=[]
    for name,case in evolution['geometries'].items():
        for branch,points in case.pop('branches').items():
            for pt in points:
                beta=pt.pop('beta')
                numerical.append(dict(geometry=name,branch=int(branch),**pt,beta1=beta[0],beta2=beta[1],beta3=beta[2]))
    dump(output/'short_time_evolution.json',evolution)
    with (output/'short_time_checks.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(numerical[0]),lineterminator='\n');writer.writeheader();writer.writerows(numerical)
    dump(output/'global_frames.json',frames)
    dump(output/'canonical_representatives.json',canonical)
    meta=dict(schema='paper06-v2',geometry='Fixed class-A homogeneous spatial metric with one time-dependent scale',
        pulse={'q':1,'N':1,'V':0,'eta':'(2/5)*exp(-t^2)','interval':[-2,2],'interpretation':'Off-shell. Equal endpoint primitives; local densities need not vanish.'},
        orbit={'q0':1,'proper_time_range':[0,4],'domain':'q>0','interpretation':'Stationary curves of the reduced action. Nil and Sol violate the omitted tracefree coframe equation.'},
        conventions={'torsion':'geometric','volume_form':'e0123','signature':'-+++','sphere':'physical radius q'},
        shape_domain='After full torsion elimination: homogeneous tracefree shape, local torus families and diagonal Nil/Sol evolution.',
        numerical_control=dict(method='DOP853',q0=1,times=[.01,.02,.04],rtol=2e-13,atol=2e-15,max_step=.002,
            scope='Twelve finite-time controls, not certified interval error bounds.'),
        python=platform.python_version(),dependencies={name:version(name) for name in ('sympy','numpy','scipy','matplotlib')},
        check_scope='Symbolic identities and finite numerical controls; the count is not a count of independent theorems.',
        sympy=s.__version__,identities=sum(len(x['identities']) for x in (geo,red,cls))+len(checks)+sum(len(x['identities']) for x in extra),
        files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.iterdir())})
    dump(output/'manifest.json',meta)
    return meta
