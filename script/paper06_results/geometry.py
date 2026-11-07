"""Cartan exterior algebra for a homogeneous class-A coframe.

Signature (-,+,+,+), positive e0123, geometric torsion.
d sigma^i = (1/2) n_ij epsilon_jkl sigma^k wedge sigma^l.
The constant symmetric structure matrix is diagonalized by an SO(3) rotation.
Derivatives D are with respect to proper time; q'=H q, eta'=etaD, b'=bD.
"""
from itertools import permutations, product
import sympy as s

SIGN=(-1,1,1,1)
q=s.Symbol('q',positive=True)
n1,n2,n3,H,eta,etaD,b,bD=s.symbols('n1 n2 n3 H eta etaD b bD',real=True)
HD,etaDD,bDD=s.symbols('HD etaDD bDD',real=True)
ns=(n1,n2,n3)
PAIRS=((q,H*q),(H,HD),(eta,etaD),(etaD,etaDD),(b,bD),(bD,bDD))
E=[{(i,):s.S.One} for i in range(4)]
CYCLIC=((1,2,3),(2,3,1),(3,1,2))

def simp(x):return s.factor(s.cancel(x))
def clean(f):return {k:simp(v) for k,v in f.items() if s.cancel(v)!=0}
def add(*forms):
    out={}
    for f in forms:
        for k,v in f.items():out[k]=out.get(k,0)+v
    return clean(out)
def scale(a,f):return clean({k:a*v for k,v in f.items()})
def wedge(f,g):
    out={}
    for k,a in f.items():
        for l,bv in g.items():
            kl=k+l
            if len(set(kl))!=len(kl):continue
            sign=(-1)**sum(kl[i]>kl[j] for i in range(len(kl)) for j in range(i+1,len(kl)))
            key=tuple(sorted(kl));out[key]=out.get(key,0)+sign*a*bv
    return clean(out)
def D(x):return sum(s.diff(x,v)*vd for v,vd in PAIRS)
def exterior(f,de):
    out={}
    for idx,val in f.items():
        out=add(out,wedge(scale(D(val),E[0]),{idx:1}))
        for j,a in enumerate(idx):
            out=add(out,scale((-1)**j*val,wedge(wedge({idx[:j]:1},de[a]),{idx[j+1:]:1})))
    return out
def tensor_from_forms(R):
    out=s.MutableDenseNDimArray.zeros(4,4,4,4)
    for a,bv in product(range(4),repeat=2):
        for (i,j),v in R[a][bv].items():
            out[a,bv,i,j]=SIGN[a]*v;out[a,bv,j,i]=-SIGN[a]*v
    return out
def ricci(R):return s.Matrix(4,4,lambda j,k:simp(sum(SIGN[a]*R[a,j,a,k] for a in range(4))))
def scalar(R):
    rc=ricci(R);return simp(sum(SIGN[a]*rc[a,a] for a in range(4)))
def internal(R,metric_weights=True):
    out=0
    for c in range(4):
        for d in range(c+1,4):
            weight=SIGN[c]*SIGN[d] if metric_weights else 1
            for a,bv,e,f in permutations(range(4)):
                out-=s.LeviCivita(a,bv,e,f)*weight*R[a,bv,c,d]*R[e,f,c,d]/8
    return simp(out)
def trace2(X,Y):return add(*[wedge(X[a][bv],Y[bv][a]) for a,bv in product(range(4),repeat=2)])
def trace3(X,Y,Z):return add(*[wedge(wedge(X[a][bv],Y[bv][c]),Z[c][a]) for a,bv,c in product(range(4),repeat=3)])
def connection(de,om):
    dom=[[exterior(om[a][bv],de) for bv in range(4)] for a in range(4)]
    torsion=[add(de[a],*[wedge(om[a][bv],E[bv]) for bv in range(4)]) for a in range(4)]
    rf=[[add(dom[a][bv],*[wedge(om[a][z],om[z][bv]) for z in range(4)]) for bv in range(4)] for a in range(4)]
    cs=add(scale(-1,trace2(om,dom)),scale(-s.Rational(2,3),trace3(om,om,om)))
    return dom,torsion,rf,cs
def class_a():
    """Construct all forms before taking any static or shell restriction."""
    de=[{}]+[add(scale(H,wedge(E[0],E[i])),scale(ns[i-1]/q,wedge(E[j],E[k]))) for i,j,k in CYCLIC]
    ell=[sum(ns)/2-ni for ni in ns]
    aa=[(li-eta)/q for li in ell]
    om=[[{} for _ in range(4)] for _ in range(4)]
    for i in range(1,4):om[0][i]=om[i][0]=scale(b,E[i])
    for i,j,k in permutations(range(1,4)):om[j][k]=scale(s.LeviCivita(i-1,j-1,k-1)*aa[i-1],E[i])
    dom,T,rf,CS=connection(de,om)
    Q=add(*[scale(SIGN[a],wedge(E[a],T[a])) for a in range(4)])
    TT=add(*[scale(SIGN[a],wedge(T[a],T[a])) for a in range(4)])
    EER=add(*[scale(SIGN[a],wedge(wedge(E[a],E[bv]),rf[a][bv])) for a,bv in product(range(4),repeat=2)])
    P=scale(-1,trace2(rf,rf))
    # Reference coframe is the same; b=H, eta=0 is its LC connection.
    lc=lambda f:clean({idx:simp(v.subs({eta:0,b:H})) for idx,v in f.items()})
    bar=[[lc(f) for f in row] for row in om]
    # LC curvature is obtained by substituting bD=HD in the full result.
    lcR=lambda f:clean({idx:simp(v.subs({eta:0,etaD:0,b:H,bD:HD})) for idx,v in f.items()})
    Rbar=[[lcR(f) for f in row] for row in rf]
    K=[[add(om[a][bv],scale(-1,bar[a][bv])) for bv in range(4)] for a in range(4)]
    DK=[[add(exterior(K[a][bv],de),*[wedge(bar[a][c],K[c][bv]) for c in range(4)],*[wedge(K[a][c],bar[c][bv]) for c in range(4)]) for bv in range(4)] for a in range(4)]
    trans=add(scale(-2,trace2(K,Rbar)),scale(-1,trace2(K,DK)),scale(-s.Rational(2,3),trace3(K,K,K)))
    return dict(de=de,omega=om,T=T,Rforms=rf,R=tensor_from_forms(rf),Q=Q,TT=TT,EER=EER,NY=exterior(Q,de),P=P,CS=CS,transgression=trans,ell=ell,reference_CS=lc(CS),reference_P=scale(-1,trace2(Rbar,Rbar)))

def forms_json(f):return {','.join(map(str,k)):str(v) for k,v in f.items()}
