"""Independent component curvature from coframe coefficients and Koszul.

C records de^a=(1/2) C^a_cd e^c wedge e^d. G^a_bc is the
connection coefficient along e^c. This module does not import public geometry.
"""
from itertools import product, combinations
import sympy as s
SIGN=(-1,1,1,1)
I3=list(product(range(4),repeat=3)); I4=list(product(range(4),repeat=4))

def lc_connection(C):
    return s.MutableDenseNDimArray([s.expand(SIGN[a]*(SIGN[c]*C[c,b,a]-SIGN[b]*C[b,a,c]-SIGN[a]*C[a,c,b])/2)
                                   for a,b,c in I3],(4,4,4))

def curvature(G,C,D):
    R=s.MutableDenseNDimArray.zeros(4,4,4,4)
    for a,b,c,d in I4:
        v=(D(G[a,b,d]) if c==0 else 0)-(D(G[a,b,c]) if d==0 else 0)
        v+=sum(G[a,e,c]*G[e,b,d]-G[a,e,d]*G[e,b,c]+G[a,b,e]*C[e,c,d] for e in range(4))
        R[a,b,c,d]=s.expand(SIGN[a]*v)
    return R

def curvature_variation(G,G1,C,C1,D):
    R=s.MutableDenseNDimArray.zeros(4,4,4,4)
    for a,b,c,d in I4:
        v=(D(G1[a,b,d]) if c==0 else 0)-(D(G1[a,b,c]) if d==0 else 0)
        v+=sum(G1[a,e,c]*G[e,b,d]+G[a,e,c]*G1[e,b,d]
               -G1[a,e,d]*G[e,b,c]-G[a,e,d]*G1[e,b,c]
               +G1[a,b,e]*C[e,c,d]+G[a,b,e]*C1[e,c,d] for e in range(4))
        R[a,b,c,d]=s.expand(SIGN[a]*v)
    return R

def coframe(n,q,H):
    C=s.MutableDenseNDimArray.zeros(4,4,4)
    for i in range(3): C[i+1,0,i+1]=H; C[i+1,i+1,0]=-H
    for i,j,k in product(range(3),repeat=3):
        C[i+1,j+1,k+1]=sum(n[i,z]*s.LeviCivita(z,j,k) for z in range(3))/q
    return C

def solve_contortion(T):
    slots=[(a,b,c) for a,b in combinations(range(4),2) for c in range(4)]
    vs=s.symbols('K0:24'); K=s.MutableDenseNDimArray.zeros(4,4,4)
    for (a,b,c),v in zip(slots,vs): K[a,b,c]=SIGN[a]*v; K[b,a,c]=-SIGN[b]*v
    eq=[K[a,d,c]-K[a,c,d]-T[a,c,d] for a in range(4) for c,d in combinations(range(4),2)]
    M,rhs=s.linear_eq_to_matrix(eq,vs)
    vals=dict(zip(vs,M.inv()*rhs))
    return s.MutableDenseNDimArray([s.sympify(K[idx]).subs(vals) for idx in I3],(4,4,4))

def ricci(R):
    return s.Matrix(4,4,lambda b,d:s.expand(sum(SIGN[a]*R[a,b,a,d] for a in range(4))))
def scalar(R):
    rc=ricci(R); return s.expand(sum(SIGN[a]*rc[a,a] for a in range(4)))
def quadratic(R):
    i1=sum(SIGN[a]*SIGN[b]*SIGN[c]*SIGN[d]*R[a,b,c,d]**2 for a,b,c,d in I4)
    i2=sum(SIGN[a]*SIGN[b]*SIGN[c]*SIGN[d]*R[a,b,c,d]*R[c,d,a,b] for a,b,c,d in I4)
    rc=ricci(R)
    ir=sum(SIGN[a]*SIGN[b]*rc[a,b]**2 for a,b in product(range(4),repeat=2))
    return tuple(s.factor(s.cancel(v)) for v in (i1,i2,ir))
