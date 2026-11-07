"""Reader-facing derivations for shape."""
import sympy as s
from .identities import Identities, simp
from .components import *

def derive_first_variation():
    """General symmetric class-A structure and five tracefree shape directions."""
    from itertools import product
    ck = Identities()
    q, N, VOL, kap = s.symbols('q N volume kappa', positive=True)
    H, HD = s.symbols('H HD', real=True)
    n0, n1, n2, n3, n4, n5 = s.symbols('n11 n22 n33 n12 n13 n23', real=True)
    n = s.Matrix([[n0, n3, n4], [n3, n1, n5], [n4, n5, n2]])
    vs = s.symbols('u0:5', real=True)
    vd = s.symbols('ud0:5', real=True)
    vdd = s.symbols('udd0:5', real=True)
    shape = lambda v: s.Matrix([[v[0], v[2], v[3]], [v[2], v[1], v[4]], [v[3], v[4], -v[0] - v[1]]])
    B = shape(vs)
    BD = shape(vd)
    BDD = shape(vdd)
    Dpairs = [(q, H * q), (H, HD)] + list(zip(vs, vd)) + list(zip(vd, vdd))
    D = lambda x: sum((s.diff(x, a) * av for a, av in Dpairs))
    chi = s.expand((s.trace(n) ** 2 / 2 - s.trace(n * n)) / 6)
    X = 2 * n * n - s.trace(n) * n + 4 * chi * s.eye(3)
    ck.zero('tracefree_matrix', s.trace(X))
    C = coframe(n, q, H)
    C1 = s.MutableDenseNDimArray.zeros(4, 4, 4)
    for i, j in product(range(3), repeat=2):
        C1[i + 1, 0, j + 1] = BD[i, j]
        C1[i + 1, j + 1, 0] = -BD[i, j]
    for i, j, k in product(range(3), repeat=3):
        C1[i + 1, j + 1, k + 1] = s.expand(sum((B[i, a] * C[a + 1, j + 1, k + 1] - C[i + 1, a + 1, k + 1] * B[a, j] - C[i + 1, j + 1, a + 1] * B[a, k] for a in range(3))))
    dn = B * n + n * B
    ck.zero('coframe_transform_equals_delta_n', [C1[i + 1, j + 1, k + 1] - sum((dn[i, a] * s.LeviCivita(a, j, k) for a in range(3))) / q for i, j, k in product(range(3), repeat=3)])
    G = lc_connection(C)
    G1 = lc_connection(C1)
    R = curvature(G, C, D)
    R1 = curvature_variation(G, G1, C, C1, D)
    rc = ricci(R)
    sc = scalar(R)
    dsc = scalar(R1)
    ck.zero('independent_full_scalar_background', sc - 6 * (HD + 2 * H ** 2 + chi / q ** 2))
    ck.zero('independent_scalar_shape_variation', dsc + 2 * s.trace(X * B) / q ** 2)
    ck.zero('independent_spatial_Ricci', rc[1:4, 1:4] - (HD + 3 * H ** 2 + 2 * chi / q ** 2) * s.eye(3) - X / q ** 2)
    ck.zero('momentum_background', rc[0, 1:4])
    ck.zero('no_shape_velocity_or_acceleration_linear_scalar', [s.diff(dsc, v) for v in (*vd, *vdd)])
    R3 = (s.trace(n) ** 2 / 2 - s.trace(n * n)) / q ** 2
    delta_R3 = (s.trace(n) * s.trace(dn) - 2 * s.trace(n * dn)) / q ** 2
    L1 = s.expand(VOL * N * q ** 3 * delta_R3 / (2 * kap ** 2))
    target = -VOL * N * q * s.trace(X * B) / kap ** 2
    ck.zero('ADM_action_equals_component_action', L1 - VOL * N * q ** 3 * dsc / (2 * kap ** 2))
    ck.zero('shape_Euler_from_full_coframe', L1 - target)
    E = s.Matrix([simp(s.diff(L1, v)) for v in vs])
    for a in range(5):
        ck.zero('shape_direction_' + str(a), E[a] + VOL * N * q * s.trace(X * B.diff(vs[a])) / kap ** 2)
    eps = s.Symbol('epsilon')
    kin = s.trace((H * s.eye(3) + eps * BD) ** 2) - s.trace(H * s.eye(3) + eps * BD) ** 2
    ck.zero('ADM_kinetic_linear_shape', s.diff(kin, eps).subs(eps, 0))
    ck.zero('GHY_trace_linear_shape', s.diff(s.trace(H * s.eye(3) + eps * BD), eps))
    ck.zero('spatial_volume_linear_shape', s.trace(B))
    ck.zero('GHY_exact_determinant_identity', s.trace(BD))
    a, b, c = s.symbols('a b c', real=True)
    sub = {n0: a, n1: b, n2: c, n3: 0, n4: 0, n5: 0}
    Xd = X.subs(sub).applyfunc(simp)
    chid = simp(chi.subs(sub))
    diff = [simp(Xd[0, 0] - Xd[1, 1]), simp(Xd[1, 1] - Xd[2, 2]), simp(Xd[2, 2] - Xd[0, 0])]
    ck.zero('pair_difference_factorization', [diff[0] - (a - b) * (a + b - c), diff[1] - (b - c) * (b + c - a), diff[2] - (c - a) * (c + a - b)])
    M = s.Matrix([[1, 1, -1], [-1, 1, 1], [1, -1, 1]])
    ck.ok('all_distinct_case_impossible_linear_system', M.det() != 0)
    ck.zero('two_equal_case_reduction', diff[1].subs(b, a) - c * (a - c))
    ck.zero('all_equal_branch', Xd.subs({b: a, c: a}))
    ck.zero('all_equal_chi', chid.subs({b: a, c: a}) - a * a / 4)
    for zero in range(3):
        vals = [a, a, a]
        vals[zero] = 0
        nn = s.diag(*vals)
        ch = (s.trace(nn) ** 2 / 2 - s.trace(nn * nn)) / 6
        ck.zero(f'flat_branch_{zero}', 2 * nn * nn - s.trace(nn) * nn + 4 * ch * s.eye(3))
        ck.zero(f'flat_chi_{zero}', ch)
        ell = [sum(vals) / 2 - v for v in vals]
        ck.zero(f'flat_sectional_{zero}', [vals[i] * ell[i] - ell[(i + 1) % 3] * ell[(i + 2) % 3] for i in range(3)])
    rows = {}
    for name, vals in {'T3': (0, 0, 0), 'Nil': (0, 0, -1), 'Sol': (0, 1, -1), 'S3': (2, 2, 2)}.items():
        ss = {**sub, a: vals[0], b: vals[1], c: vals[2]}
        xx = Xd.subs(dict(zip((a, b, c), vals)))
        ee = E.subs(sub).subs(dict(zip((a, b, c), vals))).applyfunc(simp)
        expected = {'T3': s.zeros(3), 'Nil': s.diag(-1, -1, 2) / 3, 'Sol': s.diag(-4, 2, 2) / 3, 'S3': s.zeros(3)}[name]
        ck.zero(name + '_class_A_residual', xx - expected)
        rows[name] = dict(n=vals, chi=chid.subs(dict(zip((a, b, c), vals))), q2_Xi=xx, shape_Euler=ee)
    momentum = ricci(R1)[0, 1:4].subs(sub).applyfunc(simp)
    return dict(identities=ck.names, chi=chi, q2_Xi=X, shape_matrix=B, delta_R=dsc, delta_L=L1, shape_Euler=E, momentum_linear_diagonal_frame=momentum, pair_differences=diff, classification=dict(isotropic='(a,a,a)', flat='permutations of (a,a,0)', real_scale_vacuum='flat branch only; q constant; local classification, no new compact quotient asserted'), geometries=rows, proof='All distinct eigenvalues force an invertible homogeneous linear system and hence equality at zero. Thus a pair is equal; remaining equation is c*(a-c)=0.')

def derive_shape_action():
    """Shape quadratic action and local torus constraint cone."""
    ck = Identities()
    ep = s.Symbol('epsilon', real=True)
    q0, N0, vol, kap = s.symbols('q0 N0 volume kappa', positive=True)
    u = s.symbols('u0:5', real=True)
    v = s.symbols('v0:5', real=True)
    shape = lambda a: s.Matrix([[a[0], a[2], a[3]], [a[2], a[1], a[4]], [a[3], a[4], -a[0] - a[1]]])
    B = shape(u)
    Bd = shape(v)
    basis = [B.diff(z) for z in u]
    Gram = s.Matrix(5, 5, lambda i, j: s.trace(basis[i] * basis[j]))
    ck.ok('five_shape_Gram_positive', all((x > 0 for x in Gram.eigenvals())))
    n11, n22, n33, n12, n13, n23 = s.symbols('n11 n22 n33 n12 n13 n23', real=True)
    n = s.Matrix([[n11, n12, n13], [n12, n22, n23], [n13, n23, n33]])
    chi = (s.trace(n) ** 2 / 2 - s.trace(n * n)) / 6
    X = 2 * n * n - s.trace(n) * n + 4 * chi * s.eye(3)
    dn = B * n + n * B
    d2n = B * B * n / 2 + B * n * B + n * B * B / 2
    U2 = 2 * s.trace(B * n) ** 2 + 2 * s.trace(n) * s.trace(B * B * n) - 4 * s.trace(B * n * B * n) - 4 * s.trace(B * B * n * n)
    f = lambda m: s.trace(m) ** 2 / 2 - s.trace(m * m)
    fseries = s.expand(f(n + ep * dn + ep ** 2 * d2n))
    ck.zero('potential_first_order', fseries.coeff(ep, 1) + 2 * s.trace(X * B))
    ck.zero('potential_second_order', fseries.coeff(ep, 2) - U2)
    C0 = coframe(n, s.S.One, 0)
    C1 = coframe(dn, s.S.One, 0)
    C2 = coframe(d2n, s.S.One, 0)
    G0 = lc_connection(C0)
    G1 = lc_connection(C1)
    G2 = lc_connection(C2)
    R11 = curvature(G1, C1, lambda z: 0)
    R02 = curvature_variation(G0, G2, C0, C2, lambda z: 0)
    ck.zero('independent_Koszul_second_order', scalar(R11) + scalar(R02) - U2)
    ck.zero('sphere_EH_shape_potential', U2.subs({n11: 2, n22: 2, n33: 2, n12: 0, n13: 0, n23: 0}) + 8 * s.trace(B * B))
    ck.zero('torus_shape_potential', U2.subs({n11: 0, n22: 0, n33: 0, n12: 0, n13: 0, n23: 0}))
    eminus = s.eye(3) - ep * B + ep ** 2 * B * B / 2
    dexp = 2 * ep * Bd + 2 * ep ** 2 * (Bd * B + B * Bd)
    Sm = (eminus * dexp * eminus / 2).applyfunc(s.expand)
    ck.zero('noncommuting_shape_velocity_first', Sm.applyfunc(lambda z: z.coeff(ep, 1)) - Bd)
    ck.zero('noncommuting_shape_velocity_second', Sm.applyfunc(lambda z: z.coeff(ep, 2)))
    ck.zero('volume_preserving_shape', s.trace(B))
    ck.zero('GHY_trace_shape_independent', s.trace(Bd))
    rho, rd, rdd, m, N0d = s.symbols('rho rhod rhodd lapse1 N0d', real=True)
    G = vol * q0 ** 3 / (2 * kap ** 2)
    T = s.trace(Bd * Bd) - 6 * rd ** 2
    L2 = G * T / N0
    EN2 = -G * T / N0 ** 2
    L3 = G * T * (3 * rho - m / N0) / N0
    Nfull = s.Symbol('Nfull', positive=True)
    Lbefore = G * s.exp(3 * ep * rho) * ep ** 2 * T / Nfull
    Lexpanded = s.series(Lbefore.subs(Nfull, N0 + ep * m), ep, 0, 4).removeO().expand()
    ENexpanded = s.series(s.diff(Lbefore, Nfull).subs(Nfull, N0 + ep * m), ep, 0, 3).removeO().expand()
    ck.zero('quadratic_action_from_unexpanded_lapse', Lexpanded.coeff(ep, 2) - L2)
    ck.zero('cubic_action_from_unexpanded_lapse', Lexpanded.coeff(ep, 3) - L3)
    ck.zero('constraint_varied_before_expansion', ENexpanded.coeff(ep, 2) - EN2)
    ck.zero('lapse_constraint_from_cubic_action', s.diff(L3, m) - EN2)
    ck.zero('lapse_missing_from_quadratic_action', s.diff(L2, m))
    ck.zero('scale_shape_velocity_Hessian', s.hessian(L2, [rd, *v]) - 2 * G / N0 * s.diag(s.Integer(-6), Gram))
    cone = s.trace(Bd * Bd) - 6 * rd ** 2
    ck.zero('shape_frozen_limit', cone.subs(dict.fromkeys(v, 0)) + 6 * rd ** 2)
    acc = s.symbols('acc0:5', real=True)
    dt = lambda z: s.diff(z, rd) * rdd + s.diff(z, N0) * N0d + sum((s.diff(z, a) * av for a, av in zip(v, acc)))
    Erho = -dt(s.diff(L2, rd))
    EB = s.Matrix([-dt(s.diff(L2, a)) for a in v])
    ck.zero('radius_linear_Euler', Erho - 12 * G * (rdd / N0 - rd * N0d / N0 ** 2))
    ck.zero('shape_linear_Euler', EB + 2 * G * Gram * (s.Matrix(acc) / N0 - s.Matrix(v) * N0d / N0 ** 2))
    h = s.Symbol('h', real=True, nonzero=True)
    tau, aa = s.symbols('tau a', real=True)
    Sigma = shape(s.symbols('sigma0:5', real=True))
    P = s.eye(3) / 3 + Sigma / (3 * h)
    ck.zero('Kasner_trace', s.trace(P) - 1)
    ck.zero('Kasner_square_cone_equivalence', s.trace(P * P) - 1 - (s.trace(Sigma * Sigma) - 6 * h * h) / (9 * h * h))
    F = 1 + 3 * ep * h * tau
    K = 3 * ep * h * P / F
    theta = s.trace(K)
    ck.zero('exact_Kasner_spatial_Ricci', K.diff(tau) + theta * K)
    R00 = -s.diff(theta, tau) - s.trace(K * K)
    ck.zero('exact_Kasner_temporal_Ricci', R00 + ep ** 2 * (s.trace(Sigma * Sigma) - 6 * h * h) / F ** 2)
    ck.zero('exact_congruence_family_tangent', 2 * aa * s.eye(3) + 2 * B + 6 * h * tau * P - 2 * ((aa + h * tau) * s.eye(3) + B + tau * Sigma))
    ck.zero('volume_family_tangent', s.diff(q0 * s.exp(ep * aa) * F ** s.Rational(1, 3), ep).subs(ep, 0) - q0 * (aa + h * tau))
    ck.zero('Kasner_constraint_propagates', s.diff(theta ** 2 - s.trace(K * K), tau) + 2 * theta * (theta ** 2 - s.trace(K * K)))
    ck.ok('zero_h_cone_forces_zero_shear', all((e > 0 for e in Gram.eigenvals())))
    coords = s.symbols('xcoord ycoord zcoord')
    Km = s.Matrix(3, 3, lambda i, j: s.Symbol(f'K{i}{j}'))
    momentum = s.Matrix([sum((s.diff(Km[j, i] - s.trace(Km) * s.KroneckerDelta(j, i), coords[j]) for j in range(3))) for i in range(3)])
    ck.zero('torus_momentum_spatial_homogeneity', momentum)
    return dict(identities=ck.names, shape=B, Gram=Gram, U2=s.expand(U2), shape_action='V N q^3/(2 kappa^2)*[-6H^2+tr((D B)^2)+(6chi-2tr(XB)+U2)/q^2]+O(B^3)', L2=L2, L3=L3, lapse_equation_second=EN2, linear_radius=Erho, linear_shape=EB, cone='tr(Sigma^2)=6h^2', Kasner_matrix=P, exact_metric='q0^2 exp(epsilon M0) (1+3 epsilon h tau)^(2P) exp(epsilon M0), M0=a I+B0; positive interval', exact_scale='q0 exp(epsilon a) (1+3 epsilon h tau)^(1/3)', integrability='All real homogeneous linearized directions on the cone integrate locally: h!=0 uses the exact family; h=0 forces Sigma=0 and leaves static metric moduli.', endpoints='Local/initial-data or varying-boundary comparison. Holding both metric endpoints fixed at the background removes the affine first-order coefficients.', scope='Finite homogeneous sector on a fixed torus; no full gravitational degree count or global-in-time stability claim.')
