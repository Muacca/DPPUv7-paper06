"""Reader-facing derivations for invariants."""
import sympy as s
from .identities import Identities, simp
from .components import *

def derive_quadratic(geometry_data):
    """Quadratic contractions from a general symmetric class-A coframe."""
    from itertools import product, permutations
    ck = Identities()
    q = s.Symbol('q', positive=True)
    H, HD, eta, etaD, b, bD = s.symbols('H HD eta etaD b bD', real=True)
    na, nb, nc, nd, ne, nf = s.symbols('n1 n2 n3 n12 n13 n23', real=True)
    n = s.Matrix([[na, nd, ne], [nd, nb, nf], [ne, nf, nc]])
    chi = s.expand((s.trace(n) ** 2 / 2 - s.trace(n * n)) / 6)
    X = 2 * n * n - s.trace(n) * n + 4 * chi * s.eye(3)
    J2 = s.expand(s.trace(X * X))
    D = lambda f: sum((s.diff(f, x) * v for x, v in ((q, H * q), (H, HD), (eta, etaD), (b, bD))))
    C = coframe(n, q, H)
    GLC = lc_connection(C)
    T = s.MutableDenseNDimArray.zeros(4, 4, 4)
    for i in range(1, 4):
        T[i, 0, i] = H - b
        T[i, i, 0] = b - H
    for i, j, k in permutations(range(1, 4)):
        T[i, j, k] = 2 * eta * s.LeviCivita(i - 1, j - 1, k - 1) / q
    K = solve_contortion(T)
    G = s.MutableDenseNDimArray([s.expand(GLC[idx] + K[idx]) for idx in I3], (4, 4, 4))
    ck.zero('geometric_torsion_first_Cartan', [C[a, c, d] + G[a, d, c] - G[a, c, d] - T[a, c, d] for a, c, d in I3])
    ck.zero('metric_compatible_connection', [SIGN[a] * G[a, bv, c] + SIGN[bv] * G[bv, a, c] for a, bv, c in I3])
    R = curvature(G, C, D)
    public = geometry_data
    locals_map = {str(v): v for v in (q, H, HD, eta, etaD, b, bD, na, nb, nc)}
    pubtensor = {tuple(map(int, k.split(','))): s.sympify(str(v), locals=locals_map) for k, v in public['tensor'].items()}
    diag = {nd: 0, ne: 0, nf: 0}
    ck.zero('derived_diagonal_tensor_all_256', [R[idx].subs(diag) - pubtensor.get(idx, 0) for idx in I4])
    x = bD + H * b
    y = b * b + (chi - eta * eta) / q ** 2
    z = -etaD / q
    w = 2 * b * eta / q
    i1, i2, ir = quadratic(R)
    f1 = 12 * (x * x + y * y - w * w - z * z) + 4 * J2 / q ** 4
    f2 = 12 * (x * x + y * y + 2 * w * z) + 4 * J2 / q ** 4
    fr = 12 * (x * x + x * y + y * y) + J2 / q ** 4
    ck.zero('I1_general_symmetric_n', i1 - f1)
    ck.zero('I2_general_symmetric_n', i2 - f2)
    ck.zero('Ricci_square_general_symmetric_n', ir - fr)
    ck.zero('pair_exchange_difference', i1 - i2 + 12 * (w + z) ** 2)
    rc = ricci(R)
    rct = s.zeros(4)
    rct[0, 0] = -3 * x
    rct[1:4, 1:4] = (x + 2 * y) * s.eye(3) + X / q ** 2
    ck.zero('Ricci_matrix', rc - rct)
    ck.zero('scalar_curvature', scalar(R) - 6 * (x + y))
    xx, zz, ww, y1, y2, y3 = s.symbols('x z w y1 y2 y3', real=True)
    block = s.MutableDenseNDimArray.zeros(4, 4, 4, 4)

    def entry(a, bv, c, d, v):
        block[a, bv, c, d] = v
        block[bv, a, c, d] = -v
        block[a, bv, d, c] = -v
        block[bv, a, d, c] = v
    for (i, j, k), yi in zip(((1, 2, 3), (2, 3, 1), (3, 1, 2)), (y1, y2, y3)):
        entry(0, i, 0, i, -xx)
        entry(0, i, j, k, -ww)
        entry(j, k, 0, i, zz)
        entry(j, k, j, k, yi)
    bi1, bi2, bir = quadratic(block)
    yb = (y1 + y2 + y3) / 3
    Y2 = sum(((yi - yb) ** 2 for yi in (y1, y2, y3)))
    ck.zero('bivector_I1', bi1 - (12 * (xx ** 2 + yb ** 2 - ww ** 2 - zz ** 2) + 4 * Y2))
    ck.zero('bivector_I2', bi2 - (12 * (xx ** 2 + yb ** 2 + 2 * ww * zz) + 4 * Y2))
    ck.zero('bivector_Ricci_square', bir - (12 * (xx ** 2 + xx * yb + yb ** 2) + Y2))
    lc = {eta: 0, etaD: 0, b: H, bD: HD}
    Rlc = R.applyfunc(lambda v: s.expand(s.sympify(v).subs(lc, simultaneous=True)))
    rclc = ricci(Rlc)
    scalarlc = scalar(Rlc)
    g = s.diag(*SIGN)
    W = s.MutableDenseNDimArray.zeros(4, 4, 4, 4)
    for a, bv, c, d in I4:
        W[a, bv, c, d] = s.expand(Rlc[a, bv, c, d] - (g[a, c] * rclc[bv, d] - g[a, d] * rclc[bv, c] - g[bv, c] * rclc[a, d] + g[bv, d] * rclc[a, c]) / 2 + scalarlc * (g[a, c] * g[bv, d] - g[a, d] * g[bv, c]) / 6)
    ck.zero('LC_pair_symmetry', [Rlc[a, bv, c, d] - Rlc[c, d, a, bv] for a, bv, c, d in I4])
    ck.zero('LC_Weyl_trace', ricci(W))
    W2 = quadratic(W)[0]
    ck.zero('LC_Weyl_square', W2 - 2 * J2 / q ** 4)
    ck.zero('LC_invariant_identity', (i1 - 2 * ir + scalar(R) ** 2 / 3).subs(lc, simultaneous=True) - W2)
    vv = s.Symbol('V', real=True)
    static = {H: 0, bD: 0, etaD: 0, b: -vv / 3}
    ys = vv ** 2 / 9 + (chi - eta ** 2) / q ** 2
    ws = -2 * vv * eta / (3 * q)
    ck.zero('static_I1', i1.subs(static, simultaneous=True) - 12 * (ys ** 2 - ws ** 2) - 4 * J2 / q ** 4)
    ck.zero('static_I2', i2.subs(static, simultaneous=True) - 12 * ys ** 2 - 4 * J2 / q ** 4)
    ck.zero('static_Ricci_square', ir.subs(static, simultaneous=True) - 12 * ys ** 2 - J2 / q ** 4)
    rows = {}
    for name, vals in {'T3': (0, 0, 0), 'Nil': (0, 0, -1), 'Sol': (0, 1, -1), 'S3': (2, 2, 2), 'Sol_half': (0, s.Rational(1, 2), -s.Rational(1, 2))}.items():
        sub = {**diag, na: vals[0], nb: vals[1], nc: vals[2]}
        jj = simp(J2.subs(sub))
        ch = simp(chi.subs(sub))
        ck.zero(name + '_J2', jj - {'T3': 0, 'Nil': s.Rational(2, 3), 'Sol': s.Rational(8, 3), 'S3': 0, 'Sol_half': s.Rational(1, 6)}[name])
        rows[name] = dict(n=vals, chi=ch, J2=jj, q2_Xi=X.subs(sub), LC_Weyl_square=2 * jj / q ** 4)
        if name in ('Nil', 'Sol'):
            for label, expr, coef in [('I1', i1, 4), ('I2', i2, 4), ('Ricci_square', ir, 1)]:
                shell = simp(expr.subs(sub).subs(lc, simultaneous=True).subs(HD, -H ** 2))
                shell = simp(s.expand(shell).subs(H ** 2, -ch / q ** 2))
                ck.zero(name + '_reduced_curve_' + label, shell - coef * jj / q ** 4)
        if name in ('T3', 'S3'):
            ck.zero(name + '_shape_free_curvature_information', jj)
    nil = {**diag, na: 0, nb: 0, nc: -1}
    solhalf = {**diag, na: 0, nb: s.Rational(1, 2), nc: -s.Rational(1, 2)}
    equal = {}
    for label, expr, expect in [('I1', i1, 2 / q ** 4), ('I2', i2, 2 / q ** 4), ('Ricci_square', ir, 1 / (2 * q ** 4)), ('LC_Weyl_square', W2, 1 / q ** 4)]:
        diff = simp(expr.subs(nil) - expr.subs(solhalf))
        ck.zero('equal_chi_' + label, diff - expect)
        equal[label] = diff
    return dict(identities=ck.names, definitions=dict(I1='R_abcd R^abcd', I2='R_abcd R^cdab', Ricci='Ric_bd=eta^ac R_abcd', Ricci_square='Ric_ab Ric^ab', J2='tr((q^2 Xi)^2)'), chi=chi, q2_Xi=X, J2=J2, blocks=dict(x=x, y_mean=y, z=z, w=w), invariants=dict(I1=i1, I2=i2, Ricci_square=ir), compact=dict(I1='12*(x^2+y_mean^2-w^2-z^2)+4*J2/q^4', I2='12*(x^2+y_mean^2+2*w*z)+4*J2/q^4', Ricci_square='12*(x^2+x*y_mean+y_mean^2)+J2/q^4'), pair_difference='I1-I2=-12*(w+z)^2', LC_Weyl_square=2 * J2 / q ** 4, geometries=rows, equal_chi_Nil_minus_Sol_half=equal, reduced_curve='On Nil/Sol reduced curves: I1=I2=4 J2/q^4, Ricci_square=J2/q^4. These are not full vacuum solutions.', scope='Three specified contractions on constant class A, common scale and two-amplitude geometric torsion; not a generating set for all curvature invariants.')

def derive_cubic(geometry_data):
    """Cubic contractions and the spatial Ricci inverse problem."""
    from itertools import product
    ck = Identities()
    a, b, c = s.symbols('a b c', real=True)
    q = s.Symbol('q', positive=True)
    n = s.diag(a, b, c)
    t = s.trace(n)
    chi = (t * t / 2 - s.trace(n * n)) / 6
    X = (2 * n * n - t * n + 4 * chi * s.eye(3)).applyfunc(s.expand)
    J2 = s.expand(s.trace(X * X))
    J3 = s.expand(s.trace(X ** 3))
    lam = s.Symbol('lambda')
    ck.zero('tracefree_cubic_characteristic', (lam * s.eye(3) - X).det() - (lam ** 3 - J2 * lam / 2 - J3 / 3))
    disc = s.discriminant(lam ** 3 - J2 * lam / 2 - J3 / 3, lam)
    ck.zero('real_eigenvalue_discriminant', disc - (J2 ** 3 / 2 - 3 * J3 ** 2))
    ck.zero('discriminant_eigenvalue_differences', disc - (X[0, 0] - X[1, 1]) ** 2 * (X[1, 1] - X[2, 2]) ** 2 * (X[2, 2] - X[0, 0]) ** 2)
    x, y, z, w, d1, d2 = s.symbols('x y z w d1 d2', real=True)
    ds = [d1, d2, -d1 - d2]
    j2 = sum((v * v for v in ds))
    j3 = sum((v ** 3 for v in ds))
    R = s.MutableDenseNDimArray.zeros(4, 4, 4, 4)

    def put(i, j, k, l, val):
        R[i, j, k, l] = val
        R[j, i, k, l] = -val
        R[i, j, l, k] = -val
        R[j, i, l, k] = val
    for (i, j, k), d in zip(((1, 2, 3), (2, 3, 1), (3, 1, 2)), ds):
        put(0, i, 0, i, -x)
        put(0, i, j, k, -w)
        put(j, k, 0, i, z)
        put(j, k, j, k, y - d / q ** 2)
    rc = ricci(R)
    mixed = s.diag(*SIGN) * rc
    Cric = s.expand(s.trace(mixed ** 3))
    Cr = s.S.Zero
    for aa, bb, cc, dd, ee, ff in product(range(4), repeat=6):
        Cr += SIGN[aa] * SIGN[bb] * SIGN[cc] * SIGN[dd] * SIGN[ee] * SIGN[ff] * R[aa, bb, cc, dd] * R[cc, dd, ee, ff] * R[ee, ff, aa, bb]
    Cr = s.expand(Cr)
    fr = 27 * x ** 3 + 3 * (x + 2 * y) ** 3 + 3 * (x + 2 * y) * j2 / q ** 4 + j3 / q ** 6
    fk = 24 * (x ** 3 + y ** 3 + 3 * w * z * (x + y)) + 24 * y * j2 / q ** 4 - 8 * j3 / q ** 6
    ck.zero('mixed_Ricci_cube', Cric - fr)
    ck.zero('oriented_all_index_curvature_cube', Cr - fk)
    ck.zero('spatial_Ricci_cube', s.trace((2 * chi * s.eye(3) + X) ** 3) - (24 * chi ** 3 + 6 * chi * J2 + J3))
    pub = geometry_data
    syms = {name: s.Symbol(name, real=True) for name in ('n1', 'n2', 'n3', 'H', 'b', 'bD', 'eta', 'etaD')}
    syms['q'] = q
    Rp = s.MutableDenseNDimArray.zeros(4, 4, 4, 4)
    for key, value in pub['tensor'].items():
        Rp[tuple(map(int, key.split(',')))] = s.sympify(str(value), locals=syms)
    Rcp = ricci(Rp)
    Cp = s.expand(s.trace((s.diag(*SIGN) * Rcp) ** 3))
    Kp = 0
    for aa, bb, cc, dd, ee, ff in product(range(4), repeat=6):
        Kp += SIGN[aa] * SIGN[bb] * SIGN[cc] * SIGN[dd] * SIGN[ee] * SIGN[ff] * Rp[aa, bb, cc, dd] * Rp[cc, dd, ee, ff] * Rp[ee, ff, aa, bb]
    subn = {a: syms['n1'], b: syms['n2'], c: syms['n3']}
    chip = chi.subs(subn)
    jp2 = J2.subs(subn)
    jp3 = J3.subs(subn)
    xp = syms['bD'] + syms['H'] * syms['b']
    yp = syms['b'] ** 2 + (chip - syms['eta'] ** 2) / q ** 2
    zp = -syms['etaD'] / q
    wp = 2 * syms['b'] * syms['eta'] / q
    ck.zero('derived_tensor_Ricci_cube', Cp - (27 * xp ** 3 + 3 * (xp + 2 * yp) ** 3 + 3 * (xp + 2 * yp) * jp2 / q ** 4 + jp3 / q ** 6))
    ck.zero('derived_tensor_curvature_cube', Kp - (24 * (xp ** 3 + yp ** 3 + 3 * wp * zp * (xp + yp)) + 24 * yp * jp2 / q ** 4 - 8 * jp3 / q ** 6))
    examples = {}
    for name, vals in {'pair_A': (4, 4, 1), 'pair_B': (2, 2, 3), 'Nil': (0, 0, -1), 'Sol': (0, 1, -1), 'T3': (0, 0, 0), 'S3': (2, 2, 2)}.items():
        su = dict(zip((a, b, c), vals))
        examples[name] = dict(n=vals, chi=simp(chi.subs(su)), J2=simp(J2.subs(su)), J3=simp(J3.subs(su)), X=X.subs(su))
    ck.zero('same_chi_example', examples['pair_A']['chi'] - examples['pair_B']['chi'])
    ck.zero('same_J2_example', examples['pair_A']['J2'] - examples['pair_B']['J2'])
    ck.zero('opposite_J3_example', [examples['pair_A']['J3'] + 6, examples['pair_B']['J3'] - 6])
    ck.zero('global_n_sign_ambiguity', [chi.subs({a: -a, b: -b, c: -c}, simultaneous=True) - chi, J2.subs({a: -a, b: -b, c: -c}, simultaneous=True) - J2, J3.subs({a: -a, b: -b, c: -c}, simultaneous=True) - J3])
    l1, l2, l3 = s.symbols('ell1 ell2 ell3', real=True)
    su = {a: l2 + l3, b: l1 + l3, c: l1 + l2}
    rbar = (2 * chi * s.eye(3) + X).subs(su).applyfunc(simp)
    ck.zero('Ricci_product_representation', rbar - s.diag(2 * l2 * l3, 2 * l1 * l3, 2 * l1 * l2))
    rv = [rbar[i, i] for i in range(3)]
    ck.zero('nonzero_inverse_ell_square', [rv[1] * rv[2] - 2 * rv[0] * l1 * l1, rv[2] * rv[0] - 2 * rv[1] * l2 * l2, rv[0] * rv[1] - 2 * rv[2] * l3 * l3])
    ck.zero('Ricci_determinant_nonnegative', rv[0] * rv[1] * rv[2] - 8 * (l1 * l2 * l3) ** 2)
    tt = s.Symbol('t', real=True, nonzero=True)
    deg = (tt + 1 / tt, 1 / tt, tt)
    rdeg = (2 * chi * s.eye(3) + X).subs(dict(zip((a, b, c), deg))).applyfunc(simp)
    ck.zero('degenerate_continuous_inverse_family', rdeg - s.diag(2, 0, 0))
    ck.zero('flat_continuous_inverse_family', (2 * chi * s.eye(3) + X).subs({a: tt, b: tt, c: 0}))
    for name in ('T3', 'S3'):
        ck.zero(name + '_J2_J3_zero', [examples[name]['J2'], examples[name]['J3']])
    return dict(identities=ck.names, chi=chi, J2=J2, J3=J3, characteristic='lambda^3-J2*lambda/2-J3/3', real_spectrum_condition='J2>=0, J2^3>=6*J3^2; additionally constrain class-A realization separately', cubic_definitions=dict(Ricci='tr((Ric^a_b)^3)', curvature='R_ab^cd R_cd^ef R_ef^ab; all ordered indices'), compact=dict(Ricci='27*x^3+3*(x+2*ybar)^3+3*(x+2*ybar)*J2/q^4+J3/q^6', curvature='24*(x^3+ybar^3+3*w*z*(x+ybar))+24*ybar*J2/q^4-8*J3/q^6'), examples=examples, degenerate_inverse_family=deg, inverse='With all spatial Ricci eigenvalues nonzero and positive product, ell is determined up to simultaneous sign (and eigenvalue permutation). Exactly two nonzero Ricci eigenvalues are impossible. With exactly one nonzero, ell has one zero and a fixed product of the other two: continuous ambiguity. With all zero, at most one ell is nonzero: the flat branch.', caveat='Pointwise Ricci spectrum is reconstructed, not compact quotient or local isometry class. A covariant derivative invariant is a possible later diagnostic; not computed.')
