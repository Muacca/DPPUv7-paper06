"""Reader-facing derivations for frames."""
import sympy as s
from .identities import Identities, simp
from .components import *

def derive_global_frames():
    """Periodic coframes on one fixed flat torus."""
    ck = Identities()
    x, y, z = s.symbols('x y z', real=True)
    m = s.symbols('m', integer=True)
    coords = (x, y, z)
    O = s.Matrix([[s.cos(m * z), s.sin(m * z), 0], [-s.sin(m * z), s.cos(m * z), 0], [0, 0, 1]])
    ts = lambda v: s.trigsimp(s.expand_trig(v))
    mat = lambda M: M.applyfunc(ts)

    def wedge(a, b):
        return s.Matrix([a[0] * b[1] - a[1] * b[0], a[0] * b[2] - a[2] * b[0], a[1] * b[2] - a[2] * b[1]])

    def exterior(a):
        return s.Matrix([s.diff(a[j], coords[i]) - s.diff(a[i], coords[j]) for i, j in ((0, 1), (0, 2), (1, 2))])
    sig = [O.row(i).T for i in range(3)]
    ck.zero('same_flat_metric', mat(O.T * O - s.eye(3)))
    ck.zero('oriented_coframe', ts(O.det() - 1))
    ck.zero('periodic_global_frame', mat(O.subs(z, z + 2 * s.pi) - O))
    ck.zero('structure_1', mat(exterior(sig[0]) + m * wedge(sig[1], sig[2])))
    ck.zero('structure_2', mat(exterior(sig[1]) + m * wedge(sig[2], sig[0])))
    ck.zero('structure_3', exterior(sig[2]))
    W = mat(-O.diff(z) * O.T)
    ck.zero('connection_gauge_transform', W - s.Matrix([[0, -m, 0], [m, 0, 0], [0, 0, 0]]))
    for i in range(3):
        ck.zero('torsion_' + str(i), mat(exterior(sig[i]) + sum((W[i, j] * wedge(s.Matrix([0, 0, 1]), sig[j]) for j in range(3)), s.zeros(3, 1))))
    ck.zero('coordinate_connection_constant', W.diff(z))
    ck.zero('coordinate_curvature', s.Matrix([sum((W[i, k] * W[k, j] for k in range(3))) * wedge(s.Matrix([0, 0, 1]), s.Matrix([0, 0, 1]))[r] for i in range(3) for j in range(3) for r in range(3)]))
    n = s.diag(-m, -m, 0)
    C = coframe(n, s.S.One, s.S.Zero)
    G = lc_connection(C)
    ck.zero('independent_Koszul_connection', [G[i + 1, j + 1, k + 1] - (W[i, j] if k == 2 else 0) for i in range(3) for j in range(3) for k in range(3)])
    R = curvature(G, C, lambda f: s.S.Zero)
    ck.zero('independent_full_curvature', [R[idx] for idx in I4])
    ck.zero('parallel_transport_equation', mat(O.diff(z) + W * O))
    ck.zero('z_vector_holonomy', mat(O.subs(z, 2 * s.pi) * O.subs(z, 0).T - s.eye(3)))
    ck.zero('xy_connection_integrands', s.Matrix([W[i, j] * s.diff(z, c) for i in range(3) for j in range(3) for c in (x, y)]))
    a, b, c = s.symbols('a b c', real=True)
    dalpha = mat(exterior(a * sig[0] + b * sig[1] + c * sig[2]))
    ck.zero('constant_coefficient_closed_test', mat(dalpha - (-m * a * wedge(sig[1], sig[2]) - m * b * wedge(sig[2], sig[0]))))
    mnz = s.symbols('mnz', nonzero=True, real=True)
    closed_map = mat(dalpha.jacobian(s.Matrix([a, b, c])))
    ck.ok('restricted_closed_dim_nonzero', closed_map.subs(m, mnz).rank() == 2)
    ck.ok('restricted_closed_dim_zero', closed_map.subs(m, 0).rank() == 0)
    ck.ok('structure_rank_nonzero', s.diag(-mnz, -mnz, 0).rank() == 2)
    ck.ok('structure_rank_zero', n.subs(m, 0).rank() == 0)
    recovered = mat(O.T * O)
    ck.zero('all_three_global_closed_forms', recovered - s.eye(3))
    for i in range(3):
        ck.zero('cartesian_form_closed_' + str(i), exterior(recovered.row(i).T))
    periods = s.Matrix(3, 3, lambda i, j: s.integrate(recovered[i, j], (coords[j], 0, 2 * s.pi)))
    ck.zero('de_Rham_independent_periods', periods - 2 * s.pi * s.eye(3))
    k = s.symbols('kx ky kz', integer=True)
    laplace = sum((s.diff(s.exp(s.I * sum((k[i] * coords[i] for i in range(3)))), v, 2) for v in coords))
    ck.zero('flat_Fourier_Laplacian', laplace + sum((v * v for v in k)) * s.exp(s.I * sum((k[i] * coords[i] for i in range(3)))))
    chi = (s.trace(n) ** 2 / 2 - s.trace(n * n)) / 6
    X = 2 * n * n - s.trace(n) * n + 4 * chi * s.eye(3)
    ck.zero('flat_invariants', [chi, s.trace(X * X), s.trace(X ** 3)])
    return dict(identities=ck.names, periods='x,y,z each have period 2*pi', coframe=O, structure_matrix=n, connection_dz=W, curvature='zero', metric='dx^2+dy^2+dz^2', global_condition='m integer; noninteger m fails single-valued coframe on this fixed torus', comparison=[dict(m=0, rank_n=0, constant_coefficient_closed_dimension=3, b1_T3=3, vector_holonomy='identity'), dict(m='nonzero integer', rank_n=2, constant_coefficient_closed_dimension=1, b1_T3=3, vector_holonomy='identity')], periods_matrix=periods, scope='One fixed flat T3 metric and its global oriented coframes. No spin holonomy or quotient classification; rank n belongs to the chosen homogeneous presentation.')
