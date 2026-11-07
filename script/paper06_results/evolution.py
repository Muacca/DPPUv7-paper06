"""Reader-facing derivations for evolution."""
import sympy as s
from .identities import Identities, simp
from .components import *

def derive_short_time():
    """Constraint-compatible diagonal vacuum evolution and short-time controls."""
    from itertools import product
    ck = Identities()
    ms = s.symbols('m1:4', real=True)
    hs = s.symbols('h1:4', real=True)
    hds = s.symbols('hd1:4', real=True)
    theta = sum(hs)
    f = sum(ms) ** 2 / 2 - sum((v * v for v in ms))
    ri = [2 * v * v - sum(ms) * v + f for v in ms]
    C = coframe(s.diag(*ms), s.S.One, 0)
    for i in range(3):
        C[i + 1, 0, i + 1] = hs[i]
        C[i + 1, i + 1, 0] = -hs[i]
    D = lambda e: sum((s.diff(e, v) * dv for v, dv in list(zip(hs, hds)) + [(ms[i], ms[i] * (2 * hs[i] - theta)) for i in range(3)]))
    R = curvature(lc_connection(C), C, D)
    rc = ricci(R)
    sc = scalar(R)
    ck.zero('diagonal_time_Ricci', rc[0, 0] + sum(hds) + sum((v * v for v in hs)))
    ck.zero('diagonal_spatial_Ricci', [rc[i + 1, i + 1] - hds[i] - theta * hs[i] - ri[i] for i in range(3)])
    ck.zero('momentum_and_offdiagonal_closure', [rc[i, j] for i, j in product(range(4), repeat=2) if i != j])
    constraint = s.expand(f + theta ** 2 - sum((v * v for v in hs)))
    ck.zero('Hamiltonian_constraint_from_full_tensor', 2 * (rc[0, 0] + sc / 2) - constraint)
    kin = sum((v * v for v in hs)) - theta ** 2
    for i in range(3):
        df = sum((s.diff(f, ms[j]) * ms[j] * (2 * int(i == j) - 1) for j in range(3)))
        Euler = f + df + kin - 2 * (hds[i] - sum(hds)) - 2 * theta * (hs[i] - theta)
        ck.zero('action_Euler_' + str(i), Euler + 2 * (rc[i + 1, i + 1] - sc / 2))
    q, H, s1, s2 = s.symbols('q H s1 s2', real=True)
    ss = [s1, s2, -s1 - s2]
    X = [s.expand(ri[i] - f / 3) for i in range(3)]
    states = [q, H, s1, s2, *ms]
    rhs = [q * H, -3 * H ** 2 - f / 3, -3 * H * s1 - X[0], -3 * H * s2 - X[1]] + [ms[i] * (2 * ss[i] - H) for i in range(3)]
    Lie = lambda e: s.expand(sum((s.diff(e, v) * dv for v, dv in zip(states, rhs))))
    con = s.expand(f + 6 * H ** 2 - sum((v * v for v in ss)))
    ck.zero('constraint_propagation', Lie(con) + 6 * H * con)
    ck.zero('tracefree_shear_closure', sum(X))
    ck.zero('Raychaudhuri_on_constraint', Lie(q * H) + q * sum((v * v for v in ss)) / 3 + q * con / 3)
    q0 = s.Symbol('q0', positive=True)
    h0 = s.Symbol('H0', real=True)
    ns = s.symbols('n1:4', real=True)
    fn = sum(ns) ** 2 / 2 - sum((v * v for v in ns))
    chi = fn / 6
    Xn = [s.expand(2 * v * v - sum(ns) * v + 2 * fn / 3) for v in ns]
    J2 = sum((v * v for v in Xn))
    initial = {q: q0, H: h0, s1: 0, s2: 0, **dict(zip(ms, [v / q0 for v in ns]))}
    reduceH = lambda e: s.factor(s.rem(s.Poly(s.cancel(e), h0), s.Poly(h0 * h0 + chi / q0 ** 2, h0)).as_expr())

    def jet(expr, order):
        arr = []
        for k in range(order + 1):
            arr.append(reduceH(expr.subs(initial)))
            expr = Lie(expr)
        return arr
    qjet = jet(q, 4)
    betajets = [[s.S.Zero] + jet(si, 3) for si in ss]
    ck.zero('scale_first_three_derivatives', [qjet[0] - q0, qjet[1] - q0 * h0, qjet[2], qjet[3]])
    ck.zero('scale_fourth_derivative', qjet[4] + 2 * J2 / (3 * q0 ** 3))
    for i in range(3):
        ck.zero('shape_initial_acceleration_' + str(i), betajets[i][2] + Xn[i] / q0 ** 2)
        ck.zero('shape_cubic_derivative_' + str(i), betajets[i][3] - 5 * h0 * Xn[i] / q0 ** 2)
    ck.zero('initial_constraint', con.subs(initial).subs(h0 ** 2, -chi / q0 ** 2))
    import numpy as np
    from scipy.integrate import solve_ivp
    fnrhs = s.lambdify(states, rhs, 'numpy')
    fconstraint = s.lambdify(states, con, 'numpy')
    tau = s.Symbol('tau', real=True)
    rows = {}
    for name, nn in {'Nil': (0, 0, -1), 'Sol': (0, 1, -1)}.items():
        sn = dict(zip(ns, nn))
        ch = simp(chi.subs(sn))
        xx = [simp(v.subs(sn)) for v in Xn]
        jj = simp(J2.subs(sn))
        case = dict(n=nn, chi=ch, X=xx, J2=jj, branches={})
        for sign in (-1, 1):
            hv = sign * np.sqrt(float(-ch))
            y0 = np.array([1, hv, 0, 0, *nn, 0, 0], dtype=float)

            def fun(t, y):
                return np.r_[np.asarray(fnrhs(*y[:7]), dtype=float), y[2], y[3]]
            times = np.array([0.01, 0.02, 0.04])
            sol = solve_ivp(fun, (0, float(times[-1])), y0, method='DOP853', rtol=2e-13, atol=2e-15, t_eval=times, max_step=0.002)
            ck.ok(name + str(sign) + '_numerical_success', sol.success)
            records = []
            for tnum, yv in zip(times, sol.y.T):
                beta = np.r_[yv[7:9], -yv[7] - yv[8]]
                b3 = np.array([-float(v) * tnum ** 2 / 2 + 5 * hv * float(v) * tnum ** 3 / 6 for v in xx])
                q4 = 1 + hv * tnum - float(jj) * tnum ** 4 / 36
                errb = float(np.linalg.norm(beta - b3))
                errq = float(abs(yv[0] - q4))
                cval = float(fconstraint(*yv[:7]))
                ck.ok(name + str(sign) + str(tnum) + '_constraint', abs(cval) < 1e-11)
                records.append(dict(tau_over_q0=float(tnum), q=float(yv[0]), beta=beta.tolist(), shape_error_order3=errb, scale_error_order4=errq, constraint=cval, shape_scaled_remainder=errb / tnum ** 4, scale_scaled_remainder=errq / tnum ** 5))
            for key in ('shape_scaled_remainder', 'scale_scaled_remainder'):
                vals = [v[key] for v in records]
                ck.ok(name + str(sign) + '_' + key + '_bounded_at_checked_points', max(vals) < 0.5)
            case['branches'][str(sign)] = records
        rows[name] = case
    return dict(identities=ck.names, states=states, rhs=rhs, constraint=con, propagation='D constraint=-6H constraint', q_derivatives=qjet, beta_derivatives=betajets, short_time='u=tau/q0; beta=-X*u^2/2+5*(H0*q0)*X*u^3/6+O(u^4); q/q0=1+(H0*q0)*u-J2*u^4/36+O(u^5)', geometries=rows, numerical_control='q0=1; both H signs; DOP853; times .01,.02,.04; finite checks, not rigorous error bounds.', scope='Constraint-compatible local diagonal vacuum development. Common-scale curves are comparison curves, not vacuum backgrounds. No stability or long-time conclusion.')
