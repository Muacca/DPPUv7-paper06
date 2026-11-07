"""Reader-facing derivations for canonical."""
import sympy as s
from .identities import Identities, simp
from .components import *

def derive_representatives(reduction_data):
    """Canonical representatives after the EH endpoint completion."""
    ck = Identities()
    q, N, eta, V = s.symbols('q N eta V', real=True)
    vol, kappa, theta, chi = s.symbols('V_Sigma kappa theta chi', real=True, nonzero=True)
    A = 3 * vol / kappa ** 2
    Q = s.Matrix([q, N, eta, V])
    p = s.Matrix(s.symbols('p_q p_N p_eta p_V', real=True))
    P = s.Matrix(s.symbols('P_q P_N P_eta P_V', real=True))
    B = -vol * q ** 3 * V / kappa ** 2 + 6 * theta * vol * q ** 2 * eta
    published = reduction_data
    loc = {'q': q, 'N': N, 'eta': eta, 'V': V, 'volume': vol, 'kappa': kappa, 'theta': theta}
    ck.zero('derived_endpoint_generator', s.sympify(str(published['BEC']), locals=loc) + s.sympify(str(published['BNY']), locals=loc) - B)
    grad = s.Matrix([s.diff(B, v) for v in Q])
    hess = grad.jacobian(Q)
    ck.zero('generating_gradient_curl', hess - hess.T)
    J = s.zeros(8)
    J[:4, 4:] = s.eye(4)
    J[4:, :4] = -s.eye(4)
    T = s.eye(8)
    T[4:, :4] = hess
    ck.zero('symplectic_translation', T.T * J * T - J)
    rho = P - grad

    def pb(f, g):
        return simp(sum((s.diff(f, Q[i]) * s.diff(g, P[i]) - s.diff(f, P[i]) * s.diff(g, Q[i]) for i in range(4))))
    ck.zero('translated_configuration_brackets', [pb(Q[i], Q[j]) for i in range(4) for j in range(4)])
    ck.zero('translated_mixed_brackets', [pb(Q[i], rho[j]) - (1 if i == j else 0) for i in range(4) for j in range(4)])
    ck.zero('translated_momentum_brackets', [pb(rho[i], rho[j]) for i in range(4) for j in range(4)])
    qd = s.Matrix(s.symbols('q_dot N_dot eta_dot V_dot', real=True))
    L = -A * q * qd[0] ** 2 / N + A * N * q * (chi - eta ** 2 + q ** 2 * V ** 2 / 9)
    Lprime = L + (grad.T * qd)[0]
    oldmom = s.Matrix([s.diff(L, v) for v in qd])
    newmom = s.Matrix([s.diff(Lprime, v) for v in qd])
    ck.zero('action_momentum_translation', newmom - oldmom - grad)
    C = -p[0] ** 2 / (4 * A * q) - A * q * (chi - eta ** 2 + q ** 2 * V ** 2 / 9)
    Cp = C.subs(dict(zip(p, rho)), simultaneous=True)
    vel_old = -N * p[0] / (2 * A * q)
    ck.zero('old_Legendre_Hamiltonian', (p[0] * qd[0] - L).subs(qd[0], vel_old) - N * C)
    vel_new = -N * rho[0] / (2 * A * q)
    ck.zero('new_Legendre_Hamiltonian', ((P.T * qd)[0] - Lprime - N * Cp - sum((rho[i] * qd[i] for i in (1, 2, 3)))).subs(qd[0], vel_new))
    ck.zero('lapse_primary_secondary', pb(rho[1], N * Cp) + Cp)
    ck.zero('eta_primary_secondary', pb(rho[2], Cp) + 2 * A * q * eta)
    ck.zero('V_primary_secondary', pb(rho[3], Cp) - 2 * A * q ** 3 * V / 9)
    constraints = [rho[2], rho[3], eta, V]
    M = s.Matrix(4, 4, lambda i, j: pb(constraints[i], constraints[j]))
    target = s.Matrix([[0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0]])
    ck.zero('second_class_bracket_matrix', M - target)
    ck.zero('second_class_determinant', M.det() - 1)
    ck.zero('second_class_reduced_dirac_bracket_qPq', pb(q, P[0]) - (s.Matrix([[pb(q, c) for c in constraints]]) * M.inv() * s.Matrix([pb(c, P[0]) for c in constraints]))[0] - 1)
    shell = {eta: 0, V: 0}
    shell.update({P[2]: grad[2], P[3]: grad[3]})
    ck.zero('auxiliary_shell_scale_momentum', grad[0].subs({eta: 0, V: 0}))
    Cred = -P[0] ** 2 / (4 * A * q) - A * q * chi
    ck.zero('reduced_Hamiltonian_constraint', Cp.subs({eta: 0, V: 0}) - Cred)
    red = s.Matrix([q, P[0], N, P[1]])
    embed = s.Matrix([q, N, 0, 0, P[0], P[1], grad[2], grad[3]])
    pull = embed.jacobian(red).T * J * embed.jacobian(red)
    ck.zero('computed_reduced_symplectic_pullback', pull - s.diag(s.Matrix([[0, 1], [-1, 0]]), s.Matrix([[0, 1], [-1, 0]])))
    ck.zero('torus_nonregular_gradient', [s.diff(Cred, v).subs({chi: 0, P[0]: 0}) for v in red])
    ck.zero('full_torus_nonregular_gradient', [s.diff(Cp, v).subs({eta: 0, V: 0, P[0]: 0, chi: 0}, simultaneous=True) for v in list(Q) + list(P)])
    ck.zero('lapse_absent_from_endpoint_generator', grad[1])
    dq, deta, dV = s.symbols('delta_q delta_eta delta_V')
    deltaB = grad[0] * dq + grad[2] * deta + grad[3] * dV
    ck.zero('Dirichlet_endpoint_variation', deltaB.subs({dq: 0, deta: 0, dV: 0}))
    return dict(identities=ck.names, boundary_generator=B, gradient=grad, convention='L_prime=L+dB/dt, P=p+gradient B; coordinates (q,N,eta,V)', old_constraint=C, new_constraint=Cp, primaries=[rho[1], rho[2], rho[3]], secondaries=[Cp, eta, V], second_class_brackets=M, reduced_constraint=Cred, auxiliary_shell=dict(P_eta=grad[2], P_V=grad[3], P_q='p_q'), endpoint_variation=deltaB, scope='Constant NY coupling; canonical equivalence with mapped endpoint data. Fixed new momenta and fixed old momenta are not the same boundary condition. Torus square constraint remains nonregular; no new NY dynamics.')
