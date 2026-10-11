"""Energy extension with a reciprocal metric/charge vertex.

delta=L/q is motivated by the free LLL projection. The ground-state drift
source remains a prescribed closure; this is not a complete QED effective action.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def verify():
    r, t = sp.symbols('r t', real=True)
    L, s, delta = sp.symbols('L s delta', real=True)
    S, C, at, ar = [sp.Function(n)(r) for n in ['S', 'C', 'A_t', 'A_r']]
    F, Phi, X = [sp.Function(n)(t, r) for n in ['F', 'Phi', 'X']]
    k = L-2
    dot = lambda v: sp.diff(v, t)
    prime = lambda v: sp.diff(v, r)
    h0, h1 = S*prime(F)/k, dot(F)/(k*S)
    w0 = -2*h0/C-S*prime(Phi)/(L*at)
    w1 = -2*h1/C-dot(Phi)/(L*S*ar)
    G = prime(h0)-dot(h1)-prime(C)*h0/C
    E = dot(w1)-prime(w0)
    combination = F-2*Phi/L+delta*s*X/L
    eqf = sp.diff(F, t, 2)/(S*C)-prime(S*prime(F)/C)+k*combination/C**2
    eqp = sp.diff(Phi, t, 2)/(S*ar)-prime(S*prime(Phi)/at)+L*(Phi+s*X)/C-2*L*combination/C**2
    eqx = sp.diff(X, t, 2)/S-prime(S*prime(X))+(Phi+s*X)/C+delta*combination/C**2
    checks = {}
    def zero(name, value):
        assert sp.factor(value) == 0, name
        checks[name] = 'exact_zero'
    # Independent variation of the charge action fixes the reciprocal
    # source signs before imposing the Einstein-Maxwell constraints.
    hh0, hh1, aa0, aa1 = [sp.Function(n)(t, r) for n in ['h0', 'h1', 'w0', 'w1']]
    eeff = dot(aa1)-prime(aa0)+delta*(prime(hh0/C)-dot(hh1/C))
    lag = s*(dot(X)**2/(2*S)-S*prime(X)**2/2-X*eeff)
    def el(field):
        return sp.diff(lag, field)-dot(sp.diff(lag, dot(field)))-prime(sp.diff(lag, prime(field)))
    zero('action_Maxwell_t', el(aa0)+s*prime(X))
    zero('action_Maxwell_r', el(aa1)-s*dot(X))
    zero('action_metric_t', el(hh0)-delta*s*prime(X)/C)
    zero('action_metric_r', el(hh1)+delta*s*dot(X)/C)
    zero('action_anomaly', el(X)+s*(sp.diff(X, t, 2)/S-prime(S*prime(X))+eeff))
    zero('gravitational_definition', C*G+2*Phi/L-delta*s*X/L-F+C*C*eqf/k)
    zero('electric_definition', E-(Phi+s*X)/C+eqp/L+2*eqf/k)
    zero('anomaly_reconstruction', -sp.diff(X, t, 2)/S+prime(S*prime(X))-E-delta*G/C
         +eqx-eqp/L-(2+delta)*eqf/k)
    et = -S*prime(C*G)/(2*C)+k*h0/(2*C)+at*(w0+2*h0/C)/C+delta*s*S*prime(X)/(2*L*C)
    er = -dot(G)/(2*S)+k*h1/(2*C)+ar*(w1+2*h1/C)/C+delta*s*dot(X)/(2*L*C*S)
    zero('Einstein_t', et-S*prime(C*C*eqf)/(2*k*C))
    zero('Einstein_r', er-dot(eqf)*C/(2*k*S))
    zero('extra_axial_stress_divergence', dot(-delta*s*prime(X)/(2*L))+prime(delta*s*dot(X)/(2*L)))
    v = sp.Matrix(sp.symbols('f p x', real=True))
    potential = (L*v[0]-2*v[1]+delta*s*v[2])**2/C**2+L*(v[1]+s*v[2])**2/C
    V = sp.hessian(potential, tuple(v))/2
    fields = sp.Matrix([F, Phi, X])
    K = sp.diag(L*L/(k*S*C), 1/(S*ar), L*s/S)
    B = sp.diag(L*L*S/(k*C), S/at, L*s*S)
    system = K*fields.diff(t, 2)-(B*fields.diff(r)).diff(r)+V*fields
    for i, eq in enumerate([L*L*eqf/k, eqp, L*s*eqx]):
        zero('weighted_equation_'+str(i), system[i]-eq)
    density = (dot(fields).dot(K*dot(fields))+prime(fields).dot(B*prime(fields))+fields.dot(V*fields))/2
    flux = -dot(fields).dot(B*prime(fields))
    zero('local_energy_identity', dot(density)+prime(flux)-dot(fields).dot(system))
    null = sp.Matrix([-(2+delta)*s/L, -s, 1])
    for i in range(3):
        zero('constant_kernel_'+str(i), (V*null)[i])
    return {'identities': checks, 'potential_squares': '(L F-2 Phi+delta s X)^2/C^2 + L(Phi+s X)^2/C',
            'kinetic_and_gradient': 'Same positive diagonal matrices as delta=0',
            'physical_substitution': 'delta=L/q, s=s_bare*beta_l(q), for 2<=l<q',
            'closure': 'Prescribed ground-state drift plus reciprocal projected charge vertex',
            'microscopic_ground_state_response_proved': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = {'status': 'passed_conditional_reciprocal_vertex_energy',
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'sympy_version': sp.__version__, 'checks': verify(), 'lean_certified': False}
    args.output.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
