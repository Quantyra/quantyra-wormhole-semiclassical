"""Actual projection-value budgets and independent normalization fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'cold-projection-value-upper-bound'
PREVIOUS = 'coercive-signed-moment-certificate'
OWN = ['checks/cold_projection_value_upper_bound.py', 'notes/'+STEM+'.md',
       'evidence/cold-projection-value-source-review.json']
ROWS, FIXTURES = [], []
PROTOCOL = {'precision_digits': [50, 80], 'identity_tolerance': '1e-28',
            'cold_direction_upper': '1e49', 'coefficient_sufficient': '1e-49',
            'residual_targets': ['2e-49', '2e-14', '2e21'],
            'terminal_target': '2e56', 'exact_terminal_target': '1e57',
            'E_lower': '1.125e35', 'actual_source_evaluated': False,
            'warp_cases': [dict(lam='1.7', B='1.3', period='1.4', h='.1', c='2', left='.04'),
                           dict(lam='3.2', B='2.1', period='2.3', h='.02', c='1.3', left='1')]}


def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    row = {'name': name, 'passed': bool(ok), 'detail': str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, value):
    value = S.simplify(S.sympify(value).doit())
    check(name, value == 0, value)


def prior():
    return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]


def paths():
    return sorted(set(prior()['input_sha256_after']) | set(OWN) | {'evidence/'+PREVIOUS+'.json'})


def provenance():
    p = prior()
    check('predecessor_passed', p['status'] == 'passed' and len(p['checks']) == 108 and all(c['passed'] for c in p['checks']))
    check('predecessor_immutable', p['input_sha256_before'] == p['input_sha256_after'])
    bad = [f for f, h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest() != h]
    check('predecessor_hashes', not bad, bad)
    review = json.loads((ROOT/OWN[2]).read_bytes())
    check('retained_inputs_pinned', all(f in p['input_sha256_after'] for f in review['retained_inputs']))
    check('scope_retained', review['actual_direction_upper_bound_proved'] and not review['actual_source_certificate_passed'] and not review['goal_complete'])


def budgets():
    R = S.Rational
    xi = R('30.1')*R('3e-19')+R('1e-8')
    check('seed_vertex_budget', xi < R('1.000001e-8'), xi)
    at = R('1.026')*xi+R('3e-7')*2520
    ass = R('4e-7')*xi+R('1.000002')*2520
    check('full_TT_channel', at < R('.000757'), at)
    check('full_scalar_channel', ass < 2521, ass)
    at, ass, sg = R('.000757'), R(2521), R('4e-17')
    a, b = [R('4e6'), R(4)], [R(1000), R('2e-5')]
    energy = sum((aa*sg+bb/2)*v/180 for aa, bb, v in zip(a, b, [at, ass]))
    check('elliptic_energy_budget', energy < R('.00225'), energy)
    outgoing = 15*sg/180*sum((R('1.001')*bb/13+3*aa)*v for aa, bb, v in zip(a, b, [at, ass]))
    check('outgoing_pointwise_budget', outgoing < R('2e-13'), outgoing)
    check('TT_value_factor', R(25)**2 < R(5)**2*9*3)
    check('TT_gradient_factor', R(2161)**2 < R(87)**2*R('.90')*80*3**2)
    zero('return_gradient_sum', R('.8')*420+73*25-2161)
    zero('TT_value_budget', 5*at-R('.003785'))
    zero('TT_gradient_budget', 87*at-R('.065859'))
    check('Liouville_connection_bound', R('2.24')+(1/R('.04')+R('3.25'))/2 < 17)
    def sobolev_square(value, gradient):
        return value**2+2*value*(gradient+17*value)
    check('elliptic_pointwise_budget', sobolev_square(R('.001125'), R('.00225')) < R('.0072')**2)
    check('returned_TT_pointwise_budget', sobolev_square(R('.003785'), R('.065859')) < R('.032')**2)
    check('full_adjoint_pointwise_budget', R('1.000001')+R('.0072')+R('.032')+R('2e-13') < R('1.04'))
    check('cutoff_support_geometry', R('5.98')**2*R('.04')**2 > R('2.1')**-4)
    check('projection_pointwise_budget', R('2.1')*R('1.04')/R('2.9') < R('.8'))
    exp_upper = sum(R(1, S.factorial(j)) for j in range(6))+R(1, 720)/(1-R(1, 7))
    check('elementary_exponential_upper', exp_upper < R(11, 4))
    check('cold_log_lower', R(11, 4)**69 < R(2000)**10/100)
    check('thermal_power_upper', 8**3 < 23**2)
    cold_upper = R('.8')/70*23*2000**15
    check('whole_cold_direction_upper', cold_upper < R('1e49'), cold_upper)
    E = R(PROTOCOL['E_lower'])
    subtotal = sum(R(v)/E**j for j, v in enumerate(PROTOCOL['residual_targets']))+R(PROTOCOL['terminal_target'])/E**3
    coeff = subtotal/25
    exact = R(PROTOCOL['exact_terminal_target'])/(25*E**3)
    check('inexact_sufficient_coefficient', coeff < R('2.71e-50'), coeff)
    check('inexact_cold_subtraction', R('1e49')*coeff < R('.271'))
    check('exact_sufficient_coefficient', exact < R('2.81e-50'), exact)
    check('exact_cold_subtraction', R('1e49')*exact < R('.281'))
    check('old_threshold_not_sufficient_via_new_upper', R('1e49')/R('3e43') > 1)
    check('old_three_step_terminal_requires_tightening', R('1e49')*R('1e63')/(25*E**3) > 1)
    zero('new_terminal_target_factor', R('1e63')/R('1e57')-10**6)
    FIXTURES.append({'fixture': 'exact_charged_budgets', 'cold_upper': str(cold_upper),
                     'inexact_coefficient': str(coeff), 'exact_coefficient': str(exact),
                     'actual_residual_evaluated': False})


def symbolic_identities():
    y = S.symbols('y', positive=True)
    lam, gamma, h, t, GG, alpha, nn, theta, vv = S.symbols('lambda gamma h t G alpha N theta v', positive=True)
    U = lam**R(1, 4)*S.sqrt(gamma)*t*GG**R(1, 4)
    phi = vv/U
    direct = S.sqrt(gamma)*theta*t**-4*phi/(alpha*lam**R(-1, 4)*nn)
    zero('physical_projection_normalization', S.sqrt(h)*direct-theta*t**-5*(GG/h**2)**R(-1, 4)*vv/(alpha*nn))
    tf, gf = S.Function('t')(y), S.Function('G')(y)
    Uf = lam**R(1, 4)*S.sqrt(gamma)*tf*gf**R(1, 4)
    zero('Liouville_log_derivative', S.diff(Uf, y)/Uf-S.diff(tf, y)/tf-S.diff(gf, y)/(4*gf))
    zero('physical_scaled_measure', S.sqrt(lam)*gamma*t**2*S.sqrt(GG)-U**2)
    # An independent vector polynomial checks the exact identity underlying (9).
    x, z = S.symbols('x z', real=True)
    vector = S.Matrix([x+x**2, 1-2*x+3*x**2])
    norm2 = (vector.T*vector)[0]
    zero('vector_norm_derivative', S.diff(norm2, x)-2*(vector.T*vector.diff(x))[0])
    zero('vector_interval_fundamental_theorem', norm2-norm2.subs(x, z)-S.integrate(S.diff(norm2, x), (x, z, x)))
    wc = y*(1+S.log(R(1, 100)/y))
    zero('cold_weight_monotonicity_derivative', S.diff(wc, y)-S.log(R(1, 100)/y))
    FIXTURES.append({'fixture': 'exact_scaling_and_vector_identity', 'actual_source': False})


R = S.Rational


def independent_warps():
    for case, data in enumerate(PROTOCOL['warp_cases']):
        levels = []
        for digits in PROTOCOL['precision_digits']:
            with mp.workdps(digits):
                lam, B, period, h, c, left = [mp.mpf(data[k]) for k in ['lam', 'B', 'period', 'h', 'c', 'left']]
                right = left+1
                gamma = 4*mp.pi*period*B**3
                t = lambda y: 1+h*mp.cosh(c*y)
                tp = lambda y: h*c*mp.sinh(c*y)
                U = lambda y: lam**mp.mpf('.25')*mp.sqrt(gamma)*t(y)*mp.sqrt(tp(y))
                v = lambda y: mp.matrix([mp.exp(-y)*(1+y), mp.sin(y)+y*y/10])
                f = lambda y: v(y)/U(y)
                norm2 = lambda w: sum(a*a for a in w)
                invt = lambda tt: mp.acosh((tt-1)/h)/c
                # Direct areal-volume integral vs independent Liouville-coordinate integral.
                phys = mp.sqrt(lam)*mp.quad(lambda tt: gamma*tt**2*norm2(f(invt(tt))), [t(left), t(right)])
                transformed = mp.quad(lambda y: norm2(v(y)), [left, right])
                err = abs(phys-transformed)/(1+abs(transformed))
                check(f'independent_volume_{case}_{digits}', err < mp.mpf('1e-28'), mp.nstr(err, digits))
                m = lambda y: tp(y)/t(y)+c/(2*mp.tanh(c*y))
                rho0, rho1 = left/mp.sqrt(lam), right/mp.sqrt(lam)
                radial = lambda rho: f(mp.sqrt(lam)*rho)
                physical_grad = mp.sqrt(lam)*mp.quad(lambda rho: U(mp.sqrt(lam)*rho)**2*norm2(mp.diff(radial, rho))/lam, [rho0, rho1])
                transformed_grad = mp.quad(lambda y: norm2(mp.diff(v, y)-m(y)*v(y)), [left, right])
                err = abs(physical_grad-transformed_grad)/(1+abs(transformed_grad))
                check(f'independent_gradient_{case}_{digits}', err < mp.mpf('1e-28'), mp.nstr(err, digits))
                omitted = mp.quad(lambda y: norm2(mp.diff(v, y)), [left, right])
                gap = abs(omitted-transformed_grad)
                check(f'omitted_connection_detected_{case}_{digits}', gap > mp.mpf('.01'), mp.nstr(gap, digits))
                point = (left+right)/2
                theta, alpha, nn = mp.mpf('.6'), mp.mpf('.995'), mp.mpf('3.2')
                direct = mp.sqrt(gamma)*theta*t(point)**-4*f(point)/(alpha*lam**mp.mpf('-.25')*nn)
                alternate = theta*t(point)**-5*(tp(point)**2/h**2)**mp.mpf('-.25')*v(point)/(alpha*nn*mp.sqrt(h))
                err = mp.norm(direct-alternate)/(1+mp.norm(direct))
                check(f'independent_projection_{case}_{digits}', err < mp.mpf('1e-28'), mp.nstr(err, digits))
                wrong = theta*t(point)**-4*v(point)/(alpha*nn)
                check(f'omitted_Liouville_detected_{case}_{digits}', mp.norm(wrong-direct) > mp.mpf('.001'))
                levels.append({'precision': digits, 'physical_value_norm_squared': mp.nstr(phys, digits),
                               'physical_gradient_norm_squared': mp.nstr(physical_grad, digits),
                               'projection_value': mp.nstr(mp.norm(direct), digits)})
        with mp.workdps(80):
            for key in ['physical_value_norm_squared', 'physical_gradient_norm_squared', 'projection_value']:
                a, b = mp.mpf(levels[0][key]), mp.mpf(levels[1][key])
                err = abs(a-b)/(1+abs(b))
                check(f'precision_{case}_{key}', err < mp.mpf('1e-28'), mp.nstr(err, 80))
        FIXTURES.append({'fixture': 'independent_smooth_warp', 'case': case, 'levels': levels,
                         'actual_charged_geometry': False, 'actual_source': False})


def main():
    sys.stdout.reconfigure(encoding='utf-8', newline='\n')
    sys.stderr.reconfigure(encoding='utf-8', newline='\n')
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before, error = hashes(), None
    print(json.dumps({'event': 'start', 'story': 'S167 actual cold projection value', 'input_count': len(ps)}), flush=True)
    try:
        provenance()
        budgets()
        symbolic_identities()
        independent_warps()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c['passed'] for c in ROWS)
    run = {'story': 'S167 actual cold projection value', 'created_utc': datetime.now(timezone.utc).isoformat(),
           'status': 'passed' if good else 'failed', 'exception': error,
           'python': sys.version, 'sympy': S.__version__, 'mpmath': mp.__version__,
           'protocol': PROTOCOL, 'checks': ROWS, 'fixtures': FIXTURES,
           'input_sha256_before': before, 'input_sha256_after': after, 'unchanged_inputs': before == after,
           'actual_direction_upper_bound_proved': good, 'conditional_sufficient_cold_value_certificate_proved': good,
           'actual_source_evaluated': False, 'actual_source_certificate_passed': False,
           'actual_signed_coefficient_bounded_to_cold_target': False, 'actual_signed_compatibility_zero_proved': False,
           'actual_weighted_graph_iterates_certified': False, 'cold_metric_response_proved': False,
           'nonlinear_neighborhood_preserved': False, 'nonlinear_solution_constructed': False,
           'physical_source_projection_changed': False, 'physical_stability_proved': False,
           'formal_verification': False, 'S167_complete': False, 'S168_complete': False, 'goal_complete': False}
    if not good:
        for kind, p in zip(['checker', 'note', 'review'], OWN):
            run['failed_'+kind+'_source_utf8'] = (ROOT/p).read_text(encoding='utf-8')
    dest = ROOT/'evidence'/(STEM+'.json')
    saved = json.loads(dest.read_bytes()) if dest.exists() else {'runs': []}
    saved['runs'].append(run)
    dest.write_text(json.dumps(saved, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status': run['status'], 'checks': len(ROWS), 'unchanged_inputs': before == after}), flush=True)
    return 0 if good else 1


if __name__ == '__main__':
    raise SystemExit(main())
