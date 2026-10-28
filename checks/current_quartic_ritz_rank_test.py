"""Outward arithmetic for the full-operator rank-one test; no finite-mode inverse."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
from math import factorial
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-quartic-ritz-rank-test'
RUNTIME = ROOT / 'evidence/current-coupled-mass-adjoint-runtime.json'
runtime = json.loads(RUNTIME.read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/' + STEM + '.md', 'checks/current_quartic_ritz_rank_test.py',
    'notes/current-finite-local-action-assembly.md',
    'notes/current-canonical-heat-restoration-assembly.md',
    'notes/current-geodesic-two-sided-remainder.md',
    'notes/stationary-curved-quantum-tt.md',
    'notes/current-spatial-inverse-output-test.md',
    'notes/current-quartic-resolvent-dual-test.md',
    'evidence/current-finite-local-action-assembly.json',
    'evidence/current-canonical-heat-restoration-assembly.json',
    'evidence/current-geodesic-two-sided-remainder.json',
    'evidence/current-spatial-inverse-output-test.json',
    'evidence/current-quartic-resolvent-dual-test-acceptance.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(v):
    d = v['decimal_enclosure']
    return arb(int(d['mid']), int(d['rad'])) * arb(10)**int(d['exponent'])


def enc(x):
    m, r, e = x.mid_rad_10exp(60)
    return {'ball': x.str(40), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def receipt(stem):
    obj = json.loads((ROOT / 'evidence' / (stem+'.json')).read_bytes())
    run = obj.get('runs', [obj])[-1]
    expected = ('error_estimator_test_passed_signed_action_unresolved'
                if stem == 'current-spatial-inverse-output-test' else 'passed')
    assert run['status'] == expected
    return run


def pick(run, kind):
    return next(x for x in run['numerics'] if x['kind'] == kind)


def add(*ps):
    out = [arb(0)] * max(map(len, ps))
    for p in ps:
        for j, c in enumerate(p):
            out[j] += c
    return out


def scale(p, c):
    return [c*x for x in p]


def mul(p, q):
    out = [arb(0)] * (len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out


def diff(p):
    return [j*p[j] for j in range(1, len(p))] or [arb(0)]


def shift(p, k):
    if k >= 0:
        return [arb(0)]*k+p
    assert all(x.is_zero() for x in p[:-k])
    return p[-k:]


def primitive_at(p, end):
    return sum((c*end**(j+1)/(j+1) for j, c in enumerate(p)), arb(0))


def read_poly(items):
    assert all(x['log_power'] == 0 and x['power'] >= 0 for x in items)
    out = [arb(0)]*(1+max(x['power'] for x in items))
    for x in items:
        out[x['power']] = dec(x['coefficient'])
    return out


def execute():
    ctx.prec = 512
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for pkg in runtime['packages'].values():
        for name, expected in pkg['files'].items():
            assert digest(Path(pkg['root']) / name) == expected

    local = receipt('current-finite-local-action-assembly')
    restoration = receipt('current-canonical-heat-restoration-assembly')
    receipt('current-geodesic-two-sided-remainder')
    sp = receipt('current-spatial-inverse-output-test')['result']['constants']
    endpoint = pick(local, 'actual_seed_correlated_endpoint')
    assert endpoint['exact_Q_bolt'] == '0'
    L = dec(endpoint['z_h'])
    Q = [arb(0)] * 7
    for k, value in zip([0, 2, 3, 4, 5, 6], endpoint['Q_coefficients']):
        Q[k] = dec(value)
    beta = Q[0]
    constants = pick(restoration, 'common_restoration_coefficients')
    C, B2 = dec(constants['C']), dec(constants['B2'])
    q = arb(2000)**4
    epsilon = dec(sp['epsilon'])
    eta0 = arb('.0022')
    qt = arb('.338')/(1-epsilon)
    ft = qt/(1+qt)
    assert arb('.99') < L < 1 and 0 < epsilon < arb('1.227e-74')

    mr = pick(local, 'signed_actual_local_matrix')
    assert mr['complete_scalar_and_mixed_rows']
    M = [[dec(x) for x in row] for row in mr['matrix']]
    A = [[M[i][j] for j in [0, 1]] for i in [0, 1]]
    H = [[M[i][j] for j in [2, 4]] for i in [0, 1]]
    detA = A[0][0]*A[1][1]-A[0][1]*A[1][0]
    assert A[0][0] > 0 and detA > 0
    Ai = [[A[1][1]/detA, -A[0][1]/detA],
          [-A[1][0]/detA, A[0][0]/detA]]
    B = [[sum((H[k][i]*Ai[k][l]*H[l][j]
               for k in range(2) for l in range(2)), arb(0))/beta
          for j in range(2)] for i in range(2)]
    detB = B[0][0]*B[1][1]-B[0][1]*B[1][0]
    traceB = B[0][0]+B[1][1]
    assert B[0][0] > 0 and detB > 0
    assert B[0][1].overlaps(B[1][0])
    low_eigen = detB/traceB

    known = pick(restoration, 'known_TT_matrix_with_shared_vacuum')
    assert known['basis'] == ['TT_J4', 'TT_J5', 'TT_zero_J']
    geo_pref = 2*(4*C*C/q)/(525*arb.pi()*B2*beta)
    S2_upper = 4*C*C/q
    # Positive monotone sum: first term plus its integral from 1 to infinity.
    S2_integral_test = 2*C*C*(q+2)/(q+1)**2 + 2*C*C/(q+1)
    assert S2_integral_test < S2_upper
    moment0 = Fraction(2**5*factorial(3)**4, factorial(7))
    moment1 = Fraction(2**5*factorial(3)**2*factorial(4)*factorial(2), factorial(7))
    assert moment0 == Fraction(288, 35) and moment1 == Fraction(384, 35)
    assert (moment1-moment0)/1440 == Fraction(1, 525)

    trials = {r['name']: r for r in local['numerics'] if r['kind'] == 'actual_global_trial'}
    energies = []
    for name, ki in [('TT_J4', 0), ('TT_zero_J', 2)]:
        v = read_poly(trials[name]['u'])
        if name == 'TT_J4':
            extra = [arb(0)]*6
            for k in range(2, len(Q)):
                extra[k-1] += k*Q[k]/(k-1)
            xi = scale(add(shift(Q, 3), shift(extra, 4)), 4)
            J = [arb(0)]*4+[arb(1)]
        else:
            xi, J = [arb(0)]*4+[arb(1)], [arb(0)]
        v_rebuilt = scale(add(xi, scale(diff(mul(Q, J)), arb(1)/2)), -1)
        assert all(x.contains(0) for x in add(v, scale(v_rebuilt, -1)))
        # Multiplied York differential identity, avoiding a division by Q.
        # -z*(Q^2 J')'+4 Q^2 J' = Q*(4xi-z xi').
        qqjp = mul(mul(Q, Q), diff(J))
        york = add(scale(shift(diff(qqjp), 1), -1), scale(qqjp, 4),
                   scale(mul(Q, add(scale(xi, 4), scale(shift(diff(xi), 1), -1))), -1))
        assert all(x.contains(0) for x in york)
        assert all(x.contains(0) for x in v[:3])
        av = mul(Q, diff(diff(diff(v))))
        dv = add(scale(mul(diff(Q), diff(diff(v))), arb(3)/2),
                 scale(mul(diff(diff(Q)), diff(v)), arb(1)/2))
        polynomial = mul(Q, add(scale(mul(av, av), arb(5)/8),
                                 scale(mul(av, dv), arb(3)/2), mul(dv, dv)))
        integral = primitive_at(polynomial, L)
        assert integral > 0
        egeo = geo_pref*integral
        ek = -dec(known['matrix'][ki][ki])/beta
        ei = ek+egeo
        assert ei > 0 and egeo < arb('1e-32')
        energies.append({'trial': name, 'flow_polynomial_integral': enc(integral),
                         'known_energy_diagonal': enc(ek), 'complete_geodesic_nonlocal_upper': enc(egeo),
                         'energy_diagonal_upper': enc(ei)})
    L_up = (1+epsilon)*sum((dec(e['energy_diagonal_upper']) for e in energies), arb(0))/(1-eta0)
    mu = low_eigen/L_up
    assert mu > 0
    flo = mu/(1+qt)
    residual_ideal_upper = qt*(qt+2)/(2*(1+qt).sqrt())
    W0, W1 = dec(sp['W0_complete_scalar_physical_upper']), dec(sp['W1_complete_mixed_physical_upper'])
    Gamma = lambda z: W0*(3*z+3*z*z+z**3)+W1*(2*z+z*z)
    bolt = arb('1.667e-9')
    rank_estimator = Gamma(mu/ft)

    # Same-domain alternative, with exact applications kept distinct from targets.
    ed = dec(sp['alternative_eta'])
    e0 = (eta0+epsilon)/(1-epsilon)
    r = (e0+ed)/(1-e0)
    de = ft*(ed+r**12)/(1-ed)
    rho_exact = (1+qt)*de
    S_bound = qt/(1-e0)
    A12 = lambda x: sum((x**k for k in range(12)), arb(0))
    D12 = lambda x: sum((k*x**(k-1) for k in range(1, 12)), arb(0))
    numerical_example = arb('5e-32')
    coeff_S = A12(r+numerical_example)
    coeff_K = S_bound*D12(r+numerical_example)
    aapp_example = (coeff_S+coeff_K)*numerical_example
    rho_example = (1+qt)*(de+aapp_example)
    residual_target = arb('1e-30')
    residual_room_exact = residual_target-rho_exact
    common_target_ceiling_linearized = residual_room_exact/((1+qt)*(A12(r)+S_bound*D12(r)))
    common_low, common_high = arb('5.8261192456e-32'), arb('5.8261192458e-32')
    common_report = lambda t: (1+qt)*(de+(A12(r+t)+S_bound*D12(r+t))*t)
    assert common_report(common_low) < residual_target
    assert common_report(common_high) > residual_target
    assert rho_example < residual_target and rho_exact < residual_target
    assert flo > arb('1e-5') and rank_estimator > bolt

    c = dict(q_t=qt, f_t=ft, epsilon=epsilon, beta=beta,
             complete_scalar_determinant=detA, restricted_B_determinant=detB,
             restricted_B_trace=traceB, restricted_B_lower_eigenvalue=low_eigen,
             full_angular_inverse_fourth_sum_upper=S2_upper,
             full_angular_integral_test_upper=S2_integral_test,
             complete_test_space_energy_upper=L_up,
             full_rank_one_residual_lower_bound=mu,
             full_rank_one_map_error_lower_bound=flo,
             rank_one_map_error_over_1e_minus30=flo/arb('1e-30'),
             ideal_exact_Ritz_residual_upper=residual_ideal_upper,
             rank_one_residual_to_bolt_upper_estimator_floor=rank_estimator,
             rank_one_bolt_estimator_floor_over_whole_allowance=rank_estimator/bolt,
             alternative_cutoff_eta=ed, alternative_contraction=r,
             alternative_exact_application_map_error_upper=de,
             alternative_exact_application_residual_upper=rho_exact,
             residual_room_after_exact_application_bound=residual_room_exact,
             common_application_threshold_linearized_NOT_attained=common_target_ceiling_linearized,
             common_application_threshold_certified_bracket_lower=common_low,
             common_application_threshold_certified_bracket_upper=common_high,
             proposed_common_application_error_NOT_attained=numerical_example,
             source_application_error_multiplier=coeff_S,
             K_application_error_multiplier=coeff_K,
             conditional_total_application_error=aapp_example,
             conditional_complete_residual_upper=rho_example,
             conditional_residual_screen_slack=residual_target-rho_example,
             conditional_source_application_residual_debit=(1+qt)*coeff_S*numerical_example,
             conditional_K_application_residual_debit=(1+qt)*coeff_K*numerical_example,
             alternative_output_error_exact_applications=Gamma(de/ft),
             conditional_complete_bolt_error_upper=Gamma((de+aapp_example)/ft),
             conditional_bolt_error_fraction_whole_allowance=Gamma((de+aapp_example)/ft)/bolt)
    return {'precision_bits': ctx.prec, 'constants': {k: enc(v) for k, v in c.items()},
            'restricted_scalar_matrix': [[enc(x) for x in row] for row in A],
            'restricted_mixed_matrix': [[enc(x) for x in row] for row in H],
            'variational_B_lower_matrix': [[enc(x) for x in row] for row in B],
            'full_quantum_energy_upper_rows': energies,
            'unknown': {'implemented_coefficient_error': None, 'source_operator_application_error': None,
                        'K_operator_application_error': None, 'actual_signed_G_t': None,
                        'other_coupled_entries': None, 'remaining_physical_margin': None},
            'scope': {'actual_full_operator_rank_one_residual_lower_bound': True,
                      'complete_quantum_energy_in_lower_bound': True,
                      'finite_quantum_mode_replacement': False,
                      'actual_map_error_lower_bound': True,
                      'physical_output_error_lower_bound': False,
                      'ideal_exact_Ritz_residual_upper_bound': True,
                      'implemented_complete_correction': False,
                      'numerical_alternative_application_errors_attained': False,
                      'geometric_admission_or_existence_or_stability': False,
                      'construction_obstruction': False, 'full_goal_complete': False}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/'+STEM+'.json')
    args = parser.parse_args()
    output = ROOT / args.output
    if output.exists():
        raise RuntimeError('Execution receipts are immutable; select a fresh output path.')
    protected = json.loads(BASELINE.read_bytes())
    verify_protected = lambda: all(digest(Path(p)) == h for p, h in protected.items())
    assert verify_protected()
    before = {p: digest(ROOT / p) for p in INPUTS}
    data = {'started_utc': datetime.now(timezone.utc).isoformat(),
            'previous_turn': 'no progress toward residual: requested email review',
            'input_sha256_before': before, 'python': runtime['python_sha256'],
            'ordinary_argument_utf8': (ROOT / INPUTS[0]).read_text(encoding='utf-8'),
            'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8')}
    try:
        data['result'] = execute()
        data['status'] = 'passed'
    except Exception as exc:
        data['status'] = 'failed'
        data['error'] = repr(exc)
        raise
    finally:
        data['input_sha256_after'] = {p: digest(ROOT / p) for p in INPUTS}
        data['inputs_unchanged'] = before == data['input_sha256_after']
        data['protected_paths_unchanged'] = verify_protected()
        data['protected_path_count'] = len(protected)
        data['finished_utc'] = datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
    assert data['inputs_unchanged'] and data['protected_paths_unchanged']
    for name, value in data['result']['constants'].items():
        print(name+': '+value['ball'])
    print('Receipt: '+str(output))


if __name__ == '__main__':
    main()
