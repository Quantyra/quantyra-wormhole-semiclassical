"""Filtered original-output arithmetic, not a signed quantum evaluation."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-filtered-mixed-green-residual'
sys.path.insert(0, str(ROOT / 'checks'))
from current_quartic_ritz_rank_test import arb, ctx, dec, enc

INPUTS = [
    'notes/' + STEM + '.md', 'checks/current_filtered_mixed_green_residual.py',
    'checks/current_quartic_ritz_rank_test.py',
    'notes/current-mixed-green-input-domain.md',
    'notes/current-direct-triple-return-attempt.md',
    'notes/current-global-mixed-schur-preconditioner.md',
    'notes/current-complete-bolt-energy-bound.md',
    'notes/current-signed-observation-regularization-test.md',
    'notes/current-spatial-inverse-output-test.md',
    'evidence/current-mixed-green-input-domain.json',
    'evidence/current-mixed-green-input-domain-verification.json',
    'evidence/current-global-mixed-schur-preconditioner.json',
    'evidence/current-complete-bolt-energy-bound.json',
    'evidence/current-joint-radius-composition-test.json',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-spatial-inverse-output-test.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def run(stem):
    data = json.loads((ROOT / 'evidence' / (stem + '.json')).read_bytes())
    return data.get('runs', [data])[-1]


def accepted_result(stem):
    data = run(stem)
    assert data['status'] == 'passed', (stem, data['status'])
    return data['result']


def execute():
    ctx.prec = 512
    runtime = run('current-coupled-mass-adjoint-runtime')
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for package in runtime['packages'].values():
        for name, expected in package['files'].items():
            assert digest(Path(package['root']) / name) == expected, name
    verification = run('current-mixed-green-input-domain-verification')
    assert verification['status'] == 'ordinary_input_domain_argument_and_frozen_arithmetic_reviewed'
    for name, expected in verification['artifact_sha256'].items():
        assert digest(ROOT / name) == expected, name

    c = accepted_result('current-joint-radius-composition-test')['constants']
    b = accepted_result('current-complete-bolt-energy-bound')['constants']
    e = accepted_result('current-coupled-truncation-energy-error')['parameters']
    mixed = run('current-global-mixed-schur-preconditioner')
    assert mixed['status'] == 'passed'
    spatial = run('current-spatial-inverse-output-test')
    assert spatial['status'] == 'error_estimator_test_passed_signed_action_unresolved'
    a = run('current-mixed-green-input-domain')
    assert a['status'] == 'input_domain_and_budget_test_passed_complete_signed_estimate_unresolved'
    a = a['result']['constants']

    beta, d, gamma = dec(c['beta']), dec(c['d']), dec(e['gamma'])
    alpha, pref = dec(a['alpha_correlated']), dec(a['pref_t'])
    eps = dec(spatial['result']['constants']['epsilon'])
    q = arb('.338')
    qt = q / (1 - eps)
    ft = qt / (1 + qt)
    p0 = 12 * beta - arb('1.96')
    rm = arb('1.556e-23')
    KV, Kb = dec(b['KV']), dec(b['Kb'])
    assert 0 < Kb < 1 and 0 < rm < 1 and p0 > arb('.03999999995')
    # This factor is 1/beta because dvol/rho_action=|W_action| z^-4 dz.
    Kraw = pref * KV * Kb * q.sqrt() * ft**2 / (beta * (1 - eps).sqrt())
    K = Kraw / (p0 * (1 - rm))
    r = arb('1e-110')
    astar = gamma * d / 4
    assert 0 < astar < arb('.99')
    Db = 4 / (gamma * d).sqrt()
    Ddiff = 4 * r.sqrt() / (gamma * d)
    Hb = beta * alpha * Db / 2
    Dh = beta * alpha * Ddiff / 2
    Epoint = K * Dh
    E0 = K * (Dh + rm * Hb)
    E1 = K * (Dh + rm**2 * Hb)
    budget = arb('1.667e-9')
    share = budget / 100
    implementation_contract = arb('4e30')
    Eimplemented = E0 + K * (1 + rm) * implementation_contract
    cap = (share / K - Dh - rm * Hb) / (1 + rm)
    assert 0 < K < arb('2.162e-42')
    assert 0 < E0 < arb('7.415e-12') < share
    assert 0 < Epoint < arb('1.6e-23')
    assert implementation_contract < cap
    assert Eimplemented < share
    # Interval support for the conditional precision contract is not a solve.
    relative_contract = implementation_contract / Hb
    assert arb('1.8e-23') < relative_contract < arb('1.9e-23')

    baseline = json.loads(BASELINE.read_bytes())
    assert all(digest(Path(p)) == expected for p, expected in baseline.items())
    constants = {
        'beta': beta, 'd': d, 'gamma': gamma, 'alpha_original_radius_cost': alpha,
        'pref_t': pref, 'epsilon': eps, 'p0': p0, 'r_mix': rm, 'K_V': KV,
        'K_b_complete': Kb, 'raw_Hmu_to_filtered_physical_multiplier': Kraw,
        'preconditioned_residual_to_complete_physical_multiplier': K,
        'point_trace_averaging_length': astar, 'point_resolvent_bound': Db,
        'point_to_actual_averaged_resolvent_difference_bound': Ddiff,
        'h_b_norm_upper': Hb, 'h_r_minus_h_b_norm_upper': Dh,
        'point_replacement_physical_charge': Epoint,
        'exact_two_resolvent_input_charge': E0,
        'exact_input_fraction_of_whole_bolt_not_free_margin': E0 / budget,
        'fallback_one_extra_exact_Neumann_charge': E1,
        'conditional_one_percent_allocation': share,
        'conditional_sum_of_application_error_contract': implementation_contract,
        'conditional_full_input_charge_with_contract': Eimplemented,
        'conditional_application_error_ceiling': cap,
        'conditional_contract_relative_to_h_b_upper': relative_contract,
    }
    return {
        'precision_bits': ctx.prec, 'constants': {k: enc(v) for k, v in constants.items()},
        'unknowns': ['actual_e_D', 'actual_e_P', 'signed_Zr', 'complete_quantum_inverse_applications',
                     'quantum_application_errors', 'e_other', 'remaining_coupled_margins'],
        'scope': {
            'filtered_input_comparison_proved': True,
            'unfiltered_Hmu_to_Et_comparison_proved': False,
            'actual_global_Green_profile_evaluated': False,
            'complete_signed_response_evaluated': False,
            'bolt_certificate_closed': False, 'original_defect_closed': False,
            'graph_and_other_physical_outputs_closed': False, 'mass_column_closed': False,
            'nonlinear_admission_contraction': False, 'exact_compatibility': False,
            'physical_matching': False, 'semiclassical_validity_noise': False,
            'all_sector_stability': False, 'construction_obstruction': False,
            'full_objective_achieved': False,
        },
        'protected_paths_unchanged': len(baseline),
    }


if __name__ == '__main__':
    hashes = {name: digest(ROOT / name) for name in INPUTS}
    result = {
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256_before': hashes,
        'ordinary_argument_utf8': (ROOT / INPUTS[0]).read_text(encoding='utf-8'),
        'source_utf8': Path(__file__).read_text(encoding='utf-8'),
        'runtime_python_sha256': digest(Path(sys.executable)),
    }
    try:
        result['result'] = execute()
        result['input_sha256_after'] = {name: digest(ROOT / name) for name in INPUTS}
        assert result['input_sha256_after'] == hashes
        result['status'] = 'filtered_input_bound_passed_complete_signed_response_unresolved'
    except Exception:
        result['status'] = 'failed'
        result['traceback'] = traceback.format_exc()
    result['finished_utc'] = datetime.now(timezone.utc).isoformat()
    out = ROOT / 'evidence' / (STEM + '.json')
    data = json.loads(out.read_bytes()) if out.exists() else {'runs': []}
    data['runs'].append(result)
    out.write_bytes((json.dumps(data, indent=2) + '\n').encode('utf-8'))
    constants = result.get('result', {}).get('constants', {})
    display = ['preconditioned_residual_to_complete_physical_multiplier',
               'exact_two_resolvent_input_charge', 'point_replacement_physical_charge',
               'conditional_full_input_charge_with_contract',
               'conditional_application_error_ceiling']
    print(json.dumps({'status': result['status'], 'traceback': result.get('traceback'),
                      'constants': {k: constants[k]['ball'] for k in display if k in constants}}, indent=2))
    if result['status'] == 'failed':
        raise SystemExit(1)
