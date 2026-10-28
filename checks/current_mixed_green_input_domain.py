"""Actual mixed-Green input audit; does not evaluate the remaining quantum action."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-mixed-green-input-domain'
sys.path.insert(0, str(ROOT / 'checks'))
from current_quartic_ritz_rank_test import arb, ctx, dec, enc, receipt, pick
sys.path.insert(0, str(ROOT))
from operators.signed_local_core import SignedLocalCore
import sympy as sp

INPUTS = [
    'notes/' + STEM + '.md', 'checks/current_mixed_green_input_domain.py',
    'operators/signed_local_core.py', 'checks/current_quartic_ritz_rank_test.py',
    'notes/current-signed-adjoint-trial-attempt.md',
    'notes/current-signed-observation-regularization-test.md',
    'notes/current-signed-strong-local-core.md',
    'notes/current-finite-local-action-assembly.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/current-structured-weak-response-attempt.md',
    'notes/current-spatial-remainder-regulator-test.md',
    'notes/current-quartic-resolvent-dual-test.md',
    'evidence/current-finite-local-action-assembly.json',
    'evidence/current-canonical-heat-restoration-assembly.json',
    'evidence/current-signed-adjoint-trial-attempt.json',
    'evidence/current-spatial-inverse-output-test.json',
    'evidence/current-quartic-resolvent-dual-test-acceptance.json',
    'evidence/current-quartic-resolvent-dual-test-verification.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def last_run(stem):
    data = json.loads((ROOT / 'evidence' / (stem + '.json')).read_bytes())
    return data.get('runs', [data])[-1]


def execute():
    ctx.prec = 512
    core = SignedLocalCore()
    z = core.z
    k = sp.Symbol('k')
    trace_rows = core.independent_trace_variation()
    qsub = core.seed_substitution()
    rows = {
        name: [sp.expand(z**4 * v.subs(qsub)) for v in row[1]]
        for name, row in trace_rows.items()
    }
    coefficients = [sp.expand(sum(row[j] for row in rows.values())) for j in range(5)]
    assert all(a.is_polynomial(z) and not a.has(sp.log(z)) for a in coefficients)
    assert all(a.coeff(z, 1) == 0 for a in coefficients)
    Q = qsub[core.Q[0]]
    Dr = 96 * sp.pi * core.cR - 1 / (60 * sp.pi)
    assert sp.cancel(coefficients[4] - 2 * core.C * Dr * Q**2 / core.B2) == 0

    contributions = {
        name: sp.factor(sum(a.coeff(z, 0) * k**j for j, a in enumerate(row)))
        for name, row in rows.items()
    }
    lam = core.C * core.beta / core.B2
    physical_lambda0 = lam + 11 * lam**2 / (360 * sp.pi)
    normal = sp.expand(sum(contributions.values()).subs(core.lam0, physical_lambda0))
    expected = 2 * core.beta * (k - 4) * (k + 3) * (
        Dr * lam * k * (k - 3) - (1 + 11 * lam / (180 * sp.pi)))
    assert sp.cancel(normal - expected) == 0
    assert sp.cancel(normal.subs(k, 4)) == 0
    assert sp.cancel(normal.subs(k, -3)) == 0
    print('Exact complete mixed row, curvature cancellation and first recurrence: verified', flush=True)

    local = receipt('current-finite-local-action-assembly')
    restored = receipt('current-canonical-heat-restoration-assembly')
    endpoint = pick(local, 'actual_seed_correlated_endpoint')
    params = pick(restored, 'common_restoration_coefficients')
    beta = dec(endpoint['Q_coefficients'][0])
    C, B2 = dec(params['C']), dec(params['B2'])
    lam_num = C * beta / B2
    dlam = (96 * arb.pi() * arb('1e22') - 1 / (60 * arb.pi())) * lam_num
    cv = 1 + 11 * lam_num / (180 * arb.pi())
    sigp = (3 + (9 + 4 * cv / dlam).sqrt()) / 2
    sigm = 3 - sigp
    assert arb('18427.85') < sigp < arb('18427.86')
    assert arb('-18424.86') < sigm < arb('-18424.85')
    # This isolates every noninteger root separation; 4 - (-3) = 7 is retained.
    assert arb('18423.85') < sigp - 4 < arb('18423.86')
    assert arb('18430.85') < sigp + 3 < arb('18430.86')
    assert arb('18428.85') < 4 - sigm < arb('18428.86')
    assert arb('18421.85') < -3 - sigm < arb('18421.86')
    assert arb('36852.70') < sigp - sigm < arb('36852.72')
    p5 = 16 * beta * (10 * dlam - cv)
    dp4 = 14 * beta * (4 * dlam - cv)
    assert p5 < 0 and dp4 < 0

    # Preserve the ordinary audit attached to the accepted complete quartic column.
    verification = last_run('current-quartic-resolvent-dual-test-verification')
    assert verification['status'] == 'ordinary_column_argument_audited_and_numerical_enclosure_verified'
    for name, expected_hash in verification['artifact_sha256'].items():
        assert digest(ROOT / name) == expected_hash, name
    column = last_run('current-quartic-resolvent-dual-test-acceptance')
    assert column['status'] == 'complete_column_arithmetic_executed_requires_argument_audit'
    mpsi = dec(column['result']['constants']['zero_centered_massive_column_bound'])
    assert mpsi < arb('5.859156e-32')
    trial = last_run('current-signed-adjoint-trial-attempt')['result']['trial']['multipliers']
    p1, p2, p3 = (dec(trial[name]) for name in ['p1', 'p2', 'p3'])
    spatial = receipt('current-spatial-inverse-output-test')['result']['constants']
    eps = dec(spatial['epsilon'])
    pref = dec(spatial['original_physical_prefactor'])
    qt = arb('.338') / (1 - eps)
    ft = qt / (1 + qt)
    epsi = pref * (ft * p1 + ft**2 * p2 + ft**3 * p3) * mpsi
    kinput = pref * ft**3 * (1 + qt)
    assert 0 < epsi < arb('2.811e-42')
    assert epsi / arb('1.667e-9') < arb('1.687e-33')
    assert 0 < kinput < arb('1.082e-44')

    baseline = json.loads(BASELINE.read_bytes())
    assert all(digest(Path(p)) == expected_hash for p, expected_hash in baseline.items())
    return {
        'precision_bits': ctx.prec,
        'indicial_contributions': {name: str(value) for name, value in contributions.items()},
        'complete_mixed_euler_coefficients': [str(a) for a in coefficients],
        'exact_normal_factor': str(expected),
        'exact_first_coefficient_all_zero': True,
        'exact_resonance_retained': '4 - (-3) = 7; no blanket nonresonance assumption',
        'constants': {name: enc(value) for name, value in {
            'lambda': lam_num, 'd_lambda': dlam, 'c_v': cv,
            'sigma_plus': sigp, 'sigma_minus': sigm,
            'P_C_5': p5, 'P_C_derivative_4': dp4,
            'accepted_massive_quartic_column_bound': mpsi,
            'alpha_correlated': p3, 'pref_t': pref,
            'quartic_all_three_physical_charge_upper': epsi,
            'quartic_fraction_of_whole_bolt_not_free_margin': epsi / arb('1.667e-9'),
            'full_input_error_physical_multiplier': kinput,
        }.items()},
        'ordinary_source': {
            'url': 'https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode.pdf',
            'sections_read': '4.3-4.4; Theorems 4.11-4.13; equations (4.115)-(4.117)',
            'scope': 'convergence/classification of actual analytic Euler ODE; no quantum or global numerical estimate',
        },
        'unknowns': ['a_4', 'a_plus', 'numerical_input_C5_norm', 'e_chi',
                     'retained_signed_quantum_action', 'complete_Zr', 'complete_application_errors',
                     'remaining_coupled_margins'],
        'scope': {
            'actual_mixed_green_outer_input_domain_test': True,
            'affine_h1_log_result_preserved': True,
            'quartic_contribution_to_complete_three_residual_certificate_bounded': True,
            'complete_quantum_evaluation': False, 'bolt_certificate_closed': False,
            'original_strong_source_and_graph_maps': False,
            'nonlinear_admission_contraction': False, 'exact_compatibility': False,
            'physical_matching': False, 'semiclassical_validity_noise': False,
            'all_sector_stability': False, 'construction_obstruction': False,
            'full_objective_achieved': False,
        },
        'protected_paths_unchanged': len(baseline),
    }


if __name__ == '__main__':
    hashes = {name: digest(ROOT / name) for name in INPUTS}
    run = {
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256_before': hashes,
        'ordinary_argument_utf8': (ROOT / INPUTS[0]).read_text(encoding='utf-8'),
        'source_utf8': Path(__file__).read_text(encoding='utf-8'),
        'runtime_python_sha256': digest(Path(sys.executable)),
    }
    try:
        run['result'] = execute()
        run['input_sha256_after'] = {name: digest(ROOT / name) for name in INPUTS}
        assert run['input_sha256_after'] == hashes
        run['status'] = 'input_domain_and_budget_test_passed_complete_signed_estimate_unresolved'
    except Exception:
        run['status'] = 'failed'
        run['traceback'] = traceback.format_exc()
    run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    out = ROOT / 'evidence' / (STEM + '.json')
    data = json.loads(out.read_bytes()) if out.exists() else {'runs': []}
    data['runs'].append(run)
    out.write_bytes((json.dumps(data, indent=2) + '\n').encode('utf-8'))
    print(json.dumps({'status': run['status'], 'traceback': run.get('traceback'),
                      'constants': run.get('result', {}).get('constants')}, indent=2), flush=True)
    if run['status'] == 'failed':
        raise SystemExit(1)
