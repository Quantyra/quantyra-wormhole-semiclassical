"""Full first-inverse residual estimator; no signed quantum solve is claimed."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb, ctx, dec, enc
from flint import acb

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-first-inverse-source-screen'
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS = [
    'notes/'+STEM+'.md', 'checks/current_first_inverse_source_screen.py',
    'checks/current_quartic_ritz_rank_test.py',
    'notes/current-filtered-mixed-green-residual.md',
    'notes/current-complete-bolt-energy-bound.md',
    'notes/current-global-green-trial-enclosure.md',
    'notes/current-signed-adjoint-trial-attempt.md',
    'notes/current-quartic-ritz-rank-test.md',
    'evidence/current-global-green-trial-enclosure.json',
    'evidence/current-global-green-trial-enclosure-verification.json',
    'evidence/current-filtered-mixed-green-residual.json',
    'evidence/current-signed-adjoint-trial-attempt.json',
    'evidence/current-quartic-ritz-rank-test-final.json',
    'evidence/current-quartic-ritz-rank-test-verification.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def result(name, status):
    data = json.loads((ROOT/'evidence'/(name+'.json')).read_bytes())
    run = data.get('runs', [data])[-1]
    assert run['status'] == status
    return run['result']


def execute():
    ctx.prec = 512
    green = result('current-global-green-trial-enclosure',
                   'explicit_global_input_enclosed_signed_quantum_response_unresolved')
    p = {k: dec(v) for k, v in green['trial_parameters'].items()}
    f = result('current-filtered-mixed-green-residual',
               'filtered_input_bound_passed_complete_signed_response_unresolved')['constants']
    f = {k: dec(v) for k, v in f.items()}
    trial = result('current-signed-adjoint-trial-attempt', 'unresolved')
    rank = result('current-quartic-ritz-rank-test-final', 'passed')
    ell, nu, L = p['ell'], p['nu'], p['L']
    At = dec(green['constants']['A_t'])
    SP = arb('1e-16')
    assert nu > ell >= 1 and L > 2*SP and At > 0
    psi = lambda x: acb(x).digamma().real
    trigamma = lambda x: acb(x).polygamma(1).real
    cross = (psi(nu+1)-psi(ell+1))/(nu*(nu+1)-ell*(ell+1))
    J = trigamma(ell+1)/(2*ell+1)+trigamma(nu+1)/(2*nu+1)-2*cross
    x0 = p['a']*SP/p['gamma']
    tail = 1/(27*x0**3)
    assert x0 > 1 and J > tail > 0
    Hlower = At*(p['gamma']/p['a']*(J-tail)).sqrt()/L**2
    Hupper = At*(p['gamma']/p['a']*J).sqrt()/(L-2*SP)**2
    assert Hupper/Hlower-1 < arb('1e-14')
    q, epsilon = arb('.338'), f['epsilon']
    qt = q/(1-epsilon)
    ft = qt/(1+qt)
    Ks = f['K_V']*f['K_b_complete']/p['beta']
    Kg = Ks*q.sqrt()/(1-epsilon).sqrt()
    g_upper = Kg*Hupper
    vcap = dec(trial['known_complete_certificate_debits']['psi_retained_energy_upper'])
    coefficient = p['alpha']*q
    trial_norm = coefficient*vcap
    residual_upper = g_upper+trial_norm
    Kfirst = f['pref_t']*ft**2
    Kscalar = Kfirst*q.sqrt()/(1-epsilon).sqrt()
    B_lower = dec(rank['variational_B_lower_matrix'][1][1])
    psi_Hmu = (L**5/5).sqrt()
    Cfloor = B_lower/(vcap*psi_Hmu)
    report_floor = Kfirst*Cfloor*Hlower
    bolt = arb('1.667e-9')
    current_report = Kfirst*residual_upper
    baseline_report = Kfirst*g_upper
    conditional_C_ceiling = bolt/(Kfirst*Hlower)
    assert Cfloor > arb('.3853')
    assert Hlower > arb('3.3968e48') and Hupper < arb('3.3969e48')
    assert report_floor > bolt and current_report > bolt
    assert Cfloor > conditional_C_ceiling
    assert report_floor/bolt > arb('2e13')
    assert current_report/bolt > arb('1e14')
    assert current_report >= baseline_report
    c = dict(
        full_profile_square_integral=J, mixed_profile_integral=cross,
        entire_omitted_profile_tail_upper=tail,
        actual_trial_Hmu_norm_lower=Hlower, actual_trial_Hmu_norm_upper=Hupper,
        relative_Hmu_norm_bracket_width=Hupper/Hlower-1,
        scalar_source_uniform_upper=Ks, induced_source_uniform_upper=Kg,
        complete_source_dual_upper=g_upper, nonzero_trial_coefficient=coefficient,
        complete_trial_energy_upper=trial_norm, complete_first_residual_upper=residual_upper,
        first_residual_to_original_bolt_multiplier=Kfirst,
        first_scalar_residual_to_original_bolt_multiplier=Kscalar,
        complete_quartic_variational_B_lower=B_lower,
        complete_quartic_energy_upper=vcap,
        actual_uniform_source_operator_norm_lower=Cfloor,
        zero_trial_source_report=baseline_report,
        nonzero_trial_first_error_report=current_report,
        trial_addition_to_error_report=Kfirst*trial_norm,
        current_report_improvement_required=current_report/bolt,
        best_uniform_source_report_floor=report_floor,
        best_uniform_source_report_over_allowance=report_floor/bolt,
        conditional_uniform_source_constant_ceiling=conditional_C_ceiling,
        conditional_first_residual_ceiling_all_other_debits_zero=bolt/Kfirst,
        conditional_scalar_residual_ceiling_all_other_debits_zero=bolt/Kscalar,
        retained_input_debit=dec(green['constants']['complete_input_debit']),
        known_partial_subtotal_not_free_margin=dec(green['constants']['known_partial_subtotal_not_free_margin']),
    )
    assert all(x.is_finite() for x in c.values())
    return {
        'precision_bits': ctx.prec, 'constants': {k:enc(v) for k,v in c.items()},
        'scope': {
            'actual_input_norm_enclosed': True,
            'complete_first_residual_finite_upper': True,
            'uniform_source_estimator_class_rejected': True,
            'actual_residual_lower_bound': False, 'actual_G_lower_bound': False,
            'actual_scalar_source_evaluated': False, 'signed_quantum_action_evaluated': False,
            'complete_Zr_certificate': False, 'full_bolt_or_coupled_admission': False,
            'nonlinear_existence_or_matching_or_validity_or_stability': False,
            'construction_obstruction': False, 'full_objective_complete': False,
        },
        'unknowns': ['actual_scalar_source_and_residual', 'complete_signed_quantum_applications',
                     'Zr', 'subsequent_application_errors', 'e_other', 'other_coupled_entries'],
        'argument_review': {
            'type': 'author ordinary argument review; not independent or Lean verification',
            'checks': ['paired form-core source extension', 'Legendre Wronskian endpoint constants',
                       'cutoff and original z^-4 measure bounds', 'full quantum quartic energy',
                       'full scalar variational lower bound', 'estimator floor is not actual error',
                       'same physical prefactor and two retained F filters'],
        },
    }


if __name__ == '__main__':
    out = ROOT/'evidence'/(STEM+'.json')
    if out.exists():
        raise RuntimeError('Preserve immutable execution evidence; choose a fresh output.')
    before = {p:digest(ROOT/p) for p in INPUTS}
    protected = json.loads(BASELINE.read_bytes())
    assert all(digest(Path(p)) == h for p,h in protected.items())
    record = {'started_utc':datetime.now(timezone.utc).isoformat(),
              'input_sha256_before':before, 'python_sha256':digest(Path(sys.executable)),
              'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
              'arithmetic_source_utf8':Path(__file__).read_text(encoding='utf-8')}
    try:
        record['result'] = execute()
        after = {p:digest(ROOT/p) for p in INPUTS}
        assert after == before
        assert all(digest(Path(p)) == h for p,h in protected.items())
        record.update(status='uniform_source_estimator_rejected_full_signed_response_unresolved',
                      input_sha256_after=after, protected_paths_unchanged=len(protected))
    except Exception:
        import traceback
        record.update(status='failed', traceback=traceback.format_exc())
    record['finished_utc'] = datetime.now(timezone.utc).isoformat()
    out.write_text(json.dumps(record,indent=2)+'\n', encoding='utf-8')
    constants = record.get('result',{}).get('constants',{})
    keys = ['complete_first_residual_upper','nonzero_trial_first_error_report',
            'current_report_improvement_required','actual_uniform_source_operator_norm_lower',
            'best_uniform_source_report_floor','best_uniform_source_report_over_allowance',
            'conditional_scalar_residual_ceiling_all_other_debits_zero']
    print(json.dumps({'status':record['status'], 'error':record.get('traceback'),
                      'constants':{k:constants[k]['ball'] for k in keys if k in constants}},indent=2))
    if record['status']=='failed':
        raise SystemExit(1)
