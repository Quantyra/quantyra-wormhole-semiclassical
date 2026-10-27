"""Original-energy spatial-cutoff arithmetic, not signed action evaluation."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-spatial-remainder-regulator-test'
runtime = json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/'+STEM+'.md',
    'checks/current_spatial_remainder_regulator_test.py',
    'notes/stationary-curved-quantum-tt.md',
    'notes/current-geodesic-two-sided-remainder.md',
    'notes/current-seed-full-hessian-realization.md',
    'notes/current-coupled-truncation-energy-error.md',
    'notes/current-signed-adjoint-trial-attempt.md',
    'notes/current-signed-observation-regularization-test.md',
    'notes/current-signed-weak-source-transfer.md',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-complete-polynomial-weak-transfer.json',
    'evidence/current-signed-adjoint-trial-attempt.json',
    'evidence/current-adjoint-uniform-estimator-test.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(value):
    d=value['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])


def enc(value):
    m,r,e=value.mid_rad_10exp(60)
    return {'ball':value.str(35),'decimal_enclosure':
            {'mid':str(m),'rad':str(r),'exponent':str(e)}}


def prior(stem):
    run=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert run['status']=='passed'
    return run['result']


def execute():
    ctx.prec=512
    assert digest(Path(sys.executable))==runtime['python_sha256']
    for package in runtime['packages'].values():
        for name,expected in package['files'].items():
            assert digest(Path(package['root'])/name)==expected
    original=prior('current-coupled-truncation-energy-error')
    beta,C,B2=[dec(original['parameters'][k]) for k in ['beta','C','B2']]
    CL,CD=[dec(original[k]) for k in
           ['global_shape_L2_constant','global_shape_Laplacian_constant']]
    epsilon=dec(prior('current-complete-polynomial-weak-transfer')['constants']['epsilon'])
    assert 0<epsilon<arb('1.227e-74')
    K,R0,delta=arb('1e13'),arb(1)/400,arb('1e-22')
    assert 0<delta<R0
    xE=(arb.pi()*CL*CD).sqrt()
    zE=(3*arb.pi()*CD**2/4+5*arb.pi()*CL*CD).sqrt()
    measure_factor=2*K*C/(B2*beta*(1-epsilon))
    linear=measure_factor*xE*zE
    quadratic=measure_factor*zE**2/8
    eta=linear*delta+quadratic*delta**2
    uniform=prior('current-adjoint-uniform-estimator-test')
    old=uniform['all_trial_uniform_estimator']
    necessary=dec(old['necessary_new_uniform_eta_ceiling'])
    eta_old=dec(uniform['eta_available_upper'])
    assert eta<necessary/10
    trial=json.loads((ROOT/'evidence/current-signed-adjoint-trial-attempt.json').read_bytes())['result']
    assert trial['scope']['explicit_nonzero_paired_trial_constructed']
    ell_b=dec(trial['trial']['multipliers']['p3'])
    vcap=dec(trial['known_complete_certificate_debits']['psi_retained_energy_upper'])
    K3=arb('8.1e-45')
    physical_per_phi=K3*eta
    fixed=physical_per_phi*2*ell_b*vcap
    allowance=arb('1.667e-9')-arb('6.447981645e-11')
    assert fixed<allowance
    vals={'K_envelope':K,'anchor_R0':R0,'chosen_delta':delta,
        'C_L_original':CL,'C_D_original':CD,'x_E_flow_gradient':xE,
        'z_E_flow_second_derivative':zE,'epsilon':epsilon,
        'action_measure_and_envelope_factor':measure_factor,
        'eta_linear_coefficient':linear,'eta_quadratic_coefficient':quadratic,
        'eta_delta_original_energy_operator_upper':eta,
        'necessary_old_uniform_eta_ceiling_NOT_sufficient':necessary,
        'eta_delta_over_necessary_uniform_ceiling':eta/necessary,
        'eta_delta_over_known_full_remainder_upper':eta/eta_old,
        'retained_center_triangle_upper_NOT_actual_value':eta_old+eta,
        'unchanged_old_all_trial_physical_estimator_failure_ratio':
            dec(old['improved_ratio_to_full_bolt_margin']),
        'ell_b_original':ell_b,'psi_energy_upper':vcap,
        'physical_cutoff_error_per_phi_energy_unit':physical_per_phi,
        'fixed_first_two_functional_cutoff_error_upper':fixed,
        'conditional_phi_energy_ceiling_all_unknown_debits_zero':
            (allowance-fixed)/physical_per_phi}
    return {'constants':{k:enc(v) for k,v in vals.items()},
        'actual_phi_energy_norm':None,'signed_retained_kernel_value':None,
        'signed_local_contacts':None,'signed_Zr':None,
        'scope':{
            'spatial_near_form_norm_convergence_proved_in_note':True,
            'principal_log_and_matched_contacts_retained':True,
            'all_three_weak_residual_cutoff_errors_accounted_for':True,
            'original_heat_source_decay_certified':False,
            'complete_signed_action_evaluated':False,
            'complete_regulator_and_evaluation_certificate':False,
            'new_signed_physical_margin':False,
            'actual_kappa_m_or_original_defect_certified':False,
            'construction_obstruction':False,
            'nonlinear_existence_or_stability_proved':False}}


if __name__=='__main__':
    output=ROOT/'evidence'/(STEM+'.json')
    assert not output.exists(),'Preserve immutable receipts.'
    before={p:digest(ROOT/p) for p in INPUTS}
    receipt={'started_utc':datetime.now(timezone.utc).isoformat(),
        'input_sha256_before':before,'python':sys.executable,
        'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
        'arithmetic_source_utf8':Path(__file__).read_text(encoding='utf-8')}
    try:
        receipt['result']=execute()
        receipt['status']='spatial_cutoff_arithmetic_passed_signed_evaluation_unresolved'
    except BaseException as exc:
        receipt['status']='failed'
        receipt['failure']=repr(exc)
        raise
    finally:
        receipt['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    assert before==receipt['input_sha256_after']
    print(json.dumps({'status':receipt['status'],
        'constants':{k:v['ball'] for k,v in receipt['result']['constants'].items()}},indent=2))
