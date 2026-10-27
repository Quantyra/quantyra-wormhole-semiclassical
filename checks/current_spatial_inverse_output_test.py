"""Original-output error propagation; does not evaluate the signed operator."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
STEM='current-spatial-inverse-output-test'
runtime=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0,str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS=[
    'notes/'+STEM+'.md','checks/current_spatial_inverse_output_test.py',
    'notes/current-spatial-remainder-regulator-test.md',
    'notes/current-joint-radius-composition-test.md',
    'notes/current-signed-adjoint-trial-attempt.md',
    'notes/current-complete-polynomial-weak-transfer.md',
    'notes/current-geodesic-two-sided-remainder.md',
    'notes/current-angular-resummed-quantum-operator.md',
    'checks/current_subtracted_thermal_action.py',
    'evidence/current-spatial-remainder-regulator-test.json',
    'evidence/current-joint-radius-composition-test.json',
    'evidence/current-signed-adjoint-trial-attempt.json',
    'evidence/current-complete-polynomial-weak-transfer.json',
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
    for pkg in runtime['packages'].values():
        for name,expected in pkg['files'].items():
            assert digest(Path(pkg['root'])/name)==expected
    spatial=json.loads((ROOT/'evidence/current-spatial-remainder-regulator-test.json').read_bytes())
    assert spatial['status']=='spatial_cutoff_arithmetic_passed_signed_evaluation_unresolved'
    sc=spatial['result']['constants']
    a,b=[dec(sc[k]) for k in ['eta_linear_coefficient','eta_quadratic_coefficient']]
    eps=dec(prior('current-complete-polynomial-weak-transfer')['constants']['epsilon'])
    joint=prior('current-joint-radius-composition-test')
    p_old=dec(joint['constants']['p_aff'])
    trial=json.loads((ROOT/'evidence/current-signed-adjoint-trial-attempt.json').read_bytes())['result']
    pref=dec(trial['known_complete_certificate_debits']['physical_prefactor_with_all_comparisons'])
    E32=dec(trial['known_complete_certificate_debits']['already_accepted_E32'])
    q=arb('.338'); f=q/(1+q); qt=q/(1-eps); ft=qt/(1+qt)
    deb=joint['complete_debits']
    oldW0=dec(deb['composed_L4_inverse_P_return'])+dec(deb['complete_composition_residual_return'])
    oldW1=dec(deb['complete_correlated_mixed_return'])
    W0=oldW0*(pref/p_old)*(ft/f)**3/(1-eps).sqrt()
    W1=oldW1*(pref/p_old)*(ft/f)**2/(1-eps).sqrt()
    sensitivity=3*W0+2*W1
    budget=arb('1.667e-9'); test_share=budget/100
    def eta(delta):return a*delta+b*delta**2
    def report(z):return W0*(3*z+3*z*z+z**3)+W1*(2*z+z*z)
    old_delta=arb('1e-22'); old_eta=eta(old_delta)
    assert old_eta.overlaps(dec(sc['eta_delta_original_energy_operator_upper']))
    old_error=report(old_eta/(1-old_eta))
    assert old_error>budget
    delta=arb('1e-38'); new_eta=eta(delta)
    assert 0<new_eta<1
    cutoff_error=report(new_eta/(1-new_eta))
    assert cutoff_error<test_share
    e0=(arb('.0022')+eps)/(1-eps)
    contraction=(e0+new_eta)/(1-e0)
    assert 0<contraction<1
    m=12
    z_exact_apps=(new_eta+contraction**m)/(1-new_eta)
    exact_apps_error=report(z_exact_apps)
    assert exact_apps_error<test_share
    app_target=arb('1e-30')
    conditional_error=report(z_exact_apps+app_target/ft)
    assert conditional_error<test_share
    assert E32<arb('6.131257147e-36')
    candidates={str(n):report((new_eta+contraction**n)/(1-new_eta)) for n in range(1,m+1)}
    first_whole=next(n for n in range(1,m+1) if candidates[str(n)]<budget)
    first_share=next(n for n in range(1,m+1) if candidates[str(n)]<test_share)
    vals={'original_physical_prefactor':pref,'unchanged_E32':E32,'epsilon':eps,
        'W0_complete_scalar_physical_upper':W0,'W1_complete_mixed_physical_upper':W1,
        'physical_error_sensitivity_S':sensitivity,
        'old_delta':old_delta,'old_eta':old_eta,
        'old_complete_inverse_output_error_upper':old_error,
        'old_estimator_overrun_factor':old_error/budget,
        'necessary_eta_ceiling_all_other_debits_zero':budget/sensitivity,
        'necessary_delta_ceiling_all_other_debits_zero':budget/(sensitivity*a),
        'alternative_delta':delta,'alternative_eta':new_eta,
        'alternative_cutoff_output_error_upper':cutoff_error,
        'geodesic_Et_comparison_error':e0,'Neumann_contraction_upper':contraction,
        'Neumann_relative_error_upper':contraction**m,
        'cutoff_plus_Neumann_error_exact_applications':exact_apps_error,
        'proposed_application_operator_error_NOT_attained':app_target,
        'conditional_output_error_if_application_target_met':conditional_error,
        'conditional_fraction_of_whole_bolt_allowance':conditional_error/budget}
    return {'constants':{k:enc(v) for k,v in vals.items()},'internal_Neumann_terms':m,
        'first_term_count_whole_allowance_exact_applications':first_whole,
        'first_term_count_one_percent_exact_applications':first_share,
        'term_count_diagnostic_errors_exact_applications':{k:enc(v) for k,v in candidates.items()},
        'actual_application_error':None,'actual_signed_G_norm':None,
        'signed_output_quadrature_error':None,'other_bolt_debits':None,
        'scope':{'paired_seed_form_inverse_defined':True,
            'old_cutoff_full_output_estimator_fails':True,
            'alternative_cutoff_and_inverse_truncation_error_bounded':True,
            'conditional_application_target_attained':False,
            'full_quantum_operator_numerically_applied':False,
            'complete_signed_certificate':False,'full_admission_or_stability':False,
            'construction_obstruction':False}}

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
        receipt['status']='error_estimator_test_passed_signed_action_unresolved'
    except BaseException as exc:
        receipt['status']='failed';receipt['failure']=repr(exc)
        raise
    finally:
        receipt['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        output.write_bytes((json.dumps(receipt,indent=2)+'\n').encode('utf-8'))
    assert before==receipt['input_sha256_after']
    print(json.dumps({'status':receipt['status'],
        'constants':{k:v['ball'] for k,v in receipt['result']['constants'].items()},
        'first_term_count_whole_allowance_exact_applications':receipt['result']['first_term_count_whole_allowance_exact_applications'],
        'first_term_count_one_percent_exact_applications':receipt['result']['first_term_count_one_percent_exact_applications']},indent=2))
