"""Record a concrete paired-domain trial; the signed quantum residual is unresolved.

Only known norm debits are evaluated here. There is deliberately no fake
quantum evaluator, selected-mode surrogate, or successful feasibility status.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys, json
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-signed-adjoint-trial-attempt'
OWN=['notes/'+STEM+'.md','checks/current_signed_adjoint_trial_attempt.py']
INPUTS=OWN+['evidence/'+STEM+'-review.json',
 'evidence/current-complete-polynomial-weak-transfer.json',
 'evidence/current-schur-bolt-functional.json',
 'evidence/current-localized-bolt-adjoint-sensitivity.json',
 'evidence/current-interpolated-complete-bolt-test.json',
 'evidence/current-coupled-mass-adjoint-runtime.json',
 'notes/current-coupled-adjoint-residual-audit.md',
 'notes/current-common-regulator-strong-limit-test.md',
 'notes/current-angular-resummed-quantum-operator.md',
 'notes/current-bolt-resummed-dirac-parametrix.md',
 'notes/current-physical-york-dual-assembly.md']

def sha(path):return sha256(path.read_bytes()).hexdigest()
def dec(x):
    if not isinstance(x,dict):return arb(x.strip('[]'))
    y=x['decimal_enclosure']
    return arb(y['mid']+' +/- '+y['rad'])*arb(10)**int(y['exponent'])
def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def receipt(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed'
    return r

def execute():
    ctx.prec=512
    review=json.loads((ROOT/'evidence'/(STEM+'-review.json')).read_bytes())
    assert all(sha(ROOT/p)==h for p,h in review['scientific_sha256'].items())
    rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
    assert sha(Path(sys.executable))==rt['python_sha256']
    for pkg in rt['packages'].values():
        for p,h in pkg['files'].items():assert sha(Path(pkg['root'])/p)==h
    transfer=receipt('current-complete-polynomial-weak-transfer')['result']['constants']
    shape=receipt('current-localized-bolt-adjoint-sensitivity')['result']['parameters']
    broad=next(x for x in receipt('current-schur-bolt-functional')['numerics']
        if x.get('kind')=='actual_Schur_trial_diagonal_and_source_threshold' and x.get('name')=='TT_zero_J')
    Lg=dec(broad['actual_l_diagonal_upper_over_action'])
    h=dec(shape['h']);beta=dec(shape['beta']);alpha=(1+h)/h
    q=arb(169)/500;f=q/(1+q)
    eps,dQ,dw,h3,E32,CJ,DG=[dec(transfer[k]) for k in
        ['epsilon','dQ','dw','h3','E32_complete','CJ','DG_complete']]
    assert eps>0 and eps<arb('.5')
    qt=q/(1-eps);vcap=((1+eps)*Lg/beta).sqrt()
    pref=h3*(1+dQ)*(arb('1.1e-10')+dw)*(1+eps).sqrt()
    ca=alpha*q*f*f*vcap;ctau=ca
    known=pref*(ca+ctau)+E32;budget=arb('1.667e-9')
    assert known<budget
    conditional=(budget-known)/pref
    # Fallback upper expression only: this does not evaluate the signed Z.
    Ctilde=(CJ+DG)/(1-eps).sqrt()
    Zfallback=Ctilde+qt*alpha*f*f*vcap
    fallback=pref*(ca+ctau+Zfallback)+E32
    assert fallback>budget
    point=receipt('current-interpolated-complete-bolt-test')['result']['constants']
    Knew=dec(point['complete_unfiltered_Knew'])
    amplitudes={'p1':alpha*f*f,'p2':alpha*f,'p3':alpha,'a':alpha*q*f*f}
    return {
      'outcome':'unresolved_calculation_no_new_signed_feasibility_evidence',
      'trial':{'TT_profile':'xi=z^4, J=0, radial_b=-z^4',
        'physical_components':'nu=xi/2, T_shape=-xi',
        'multipliers':{k:enc(v) for k,v in amplitudes.items()},
        'q_used_as_trial_parameter_not_measured_eigenvalue':True},
      'known_complete_certificate_debits':{k:enc(v) for k,v in {
        'Lg_actual_full_diagonal_upper':Lg,'psi_retained_energy_upper':vcap,
        'physical_prefactor_with_all_comparisons':pref,'trial_a_norm_upper':ca,
        'Riesz_residual_upper':ctau,'trial_a_physical_debit':pref*ca,
        'Riesz_physical_debit':pref*ctau,'already_accepted_E32':E32,
        'known_partial_total':known,'known_partial_fraction':known/budget}.items()},
      'conditional_targets_not_available_margins':{k:enc(v) for k,v in {
        'bolt_after_known_debits_if_all_others_zero':budget-known,
        'signed_Z_plus_eZ_if_all_others_zero':conditional,
        'direct_affine_full_residual_cap_if_whole_margin_free':budget/Knew}.items()},
      'reused_fallback_not_signed_evaluation':{k:enc(v) for k,v in {
        'Z_upper_expression':Zfallback,'complete_report':fallback,
        'overrun':fallback/budget}.items()},
      'unknowns':{k:None for k in ['retained_quantum_action_on_trial_and_Green_return',
        'signed_transferred_residual_Z','weak_action_core_limit_error',
        'retained_solve_quadrature_and_Riesz_errors_eZ','other_bolt_debits',
        'original_Y_defect','full_admission_comparison','actual_kappa_m','all_other_coupled_entries']},
      'first_missing_quantitative_implication':
        'A certified E_t* enclosure of the common-restored full quantum action on chi_N-alpha*z^4, then the correlated three-residual transfer. Density defines the form but supplies no numerical core-limit modulus or orbit-kernel evaluator.',
      'scope':{'explicit_nonzero_paired_trial_constructed':True,
        'retained_signed_quantum_response_numerically_evaluated':False,
        'full_residual_enclosed':False,'new_consequential_feasibility_bound':False,
        'requested_milestone_complete':False,'bolt_admission_closed':False,
        'full_admission_closed':False,'construction_obstruction':False,'goal_complete':False}}

if __name__=='__main__':
    before={p:sha(ROOT/p) for p in INPUTS}
    r={'started_utc':datetime.now(timezone.utc).isoformat(),
       'input_sha256_before':before,'frozen_sources_utf8':{p:(ROOT/p).read_text() for p in OWN}}
    r['result']=execute()
    r['input_sha256_after']={p:sha(ROOT/p) for p in INPUTS}
    assert r['input_sha256_after']==before
    r['status']='unresolved'
    r['arithmetic_execution']='completed; not a quantum-response evaluation or feasibility pass'
    r['finished_utc']=datetime.now(timezone.utc).isoformat()
    suffix='-replay' if '--replay' in sys.argv else ''
    path=ROOT/'evidence'/(STEM+suffix+'.json')
    assert not path.exists()
    path.write_bytes((json.dumps(r,indent=2)+'\n').encode())
    print(json.dumps({'status':r['status'],'scope':r['result']['scope'],
        'known_partial_debit':r['result']['known_complete_certificate_debits']['known_partial_total']['ball'],
        'conditional_Z_allowance':r['result']['conditional_targets_not_available_margins']['signed_Z_plus_eZ_if_all_others_zero']['ball'],
        'fallback_overrun':r['result']['reused_fallback_not_signed_evaluation']['overrun']['ball']},indent=2))
