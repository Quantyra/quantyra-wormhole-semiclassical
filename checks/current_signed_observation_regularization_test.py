"""Complete-map observation error only; NOT a signed quantum evaluation."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-signed-observation-regularization-test'
OWN=['notes/'+STEM+'.md','checks/current_signed_observation_regularization_test.py']
INPUTS=OWN+['evidence/'+STEM+'-review.json',
 'evidence/current-complete-polynomial-weak-transfer.json',
 'evidence/current-joint-radius-composition-test.json',
 'evidence/current-complete-bolt-energy-bound.json',
 'evidence/current-coupled-truncation-energy-error.json',
 'evidence/current-localized-bolt-adjoint-sensitivity.json',
 'evidence/current-signed-adjoint-trial-attempt.json',
 'evidence/current-coupled-mass-adjoint-runtime.json',
 'notes/current-global-mixed-schur-preconditioner.md',
 'notes/current-joint-radius-composition-test.md',
 'notes/current-complete-bolt-energy-bound.md',
 'notes/current-complete-polynomial-weak-transfer.md',
 'notes/current-signed-adjoint-trial-attempt.md',
 'notes/current-common-regulator-strong-limit-test.md']

def sha(p):return sha256(p.read_bytes()).hexdigest()
def dec(x):
    if not isinstance(x,dict):return arb(x.strip('[]'))
    d=x['decimal_enclosure']
    return arb(d['mid']+' +/- '+d['rad'])*arb(10)**int(d['exponent'])
def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def receipt(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed'
    return r['result']

def execute():
    ctx.prec=512
    review=json.loads((ROOT/'evidence'/(STEM+'-review.json')).read_bytes())
    assert all(sha(ROOT/p)==h for p,h in review['scientific_sha256'].items())
    rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
    assert sha(Path(sys.executable))==rt['python_sha256']
    for pkg in rt['packages'].values():
        for p,h in pkg['files'].items():assert sha(Path(pkg['root'])/p)==h
    transfer=receipt('current-complete-polynomial-weak-transfer')['constants']
    composition=receipt('current-joint-radius-composition-test')['constants']
    mixed=receipt('current-complete-bolt-energy-bound')['constants']
    energy=receipt('current-coupled-truncation-energy-error')['parameters']
    local=receipt('current-localized-bolt-adjoint-sensitivity')['parameters']
    trial=json.loads((ROOT/'evidence/current-signed-adjoint-trial-attempt.json').read_bytes())
    assert trial['status']=='unresolved'
    ledger=trial['result']['known_complete_certificate_debits']
    eps=dec(transfer['epsilon']);gamma=dec(energy['gamma']);h=dec(local['h'])
    beta=dec(composition['beta']);d=dec(composition['d'])
    ellb=(1+h)/h;q=arb('.338');qt=q/(1-eps);ft=qt/(1+qt)
    U2=dec(composition['U2']);R0=dec(composition['R0_complete_residual'])
    S0=1/arb('23.99').sqrt();S2=(beta/(6*d)).sqrt()
    p0=12*beta-(arb(7)/5)**2;rmix=arb('1.556e-23')
    assert p0>arb('.03999999995') and rmix<1 and eps<arb('.5')
    D0=U2+S2*S0*R0/beta
    D1=dec(mixed['KV'])*dec(mixed['Kb'])*q.sqrt()/(2*d*p0*(1-rmix))
    Dt=(D0*ft**3+D1*ft**2)/(1-eps).sqrt()
    Kreg=4*ellb*Dt/gamma
    pref=dec(ledger['physical_prefactor_with_all_comparisons'])
    fixed_r=arb('1e-110');budget=arb('1.667e-9');allocation=budget/100
    Ereg=pref*Kreg*fixed_r.sqrt()
    r_cap=(allocation/(pref*Kreg))**2
    assert fixed_r<arb('.01') and 10*fixed_r<gamma
    # Recover the actual coefficients with the correlated seed identities.
    C=dec(energy['C']);u=1/C;eps_seed=1-u-6*beta;vh=1+h
    tau=vh*h*h*(beta*(1+2/vh+3/vh**2)+u/vh**2+eps_seed/vh**4)
    qpp_guard=60*abs(beta)+26*u+82*abs(eps_seed)+6*abs(tau)
    assert qpp_guard<20
    rho_lo=2*(fixed_r/(2*gamma+10*fixed_r)).sqrt()
    rho_hi=2*(fixed_r/(2*gamma-10*fixed_r)).sqrt()
    known=dec(ledger['known_partial_total'])+Ereg
    conditional=(budget-known)/pref
    return {'outcome':'unresolved_signed_response_with_quantified_observation_replacement',
      'constants':{k:enc(v) for k,v in {
        'D0_complete_scalar_output_Laplacian':D0,
        'D1_complete_mixed_output_Laplacian':D1,
        'Dt_complete_retained_output_Laplacian':Dt,
        'Kreg_original_bolt_map':Kreg,'pref_t':pref,'fixed_averaging_r':fixed_r,
        'Ereg_original_P32_bolt_debit':Ereg,'fraction_of_whole_bolt_allowance':Ereg/budget,
        'proposed_one_percent_allocation':allocation,'estimator_r_cap_for_one_percent':r_cap,
        'proper_support_width_lower':rho_lo,'proper_support_width_upper':rho_hi,
        'global_Q_second_derivative_guard':qpp_guard,
        'known_partial_total_with_observation_error':known,
        'conditional_Zr_plus_error_if_other_debits_zero':conditional}.items()},
      'chosen_r_fits_proposed_allocation':bool(Ereg<allocation),
      'unknowns':{k:None for k in ['signed_transfer_Zr','retained_common_quantum_action',
        'UV_restoration_and_quadrature_error','inverse_and_Riesz_error',
        'discretization_and_arithmetic_error','other_bolt_debits','all_other_admission_outputs',
        'original_Y_defect','actual_mass_column','other_coupled_columns']},
      'scope':{'observation_difference_bounded_for_complete_structured_map':True,
        'physical_norm_changed':False,'P32_transport_refined':False,
        'raw_heat_cutoff_strong_convergence_proved':False,
        'signed_quantum_response_evaluated':False,'requested_milestone_complete':False,
        'bolt_allowance_closed':False,'full_admission_closed':False,
        'construction_obstruction':False,'goal_complete':False}}

if __name__=='__main__':
    before={p:sha(ROOT/p) for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),
      'input_sha256_before':before,
      'frozen_sources_utf8':{p:(ROOT/p).read_text() for p in OWN}}
    run['result']=execute()
    run['input_sha256_after']={p:sha(ROOT/p) for p in INPUTS}
    assert before==run['input_sha256_after']
    run['status']='unresolved'
    run['arithmetic_execution']='completed_supporting_observation_test_only'
    run['finished_utc']=datetime.now(timezone.utc).isoformat()
    suffix='-replay' if '--replay' in sys.argv else ''
    path=ROOT/'evidence'/(STEM+suffix+'.json');assert not path.exists()
    path.write_bytes((json.dumps(run,indent=2)+'\n').encode())
    print(json.dumps({'status':run['status'],
      'constants':{k:v['ball'] for k,v in run['result']['constants'].items()},
      'chosen_r_fits_proposed_allocation':run['result']['chosen_r_fits_proposed_allocation'],
      'scope':run['result']['scope']},indent=2))
