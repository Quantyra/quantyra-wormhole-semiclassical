"""Exact algebra and Arb support; does not apply any physical inverse."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb, ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-schur-polynomial-return-test'
OWN=['notes/'+STEM+'.md', 'checks/current_schur_polynomial_return_test.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-schur-form-preconditioner-test.md',
 'notes/current-schur-residual-phase-domain.md',
 'notes/current-schur-physical-output-test.md',
 'evidence/current-schur-physical-output-test.json',
 'notes/current-bolt-form-commutator-test.md',
 'evidence/current-bolt-form-commutator-test.json',
 'notes/current-seed-full-hessian-realization.md']
CHECKS=[]

def check(name,ok):
    CHECKS.append({'name':name,'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'),flush=True)
    if not ok: raise AssertionError(name)

def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}

def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    for p,h in review['scientific_sha256'].items():
        check('frozen_'+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)
    for p in (INPUTS[6],INPUTS[8]):
        check('prior_passed_'+p,json.loads((ROOT/p).read_bytes())['runs'][-1]['status']=='passed')
    t,x=sp.symbols('t x',real=True)
    q=sp.Rational(169,500); N=32; m=N-3
    ch0,ch1=sp.Poly(1,x),sp.Poly(x,x)
    for j in range(1,m):ch0,ch1=ch1,sp.Poly(2*x*ch1.as_expr()-ch0.as_expr(),x)
    check('exact_chebyshev_recurrence',ch1==sp.Poly(sp.chebyshevt(m,x),x))
    den=ch1.eval(-1-2/q)
    rp=sp.Poly(-t**3*ch1.as_expr().subs(x,2*t/q-1)/den,t)
    pip,rem=sp.div(sp.Poly(1-rp.as_expr(),t),sp.Poly(1+t,t))
    check('residual_normalization_at_minus_one',rp.eval(-1)==1)
    check('candidate_is_exact_degree_31_polynomial',rem.is_zero and pip.degree()==31)
    check('candidate_exact_residual',sp.Poly(1-(1+t)*pip.as_expr()-rp.as_expr(),t).is_zero)
    check('three_required_residual_factors',all(rp.nth(j)==0 for j in range(3)) and rp.nth(3)!=0)
    check('whole_interval_residual_envelope_attained_at_q',abs(rp.eval(q))==q**3/abs(den))
    ff=t/(1+t)
    ee=-t
    for p in (2,3):
        lhs=ee**p
        rhs=(-1)**p*ff**p*(1-ff)**(-p-1)/(1+t)
        check('signed_full_power_identity_p'+str(p),sp.cancel(lhs-rhs)==0)
    check('odd_power_previous_unsigned_identity_has_wrong_sign',sp.cancel(ee**3-ff**3*(1-ff)**(-4)/(1+t))!=0)
    err=rp.as_expr()/(1+t)
    fact=-ff**3*(1+t)**2*ch1.as_expr().subs(x,2*t/q-1)/den
    check('complete_physical_error_factorization',sp.cancel(err-fact)==0)
    deffact=-ff**3*(1+t)**3*ch1.as_expr().subs(x,2*t/q-1)/den
    check('original_source_defect_factorization',sp.cancel(rp.as_expr()-deffact)==0)
    qn=arb(169)/500; x0=1+2/qn
    c0,c1=arb(1),x0
    for j in range(1,m):c0,c1=c1,2*x0*c1-c0
    D=c1; dr=arb(str(abs(den.p)))/arb(str(den.q))
    check('independent_arb_recurrence_contains_exact_denominator',D.overlaps(dr))
    check('denominator_positive',D>1)
    delta=qn**3/D; h2=(1+qn)**2/D; h3=(1+qn)**3/D
    pref=arb('1.1e-10')*h3
    oldpref=(1+qn)**5*arb('1.1e-10')*qn**29/(1-qn)
    geometric_exact_pref=arb('1.1e-10')*(1+qn)**3*qn**29
    ceilings={
      'original_defect_Cdef3_times_IY':arb('.5')/h3,
      'graph_Cg3_if_full_remaining_margin':arb('5.34235e-12')/pref,
      'physical_CX3_if_5e_minus9_margin':arb('5e-9')/pref,
      'feedback_Cg3_times_Ifb_if_other_debits_zero':arb('.0012069')/h2,
      'mass_Cmass3_times_Im_if_trial_and_other_errors_zero':arb('5e-5')/h2}
    check('complete_form_residual_below_half',delta<arb('.5'))
    check('strict_improvement_over_geometric_same_output_maps',pref<geometric_exact_pref and geometric_exact_pref<oldpref)
    check('positive_conditional_ceilings',all(v>0 for v in ceilings.values()))
    check('conditional_graph_budget_identity',(pref*ceilings['graph_Cg3_if_full_remaining_margin']).contains(arb('5.34235e-12')))
    check('conditional_defect_budget_identity',(h3*ceilings['original_defect_Cdef3_times_IY']).contains(arb('.5')))
    fast=arb('.028913')+arb('1e-11')*arb('1.0706e10')
    check('two_fast_screen_unchanged',fast<arb('.1599164'))
    return {'precision_bits':ctx.prec,'degree':N,'chebyshev_degree':m,
      'q_exact':str(q),'polynomial_pi31_exact_sha256':sha256(str(pip.as_expr()).encode()).hexdigest(),
      'denominator_exact':str(abs(den)),'denominator_D':enc(D),
      'complete_form_defect_upper':enc(delta),'output_factor_h2':enc(h2),'source_defect_factor_h3':enc(h3),
      'affine_error_prefactor':enc(pref),'previous_geometric_prefactor':enc(oldpref),
      'exact_signed_geometric_prefactor':enc(geometric_exact_pref),
      'improvement_over_recorded_geometric':enc(oldpref/pref),
      'improvement_over_signed_geometric':enc(geometric_exact_pref/pref),
      'conditional_ceilings':{k:enc(v) for k,v in ceilings.items()},
      'unchanged_two_fast_screen':enc(fast),
      'retained_joint_budgets':[
        'Y_other+4.084e-22 A_o+Y_phase <= 5.34235e-12',
        'eta_other+.0195239 gamma_E+3.712e-12 A_o J_E+eta_phase < .0012069'],
      'odd_power_erratum':'General power formula requires (-1)^p; previous absolute bounds unchanged.',
      'scope':{'actual_paired_domain_retained':True,'complete_common_restored_operator_used':True,
        'complete_form_residual_certified':True,'three_residual_factors_retained':True,
        'numerical_quantum_inverse_applied':False,'original_strong_delta_certified':False,
        'output_constants_numerically_bounded':False,'actual_kappa_m_computed':False,
        'signed_coupled_columns_computed':False,'nonlinear_admission_contraction':False,
        'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
        'all_sector_stability':False,'actual_controlled_solution':False,
        'construction_obstruction':False,'full_goal_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
      'input_sha256_before':hashes,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      'runtime':{'python':sys.version,'python_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute()
        run['input_sha256_after']={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as out:out.write(summary.encode())
    print(summary,flush=True)
    if run['status']!='passed':raise SystemExit(1)
