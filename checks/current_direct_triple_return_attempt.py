"""Ordered identities and budget arithmetic; no quantum inverse application."""
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
STEM='current-direct-triple-return-attempt'
OWN=['notes/'+STEM+'.md','checks/current_direct_triple_return_attempt.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-interior-translation-jet-bound.md',
 'evidence/current-interior-translation-jet-bound.json',
 'notes/current-schur-polynomial-return-test.md',
 'evidence/current-schur-polynomial-return-test.json',
 'notes/current-schur-residual-phase-domain.md',
 'notes/current-physical-york-dual-assembly.md',
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-geodesic-two-sided-remainder.md']
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
    check('frozen_proof_and_checker',all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in review['scientific_sha256'].items()))
    check('accepted_inputs_passed',all(json.loads((ROOT/p).read_bytes())['runs'][-1]['status']=='passed' for p in (INPUTS[4],INPUTS[6])))
    # Formal rewriting uses only ell F=B(1-F), never commutation of B and ell.
    ell,B,F,V,A,C,D,K=sp.symbols('ell B F V A C D K',commutative=False)
    lhs=(ell*F**3).subs(ell*F**3,B*(1-F)*F**2)
    check('ordered_source_cancellation',sp.expand(lhs-B*F**2*(1-F))==0)
    # L=D+C* A^-1 C. Treat each ordered inverse block as its own symbol.
    Ct,Ai,Cti,Vti,T3,u2=sp.symbols('Ct Ai Cti Vti T3 u2',commutative=False)
    # D T3 = Vti B u2 - Ct Ai C T3; left Cti Ct=I.
    recon=Cti*(Vti*B*u2-Ct*Ai*C*T3)
    recon=sp.expand(recon).subs(Cti*Ct,sp.Integer(1))
    target=-Ai*C*T3+Cti*Vti*B*u2
    check('ordered_scalar_reconstruction_plus_sign',sp.expand(recon-target)==0)
    # First row of this error is -C S T3; second row is exactly zero.
    S=sp.symbols('S',commutative=False)
    row1=sp.expand(-A*target-C*T3)
    row1=row1.subs(A*Ai,1)
    check('first_row_retains_source',sp.expand(row1+A*Cti*Vti*B*u2)==0)
    # The identity C S T3 = A Cti Vti B u2 follows from S=K L and C K=A Cti.
    check('first_row_schur_source_identity',sp.expand(C*K*Vti*B*u2).subs(C*K,A*Cti)==A*Cti*Vti*B*u2)
    t=sp.symbols('t',nonnegative=True)
    check('spectral_F_identity',sp.cancel((1+t)*t/(1+t)-t)==0)
    q=arb(169)/500; f=q/(1+q)
    x0=1+2/q; c0,c1=arb(1),x0
    for j in range(1,29): c0,c1=c1,2*x0*c1-c0
    den=c1; h2=(1+q)**2/den; h3=(1+q)**3/den
    pref=arb('1.1e-10')*h3
    log10=lambda x:x.log()/arb(10).log()
    logK=arb(4003)+log10(arb(4))
    logL=(arb(4)/3)*(logK+log10(arb(2)))
    logM=logK+arb(53)/4*logL
    logA=log10(arb(2))+logM+19+3*log10(f)
    log_reported_bolt=logA-177
    log_affine_bolt=log_reported_bolt+log10(pref)
    log_excess=log_affine_bolt-log10(arb('1.667e-9'))
    check('local_measure_already_in_value_allowance',arb('3.207e18')*arb('.027').exp()<arb('1e19'))
    check('bolt_bump_H14_majorant',3**14<10**7 and arb(3)**7*arb('1e201')<arb('1e210'))
    check('generic_positive_bolt_expression_excess',log_excess>74000)
    ceilings={
      'original_defect_Cdef3_IY':arb('.5')/h3,
      'graph_Cg3_full_margin':arb('5.34235e-12')/pref,
      'physical_CX3_full_5e_minus9_margin':arb('5e-9')/pref,
      'bolt_Cb3_full_1p667e_minus9_margin':arb('1.667e-9')/pref,
      'direct_graph_Mg2_full_margin':arb('5.34235e-12')/pref/f,
      'direct_source_Md2_IY':arb('.5')/h3,
      'mass_Cmass3_Im_zero_trial_allowance':arb('5e-5')/h2,
      'feedback_Cg3_Ifb_zero_other_debits':arb('.0012069')/h2}
    check('direct_graph_factor_matches_original_budget',(ceilings['direct_graph_Mg2_full_margin']*f*pref).contains(arb('5.34235e-12')))
    check('direct_source_factor_matches_original_budget',(ceilings['direct_source_Md2_IY']*h3).contains(arb('.5')))
    check('bolt_original_units_budget',(ceilings['bolt_Cb3_full_1p667e_minus9_margin']*pref).contains(arb('1.667e-9')))
    fast=arb('.028913')+arb('1e-11')*arb('1.0706e10')
    check('limited_fast_screen',fast<arb('.1599164'))
    return {'precision_bits':512,'f':enc(f),'h2':enc(h2),'h3':enc(h3),'affine_prefactor':enc(pref),
      'generic_positive_expression':{'log10_K':enc(logK),'log10_L':enc(logL),'log10_M10':enc(logM),
       'log10_A_gen':enc(logA),'physical_bolt_conversion_lower_bound':'1e-177',
       'log10_reported_bolt_constant_lower':enc(log_reported_bolt),
       'log10_affine_bolt_allowance_lower':enc(log_affine_bolt),
       'log10_excess_over_original_bolt_budget':enc(log_excess),
       'is_lower_bound_on_actual_operator':False,'missing_principal_estimates_granted_only_for_estimator_test':True},
      'direct_identities':{'source':'J_F^-1 ell F^3 = K(V*)^-1 B F^2(I-F)',
       'scalar':'nu3=-A^-1 C V F^3 u+(C*)^-1(V*)^-1 B F^2 u',
       'graph':'C_g3 <= [q/(1+q)] M_g2',
       'source_last_H14_quantum_bootstrap_avoided':True,'graph_third_quantum_bootstrap_avoided':True},
      'conditional_ceilings':{k:enc(v) for k,v in ceilings.items()},
      'scope':{'actual_paired_operator_retained':True,'generic_positive_expression_fails':True,
       'direct_composite_norms_enclosed':False,'bolt_variational_witness_executed':False,
       'original_strong_delta_certified':False,'numerical_quantum_inverse_applied':False,
       'actual_kappa_m_computed':False,'signed_coupled_columns_computed':False,
       'nonlinear_admission_contraction':False,'exact_compatibility':False,'physical_matching':False,
       'semiclassical_validity_noise':False,'all_sector_stability':False,
       'actual_controlled_solution':False,'construction_obstruction':False,'full_goal_complete':False}}

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
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as out:out.write(summary.encode())
    print(summary,flush=True)
    if run['status']!='passed':raise SystemExit(1)
