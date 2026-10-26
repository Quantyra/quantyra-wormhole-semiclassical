"""Exact complete-output factorization; its remaining joint norm is not computed."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb, ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from operators.signed_local_core import SignedLocalCore
STEM='current-complete-source-square-factorization'
OWN=['notes/'+STEM+'.md','checks/current_complete_source_square_factorization.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['operators/signed_local_core.py','notes/current-seed-full-hessian-realization.md',
 'notes/current-schur-residual-phase-domain.md','notes/current-schur-polynomial-return-test.md',
 'notes/current-global-mixed-schur-preconditioner.md','evidence/current-global-mixed-schur-preconditioner.json',
 'notes/current-weighted-bolt-source-return.md','evidence/current-weighted-bolt-source-return.json',
 'evidence/current-localized-bolt-adjoint-sensitivity.json','notes/current-direct-triple-return-attempt.md']
CHECKS=[]
def sh(p):return sha256(p.read_bytes()).hexdigest()
def check(name,ok):
    CHECKS.append({'name':name,'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'),flush=True)
    if not ok:raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
    if isinstance(x,dict):
        e=x['decimal_enclosure'];return arb(e['mid']+' +/- '+e['rad'])*arb(10)**int(e['exponent'])
    return arb(x.strip('[]'))
def prior(stem):
    r=json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]
    check('accepted_'+stem,r['status']=='passed')
    return r
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_checker_frozen',all(sh(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
    old=prior('current-weighted-bolt-source-return')
    check('accepted_output_proof_unchanged',sh(ROOT/'notes/current-weighted-bolt-source-return.md')==old['input_sha256_before']['notes/current-weighted-bolt-source-return.md'])
    src=prior('current-localized-bolt-adjoint-sensitivity')
    mixed=prior('current-global-mixed-schur-preconditioner')
    check('actual_mixed_remainder_retained',dec(mixed['result']['numerics']['complete_mixed_preconditioned_remainder_upper'])<arb('1.556e-23'))
    # Exact affine density variation, including off-shell contacts.
    t,T,v,w,lap=sp.symbols('t T V w DeltaT',real=True)
    wt=sp.exp(-2*t*T)*(w-v-2*t*lap)+v
    second=sp.diff(sp.exp(2*t*T)*wt**2,t,2).subs(t,0)
    target=8*(-lap+(v-w)*T)**2+4*(2*v*w-w*w)*T*T
    check('complete_affine_square_density_identity',sp.expand(second-target)==0)
    # Check ORIGINAL local action units, without constructing a matrix of modes.
    c=SignedLocalCore();z=c.z;nu,tt=c.fields
    cv=c.curvature.subs({a:0 for a in nu})
    sq=8*sp.pi*c.cR/(c.C*c.B2)*sp.exp(2*tt[0])*(2+c.C*sp.exp(-2*tt[0])*cv+12*c.C*c.beta/z**2)**2
    jets=sp.symbols('v0:3');alpha=sp.symbols('alpha')
    variation={tt[i]:alpha*jets[i] for i in range(3)}
    sq2=sp.diff(sq.subs(variation),alpha,2).subs(alpha,0)
    cv0=cv.subs({tt[i]:0 for i in range(3)})
    dlap=(c.Q[0]*(jets[2]-jets[1])+z*c.Q[1]*jets[1])/z**2
    V0=2*z*z/c.C+12*c.beta;wb=z*z*cv0+V0
    JT=-z*z*dlap-z*z*cv0*jets[0]
    dc=96*sp.pi*c.cR*c.C/c.B2
    normalized=((2*dc/3)*JT**2+(dc/3)*(2*V0*wb-wb**2)*jets[0]**2)/z**4
    check('square_coefficient_in_original_local_action_units',sp.expand(sq2-normalized)==0)
    Z,R,X=sp.symbols('Z R X',commutative=False)
    remainder=sp.expand(Z*(1-X)-(1-Z*(R+X)))
    check('noncommuting_inverse_remainder_identity',sp.expand(remainder.subs(Z*R,1-Z))==0)
    # Parameters and output conversions remain the original ones.
    pars=src['result']['parameters']
    h,gamma,beta,C,B2=[dec(pars[k]) for k in ['h','gamma','beta','C','B2']]
    g=C*beta/B2;dc=96*arb.pi()*arb('1e22')*C/B2
    d=dc-C/(60*arb.pi()*B2);rm=arb('1.556e-23')
    # The exact positive difference is much smaller than the shared input radius.
    check('distinguish_square_and_trace_diffusion',C/(60*arb.pi()*B2)>0 and arb('1e-8')<d and dc<arb('2e-8'))
    Vmax=2/C+12*beta;wm=arb('1.1e-27')*beta
    tm=beta*(arb('1.1e-27')+12*11*g/(360*arb.pi()))
    mm=2*arb('2.001')*wm+wm*wm
    check('actual_off_shell_guards',Vmax<arb('2.001') and 0<wm<tm<arb('2e-28'))
    q=arb(169)/500;f=q/(1+q);ellb=(1+h)/h
    lp=ellb*(arb('2.5')*h)**(arb(7)/5)/gamma.sqrt()
    A0=arb('17.6')*f**3*h**(-arb(7)/5)
    R2=lp*arb('2e6')*dc/3*(rm+26*tm)/(1-rm)
    R0=lp*arb('2e6')*dc/6*26*mm/(1-rm)
    r=old['result'];pref=dec(r['affine_prefactor']);Cref=dec(r['reference_three_residual_map']['C_b3_ref_upper'])
    reference=dec(r['reference_three_residual_map']['original_affine_bolt_debit_upper'])
    check('same_original_point_conversion',(lp*3000).overlaps(dec(r['actual_output_factors']['original_physical_bolt_factor_from_weighted_L2'])))
    check('same_reference_debit',(pref*Cref).overlaps(reference))
    debit0=pref*R0*A0;alpha2=pref*R2;budget=arb('1.667e-9')
    remaining=budget-reference-debit0
    check('complete_off_shell_square_debit',debit0<arb('1e-18'))
    check('conditional_remaining_margin_positive',remaining>0)
    Jceiling=remaining/pref;A2ceiling=remaining/alpha2
    square_slice=3*Jceiling/dc
    check('whole_budget_joint_slice',(pref*Jceiling+reference+debit0).overlaps(budget))
    check('whole_budget_derivative_slice',(alpha2*A2ceiling+reference+debit0).overlaps(budget))
    # Complete correction: estimate the opposite-weight adjoints before norms.
    s=arb(7)/5;mD=1-d*s*s
    check('resolvent_difference_contraction_hypothesis',mD>arb('.5'))
    gD=1/(2*(d*mD).sqrt())+s/mD
    KJL=1/d+2*gD+(arb('2.001')+54*beta+tm)/mD
    t1=arb('1.4e-26')*beta**(arb(3)/2);t2=arb('8.6e-25')*beta**2
    Dt=1900*tm+132*t1+26*t2
    EJ=2*d*Dt+C/(60*arb.pi()*B2)*arb('1e12')+40*beta*arb('4e-27')*26/(60*arb.pi())
    KJP=1900+132+26*(arb('2.001')+54*beta+tm)
    Kt=tm*KJP+132*t1+26*(t2+2*t1)
    Ksum=KJL*EJ/2+Kt
    Rmix0=lp*arb('2e6')*dc/3*Ksum/(1-rm)
    debitmix=pref*Rmix0*A0
    revised_remaining=remaining-debitmix
    check('entire_mixed_square_return_debit',debitmix<arb('1e-9'))
    check('joint_main_margin_after_actual_remainder_bound',revised_remaining>0)
    revised_J=revised_remaining/pref
    check('complete_revised_budget_identity',(pref*revised_J+debitmix+reference+debit0).overlaps(budget))
    return {'precision_bits':512,'parameters':{'dc':enc(dc),'d':enc(d),'t0_guard':enc(tm),'w_guard':enc(wm),'m_guard':enc(mm)},
      'known_original_maps':{'A0_upper':enc(A0),'physical_point_conversion_without_numerator':enc(lp),
       'R2_coefficient':enc(R2),'R0_coefficient':enc(R0),'off_shell_square_affine_debit_upper':enc(debit0),
       'reference_affine_debit_upper':enc(reference),
       'mD':enc(mD),'J_adjoint_LD_inverse_upper':enc(KJL),'entire_E_P_inverse_upper':enc(EJ),
       'adjoint_combined_remainder_Ksum_upper':enc(Ksum),'mixed_correction_A0_coefficient':enc(Rmix0),
       'complete_mixed_square_affine_debit_upper':enc(debitmix)},
      'complete_budget':{'inequality':'p_aff J_b + debit_mix + debit0 + debit_ref <= B_remaining',
       'p_aff':enc(pref),'alpha2':enc(alpha2),'full_margin_only':enc(budget),'remaining_after_known_debits':enc(remaining),
       'J_b_actual':None,'A2_actual':None,'M_rem_actual':None,
       'conditional_J_b_ceiling_A2_zero':enc(Jceiling),'conditional_A2_ceiling_J_b_zero':enc(A2ceiling),
       'initial_filtered_square_sensitivity_before_mixed_debit':enc(square_slice),
       'revised_remaining_after_all_proved_square_debits':enc(revised_remaining),
       'revised_conditional_J_b_ceiling':enc(revised_J),
       'revised_filtered_square_sensitivity_other_main_zero':enc(3*revised_J/dc),
       'mixed_remainder_fraction_of_full_margin':enc(debitmix/budget),
       'A2_removed_from_executed_remainder_bound':True,
       'sensitivity_slices_not_simultaneous_allocations':True},
      'scope':{'exact_actual_domain_square_inverse_factorization':True,'off_shell_density_debit_bounded':True,
       'entire_mixed_square_correction_bounded_without_A2':True,
       'full_remaining_source_norm_computed':False,'joint_main_norm_computed':False,'A2_computed':False,
       'full_Cb3_certified':False,'source_estimator_failure_proved':False,'construction_obstruction':False,
       'actual_quantum_inverse_applied':False,'original_delta_certified':False,'all_coupled_entries_computed':False,
       'actual_kappa_m_computed':False,'nonlinear_admission_contraction':False,'exact_compatibility':False,
       'physical_matching':False,'semiclassical_validity_noise':False,'all_sector_stability':False,'full_goal_complete':False}}
if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sh(ROOT/p) for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      'runtime':{'python':sys.version,'python_sha256':sh(Path(sys.executable)),'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={p:sh(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8') as log:log.write(summary)
    print(summary,flush=True)
    sys.exit(0 if run['status']=='passed' else 1)
