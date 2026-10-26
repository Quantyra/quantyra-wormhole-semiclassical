"""Supporting arithmetic for the actual strong quotient reconstruction test.

No retained signed center, physical correction or mass-column certificate.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb,factorial
import sys,json,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1];STEM='current-strong-quotient-reconstruction-test'
OWN=['notes/'+STEM+'.md','checks/current_strong_quotient_reconstruction_test.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-finite-local-action-assembly.md',
 'notes/current-canonical-interior-schur-bound.md',
 'notes/current-coupled-truncation-energy-error.md',
 'evidence/current-coupled-truncation-energy-error.json',
 'notes/current-global-far-return-energy-bound.md',
 'evidence/current-global-mass-conversion.json',
 'evidence/current-adaptive-dirac-trace-remainder.json']
CHECKS=[]
def check(name,ok):
    assert name not in {x['name'] for x in CHECKS}
    CHECKS.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(55)
    return {'ball':str(x),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
    d=x['decimal_enclosure'];return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def prior(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed';return r['result']
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_frozen',sha256((ROOT/OWN[0]).read_bytes()).hexdigest()==review['scientific_sha256'][OWN[0]])
    p=prior('current-global-mass-conversion')['parameters']
    pp=prior('current-adaptive-dirac-trace-remainder')['actual_parameters']
    beta,u,eps,tau=[dec(p[k]) for k in ['beta','u','eps','tau']]
    C,B2=[dec(pp[k]) for k in ['C','B2']]
    DR=96*arb.pi()*arb('1e22')-1/(60*arb.pi())
    rho=6*DR*C/B2;sig=arb('4.005');lo=arb('.14');hi=arb('.16')
    qlo=beta-(4*beta+2*u+tau)*hi**3-2*abs(eps)*hi**5
    check('actual_positive_parameters',beta>0 and C>0 and B2>0 and u>0 and tau>0)
    check('discarded_polynomial_group_nonnegative',3*beta+u+eps-abs(eps)*hi**2>0)
    check('whole_annulus_Q_lower',qlo>arb('.16'))
    check('complete_scalar_fourth_coefficient',rho>arb('1e-7') and rho<arb('1.1e-7'))
    cp=(rho*beta).sqrt()*arb('.16')*lo**((13-2*sig)/2)
    check('uniform_asymptotic_ratio_coefficient',cp>arb('1.5e-7'))
    z=S.symbols('z',positive=True);sigma=S.symbols('sigma',real=True)
    rs,cr,cs,b2=S.symbols('rho c_R C B2',positive=True)
    q=S.Function('Q')(z);s=S.Function('s')(z)
    row=z**(4-sigma)*rs*q**2*S.diff(s,z,4)
    for _ in range(3):row=S.expand(z*S.diff(row,z))
    check('original_log_graph_highest_derivative',S.simplify(row.coeff(S.diff(s,z,7))-rs*q**2*z**(7-sigma))==0)
    ar=8*S.pi*cr/(cs*b2)
    check('local_action_complete_scalar_principal',S.simplify(72*ar*cs**2-cs/(10*S.pi*b2)-6*(96*S.pi*cr-1/(60*S.pi))*cs/b2)==0)
    er=prior('current-coupled-truncation-energy-error')
    e=dec(er['complete_energy_defect']);ef=dec(er['complete_source_dual_error'])
    h=arb('1.338').sqrt();Kt=7/(1-7*e)
    se=e*(1+7*h+(h+e)*Kt*(1+7*h))
    xe=(ef+e*arb('1.1e-10'))/(1-e)
    check('quotient_energy_error',se<arb('9e-24'))
    check('defined_linear_solutions_energy_error',xe<arb('1.123e-35'))
    check('approximate_quotient_inverse_upper',1/(1-e)<arb('1.0000001'))
    ca=arb('5.34235e-12')/arb('4.084e-22')
    cf=arb('.0012069')/arb('3.712e-12')
    cg=arb('.0012069')/arb('.0195239')
    check('conditional_affine_amplifier_ceiling',ca>arb('1.30e10') and ca<arb('1.31e10'))
    check('conditional_feedback_product_ceiling',cf>arb('3.25e8') and cf<arb('3.26e8'))
    check('conditional_energy_comparison_ceiling',cg>arb('.0618') and cg<arb('.06182'))
    fast=arb('.028913')+arb('1e-11')*arb('1.0706e10')
    check('conditional_two_fast_screen',fast<arb('.1599164'))
    return {'precision_bits':ctx.prec,'actual_parameters':{'beta':enc(beta),'C':enc(C),'B2':enc(B2)},
     'complete_scalar_fourth_coefficient_rho':enc(rho),'whole_annulus_Q_lower':enc(qlo),
     'asymptotic_graph_to_energy_ratio_coefficient_lower':enc(cp),
     'ratio_scope':'G_J/E = K^5 (c_psi+o(1)); this is not a finite-K interval test.',
     'exact_quotient_energy_inverse_upper':1,'approximate_quotient_energy_inverse_upper':enc(1/(1-e)),
     'quotient_energy_error_upper_expression':enc(se),'defined_linear_solutions_energy_error_upper_expression':enc(xe),
     'conditional_ceilings':{'affine_A_o':enc(ca),'feedback_A_o_times_J_E':enc(cf),'gamma_E':enc(cg)},
     'two_fast_conditional_screen':enc(fast),'two_fast_screen_reserve':enc(arb('.1599164')-fast),
     'scope':{'ordinary_actual_domain_transition_counterexample':True,
       'energy_only_transition_source_estimator_has_no_finite_constant':True,
       'actual_paired_mass_log_domain_preserved':True,'common_regulator_complete_spectrum_preserved':True,
       'signed_center_computed':False,'strong_quotient_parametrix_computed':False,'actual_kappa_m_computed':False,
       'complete_coupled_bound_certified':False,'nonlinear_admission':False,'nonlinear_contraction':False,
       'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
       'all_sector_stability':False,'actual_controlled_solution':False,'construction_obstruction':False,
       'full_goal_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':S.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'supporting_checks':len(CHECKS),
       'result':run.get('result'),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
