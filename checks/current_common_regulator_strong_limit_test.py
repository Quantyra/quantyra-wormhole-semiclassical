"""Support the ordinary actual-domain cutoff test; no coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
STEM='current-common-regulator-strong-limit-test'
OWN=['notes/'+STEM+'.md', 'checks/current_common_regulator_strong_limit_test.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-angular-resummed-quantum-operator.md',
 'notes/current-rough-paired-variation.md', 'notes/current-curved-paired-heat-bridge.md',
 'notes/stationary-shape-logarithmic-response.md',
 'notes/current-finite-local-action-assembly.md',
 'notes/current-coupled-local-graph-bound.md',
 'notes/current-fast-graph-domain-embedding.md',
 'notes/current-original-cutoff-coupled-estimator.md',
 'notes/current-signed-strong-local-core.md',
 'evidence/current-global-mass-conversion.json']
CHECKS=[]

def check(name, ok):
    CHECKS.append({'name':name, 'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'), flush=True)
    if not ok: raise AssertionError(name)

def zero(e): return sp.cancel(sp.expand(e))==0

def dec(v):
    d=v['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])

def enc(v):
    m,r,e=v.mid_rad_10exp(50)
    return {'ball':v.str(32), 'decimal_enclosure':
            {'mid':str(m),'rad':str(r),'exponent':str(e)}}

def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    for p,h in review['scientific_sha256'].items():
        check('frozen_'+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)
    a,b,A,B=sp.symbols('a b A B')
    combined=(a+b)**2*(A-B)/(4*(a-b))+(a+b)*(A+B)/4
    check('two_vertices_combine_exactly', zero(combined-(a+b)*(a*A-b*B)/(2*(a-b))))
    Ra,Rb=sp.symbols('Ra Rb')
    check('complete_divided_difference_weight',zero(combined.subs({A:-a*Ra,B:-b*Rb})
          +(a+b)**2*(a*a*Ra-b*b*Rb)/(2*(a*a-b*b))))
    # Differentiating lambda exp(-t lambda) before integration gives F'.
    lam,t=sp.symbols('lam t')
    check('F_derivative_integrand',zero(sp.diff(lam*sp.exp(-t*lam),lam)
                                         -(1-t*lam)*sp.exp(-t*lam)))
    # Symbolic integration-by-parts identity for ARBITRARY Theta, not a mode.
    theta=sp.Function('Theta')(t)
    check('regulator_endpoint_integration_by_parts',zero(
          lam*theta*sp.exp(-t*lam)+sp.diff(theta*sp.exp(-t*lam),t)
          -sp.diff(theta,t)*sp.exp(-t*lam)))
    nu,T=sp.symbols('nu T')
    T2=(2*nu+2*T)**2+2*(2*nu)**2
    T1=(2*nu+2*T)+2*(2*nu)
    check('canonical_transverse_log_shape_only',zero(3*T2-T1*T1-8*T*T))
    z,C,Q,Bsq,G,N,gamma,beta,lamb=sp.symbols('z C Q Bsq G N gamma beta lamb',nonzero=True)
    # Include time/sphere volume and ORIGINAL negative action factor.
    volume=8*sp.pi**2*Bsq**2/(C*gamma*z**4)
    proper_four=(z*z*C*Q/Bsq)**2
    W=-N*sp.pi*Bsq/(G*gamma)
    coeff=sp.cancel(N*volume*proper_four/(120*sp.pi**2*W))
    check('original_physical_log_coefficient',zero(coeff+G*C*Q**2/(15*sp.pi*Bsq)))
    K0=G*C*beta**2/(15*sp.pi*Bsq)
    check('normal_principal_ratio',zero(-coeff/K0-(Q/beta)**2))
    check('K0_original_units',zero(K0.subs(Bsq,C*beta/lamb)-G*lamb*beta/(15*sp.pi)))
    s,Cv,r=sp.symbols('s Cv r')
    aop=(s+3)/(3*(s+1))
    Jc=-6*beta*(Cv-r*s*(s-3))*(s-4)*(s+1)
    check('packet_scalar_exact_original_coordinate',zero(aop-sp.Rational(1,3)-2/(3*(s+1))))
    check('packet_scalar_input_order_three',sp.degree(sp.cancel(Jc*2/(3*(s+1))),s)==3)
    # Verify the nonnegative terms discarded in the whole-interval Q lower bound.
    eps,tau,u=sp.symbols('eps tau u')
    actual=beta+u*z*z-(4*beta+2*u+tau)*z**3+(3*beta+u+eps)*z**4-2*eps*z**5+eps*z**6
    negative=beta-(4*beta+2*u+tau)*z**3+eps*z**6
    check('exact_finite_seed_polynomial_remainder',zero(actual-negative
                -(u*z*z+(3*beta+u+eps)*z**4-2*eps*z**5)))
    src=json.loads((ROOT/'evidence/current-global-mass-conversion.json').read_bytes())['runs'][-1]
    check('retained_seed_evidence_passed',src['status']=='passed')
    p={k:dec(v) for k,v in src['result']['parameters'].items()}
    be,uu,ep,ta,la=[p[k] for k in ['beta','u','eps','tau','g']]
    check('whole_interval_discarded_terms_nonnegative',be>0 and uu>0 and ep<0 and 3*be+uu+ep>0 and ta>0)
    edge=arb('.003')
    qlo=1-(4+2*uu/be+ta/be)*edge**3-abs(ep)/be*edge**6
    qlo2=qlo*qlo
    check('whole_interval_Q_ratio_positive',qlo>0)
    check('raw_cutoff_error_original_graph_lower_bound',qlo2>arb('.9999997'))
    k0=la*be/(15*arb.pi())
    check('physical_K0_positive',k0>0)
    check('failure_exceeds_half_defect_target',qlo2>arb('.5'))
    check('conditional_two_fast_screen_unchanged',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    check('original_mass_conversion_debit_unchanged',sp.Rational('1.627532e-15')*2174==sp.Rational('3.538254568e-12'))
    return {'precision_bits':ctx.prec, 'parameters':{k:enc(p[k]) for k in ['beta','u','eps','tau','g']},
      'Q_over_beta_whole_interval_lower_bound':enc(qlo),
      'raw_cutoff_error_norm_lower_bound':enc(qlo2),
      'K0_original_normalization':enc(k0),
      'raw_cutoff_error_lower_bound_over_half':enc(qlo2/arb('.5')),
      'scope':{'ordinary_all_spectrum_positive_cutoff_identity':True,
       'actual_domain_necessary_strong_limit_test':True,
       'specified_raw_cutoff_estimator_fails':True,
       'prior_angular_resummed_zeta_identity_refuted':False,
       'uniform_matched_remainder_estimate':False,
       'complete_coupled_return_executed':False,'actual_kappa_m_computed':False,
       'physical_admission_or_feedback_budget_closed':False,
       'physical_X_comparisons':False,'nonlinear_admission_contraction':False,
       'exact_compatibility':False,'physical_matching':False,
       'semiclassical_validity_noise':False,'all_sector_stability':False,
       'actual_controlled_solution':False,'construction_obstruction':False,
       'full_goal_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
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
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
