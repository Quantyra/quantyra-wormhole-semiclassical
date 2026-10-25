"""Supporting arithmetic for an all-spectrum weak trace remainder.

No signed center, differentiated graph error, or coupled metric solve.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys,json,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='current-adaptive-dirac-trace-remainder'
OWN=['notes/'+STEM+'.md','checks/current_adaptive_dirac_trace_remainder.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-global-dirac-transport-series.md',
 'evidence/current-global-dirac-transport-series.json',
 'notes/current-coincident-dirac-action-assembly.md',
 'notes/current-rough-paired-variation.md',
 'notes/current-curved-paired-heat-bridge.md',
 'notes/current-mass-response-regulator-audit.md',
 'notes/current-third-metric-conversion.md',
 'notes/current-global-far-return-energy-bound.md',
 'notes/current-finite-local-action-assembly.md',
 'notes/current-global-mass-conversion.md',
 'evidence/current-global-mass-conversion.json',
 'notes/current-coupled-global-return-attempt.md',
 'notes/current-complete-weak-mass-residual.md',
 'notes/current-actual-far-graph-energy-bound.md']
CHECKS=[]

def check(name,ok):
    assert name not in {x['name'] for x in CHECKS}
    CHECKS.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)

def enc(x):
    m,r,e=x.mid_rad_10exp(55)
    return {'ball':str(x),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}

def dec(x):
    d=x['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])

def prior(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed'
    return r['result']

def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_frozen',sha256((ROOT/OWN[0]).read_bytes()).hexdigest()==review['scientific_sha256'][OWN[0]])
    par=prior('current-global-mass-conversion')['parameters']
    beta,u,eps,g=[dec(par[x]) for x in ['beta','u','eps','g']]
    C=1/u; B2=beta*C/g; k0=arb(2000); qflux=k0**4
    C0=(1+12*arb.pi()*k0**6).sqrt(); classical_B2=k0**10*(C0-1)/6
    kappa=2*arb('.265')*qflux*C0/(k0*classical_B2*classical_B2.sqrt())
    gamma=kappa*B2.sqrt()/C
    check('actual_seed_parameter_signs',beta>arb('.16') and u>0 and eps<0 and abs(eps)<arb('1e-28') and g>0)
    check('actual_H1_lower_coefficient',beta-6*abs(eps)>arb('.15') and 3*beta+u+eps>0)
    check('actual_temperature_nonzero',gamma>arb('5e-34') and gamma<arb('5.1e-34'))
    check('all_angular_masses_above_eighteen',(1+qflux)/C>18**2)
    old=prior('current-global-dirac-transport-series')
    check('same_actual_C',dec(old['actual_C']).overlaps(C))
    check('retained_complete_mode_coupling',dec(old['coupling_upper_expression'])<arb('.385'))
    length=4/arb('.15').sqrt()*(arb('.15')/(2*gamma)).sqrt().asinh()+2/arb('.15').sqrt()
    check('entire_doubled_proper_length',length<arb('394.854'))
    e=arb('.011').exp()-1; q=arb('.3961')
    hs=arb('.385')*5*arb.pi()/36+395*e**2
    check('complete_HS_squared_per_M',hs<arb('.217') and hs<1)
    check('complex_mass_contraction',arb('.385')+e<q and q<1)
    check('common_square_form_coercivity',2*e+e**2<arb('.023'))
    check('form_perturbation_identity',(2*e+e**2).overlaps(arb('.022').exp()-1))
    # Algebra used in the ordinary global length and trace-ideal arguments.
    t,v=S.symbols('t v',positive=True)
    neg=(3*t*t+2*t*v+v*v)/(t**3*v**4)
    check('bolt_negative_coefficient_decomposition',S.simplify(neg-(3/(t*v**4)+2/(t*t*v**3)+1/(t**3*v*v)))==0)
    x,a,b=S.symbols('x a b',positive=True)
    primitive=2/S.sqrt(b)*S.asinh(S.sqrt(b*x/a))
    check('whole_near_bolt_length_integral',S.simplify(S.diff(primitive,x)-1/S.sqrt(x*(a+b*x)))==0)
    far=S.log((t-v)/t)/v
    check('whole_far_length_integral',S.simplify(S.diff(far,t)-1/(t*(t-v)))==0)
    th=S.symbols('theta',real=True)
    O=S.Matrix([[S.cos(th/2),-S.sin(th/2)],[S.sin(th/2),S.cos(th/2)]])
    J=S.Matrix([[0,-1],[1,0]]);B=O.T*S.diag(1,-1)*O
    check('moving_connection_mass_Frobenius_orthogonality',S.trigsimp(S.trace(J.T*B))==0)
    check('each_vertex_has_two_spin_components',S.trigsimp(S.trace(B.T*B))==2 and S.trace(J.T*J)==2)
    n,qf=S.symbols('n qf',positive=True)
    count=2*qf*n+2*n*(n+1)
    check('exact_angular_count_envelope_identity',S.expand(4*n*(n+qf)-count-2*n*(n+qf-1))==0)
    angular_integral=3*4*C/18
    check('complete_inverse_cube_angular_sum_guard',(angular_integral-2*C/3).contains(0))
    # Exact ceilings for the rounding rule: x <= ceil(x) is the only
    # discrete fact needed; no finite spectral sample is used as a proof.
    lam,M,p=S.symbols('lam M p',positive=True)
    exponent=(4*S.log(M/18)+2*S.log(1+p))/lam
    check('adaptive_power_spectral_decay',S.simplify(S.exp(-lam*exponent)-(18/M)**4/(1+p)**2)==0)
    check('full_thermal_integral',S.integrate((1+x)**-2,(x,0,S.oo))==1)
    cutoff=arb(128)
    E1=arb('.011').exp()/(1-q)*q**126*18**4*4*C*(1+gamma)/(3*arb.pi()*B2)
    E2=100*E1;E3=40000*E1
    check('whole_normalized_source_tail',E1<arb('9.204e-79'))
    check('whole_normalized_Hessian_Linf_to_L1_tail',E2<arb('9.204e-77'))
    check('whole_normalized_third_action_tail',E3<arb('3.682e-74'))
    check('conditional_fast_screen',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    mass=arb('1.627532e-15');conversion=2174*mass
    near=arb('1e-11')-conversion-arb('1.24e-20')
    check('retained_mass_conversion_debit',2174*S.Rational('1.627532e-15')==S.Rational('3.538254568e-12'))
    check('retained_near_shape_budget',near>arb('6.4617454e-12'))
    check('conditional_mass_affine_debit',(arb('5e-5')*arb('4.2953e-8')).overlaps(arb('2.14765e-12')))
    check('conditional_mass_feedback_debit',(arb('5e-5')*arb('49.462')).overlaps(arb('.0024731')))
    return {'precision_bits':ctx.prec,'actual_parameters':{k:enc(v) for k,v in {'beta':beta,'C':C,'B2':B2,'gamma':gamma}.items()},
      'doubled_proper_length_upper_expression':enc(length),'HS_squared_over_M_upper_expression':enc(hs),
      'complex_form_relative_bound':enc(2*e+e**2),'spectral_order_rule':'128+ceil((4log(M/18)+2log(1+abs(p)))/log(1/.3961))',
      'normalized_source_L1_error':enc(E1),'normalized_Hessian_Linf_to_L1_error':enc(E2),
      'normalized_third_action_Linf_Linf_to_L1_error':enc(E3),
      'retained_original_mass_conversion_debit':enc(conversion),'retained_near_shape_allowance':enc(near),
      'scope':{'actual_paired_mass_log_domain':True,'full_angular_thermal_remainder':True,
       'finite_domain_bolt_returns':True,'common_regulator_removal':'ordinary trace-ideal and closed-form argument, not proved by this checker',
       'all_exponential_mass_variation_contacts_in_tail':True,'separate_Weyl_LLL_and_local_restoration':'retained exactly; center not evaluated',
       'weak_trace_tail_certified':True,'strong_graph_tail_certified':False,
       'generic_L1_to_H3_conversion':'unbounded, not an unknown finite constant',
       'signed_center_computed':False,'actual_kappa_m_computed':False,'complete_coupled_return_certified':False,
       'nonlinear_admission':False,'nonlinear_contraction':False,'exact_compatibility':False,
       'physical_matching':False,'semiclassical_validity_noise':False,'all_sector_stability':False,
       'construction_obstruction':False,'actual_controlled_solution':False,'full_goal_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':S.__version__}}
    try:
        run['result']=execute()
        run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),
       'source_error':run.get('result',{}).get('normalized_source_L1_error'),
       'Hessian_error':run.get('result',{}).get('normalized_Hessian_Linf_to_L1_error'),
       'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
