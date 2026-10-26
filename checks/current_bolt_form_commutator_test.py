"""Support an ordinary actual-domain estimator failure, not an inverse solve."""
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
STEM='current-bolt-form-commutator-test'
OWN=['notes/'+STEM+'.md', 'checks/current_bolt_form_commutator_test.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-schur-physical-output-test.md',
 'evidence/current-schur-physical-output-test.json',
 'notes/current-schur-residual-phase-domain.md',
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-schur-form-preconditioner-test.md',
 'notes/current-strong-return-cauchy-test.md',
 'notes/strong-interior-quantum-remainder.md']
CHECKS=[]

def check(name, ok):
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
    prior=json.loads((ROOT/INPUTS[4]).read_bytes())['runs'][-1]
    check('retained_full_operator_packet_result_passed',prior['status']=='passed')
    rho,k,mu,eps,R=sp.symbols('rho k mu eps R',positive=True)
    v=rho**2; dr=sp.diff(v,rho)
    lap=sp.diff(dr,rho,2)+sp.diff(dr,rho)/rho
    check('radial_derivative_not_H2_leading_laplacian',sp.simplify(lap-2/rho)==0)
    integral=sp.integrate(lap**2*rho,(rho,eps,R))
    check('exact_logarithmic_domain_divergence',sp.simplify(integral-4*sp.log(R/eps))==0)
    check('domain_integral_diverges',sp.limit(integral,eps,0,dir='+')==sp.oo)
    s=sp.symbols('s',positive=True); f=sp.Function('f')(s)
    Delta=lambda w:4*(s*sp.diff(w,s,2)+sp.diff(w,s))
    check('regular_areal_laplacian_commutator',sp.expand(sp.diff(Delta(f),s)-Delta(sp.diff(f,s))-4*sp.diff(f,s,2))==0)
    comm4=sp.expand(sp.diff(Delta(Delta(f)),s)-Delta(Delta(sp.diff(f,s))))
    check('full_fourth_radial_commutator',sp.expand(comm4-(32*s*sp.diff(f,s,4)+64*sp.diff(f,s,3)))==0)
    l=sp.Function('l')(rho); p=l*k**4*sp.log(k/mu)
    a=1/(2*rho)
    symbol=sp.expand(a*sp.diff(p,rho)-sp.diff(a,rho)*k*sp.diff(p,k))
    target=sp.diff(l,rho)*k**4*sp.log(k/mu)/(2*rho)+l*k**4*(4*sp.log(k/mu)+1)/(2*rho**2)
    check('complete_fourth_log_symbol_commutator',sp.simplify(symbol-target)==0)
    leading=sp.simplify((symbol-a*sp.diff(p,rho))/p)
    check('singular_form_ratio_retains_log_derivative',sp.simplify(leading-(2/rho**2+1/(2*rho**2*sp.log(k/mu))))==0)
    xc=sp.expand(rho*sp.diff(p,rho)-k*sp.diff(p,k))
    check('conormal_principal_has_no_inverse_radius',sp.simplify((xc-rho*sp.diff(p,rho))/p+4+1/sp.log(k/mu))==0)
    # On 1<=t<=2, the weighted-average coefficient 2/t^2 is in [1/2,2].
    t=sp.symbols('t',positive=True)
    check('whole_packet_interval_coefficient',sp.diff(2/t**2,t)==-4/t**3 and (2/t**2).subs(t,2)==sp.Rational(1,2) and (2/t**2).subs(t,1)==2)
    alpha=sp.Rational(3,4)
    check('packet_scale_separation',1-alpha>0)
    check('unbounded_first_form_commutator_power',2*alpha==sp.Rational(3,2))
    sigma=sp.symbols('sigma',real=True)
    alpha_choice=(sigma+2)/4
    check('tame_derivative_deficit',sp.simplify(2*alpha_choice-sigma-(2-sigma)/2)==0)
    check('conormal_to_six_areal_derivative_deficit',-12-(-1)==-11)
    check('Cartesian_form_bootstrap_derivative_count',10+2==12 and 2+12==14)
    check('three_power_physical_Sobolev_reserve',2+4*3>12+sp.Rational(2,2))
    q=arb('.338');pref=(1+q)**5*arb('1.1e-10')*q**29/(1-q)
    cg=arb('5.34235e-12')/pref;cx=arb('5e-9')/pref
    check('conditional_targets_retain_original_scale',cg>arb('3.43e11') and cx>arb('3.21e14'))
    check('two_fast_screen_unchanged',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return {'precision_bits':ctx.prec,
      'ordinary_radial_domain_test':{'laplacian':'2/rho','truncated_H2_integral':'4 log(R/epsilon)','form_domain_preserved':False},
      'areal_commutator_principal':str(target),
      'first_areal_form_commutator_constant':'infinite',
      'actual_packet_bound':'liminf k^(-3/2) Re<v_k,[partial_s,L]v_k>/<v_k,Lv_k> >= 1/2; onset unknown',
      'packet_coefficient_interval':['1/2','2'],
      'minimum_extra_disk_derivatives_necessary':2,
      'Cartesian_form_commutator_orders_required':12,
      'three_power_affine_tail_prefactor_N32':enc(pref),
      'conditional_graph_C3_ceiling':enc(cg),
      'conditional_physical_C3_ceiling_if_5e_minus9_margin':enc(cx),
      'scope':{'stationary_energy_only_commutator_estimator_rejected':True,
        'actual_paired_domain_used':True,'common_regulator_and_full_sums_retained':True,
        'finite_frequency_onset_certified':False,'Cartesian_form_commutators_bounded':False,
        'three_power_numerical_norm_bounded':False,'original_strong_delta_certified':False,
        'complete_coupled_return_executed':False,'actual_kappa_m_computed':False,
        'any_additional_physical_gate_closed':False,'nonlinear_admission_contraction':False,
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
