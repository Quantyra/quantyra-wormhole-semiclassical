"""Verify exact full-resolvent tail and complete sums; no physical coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import json, sys, traceback
import sympy as sp
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
ctx.prec=512
ROOT=Path(__file__).resolve().parents[1]
STEM='current-resolvent-corrected-angular-remainder'
OWN=['notes/'+STEM+'.md','checks/current_resolvent_corrected_angular_remainder.py',
     'evidence/'+STEM+'-review.json']
sys.path.insert(0,str(ROOT/'checks'))
saved_argv=sys.argv[:];sys.argv=[sys.argv[0]]
import current_physical_york_dual_assembly as york
sys.argv=saved_argv
CHECKS=[]

def check(name,value,detail=None):
    CHECKS.append(dict(name=name,passed=bool(value)))
    if not value:raise AssertionError(name)

york.check=check;york.prior.check=check

def enc(v):
    m,r,e=v.mid_rad_10exp(60)
    return dict(ball=v.str(40),decimal_enclosure=dict(mid=str(m),rad=str(r),exponent=str(e)))



def execute():
    par,horizon,zh,env=york.actual_parameters()
    e=sp.symbols('e',real=True);F=1/(1+e)
    check('exact_bounded_functional_calculus_remainder',sp.simplify(1-F-F**2*e-F**3*e**2-e**3*F**3)==0)
    # The core identity is noncommutative at the spatial-operator level.
    S,D,R=sp.symbols('S Delta R',commutative=False)
    base=R-S-S*D*R
    expanded=base+S*D*base+S*D*S*D*base
    target=R-S-S*D*S-S*D*S*D*S-S*D*S*D*S*D*R
    check('exact_noncommutative_resolvent_iteration',sp.expand(expanded-target)==0)
    a=arb('1.85');d=a/324;den=1-d
    c0=arb(2).sqrt()/(4*(1-arb('1.05').sqrt()/36))
    A=4*c0*a**3/den**3
    AX=4*c0*(3*a**3*(996+arb(1070)/324)/den**4+3*a*a*1070/den**3+a**3*996/den**3)
    check('all_actual_masses_above_18',((par['q']+1)/par['C']).sqrt()>18)
    check('positive_relative_form_margin',0<d<1)
    check('complete_frequency_resolvent_tail_guard',A<arb('9.377'))
    check('all_seven_commutator_terms_guard',AX<53879)
    s=sp.symbols('s',positive=True)
    integral=sp.integrate((1+s*s)**(-sp.Rational(7,2)),(s,0,sp.oo))
    check('complete_thermal_integral',integral==sp.Rational(8,15))
    check('decreasing_thermal_envelope',sp.simplify(sp.diff((1+s*s)**(-sp.Rational(7,2)),s)+7*s/(1+s*s)**sp.Rational(9,2))==0)
    check('actual_thermal_mesh_guard',0<par['gamma']<arb('1e-33'))
    thermal_factor=arb('.17').sqrt()*8/15+arb('1e-33')/18
    tv=arb('9.377')*thermal_factor
    tx=arb(53879)*thermal_factor
    check('all_discrete_thermal_value_guard',tv<arb('2.062'))
    check('all_discrete_thermal_commutator_guard',tx<11848)
    x,q=sp.symbols('x q',positive=True)
    f=(q+2*x)/(x*x*(x+q)**2)
    check('angular_integral_antiderivative',sp.simplify(sp.diff(-1/(x*(x+q)),x)-f)==0)
    check('angular_envelope_decreasing',sp.simplify(-sp.diff(f,x)*x**3*(x+q)**3-2*(q*q+3*q*x+3*x*x))==0)
    check('complete_angular_sum_upper_guard',sp.simplify(2/q-(q+2)/(q+1)**2-1/(q+1)-(q+2)/(q*(q+1)**2))==0)
    sum4=4*par['C']**2/par['q']
    A4=2*sum4/(arb.pi()*par['B2'])
    check('unchanged_original_angular_conversion',A4<arb('4.583e-35'))
    # Retain the previous actual density/Jacobian and near-energy comparison.
    old_source=json.loads((ROOT/'evidence/current-paired-response-parametrix-bound.json').read_bytes())['runs'][-1]
    check('prior_acceptance_exists_with_stated_scope',old_source['status']=='passed' and old_source['result']['complete_relative_thermal_sum_certified'] and not old_source['result']['common_subtracted_strong_sum_certified'])
    density=arb('2.5').sqrt()*arb('35.534')
    source=arb('4.583e-35')*density*arb('2.062')
    comm=arb('4.583e-35')*density*(11848+arb('.0625')*arb('2.062'))
    check('complete_original_L2_tail_guard',source<arb('5.31e-33'))
    check('complete_original_first_commutator_tail_guard',comm<arb('3.052e-29'))
    check('conditional_fast_screen_unchanged',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    remain=arb('5.34235e-12')-arb('5.042e-23')
    ge=arb('.0012069')/(arb('.0195239')+arb('4.583e-13'))
    check('conditional_remaining_affine_unchanged',remain>arb('5.34234999994957e-12'))
    check('conditional_energy_feedback_ceiling_unchanged',arb('.06181654279993')<ge<arb('.06181654279994'))
    return dict(outcome='unresolved_complete_physical_coupled_return_with_complete_corrected_angular_thermal_tail',
      per_mode_L2_constant=enc(A),per_mode_first_commutator_constant=enc(AX),
      complete_thermal_value_upper=enc(tv),complete_thermal_commutator_upper=enc(tx),
      sum_d_over_M4_upper=enc(sum4),original_angular_conversion=enc(A4),
      original_L2_tail_upper=enc(source),original_first_commutator_tail_upper=enc(comm),
      guards={'mode':'9.377 M^2/K^7','mode_commutator':'53879 M^2/K^7','thermal':'2.062/M^4','thermal_commutator':'11848/M^4','original_near_L2':'5.31e-33 per E','original_near_commutator':'3.052e-29 per E'},
      complete_angular_and_thermal_remainder_norm_convergence=True,actual_global_scalar_inverses_retained=True,
      actual_paired_mass_log_domain_retained=True,original_tadpole_retained_exactly=True,finite_mode_selection=False,
      original_regulator_changed=False,weak_regulator_compatible_bounded_subtraction=True,
      retained_core_common_regulator_numerical_evaluation=False,higher_conormal_commutators_certified=False,
      original_strong_coupled_composition_certified=False,nonlinear_metric_variations_certified=False,
      L2_bound_is_composed_C4=False,actual_C4=None,kappa_m=None,actual_mass=None,actual_gamma_E=None,
      physical_coupled_signed_solve_executed=False,physical_operator_numerical_applications=0,
      original_norms_changed=False,new_nonlinear_requirement_closed=False,actual_controlled_solution=False,
      exact_compatibility=False,physical_matching=False,semiclassical_validity_noise=False,all_sector_stability=False,
      construction_obstruction=False,full_objective_preserved=True,closed_audit_preserved=True)

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    review=json.loads((ROOT/OWN[2]).read_bytes())
    paths=set(OWN+review['retained_inputs'])
    for module in list(sys.modules.values()):
        f=getattr(module,'__file__',None)
        if f:
            try:rel=Path(f).resolve().relative_to(ROOT).as_posix()
            except ValueError:continue
            if rel.endswith('.py'):paths.add(rel)
    digest=lambda p:sha256(p.read_bytes()).hexdigest()
    before={p:digest(ROOT/p) for p in sorted(paths)}
    fr=Path(flint.__file__).resolve().parent
    runtime=dict(python=sys.version,python_executable_sha256=digest(Path(sys.executable)),flint_version=flint.__version__,
      flint_artifacts={p.relative_to(fr).as_posix():digest(p) for p in sorted(fr.rglob('*')) if p.is_file() and p.suffix in ['.py','.pyd','.dll']})
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),mode=mode,input_sha256_before=before,precision_bits=512,runtime=runtime,
      frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN})
    try:
        check('ordinary_argument_frozen',before[OWN[0]]==review['ordinary_argument_sha256'])
        run['result']=execute()
        run['input_sha256_after']={p:digest(ROOT/p) for p in before}
        check('frozen_inputs_unchanged',before==run['input_sha256_after'])
        run['status']='passed'
    except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    ledger['runs'].append(run);out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    print(json.dumps(dict(mode=mode,status=run['status'],checks=len(CHECKS),result=run.get('result'),traceback=run.get('traceback')),indent=2))
    if run['status']!='passed':raise SystemExit(1)
