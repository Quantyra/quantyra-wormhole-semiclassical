"""Verify signed auxiliary asymptote and coupled estimator test; no physical solve."""
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
STEM='current-relative-angular-coupled-test'
OWN=['notes/'+STEM+'.md','checks/current_relative_angular_coupled_test.py',
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
    co={k:env[v] for k,v in zip(york.powers,york.qs)}
    def deriv(j):
        return {k-j:v*factorial(k)/factorial(k-j) for k,v in co.items() if k>=j}
    def bern_interval(d,lo,hi):
        n=max(d);a=arb(lo);width=arb(hi)-a
        shifted={i:sum((d[k]*comb(k,i)*a**(k-i) for k in d if k>=i),arb(0)) for i in range(n+1)}
        return [sum((shifted[i]*width**i*comb(k,i)/comb(n,i) for i in range(k+1)),arb(0)) for k in range(n+1)]
    qb=bern_interval(co,'.02','.03')
    qp=bern_interval(deriv(1),'.02','.03')
    qpp=bern_interval(deriv(2),'.02','.03')
    check('whole_actual_test_interval_Q_lower',all(v>arb('.16') for v in qb))
    check('whole_actual_test_interval_Q_upper',all(v<arb('.17') for v in qb))
    check('whole_actual_test_interval_Qprime',all(abs(v)<arb('.0018') for v in qp))
    check('whole_actual_test_interval_Qsecond_negative',all(v<-arb('.0775') for v in qpp))
    x=sp.symbols('x',real=True)
    w=sp.Function('w')(x); rr=sp.sqrt(1+w*w)
    B1=sp.diff(w,x)/(2*rr**3)
    tau2=B1**2/(2*rr)-w*sp.diff(B1,x)/(2*rr**2)
    alpha2=B1**2/2+sp.diff(B1,x)/(2*w*rr)
    v2=-(w*w+2)/rr**2*alpha2+tau2/rr
    check('relative_potential_subprincipal',sp.simplify(v2+B1**2/2+sp.diff(B1,x)/(w*rr))==0)
    d2=-tau2+w*w*v2/rr**3
    stated=-(1+2*w*w)*B1**2/(2*rr**3)+w*(w*w-1)*sp.diff(B1,x)/(2*rr**4)
    check('complete_tadpole_plus_bubble_subprincipal',sp.simplify(d2-stated)==0)
    ww,ll,lp=sp.symbols('w ell ellprime',real=True,positive=True)
    d2poly=sp.simplify(stated.subs({sp.diff(w,x,2):(lp+ll**2)*w,sp.diff(w,x):ll*w,w:ww},simultaneous=True))
    # Substitute in stages to ensure all functions are converted.
    d2poly=sp.simplify(stated.subs({sp.diff(w,x,2):(lp+ll**2)*w,sp.diff(w,x):ll*w}).subs(w,ww))
    a=-(1+2*ww**2)*ww**2/(8*(1+ww**2)**sp.Rational(9,2))-3*ww**4*(ww**2-1)/(4*(1+ww**2)**sp.Rational(9,2))
    b=ww**2*(ww**2-1)/(4*(1+ww**2)**sp.Rational(7,2))
    check('complete_integrand_decomposition',sp.simplify(d2poly-ll**2*a-(lp+ll**2)*b)==0)
    ia=sp.integrate(a,(ww,0,sp.oo));ib=sp.integrate(b,(ww,0,sp.oo))
    check('entire_scaled_thermal_first_integral',ia==-sp.Rational(37,420))
    check('entire_scaled_thermal_second_integral',ib==sp.Rational(1,60))
    check('full_thermal_local_coefficient',sp.expand(ll**2*ia+(lp+ll**2)*ib-lp/60+ll**2/14)==0)
    z=sp.symbols('z',positive=True);Q=sp.Function('Q')(z);f=sp.sqrt(Q)
    ell=-sp.diff(Q,z)/(2*f)
    ellprime=f*sp.diff(ell,z)
    caux=sp.simplify(f*(ellprime/60-ell**2/14))
    formula=-f*sp.diff(Q,z,2)/120-23*sp.diff(Q,z)**2/(1680*f)
    check('actual_geometry_conversion',sp.simplify(caux-formula)==0)
    source=sp.simplify(-2*caux/f)
    check('original_signed_source_density',sp.simplify(source-sp.diff(Q,z,2)/60-23*sp.diff(Q,z)**2/(840*Q))==0)
    cmin=arb('.16').sqrt()*arb('.0775')/120-23*arb('.0018')**2/(1680*arb('.16').sqrt())
    jmin=(arb('.0775')/60-23*arb('.0018')**2/(840*arb('.16')))/(arb.pi()*par['B2'])
    check('nonzero_auxiliary_thermal_limit_whole_interval',cmin>arb('2.58e-4'))
    check('nonzero_original_source_limit_whole_interval',jmin>arb('4.8e-47'))
    n,N,q,C=sp.symbols('n N q C',positive=True,integer=True)
    d=2*(q+2*n);mass2=n*(n+q)/C
    dn=sp.summation(d,(n,1,N))
    check('exact_all_angular_weight',sp.expand(dn-2*q*N-2*N*(N+1))==0)
    check('angular_weight_diverges',sp.limit(dn,N,sp.oo)==sp.oo)
    check('Mminus2_residual_is_borderline',sp.limit(n*d/mass2,n,sp.oo)==4*C)
    check('Mminus4_summable_decay',sp.limit(n**3*d/mass2**2,n,sp.oo)==4*C*C)
    # Established bound 4 C^2/q is conservative: compare decreasing F to integral.
    xx=sp.symbols('xx',positive=True)
    F=(q+2*xx)/(xx**2*(xx+q)**2)
    check('mass_fourth_integral_antiderivative',sp.simplify(sp.diff(-1/(xx*(xx+q)),xx)-F)==0)
    check('mass_fourth_decreasing_numerator',sp.simplify(-sp.diff(F,xx)*xx**3*(xx+q)**3-2*(q*q+3*q*xx+3*xx*xx))==0)
    # Sum_{n>=1}F(n) <= F(1)+int_1^infty F < 2/q.
    check('mass_fourth_complete_sum_guard',sp.simplify(2/q-(q+2)/(q+1)**2-1/(q+1)-(q+2)/(q*(q+1)**2))==0)
    sum4=4*par['C']**2/par['q']
    A4=2*sum4/(arb.pi()*par['B2'])
    check('original_composed_angular_factor',A4<arb('4.583e-35'))
    rem=arb('5.34235e-12')-arb('5.042e-23')
    actual_ceiling=rem/(A4*arb('1.1e-10'))
    affine_guard=arb('4.583e-35')*arb('1e33')*arb('1.1e-10')
    check('alternative_affine_guard_exact',sp.Rational('4.583e-35')*sp.Rational('1e33')*sp.Rational('1.1e-10')==sp.Rational('5.0413e-12'))
    check('alternative_other_affine_reserve',rem-affine_guard>arb('3.0104999994957e-13'))
    ge=arb('.0012069')/(arb('.0195239')+arb('4.583e-13')+arb('.04583'))
    check('alternative_feedback_sufficient_guard',arb('.01846')<ge)
    check('conditional_fast_scope_numerical_screen',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    check('actual_temperature_is_positive_and_fixed',0<par['gamma']<arb('1e-33'))
    return dict(outcome='unresolved_physical_coupled_estimate_with_verified_separate_signed_angular_estimator_failure',
      signed_auxiliary_limit='-sqrt(Q)*Qsecond/120-23*Qprime^2/(1680*sqrt(Q))',
      coefficient_lower_expression=enc(cmin),original_source_magnitude_lower_expression=enc(jmin),
      whole_interval=['.02','.03'],all_thermal_integrals=[str(ia),str(ib)],
      angular_weight='2qN+2N(N+1)',actual_sum4_upper=enc(sum4),original_composed_A4=enc(A4),
      conditional_C4_sole_affine_ceiling=enc(actual_ceiling),
      conditional_C4_selected_allowance='1e33',conditional_affine_remainder_guard=enc(affine_guard),
      conditional_other_affine_reserve=enc(rem-affine_guard),conditional_gamma_E_ceiling=enc(ge),
      actual_C4=None,actual_gamma_E=None,kappa_m=None,actual_mass=None,
      signed_separate_piece_estimator_failure=True,actual_coupled_weak_inverse_does_not_remove_this_failure=True,
      complete_thermal_asymptotic_executed=True,complete_angular_asymptotic_executed=True,
      physical_common_regulator_limit_executed=False,common_subtracted_strong_sum_certified=False,
      physical_coupled_signed_solve_executed=False,physical_operator_numerical_applications=0,
      actual_paired_domain_retained=True,finite_mode_selection=False,all_coincidence_and_companion_terms_required=True,
      actual_construction_obstruction=False,nonlinear_admission=False,nonlinear_contraction=False,
      exact_compatibility=False,physical_matching=False,semiclassical_validity_noise=False,all_sector_stability=False,
      full_objective_preserved=True,closed_audit_preserved=True)

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
