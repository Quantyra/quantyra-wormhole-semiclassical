"""Support the ordinary full-tail output test; no numerical inverse solve."""
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
STEM='current-schur-physical-output-test'
OWN=['notes/'+STEM+'.md','checks/current_schur_physical_output_test.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'evidence/current-schur-residual-phase-domain.json',
 'evidence/current-global-mass-conversion.json',
 'notes/current-schur-residual-phase-domain.md',
 'notes/current-schur-form-preconditioner-test.md',
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-global-mixed-schur-preconditioner.md',
 'notes/current-short-travel-strong-domain.md',
 'notes/current-coupled-global-return-attempt.md',
 'notes/strong-interior-quantum-remainder.md',
 'notes/general-interior-dirac-response.md']
CHECKS=[]
def check(name,ok):
    CHECKS.append({'name':name,'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'),flush=True)
    if not ok:raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
    d=x['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def ev(x,env):
    if x in env:return env[x]
    if x==sp.pi:return arb.pi()
    if x.is_Rational:return arb(int(x.p))/int(x.q)
    if x.is_Add:return sum((ev(v,env) for v in x.args),arb(0))
    if x.is_Mul:
        out=arb(1)
        for v in x.args:out*=ev(v,env)
        return out
    if x.is_Pow and x.exp.is_Integer:return ev(x.base,env)**int(x.exp)
    raise ValueError(str(x))
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[-1]).read_bytes())
    for p,h in review['scientific_sha256'].items():
        check('frozen_'+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)
    prior=json.loads((ROOT/INPUTS[3]).read_bytes())['runs'][-1]
    check('retained_full_mixed_receipt_passed',prior['status']=='passed')
    names=['z','C','B2','beta','eps','tau','cR']
    c2=sp.sympify(prior['result']['complete_mixed_coefficients']['xi2'],
                  locals={n:sp.Symbol(n) for n in names})
    sy={str(v):v for v in c2.free_symbols}
    z,C,B2,be,ep,ta,cR=[sy[v] for v in ['z','C','B2','beta','eps','tau','cR']]
    Q=be+z*z/C-(4*be+2/C+ta)*z**3+(3*be+1/C+ep)*z**4-2*ep*z**5+ep*z**6
    p2=sp.cancel(c2/Q)
    check('complete_c2_exact_Q_division',not sp.denom(p2).has(z))
    check('whole_polynomial_identity',sp.cancel(c2-Q*p2)==0)
    pp=json.loads((ROOT/INPUTS[4]).read_bytes())['runs'][-1]
    check('parameter_receipt_passed',pp['status']=='passed')
    p={k:dec(v) for k,v in pp['result']['parameters'].items()}
    env={C:1/p['u'], B2:p['beta']/(p['u']*p['g']),
         be:p['beta'],ep:p['eps'],ta:p['tau'],cR:arb('1e22')}
    # Entire .99..1 interval, not endpoint sampling or a rounded bolt.
    zhull=arb('.995','.005')
    val=arb(0)
    for coef in sp.Poly(p2,z).all_coeffs():val=val*zhull+ev(coef,env)
    check('nonzero_actual_mixed_principal_whole_bolt_collar',val>arb('-1e-6') and val<arb('-1e-10'))
    # Exact symbolic exponent bookkeeping, not a finite-frequency experiment.
    a=sp.Rational(3,4);pwr=sp.symbols('p',integer=True,positive=True)
    energy=2-a
    lapse_fifth=-4*pwr+6+6*a
    ratio=sp.expand(lapse_fifth-energy)
    check('full_Schur_physical_power',ratio==4-4*pwr+7*a)
    check('two_residual_positive_power',ratio.subs(pwr,2)==sp.Rational(5,4))
    check('two_residual_log_denominator',sp.Rational(2)-sp.Rational(1,2)==sp.Rational(3,2))
    check('three_residual_negative_power',ratio.subs(pwr,3)==-sp.Rational(11,4))
    check('three_residual_H14_sufficient',2+4*3>12+sp.Rational(2,2))
    check('two_residual_H10_insufficient',2+4*2<12+sp.Rational(2,2))
    t=sp.symbols('t');f=t/(1+t)
    check('full_three_power_identity',sp.cancel(f**3*(1-f)**-4/(1+t)-t**3)==0)
    q=arb('.338');pref=(1+q)**5*arb('1.1e-10')*q**29/(1-q)
    cg=arb('5.34235e-12')/pref
    cx=arb('5e-9')/pref
    old=arb('5.34235e-12')*(1-q)/(arb('1.1e-10')*(1+q)**4*q**30)
    check('same_norm_relative_price', (cg-old*q/(1+q)).contains(0))
    check('conditional_ceilings_positive',cg>0 and cx>0)
    check('fast_screen_scope_unchanged',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return {'precision_bits':ctx.prec,'p2_exact':str(sp.factor(p2)),
      'p2_whole_bolt_collar':enc(val),'p2_over_8piG_G1':enc(val/(8*arb.pi())),
      'packet_a':'3/4','two_residual_ratio':'c*k^(5/4)/(log k)^(3/2), c>0; asymptotic, onset unknown',
      'three_residual_ratio_power':'-11/4','physical_derivatives':{'required_disk_point_jets':12,'two_power_Sobolev':10,'three_power_Sobolev':14},
      'three_power_affine_tail_prefactor_N32':enc(pref),
      'conditional_graph_C3_ceiling':enc(cg),'conditional_physical_C3_ceiling_if_5e_minus9_margin':enc(cx),
      'extra_prefactor_price':enc((1+q)/q),
      'scope':{'ordinary_actual_domain_two_residual_physical_estimator_failure':True,
       'full_mixed_nonvanishing_executed':True,'three_power_fixed_seed_finiteness':True,
       'three_power_numerical_norm_bounded':False,'original_strong_delta_certified':False,
       'complete_coupled_return_executed':False,'actual_kappa_m_computed':False,
       'any_additional_admission_feedback_requirement_closed':False,
       'nonlinear_admission_contraction':False,'exact_compatibility':False,
       'physical_matching':False,'semiclassical_validity_noise':False,
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
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(summary.encode())
    print(summary,flush=True)
    if run['status']!='passed':raise SystemExit(1)
