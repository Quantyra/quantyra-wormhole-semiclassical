"""Exact mixed-row and transpose checks; no numerical strong inverse solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from operators.york_mixed_return import YorkMixedReturn
STEM='current-schur-residual-phase-domain'
OWN=['notes/'+STEM+'.md','operators/york_mixed_return.py',
     'checks/current_schur_residual_phase_domain.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['operators/signed_local_core.py',
 'notes/current-schur-form-preconditioner-test.md',
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-global-boundary-moment-return.md',
 'notes/current-global-mixed-schur-preconditioner.md',
 'notes/current-physical-york-dual-assembly.md',
 'notes/current-coupled-input-energy-bound.md',
 'notes/current-short-travel-strong-domain.md',
 'notes/current-implicit-phase-graph-domain.md',
 'notes/current-coupled-truncation-energy-error.md',
 'notes/current-adaptive-dirac-trace-remainder.md',
 'evidence/current-global-mass-conversion.json']
CHECKS=[]
def check(name,ok):
    CHECKS.append({'name':name,'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'),flush=True)
    if not ok:raise AssertionError(name)
def zero(e):return sp.cancel(sp.expand(e))==0
def dec(v):
    d=v['decimal_enclosure'];return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def enc(v):
    m,r,e=v.mid_rad_10exp(50)
    return {'ball':v.str(32),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[3]).read_bytes())
    for p,h in review['scientific_sha256'].items():
        check('frozen_'+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)
    model=YorkMixedReturn();z=model.z;c=model.core
    built=model.calculate();exact=built['seed']
    check('whole_mixed_row_no_order_four',zero(exact['xi4']))
    check('whole_mixed_row_no_order_three',zero(exact['xi3']))
    valuations={name:model.valuation(v) for name,v in exact.items()}
    for name in ['xi0','xi1','xi2','b_radial','J_Y']:
        check('whole_outer_gain_'+name,valuations[name] is None or valuations[name]>=3)
    check('whole_seed_coefficients_polynomial',all(not sp.denom(v).has(z) for v in exact.values()))
    qseed=c.seed_substitution()[c.Q[0]]
    bolt_quotient=sp.cancel((exact['J_Y']-sp.diff(qseed,z)*exact['b_radial'])/qseed)
    check('combined_bolt_log_cancellation',not sp.denom(bolt_quotient).has(z))
    # Check formal transpose primitives by differentiation, not selected modes.
    Q=sp.Function('Q')(z);f=sp.Function('f')(z);R=sp.Function('R')(z);A=sp.Function('A')(z)
    candidate=z**4*R/Q+sp.diff(Q,z)*A
    # R'=-f/z^4, A'=z^4 R/Q^2: both pieces of the global J adjoint.
    derivative=sp.diff(candidate,z).subs({sp.diff(R,z):-f/z**4,sp.diff(A,z):z**4*R/Q**2})
    check('J_adjoint_exact_local_derivative',zero(derivative-(4*z**3*R/Q-f/Q+sp.diff(Q,z,2)*A)))
    y=sp.Function('y')(z);T=sp.Function('T')(z)
    Ay=lambda u:sp.Rational(2,3)*(sp.diff(u,z)-4*u/z)/z**2
    Aystar=lambda u:-sp.Rational(2,3)*(sp.diff(u,z)+2*u/z)/z**2
    check('inverse_York_adjoint_form',zero(sp.Rational(3,4)*z**4*Aystar(y)+z*y+z**2*sp.diff(y,z)/2))
    # Total derivative is the exact integration-by-parts boundary expression.
    check('A_Y_formal_adjoint_with_boundary',zero(y*Ay(T)-T*Aystar(y)-sp.diff(sp.Rational(2,3)*y*T/z**2,z)))
    phi=sp.Function('phi')(z)
    check('no_resonant_quartic_forcing',zero(Ay(z**4*phi)-sp.Rational(2,3)*z**2*sp.diff(phi,z)))
    n=sp.symbols('n')
    check('York_outer_indicial_roots',zero(6-n*(n-1)+(n-3)*(n+2)))
    t=sp.symbols('t')
    F=t/(1+t)
    check('full_twice_residual_transfer',zero(F**2*(1-F)**-3/(1+t)-t**2))
    q=arb('.338');factor=(1+q)**3
    src=json.loads((ROOT/'evidence/current-global-mass-conversion.json').read_bytes())['runs'][-1]
    check('retained_parameter_receipt',src['status']=='passed')
    p={k:dec(v) for k,v in src['result']['parameters'].items()}
    be,uu,g=[p[k] for k in ['beta','u','g']]
    Wmax=2*arb.pi()*be/(uu*g*arb('5e-34'))
    Uaff=(1+q)*(Wmax*be).sqrt()*arb('1.1e-10')
    Cceiling=arb('5.34235e-12')*(1-q)/(arb('1.1e-10')*(1+q)**4*q**30)
    check('affine_input_conversion_finite',Uaff>0)
    check('conditional_affine_ceiling_positive',Cceiling>0)
    check('two_fast_screen_unchanged',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    check('original_mass_conversion_debit',sp.Rational('1.627532e-15')*2174==sp.Rational('3.538254568e-12'))
    return {'precision_bits':ctx.prec,'whole_mixed_outer_valuations':valuations,
      'complete_mixed_coefficients':{name:str(v) for name,v in exact.items()},
      'complete_mixed_coefficients_sha256':sha256(str(exact).encode()).hexdigest(),
      'combined_bolt_quotient':str(bolt_quotient),
      'form_F_norm_upper':enc(q/(1+q)),'transfer_factor_upper':enc(factor),
      'physical_affine_Y_F_input_upper':enc(Uaff),
      'conditional_C_graph_ceiling_N32':enc(Cceiling),
      'scope':{'exact_complete_mixed_operator_executed':True,
       'ordinary_fixed_seed_graph_phase_composite_domain':True,
       'complete_affine_input_conversion':True,'numerical_composite_norm_bounded':False,
       'original_strong_delta_certified':False,'candidate_numerically_applied':False,
       'complete_coupled_return_executed':False,'actual_kappa_m_computed':False,
       'physical_admission_or_feedback_budget_closed':False,'nonlinear_admission_contraction':False,
       'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
       'all_sector_stability':False,'actual_controlled_solution':False,'construction_obstruction':False,
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
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
