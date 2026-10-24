"""Support the ordinary bolt clock-gap counterexample; no quantum response bound."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-bolt-frequency-domain-test'
OWN=['notes/'+STEM+'.md','checks/current_bolt_frequency_domain_test.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def execute():
 pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 z,be,u,ep,ta=S.symbols('z beta u ep tau',real=True)
 Q=be+u*z**2+(-4*be-2*u-ta)*z**3+(3*be+u+ep)*z**4-2*ep*z**5+ep*z**6
 expected=2*u+(-24*be-12*u-6*ta)*z+(36*be+12*u+12*ep)*z**2-40*ep*z**3+30*ep*z**4
 check('exact_seed_second_derivative',S.expand(S.diff(Q,z,2)-expected)==0)
 x,gam=S.symbols('x gamma',positive=True)
 check('Taylor_relative_error_in_layer',S.simplify((6*x*x)/(gam*x)).subs(x,gam/6)==1)
 check('positive_bolt_margin',zh>arb('.999') and zh<1 and h>0)
 gamma=pars['gamma'];M=((1+pars['q'])/pars['C']).sqrt();B=arb('.001')
 Q2=60*pars['beta']+26*pars['u']+82*pars['eps'].abs_upper()+6*pars['tau'].abs_upper()
 check('global_Q_second_derivative_guard',Q2<arb('10.000000000326') and Q2<12)
 delta=gamma/6;p0=gamma/2;eps=p0/100;V0=M*B.exp()*gamma/arb(2).sqrt()
 length=2*arb.pi()/eps;start=2/gamma*(2*V0/eps).log();T=2*arb.pi()/gamma
 check('entire_bolt_layer_is_in_domain',delta>0 and delta<zh and delta<arb('8.430e-35'))
 check('actual_gamma_M1_boxes',gamma>arb('5.057e-34') and gamma<arb('5.059e-34') and M>18 and M<19)
 t=S.symbols('t',real=True);chi=S.sqrt(2)*S.sin(S.pi*t)
 check('packet_exact_norm',S.integrate(chi*chi,(t,0,1))==1)
 check('packet_exact_derivative_norm',S.integrate(S.diff(chi,t)**2,(t,0,1))==S.pi**2)
 check('packet_extension_zero_traces',chi.subs(t,0)==0 and chi.subs(t,1)==0)
 s1=S.Matrix([[0,1],[1,0]]);s2=S.Matrix([[0,-S.I],[S.I,0]]);s3=S.diag(1,-1)
 check('transfer_to_clock_potential_identity',S.I*s1*s3==s2)
 ev=S.Matrix([1,S.I])/S.sqrt(2)
 check('packet_clock_eigenvector',s2*ev==ev and (ev.conjugate().T*ev)[0]==1)
 ee,vv,gg=S.symbols('epsilon V0 gamma',positive=True)
 RR=2/gg*S.log(2*vv/ee);LL=2*S.pi/ee
 check('exact_packet_residual_budget',S.simplify(S.pi/LL+vv*S.exp(-gg*RR/2)-ee)==0)
 check('positive_optical_packet_start',start>0)
 check('packet_scales',length<arb('2.485e36') and start<arb('3.377e34') and gamma*start<arb('17.078'))
 check('optical_scaled_length_identity',(gamma*length).overlaps(400*arb.pi()))
 check('small_clock_residual',eps<arb('2.529e-36') and eps/(arb('.01')*M)<arb('1.402e-35'))
 check('global_clock_resolvent_conditional_scale',1/p0>arb('3.954e33'))
 check('unsupported_exponent_size',arb('.01')*M*T>arb('2.242e33'))
 check('remaining_C3_ceiling',S.Rational('.159916396')/S.Rational('1e-11')==15991639600)
 check('conditional_return_budget',S.Rational('.01')/S.Rational('1e-11')==10**9 and S.Rational('.159916396')-S.Rational('.01')==S.Rational('.149916396'))
 check('original_ball_allocation',S.Rational('5e-9')+S.Rational('1e-8')/4+S.Rational('5e-17')*S.Rational('5e7')==S.Rational('1e-8'))
 check('original_contraction_allocation',S.Rational(1,4)+S.Rational('1e-8')*S.Rational('5e7')==S.Rational(3,4))
 vals=dict(gamma=gamma,M1=M,Q2_guard=Q2,bolt_layer_delta=delta,p0=p0,packet_residual_target=eps,packet_length=length,packet_start=start,optical_scaled_start=gamma*start,optical_scaled_length=gamma*length,potential_bound_at_start=V0*(-gamma*start/2).exp(),derivative_bound=arb.pi()/length,nominal_gap=arb('.01')*M,residual_over_nominal_gap=eps/(arb('.01')*M),period=T,unsupported_thermal_exponent=arb('.01')*M*T,clock_resolvent_at_first_frequency=1/p0)
 return dict(numerics={k:old.encoded(v) for k,v in vals.items()},continuum_norm_gap_counterexample_executed=True,all_angular_masses_covered_by_argument=True,finite_physical_mode_approximation=False,spectral_conclusion_conditional_on_self_adjoint_clock_realization=True,physical_causal_realization_selected=False,Euclidean_mass_bound_retracted=False,global_gap_method_obstructed=True,local_thermal_response_obstructed=False,local_thermal_bound_known=False,remaining_seed_nonlocal_target='.159916396',complete_C3_ceiling_if_eta_zero='15991639600',proposed_return_C3_target='1e9',return_C3_known=False,proposed_return_calculation_executed=False,common_strong_regulator_limit_proved=False,full_signed_projection_executed=False,actual_curved_thermal_remainder_certified=False,uniform_complete_metric_variation_certified=False,original_X_graph_comparison_known=False,physical_domain_changed=False,physical_source_changed=False,normalization_changed=False,actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False)

def main():
 ctx.prec=512;review=json.loads((ROOT/OWN[2]).read_bytes());paths=set(OWN+review['retained_inputs'])
 for mod in list(sys.modules.values()):
  f=getattr(mod,'__file__',None)
  if f:
   try:p=Path(f).resolve().relative_to(ROOT).as_posix()
   except ValueError:continue
   if p.endswith('.py'):paths.add(p)
 before={p:digest(p) for p in sorted(paths)};out=ROOT/'evidence'/f'{STEM}.json'
 archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256']);run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 result=run.get('result',{})
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),clock_gap_test=result.get('continuum_norm_gap_counterexample_executed')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
