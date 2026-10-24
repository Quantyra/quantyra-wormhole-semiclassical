"""Uniform principal-coefficient test and correction of the selected projection status."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
from functools import lru_cache
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-uniform-principal-coefficient'
OWN=['notes/'+STEM+'.md','checks/current_uniform_principal_coefficient.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def execute():
 ctx.prec=512
 prior=json.loads((ROOT/'evidence/current-signed-resolvent-projection.json').read_bytes())['runs'][-1]
 independent=json.loads((ROOT/'evidence/current-signed-resolvent-projection-closeout-verification.json').read_bytes())['attempts'][-1]
 check('prior_selected_screen_executed',prior['status']=='passed' and prior['complete_signed_projection_evaluated'] and prior['full_projection_milestone_complete'])
 check('sign_was_not_measured',not prior['signed_projection_value_measured'])
 check('prior_independent_verification_passed',independent['status']=='passed')
 check('prior_all_recorded_inputs_unchanged',all(digest(f)==h for f,h in prior['input_sha256_after'].items()))
 certs=[v for v in prior['numerics'] if v['kind']=='whole_signed_projection_certificate'];totals={}
 for row in certs:
  total=sum((old.decode(v) for v in row['errors'].values()),arb(0));totals[str(row['N'])]=old.encoded(total)
  check('complete_error_inventory_'+str(row['N']),len(row['errors'])==13 and total.overlaps(old.decode(row['E_rest'])))
  check('sufficient_screen_'+str(row['N']),total<arb('1.287e-40') and row['P_hat_rest']=='0' and row['complete_signed_interval'])
 check('retained_rounding_guard',old.decode(totals['11'])<arb('2.123e-42'))
 pars,h,zh,env=old.scalar.prev.york.actual_parameters();sig=arb('4.005');a=arb('.01');A=arb('.001');beta=pars['beta'];u=pars['u'];ep=pars['eps'];tau=pars['tau']
 q=[beta,arb(0),u,-4*beta-2*u-tau,3*beta+u+ep,-2*ep,ep];q=[v/beta for v in q]
 poly=[sum((q[i]*q[k-i] for i in range(max(0,k-6),min(6,k)+1)),arb(0)) for k in range(13)]
 R=[sum((abs(poly[k])*k**j*a**k for k in range(1,13)),arb(0)) for j in range(4)]
 B=[A*3**j*a**3/6 for j in range(4)];ex=(2*B[0]).exp()
 E=[ex-1,2*ex*B[1],ex*(2*B[2]+4*B[1]**2),ex*(2*B[3]+12*B[1]*B[2]+8*B[1]**3)]
 G=[2*B[0]**2*ex,2*(ex-1)*B[1],2*(ex-1)*B[2]+4*ex*B[1]**2,2*(ex-1)*B[3]+12*ex*B[1]*B[2]+8*ex*B[1]**3]
 H=[sum((comb(j,k)*R[k]*E[j-k] for k in range(j+1)),arb(0))+G[j] for j in range(4)]
 weights=[100,51,12,1];cq=arb('2.328');eta=cq*sum((weights[j]*H[j] for j in range(4)),arb(0))
 V=[arb(1801)/arb(2).sqrt()*(sig+1)**j*a**sig for j in range(4)]
 AA=[sum((comb(j,k)*R[k]*([ex]+E[1:])[j-k] for k in range(j+1)),arb(0))+E[j] for j in range(4)]
 JJ=[2*sum((comb(j,k)*AA[k]*V[j-k] for k in range(j+1)),arb(0)) for j in range(4)]
 derivative=cq*sum((weights[j]*JJ[j] for j in range(4)),arb(0))
 x=S.symbols('x');bb=S.Function('b')(x);rr=S.symbols('rho');vv=S.symbols('v')
 check('principal_repartition_exact',S.simplify((1+rr)*(S.exp(-2*bb)-1)+2*bb-(rr*(S.exp(-2*bb)-1)+S.exp(-2*bb)-1+2*bb))==0)
 check('coefficient_derivative_exact',S.simplify(S.diff((1+rr)*(S.exp(-2*x)-1)+2*x,x)+2*((1+rr)*S.exp(-2*x)-1))==0)
 z=S.symbols('z',positive=True);f=S.Function('f')(z);expr=f
 for j in range(4):
  # For a cubic controlled by its third derivative, all Euler coefficients
  # through three derivatives are nonnegative and their cubic action is 3^j.
  if j:expr=S.expand(z*S.diff(expr,z))
  cube=expr.subs(f,z**3/6).doit()
  check('mass_phase_Euler_bound_'+str(j),S.simplify(cube-3**j*z**3/6)==0)
 cutoff=[1,3,10,60]
 generated=[sum(comb(3,l)*comb(l,j)*cutoff[l-j] for l in range(j,4)) for j in range(4)]
 check('cutoff_multiplier_weights',generated==weights)
 check('principal_family_guard',eta>0 and eta<arb('2e-11'))
 check('fast_coefficient_orientation_guard',derivative>0 and derivative<arb('8.39e-7'))
 check('phase_suppression_not_normalization',B[0]<arb('1.667e-10') and R[0]<arb('8.061e-6'))
 remaining=arb('.030416395999')-arb('2e-11');target=remaining/arb('1e-11')
 check('remaining_budget_positive',remaining>arb('.0304163959789'))
 check('other_C3_threshold',target>arb('3.0416395978e9') and target<arb('3.0416395980e9'))
 oldbolt=json.loads((ROOT/'evidence/current-localized-bolt-adjoint-sensitivity.json').read_bytes())['runs'][-1]
 check('localized_bolt_result_retained',oldbolt['status']=='passed' and all(digest(f)==h for f,h in oldbolt['input_sha256_after'].items()))
 bolt=oldbolt['result']['necessary_projection_ceiling']
 return dict(prior_selected_screen_inputs=len(prior['input_sha256_after']),prior_screen_totals=totals,prior_selected_screen_revalidated=True,complete_selected_signed_interval_certified=True,projection_newly_evaluated=False,sign_measured=False,prior_projection_status_correction=True,old_uniform_bolt_projection_ceiling=bolt,curvature_coefficient_caps=[old.encoded(v) for v in R],background_Euler_caps=[old.encoded(v) for v in B],exponential_quadratic_caps=[old.encoded(v) for v in G],complete_principal_defect_caps=[old.encoded(v) for v in H],principal_defect_upper=old.encoded(eta),coefficient_derivative_orientation_upper=old.encoded(derivative),principal_debit='2e-11',remaining_budget=old.encoded(remaining),best_case_other_C3_ceiling=old.encoded(target),equal_split_seed_ceiling=old.encoded(remaining/2),equal_split_other_C3_ceiling=old.encoded(target/2),principal_family_theorem_executed=True,full_third_response_bound=False,full_actual_common_strong_regulator_limit=False,physical_domain_changed=False,normalization_changed=False,global_admission=False,actual_controlled_solution=False,nonlinear_closure=False,exact_compatibility=False,physical_matching=False,semiclassical_validity_noise=False,all_sector_stability=False)

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
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),principal_defect=result.get('principal_defect_upper'),coefficient_derivative=result.get('coefficient_derivative_orientation_upper'),screen_revalidated=result.get('prior_selected_screen_revalidated')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
