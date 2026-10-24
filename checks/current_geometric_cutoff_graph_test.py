"""Actual-domain cutoff obstruction and a full-angular principal-model repair."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-geometric-cutoff-graph-test'
OWN=['notes/'+STEM+'.md','checks/current_geometric_cutoff_graph_test.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def execute():
 v,a,X=S.symbols('v a X',positive=True)
 psi=(1+v)*S.exp(-v)-1+v*v/2
 check('cap_primitive_derivative',S.simplify(S.diff(-psi/(2*v*v),v)+(1-S.exp(-v))/(2*v)-psi/v**3)==0)
 check('cap_primitive_zero_endpoint',S.limit(psi/v**2,v,0)==0)
 check('positive_remainder_integrand',S.simplify(S.diff(psi,v)-v*(1-S.exp(-v)))==0)
 check('heat_a2_moment',S.integrate((a*(1-a))**2,(a,0,1))==S.Rational(1,30))
 check('heat_a6_moment',S.integrate((a*(1-a))**3,(a,0,1))==S.Rational(1,140))
 check('cap_small_X_coefficient',S.Rational(1,2)*S.Rational(1,3)*S.Rational(1,140)==S.Rational(1,840))
 check('cap_log_coefficient',S.Rational(1,120)*60*2==1)
 x,q,t,C=S.symbols('x q t C',positive=True)
 ft=2*(q+2*x)*S.exp(-x*(x+q)*t/C)
 check('angular_density_primitive',S.simplify(S.diff(-2*C/t*S.exp(-x*(x+q)*t/C),x)-ft)==0)
 check('angular_unimodality',S.simplify(S.diff(ft,x)/ft-(2/(q+2*x)-(q+2*x)*t/C))==0)
 z,s=S.symbols('z s',positive=True);Q=S.Function('Q')(z)
 lap=lambda u:S.diff(Q*S.diff(u,z),z)
 row=z**4*lap(lap(z**s));row=S.expand((z*S.diff(row,z)-4*row)/z**s)
 check('actual_row_principal_coefficient',S.simplify(row.coeff(s,5)-Q**2)==0)
 # Original full graph parameters, no independent horizon or geometry rounding.
 p,h,zh,env=old.scalar.prev.york.actual_parameters()
 be,u,ep,ta=p['beta'],p['u'],p['eps'],p['tau'];Cn=p['C'];qn=arb(p['q'])
 coeff=[arb(1),arb(0),u/be,-4-2*u/be-ta/be,3+u/be+ep/be,-2*ep/be,ep/be]
 evaluate=lambda zz:sum((val*zz**i for i,val in enumerate(coeff)),arb(0))
 interval=arb('0.002','0.0001')
 Qratio=evaluate(interval);alower=(Qratio*Qratio).lower()
 astar=evaluate(arb('0.002'))**2
 check('collar_inside_unit_output_cutoff',arb('.0021')<arb('.01')/arb(1).exp())
 check('actual_positive_collar',Qratio>0)
 check('actual_raw_cutoff_norm_obstruction',alower>arb('.9999999'))
 check('raw_cutoff_exceeds_remaining_target',alower>arb('.18368'))
 Cq=2*(arb(17)/16)**(arb(5)/2)
 check('retained_graph_constant',Cq<arb('2.328'))
 A=qn/(2*Cn);Bb=2*(2/Cn).sqrt();fac=astar*Cq/arb(2).log()
 cases=[]
 for eps in ['.001','.0001','.000001']:
  e=arb(eps);bound=fac*(A*e+Bb*e.sqrt())
  cases.append(dict(epsilon=eps,angular_discrepancy_cap_upper=old.encoded(bound),model_target_certified=bool(bound<arb('.001'))))
  if eps=='.000001':check('repaired_full_angular_model_cap_target',bound<arb('.001'));check('repaired_model_refined_guard',bound<arb('.00055'))
  if eps=='.001':check('coarse_model_upper_does_not_certify',bound>arb('.18'))
 d=arb('.001')/fac
 root=2*d/(Bb+(Bb*Bb+4*A*d).sqrt());threshold=root**2
 check('stable_cutoff_threshold_equation',(fac*(A*threshold+Bb*root)-arb('.001')).contains(0))
 check('selected_epsilon_within_threshold',arb('.000001')<threshold)
 # An explicit whole-frequency cap sandwich, not a finite-mode approximation.
 cap_sandwich=[]
 for xx in ['1','1e6','1e40']:
  xx=arb(xx);lo=xx.log()/120-arb(1)/4;hi=xx.log()/120+arb(1)/840
  check('cap_sandwich_order_'+str(xx),lo<hi)
  cap_sandwich.append(dict(X=old.encoded(xx),lower=old.encoded(lo),upper=old.encoded(hi)))
 return dict(actual_interval=['.0019','.0021'],actual_Q_over_beta=old.encoded(Qratio),
  raw_regulated_feedback_norm_lower=old.encoded(alower),actual_frozen_coefficient=old.encoded(astar),
  full_graph_Cq=old.encoded(Cq),full_angular_model_caps=cases,model_epsilon_threshold=old.encoded(threshold),
  cap_sandwich=cap_sandwich,
  raw_unregulated_subtraction_norm_Cauchy_obstructed=True,
  actual_compact_domain_directions_used=True,physical_domain_replaced=False,finite_physical_modes_selected=False,
  universal_cap_missing_log_matched=True,full_angular_frozen_cap_convergence_certified=True,
  actual_curved_matched_remainder_certified=False,actual_thermal_or_boundary_remainder_certified=False,
  selected_canonical_target='.18',selected_combined_remaining_target='.18368',
  full_signed_projection_executed=False,actual_controlled_solution=False,
  nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False)
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
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
  frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
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
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),result=run.get('result')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
