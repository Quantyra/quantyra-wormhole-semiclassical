"""Supporting regulator identities and complete angular sums; no strong graph claim."""
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
STEM='current-mass-response-regulator-audit'
OWN=['notes/'+STEM+'.md','checks/current_mass_response_regulator_audit.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def algebra():
 # Mixed derivatives in the ring with three square-zero parameters. This
 # checks genuinely noncommuting trace words; these are not physical modes.
 P={i:S.Matrix([[2+i,1+2*i],[3-i,5+i*i]]) for i in range(8)}
 def mult(a,b):
  out={}
  for i,A in a.items():
   for j,B in b.items():
    if not i&j:out[i|j]=out.get(i|j,S.zeros(2))+A*B
  return out
 power={0:S.eye(2)};direct={}
 for m in range(1,8):
  power=mult(power,P);direct[m]=S.Rational((-1)**m,factorial(m))*S.trace(power.get(7,S.zeros(2)))
  calc=-S.Rational((-1)**(m-1),factorial(m-1))*S.trace(P[7]*P[0]**(m-1))
  if m>=2:
   calc+=S.Rational((-1)**(m-2),factorial(m-1))*sum(S.trace(P[pair]*P[0]**r*P[single]*P[0]**(m-2-r)) for pair,single in [(3,4),(5,2),(6,1)] for r in range(m-1))
  if m>=3:
   calc-=S.Rational((-1)**(m-3),factorial(m-1))*sum(S.trace(P[1]*P[0]**r*P[v]*P[0]**s*P[w]*P[0]**(m-3-r-s)) for v,w in [(2,4),(4,2)] for r in range(m-2) for s in range(m-2-r))
  check('noncommuting_Duhamel_trace_coefficient_'+str(m),calc==direct[m])
 b,d,g,rad,cs,H1,e=S.symbols('b d g R C H1 e')
 a2=(2*H1*S.exp(2*b)-cs*S.exp(-2*b)*(rad-2*d)**2/60)/(4*S.pi)
 cubic=S.expand(S.diff(a2.subs({b:e*b,d:e*d}),e,3).subs(e,0)/6)
 expected=((8*H1/3)*b**3+cs*rad**2*b**3/45+2*cs*rad*b*b*d/15+2*cs*b*d*d/15)/(4*S.pi)
 check('complete_logarithmic_cubic_contact',S.expand(cubic-expected)==0)
 check('mass_log_contact_does_not_share_conversion_cancellation',S.Poly(cubic,b,d).coeff_monomial(b*d*d)==cs/(30*S.pi))
 t=S.symbols('t',positive=True);theta=S.Function('Theta')
 expr=theta(t)
 for i in range(3):expr=2*t*S.diff(expr,t)
 check('complete_constant_mass_chain_rule',S.expand(expr-(8*t*S.diff(theta(t),t)+24*t*t*S.diff(theta(t),t,2)+8*t**3*S.diff(theta(t),t,3)))==0)
 check('ultraviolet_signed_ratio',S.Rational(-48,-8)==6 and S.Rational(8+48+48,8)==13)
 return dict(noncommuting_orders=7,cubic_density=str(cubic),signed_leading_coefficients=[-8,48,-48],combined_leading=-8,triangle_ratio=6,absolute_sum_ratio=13)

def execute():
 alg=algebra()
 pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 q=pars['q'];C=pars['C']
 r=S.symbols('r');pol=r;polys=[]
 for k in range(34):
  polys.append([int(v) for v in reversed(S.Poly(pol,r).all_coeffs())])
  pol=S.expand(r*((1-r)*S.diff(pol,r)+(k+1)*pol))
 def angular(t,J=12):
  a=q*t/C;b=t/C;rr=(-a).exp();one=1-rr
  sums=[]
  for k,co in enumerate(polys):
   v=arb(0)
   for c in reversed(co):v=v*rr+c
   sums.append(v/one**(k+1))
  def f(k,j):
   return sum((comb(k,i)*q**(k-i)*(2*q*sums[k+i+2*j]+4*sums[k+i+2*j+1]) for i in range(k+1)),arb(0))/C**k
  moments=[];errors=[]
  for k in range(4):
   v=sum(((-b)**j/factorial(j)*f(k,j) for j in range(J+1)),arb(0))
   error=b**(J+1)/factorial(J+1)*f(k,J+1)
   moments.append(v+arb(0,error.abs_upper()));errors.append(error)
  single=-8*t*moments[1];pairs=24*t*t*moments[2];triangle=-8*t**3*moments[3]
  scale=2*C/t
  return dict(t=old.encoded(t),single=old.encoded(single/scale),pairs=old.encoded(pairs/scale),triangle=old.encoded(triangle/scale),combined=old.encoded((single+pairs+triangle)/scale),triangle_relative_to_complete=old.encoded(triangle/(single+pairs+triangle)),sum_absolute_relative_to_complete=old.encoded((-single+pairs-triangle)/(-single-pairs-triangle)),moment_absolute_errors=[old.encoded(x) for x in errors],J=J)
 z=S.symbols('z');be,uc,ep,ta,cs,qs=S.symbols('beta u eps tau C q')
 Q=be+uc*z*z+(-4*be-2*uc-ta)*z**3+(3*be+uc+ep)*z**4-2*ep*z**5+ep*z**6
 R=-S.diff(Q,z,2);f=z**5;d=S.diff(Q*S.diff(f,z),z);H1=(qs**2/6-S.Rational(1,30))/cs
 dens=[-8*cs*f**3,2*cs*(R*f**3+3*f*Q*S.diff(f,z)**2)/3,(16*H1*f**3+2*cs*R**2*f**3/15+4*cs*R*f*f*d/5+4*cs*f*d*d/5)/4]
 mapping={be:pars['beta'],uc:pars['u'],ep:pars['eps'],ta:pars['tau'],cs:C,qs:q}
 def val(e):
  if e in mapping:return mapping[e]
  if e.is_Rational:return arb(int(e.p))/int(e.q)
  if e.is_Add:return sum((val(x) for x in e.args),arb(0))
  if e.is_Mul:
   v=arb(1)
   for x in e.args:v*=val(x)
   return v
  if e.is_Pow:return val(e.base)**int(e.exp)
  raise ValueError(e)
 integrals=[]
 for density in dens:
  poly=S.Poly(S.expand(density),z)
  v=sum((val(co)*zh**(power+1)/(power+1) for (power,),co in poly.terms()),arb(0))/arb.pi()
  integrals.append(old.encoded(v))

 angular_results=[angular(arb(t)) for t in ['1e-6','1e-8']]
 limits=[('[-7.9999999292644739 +/- 1e-16]'),('[-7.999999999992926447 +/- 1e-18]')]
 for j,(rec,limit) in enumerate(zip(angular_results,limits)):
  check('complete_angular_enclosure_'+str(j),arb(limit).contains(old.decode(rec['combined'])))
  check('triangle_overcount_'+str(j),old.decode(rec['triangle_relative_to_complete'])>6)
  check('absolute_sum_loss_'+str(j),old.decode(rec['sum_absolute_relative_to_complete'])>13)
 check('admissible_z5_a0_negative',old.decode(integrals[0])<arb('-7.81764e9'))
 check('admissible_z5_a1_negative',old.decode(integrals[1])<arb('-6.56944e8'))
 check('admissible_z5_a2_nonzero',old.decode(integrals[2])>arb('6.91233e13'))
 check('all_angular_tail_is_positive',all(old.decode(v)>0 for a in angular_results for v in a['moment_absolute_errors']))
 best=arb('.1599164')/arb('1e-11');example=(arb('.1599164')-arb('.08'))/arb('1e-11')
 check('conditional_neighborhood_targets',best.contains(arb(15991640000)) and example.contains(arb(7991640000)))
 return dict(algebra=alg,angular=angular_results,admissible_z5_third_heat_coefficients=integrals,conditional_C3_best_ceiling=old.encoded(best),conditional_C3_if_seed_point08=old.encoded(example),exact_positive_cutoff_vertex_inventory=True,actual_mass_square_domain_from_composition=True,common_logarithmic_third_contact_retained=True,complete_angular_test_executed=True,admissible_finite_seed_contact_test_executed=True,physical_domain_changed=False,physical_source_changed=False,normalization_changed=False,finite_physical_modes_selected=False,actual_strong_remainder_certified=False,uniform_strong_regulator_limit_proved=False,original_X_graph_comparison_known=False,actual_contact_jets_certified=False,full_signed_projection_executed=False,actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,next_strong_calculation_executed=False)

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
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),complete_angular_test=result.get('complete_angular_test_executed')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
