"""Quantified fast-sector graph admission; actual strong data remain unknown."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-fast-graph-domain-embedding'
OWN=['notes/'+STEM+'.md','checks/current_fast_graph_domain_embedding.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def execute():
 pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 sigma=arb('4.005');delta=arb('.005');a=arb('.01');E=arb('1e-11')
 be=pars['beta'];r=(96*arb.pi()*arb('1e22')-1/(60*arb.pi()))*pars['g']
 Cv=1+11*pars['g']/(180*arb.pi());A=Cv-r*sigma*(sigma-3)
 check('positive_scalar_A',A>arb('.999999'))
 candidates=[2*arb(2).sqrt()/(6*be*A*delta*(sigma+1)),200*arb(2).sqrt()/(6*be*A),
             (1+arb('1e-4'))**(arb(3)/2)/(12*be*(A*r).sqrt())]
 Ctheta=max(v.upper() for v in candidates);Cb=arb(1801);aa=(sigma+3)/(3*(sigma+1));Cnu=Ctheta+aa*Cb
 ck=[]
 for k,expected in enumerate([63,7,3,3,7,63]):
  v=S.Rational(factorial(2*k)*factorial(10-2*k),2*4**5*factorial(k)*factorial(5-k)*factorial(5))
  check('exact_sobolev_constant_'+str(k),v==S.Rational(expected,512));ck.append((arb(expected)/512).sqrt())
 z,s,d=S.symbols('z s d');q=S.symbols('q0:7');Q=sum(q[i]*z**i for i in range(7))
 values=[be,arb(0),pars['u'],-4*be-2*pars['u']-pars['tau'],3*be+pars['u']+pars['eps'],-2*pars['eps'],pars['eps']]
 mapping=dict(zip(q,values))
 def val(v):
  if v in mapping:return mapping[v]
  if v.is_Rational:return arb(int(v.p))/int(v.q)
  if v.is_Add:return sum((val(w) for w in v.args),arb(0))
  if v.is_Mul:
   out=arb(1)
   for w in v.args:out*=val(w)
   return out
  if v.is_Pow:return val(v.base)**int(v.exp)
  raise ValueError(v)
 def bern(expr,right):
  poly=S.Poly(S.expand(expr),z)
  if poly.is_zero:return [arb(0)]
  n=poly.degree();power=[val(poly.nth(i))*right**i for i in range(n+1)]
  return [sum((power[i]*comb(k,i)/comb(n,i) for i in range(k+1)),arb(0)) for k in range(n+1)]
 qcoll=bern(Q,a);qmin=min(v.lower() for v in qcoll);qmax=max(v.upper() for v in qcoll)
 qwhole=bern(Q,zh);qg=max(v.upper() for v in qwhole)
 check('actual_collar_positive',qmin>arb('.16'))
 check('global_Q_guard',qg<arb('.2'))
 check('outer_proper_collar_inside_phase',a/qmin.sqrt()<arb('.1'))
 P=S.Integer(1);symbols=[];outer=[];global3=[];bounds=[]
 for j in range(6):
  symbols.append(str(P));poly=S.Poly(S.expand(P.subs(s,d+S.Rational(801,200))),d)
  local=arb(0);glob=arb(0);items=[]
  for k in range(j+1):
   expr=poly.nth(k);outb=bern(expr,a);mb=max(v.abs_upper() for v in outb)
   local+=mb*ck[k]
   item=dict(k=k,symbol=str(expr),outer_coefficient_upper=old.encoded(mb))
   if j<=3:
    gb=max(v.abs_upper() for v in bern(expr,zh));glob+=gb*ck[k];item['global_coefficient_upper']=old.encoded(gb)
   items.append(item)
  if j%2:local*=qmax.sqrt();glob*=arb('.2').sqrt()
  T=qmax**(arb(4-j)/2) if j<=4 else 1/qmin.sqrt()
  bound=local*a**delta*T/(factorial(j)*arb(100)**j*(1+(qmin.sqrt()/a).log()))
  outer.append(bound)
  if j<=3:global3.append(glob)
  bounds.append(dict(order=j,Lj_outer=old.encoded(local),outer_phase_multiplier=old.encoded(bound),terms=items))
  if j<5:
   core=(s-j)*P+z*S.diff(P,z)
   P=S.expand(core if j%2==0 else Q*core+z*S.diff(Q,z)*P/2)
 A6=max(v.upper() for v in outer);W6=ck[0]+arb('.2').sqrt()*(sigma*ck[0]+ck[1])
 check('value_controls_outer_phase_guard',all(outer[0].lower()>v.upper() for v in outer[1:]))
 components=[]
 for name,C in [('b',Cb),('theta',Ctheta),('nu',Cnu)]:
  outer_bound=A6*C*E;w1_bound=W6*C*E
  check('fast_outer_admission_'+name,outer_bound<arb('1e-8'))
  check('fast_W1_admission_'+name,w1_bound<arb('1e-3'))
  components.append(dict(component=name,H6_constant=old.encoded(C),conditional_outer_phase_upper=old.encoded(outer_bound),
   conditional_W1_upper=old.encoded(w1_bound),conditional_W3_jets=[old.encoded(v*C*E) for v in global3],
   remaining_asym_outer_allowance=old.encoded(arb('1e-8')-outer_bound),remaining_asym_W1_allowance=old.encoded(arb('1e-3')-w1_bound)))
 Eout=arb('1e-8')/(A6*Cnu);Ew1=arb('1e-3')/(W6*Cnu)
 check('selected_E_below_outer_threshold',E<Eout)
 check('selected_E_below_multiplier_threshold',E<Ew1)
 return dict(sigma='4.005',data_target='1e-11',actual_data_target_certified=False,
  scalar_coercivity_A=old.encoded(A),scalar_candidates=[old.encoded(v) for v in candidates],
  triangular_multiplier=old.encoded(aa),Q_outer_lower=old.encoded(qmin),Q_outer_upper=old.encoded(qmax),Q_global_upper=old.encoded(qg),
  Sobolev_constants=[old.encoded(v) for v in ck],proper_derivative_polynomials=symbols,proper_derivative_bounds=bounds,
  outer_phase_embedding=old.encoded(A6),W1_embedding=old.encoded(W6),components=components,
  sufficient_fast_E_outer=old.encoded(Eout),sufficient_fast_E_W1=old.encoded(Ew1),
  uniform_fast_graph_coefficient_admission_proved=True,weak_C2_action_applies_to_fast_graph_family=True,
  actual_mass_slow_allowances_certified=False,full_X_ball_admission_proved=False,
  actual_matched_curved_remainder_certified=False,physical_source_or_norm_changed=False,physical_boundary_changed=False,
  full_signed_projection_executed=False,actual_controlled_solution=False,nonlinear_preserved_ball=False,exact_compatibility=False,
  physical_matching=False,validity_noise=False,all_sector_stability=False)
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
 result=run.get('result',{})
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),summary={k:result[k] for k in ['scalar_candidates','outer_phase_embedding','W1_embedding','components','sufficient_fast_E_outer','sufficient_fast_E_W1'] if k in result}),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
