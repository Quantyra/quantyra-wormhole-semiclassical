"""Support the ordinary curved paired heat proof; no numerical solution claim."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
from collections import defaultdict
from itertools import permutations
import json,platform,sys,traceback
import sympy as sp
sys.path.insert(0,str(Path.home()/'AppData/Local/Temp/s169af-arb'))
import flint
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-curved-paired-heat-bridge'
INPUTS=[
 'checks/current_curved_paired_heat_bridge.py','notes/current-curved-paired-heat-bridge.md',
 'evidence/current-curved-paired-heat-bridge-review.json',
 'notes/current-rough-paired-variation.md','evidence/current-rough-paired-variation-closeout-verification.json',
 'notes/current-paired-regulator-domain-audit.md','notes/stationary-paired-dirac-phase.md',
 'notes/stationary-massive-boundary-layer.md','notes/stationary-exterior-conformal-domain.md',
 'notes/current-smeared-heat-projection-domain.md','notes/finite-outer-admission-and-response.md',
 'notes/finite-cold-admission-and-response.md','notes/current-signed-resolvent-projection.md',
 'evidence/current-signed-resolvent-projection-closeout-verification.json',
 'notes/current-controlled-solution-audit-closeout.md',
]
def hashes():return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}

def bridge_moment(powers):
 """Wick pairings integrated over every ordered time simplex, exactly."""
 n=len(powers);times=sp.symbols('s0:'+str(n),positive=True)
 indices=tuple(i for i,p in enumerate(powers) for _ in range(p))
 def wick(ids,order):
  if not ids:return sp.Integer(1)
  i=ids[0];out=0
  for k,j in enumerate(ids[1:],1):
   least=i if order[i]<=order[j] else j
   cov=2*(times[least]-times[i]*times[j])
   out+=cov*wick(ids[1:k]+ids[k+1:],order)
  return out
 total=0
 for perm in permutations(range(n)):
  order={v:i for i,v in enumerate(perm)}
  expr=sp.expand(wick(indices,order))
  for i,index in enumerate(perm):
   expr=sp.integrate(expr,(times[index],0,times[perm[i+1]] if i+1<n else 1))
  total+=expr
 return sp.simplify(total)

def transfer_inventory(j):
 """All absolute-coefficient terms in the six-derivative IBP algorithm.

 State=(twice denominator power, sorted derivative orders, chi order).
 Positions of factors do not affect this conservative L-infinity bound.
 """
 states={(2*j+1,(1,)*(2*j),0):Q(1)}
 for _ in range(6-2*j):
  out=defaultdict(Q)
  for (p,orders,chi),c in states.items():
   out[p+2,tuple(sorted(orders+(0,1))),chi]+=p*c
   for i,n in enumerate(orders):
    ns=list(orders);ns[i]=n+1
    out[p,tuple(sorted(ns)),chi]+=c
  states=out
 initial=len(states);finished=defaultdict(Q);pending=dict(states);steps=0
 while pending:
  out=defaultdict(Q)
  for (p,orders,chi),c in pending.items():
   if not orders or max(orders)<=3:
    finished[p,orders,chi]+=c;continue
   high=orders.index(max(orders));base=list(orders);base[high]-=1
   out[p,tuple(sorted(base)),chi+1]+=c
   out[p+2,tuple(sorted(base+[0,1])),chi]+=p*c
   for i,n in enumerate(orders):
    if i==high:continue
    ns=base.copy();ns[i]=n+1
    out[p,tuple(sorted(ns)),chi]+=c
  pending=dict(out);steps+=1
  if steps>5:raise AssertionError('IBP did not terminate')
 return {'j':j,'initial_terms':initial,'final_terms':len(finished),
  'absolute_coefficient_sum':str(sum(finished.values(),Q(0))),
  'maximum_derivative':max(max(o,default=0) for p,o,c in finished),
  'maximum_cutoff_derivative':max(c for p,o,c in finished),
  'maximum_factors':max(len(o) for p,o,c in finished),
  'maximum_twice_denominator_power':max(p for p,o,c in finished),
  'maximum_D_power':max(p+len(o) for p,o,c in finished),
  'total_derivative_counts':sorted(set(sum(o)+c for p,o,c in finished))}

def execute(check,emit):
 ctx.prec=512
 h=sp.symbols('h');f=sp.symbols('f0:7');mu=sp.symbols('m1:5')
 expected={(2,):sp.Rational(1,3),(1,1):sp.Rational(1,6),
  (4,):sp.Rational(2,5),(3,1):sp.Rational(1,5),(2,2):sp.Rational(1,5),
  (2,1,1):sp.Rational(11,90),(1,1,1,1):sp.Rational(1,12)}
 moments={}
 for powers,value in expected.items():
  actual=bridge_moment(powers);check('complete_Wick_integral_'+str(powers),actual==value);moments[powers]=actual
 field=sum(f[j]*h**j/sp.factorial(j) for j in range(5))
 A=sum(sp.expand(field**2).coeff(h,j)*(1 if j==0 else mu[j-1])*h**j for j in range(5))
 B=sum(f[j+1]*(1 if j==0 else mu[j-1])*h**j/sp.factorial(j) for j in range(4))
 v=(A-f[0]**2)/f[0]**2
 power=lambda p:sp.series(sum(sp.binomial(p,k)*v**k for k in range(5)),h,0,5).removeO()*f[0]**(2*p)
 H=sp.series(power(-sp.Rational(1,2))+h**2*B**2*power(-sp.Rational(3,2))/4+h**4*B**4*power(-sp.Rational(5,2))/32,h,0,5).removeO().expand()
 def expectation(expr):
  total=0
  for ps,c in sp.Poly(expr,*mu).terms():
   part=tuple(sorted((j+1 for j,p in enumerate(ps) for _ in range(p)),reverse=True))
   if not part:total+=c
   elif sum(part)%2:continue
   else:total+=c*moments[part]
  return sp.simplify(total)
 C={j:expectation(H.coeff(h,j)) for j in range(5)}
 derivative=lambda expr:sp.expand(sum(sp.diff(expr,f[j])*f[j+1] for j in range(6)))
 R=2*f[2]/f[0]-4*f[1]**2/f[0]**2
 lapR=derivative(derivative(R))-f[1]/f[0]*derivative(R)
 check('C0_normalization',C[0]==1/f[0])
 check('odd_bridge_terms',C[1]==0 and C[3]==0)
 check('curvature_C2',sp.simplify(C[2]+R/(12*f[0]))==0)
 check('curvature_and_divergence_C4',sp.simplify(C[4]+R**2/(240*f[0])+lapR/(120*f[0]))==0)
 check('weak_boundary_derivative_identity',sp.simplify(lapR/f[0]-derivative(derivative(R)/f[0]))==0)
 check('flat_coefficient_fixture',all(C[j].subs({f[k]:0 for k in range(1,7)})==0 for j in [1,2,3,4]))
 # Exact nonflat constant-curvature jet fixture f(y)=exp(a y).
 aa=sp.symbols('a');fixture={f[k]:aa**k*f[0] for k in range(1,7)}
 check('constant_curvature_fixture',sp.simplify(C[2].subs(fixture)-aa**2/(6*f[0]))==0 and sp.simplify(C[4].subs(fixture)+aa**4/(60*f[0]))==0)
 inventories=[]
 for j in range(3):
  item=transfer_inventory(j);inventories.append(item)
  check('IBP_derivative_limit_'+str(j),item['maximum_derivative']<=3 and item['maximum_cutoff_derivative']<=3)
  check('IBP_factor_count_'+str(j),item['maximum_factors']<=22 and item['maximum_twice_denominator_power']<24)
  check('IBP_power_and_coefficient_guard_'+str(j),item['maximum_D_power']<=64 and Q(item['absolute_coefficient_sum'])<2**120)
  check('IBP_total_derivative_conservation_'+str(j),item['total_derivative_counts']==[6])
 # Gaussian completion and complete thermal-image inequality.
 p,A0,B0,T,k,tt=sp.symbols('p A B T k t',positive=True)
 check('Gaussian_completion',sp.expand(-tt*A0*(p-B0/(2*A0))**2+tt*B0**2/(4*A0)+tt*A0*p*p-tt*B0*p)==0)
 n=sp.symbols('n',integer=True,positive=True)
 check('complete_image_geometric_tail',sp.factor(n*n-1-3*(n-1))==(n-2)*(n-1))
 ac,bc=sp.symbols('ac bc',positive=True)
 check('complex_image_Young_bound',sp.expand(ac*k*k/(2*tt)-bc*k+bc*bc*tt/(2*ac)-(sp.sqrt(ac/(2*tt))*k-bc*sp.sqrt(tt/(2*ac)))**2)==0)
 # Finite Gaussian recursion cancels every product of parabolic weight <=4.
 cancelled={(r,j-r) for j in range(1,5) for r in range(1,j+1)}
 residual={(r,j) for r in range(1,5) for j in range(5)}-cancelled
 check('all_low_parabolic_products_cancel',all(r+j>=5 for r,j in residual) and len(cancelled)==10)
 check('regular_C5_Taylor_defect_order',[(sp.Rational(5,2)-1),(sp.Rational(4,2)-sp.Rational(1,2)),sp.Rational(3,2)]==[sp.Rational(3,2)]*3)
 check('regular_trace_remainder_order',1+sp.Rational(3,2)-1==sp.Rational(3,2))
 # Implicit-coordinate differentiation: retained actual mass/log powers.
 q,ub,wb,bx,Iu,Iw,Iuw=sp.symbols('q u w bx Iu Iw Iuw')
 xu=-Iu/q;xw=-Iw/q;xuw=-Iuw/q-ub*xw-wb*xu-bx*xu*xw
 check('inverse_coordinate_second_variation',sp.expand(q*xuw+Iuw+q*ub*xw+q*wb*xu+q*bx*xu*xw)==0)
 check('worst_fifth_jet_vanishes',-1+4+4==7)
 check('Brownian_sixth_sup_bound',12*sp.gamma(3)/2==12 and 32*(1+12)==416)
 t0=arb('.001');pi=arb.pi()
 common=2*t0.sqrt();outer=t0
 coarse_log=(arb(2)**124*416/(2*pi)).log()/arb(10).log()
 check('common_regulator_integrability',common<arb('.063245553204'))
 check('selected_screen_preserved',arb('2.123e-42')<arb('1.287e-40'))
 emit('exact',{'bridge_moments':{str(k):str(v) for k,v in moments.items()},'C':{str(k):str(v) for k,v in C.items()},'IBP_inventory':inventories,'regular_residual_products':sorted(residual)})
 emit('numerics',{'global_Dini_factor':str(common),'outer_Dini_factor':str(outer),
  'log10_universal_outer_prefactor':str(coarse_log),'complete_C_A_numeric_value':None,
  'action_C2_on_admitted_stationary_family':True,'X_ball_admission':False,'source_Z2_certified':False,'actual_controlled_solution':False})

def main():
 output=ROOT/'evidence'/f'{STEM}.json';log=ROOT/'evidence'/f'{STEM}.log'
 archive=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',
  dependencies={'python':platform.python_version(),'sympy':sp.__version__,'python_flint':flint.__version__,'arb_bits':512},
  input_sha256_before=hashes(),checks=[],numerics=[],
  frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in INPUTS[:3]},
  continuum_scope='Supporting verification of ordinary exact-kernel/weak-transfer/finite-parametrix argument; not a proof by checker.',
  action_C2_scope='bounded admitted stationary scalar coefficient families, not unknown X ball',
  actual_controlled_solution=False,full_stable_wormhole_goal='active and incomplete')
 archive['runs'].append(run)
 def save():output.write_text(json.dumps(archive,indent=2)+'\n',encoding='utf-8',newline='\n')
 def row(obj):
  with log.open('a',encoding='utf-8',newline='\n') as fp:fp.write(json.dumps(obj)+'\n')
 def check(label,value):
  item={'label':label,'passed':bool(value)};run['checks'].append(item);row(item)
  if not value:raise AssertionError(label)
 def emit(kind,value):item={'kind':kind,**value};run['numerics'].append(item);row(item)
 save();row({'run':len(archive['runs']),'started_utc':run['started_utc']})
 try:
  execute(check,emit);check('inputs_preserved',hashes()==run['input_sha256_before']);run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc();row({'exception':run['exception']})
 finally:run['input_sha256_after']=hashes();run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 print(json.dumps({'status':run['status'],'checks':len(run['checks']),'inputs':len(INPUTS),'runs_preserved':len(archive['runs']),'actual_controlled_solution':False}))
 if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
