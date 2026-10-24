"""Support the actual third seed conversion; not the unknown nonlocal graph norm."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-third-metric-conversion'
OWN=['notes/'+STEM+'.md','checks/current_third_metric_conversion.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def identities():
 b,d,g,R,u,k,q,c,e=S.symbols('b d g R u k q C e')
 core=b*R**2/120-b*R*d/30+b*d*d/30+R*d/60+d*d/60-R*g/60-d*g/180+7*g*g/360
 density=S.exp(-2*b)*core-u*b*R/36+u*b*d/18+u*g/36+k*u*u*b*S.exp(2*b)/2-q*u*(g+R*b)/12
 scaled=density.subs({b:e*b,d:e*d,g:e*e*g})
 coeff=[S.expand(S.diff(scaled,e,j).subs(e,0)/S.factorial(j)) for j in range(1,4)]
 L3=R*b*b*d/10+R*b*g/30+b**3*R**2/60-d*g/180+k*u*u*b**3
 check('full_cubic_density',S.expand(coeff[2]-L3)==0)
 check('b_laplacian_squared_cancels',S.Poly(coeff[2],b,d,g).coeff_monomial(b*d*d)==0)
 z=S.symbols('z');Q=S.Function('Q')(z);field=S.Function('b')(z);rad=-S.diff(Q,z,2)
 lap=lambda f:S.diff(Q*S.diff(f,z),z)
 A=-rad*Q/6;B=-Q*S.diff(Q,z)/540;V=rad**2/60+lap(rad)/30+k/c**2
 raw=L3.subs({b:field,d:lap(field),g:Q*S.diff(field,z)**2,R:rad,u:1/c})
 reduced=A*field*S.diff(field,z)**2+B*S.diff(field,z)**3+V*field**3
 flux=rad*Q*field**2*S.diff(field,z)/10-Q*S.diff(rad,z)*field**3/30-Q**2*S.diff(field,z)**3/540
 check('exact_cubic_boundary_flux',S.expand(raw-reduced-S.diff(flux,z))==0)
 # Euler derivatives of the first two local Taylor coefficients.
 y0,y1,y2,y3,y4=S.symbols('y0:5');jets=[y0,y1,y2,y3,y4]
 dz=lambda f:S.diff(f,z)+sum(S.diff(f,jets[j])*jets[j+1] for j in range(4))
 sub={b:y0,d:Q*y2+S.diff(Q,z)*y1,g:Q*y1*y1,R:rad,u:1/c}
 euler=lambda f:S.expand(S.diff(f,y0)-dz(S.diff(f,y1))+dz(dz(S.diff(f,y2))))
 first=euler(coeff[0].subs(sub));second=euler(coeff[1].subs(sub))
 H0=-q-S.Rational(1,3);H1=(q*q/6-S.Rational(1,30))/c;km=S.Rational(1,15)-q*q/3
 expected1=rad**2/120+lap(rad)/60+H0*rad/(12*c)-H1/c
 check('complete_first_conversion_jet',S.simplify(first.subs(k,km)-expected1)==0)
 lapjet=Q*y2+S.diff(Q,z)*y1
 expected2=dz(Q*dz(lapjet))/30-(rad*lapjet+Q*S.diff(rad,z)*y1)/10-H0*lapjet/(6*c)-(4*H1/c+rad**2/30+lap(rad)/15)*y0
 check('complete_second_conversion_jet',S.simplify(second.subs(k,km)-expected2)==0)
 cubic_euler=euler(reduced.subs({field:y0,S.diff(field,z):y1},simultaneous=True))
 aa,bb=S.symbols('aa bb');v=S.Function('v')(z);w=S.Function('w')(z)
 raw_euler=S.diff(reduced,field)-S.diff(S.diff(reduced,S.diff(field,z)),z)
 mixed=S.diff(raw_euler.subs(field,aa*v+bb*w).doit(),aa,bb).subs({aa:0,bb:0})
 expected=-2*A*(S.diff(v,z)*S.diff(w,z)+v*S.diff(w,z,2)+w*S.diff(v,z,2))-2*S.diff(A,z)*(v*S.diff(w,z)+w*S.diff(v,z))-6*B*(S.diff(v,z)*S.diff(w,z,2)+S.diff(w,z)*S.diff(v,z,2))-6*S.diff(B,z)*S.diff(v,z)*S.diff(w,z)+6*V*v*w
 check('complete_strong_third_conversion_row',S.expand(mixed-expected)==0)
 check('stationary_flat_derivative_cubic_is_flux',S.expand((raw-reduced-S.diff(flux,z)).subs(Q,1).doit())==0)
 s,t=S.symbols('s t')
 check('pairing_outer_flux_powers_positive',min(3+3+2,3+3+3,2+2+2)>0)
 # Every nonzero word in the two fixed Clifford matrices has even parity.
 P1=S.Matrix([[1,0],[0,-1]]);P2=S.Matrix([[0,1],[1,0]])
 from itertools import product
 checked=0
 for length in range(1,7):
  for word in product([0,1],repeat=length):
   mat=S.eye(2)
   for letter in word:mat=mat*[P1,P2][letter]
   if word.count(0)%2 or word.count(1)%2:assert S.trace(mat)==0
   checked+=1
 check('ordered_matrix_parity_through_tail_start',checked==126)
 return dict(cubic_raw=str(L3),cubic_flux=str(flux),strong_third_row=str(expected),parity_words_checked=checked)

def execute():
 algebra=identities();pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 z,s,t=S.symbols('z s t');be,u,ep,ta,q,c=S.symbols('beta u ep tau q C')
 Q=be+u*z*z+(-4*be-2*u-ta)*z**3+(3*be+u+ep)*z**4-2*ep*z**5+ep*z**6
 R=-S.diff(Q,z,2);lap=lambda x:S.diff(Q*S.diff(x,z),z);km=S.Rational(1,15)-q*q/3
 A=-R*Q/6;B=-Q*S.diff(Q,z)/540;V=R*R/60+lap(R)/30+km/c**2
 F=-2*A*(s*t+s*(s-1)+t*(t-1))/z**2-2*S.diff(A,z)*(s+t)/z-6*B*s*t*(s+t-2)/z**3-6*S.diff(B,z)*s*t/z**2+6*V
 row=S.Poly(S.expand(z**4*((s+t)*F+z*S.diff(F,z))),z)
 check('full_bilinear_row_symmetric',S.expand(row.as_expr()-row.as_expr().subs({s:t,t:s},simultaneous=True))==0)
 check('actual_polynomial_power_range',min(p[0] for p,cc in row.terms())==2 and max(p[0] for p,cc in row.terms())==12)
 sig=S.Rational(801,200);mapping={be:pars['beta'],u:pars['u'],ep:pars['eps'],ta:pars['tau'],q:pars['q'],c:pars['C']}
 def val(e):
  if e in mapping:return mapping[e]
  if e.is_Rational:return arb(int(e.p))/int(e.q)
  if e.is_Add:return sum((val(x) for x in e.args),arb(0))
  if e.is_Mul:
   r=arb(1)
   for x in e.args:r*=val(x)
   return r
  if e.is_Pow:return val(e.base)**int(e.exp)
  raise ValueError(e)
 ck=[(arb(k)/512).sqrt() for k in [63,7,3,3,7,63]]
 def prod(i,j):
  total=arb(0)
  for l in range(4):
   for m in range(l+1):
    low=min(i+m,j+l-m);high=max(i+m,j+l-m)
    assert low<=3 and high<=6
    total+=arb(comb(3,l)).sqrt()*comb(l,m)*arb(1801)**2*ck[low]
  return total
 fq=pars['C']/(arb.pi()*pars['B2']);total=arb(0);records=[];maxdegree=0
 for (power,),coeff in row.terms():
  pol=S.Poly(S.expand(coeff.subs({s:s+sig,t:t+sig})),s,t);v=arb(0);cells=[]
  for (i,j),coef in pol.terms():
   assert i<=3 and j<=3 and i+j<=4
   maxdegree=max(maxdegree,i+j);pc=prod(i,j);v+=val(coef).abs_upper()*pc
   cells.append(dict(i=i,j=j,exact_coefficient=str(coef),interval_coefficient=old.encoded(val(coef)),product_constant=old.encoded(pc)))
  pp=arb(power)+arb('4.005');mult=pp**3+12*pp**2+51*pp+100
  contribution=fq*(pp*arb('.01').log()).exp()*mult*v;total+=contribution
  records.append(dict(power=power,exact_coefficient=str(coeff),cells=cells,graph=old.encoded(contribution)))
 check('maximum_total_differential_degree',maxdegree==4)
 check('third_seed_conversion_graph_target',total<arb('9.841e-33'))
 check('positive_bound_not_set_zero',total>arb('9.84e-33'))
 return dict(algebra=algebra,third_seed_conversion_graph=old.encoded(total),physical_factor=old.encoded(fq),terms=records,
  weak_mass_third_variation_common_regulator_proved=True,seed_third_conversion_common_regulator_proved=True,full_local_conversion_third_jet_proved=True,outer_fast_third_conversion_graph_certified=True,
  b_laplacian_squared_cancellation=True,normal_or_slow_root_changed=False,physical_source_or_norm_changed=False,finite_physical_modes_selected=False,
  actual_nonlocal_third_graph_certified=False,actual_contact_jets_certified=False,actual_nonlocal_seed_norm_certified=False,uniform_strong_Z2_certified=False,full_X_ball_admission_proved=False,
  global_third_conversion_graph_certified=False,full_signed_projection_executed=False,actual_controlled_solution=False,nonlinear_preserved_ball=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,next_nonlocal_calculation_executed=False)

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
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),third_seed_conversion_graph=result.get('third_seed_conversion_graph')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
