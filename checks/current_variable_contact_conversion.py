"""Support the ordinary seed common-regulator conversion, not the open strong remainder."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
import current_curved_paired_heat_bridge as bridge
STEM='current-variable-contact-conversion'
OWN=['notes/'+STEM+'.md','checks/current_variable_contact_conversion.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(stem):return json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]['result']

def insertion_coefficients():
 h=S.symbols('h');f=S.symbols('f0:7',positive=True);v=S.symbols('v0:7');mu=S.symbols('m1:5')
 field=sum(f[j]*h**j/S.factorial(j) for j in range(5))
 A=sum(S.expand(field**2).coeff(h,j)*(1 if j==0 else mu[j-1])*h**j for j in range(5))
 B=sum(f[j+1]*(1 if j==0 else mu[j-1])*h**j/S.factorial(j) for j in range(4))
 w=(A-f[0]**2)/f[0]**2
 power=lambda p:S.series(sum(S.binomial(p,k)*w**k for k in range(5)),h,0,5).removeO()*f[0]**(2*p)
 H=S.series(power(-S.Rational(1,2))+h*h*B**2*power(-S.Rational(3,2))/4+h**4*B**4*power(-S.Rational(5,2))/32,h,0,5).removeO().expand()
 V=sum(v[j]*(1 if j==0 else mu[j-1])*h**j/S.factorial(j) for j in range(5))
 moments={};coeff=[]
 for j in range(5):
  ans=0
  for ps,c in S.Poly(S.expand(H*V).coeff(h,j),*mu).terms():
   part=tuple(sorted((a+1 for a,b in enumerate(ps) for _ in range(b)),reverse=True))
   if not part:ans+=c
   elif sum(part)%2:continue
   else:
    if part not in moments:moments[part]=bridge.bridge_moment(part)
    ans+=c*moments[part]
  coeff.append(S.factor(ans))
 D=lambda e:S.expand(sum(S.diff(e,f[j])*f[j+1]+S.diff(e,v[j])*v[j+1] for j in range(6)))
 R=2*f[2]/f[0]-4*f[1]**2/f[0]**2
 lap=lambda e:D(D(e))-f[1]/f[0]*D(e)
 a2=-R**2/120-lap(R)/60
 target=[2*v[0]/f[0],(-R*v[0]/6+lap(v[0])/3)/f[0],(a2*v[0]-(R*lap(v[0])+D(R)*v[1])/60+lap(lap(v[0]))/30)/f[0]]
 for j in range(3):check('exact_insertion_heat_coefficient_'+str(j),S.simplify(2*coeff[2*j]-target[j])==0)
 check('odd_insertion_terms_vanish',coeff[1]==0 and coeff[3]==0)
 return dict(coefficients=[str(vv) for vv in coeff],bridge_moments={str(k):str(vv) for k,vv in moments.items()})

def execute():
 exact=insertion_coefficients()
 z,s=S.symbols('z s');be,u,ep,ta,q,C=S.symbols('beta u eps tau q C')
 Q=be+u*z**2+(-4*be-2*u-ta)*z**3+(3*be+u+ep)*z**4-2*ep*z**5+ep*z**6
 der=lambda v:S.diff(v,z)+s*v/z
 lap=lambda v:der(Q*der(v))
 R=-S.diff(Q,z,2);DR=S.diff(Q*S.diff(R,z),z);H0=-q-S.Rational(1,3);H1=(q*q/6-S.Rational(1,30))/C
 H=lap(lap(S.Integer(1)))/30-(R*lap(S.Integer(1))+Q*S.diff(R,z)*s/z)/10-H0*lap(S.Integer(1))/(6*C)-(4*H1/C+R**2/30+DR/15)
 A=S.expand(z**4*H);T=S.expand(z*S.diff(A,z)+(s-4)*A)
 normal=S.Poly(T,z).nth(0)
 check('exact_conversion_normal',S.simplify(normal-be**2*(s-4)*s*(s-1)*(s-2)*(s-3)/30)==0)
 # Independently collect the finite common cutoff coefficients in (3).
 e,rr,dr,uv,gg,dd,h0,h1,cc=S.symbols('e R DR uv gg dd H0 H1 C')
 geometric2=(rr**2*uv-4*rr*gg+2*dr*uv+2*dd)/60
 insert2=(-rr**2/120-dr/60)*uv+rr*gg/60+dd/30
 finite=S.expand(2*cc*geometric2-h0*gg/3+4*h1*uv-4*cc*insert2+4*h1*uv)
 expected=(8*h1+cc*rr**2/15+2*cc*dr/15)*uv-(h0/3+cc*rr/5)*gg-cc*dd/15
 check('complete_finite_conversion_bilinear',S.expand(finite-expected)==0)
 # Whole-frequency geometric/mass comparison identity: differentiate the
 # absolutely convergent Laplace integral in q, integrate from zero.
 x,m=S.symbols('x m',positive=True)
 primitive=m*(S.log(1+x/m)-x/m+x*x/(2*m*m))
 check('all_mass_Laplace_identity',S.simplify(S.diff(primitive,x)-(x/m-x/(m+x)))==0)
 check('Laplace_identity_anchor',primitive.subs(x,0)==0)
 alpha=S.symbols('alpha');a=alpha*(1-alpha)
 for j,want in [(1,S.Rational(1,6)),(2,S.Rational(1,30)),(3,S.Rational(1,140))]:
  check('complete_parameter_moment_'+str(j),S.integrate(a**j,(alpha,0,1))==want)
 check('retained_normal_finite_contact',-S.Rational(31,30)-S.Rational(1,6)==-S.Rational(6,5))
 pars,hh,zh,env=old.scalar.prev.york.actual_parameters()
 mapping={be:pars['beta'],u:pars['u'],ep:pars['eps'],ta:pars['tau'],q:arb(pars['q']),C:pars['C']}
 def val(expr):
  if expr in mapping:return mapping[expr]
  if expr.is_Rational:return arb(int(expr.p))/int(expr.q)
  if expr.is_Add:return sum((val(v) for v in expr.args),arb(0))
  if expr.is_Mul:
   out=arb(1)
   for v in expr.args:out*=val(v)
   return out
  if expr.is_Pow:return val(expr.base)**int(expr.exp)
  raise ValueError(expr)
 prior=read('current-local-quantum-contact-graph')
 caps=[old.decode(v) for v in prior['shape_multiplier_constants']]
 items=[];eta=arb(0);t=S.symbols('t');fac=pars['C']/(arb.pi()*pars['B2'])
 for (p,),coeff in S.Poly(T-normal,z).terms():
  poly=S.Poly(S.expand(coeff.subs(s,t+S.Rational(801,200))),t)
  check('weighted_order_'+str(p),p>=2 and poly.degree()<=5)
  bound=sum((val(poly.nth(k)).abs_upper()*caps[k] for k in range(poly.degree()+1)),arb(0))
  term=fac*arb('.01')**p*(p**3+12*p*p+51*p+100)*bound
  eta+=term;items.append(dict(power=p,symbol=str(coeff),multiplier_bound=old.encoded(bound),value=old.encoded(term)))
 check('conversion_graph_allocation',eta<arb('8.36e-5'))
 remain=arb('.18')-eta
 check('remaining_mass_kernel_allocation',remain>arb('.1799164'))
 check('conversion_fraction_of_canonical_allocation',eta/arb('.18')<arb('.00047'))
 poly=S.Poly(S.expand(R*R),z)
 intR2=sum((val(poly.nth(i))*zh**(i+1)/(i+1) for i in range(poly.degree()+1)),arb(0))
 const=-2/pars['B2']/(4*arb.pi())*(8*val(H1)*zh+pars['C']*intR2/15)
 constold=-4/(arb.pi()*pars['B2'])*(zh*val(H1)+pars['C']*intR2/120)
 check('constant_probe_recovered',const.overlaps(constold))
 check('constant_probe_signed_enclosure',const>arb('-1.319287e-28') and const<arb('-1.319286e-28'))
 return dict(exact_insertion=exact,conversion_normal=str(S.factor(normal)),conversion_row=str(T),
  conversion_graph_feedback_upper=old.encoded(eta),remaining_canonical_mass_kernel_target=old.encoded(remain),
  selected_remaining_mass_kernel_target='.1799164',conversion_fraction=old.encoded(eta/arb('.18')),
  constant_probe_normalized=old.encoded(const),terms=items,shape_multiplier_constants=prior['shape_multiplier_constants'],
  seed_first_two_conversion_jets_proved=True,actual_paired_domain_retained=True,common_regulator_removed_for_seed_jets=True,
  complete_local_conversion_remainder_bound=True,retained_geometric_normal_confirmed=True,normal_or_slow_root_changed=False,
  actual_nonlocal_remainder_certified=False,full_strong_regulator_limit_proved=False,uniform_nonlinear_conversion_proved=False,
  actual_data_target_certified=False,full_X_ball_admission_proved=False,full_signed_projection_executed=False,
  actual_controlled_solution=False,nonlinear_preserved_ball=False,exact_compatibility=False,physical_matching=False,
  validity_noise=False,all_sector_stability=False,physical_source_or_norm_changed=False,finite_physical_modes_selected=False)

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
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),summary={k:result[k] for k in ['conversion_graph_feedback_upper','remaining_canonical_mass_kernel_target','constant_probe_normalized'] if k in result}),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
