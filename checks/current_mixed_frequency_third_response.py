"""Support the complete stationary mixed-frequency loop in the retained comparison."""
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
STEM='current-mixed-frequency-third-response'
OWN=['notes/'+STEM+'.md','checks/current_mixed_frequency_third_response.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def algebra():
 a,b,l,y,M=S.symbols('a b l y M',real=True);c=-a-b
 # Vertex labels: D2 at 0, D3 at a, D1 at a+b.
 D1=(l+a+b)**2+y*y+M*M;D2=l*l+y*y+M*M;D3=(l+a)**2+y*y+M*M
 check('collinear_triangle_partial_fraction',S.factor(1+D1/(b*c)+D2/(c*a)+D3/(a*b))==0)
 sx=S.Matrix([[0,1],[1,0]]);sy=S.Matrix([[0,-S.I],[S.I,0]]);eye=S.eye(2)
 shifts=[0,a,a+b];num=[M*eye-S.I*(l+r)*sx-S.I*y*sy for r in shifts]
 dots=sum((l+shifts[i])*(l+shifts[j])+y*y for i,j in [(0,1),(0,2),(1,2)])
 check('complete_Dirac_triangle_trace',S.expand(S.trace(num[0]*num[1]*num[2])-2*M*(M*M-dots))==0)
 check('complete_Dirac_bubble_trace',S.expand(S.trace(num[0]*num[1])-2*(M*M-l*(l+a)-y*y))==0)
 alpha,t,pp=S.symbols('alpha t p',positive=True);aa=alpha*(1-alpha)
 # The integrated total derivative supplies the bubble identity.
 remainder=S.diff((2*alpha-1)*S.log(1+aa*pp**2/M**2),alpha)-2*S.log(1+aa*pp**2/M**2)+(2*alpha-1)**2*pp**2/(M*M+aa*pp**2)
 check('bubble_parameter_integration_by_parts',S.simplify(remainder)==0)
 R3=S.log(1+t)-t+t*t/2;R4=R3-t**3/3
 check('R3_global_guard_derivative',S.simplify(S.diff(t**3/(2*(1+t))-R3,t)-t*t/(2*(1+t)**2))==0)
 lam=S.symbols('lambda',positive=True)
 fp=pp*(S.log(1+lam*pp**2)-lam*pp**2+(lam*pp**2)**2/2-(lam*pp**2)**3/3)
 expected=-(lam*pp**2)**4/pp**2*(126+180*lam*pp**2+70*(lam*pp**2)**2)/(1+lam*pp**2)**3
 check('third_derivative_retains_all_terms',S.simplify(S.diff(fp,pp,3)-expected)==0)
 check('third_derivative_positive_majorant',S.expand(126*(1+t)**2-(126+180*t+70*t*t))==72*t+56*t*t)
 mass,AA=S.symbols('mass A',positive=True)
 primitive=S.log(1+AA/mass**2)/AA+1/(mass**2+AA)
 check('entire_angular_count_integral',S.simplify(-S.diff(primitive,mass)-2*(2*mass**2+AA)/(mass*(mass**2+AA)**2))==0)
 check('logarithm_count_guard',S.simplify(S.diff(S.log(1+t)-t/(1+t),t)-t/(1+t)**2)==0)
 for n in range(2,13):
  In=S.factorial(n)**2/S.factorial(2*n+1)
  Inext=S.factorial(n+1)**2/S.factorial(2*n+3)
  check('beta_moment_recurrence_'+str(n),S.cancel(Inext/In-S.Rational(n+1,2*(2*n+3)))==0)
 check('I2_I3',S.integrate(aa**2,(alpha,0,1))==S.Rational(1,30) and S.integrate(aa**3,(alpha,0,1))==S.Rational(1,140))
 C,kappa,z=S.symbols('C kappa z',real=True)
 K=C*pp**4*(S.log(pp)+kappa)/(30*S.pi)
 # Regulator finite part: -a2'' log epsilon/2, with epsilon -> exp(2c)epsilon.
 cc=S.symbols('c',real=True);a2h=-C*S.exp(-2*cc)*pp**4/(30*S.pi)
 scale_contact=-cc*a2h
 check('same_regulator_soft_contact',S.simplify(S.diff(scale_contact,cc).subs(cc,0)-C*pp**4/(30*S.pi))==0)
 check('quartic_restoration_soft_identity',S.simplify(2*C*kappa*pp**4/(30*S.pi)-pp*S.diff(C*kappa*pp**4/(30*S.pi),pp)+C*pp**4/(30*S.pi)-C*(1-2*kappa)*pp**4/(30*S.pi))==0)
 u,v,r=S.symbols('u v r');arg=a*(u+r-1)+b*(v+r-1)
 for degree in [3,5,7]:
  rhs=S.integrate(S.diff(z**degree,z,3).subs(z,arg),(u,0,1),(v,0,1),(r,0,1))/2
  check('divided_difference_polynomial_'+str(degree),S.factor(rhs-(a**degree+b**degree+c**degree)/(a*b*c))==0)
 # Convexity rather than a maximum-only allocation is needed for 1792.
 F=(1+t*t)**4*S.log(2+t)
 target=8*(1+t*t)**2*(1+7*t*t)*S.log(2+t)+16*t*(1+t*t)**3/(2+t)-(1+t*t)**4/(2+t)**2
 check('convex_weight_second_derivative',S.simplify(S.diff(F,t,2)-target)==0)
 check('convex_weight_margin',8*arb(2).log()>1)
 check('frequency_tame_constant',7*256==1792 and arb(1792)/arb(2).sqrt()<1280)
 check('global_symbol_constant',arb(74)/(5*arb.pi())+arb(33)/(10*arb.pi()*arb(2).log())<7)
 nn,qq=S.symbols('n q',positive=True)
 for rr in range(2,8):
  expanded=0
  for j in range(2,rr+1):
   kk=rr-j;co=S.Rational(j-1,rr-1)*(-1)**kk*S.binomial(rr+kk-2,kk)/qq**(rr+kk-1)
   expanded+=co*(nn**(-j)+(-1)**(j+1)*(nn+qq)**(-j))
  check('entire_moment_partial_fraction_'+str(rr),S.factor(expanded-(2*nn+qq)/(nn**rr*(nn+qq)**rr))==0)
 return dict(triangle_trace_retained=True,both_triangle_orientations_and_three_pair_contacts_retained=True,finite_regulator_scale_contact_retained=True,all_frequency_divided_difference_identity=True)

def execute():
 alg=algebra()
 ctx.prec=512;pars,h,zh,env=old.scalar.prev.york.actual_parameters();C=pars['C'];q=pars['q'];B2=pars['B2'];beta=pars['beta'];MM2=(q+1)/C
 b02=arb(2000)**10*((1+12*arb.pi()*arb(2000)**6).sqrt()-1)/6
 kappa=(arb(q+1).digamma()-arb.const_euler())/2-(B2/b02).log()/2
 assert 0<kappa<16
 def moment(r):
  ans=arb(0)
  for j in range(2,r+1):
   k=r-j
   aa=arb(j-1)/(r-1)*(-1)**k*comb(r+k-2,k)/arb(q)**(r+k-1)
   ans+=aa*(arb(j).zeta()+(-1)**(j+1)*arb(j).zeta(arb(q+1)))
  return 2*C**r*ans
 mom={r:moment(r) for r in range(2,12)}
 e2,e3=S.symbols('e2 e3');powers=[S.Integer(3),S.Integer(0),-2*e2]
 for k in range(3,28):powers.append(S.expand(-e2*powers[k-2]+e3*powers[k-3]))
 def eval_poly(pp,a,b):
  return sum((arb(int(co.p))/int(co.q)*a**i*b**j for (i,j),co in S.Poly(pp,e2,e3).terms()),arb(0))
 polys={n:S.cancel(powers[2*n]-S.Rational(2*n,2*n+3)*powers[2*n+3]/e3).expand() for n in range(3,13)}
 def probe(a,b,J=12):
  cc=-a-b;ee2=a*b+a*cc+b*cc;ee3=a*b*cc;P=max(abs(a),abs(b),abs(cc));value=arb(0)
  for n in range(3,J+1):
   In=arb(factorial(n)**2)/factorial(2*n+1)
   value+=(-1)**(n+1)/arb(n)*In*mom[n-1]*eval_poly(polys[n],ee2,ee3)/(2*arb.pi())
  r=P**2/(4*MM2);tail=5*mom[2]*P**6/(280*arb.pi())*r**(J-2)*((J+2)**2/(1-r)+2*(J+2)*r/(1-r)**2+r*(1+r)/(1-r)**3)
  contact=C*(1-2*kappa)/(30*arb.pi())*ee2**2
  return dict(raw_subtracted_enclosure=old.encoded(value+arb(0,tail.abs_upper())),fourth_order_block_enclosure=old.encoded(value+contact+arb(0,tail.abs_upper())),moment_degree=J,raw_subtracted_center=old.encoded(value),entire_series_tail=old.encoded(tail),fixed_fourth_contact=old.encoded(contact),fourth_order_block_center=old.encoded(value+contact),original_normalized_center=old.encoded(-2*(value+contact)/B2),original_normalized_tail=old.encoded(2*tail/B2))
 probes={'equal_inputs':probe(arb(1),arb(1)),'zero_frequency_input':probe(arb(1),arb(-1)),'near_canceling_inputs':probe(arb(1),-1+arb(1)/2**20)}
 base=2560*C/B2;ceiling=arb('.159916395999')/arb('1e-11');tr=ceiling/base
 ans=dict(kappa4=old.encoded(kappa),first_mass_squared=old.encoded(MM2),all_angular_second_inverse_moment=old.encoded(mom[2]),normalized_Sobolev_prefactor=old.encoded(base),best_case_weighted_transfer_ceiling=old.encoded(tr),equal_split_weighted_transfer_ceiling=old.encoded(tr/2),probes=probes)
 check('first_mass_strictly_above_18',MM2>324)
 check('actual_fixed_fourth_coefficient',kappa>0 and kappa<16)
 check('positive_complete_angular_moments',all(val>0 for val in mom.values()))
 check('moment_second_exact_form',mom[2].overlaps(2*C*C/arb(q)*(arb(2).zeta()-arb(2).zeta(arb(q+1)))))
 check('diagnostic_series_tail_equal',old.decode(probes['equal_inputs']['entire_series_tail'])<arb('2.640e-15'))
 check('diagnostic_series_tail_zero',old.decode(probes['zero_frequency_input']['entire_series_tail'])<arb('3.924e-23'))
 check('diagnostic_series_tail_near_zero',old.decode(probes['near_canceling_inputs']['entire_series_tail'])<arb('3.924e-23'))
 check('signed_complete_spectrum_diagnostics',all(old.decode(row['fourth_order_block_enclosure'])<0 for row in probes.values()))
 check('original_normalized_tame_prefactor',base>arb('1.5e-29') and base<arb('1.500000001e-29'))
 check('best_case_conditional_transfer_threshold',tr>arb('1.066109306e39') and tr<arb('1.066109307e39'))
 check('equal_split_conditional_transfer_threshold',tr/2>arb('5.33054653e38') and tr/2<arb('5.33054654e38'))
 return dict(**ans,algebra=alg,all_frequency_comparison_symbol_bound='7 C P^4 log(2+P)',comparison_tame_Sobolev_constant='1280 C',entire_angular_sums_retained=True,complete_fourth_contact_retained=True,lower_Taylor_contacts_retained_separately=True,soft_and_opposing_frequency_denominators_canceled=True,flat_comparison_finite_part_identified=True,actual_curved_strong_regulator_limit_proved=False,actual_weighted_graph_transfer_known=False,sharp_input_indicator_has_H8_regularization=False,physical_domain_changed=False,physical_source_changed=False,normalization_changed=False,original_graph_changed=False,full_signed_projection_executed=False,diagnostic_values_are_bolt_projection=False,finite_physical_modes_selected=False,actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,next_actual_graph_transfer_executed=False)

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
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),normalized_prefactor=result.get('normalized_Sobolev_prefactor'),transfer_threshold=result.get('best_case_weighted_transfer_ceiling')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
