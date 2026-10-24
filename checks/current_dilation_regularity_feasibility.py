"""Supporting constants/identities for the ordinary dilation regularity test."""
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
STEM='current-dilation-regularity-feasibility'
OWN=['notes/'+STEM+'.md','checks/current_dilation_regularity_feasibility.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def execute():
 pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 t,z=S.symbols('t z');be,u,ep,ta=S.symbols('beta u eps tau')
 Q=be+u*z*z+(-4*be-2*u-ta)*z**3+(3*be+u+ep)*z**4-2*ep*z**5+ep*z**6
 chi=1-S.integrate(t**8*(1-t)**8,(t,0,t))/S.beta(9,9).expand(func=True)
 check('taper_endpoint_values',chi.subs(t,0)==1 and chi.subs(t,1)==0)
 check('taper_eight_endpoint_jets',all(S.diff(chi,t,j).subs(t,e)==0 for j in range(1,9) for e in [0,1]))
 mapping={be:pars['beta'],u:pars['u'],ep:pars['eps'],ta:pars['tau']}
 def val(e):
  if e in mapping:return mapping[e]
  if e.is_Rational:return arb(int(e.p))/int(e.q)
  if e.is_Add:return sum((val(v) for v in e.args),arb(0))
  if e.is_Mul:
   v=arb(1)
   for a in e.args:v*=val(a)
   return v
  if e.is_Pow:return val(e.base)**int(e.exp)
  raise ValueError(e)
 def bound(poly,N=64):
  n=poly.degree();coeff=[val(poly.nth(i)) for i in range(n+1)];bb=arb(0)
  for cell in range(N):
   left=arb(cell)/N;width=arb(1)/N
   pw=[width**i*sum((coeff[k]*comb(k,i)*left**(k-i) for k in range(i,n+1)),arb(0)) for i in range(n+1)]
   bern=[sum((pw[i]*comb(j,i)/comb(n,i) for i in range(j+1)),arb(0)) for j in range(n+1)]
   bb=max(bb,max(v.abs_upper() for v in bern))
  return bb
 delta=S.expand((Q/be-1).subs(z,(4+t)/100)*chi)
 inner=S.Poly(Q/be-1,z);guards=[];inner_guards=[];cap_guards=[]
 for j in range(6):
  cap=bound(S.Poly(delta,t));inn=sum((val(inner.nth(k)).abs_upper()*arb('.04')**k*k**j for k in range(1,7)),arb(0))
  cap_guards.append(cap);inner_guards.append(inn);guards.append(max(cap,inn))
  delta=S.expand((4+t)*S.diff(delta,t))
 radii={r:sum((guards[j]*arb(r)**j/factorial(j) for j in range(6)),arb(0)) for r in ['.20','.25']}
 check('whole_line_radius_point20',radii['.20']<arb('.104'))
 check('whole_line_radius_point25',radii['.25']<arb('.285'))
 r=arb('.2');d=arb('.104');h1=1/(1-d).sqrt()-1;h2=d/(1-d);eta=(h1*h1+h2*h2).sqrt()
 ratio=eta/r.cos();kernel_ratio=4*ratio;decay=(-r).exp()*r.cos()/2
 inv=r.exp()/(r.cos()*(1-ratio));kernel=1/(1-kernel_ratio)
 check('complex_L2_neumann_ratio',ratio<1)
 check('complex_kernel_neumann_ratio',kernel_ratio<1)
 check('uniform_decay_guard',decay>arb('.40'))
 check('inverse_constant_guard',inv<arb('1.45'))
 check('kernel_constant_guard',kernel<arb('2.13'))
 n,q,c=S.symbols('n q C',positive=True)
 check('exact_divergent_angular_summand',S.simplify(2*(q+2*n)*c/(n*(n+q))-2*c*(1/n+1/(n+q)))==0)
 CA=arb(factorial(5))/r**5;CV=arb(factorial(4))/r**4
 CVdensity=sum((arb(comb(4,j)*factorial(j))/r**j for j in range(5)),arb(0))
 exact_r=S.Rational(1,5)
 check('cauchy_factors',S.factorial(5)/exact_r**5==375000 and S.factorial(4)/exact_r**4==15000 and sum(S.binomial(4,j)*S.factorial(j)/exact_r**j for j in range(5))==18321 and CA.contains(375000) and CV.contains(15000) and CVdensity.contains(18321))
 targetA=arb('3.18e27')/CA;targetV=arb('2.66e30')/CV;targetVdensity=arb('2.66e30')/CVdensity
 check('conditional_contact_target_A',S.Rational('3.18e27')/375000==S.Rational('8.48e21') and targetA.contains(arb('8.48e21')))
 check('conditional_contact_target_V',arb('1.773333333333e26')<targetV)
 K=pars['C']*pars['beta']**2/(15*arb.pi()*pars['B2'])
 L=10*arb.pi()*(pars['Fl']+23*pars['g']/(180*arb.pi()))/pars['g']
 C0=(1+12*arb.pi()*arb(2000)**6).sqrt();b02=arb(2000)**10*(C0-1)/6
 kappa=-arb(6)/5+(b02*pars['g']).log()/2
 check('kappa_bound',kappa.abs_upper()<11)
 check('cotangent_guard',(arb.pi()/2)*(arb.pi()*arb('.005')/2).cot()<201)
 check('digamma_real_guard',(arb('3.005')).digamma().abs_upper()<1)
 check('logarithm_dominates_constant',319*arb(2).log()>216)
 CB=10*K*(3840+L/arb(2).log())
 check('all_frequency_symbol_constant',CB<2)
 ss=S.symbols('s');xx=S.symbols('x');sig=S.Rational(801,200)
 psi8=S.diff(S.exp(-xx*xx),xx,8).subs(xx,0)
 check('bump_eighth_jet',psi8==1680)
 # partial_z^3 D^5 on z^s is z^(s-3) s^6(s-1)(s-2).
 exponent=ss;poly=ss**5
 for j in range(3):poly*=exponent;exponent-=1
 check('exact_dilation_chain_rule',S.expand(poly-ss**6*(ss-1)*(ss-2))==0)
 check('degree_eight_leading_term',S.Poly(poly.subs(ss,sig+xx),xx).degree()==8 and S.Poly(poly.subs(ss,sig+xx),xx).LC()==1)
 rr,cv,beta=S.symbols('r cv beta')
 J=-6*beta*(cv-rr*ss*(ss-3))*(ss-4)*(ss+1);aop=(ss+3)/(3*(ss+1))
 coupled=S.cancel(J*aop)
 check('coupled_scalar_degree_four',S.Poly(coupled,ss).degree()==4)
 check('normal_H3_squared_scaling',-15-1+6+10==0)
 check('scalar_H3_scaling',-S.Rational(15,2)+4+3-S.Rational(1,2)==-1)
 check('pointwise_eighth_growth',-S.Rational(15,2)+8==S.Rational(1,2))
 check('L2_eighth_scaling',-S.Rational(15,2)+8-S.Rational(1,2)==0)
 check('C5_decay',-S.Rational(15,2)+5<0)
 z0=arb('.002');xwidth=arb(1)/100
 check('packet_support_inside_actual_collar',z0*(-xwidth).exp()>arb('.0019') and z0*xwidth.exp()<arb('.0021'))
 return dict(inner_guards=[old.encoded(v) for v in inner_guards],cap_guards=[old.encoded(v) for v in cap_guards],whole_line_guards=[old.encoded(v) for v in guards],relative_deviation={r:old.encoded(v) for r,v in radii.items()},
  selected_radius='.20',outward_deviation_guard='.104',L2_neumann_ratio=old.encoded(ratio),kernel_neumann_ratio=old.encoded(kernel_ratio),decay=old.encoded(decay),L2_inverse_constant=old.encoded(inv),kernel_constant=old.encoded(kernel),
  C_A=old.encoded(CA),C_V=old.encoded(CV),C_V_density=old.encoded(CVdensity),conditional_A_amplitude_target=old.encoded(targetA),conditional_V_amplitude_target=old.encoded(targetV),conditional_V_density_target=old.encoded(targetVdensity),
  actual_normal_symbol_upper_constant=old.encoded(CB),kappa=old.encoded(kappa),bump_eighth_jet=str(psi8),dilation_derivative_polynomial=str(S.expand(poly)),coupled_scalar_polynomial=str(coupled),
  seed_complex_first_order_domain_proved=True,seed_complex_resolvent_proved=True,bare_angular_estimator_diverges=True,uniform_five_dilation_jets_in_W3_infinity_obstructed_on_fast_graph=True,
  actual_operator_obstructed=False,actual_subtracted_action_diverges=False,complex_renormalized_response_certified=False,actual_contact_jets_certified=False,actual_nonlocal_remainder_certified=False,actual_tame_nonlinear_bound_certified=False,
  full_X_ball_admission_proved=False,full_signed_projection_executed=False,actual_controlled_solution=False,nonlinear_preserved_ball=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,
  physical_source_or_norm_changed=False,normal_or_slow_root_changed=False,finite_physical_modes_selected=False,next_tame_calculation_executed=False)

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
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),summary={k:result[k] for k in ['relative_deviation','L2_inverse_constant','kernel_constant','actual_normal_symbol_upper_constant'] if k in result}),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
