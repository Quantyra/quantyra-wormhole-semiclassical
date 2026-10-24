"""Entire local collar path bound; no global heat or signed-source certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import factorial
import sys,json,traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
from flint import arb,ctx
import current_heat_remainder_transfer_audit as previous
bh=previous.previous;gy=bh.gy
root=Path(__file__).resolve().parents[1]
STEM='current-local-boundary-heat-remainder';PREVIOUS=previous.STEM
OWN=['notes/'+STEM+'.md','checks/current_local_boundary_heat_remainder.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];encoded=bh.encoded
def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
gy.check=check;bh.prior.scalar.prev.prior.check=check
def digest(p):return sha256((root/p).read_bytes()).hexdigest()
def aa(x):return arb(abs(x).upper())
def ar(x):
 x=sp.Rational(x);return arb(int(x.p))/int(x.q)
def mul(a,b,N):
 out=[arb(0) for _ in range(N+1)]
 for i,x in enumerate(a[:N+1]):
  if x==0:continue
  for j,y in enumerate(b[:N+1-i]):
   if y!=0:out[i+j]+=x*y
 return out
def compose(a,z,N):
 out=[arb(0) for _ in range(N+1)];power=[arb(1)]+[arb(0)]*N
 for x in a:
  out=[s+x*p for s,p in zip(out,power)];power=mul(power,z,N)
 return out
def tail(M,k,p,N=64):
 x=sp.symbols('x');s=x**(p-k)*sp.diff(x**(N+1)/(1-x),x,p)
 return arb(str(M))/(arb(3)/4)**k*ar(s.subs(x,sp.Rational(2,3)))
def norm(coeff,k,p,M,odd=False,N=64):
 out=arb(0)
 for n in range(k,N+1):
  if odd and n%2==0:continue
  out+=aa(coeff[n])*(factorial(n)//factorial(n-p))*(arb(1)/2)**(n-k)
 return out+tail(M,k,p,N)

def collar():
 p,h,zh,env=gy.actual_parameters();b=p['beta'].sqrt()
 Q=[p['beta'],arb(0),p['u'],-4*p['beta']-2*p['u']-p['tau'],3*p['beta']+p['u']+p['eps'],-2*p['eps'],p['eps']]
 qp=[i*Q[i]/2 for i in range(1,7)];qpp=[i*qp[i] for i in range(1,6)]
 disk=arb('.4');qr=sum((aa(c)*disk**i for i,c in enumerate(qp)),arb(0));ql=sum((aa(c)*disk**i for i,c in enumerate(qpp)),arb(0))
 check('actual_complex_ODE_coefficient_guards',qr<arb('.225') and ql<arb('1.282') and arb('.408')<b<arb('.409'))
 R=arb('.75');check('analytic_Picard_ball_preserved',arb('.409')*R+arb('.225')*R*R/2<disk)
 check('analytic_Picard_contraction',arb('1.282')*R*R/2<1)
 check('complex_reciprocal_has_no_pole',arb('.408')-R*arb('.225')>arb('.239'))
 check('complex_full_function_guards',arb('.409')/arb('.239')<arb('1.72') and arb('.409')*arb('.225')/arb('.239')**2<arb('1.62'))
 check('complex_metric_variation_guard',arb('.4')**4*(arb('1.72')+arb('1.62')*R)<arb('.08'))
 N=64;z=[arb(0)]*(N+3);z[1]=b
 for n in range(N+1):z[n+2]=compose(qp,z,n)[n]/((n+1)*(n+2))
 zp=[(i+1)*z[i+1] for i in range(N+1)];f=[arb(0)]*(N+1);f[0]=arb(1)
 for n in range(1,N+1):f[n]=-sum((zp[k]*f[n-k]/b for k in range(1,n+1)),arb(0))
 z2=mul(z,z,N);v=[-x for x in mul(z2,z2,N)];V0=[arb(0)]+[v[i-1]/i for i in range(1,N+1)]
 fp=[(i+1)*f[i+1] for i in range(N)]+[arb(0)]
 df=[-a-bb for a,bb in zip(mul(v,f,N),mul(fp,V0,N))]
 check('exact_vanishing_low_odd_jets',f[1]==0 and all(df[i]==0 for i in [0,1,2,3,5]))
 guards=[]
 def guard(name,value,target):
  ok=bool(value<arb(str(target)))
  ROWS.append(dict(name=name,passed=ok,detail=str(value)))
  guards.append(dict(name=name,bound=encoded(value),guard=str(target),passed=ok))
 guard('collar_f_minus_one',sum((aa(f[n])*(arb('.5')**n) for n in range(1,N+1)),arb(0))+tail('1.72',0,0),'.05')
 guard('collar_v_fourth_order',norm(z,1,0,'.4')**4,'.04')
 guard('cap_v_endpoint_first_derivative',norm(v,4,1,'.0256')*arb('.5')**3,'.03')
 guard('cap_v_endpoint_second_derivative',norm(v,4,2,'.0256')*arb('.5')**2,'.3')
 guard('collar_f_second_jet',2*aa(f[2]),'3e-11')
 guard('collar_f_prime_after_linear',norm(f,3,1,'1.72'),'.6')
 guard('collar_f_second_after_constant',norm(f,3,2,'1.72'),'1.2')
 guard('collar_f_third',norm(f,3,3,'1.72'),'2')
 for k,target in enumerate(['.03','.13','.4']):guard('collar_metric_variation_derivative_'+str(k),norm(df,4,k,'.08'),target)
 for k,target in enumerate(['.15','.45','.9','1']):guard('collar_f_odd_derivative_'+str(k),norm(f,3,k,'1.72',odd=True),target)
 for k,target in enumerate(['.003','.02','.11']):guard('collar_metric_odd_derivative_'+str(k),norm(df,7,k,'.08',odd=True),target)
 NUMERICS.append(dict(kind='collar_guard_attempt',guards=guards))
 check('all_full_collar_guards',all(g['passed'] for g in guards))
 r=arb('.5');eps=arb('3e-11')
 check('path_G_guard',eps+arb('.6')*r<arb('.31'))
 check('path_time_G_guard',eps+arb('1.2')*r<arb('.61'))
 check('shifted_metric_vertex_guards',
  arb('.03')+arb('.008')*(arb('.6')*r**3+eps*r**2)<arb('.031') and
  arb('.13')+arb('.008')*(eps*r**2+arb('1.2')*r**3)<arb('.132') and
  arb('.4')+arb('.008')*arb('2')*r**3<arb('.403'))
 NUMERICS.append(dict(kind='full_analytic_collar_certificate',degree=N,analytic_radius='3/4',physical_collar_radius='1/2',
  actual_beta=encoded(p['beta']),actual_gamma=encoded(p['gamma']),Q_coefficients=[encoded(x) for x in Q],
  z_coefficients=[encoded(x) for x in z],f_coefficients=[encoded(x) for x in f],metric_coefficients=[encoded(x) for x in df],
  analytic_sup_guards=dict(z='.4',f='1.72',metric_variation='.08'),full_tail_guards=guards))
 return p

def cap_guards():
 w0=(arb('3e-11')*arb('.5')+arb('.6')*arb('.5')**2)/arb('.95')
 w1=(arb('3e-11')+arb('1.2')*arb('.5'))/arb('.95')+w0*w0
 w2=arb(2)/arb('.95')+3*w0*w1+w0**3
 check('cap_endpoint_log_jets',w0<arb('.16') and w1<arb('.66') and w2<arb('2.43'))
 check('cap_log_derivative_guards',arb('.16')+arb('.66')/2+arb('2.43')/8<arb('.8') and
  arb('3.75')*arb('.8')+arb('.66')+arb('2.43')/2<5 and
  24*arb('.8')+arb('7.5')*(arb('.66')+arb('2.43')/2)+arb('2.43')<36)
 check('cap_relative_warping_guards',arb('.8')<2 and 5+arb('.8')**2<10 and 36+3*arb('.8')*5+arb('.8')**3<100)
 check('cap_value_guards',arb('.95')*(-arb('.4')).exp()>arb('.6') and arb('1.06')*arb('.4').exp()<arb('1.6'))
 check('cap_metric_test_guards',arb('.0025')+arb('.03')/2+arb('.3')/8<arb('.055')+arb('1e-20') and
  arb('3.75')*arb('.055')+arb('.03')+arb('.3')/2<10 and
  24*arb('.055')+arb('7.5')*(arb('.03')+arb('.3')/2)+arb('.3')<100)
 l,t=sp.symbols('l t');A=2+2*l+t*(16+14*l);B=l+t*(1+7*l);C=10*l+11*l*l+t*(16+206*l+145*l*l)
 gap=sp.Poly(sp.expand(33*(1+l)**2-(A+B*A+C).subs(t,sp.Rational(1,1000))),l)
 check('complete_cap_path_derivative_bound',all(c>0 for c in gap.all_coeffs()))
 check('excursion_time_bound_increases',arb(7)/8/(16*arb('.001'))>4)
 moment=arb(0)
 for k in range(9):
  I=lambda j:arb(2).sqrt()**(j-1)*(arb(j+1)/2).gamma()
  moment+=int(sp.binomial(8,k))*((k+1)*I(k)+I(k+1)/15)
 moment*=2*(arb(1)/450).exp()
 check('complete_thermal_path_moment',moment<arb('1e10'))
 NUMERICS.append(dict(kind='complete_auxiliary_cap_guards',log_endpoint_bounds=[encoded(x) for x in [w0,w1,w2]],
  global_f_range=['.6','1.6'],relative_f_derivative_bounds=['2','10','100'],global_test_bounds=['1','10','100'],
  complete_path_derivative_coefficient='33',thermal_path_moment_bound=encoded(moment),
  actual_disk_comparison_proved=False))

def H_guards():
 F,G,t=sp.symbols('F G t',positive=True);H=sp.exp(t*G*G/(4*F))/sp.sqrt(F)
 specs=[(1,0,0,1,'.7'),(0,1,0,t*G,'.7'),(2,0,0,1,'1.2'),(1,1,0,t*G,'1.2'),(0,2,0,t,'.7'),
  (3,0,0,1,'3.5'),(2,1,0,t*G,'3.5'),(1,2,0,t,'1.2'),(0,3,0,t*t*G,'1.5'),
  (1,0,1,G*G,'1.2'),(0,1,1,G,'.7'),(2,0,1,G*G,'3.5'),(1,1,1,G,'1.2'),(0,2,1,1,'.7')]
 fmin=arb('.9025');fmax=arb('1.1236');gmax=arb('.155');tmax=arb('.001')
 hmax=(tmax*gmax*gmax/(4*fmin)).exp()/fmin.sqrt();details=[]
 for nf,ng,nt,factor,guard in specs:
  expr=sp.expand(sp.simplify(sp.diff(H,F,nf,G,ng,t,nt)/H/factor))
  upper=arb(0)
  for term in sp.Add.make_args(expr):
   powers=term.as_powers_dict();a=int(powers.get(F,0));b=int(powers.get(G,0));c=int(powers.get(t,0))
   assert b>=0 and c>=0
   coeff=sp.simplify(term/(F**a*G**b*t**c));assert coeff.is_Rational
   upper+=ar(abs(coeff))*(fmin if a<0 else fmax)**a*gmax**b*tmax**c
  upper*=hmax;name='H_derivative_'+str(nf)+str(ng)+str(nt)
  check(name,upper<arb(guard),str(upper));details.append(dict(derivatives=[nf,ng,nt],removed_factor=str(factor),bound=encoded(upper),guard=guard))
 NUMERICS.append(dict(kind='differentiated_path_function_bounds',rows=details))

def polynomial_bound(p):
 V,t=sp.symbols('V t',nonnegative=True);q=lambda x:sp.Rational(str(x))
 M=q('1.06');G=q('.31')*V
 dx=[2*M*q('.30')*V**3,q('.90')*V**2]
 tx=[2*M*q('.31')*V**2,q('.61')*V]
 dtx=[2*(M*q('.90')*V**3+q('.30')*q('.31')*V**5),q('1.8')*V**2]
 sx=[2*M*q('.031')*V**4,q('.132')*V**3]
 dsx=[2*(M*q('.0132')+q('.30')*q('.031'))*V**7,q('.0544')*V**6]
 tsx=[2*(M*q('.132')*V**4+q('.31')*q('.031')*V**6),q('.403')*V**3]
 dtsx=[2*(M*q('.0544')+q('.30')*q('.132')+q('.90')*q('.031'))*V**7+2*q('.31')*q('.0132')*V**9,q('.236')*V**6]
 tab={(1,0,0):q('.7'),(0,1,0):q('.7')*t*G,(2,0,0):q('1.2'),(1,1,0):q('1.2')*t*G,(0,2,0):q('.7')*t,
  (3,0,0):q('3.5'),(2,1,0):q('3.5')*t*G,(1,2,0):q('1.2')*t,(0,3,0):q('1.5')*t*t*G,
  (1,0,1):q('1.2')*G*G,(0,1,1):q('.7')*G,(2,0,1):q('3.5')*G*G,(1,1,1):q('1.2')*G,(0,2,1):q('.7')}
 def hd(I,j=None,k=None,nt=0):
  nf,ng=I
  for index in [j,k]:
   if index==0:nf+=1
   elif index==1:ng+=1
  return tab[nf,ng,nt]
 def delta(I):return sum(hd(I,j)*dx[j] for j in range(2))
 def timepoint(I):return sum(hd(I,j)*tx[j] for j in range(2))+t*hd(I,nt=1)
 def timedelta(I):
  return sum(hd(I,j)*dtx[j]+t*hd(I,j,nt=1)*dx[j] for j in range(2))+sum(hd(I,j,k)*dx[k]*tx[j] for j in range(2) for k in range(2))
 first=[(1,0),(0,1)];vv=q('.04')*V**4
 val=vv*delta((0,0))+sum(hd(I)*dsx[j]+delta(I)*sx[j] for j,I in enumerate(first))
 dtval=vv*timedelta((0,0))+sum(timepoint(I)*dsx[j]+hd(I)*dtsx[j]+timedelta(I)*sx[j]+delta(I)*tsx[j] for j,I in enumerate(first))
 poly=sp.Poly(sp.expand(val+dtval),V,t);terms=[];bound=arb(0);t0=arb('.001')
 for (power,ntime),coef in poly.terms():
  weight=power+2*ntime;check('nonnegative_path_monomial_'+str(power)+'_'+str(ntime),coef>0 and weight>=7)
  moment=arb(2**(power+1)-1)/(power+1)*2*(arb(power+3)/2).gamma()
  contribution=p['beta'].sqrt()/arb.pi()*ar(coef)*moment*t0.sqrt()**(weight-7)
  bound+=contribution;terms.append(dict(V_power=power,t_power=ntime,coefficient=str(coef),parabolic_weight=weight,moment_factor=encoded(moment),D_allowance=encoded(contribution)))
 check('seventh_bridge_moment_factor',arb(255)/8*2*arb(5).gamma()==1530)
 check('fifth_bridge_moment_factor',arb(63)/6*2*arb(4).gamma()==126)
 T=2*arb.pi()*p['beta'].sqrt()/p['gamma'];check('original_period_guard',T>arb('1e33'))
 thermal=arb('1e26')*11**6*2**101*factorial(100)*arb('1e67')**-94
 check('complete_thermal_image_tail',thermal<arb('1e-100'))
 c=arb(7)/8
 exited=220*p['beta'].sqrt()/arb.pi()*(t0-c/(16*t0)).exp()*(arb(25)/(64*t0**4)+arb(95)/(16*c*t0**3)+8/(c*c*t0*t0))
 check('complete_auxiliary_excursion_allowance',exited<arb('1e-8'),str(exited))
 total=bound+thermal+exited;check('prespecified_local_boundary_budget',total<600,str(total))
 guard=next(n for n in range(1,601) if arb(n)>total)
 remaining=arb('1777.7385')-guard;check('positive_other_remainder_capacity',remaining>0)
 cap=arb('2.5e-43')*(arb('1.2e-8')*guard*arb('.001')**2/2)/arb('1.0666666e-11')
 NUMERICS.append(dict(kind='whole_local_path_bound',polynomial=str(poly.as_expr()),terms=terms,
  original_period=encoded(T),zero_winding_D_bound=encoded(bound),all_nonzero_winding_D_bound=encoded(thermal),
  all_auxiliary_excursions_D_bound=encoded(exited),complete_capped_boundary_parametrix_bounded=True,
  local_D_bound=encoded(total),local_D_rounded_up_integer=guard,other_D_sufficient=encoded(remaining),
  local_B_bound=str(arb('1.2e-8')*guard),normalized_local_cap_allowance=encoded(cap),
  complete_actual_domain_remainder_bounded=False,excursion_and_propagation_bound_known=False,
  full_heat_cap_certified=False,complete_signed_projection_evaluated=False))

def main():
 ctx.prec=512
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 review=json.loads((root/OWN[2]).read_bytes());paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=bh.prior.scalar.prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,dependencies=deps,invocation=sys.argv,
  frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'])
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
 def save():out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 save()
 try:
  check('previous_scoped_execution_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('interval_library_pinned',bh.prior.scalar.prev.prior.flint.__version__=='0.9.0')
  p=collar();cap_guards();H_guards();polynomial_bound(p)
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=bh.prior.scalar.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  continuum_argument_scope='Ordinary local paired path/analytic-tail argument. Checkers verify explicit guards and integrals, not global propagation.',
  complete_actual_domain_remainder_bounded=False,excursion_and_propagation_bound_known=False,
  C2_action_on_rough_base_family_proved=False,complete_signed_projection_evaluated=False,full_heat_cap_certified=False,
  actual_source_failure_proved=False,total_source_improvement_certified=False,nonlinear_solution_constructed=False,
  physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
 save()
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
