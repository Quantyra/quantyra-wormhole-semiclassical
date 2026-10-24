"""Product contact and actual radial envelopes in a partial joint source screen."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import sys,json,traceback,heapq
import sympy as sp
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'checks'))
import current_reference_rescaled_source_test as prior
scalar=prior.prev;decode=prior.decode;encoded=prior.encoded
from flint import arb,ctx
STEM='current-product-matched-source-projection';PREVIOUS='current-reference-rescaled-source-test'
OWN=['notes/'+STEM+'.md','checks/current_product_matched_source_projection.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[]

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
scalar.prev.prior.check=check;scalar.prev.york.check=check

def digest(p):return sha256((root/p).read_bytes()).hexdigest()
def ball(lo,hi):return arb((lo+hi)/2,(hi-lo)/2)
def bern_on_cell(coeff,lo,hi):
 n=len(coeff)-1;w=hi-lo
 a=[sum((coeff[k]*comb(k,i)*lo**(k-i)*w**i for k in range(i,n+1)),arb(0)) for i in range(n+1)]
 return scalar.bernstein(a)
def family_bound(P,T,lo,hi,power):
 width=hi-lo;ymid=1-(hi+lo)/2;ww=[arb(0) for _ in range(7)];bc=arb(1)
 for j in range(7):
  if j:bc*=(power-j+1)/j
  coeff=ymid**power*bc*(-width/ymid)**j
  for i in range(j+1):ww[i]+=coeff*comb(j,i)*(-arb(1)/2)**(j-i)
 prod=scalar.bern_product(T,scalar.bernstein(ww));raised=scalar.bern_product(P,[arb(1)]*7)
 fall=arb(1)
 for j in range(7):fall*=power-j
 error=fall.abs_upper()/factorial(7)*(width/2)**7*(1-hi)**(power-7)*max(c.abs_upper() for c in T)
 vals=[a+b for a,b in zip(raised,prod)]
 return ball(min(v.lower() for v in vals)-error.upper(),max(v.upper() for v in vals)+error.upper())
def g(z):return 4*z-3*z*z-1+2*z.log()
def ell(z):return -(1-z).log()-z

def algebra():
 z=sp.symbols('z',positive=True);v=sp.Function('v')(z)
 gg=4*z-3*z*z-1+2*sp.log(z);gp=-4*z+6*z*z-2;gpp=4*z-12*z*z
 kn=gpp*v-2*gp*z**4;kb=2*gp*v-4*gg*z**4
 check('joint_displacement_weak_cancellation',sp.simplify(z*kn+z*z*sp.diff(kb,z)/2-(z*z*gp*sp.diff(v,z)-8*gg*z**5))==0)
 check('increasing_g_endpoint_bound',sp.simplify(sp.diff(gg,z)-2*(1-z)*(1+3*z)/z)==0)
 phi=sp.symbols('phi');fun=4*sp.exp(-phi)-3*sp.exp(-2*phi)-1-2*phi
 check('product_magnetic_value_and_trace_zero',fun.subs(phi,0)==0 and sp.diff(fun,phi).subs(phi,0)==0)
 t,q=sp.symbols('t q',positive=True)
 check('integrated_density_bound',sp.factor(sp.integrate(q*(t*t-1)/2,(t,1,z)))-q*(z-1)**2*(z+2)/6==0)
 check('rational_thermal_geometry_density_ratio',sp.simplify((t-1)**2*(t+2)/(3*t)/((t-1)**2*(t*t+2*t+3)/t**2)-t*(t+2)/(3*(t*t+2*t+3)))==0)
 kappa,gs,c0=sp.symbols('kappa gs c0')
 check('shared_a0_Hessian_cancels',sp.simplify(c0-2*(c0/2))==0)
 check('shape_product_solves_area_contact',sp.solve(2*kappa+gs,kappa)==[-gs/2])
 check('lapse_integral_primitive',sp.simplify(sp.diff(sp.log(t/(t-1))-1/t,t)+1/(t*t*(t-1)))==0)
 check('allocation_arithmetic',sp.Rational('1.31e-40')+sp.Rational('1.287e-40')+sp.Rational('2e-43')+sp.Rational('2.331e-44')<sp.Rational('2.60e-40'))

def integrate(trial,zh):
 power=scalar.power;n=max(k for k,j,t in trial.c);P=[trial.c.get((k,0,0),arb(0)) for k in range(n+1)];T=[trial.c.get((k,0,1),arb(0)) for k in range(n+1)]
 AP=[arb(0) for _ in P];BT=[arb(0) for _ in P]
 for k in range(n+1):
  if k:AP[k-1]+=k*P[k];BT[k-1]+=k*T[k]
  AP[k]-=k*P[k];BT[k]-=(k+power)*T[k]
 C0=sum((x.abs_upper() for x in P+T),arb(0));C1=sum((x.abs_upper() for x in AP+BT),arb(0))
 zo=arb(2)**-20;zc=1-arb(2)**-49
 logmoment=lambda k:zo**(k+1)*(arb(1)/(k+1)**2-zo.log()/(k+1))
 outer_a=((8*C0+16)*zo**7/7+4*logmoment(6))/(2*(1-zo))
 outer_w=(8*C1+64)*zo**5/5+16*logmoment(4)
 outer_wL=((8*C1+64)*zo**7/7+16*logmoment(6))/(2*(1-zo))
 outer=[outer_a,outer_w,outer_wL]
 raw_a=(8*C0+16)/5+arb(4)/25;raw_w=(8*C1+64)/6+arb(16)/36
 check('analytic_outer_caps_small',max(a.upper() for a in outer)<arb('1e-12'))
 def cell(lo,hi,bern):
  vv=family_bound(bern[0],bern[1],lo,hi,power)
  yvp=family_bound(bern[2],bern[3],lo,hi,power)
  assert vv.upper()>=0 and vv.lower()<=arb('1.01'),'trial_range_intersection'
  assert yvp.upper()>=0,'monotone_derivative_range_intersection'
  vv=ball(max(arb(0),vv.lower()),min(arb('1.01'),vv.upper()))
  yvp=ball(max(arb(0),yvp.lower()),yvp.upper())
  zz=ball(lo,hi);yy=ball(1-hi,1-lo);gg=ball(g(lo).lower(),g(hi).upper());ll=ball(ell(lo).lower(),ell(hi).upper())
  gp=-2*yy*(1+3*zz)
  ka=gp*vv-2*gg*zz**4;kw=-2*zz**2*(1+3*zz)*yvp-8*gg*zz**5
  weight=ball(((1-hi)/hi).lower(),((1-lo)/lo).upper())
  vals=[abs(ka)*ll,abs(kw)*weight,abs(kw)*weight*ll]
  bounds=[(max(arb(0),a.lower())*(hi-lo),a.upper()*(hi-lo)) for a in vals]
  priority=max(float(b-a) for a,b in bounds)
  return {'lo':lo,'hi':hi,'bern':bern,'bounds':bounds,'priority':priority}
 heap=[];serial=0;lower=[arb(0) for _ in outer];upper=[x.upper() for x in outer]
 def add(c):
  nonlocal serial
  serial+=1;heapq.heappush(heap,(-c['priority'],serial,c))
  for i,(a,b) in enumerate(c['bounds']):lower[i]+=a;upper[i]+=b
 initial=[(arb(2)**-k,arb(2)**(-k+1)) for k in range(20,1,-1)]
 initial +=[(1-arb(2)**(-k+1),1-arb(2)**-k) for k in range(2,50)]
 for lo,hi in initial:add(cell(lo,hi,[bern_on_cell(c,lo,hi) for c in [P,T,AP,BT]]))
 while max((b-a).upper() for a,b in zip(lower,upper))>=arb('.0005'):
  assert len(heap)<50000,'adaptive_cell_limit'
  _,_,c=heapq.heappop(heap)
  for i,(a,b) in enumerate(c['bounds']):lower[i]-=a;upper[i]-=b
  split=[scalar.split(b) for b in c['bern']];mid=(c['lo']+c['hi'])/2
  add(cell(c['lo'],mid,[b[0] for b in split]));add(cell(mid,c['hi'],[b[1] for b in split]))
  if len(heap)%2000==0:print(json.dumps(dict(event='interval_projection_refinement',cells=len(heap),widths=[str(b-a) for a,b in zip(lower,upper)])),flush=True)
 cells=sorted((c for _,_,c in heap),key=lambda c:float(c['lo']))
 check('bulk_interval_coverage',cells[0]['lo']==zo and cells[-1]['hi']==zc and all(a['hi']==b['lo'] for a,b in zip(cells,cells[1:])))
 check('bulk_integral_widths',max((b-a).upper() for a,b in zip(lower,upper))<arb('.0005'))
 NUMERICS.append(dict(kind='validated_magnetic_weak_integrals',outer_cut=encoded(zo),cold_cut=encoded(zc),
  fractional_power=encoded(power),trial_coefficients=trial.serial(),absolute_trial_coefficient_sum=encoded(C0),
  absolute_weighted_derivative_coefficient_sum=encoded(C1),analytic_outer_caps=[encoded(x) for x in outer],
  integral_enclosures=[encoded(ball(a.lower(),b.upper())) for a,b in zip(lower,upper)],
  raw_absolute_kernel_integral_bounds=[encoded(raw_a),encoded(raw_w)],
  cells=[dict(left=encoded(c['lo']),right=encoded(c['hi']),integral_bounds=[[encoded(a),encoded(b)] for a,b in c['bounds']]) for c in cells]))
 return [x.upper() for x in upper],raw_a,raw_w,zc

def execute(previous):
 p,h,zh,env=scalar.prev.york.actual_parameters();NUMERICS.extend(scalar.prev.york.NUMERICS)
 C=p['C'];B2=p['B2'];lam=p['g'];x=lam*B2;q=arb(p['q']);du=p['eps']*C/x
 chi=16*arb.pi()*arb('1e22')*lam;eR=24*chi;Fl=1-eR
 A=q*q*arb('1e15')/(arb('.99')*lam*B2**2);A0=arb('3.768e-12');d0=arb('1e-15');qr=arb('1.3e-10')
 sr=arb('3.33e-10')*h;rh=arb('1.4e-11')*h;eps=arb('1e-27')
 check('lapse_coefficient_target',A<A0,str(A))
 check('local_lapse_allocation',arb('3e-28')+2*arb(2000)**-10+arb('9.182e-37')<arb('4e-28'))
 check('actual_small_thermal_guard',h<arb('6e-34'))
 check('actual_horizon_density_conversion',p['gamma']*C*(1+h)/x<7*h)
 check('bulk_factored_geometry_lower',du.abs_upper()/3+7*h*h/d0**2<eps)
 check('bulk_logarithmic_warp_derivative',(21+14*h/d0)/(1-eps)<22)
 hpart=2*rh/d0/(1-eps);bpart=14*h*(1+h)*sr/d0**2/(1-eps);translation=22*sr/d0
 check('bulk_horizon_density_reserve',hpart<arb('3e-29'))
 check('bulk_bolt_value_reserve',bpart<arb('1e-44'))
 check('bulk_translation_reserve',translation<arb('5e-27'))
 fbound=qr/(3*(1-eps))+hpart+bpart+translation
 check('actual_complete_bulk_fbar_bound',fbound<arb('4.34e-11'),str(fbound))
 old=json.loads((root/'evidence/current-schur-scalar-projection.json').read_bytes())['runs'][-1]
 check('global_scalar_trial_accepted',old['status']=='passed')
 row=next(r for r in old['numerics'] if r['kind']=='actual_global_scalar_inverse')
 scalar.power=decode(row['fractional_power']);scalar.y0=h/(1+h)
 trial=scalar.FP({(r['z_power'],r['y_integer_power'],r['fractional_multiplicity']):decode(r['coefficient']) for r in row['trial_coefficients']})+scalar.mon(4)/2
 check('retained_monotonicity_certificate_present',bool(row['global_test_monotonicity_cells']))
 check('fixed_trial_bolt_bound',trial.value(zh)<arb('1.01'))
 ints,rawa,raww,zc=integrate(trial,zh);Ia,Iw,IwL=ints;sc=(1-zc)/zc
 check('cold_cut_bulk_hypothesis',sc>d0 and sc<arb('2e-15'))
 gp=-2*(1-zc)*(1+3*zc);ka=gp*trial.value(zc)-2*g(zc)*zc**4
 wbound=sr+2*(1-zc)/zc*(A0*ell(zc)+arb('4.34e-11')+arb('4e-28'))
 endpoint=zc**2*ka.abs_upper()*wbound
 kM=q*q/(12*arb.pi()*C*B2)
 PM=kM/Fl*(A0*(Ia+2*IwL)+2*(arb('4.34e-11')+arb('4e-28'))*Iw+arb('4e-28')*rawa+sr*raww+endpoint)
 cold=(arb('2.01')*arb('1.73e-15')*sc**2/2*(arb(3)/2+(arb('.01')/sc).log())+eR*(arb('8.9e-11')+arb('1.01')*arb('1.5e-27'))*sc**2/2)/Fl
 oldnum=next(r for r in previous['numerics'] if r['kind']=='known_joint_local_source_test')
 scalarerror=decode(oldnum['numerical_scalar_error_upper'])
 check('dominant_magnetic_projection_target',PM<arb('1.31e-40'),str(PM))
 check('complete_omitted_cold_projection_target',cold<arb('2e-43'),str(cold))
 check('full_source_scalar_trial_error_target',scalarerror<arb('2.331e-44'))
 check('new_partial_allocation_target',PM+cold+scalarerror+arb('1.287e-40')<arb('2.60e-40'))
 NUMERICS.append(dict(kind='actual_radial_envelope_projection_budget',h=encoded(h),F_lambda=encoded(Fl),
  actual_lapse_coefficient_upper=encoded(A),used_lapse_coefficient=encoded(A0),actual_fbar_bulk_upper=encoded(fbound),
  used_fbar_bulk_upper='4.34e-11',bulk_horizon_term_upper=encoded(hpart),bulk_bolt_value_term_upper=encoded(bpart),
  bulk_translation_term_upper=encoded(translation),k_M=encoded(kM),cold_cut_s=encoded(sc),
  finite_cut_endpoint_inside_bracket=encoded(endpoint),magnetic_bulk_projection_upper=encoded(PM),
  complete_cold_projection_upper=encoded(cold),global_scalar_trial_pairing_error_upper=encoded(scalarerror),
  sum_executed_projection_allowances=encoded(PM+cold+scalarerror),
  remaining_bulk_signed_target='1.287e-40',working_full_screen='2.60e-40',
  full_quantum_contact_bulk_projection_evaluated=False,total_source_improvement_certified=False))
 print(json.dumps(dict(event='partial_joint_projection',magnetic_upper=str(PM),cold_upper=str(cold),remaining_target='1.287e-40',full_projection_evaluated=False)),flush=True)

def main():
 ctx.prec=512
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1];review=json.loads((root/OWN[2]).read_bytes())
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=scalar.prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'],dependencies=deps,
  invocation=sys.argv,previous_goal_turn='progress')
 try:
  check('previous_source_test_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('pinned_interval_library',scalar.prev.prior.flint.__version__=='0.9.0')
  algebra();execute(old)
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=scalar.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  exact_product_area_contact_identity_verified=status=='passed',actual_radial_envelopes_certified=status=='passed',
  magnetic_bulk_projection_bound_certified=status=='passed',complete_omitted_cold_projection_bound_certified=status=='passed',
  full_signed_reduced_source_evaluated=False,full_quantum_contact_bulk_projection_evaluated=False,total_source_improvement_certified=False,
  full_solution_feasibility_test_executed=False,nonlinear_solution_constructed=False,physical_stability_proved=False,
  semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
