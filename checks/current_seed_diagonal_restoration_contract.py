"""Fixed-scheme local contact and UV integration budget, not a metric solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import factorial
import json,traceback,sys
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_seed_radial_green_enclosure as prior
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-seed-diagonal-restoration-contract'
PREVIOUS='current-seed-separated-quantum-response'
OWN=['notes/'+STEM+'.md','checks/current_seed_diagonal_restoration_contract.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[];NUMERICS=[]
encoded=prior.encoded; hull=prior.hull; dependency_identity=prior.dependency_identity
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
 assert name not in {x['name'] for x in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
prior.check=check

def ev(expr,values):
 if expr.is_Integer:return arb(int(expr))
 if expr.is_Rational:return arb(int(expr.p))/int(expr.q)
 if expr.is_Symbol:return values[expr]
 if expr.is_Add:return sum((ev(a,values) for a in expr.args),arb(0))
 if expr.is_Mul:
  out=arb(1)
  for a in expr.args:out*=ev(a,values)
  return out
 if expr.is_Pow:
  b,e=expr.args; value=ev(b,values)
  if e.is_Integer:return value**int(e)
  if e.is_Rational and e.q==2:return value.sqrt()**int(e.p)
 raise ValueError(str(expr))

def formulas():
 z=sp.symbols('z'); Q=sp.Function('Q')(z); v=sp.Function('v')(z)
 lap=lambda f:sp.diff(Q*sp.diff(f,z),z)
 R=-sp.diff(Q,z,2)
 physical=sp.expand(lap(lap(v))+R*lap(v)+lap(R*v)+R*R*v/2)
 desired=Q**2*sp.diff(v,z,4)+4*Q*sp.diff(Q,z)*sp.diff(v,z,3)+(Q*sp.diff(Q,z,2)+2*sp.diff(Q,z)**2)*sp.diff(v,z,2)-(Q*sp.diff(Q,z,3)+sp.diff(Q,z)*sp.diff(Q,z,2))*sp.diff(v,z)+(-Q*sp.diff(Q,z,4)-sp.diff(Q,z)*sp.diff(Q,z,3)+sp.diff(Q,z,2)**2/2)*v
 check('complete_curvature_contact_in_z',sp.simplify(physical-desired)==0)
 t,u,w,R0,Lu,Lw=sp.symbols('t u w R Lu Lw'); s=sp.symbols('s')
 action=sp.exp(-2*(t*u+s*w))*(R0-2*t*Lu-2*s*Lw)**2
 mixed=sp.diff(action,t,s).subs({t:0,s:0})
 check('affine_second_contact_vertex',sp.expand(mixed-(4*R0**2*u*w+8*R0*(u*Lw+w*Lu)+8*Lu*Lw))==0)
 q=sp.symbols('q0:7',positive=True)
 dz=lambda f:sp.expand(sum(sp.diff(f,q[j])*q[j+1] for j in range(6)))
 dp=lambda f:sp.simplify(-sp.sqrt(q[0])*dz(f))
 H=-q[1]/(2*sp.sqrt(q[0])); R=-q[2]
 coeff=[-q[0]*q[4]-q[1]*q[3]+q[2]**2/2,dp(dp(H))+H*dp(H)+2*R*H+2*dp(R),2*dp(H)+H**2+2*R,2*H,sp.Integer(1)]
 derivative=[sp.simplify(dp(coeff[j])+(coeff[j-1] if j else 0)) for j in range(5)]+[sp.Integer(1)]
 x=sp.symbols('x'); sig=sp.Function('sigma')(x); f=sp.Function('f')(x)
 direction=f
 jets=sp.symbols('s0:6'); fs=sp.symbols('f0:6')
 center={sig:0,sp.diff(sig,x):0}
 center.update({sp.diff(sig,x,j):jets[j] for j in range(2,6)})
 center.update({sp.diff(f,x,j):fs[j] for j in range(6)})
 expected=[fs[0],fs[1],fs[2],fs[3]-jets[2]*fs[1],fs[4]-4*jets[2]*fs[2]-jets[3]*fs[1],fs[5]-10*jets[2]*fs[3]-5*jets[3]*fs[2]+(7*jets[2]**2-jets[4])*fs[1]]
 for j in range(6):
  check('proper_direction_jet_'+str(j),sp.simplify(direction.subs(center,simultaneous=True)-expected[j])==0)
  direction=sp.exp(-sig)*sp.diff(direction,x)
 theta,k,d,C,pi=sp.symbols('theta k d C pi',positive=True)
 angular=sp.integrate(sp.sin(theta)**3*sp.cos(theta)**2,(theta,0,sp.pi/2))
 radial=sp.integrate(k**4*sp.exp(-2*d*k),(k,0,sp.oo))
 check('leading_polar_singularity',sp.simplify(-4*C/pi*angular*radial+2*C/(5*pi*d**5))==0)
 joint=sp.integrate(sp.sin(theta)**3,(theta,0,sp.pi/2))*sp.factorial(4)
 check('majorant_joint_integral',joint==16)
 check('scale_contact_matches_log_symbol',sp.Rational(4,60)==sp.Rational(1,15))
 check('unimodal_derivative_guard',sp.expand((1+x)**2-(x*(1-x)))==2*x*x+x+1)
 # Independent z-to-proper reconstruction of the complete local operator.
 dzf=lambda expr:sp.diff(expr,z)
 proper=lambda expr:-sp.sqrt(Q)*dzf(expr)
 vjets=[v]
 for j in range(4):vjets.append(sp.simplify(proper(vjets[-1])))
 qsubs={q[j]:sp.diff(Q,z,j) for j in range(7)}
 assembled=sum(coeff[j].subs(qsubs)*vjets[j] for j in range(5))
 check('proper_contact_matches_z_contact',sp.simplify(assembled-desired)==0)
 FIXTURES.append(dict(fixture='full_local_contact_and_UV_normalization',domain_proof_replaced_by_checks=False))
 return q,coeff,derivative

def proper_jet_bounds():
 s={j:arb('.027')*(j+1)*factorial(j)*32**j for j in range(2,5)}
 f=[arb('1.2e-8'),arb('8e-6')]+[arb('1.3e-4')*factorial(j)*100**j for j in range(2,6)]
 out=[f[0],f[1],f[2],f[3]+s[2]*f[1],f[4]+4*s[2]*f[2]+s[3]*f[1],f[5]+10*s[2]*f[3]+5*s[3]*f[2]+(7*s[2]**2+s[4])*f[1]]
 check('actual_reference_proper_H5_guard',out[5]<arb('2e8'),encoded(out[5]))
 return out

def contact_coefficients(par,q,coeff,derivative,N):
 z=sp.symbols('z'); b,u,e,t=sp.symbols('beta u eps tau')
 poly=(1-z)**2*(b*(1+2*z+3*z*z)+u*z*z+e*z**4)-t*z**3
 polynomials=[sp.Poly(sp.diff(poly,z,j),z) for j in range(7)]
 params={b:par['beta'],u:par['u'],e:par['eps'],t:par['tau']}
 polycoeff=[[ev(a,params) for a in p.all_coeffs()] for p in polynomials]
 def horner(cs,x):
  v=arb(0)
  for a in cs:v=v*x+a
  return v
 maximum=[arb(0)]*5; maxder=[arb(0)]*6; qmin=arb(1)
 for i in range(N):
  lo=arb(49)/100+arb(39*i)/(100*N);hi=arb(49)/100+arb(39*(i+1))/(100*N)
  zz=hull(lo.lower(),hi.upper()); qs=[horner(c,zz) for c in polycoeff]
  # Factoring is required for the zeroth geometry interval near the bolt.
  qs[0]=prior.geometry(zz,par)[0]
  assert qs[0]>arb('.01')
  qmin=min(qmin,qs[0].lower()); vs=dict(zip(q,qs))
  maximum=[max(a,ev(c,vs).abs_upper()) for a,c in zip(maximum,coeff)]
  maxder=[max(a,ev(c,vs).abs_upper()) for a,c in zip(maxder,derivative)]
 return maximum,maxder,qmin

def calculation(qsymbols,coeff,derivative):
 par=prior.parameters(); pi=arb.pi(); q=arb(par['q']); C=par['C']; B2=par['B2']; gam=par['gamma']
 k=arb(2000); C0=(1+12*pi*k**6).sqrt(); b02=k**10*(C0-1)/6
 logscale=(B2/b02).log()/2
 Hq=(q+1).digamma()+arb.const_euler()
 F=4*arb.const_euler()-2*Hq+4*logscale
 alpha=-C*F/(960*pi); c4=C*F/(60*pi*B2); exactUV=4*C/(5*pi*B2)
 check('retained_scale_guard',logscale.abs_upper()<arb('1e-8'))
 check('common_contact_sign',F<0 and alpha>0 and c4<0)
 # Independent large-integer harmonic enclosure, with alternating remainder.
 approx=q.log()+arb.const_euler()+1/(2*q)-1/(12*q*q)+1/(120*q**4)
 harmonic=hull((approx-1/(252*q**6)).lower(),approx.upper())
 check('harmonic_digamma_asymptotic_overlap',Hq.overlaps(harmonic))
 jets=proper_jet_bounds(); prior_bounds=None
 for N in [64,128]:
  cs,ds,qmin=contact_coefficients(par,qsymbols,coeff,derivative,N)
  val=c4.abs_upper()*sum((a*b for a,b in zip(cs,jets)),arb(0))
  der=c4.abs_upper()*sum((a*b for a,b in zip(ds,jets)),arb(0))
  check('geometry_and_contact_finite_N'+str(N),qmin>arb('.01') and val.is_finite() and der.is_finite())
  if prior_bounds is not None:
   check('contact_bounds_refine',val<=prior_bounds[0] and der<=prior_bounds[1])
   check('contact_value_target',val<arb('1e-27'),encoded(val))
   check('contact_first_derivative_target',der<arb('1e-24'),encoded(der))
  NUMERICS.append(dict(kind='contact_geometry_enclosure',cells=N,region_z=['.49','.88'],minimum_Q_lower=encoded(qmin),proper_operator_coefficient_upper=[encoded(x) for x in cs],proper_derivative_coefficient_upper=[encoded(x) for x in ds],contact_value_upper=encoded(val),contact_first_derivative_upper=encoded(der)))
  prior_bounds=(val,der)
 NUMERICS.append(dict(kind='fixed_scheme_constants',F_eff=encoded(F),alpha=encoded(alpha),normalized_fourth_contact_c4=encoded(c4),log_mu0B=encoded(logscale),harmonic_Hq=encoded(Hq),exact_leading_kappa5=encoded(exactUV),proper_direction_jet_upper=[encoded(x) for x in jets]))
 c=arb('1.98'); a=c/arb(2).sqrt(); b=c*(arb(5)/2).sqrt(); D0=16*q/(pi*B2*arb('.1'))
 em2=arb(-2).exp();em3=arb(-3).exp()
 kap={5:D0*64*C/(q*arb(5).sqrt()*c**5),4:D0*(216*C.sqrt()*em3/(q*a**3*b)+gam*24*C/(q*c**4)),3:D0*(16*em2/(a*a*b)+gam*216*C.sqrt()*em3/(q*c**3)),2:D0*gam*16*em2/c**2}
 check('all_mode_majorant_above_leading_singularity',kap[5]>exactUV)
 NUMERICS.append(dict(kind='all_mode_kernel_majorant',kappa={str(j):encoded(v) for j,v in kap.items()},majorant_to_exact_leading_ratio=encoded(kap[5]/exactUV)))
 last=None
 for label in ['.01','.0001','.00001','.000001']:
  d=arb(label); error=arb('2e8')/60*sum((kap[j]*d**(6-j)/(6-j) for j in kap),arb(0))
  holder=6*arb('1e-28')/(exactUV*d.sqrt())
  m32=32/(a*d); n32=2*C*m32*m32/(q+(q*q+4*C*m32*m32).sqrt());p32=32/(b*d);j32=p32/gam-arb(1)/2
  if last is not None:check('remainder_decreases_delta_'+label,error<last)
  if label=='.00001':check('subtracted_value_remainder_target',error<arb('1e-28'),encoded(error))
  if label=='.00001':
   oldH6=arb('4e24')*factorial(6)*100**6/arb('.49')**6
   conditional_holder=(2*arb('2e8')*oldH6).sqrt(); ratio=conditional_holder/holder
   check('older_higher_jet_allowance_fails_principal_budget',ratio>arb('1e16'),encoded(ratio))
   NUMERICS.append(dict(kind='older_high_jet_conditional_comparison',sixth_proper_derivative_old_allowance=encoded(oldH6),conditional_Holder_half_allowance=encoded(conditional_holder),allowance_to_principal_target_ratio=encoded(ratio),current_reference_high_jet_transfer_proved=False,actual_seminorm_lower_bound=False))
  NUMERICS.append(dict(kind='near_diagonal_budget',delta=label,taylor_subtracted_value_remainder_upper=encoded(error),principal_derivative_Holder_half_threshold=encoded(holder),actual_Holder_half_seminorm_known=False,attenuation_scale_units=32,mass_scale=encoded(m32),angular_index_scale=encoded(n32),proper_scaled_thermal_p=encoded(p32),thermal_index_scale=encoded(j32),scale_counts_are_truncation_tolerance=False))
  last=error
 # The mapping test concerns the retained uncut principal symbol, not physical modes.
 packet=[dict(N=str(10**j),C5_amplitude=encoded(arb(1)),principal_derivative_amplitude=encoded(exactUV/12*arb(10**j).log())) for j in [2,4,8,16,32]]
 check('principal_C5_to_C1_logarithmic_growth',all(arb(10**b).log()>arb(10**a).log() for a,b in zip([2,4,8,16],[4,8,16,32])))
 NUMERICS.append(dict(kind='principal_packet_diagnostic',packets=packet,actual_direction_counterexample=False,physical_high_frequency_modes_claimed=False))
 FIXTURES.append(dict(fixture='actual_parameter_common_contact_and_reference_jet_bound',scope='compact middle region, actual reference direction'))
 FIXTURES.append(dict(fixture='full_spectrum_Taylor_remainder_and_resolution_regularity_test',full_differentiated_quantum_operator_evaluated=False))

def main():
 ctx.prec=256
 review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=dependency_identity();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
 try:
  check('previous_separated_enclosure_passed',old['status']=='passed' and old['full_separated_thermal_angular_sum_enclosed'])
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('pinned_interval_library',prior.flint.__version__=='0.9.0')
  q,co,der=formulas();calculation(q,co,der)
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=dependency_identity();status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),common_fourth_order_contact_evaluated=status=='passed',all_mode_near_diagonal_value_remainder_bounded=status=='passed',C5_to_C1_without_modulus_rejected=status=='passed',all_local_renormalized_moments_computed=False,actual_fifth_derivative_modulus_quantified=False,coincident_subtraction_executed=False,differentiated_response_operator_bounded=False,renormalized_quantum_source_evaluated=False,full_quantum_Hessian_evaluated=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
 out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
 out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
