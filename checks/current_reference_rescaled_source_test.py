"""Changed physical reference correction and partial joint source screen."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
import sympy as sp
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'checks'))
import current_schur_scalar_projection as prev
import current_seed_cold_action_defect as cold
from flint import arb,ctx
STEM='current-reference-rescaled-source-test';PREVIOUS='current-schur-scalar-projection'
OWN=['notes/'+STEM+'.md','checks/current_reference_rescaled_source_test.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];encoded=prev.encoded;decode=prev.prev.endpoint.raw.decode

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
prev.prev.prior.check=check;prev.prev.york.check=check;cold.check=check

def digest(p):return sha256((root/p).read_bytes()).hexdigest()

def algebra():
 cold.identities()
 e,j,h,f=sp.symbols('epsilon J2 h f');F=1-e
 # I0''h=-f, with H=F I0''+J''. This calculation retains the seed f.
 old=e*f+j*h;new=f+(-F*f+j*h)/F
 check('full_two_row_reference_residual_identity',sp.simplify(new-j*h/F)==0)
 check('seed_subtraction_in_both_rows',sp.simplify(new-(old-e*f)/F)==0)
 chi,du,r,S=sp.symbols('chi du r S')
 check('principal_full_stress_cancellation',sp.expand(-2*chi*(-12+du*r)*S-24*chi*S+2*chi*du*r*S)==0)
 z=sp.symbols('z',positive=True)
 check('global_remaining_principal_moment',sp.integrate(4*z**5+6*z**6,(z,0,1))==sp.Rational(32,21))
 # A two-dimensional invertible saddle fixture checks the SAME Newton correction.
 H0=sp.diag(2,-3);J=sp.Matrix([[sp.Rational(1,5),sp.Rational(1,7)],[sp.Rational(1,7),sp.Rational(1,11)]])
 ff=sp.Matrix([2,3]);F=sp.Rational(9,10);hh=-H0.inv()*ff;HH=F*H0+J
 ro=ff+HH*hh;hn=hh/F;rn=ff+HH*hn
 check('physical_linear_Newton_correction_unchanged',hn-HH.inv()*rn==hh-HH.inv()*ro)
 check('physical_trial_actually_changes',hn!=hh)

def execute(old):
 p,h,zh,env=prev.prev.york.actual_parameters();NUMERICS.extend(prev.prev.york.NUMERICS)
 C=p['C'];B2=p['B2'];lam=p['g'];x=lam*B2;z=prev.prev.local.monomial
 du=p['eps']*C/x;ee=1/(60*arb.pi());eta=ee*lam;DR=96*arb.pi()*arb('1e22');dl=(DR-ee)*lam
 chi=16*arb.pi()*arb('1e22')*lam;eR=24*chi;Fl=1-eR
 beta=x/C;v=1+h;thermal=v*h*h*(beta*(1+2/v+3/v**2)+1/(C*v**2)+p['eps']/v**4)
 fa=4+2/x+thermal*C/x;fb=3+1/x+du
 Fn=z(-2)+1/x-fa*z(1)+fb*z(2)-2*du*z(3)+du*z(4)
 Dt=lambda f:-z(2)*f.derivative()
 Lap=lambda f:z(4)*(Fn*f.derivative()).derivative()
 r=4*z(5)-6*z(6);rn=-12+du*r
 tn=(3+1/x+du)*z(4)-du*(5*z(5)-arb(9)/2*z(6))
 cn=fa*z(3)-2*fb*z(4)+arb(20)/3*du*z(5)-5*du*z(6)
 P=7*cn**2+arb(22)/3*tn**2-arb(11)/72*(rn**2-144)
 Ps=2*(6+du)**2+arb(22)/3*(6+du)/x+arb(22)/3/x**2+22
 S=ee*arb(p['q'])**2/(lam*B2*B2)
 tau=S*(-10*z(4)-20*z(5)+30*z(6))+eta*(P+Ps*r/2)-dl*du*Lap(r)
 K=dl*Lap(tau)+eta*(arb(11)/36*rn-arb(7)/3*cn+arb(11)/3*tn)*tau
 On=rn/2-2*tn
 shape_poly=chi*(x/C)*(2*On*tau+2*Fn*Dt(Dt(tau))+2*(Dt(Fn)+4*z(1)*Fn)*Dt(tau))
 BT=shape_poly.integral(zh)
 scalar=next(row for row in old['numerics'] if row['kind']=='actual_canonical_scalar_projection')
 sr=next(row for row in old['numerics'] if row['kind']=='actual_global_scalar_inverse')
 prev.power=decode(sr['fractional_power']);prev.y0=h/(1+h);prev.moments.clear()
 trial=prev.FP({(c['z_power'],c['y_integer_power'],c['fractional_multiplicity']):decode(c['coefficient']) for c in sr['trial_coefficients']})
 vt=trial+prev.mon(4)/2
 pairing=lambda poly:(x/C)*sum((c*(prev.mon(k-4)*vt).integral() for (k,l),c in poly.c.items()),arb(0))
 IK=pairing(K);Iseed=pairing(tau);PT=(pairing(K-eR*tau)-BT)/Fl
 check('same_previous_scalar_known_pairing',IK.overlaps(decode(scalar['known_areal_pairing'])))
 check('coefficient_preserving_seed_subtraction',PT.overlaps((IK-eR*Iseed-BT)/Fl))
 tmajor=prev.absolute(tau)
 check('seed_trace_norm_outer_integrability',all(k>=4 and j==t==0 for k,j,t in tmajor.c))
 seed_norm=(x/C)*(prev.mon(-4)*tmajor**2).integral().sqrt()
 inverse_error=(arb('2.3300059e-44')+eR*seed_norm*arb('3.6891662e-12'))/Fl
 known_allow=PT.abs_upper()+inverse_error.upper()
 principal=4*chi*du.abs_upper()*(x/C)*arb('2.01e-7')*arb(32)/21/Fl
 bolt_change=eR/Fl*arb('3.33e-10');bolt=arb('3.33e-10')/Fl
 old_ceiling=arb('2.6035780e-40');working=arb('2.60e-40')
 budget=arb('2e-9')-bolt
 safe_old_threshold=old_ceiling*arb('1.66699e-9')/arb('1.667e-9')
 contacts=(decode(scalar['unknown_areal_upper'])+decode(scalar['unknown_density_transport_endpoint_upper']))/Fl
 check('reference_factor_positive',Fl>arb('.999999988'))
 check('physical_bolt_change_target',bolt_change<arb('4e-18'),str(bolt_change))
 check('remaining_physical_bolt_budget',budget>arb('1.66699e-9'))
 check('working_screen_conservative_after_trial_change',working<safe_old_threshold)
 check('known_joint_trace_component_target',known_allow<arb('1e-42'),str(known_allow))
 check('remaining_principal_stress_target',principal<arb('1e-43'),str(principal))
 check('joint_known_allocation',known_allow+principal<arb('1.1e-42'))
 check('remaining_signed_allocation',sp.Rational('2.589e-40')+sp.Rational('1.1e-42')==sp.Rational('2.60e-40'))
 check('inherited_contact_allowance_still_exceeds_target',contacts>22*working)
 NUMERICS.append(dict(kind='physical_reference_prefactor_test',F_lambda=encoded(Fl),epsilon_R=encoded(eR),
  old_first_bolt_upper=encoded(arb('3.33e-10')),physical_trial_bolt_change_upper=encoded(bolt_change),
  new_first_bolt_upper=encoded(bolt),new_remaining_bolt_budget_lower=encoded(budget),
  conservative_working_screen=encoded(working),scaled_old_source_screen=encoded(safe_old_threshold),
  actual_first_correction_changed=True,physical_Hessian_changed=False,action_or_normalization_changed=False))
 NUMERICS.append(dict(kind='known_joint_local_source_test',x_over_C=encoded(x/C),zh=encoded(zh),
  known_scalar_old=encoded(IK),known_R2_trace_shape=encoded(BT),old_known_difference=encoded(IK-BT),
  seed_scalar_pairing=encoded(Iseed),rescaled_known_trace_component=encoded(PT),
  seed_trace_mass_norm_upper=encoded(seed_norm),numerical_scalar_error_upper=encoded(inverse_error),
  known_trace_component_absolute_upper=encoded(known_allow),remaining_principal_stress_upper=encoded(principal),
  sum_closed_components_upper=encoded(known_allow+principal),
  remaining_signed_source_target='2.589e-40',
  inherited_scalar_contact_upper_allowance=encoded(contacts),contact_allowance_target_factor=encoded(contacts/working),
  trial=vt.serial(),fractional_power=encoded(prev.power),
  scaled_scalar_integrand_coefficients=prev.convert((x/C)*(K-eR*tau)).serial(),
  shape_integrand_coefficients=prev.convert(shape_poly).serial(),
  principal_coefficient=encoded(4*chi*du.abs_upper()*(x/C)*arb('2.01e-7')/Fl),
  complete_shape_source_evaluated=False,complete_signed_reduced_source_evaluated=False))
 print(json.dumps(dict(event='partial_joint_test',known_trace_upper=str(known_allow),principal_upper=str(principal),
  bolt_change_upper=str(bolt_change),complete_projection_evaluated=False)),flush=True)

def main():
 ctx.prec=512
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 review=json.loads((root/OWN[2]).read_bytes())
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=prev.prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'],dependencies=deps,
  invocation=sys.argv,previous_goal_turn='progress')
 try:
  check('previous_global_scalar_test_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('pinned_interval_library',prev.prev.prior.flint.__version__=='0.9.0')
  algebra();execute(old)
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=prev.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  physical_reference_prefactor_subtraction_verified=status=='passed',known_local_joint_projection_certified=status=='passed',
  remaining_principal_stress_coefficient_certified=status=='passed',complete_signed_reduced_source_evaluated=False,
  total_source_improvement_certified=False,full_solution_feasibility_test_executed=False,
  nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,
  component_complete=False,goal_complete=False,formal_verification=False)
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,
              dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
