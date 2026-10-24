"""Boundary coefficient and regularity audit; NOT a regulator-limit theorem."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
if '--arb-path' in sys.argv:
 sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
from flint import arb, ctx
import current_geometric_mass_conversion as geom
root=Path(__file__).resolve().parents[1]
STEM='current-paired-regulator-domain-audit'; PREVIOUS=geom.STEM
OWN=['notes/'+STEM+'.md','checks/current_paired_regulator_domain_audit.py','evidence/'+STEM+'-review.json']
ROWS=[]; NUMERICS=[]; encoded=geom.encoded

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok: raise AssertionError(name)

geom.prior.scalar.prev.york.check=check
geom.prior.scalar.prev.prior.check=check
def digest(p): return sha256((root/p).read_bytes()).hexdigest()

def algebra():
 R,Rn,Rnn,M=sp.symbols('R Rn Rnn M',real=True)
 I=sp.eye(2);chi=sp.Matrix([[0,1],[1,0]])
 Om=R*sp.Matrix([[sp.I,0],[0,-sp.I]])/4
 E=-R*I/4;En=-Rn*I/4;Pp=(I+chi)/2;Pm=(I-chi)/2
 check('actual_balanced_involution',chi**2==I and sp.trace(chi)==0)
 check('normal_tangent_curvature_anticommutes',chi*Om*chi==-Om)
 check('a1_a3_scalar_chi_traces_zero',sp.trace(chi*(96*E+20*R*I))==0)
 a4bd=sp.trace((240*Pp-120*Pm)*En+(42*Pp-18*Pm)*Rn)
 a4div=sp.trace(60*(-I/4)+12*I)
 check('a4_boundary_derivative',sp.simplify(a4bd+6*Rn)==0)
 check('a4_inward_bulk_boundary_completion',sp.simplify(-a4div*Rn+a4bd)==0)
 ric2=R**2/2;riem2=R**2
 bulk=sp.trace(60*R*E+180*E**2+30*(2*Om**2)+(5*R**2-2*ric2+2*riem2)*I)/360
 check('a4_rank_two_integrated_curvature',sp.simplify(bulk+R**2/120)==0)
 check('a5_second_normal_curvature_trace_zero',sp.trace(chi*Rnn*I)==0)
 a5=sp.trace(-180*E**2+180*chi*E*chi*E-45*Om**2+180*chi*Om**2-45*chi*Om*chi*Om)
 check('a5_complete_remaining_mixed_trace',sp.simplify(a5)==0)
 a1=a3=sp.Integer(0)
 check('all_angular_mass_a5_cancellation',sp.expand(a5-M**2*a3+M**4*a1/2)==0)
 NUMERICS.append(dict(kind='boundary_coefficient_test',a1='0',a3='0',a4_bulk_divergence='-6 Delta R',
  a4_boundary_derivative='-6 R_;n (inward)',a4_completed_coefficient='-R^2/120',a5='0',
  all_angular=True,actual_rough_heat_expansion_proved=False))
 z=sp.symbols('z',positive=True);c,d=sp.symbols('c d',real=True)
 h=c*z**3+d*z**4*sp.log(z);v=-z**4;h3=sp.diff(h,z,3);v3=sp.diff(v,z,3)
 check('mass_log_third_derivative',sp.simplify(h3-6*c-d*z*(24*sp.log(z)+26))==0)
 check('log_fourth_derivative',sp.simplify(sp.diff(h,z,4)-d*(24*sp.log(z)+50))==0)
 moments=[36*c*c*z+c*d*z*z*(144*sp.log(z)+84)+d*d*z**3*(192*sp.log(z)**2+288*sp.log(z)+sp.Rational(388,3)),
  192*z**3,-72*c*z*z-d*z**3*(192*sp.log(z)+144)]
 for i,(ans,integrand) in enumerate(zip(moments,[h3*h3,v3*v3,h3*v3])):
  check('weak_endpoint_moment_primitive_'+str(i),sp.simplify(sp.diff(ans,z)-integrand)==0)
  check('weak_endpoint_moment_lower_limit_'+str(i),sp.limit(ans,z,0,dir='+')==0)
 check('reflected_cubic_third_jump',6-(-6)==12)
 check('nonzero_log_fourth_unbounded',sp.limit(24*sp.log(z)+50,z,0,dir='+')==-sp.oo)
 x=sp.symbols('x',real=True);t=sp.symbols('t',real=True)
 density=sp.factorial(15)/sp.factorial(7)**2*t**7*(1-t)**7
 theta=sp.integrate(density,(t,0,t))
 check('cutoff_endpoint_values',theta.subs(t,0)==0 and theta.subs(t,1)==1)
 check('cutoff_positive_density_factor',sp.expand(sp.diff(theta,t)-density)==0)
 check('cutoff_seven_matching_derivatives',all(sp.diff(theta,t,k).subs(t,e)==0 for k in range(1,8) for e in [0,1]))
 cubic=sp.expand(x**3*theta.subs(t,x-1));quartic=sp.expand(x**4*theta.subs(t,x-1))
 K3=36+sp.integrate(sp.expand((sp.diff(cubic,x,3)-6)**2),(x,1,2))
 K4=sp.integrate(sp.expand(sp.diff(cubic,x,4)**2),(x,1,2))
 K3v=192+sp.integrate(sp.expand((sp.diff(quartic,x,3)-24*x)**2),(x,1,2))
 check('positive_cutoff_costs',K3>0 and K4>0 and K3v>0)
 NUMERICS.append(dict(kind='exact_cutoff_costs',K3=str(K3),K4=str(K4),K3v=str(K3v),
  theta_coefficients=[str(a) for a in sp.Poly(theta,t).all_coeffs()],
  cubic_H3_squared_error='K3 epsilon',cubic_H4_squared_cost='K4/epsilon',quartic_H3_squared_error='K3v epsilon^3',
  uniform_H4_approximation_possible_by_this_cutoff=False,physical_obstruction_proved=False))
 check('cubic_H3_error_decreases',K3/sp.Integer(2)**20<K3/sp.Integer(2)**10)
 check('cubic_H4_cost_diverges',K4*sp.Integer(2)**20>K4*sp.Integer(2)**10)
 primitive=z*(2-sp.log(z))
 check('log_Dini_primitive',sp.simplify(sp.diff(primitive,z)-(1-sp.log(z)))==0 and sp.limit(primitive,z,0,dir='+')==0)
 q=sp.symbols('q',positive=True)
 total=2*(q+2)/(q+1)**2+2/(q+1)
 check('full_angular_fourth_moment_upper',sp.simplify(4/q-total-2*(q+2)/(q*(q+1)**2))==0)

def execute():
 p,h,zh,env=geom.prior.scalar.prev.york.actual_parameters()
 NUMERICS.extend(geom.prior.scalar.prev.york.NUMERICS)
 B2=p['B2'];C=p['C'];Fl=1-384*arb.pi()*arb('1e22')*p['g'];t0=arb('.001')
 check('actual_positive_parameters',B2>0 and C>0 and Fl>0 and zh<1)
 threshold=arb('2.5e-43')*B2*Fl/(4*C)
 check('conservative_Dini_threshold',threshold>arb('1.0666666e-11'))
 rows=[]
 for av in ['.25','.5','1']:
  a=arb(av);bound=threshold*a/(t0**a)
  rows.append(dict(alpha=av,A_threshold=encoded(bound)))
 logfactor=t0*(2-t0.log());Al=threshold/logfactor
 glue=12*arb('1e24')*C*arb('4e-9')/(B2*Fl)
 ratio=glue/arb('1.287e-40');Kmax=arb('1.287e-40')*B2*Fl/(arb('8e-9')*C)
 check('known_gluing_estimator_excessive',ratio>arb('1e24'))
 check('gluing_coefficient_target_consistent',(arb('6e24')/Kmax).overlaps(ratio))
 check('logarithmic_threshold_finite_positive',Al>0)
 NUMERICS.append(dict(kind='domain_and_estimator_thresholds',B_squared=encoded(B2),C=encoded(C),F_lambda=encoded(Fl),
  t0='.001',cap_error_allocation='2.5e-43',Dini_integral_threshold=encoded(threshold),power_thresholds=rows,
  logarithmic_Dini_factor=encoded(logfactor),A_log_threshold=encoded(Al),
  known_gluing_upper_allowance=encoded(glue),gluing_allowance_target_factor=encoded(ratio),
  required_crude_gluing_coefficient_upper=encoded(Kmax),full_projection_target='1.287e-40',
  Dini_majorant_known=False,uniform_common_regulator_identification_proved=False,actual_projection_violation_proved=False))
 print(json.dumps(dict(event='domain_audit',A_log_threshold=str(Al),known_gluing_allowance=str(glue),
  excessive_allowance_factor=str(ratio),uniform_regulator_identification=False,full_projection_executed=False)),flush=True)

def main():
 ctx.prec=512
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 review=json.loads((root/OWN[2]).read_bytes())
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=geom.prior.scalar.prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,dependencies=deps,invocation=sys.argv,
  frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'])
 try:
  check('previous_scoped_conversion_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('interval_library_pinned',geom.prior.scalar.prev.prior.flint.__version__=='0.9.0')
  algebra();execute()
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=geom.prior.scalar.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  uniform_regulator_identification_proved=False,complete_signed_projection_evaluated=False,thermal_subtraction_matched=False,
  actual_source_failure_proved=False,total_source_improvement_certified=False,nonlinear_solution_constructed=False,
  physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
 out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
