"""Physical bolt source, singular coefficients and pairing contract, not a solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_seed_radial_green_enclosure as prior
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1];STEM='current-bolt-adjoint-domain'
PREVIOUS='current-reference-regularity-feasibility'
OWN=['notes/'+STEM+'.md','checks/current_bolt_adjoint_domain.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[];NUMERICS=[];encoded=prior.encoded
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
prior.check=check

def identities():
 r=sp.symbols('r',positive=True);lap=lambda f:sp.diff(r*sp.diff(f,r),r)/r
 phi=r*r*sp.log(r)
 check('scalar_fundamental_Laplacian',sp.simplify(lap(phi)-(4*sp.log(r)+4))==0)
 check('scalar_fundamental_bolt_flux',sp.limit(r*sp.diff(lap(phi),r),r,0,dir='+')==4)
 check('scalar_fundamental_punctured_equation',sp.simplify(lap(lap(phi)))==0)
 cut=(1-r*r)**4*phi;force=sp.expand(lap(lap(cut)))
 for i,v in enumerate([sp.Integer(1),1+3*r*r-2*r**4,r*r]):
  weak=sp.integrate(sp.expand(r*(lap(cut)*lap(v)-force*v)),(r,0,1))
  check('cutoff_distribution_weak_identity_'+str(i),sp.simplify(weak-4*v.subs(r,0))==0,str(weak))
 check('smooth_bolt_flux_zero',sp.limit(r*sp.diff(lap(1+2*r*r+3*r**4),r),r,0,dir='+')==0)
 k=sp.symbols('k',positive=True)
 check('two_dimensional_point_dual_integral',sp.integrate(2*sp.pi*k/(1+k*k)**2,(k,0,sp.oo))==sp.pi)
 check('J0_quadratic_coefficient',sp.expand(sp.series(sp.besselj(0,r),r,0,4).removeO())==1-r*r/4)
 # The logarithmic moment producing the radial TT leading coefficient.
 L=sp.symbols('L',positive=True)
 check('TT_inverse_logarithmic_moment',sp.simplify(sp.diff(sp.log(sp.log(k)),k)-1/(k*sp.log(k)))==0)
 x,y,xi=sp.symbols('x y xi',real=True);p=sp.Matrix([x,y,0,0]);norm=x*x+y*y
 Pi=sp.eye(4)-p*p.T/norm;orb=-2*xi*(sp.eye(2)-sp.Matrix([x,y])*sp.Matrix([[x,y]])/norm)
 K=sp.diag(0,0,xi,xi)
 for i in range(2):
  for j in range(2):K[i,j]=orb[i,j]
 check('bolt_TT_lift_transverse',sp.simplify(K*p)==sp.zeros(4,1))
 check('bolt_TT_lift_trace',sp.simplify(sp.trace(K))==0)
 t2=sp.simplify(sp.trace((K*Pi)**2));t1=sp.simplify(sp.trace(K*Pi))
 check('bolt_TT_lift_square',t2==6*xi*xi and t1==0)
 coeff=sp.simplify((3*t2-t1*t1)/(960*sp.pi**2*xi**2))
 check('TT_logarithmic_per_volume_coefficient',coeff==3/(160*sp.pi**2))
 N,C,R,G,D,h=sp.symbols('N C R G D h',positive=True)
 mu=8*sp.pi**2*N*R**4/C;scale=C*C/R**4;ks=6*sp.pi*N*D*C/G;kt=3*N*C/20;ell=(1+h)/h
 check('scalar_measure_normalization',sp.simplify(mu*scale*6*D/(8*sp.pi*G)-ks)==0)
 check('TT_measure_normalization',sp.simplify(mu*scale*coeff-kt)==0)
 cs=-ell/(4*ks);ct=-ell/(8*kt)
 check('scalar_point_source_coefficient',sp.simplify(-4*ks*cs-ell)==0)
 check('TT_point_source_coefficient',sp.simplify(-4*kt*ct-ell/2)==0)
 check('folded_scalar_source_invariant',sp.simplify(2*cs.subs(N,2)-cs.subs(N,1))==0)
 check('folded_TT_source_invariant',sp.simplify(2*ct.subs(N,2)-ct.subs(N,1))==0)
 v,u,e,b=sp.symbols('v u e b',positive=True)
 ext=(v-1)**2*(b+2*b/v+(3*b+u)/v**2+e/v**4)
 P=v*ext;temperature=sp.diff(P,v)/(2*v)
 kval=(v+1)*(u+3*b*(v*v+1))/(2*v**3)-e*(v-3)/(2*v**5)
 check('exact_fixed_temperature_polynomial',sp.simplify(temperature-(v-1)*kval)==0)
 check('temperature_bolt_derivative_one',sp.simplify(kval.subs(v,1).subs(b,(1-u-e)/6))==1)
 B0,V0=sp.symbols('B0 V0',positive=True)
 radius=B0+B0*h*sp.exp(V0)
 check('physical_bolt_radius_functional',sp.simplify(sp.diff(radius,V0).subs(V0,0)/(B0*h))==1 and sp.simplify(ell*h/(1+h))==1)
 FIXTURES.append(dict(fixture='physical_point_source_measure_and_singular_profiles',full_adjoint_computed=False))

def current_numbers():
 p=prior.parameters();b,u,e,gam=p['beta'],p['u'],p['eps'],p['gamma']
 def K(v):return (v+1)*(u+3*b*(v*v+1))/(2*v**3)-e*(v-3)/(2*v**5)
 def Kp(v):
  num=(v+1)*(u+3*b*(v*v+1));dn=u+3*b*(v*v+1)+(v+1)*6*b*v
  return dn/(2*v**3)-3*num/(2*v**4)-e/(2*v**5)+5*e*(v-3)/(2*v**6)
 lower=gam.lower()*(1-arb('1e-28'));upper=gam.upper()*(1+arb('1e-28'))
 fl=lower*K(1+lower)-gam;fu=upper*K(1+upper)-gam
 check('current_temperature_root_lower',fl<0,encoded(fl))
 check('current_temperature_root_upper',fu>0,encoded(fu))
 hb=prior.hull(lower,upper);v=1+hb;hmax=arb(2000)**-10
 derivative=K(v)+hb*Kp(v)
 check('current_temperature_root_monotonicity',derivative>arb('.999999'),encoded(derivative))
 check('retained_temperature_h_guards',hb>hmax/8 and hb<hmax)
 tau=v*hb**2*(b*(1+2/v+3/v**2)+u/v**2+e/v**4)
 check('retained_thermal_polynomial_guard',tau>0 and tau<arb('1e-65'),encoded(tau))
 ell=(1+hb)/hb;DR=96*arb.pi()*arb('1e22')-1/(60*arb.pi())
 NUMERICS.append(dict(kind='actual_current_bolt_normalization',h=encoded(hb),gamma=encoded(gam),h_over_gamma=encoded(hb/gam),tau=encoded(tau),physical_radius_factor=encoded(ell),D_R_in_G_equals_one_units=encoded(DR),root_lower_sign=encoded(fl),root_upper_sign=encoded(fu)))
 for N in [1,2]:
  ks=6*arb.pi()*N*DR*p['C'];kt=3*N*p['C']/20;cs=-ell/(4*ks);ct=-ell/(8*kt)
  check('scalar_singularity_nonzero_'+str(N),cs<arb('-1e-5') and cs>arb('-.01'),encoded(cs))
  check('TT_singularity_nonzero_'+str(N),ct<arb('-1e20') and ct>arb('-1e25'),encoded(ct))
  check('normalized_scalar_flux_'+str(N),(-4*ks*cs/ell).contains(1))
  check('normalized_TT_source_'+str(N),(-8*kt*ct/ell).contains(1))
  NUMERICS.append(dict(kind='actual_adjoint_leading_coefficients',identical_exteriors=N,K_scalar=encoded(ks),K_TT_log=encoded(kt),scalar_rho2_log_coefficient=encoded(cs),TT_rho2_loglog_coefficient=encoded(ct),global_adjoint_or_inverse_norm=False,finite_radius_asymptotic_error_bound=False))
 # Missing 1/log(k) after a leading TT subtraction cannot be an L2 residual.
 values=[arb(10)**(2*j)/((arb(2)*arb(10)**j).log())**2 for j in [3,6,12]]
 check('TT_logarithmic_residual_annular_growth',values[1]>values[0] and values[2]>values[1])
 NUMERICS.append(dict(kind='TT_leading_subtraction_L2_diagnostic',annular_lower_bounds=[encoded(v) for v in values],ordinary_fourth_contact_may_not_be_dropped=True,actual_residual_norm_computed=False))
 FIXTURES.append(dict(fixture='current_temperature_and_adjoint_coefficients',geometry_and_physical_normalization_unchanged=True))

def pairing_fixture():
 H=sp.Matrix([[-7,2],[2,5]]);r=sp.Matrix([3,-11]);ell=sp.Matrix([5,sp.Rational(5,2)])
 psi=H.inv()*ell;solution=-H.inv()*r
 check('coupled_indefinite_pairing_identity',(ell.T*solution+psi.T*r)[0]==0)
 app=psi+sp.Matrix([sp.Rational(1,10),sp.Rational(-1,7)]);defect=ell-H*app
 error=abs((ell.T*solution+app.T*r)[0]);dual=lambda v:v[0]**2/7+v[1]**2/5
 check('coupled_adjoint_residual_error_certificate',error**2<=dual(r)*dual(defect))
 check('unit_action_form_inverse_for_fixture',(solution[0]**2*7+solution[1]**2*5)<=dual(r))
 FIXTURES.append(dict(fixture='independent_finite_dimensional_pairing_and_error_identity',actual_quantum_pairing=False))

def main():
 ctx.prec=256
 review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=prior.dependency_identity();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
 try:
  check('previous_regularity_execution_passed',old['status']=='passed' and old['principal_derivative_quadrature_target_certified'])
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('pinned_interval_library',prior.flint.__version__=='0.9.0')
  identities();current_numbers();pairing_fixture()
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=prior.dependency_identity();status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),physical_bolt_adjoint_source_identified=status=='passed',current_bolt_scalar_and_TT_leading_coefficients_evaluated=status=='passed',smooth_homogeneous_bolt_adjoint_rejected=status=='passed',TT_leading_subtraction_alone_L2_residual_rejected=status=='passed',full_bolt_adjoint_computed=False,full_defect_bolt_pairing_evaluated=False,full_adjoint_residual_norm_bounded=False,physical_bolt_functional_dual_norm_quantified=False,all_local_renormalized_moments_computed=False,coincident_subtraction_executed=False,differentiated_response_operator_bounded=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
 out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
 out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
