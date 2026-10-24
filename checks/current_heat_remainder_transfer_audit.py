"""Boundary/bolt method audit, NOT a full heat remainder certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
from flint import arb,ctx
import current_smeared_heat_projection_domain as previous
root=Path(__file__).resolve().parents[1]
STEM='current-heat-remainder-transfer-audit';PREVIOUS=previous.STEM
OWN=['notes/'+STEM+'.md','checks/current_heat_remainder_transfer_audit.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[]
def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
previous.gy.check=check;previous.prior.scalar.prev.prior.check=check
def digest(p):return sha256((root/p).read_bytes()).hexdigest()
def trunc(p,Y,n=8):return sp.Poly(sp.expand(p),Y).as_dict()
def polycut(p,Y,n):
 return sp.Add(*(val*Y**key[0] for key,val in sp.Poly(sp.expand(p),Y).as_dict().items() if key[0]<=n))

def algebra():
 Y,z=sp.symbols('Y z');b,u,c3,c4,c5,c6=sp.symbols('b u c3 c4 c5 c6',real=True)
 Q=b*b+u*z*z+c3*z**3+c4*z**4+c5*z**5+c6*z**6
 # z''=Q'(z)/2; retain exact independent seed coefficients.
 zz=b*Y
 for n in range(7):
  rhs=polycut(sp.diff(Q,z).subs(z,zz)/2,Y,n).coeff(Y,n)
  zz+=sp.expand(rhs)*Y**(n+2)/((n+1)*(n+2))
 zz=sp.expand(zz);zp=sp.diff(zz,Y)
 check('outer_ODE_jets_through_order_six',polycut(sp.diff(zz,Y,2)-sp.diff(Q,z).subs(z,zz)/2,Y,6)==0)
 fc=[sp.Integer(1)]
 for n in range(1,8):fc.append(sp.expand(-sum(zp.coeff(Y,k)*fc[n-k]/b for k in range(1,n+1))))
 f=sum(fc[n]*Y**n for n in range(8))
 check('reciprocal_warping_jets',polycut(f*zp-b,Y,7)==0)
 v=-polycut(zz**4,Y,7);V=sp.integrate(v,Y)
 df=polycut(-v*f-sp.diff(f,Y)*V,Y,7).expand()
 expected=b**4*Y**4-u*b**4*Y**6/30-3*c3*b**5*Y**7/10
 check('coordinate_correct_conformal_variation',sp.expand(df-expected)==0)
 for n in [0,1,2,3,5]:check('vanishing_variation_jet_'+str(n),df.coeff(Y,n)==0)
 check('first_odd_variation_weight_seven',df.coeff(Y,7)==-3*c3*b**5/10)
 check('observation_shift_order_five',all(V.coeff(Y,k)==0 for k in range(5)))
 check('boundary_weight_count',sp.Rational(7,2)-1==sp.Rational(5,2))
 check('integrated_boundary_power',sp.Rational(5,2)+sp.Rational(1,2)==3)
 rho,g,q2,q3=sp.symbols('rho gamma Q2 Q3',positive=True)
 ss=g*rho**2/2+g*q2*rho**4/48+(g*q2*q2/2880-g*g*q3/480)*rho**6
 residual=sp.diff(ss,rho,2)-(g+q2*ss/2-q3*ss*ss/4)
 check('bolt_ODE_through_order_four',polycut(residual,rho,4)==0)
 A=sp.diff(ss,rho)/g
 check('bolt_warping_series',sp.expand(A-rho-q2*rho**3/12-(q2*q2/480-g*q3/80)*rho**5)==0)
 check('bolt_cartesian_even_regular',sp.limit(A/rho,rho,0)==1 and sp.expand(A.subs(rho,-rho)+A)==0)
 check('radial_log_derivative_diverges',sp.limit(-rho*sp.diff(A,rho)/A,rho,0)==-1)
 check('divergence_independent_of_temperature',sp.diff(sp.limit(-rho*sp.diff(A,rho)/A,rho,0),g)==0)
 j=sp.symbols('j',integer=True,nonnegative=True);m=j+sp.Rational(1,2)
 check('all_thermal_minus_potential',sp.expand(m*m-m-(j*j-sp.Rational(1,4)))==0)
 check('all_thermal_plus_potential',sp.expand(m*m+m-((j+1)**2-sp.Rational(1,4)))==0)
 for shift in [0,1]:
  alpha=j+shift+sp.Rational(1,2)
  check('regular_indicial_family_'+str(shift),sp.expand(alpha*(alpha-1)-((j+shift)**2-sp.Rational(1,4)))==0)
 for label,psi in [('regular',sp.sqrt(rho)),('log',sp.sqrt(rho)*sp.log(rho))]:
  check('critical_channel_'+label,sp.simplify(-sp.diff(psi,rho,2)-psi/(4*rho**2))==0)
 NUMERICS.append(dict(kind='exact_seam_and_bolt_identities',
  conformal_warping_variation=str(df),first_odd_weight=7,
  boundary_diagonal_power='5/2',integrated_power='3',
  bolt_s=str(ss),bolt_A=str(A),radial_log_derivative_residue='-1',
  complete_thermal_potentials=['j^2-1/4','(j+1)^2-1/4'],critical_log_branch_allowed=False,
  scope='Exact algebra supports the ordinary paired-domain/parametrix argument; not a continuum remainder proof.'))

def execute():
 p,h,zh,env=previous.gy.actual_parameters()
 b=p['beta'].sqrt();c3=-4*p['beta']-2*p['u']-p['tau']
 d7=-3*c3*b**5/10
 check('actual_first_odd_variation_nonzero',d7>0)
 Dguard=arb('1777.7385');hguard=arb('1.2e-8');Bguard=arb('0.000021332863')
 check('intrinsic_remainder_guard_is_sufficient',Dguard*hguard<Bguard)
 required=Bguard/hguard
 check('guard_rounded_down',Dguard<required)
 # Retained total heat budget uses the original physical normalization.
 Dtarget=arb('1.0666666e-11');U6=arb('2.3454197e-13');t0=arb('.001')
 cap=U6*t0+hguard*Dguard*t0*t0/2
 check('full_heat_budget_if_intrinsic_guard_proved',cap<Dtarget)
 NUMERICS.append(dict(kind='actual_geometry_and_conditional_target',
  actual_beta=previous.encoded(p['beta']),actual_gamma=previous.encoded(p['gamma']),
  actual_c3=previous.encoded(c3),first_odd_variation_coefficient=previous.encoded(d7),
  D_heat_threshold=previous.encoded(required),D_heat_sufficient='1777.7385',
  global_h_bound='1.2e-8',B_heat_sufficient='0.000021332863',
  conditional_heat_budget=previous.encoded(cap),retained_heat_budget='1.0666666e-11',
  D_heat_numerically_bounded=False,B_heat_numerically_bounded=False,
  intrinsic_bound_attainability='UNKNOWN',direct_bolt_cylinder_transfer='fails bounded-log-derivative hypothesis'))

def main():
 ctx.prec=512
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 review=json.loads((root/OWN[2]).read_bytes())
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=previous.prior.scalar.prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,dependencies=deps,invocation=sys.argv,
  frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'])
 try:
  check('previous_scoped_execution_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('interval_library_pinned',previous.prior.scalar.prev.prior.flint.__version__=='0.9.0')
  algebra();execute()
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=previous.prior.scalar.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  continuum_argument_scope='Ordinary paired symmetry and smooth physical-half parametrix argument; tests verify exact jets and conditional budgets only.',
  C2_action_on_rough_base_family_proved=False,complete_signed_projection_evaluated=False,full_heat_cap_certified=False,
  B_heat_numerically_bounded=False,D_heat_numerically_bounded=False,
  actual_source_failure_proved=False,total_source_improvement_certified=False,nonlinear_solution_constructed=False,
  physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
 out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
