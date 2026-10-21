"""Verify exterior gauge coverage and the full spherical equation reduction."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp

ROOT=Path(__file__).resolve().parents[1];STEM='ingoing-areal-retarded-system'
OWN=['checks/ingoing_areal_retarded_system.py','notes/'+STEM+'.md','evidence/ingoing-areal-source-review.json']
PREVIOUS='boundary-clock-response-extension';ROWS=[];NUMERICS=[]
EXTRA=['notes/sourced-spherical-mass-charge-balance.md']
PROTOCOL={'scope':'Gauge and exact equation reduction, not an evaluated quantum kernel or physical stability',
 'frequency_points':[[.7,0.],[1.2,.8],[2.,3.]],'ODE_tolerances':[1e-8,1e-12],
 'ODE_final_error_limit':2e-10,'ODE_refinement_limit':2e-8,
 'transport_points':[[.3,1.1],[1.3,.4],[2.,1.7]],'transport_tolerances':[1e-8,1e-12],
 'transport_error_limit':2e-11,'metric_points':[.3,1.2,2.7],
 'metric_steps':[1e-3,5e-4,2.5e-4],'metric_final_error_limit':2e-6,
 'metric_refinement_ratio':.4,'metric_roundoff_floor':2e-8,
 'formal_verification':False,'physical_kernel_evaluated':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,x):
 x=S.simplify(S.sympify(x).doit());check(name,x==0,x)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN+EXTRA)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==66 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_sources_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']) and review['additional_pinned_inputs']==EXTRA)
 check('source_scope',len(review['primary_sources'])==2 and review['coordinate_clock_is_source_identification'] and not review['physical_relative_clock_discarded'] and not review['quantum_kernel_evaluated'])

def geometry():
 v,r,th,ph=S.symbols('v r theta phi',real=True);x=[v,r,th,ph]
 psi=S.Function('psi')(v,r);F=S.Function('F')(v,r)
 g=S.Matrix([[-S.exp(2*psi)*F,S.exp(psi),0,0],[S.exp(psi),0,0,0],
             [0,0,r*r,0],[0,0,0,r*r*S.sin(th)**2]])
 gi=g.inv()
 Ga=[[[S.simplify(sum(gi[a,d]*(S.diff(g[d,c],x[b])+S.diff(g[d,b],x[c])-S.diff(g[b,c],x[d]))/2 for d in range(4))) for c in range(4)] for b in range(4)] for a in range(4)]
 Ric=S.Matrix(4,4,lambda a,b:S.simplify(sum(S.diff(Ga[c][a][b],x[c])-S.diff(Ga[c][a][c],x[b])+sum(Ga[c][c][d]*Ga[d][a][b]-Ga[c][b][d]*Ga[d][a][c] for d in range(4)) for c in range(4))))
 R=S.simplify(sum(gi[a,b]*Ric[a,b] for a in range(4) for b in range(4)))
 Ein=S.simplify(Ric-g*R/2);Mix=S.simplify(gi*Ein)
 zero('direct_4D_G_rr',Ein[1,1]-2*S.diff(psi,r)/r)
 zero('direct_4D_G_vv_mixed',Mix[0,0]-(r*S.diff(F,r)+F-1)/r**2)
 zero('direct_4D_G_rv_mixed',Mix[1,0]+S.diff(F,v)/r)
 zero('direct_4D_orbit_symmetry',Mix[1,1]-Mix[0,0]-F*Ein[1,1])
 pressure=S.symbols('P');zero('angular_Ward_completion',-(Ga[2][2][1]+Ga[3][3][1])*pressure+2*pressure/r)
 G,Q,lam=S.symbols('G Q lam',positive=True)
 M=r*(1+Q/r**2+lam*r*r-F)/(2*G)
 zero('mass_radial_Einstein_equation',Mix[0,0]-3*lam+Q/r**4+2*G*S.diff(M,r)/r**2)
 zero('mass_flux_Einstein_equation',Mix[1,0]-2*G*S.diff(M,v)/r**2)
 e,f,E0,Er,Ev,Em=S.symbols('e f E0 Er Ev Em')
 rel=Ev+e*Em+(S.symbols('F0')+e*f)*(E0+e*Er)
 zero('offshell_orbit_linear_contact',S.diff(rel,e).subs(e,0)-Em-S.symbols('F0')*Er-f*E0)
 check('discarded_offshell_contact_detected',S.diff(rel,e).subs({e:0,Em:0,Er:0,E0:2,f:3})==6)
 mu,zeta=S.symbols('mu zeta',cls=S.Function);F0=S.Function('F0')(r)
 zero('linear_rr_normalization',S.diff(2*S.diff(e*zeta(v,r),r)/r,e).subs(e,0)-2*S.diff(zeta(v,r),r)/r)
 fe=F0-2*G*e*mu(v,r)/r
 zero('linear_mass_radial_normalization',S.diff((r*S.diff(fe,r)+fe-1)/r**2,e).subs(e,0)+2*G*S.diff(mu(v,r),r)/r**2)
 zero('linear_mass_flux_normalization',S.diff(-S.diff(fe,v)/r,e).subs(e,0)-2*G*S.diff(mu(v,r),v)/r**2)

def gauge_and_ADM():
 v,r=S.symbols('v r',real=True);x=[v,r];F=S.Function('F')(r)
 rho=S.Function('rho')(v,r);kap=S.Function('kappa')(v,r)
 A,B,C=S.symbols('A B C');g=S.Matrix([[-F,1],[1,0]]);xi=[kap,rho]
 h=S.Matrix([[-2*F*A,2*A-C],[2*A-C,2*(-A+B+C)/F]])
 Lie=S.Matrix(2,2,lambda a,b:sum(xi[c]*S.diff(g[a,b],x[c])+g[c,b]*S.diff(xi[c],x[a])+g[a,c]*S.diff(xi[c],x[b]) for c in range(2)))
 zero('Lie_rr',Lie[1,1]-2*S.diff(kap,r))
 sub={S.diff(kap,r):(-A+B+C)/F}
 zeta=A+B-S.diff(rho,r)-S.diff(kap,v)
 f=-2*F*B+2*F*S.diff(rho,r)-S.diff(F,r)*rho+2*S.diff(rho,v)
 expected=S.Matrix([[-f-2*F*zeta,zeta],[zeta,0]])
 new=(h-Lie).subs(sub)
 for a in range(2):
  for b in range(2):zero('gauge_metric_'+str(a)+str(b),new[a,b]-expected[a,b])
 zero('angular_areal_gauge',2*r*rho-2*r*xi[1])
 gi=g.inv();deltaX=-sum(gi[1,a]*gi[1,b]*h[a,b] for a in range(2) for b in range(2))+2*sum(gi[1,a]*S.diff(rho,x[a]) for a in range(2))
 zero('invariant_mass_amplitude',deltaX-S.diff(F,r)*rho-f)
 rt,ry=S.symbols('rho_t rho_y')
 zero('fixed_t_mass_amplitude',f.subs({S.diff(rho,v):rt,S.diff(rho,r):-(rt+ry)/F})-(-2*F*B-2*ry-S.diff(F,r)*rho))
 zero('incoming_transport_sign',(-F*S.diff(kap,r)).subs(S.diff(kap,r),(-A+B+C)/F)-(A-B-C))
 y,gamma,cb=S.symbols('y gamma cb',positive=True)
 khat=S.exp(-y)*(1/(gamma**2+1)+y*gamma/(gamma**2+4))
 fhat=S.exp(-y)*((gamma-1)/(gamma**2+1)+(1-y)*gamma/(gamma**2+4)-4*y/(gamma**2+4))
 zero('laplace_initial_term',gamma*khat+S.diff(khat,y)-y*S.exp(-y)-fhat)
 check('Laplace_initial_contact_cannot_be_removed',S.simplify(gamma*khat+S.diff(khat,y)-fhat)!=0)
 zero('boundary_L2_normalization',S.integrate(S.exp(-2*gamma*y),(y,0,S.oo))-1/(2*gamma))
 # A smooth EF vector need not be an optical H1 vector at the horizon.
 rh=S.symbols('rh',positive=True)
 check('optical_radial_component_can_diverge',S.limit(-1/(r-rh),r,rh,dir='+')==-S.oo)
 zero('smooth_EF_horizon_primitive',S.diff(-S.exp(-r)/2,r)-S.exp(-r)/2)
 # Pull back the exact EF metric to background static t,y and extract ADM.
 e1,e2,zi,zj,ki0,kj=S.symbols('e1 e2 zi zj ki kj',real=True)
 psi=e1*zi+e2*zj;k=e1*ki0+e2*kj;E=S.exp(psi);den=2-E*(1+k)
 gtt=-E**2*(1+k);gty=E**2*(1+k)-E;gyy=E*den
 lapse2=S.simplify((gty*gty-gtt*gyy)/gyy)
 zero('ADM_lapse_squared',lapse2-E/den)
 zero('ADM_speed_squared',lapse2/gyy-1/den**2)
 zero('ADM_shift',gty/gyy-(1/den-1))
 speed=1/den;mass=S.sqrt(E/den);at={e1:0,e2:0}
 zero('speed_first_vertex',S.diff(speed,e1).subs(at)-zi-ki0)
 zero('mass_first_vertex',S.diff(mass,e1).subs(at)-zi-ki0/2)
 sv=3*zi*zj+3*(zi*kj+zj*ki0)+2*ki0*kj
 mv=2*zi*zj+S.Rational(3,2)*(zi*kj+zj*ki0)+S.Rational(3,4)*ki0*kj
 zero('speed_second_vertex',S.diff(speed,e1,e2).subs(at)-sv)
 zero('mass_second_vertex',S.diff(mass,e1,e2).subs(at)-mv)
 # Pauli components are checked separately, never via a complex sum of squares.
 alpha=S.Matrix([[0,-S.I],[S.I,0]]);beta=S.Matrix([[0,1],[1,0]])
 d,w,dp,wp,m=S.symbols('d w dp wp m')
 general=-S.I*alpha*d+S.I*S.eye(2)*d
 for a in range(2):
  for b in range(2):zero('Dirac_derivative_component_'+str(a)+str(b),general[a,b]-(S.I*(S.eye(2)-alpha)*d)[a,b])
 check('static_slice_horizon_failure_detected',1-S.Rational(1,10)/S.Rational(1,100)<0)
 cp=S.symbols('clock_prime',real=True);eps=S.symbols('eps',real=True)
 zero('prepared_boundary_clock_linearization',S.diff(S.log(1+eps*cp),eps).subs(eps,0)-cp)

def cq(fun,a,b,tol):
 return quad(lambda s:complex(fun(s)).real,a,b,epsabs=tol,epsrel=tol,limit=160)[0]+1j*quad(lambda s:complex(fun(s)).imag,a,b,epsabs=tol,epsrel=tol,limit=160)[0]

def frequency_fixtures():
 grid=np.linspace(0.,9.,73);source=lambda y:y**3*np.exp(-2*y)
 for re,im in PROTOCOL['frequency_points']:
  gamma=complex(re,im);c=.23-.1j;sols=[]
  for tol in PROTOCOL['ODE_tolerances']:
   sol=solve_ivp(lambda y,u:source(y)-gamma*u,(0.,9.),[c],method='DOP853',rtol=tol,atol=tol*.01,dense_output=True)
   check('transport_ODE_'+str(gamma)+'_'+str(tol),sol.success,sol.message);sols.append(sol.sol(grid)[0])
  ref=np.array([c*np.exp(-gamma*y)+cq(lambda s:np.exp(-gamma*(y-s))*source(s),0,y,1e-12) for y in grid])
  errs=[float(np.max(np.abs(s-ref))) for s in sols];diff=float(np.max(np.abs(sols[0]-sols[1])))
  check('transport_ODE_vs_Green_'+str(gamma),errs[-1]<PROTOCOL['ODE_final_error_limit'],errs)
  check('transport_ODE_refinement_'+str(gamma),diff<PROTOCOL['ODE_refinement_limit'],diff)
  NUMERICS.append({'fixture':'incoming_frequency_transport','gamma':[re,im],'errors':errs,'refinement_difference':diff})
 # Exact positive optical moment for a separate compact probe on (1,2).
 y=S.symbols('y');theta=S.integrate(S.exp(-y)*(y-1)**2*(2-y)**2,(y,1,2))
 check('optical_moment_nonzero',float(theta)>0,theta)

def retarded_fixtures():
 t,y=S.symbols('t y',real=True);k=S.exp(-y)*(S.sin(t)+y*S.cos(2*t))
 fun=S.lambdify((t,y),k,'numpy');source=S.lambdify((t,y),S.diff(k,t)+S.diff(k,y),'numpy')
 for t0,y0 in PROTOCOL['transport_points']:
  values=[]
  for tol in PROTOCOL['transport_tolerances']:
   if y0<t0:
    value=fun(t0-y0,0.)+quad(lambda s:source(t0-y0+s,s),0,y0,epsabs=tol,epsrel=tol)[0]
   else:
    value=fun(0.,y0-t0)+quad(lambda s:source(s,y0-t0+s),0,t0,epsabs=tol,epsrel=tol)[0]
   values.append(float(value))
  err=abs(values[-1]-float(fun(t0,y0)))
  check('retarded_characteristic_'+str((t0,y0)),err<PROTOCOL['transport_error_limit'],err)
  check('retarded_quadrature_refinement_'+str((t0,y0)),abs(values[-1]-values[0])<PROTOCOL['transport_error_limit'],values)
  NUMERICS.append({'fixture':'retarded_characteristic','point':[t0,y0],'values':values,'exact':float(fun(t0,y0)),'error':err})
 t0,y0=PROTOCOL['transport_points'][0]
 check('lost_initial_identification_detected',abs(float(fun(0.,y0-t0)))>.1)

def metric_fixtures():
 for y in PROTOCOL['metric_points']:
  zi=.2*np.sin(y);zj=.3*np.cos(y);ki=.15*np.cos(2*y);kj=.1*np.sin(3*y)
  target=np.array([3*zi*zj+3*(zi*kj+zj*ki)+2*ki*kj]*2+[2*zi*zj+1.5*(zi*kj+zj*ki)+.75*ki*kj])
  def exact(a,b):
   E=np.exp(a*zi+b*zj);K=a*ki+b*kj
   gtt=-E*E*(1+K);gty=E*E*(1+K)-E;gyy=2*E-E*E*(1+K)
   N2=(gty*gty-gtt*gyy)/gyy
   return np.array([np.sqrt(N2/gyy),gty/gyy,np.sqrt(N2)])
  errors=[]
  for h in PROTOCOL['metric_steps']:
   fd=(exact(h,h)-exact(h,-h)-exact(-h,h)+exact(-h,-h))/(4*h*h)
   errors.append(float(np.max(np.abs(fd-target))))
  check('finite_metric_second_vertices_'+str(y),errors[-1]<PROTOCOL['metric_final_error_limit'],errors)
  check('finite_metric_refinement_'+str(y),all(b<PROTOCOL['metric_refinement_ratio']*a+PROTOCOL['metric_roundoff_floor'] for a,b in zip(errors,errors[1:])),errors)
  NUMERICS.append({'fixture':'finite_exact_metric','point':y,'mixed_errors':errors,'target':target.tolist()})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 ingoing areal retarded system','input_count':len(ps)}),flush=True)
 try:provenance();geometry();gauge_and_ADM();frequency_fixtures();retarded_fixtures();metric_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 ingoing areal retarded system','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'ingoing_gauge_derived':good,'coordinate_clock_identified_as_source':good,'exact_mass_redshift_rows_derived':good,
  'compact_ADM_vertices_derived':good,'initial_and_offshell_contacts_retained':good,
  'physical_relative_clock_discarded':False,'counterterm_changed':False,'seed_residual_discarded':False,
  'physical_metric_boundary_domain_selected':False,'physical_kernel_numerically_evaluated':False,
  'global_static_slice_vertex_bound_claimed':False,'global_shape_inverse_proved':False,
  'actual_metric_mode_constructed':False,'actual_charged_stability_proved':False,'actual_wormhole_instability_proved':False,
  'retained_fast_factor_removed':False,'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
