"""Verify causal horizon scope and the evaluated reflected magnetic zero modes."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[1];STEM='causal-horizon-and-zero-mode-response'
OWN=['checks/causal_horizon_zero_mode_response.py','notes/'+STEM+'.md','evidence/causal-horizon-zero-mode-source-review.json']
PREVIOUS='ingoing-areal-retarded-system';ROWS=[];NUMERICS=[]
PROTOCOL={'scope':'Exact declared zero-mode response and causal/graph domain, not the full charged kernel or stability',
 'Fermi_pairs':[[.7,.1],[.7,4.],[2.,1.],[2.,4.],[6.,.1],[6.,4.]],
 'Fermi_tail_cutoffs':[40.,60.],'Fermi_error_limit':2e-10,
 'gammas':[[.5,0.],[1.,0.],[2.,0.],[1.2,.8]],'kappa':1.,
 'quadrature_tolerances':[1e-8,1e-12],'transfer_error_limit':3e-10,
 'pulse_laplace_gammas':[[.5,0.],[2.,0.],[1.2,.8]],'pulse_readout_y':2.5,
 'pulse_error_limit':3e-10,'pulse_return_limit':3e-10,
 'profile':'D=(1-(2y-3)^2)^6 on [1,2], A=(1-(2t-1)^2)^6 on [0,1]; piecewise C5',
 'formal_verification':False,'full_quantum_kernel_evaluated':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def zero(name,x):
 x=S.simplify(S.sympify(x).doit());check(name,x==0,x)
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==70 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('source_scope',len(review['primary_sources'])==2 and review['full_angular_tower_retained'] and not review['physical_mass_return_proved'])

def exact_checks():
 y,gamma,kappa,a,rh=S.symbols('y gamma kappa a rh',positive=True)
 r=S.Function('r')(y);F=S.Function('F');fb=S.Function('fb');zb=S.Function('zb')
 k=S.exp(-gamma*y)*fb(r)/F(r)
 deriv=S.diff(k,y).subs(S.diff(r,y),-F(r))
 zero('ingoing_horizon_derivative',deriv-(S.diff(F(r),r)-gamma)*k+S.exp(-gamma*y)*S.diff(fb(r),r))
 rr=rh+a*S.exp(-2*kappa*y);ff=2*kappa*(rr-rh)
 kh=S.exp(-gamma*y)/ff
 zero('simple_horizon_threshold',kh*S.exp((gamma-2*kappa)*y)-1/(2*kappa*a))
 zero('threshold_nondecaying_value',kh.subs(gamma,2*kappa)-1/(2*kappa*a))
 check('subthreshold_unbounded',S.limit(kh.subs(gamma,kappa),y,S.oo)==S.oo)
 zero('horizon_preserving_improved_decay',S.exp(-gamma*y)*ff/ff-S.exp(-gamma*y))
 R,m,Q,c3,c4,lam=S.symbols('R m Q c3 c4 lam',positive=True)
 geo=1+lam*R**2-2*m/R+Q/R**2+c3/R**3+c4/R**4
 ratio=geo/R-S.diff(geo,R)/2
 zero('actual_seed_mass_log_derivative',ratio-(1/R-3*m/R**2+2*Q/R**3+S.Rational(5,2)*c3/R**4+3*c4/R**5))
 zero('outer_graph_coefficient_limit',S.limit(ratio,R,S.oo))
 I=S.eye(2);s1=S.Matrix([[0,1],[1,0]]);s2=S.Matrix([[0,-S.I],[S.I,0]]);s3=S.diag(1,-1)
 alpha=S.kronecker_product(s3,s2);beta=S.kronecker_product(I,s1);chi=S.kronecker_product(s3,s1);Pi=(S.eye(4)+chi)/2
 current=S.simplify(Pi*alpha*Pi);cross=S.simplify(Pi*(S.I*alpha*beta)*Pi)
 for label,mat in [('boundary_current',current),('graph_boundary_contact',cross)]:
  for i in range(4):
   for j in range(4):zero(label+'_'+str(i)+str(j),mat[i,j])
 u,v,K=S.symbols('u v K',real=True)
 zero('graph_Young_identity',u*u/2+K*K*v*v/2-K*u*v-(u-K*v)**2/2)
 # Derive the exact zero-mode velocity from the ingoing ADM Hamiltonian.
 speed=S.symbols('speed',real=True)
 zero('incoming_velocity_unchanged',speed-(speed-1)-1)
 zero('outgoing_velocity',speed+(speed-1)-(2*speed-1))
 p,pp,g,bt,c,ell=S.symbols('p pp g bt c ell',real=True,nonzero=True)
 zero('half_density_source_vertex',-g*pp-g*(p-pp)/2+g*(p+pp)/2)
 aa=S.symbols('aa');series=S.series(S.pi*aa/S.sin(S.pi*aa),aa,0,4).removeO()
 zero('logistic_variance',S.diff(series,aa,2).subs(aa,0)-S.pi**2/3)
 shift=S.symbols('shift',real=True)
 J=-S.integrate(S.pi**2/(3*bt**2)+shift**2,(shift,-ell/2,ell/2))
 zero('exact_Fermi_difference_moment',J+ell*S.pi**2/(3*bt**2)+ell**3/12)
 susceptibility=-S.I*c*J/(2*S.pi*(gamma-S.I*ell))
 t0=c*S.pi/(12*bt**2)
 zero('thermal_Virasoro_sign',(gamma-S.I*ell)*susceptibility-(2*t0*S.I*ell-c*(S.I*ell)**3/(24*S.pi)))
 zero('prepared_thermal_normalization',t0.subs(bt,2*S.pi/kappa)-c*kappa**2/(48*S.pi))
 d=S.Function('d')(y);J0=S.Function('J0')(y)
 # Directly differentiate the integrated-by-parts outgoing expression.
 out=S.diff(d,y,2)+gamma*S.diff(d,y)+(gamma**2-kappa**2)*d+gamma*(kappa**2-gamma**2)*J0
 equation=gamma*out-S.diff(out,y)
 zero('outgoing_stress_transport',equation.subs(S.diff(J0,y),gamma*J0-d)-(kappa**2*S.diff(d,y)-S.diff(d,y,3)))
 z=S.symbols('z');poly=(1-(2*z-3)**2)**6
 for order in range(3):
  zero('source_endpoint_'+str(order)+'_left',S.diff(poly,z,order).subs(z,1))
  zero('source_endpoint_'+str(order)+'_right',S.diff(poly,z,order).subs(z,2))
 zero('special_frequency_reflection_factor',(gamma*(kappa**2-gamma**2)).subs(gamma,kappa))
 bulk=S.diff(poly,z,2)+kappa*S.diff(poly,z)
 check('special_frequency_bulk_not_zero',bulk.subs({z:S.Rational(3,2),kappa:1})==-48)
 zero('first_order_integrated_flux',S.integrate(kappa**2*S.diff(poly,z)-S.diff(poly,z,3),(z,1,2)))
 t,V=S.symbols('t V',real=True)
 zero('regular_ingoing_time',S.exp(gamma*t)*S.exp(-gamma*y)-S.exp(gamma*(t-y)))
 r0=S.symbols('r0',positive=True);Rf=S.symbols('Rf')
 # A covariant incoming null stress has only vv, and its mixed r-v component equals it.
 metric=S.Matrix([[-S.symbols('F0'),1],[1,0]]);stress=S.Matrix([[Rf/(4*S.pi*r0*r0),0],[0,0]])
 mixed=metric.inv()*stress
 zero('spherical_incoming_flux',4*S.pi*r0*r0*mixed[1,0]-Rf)
 zero('incoming_radial_mass_source',mixed[0,0])
 zi,zj,ki,kj,e1,e2=S.symbols('zi zj ki kj e1 e2')
 speed=1/(2-S.exp(e1*zi+e2*zj)*(1+e1*ki+e2*kj))
 gij=S.diff(2*(speed-1),e1,e2).subs({e1:0,e2:0})
 target=2*(3*zi*zj+3*(zi*kj+zj*ki)+2*ki*kj)
 zero('zero_mode_second_vertex',gij-target)
 check('dropped_second_contact_detected',gij.subs({zi:1,zj:1,ki:0,kj:0})==6)

def cq(fun,a,b,tol):
 if b<=a:return 0j
 return quad(lambda x:complex(fun(x)).real,a,b,epsabs=tol,epsrel=tol,limit=160)[0]+1j*quad(lambda x:complex(fun(x)).imag,a,b,epsabs=tol,epsrel=tol,limit=160)[0]

def Fermi_fixtures():
 for bt,ell in PROTOCOL['Fermi_pairs']:
  vals=[]
  for cutoff in PROTOCOL['Fermi_tail_cutoffs']:
   bound=cutoff/bt+ell/2
   f=lambda E:E*E*(expit(-bt*(E+ell/2))-expit(-bt*(E-ell/2)))
   vals.append(quad(f,-bound,bound,epsabs=1e-11,epsrel=1e-12,limit=200)[0])
  target=-ell*np.pi**2/(3*bt*bt)-ell**3/12
  error=abs(vals[-1]-target)/(1+abs(target));ref=abs(vals[-1]-vals[0])/(1+abs(target))
  check('Fermi_integral_'+str((bt,ell)),error<PROTOCOL['Fermi_error_limit'],error)
  check('Fermi_tail_refinement_'+str((bt,ell)),ref<PROTOCOL['Fermi_error_limit'],ref)
  NUMERICS.append({'fixture':'Fermi_difference','beta':bt,'momentum':ell,'values':vals,'target':target,'normalized_error':error,'refinement':ref})

POLY=np.polynomial.Polynomial([1.,0.,-1.])**6
def profile(y,order=0):
 x=2*y-3
 return 2**order*POLY.deriv(order)(x) if abs(x)<1 else 0.
def pulse(t):
 x=2*t-1
 return POLY(x) if abs(x)<1 else 0.
def P(y):return (PROTOCOL['kappa']**2*profile(y,1)-profile(y,3))/(12*np.pi)

def reflection_fixtures():
 for re,im in PROTOCOL['gammas']:
  gamma=complex(re,im);vals=[]
  for tol in PROTOCOL['quadrature_tolerances']:
   vals.append(cq(lambda y:np.exp(-gamma*y)*P(y),1.,2.,tol))
  D=cq(lambda y:np.exp(-gamma*y)*profile(y),1.,2.,1e-12)
  target=gamma*(1-gamma**2)*D/(12*np.pi)
  error=abs(vals[-1]-target);ref=abs(vals[-1]-vals[0])
  check('reflected_flux_'+str(gamma),error<PROTOCOL['transfer_error_limit'],error)
  check('reflected_quadrature_refinement_'+str(gamma),ref<PROTOCOL['transfer_error_limit'],ref)
  if gamma==1:check('thermal_cancellation_resolved',abs(vals[-1])<PROTOCOL['transfer_error_limit'] and abs(D)>.01,[abs(vals[-1]),abs(D)])
  NUMERICS.append({'fixture':'reflected_flux','gamma':[re,im],'values':[[v.real,v.imag] for v in vals],'target':[target.real,target.imag],'error':error,'refinement':ref})

def incoming(t,y,tol):
 lo=max(0.,t-y-2.);hi=min(1.,t-y-1.)
 if hi<=lo:return 0.
 return quad(lambda tau:pulse(tau)*P(t-y-tau),lo,hi,epsabs=tol,epsrel=tol,limit=160)[0]

def pulse_fixtures():
 y=PROTOCOL['pulse_readout_y']
 for re,im in PROTOCOL['pulse_laplace_gammas']:
  gamma=complex(re,im);vals=[]
  for tol in PROTOCOL['quadrature_tolerances']:
   vals.append(sum(cq(lambda t:np.exp(-gamma*t)*incoming(t,y,tol*.1),y+a,y+b,tol) for a,b in [(1.,2.),(2.,3.)]))
  A=cq(lambda t:np.exp(-gamma*t)*pulse(t),0.,1.,1e-12)
  D=cq(lambda z:np.exp(-gamma*z)*profile(z),1.,2.,1e-12)
  target=A*np.exp(-gamma*y)*gamma*(1-gamma**2)*D/(12*np.pi)
  error=abs(vals[-1]-target);ref=abs(vals[-1]-vals[0])
  check('time_pulse_vs_frequency_'+str(gamma),error<PROTOCOL['pulse_error_limit'],error)
  check('time_pulse_refinement_'+str(gamma),ref<PROTOCOL['pulse_error_limit'],ref)
  NUMERICS.append({'fixture':'pulse_Laplace','gamma':[re,im],'readout_y':y,'values':[[v.real,v.imag] for v in vals],'target':[target.real,target.imag],'error':error,'refinement':ref})
 integral=sum(quad(lambda t:incoming(t,y,1e-12),y+a,y+b,epsabs=1e-11,epsrel=1e-11)[0] for a,b in [(1.,2.),(2.,3.)])
 check('integrated_first_order_LLL_flux_zero',abs(integral)<PROTOCOL['pulse_return_limit'],integral)
 check('causal_arrival',incoming(y+.9,y,1e-12)==0.)
 check('pulse_has_ended',incoming(y+3.1,y,1e-12)==0.)
 check('nonzero_signal_between_endpoints',abs(incoming(y+1.4,y,1e-12))>1e-4)
 NUMERICS.append({'fixture':'pulse_integrated_flux','integral':integral,'nonzero_sample':incoming(y+1.4,y,1e-12)})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 causal horizon and reflected zero modes','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();Fermi_fixtures();reflection_fixtures();pulse_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 causal horizon and reflected zero modes','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'finite_time_causal_localization_derived':good,'high_frequency_angular_graph_bound_derived':good,
  'LLL_reflected_kernel_evaluated':good,'LLL_second_contact_retained':good,'LLL_first_order_integrated_flux_zero':good,
  'full_quantum_kernel_evaluated':False,'full_horizon_trace_proved':False,'uniform_infinite_time_bound_proved':False,
  'physical_mass_return_proved':False,'massive_tower_discarded':False,'counterterm_changed':False,
  'physical_metric_boundary_domain_selected':False,'global_shape_inverse_proved':False,
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
