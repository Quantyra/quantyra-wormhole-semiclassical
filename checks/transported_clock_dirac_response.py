"""Transported clock response: canonical geometry, state evolution, and contacts."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1];STEM='transported-clock-dirac-response'
OWN=['checks/transported_clock_dirac_response.py','notes/'+STEM+'.md','evidence/transported-clock-source-review.json']
PREVIOUS='boundary-clock-curvature-reconstruction';ROWS=[];NUMERICS=[]
PROTOCOL={'scope':'Coordinate clock/frame response in the transported state and domain, not the physical shape remainder',
 'matrix_seed':1729,'energies':[-3.,-1.,-1.,.4,2.2,4.], 'inverse_temperatures':[.5,2.,5.],
 'time_profile':'256 t^4(1-t)^4 on [0,1], zero outside',
 'ODE_rtol':[1e-7,1e-10],'ODE_atol':[1e-9,1e-12],'ODE_max_step':[.04,.02],
 'coarse_error_tolerance':2e-6,'fine_error_tolerance':2e-9,'refinement_ratio':.1,'refinement_floor':1e-11,
 'frequency_fixtures':[[1.,0.],[.4,.8],[2.,.2]],'frequency_tolerance':2e-12,
 'pullback_amplitudes':[.02,.01,.005],'pullback_refinement_ratio':.35,'pullback_final_tolerance':2e-3,
 'pullback_fixtures':[[.2,.3,.7],[.5,1.,.7],[.35,.6,1.3],[.7,1.2,1.3]],
 'physical_metric_modes_computed':False,'actual_shape_kernel_evaluated':False,
 'instantaneous_rethermalization_adopted':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def zero(name,x):
 x=S.simplify(S.sympify(x).doit());check(name,x==0,x)
def matrixzero(name,x):
 entries=[S.simplify(v) for v in x]
 check(name,all(v==0 for v in entries),entries)
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})
def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==71 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 check('domain_and_state_pinned',all('notes/'+n+'.md' in p['input_sha256_after'] for n in ['coupled-fermion-domain','thermal-shape-response-and-subtractions','neutral-spherical-current-decoupling']))
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('source_scope',len(rev['primary_sources'])==2 and not rev['instantaneous_rethermalization_adopted'] and not rev['curvature_free_column_identified_with_pure_gauge'])

def exact_checks():
 t,y,theta,phi=S.symbols('t y theta phi',real=True);I=S.I
 alpha=S.Matrix([[0,-I],[I,0]]);beta=S.Matrix([[0,1],[1,0]]);ident=S.eye(2)
 a=S.Function('a')(t,y);b=S.Function('b')(t,y);m=S.Function('m')(y)
 psi=S.Matrix([S.Function('u')(y),S.Function('v')(y)])
 H=lambda q:-I*alpha*S.diff(q,y)+beta*m*q
 Qa=lambda q:a*H(q)-I*alpha*S.diff(a,y)*q/2
 Pb=lambda q:-I*(b*S.diff(q,y)+S.diff(b,y)*q/2)
 Gen=lambda q:Qa(q)-Pb(q)
 vertex=S.diff(Gen(psi),t)-I*(Gen(H(psi))-H(Gen(psi)))
 direct=-I*alpha*((S.diff(a,t)-S.diff(b,y))*S.diff(psi,y)+(S.diff(a,t,y)-S.diff(b,y,2))*psi/2)
 direct+=beta*(S.diff(a,t)*m+b*S.diff(m,y))*psi
 direct+=I*((S.diff(b,t)-S.diff(a,y))*S.diff(psi,y)+(S.diff(b,t,y)-S.diff(a,y,2))*psi/2)
 matrixzero('local_clock_vertex_from_differential_commutator',vertex-direct)
 matrixzero('time_generator_anticommutator',Qa(psi)-(a*H(psi)+H(a*psi))/2)
 chi=S.Function('chi')(y);Q=lambda q:chi*H(q)-I*alpha*S.diff(chi,y)*q/2
 matrixzero('mass_independent_time_commutator',Q(H(psi))-H(Q(psi))-S.diff(chi,y)*S.diff(psi,y)-S.diff(chi,y,2)*psi/2)
 pullback=-I*a*H(psi)+b*S.diff(psi,y)+S.diff(b,y)*psi/2-S.diff(a,y)*alpha*psi/2
 matrixzero('current_half_density_generator',pullback+I*Gen(psi))
 check('omitted_half_density_detected',S.simplify((pullback-(-I*a*H(psi)+b*S.diff(psi,y)))[0])!=0)
 tau3=S.diag(1,-1);af=S.kronecker_product(tau3,alpha);bf=S.kronecker_product(S.eye(2),beta)
 projector=(S.eye(4)+S.kronecker_product(tau3,beta))/2
 matrixzero('actual_paired_projector_idempotent',projector*projector-projector)
 matrixzero('actual_paired_projector_hermitian',projector.H-projector)
 matrixzero('actual_paired_boundary_current_zero',projector*af*projector)
 check('actual_projector_maximal_rank',projector.rank()==2)
 matrixzero('scalar_boundary_half_density_preserves_projector',(2*S.eye(4))*projector-projector*(2*S.eye(4)))
 # Coordinate covariance of all four mixed equation components, including off-shell shift contact.
 F=S.Function('F')(y);r=S.Function('r')(y)
 bt=S.Function('bt')(y);br=S.Function('br')(y);bo=S.Function('bo')(y)
 E=S.diag(bt,br,bo,bo);xi=[a,b,0,0];coords=[t,y,theta,phi];Lie=S.zeros(4)
 for i in range(4):
  for j in range(4):
   Lie[i,j]=sum(xi[k]*S.diff(E[i,j],coords[k])-E[k,j]*S.diff(xi[i],coords[k])+E[i,k]*S.diff(xi[k],coords[j]) for k in range(4))
 for i,z in enumerate([bt,br,bo,bo]):zero('mixed_residual_diagonal_'+str(i),Lie[i,i]-b*S.diff(z,y))
 zero('mixed_residual_yt',Lie[1,0]-(br-bt)*S.diff(b,t))
 zero('mixed_residual_ty',Lie[0,1]-(bt-br)*S.diff(a,y))
 zero('off_shell_shift_contact',Lie[1,0]+Lie[0,1]-(br-bt)*(S.diff(b,t)-S.diff(a,y)))
 w=F*r*r;div=S.diff(a,t)+S.diff(b,y)+S.diff(w,y)*b/w
 for name,row in [('trace',bt+br+2*bo),('orbit',bt+br)]:
  zero('fixed_density_'+name,w*(b*S.diff(row,y)+div*row)-S.diff(a*w*row,t)-S.diff(b*w*row,y))
 flux=S.symbols('flux');mag=S.zeros(4);mag[2,3]=flux*S.sin(theta)/2;mag[3,2]=-mag[2,3]
 LieM=S.zeros(4)
 for i in range(4):
  for j in range(4):
   LieM[i,j]=sum(xi[k]*S.diff(mag[i,j],coords[k])+mag[k,j]*S.diff(xi[k],coords[i])+mag[i,k]*S.diff(xi[k],coords[j]) for k in range(4))
 matrixzero('fixed_magnetic_flux_under_clock_transport',LieM)
 g,Delta,fa,fb,Gab=S.symbols('g Delta fa fb Gab')
 V=(g-I*Delta)*Gab;dC=-I*Gab*(fb-fa)
 zero('retarded_denominator_cancellation',(g-I*Delta)*dC+I*V*(fb-fa))
 zero('degenerate_thermal_block',dC.subs(fb,fa))
 hb=S.diag(-2,-1,-1,3);Cb=S.diag(S.Rational(4,5),S.Rational(2,3),S.Rational(2,3),S.Rational(1,5))
 G=S.Matrix([[1,2,I,0],[2,0,1,-I],[-I,1,3,2],[0,I,2,-1]])
 O=S.Matrix([[2,I,1,0],[-I,1,0,2],[1,0,-1,I],[0,2,-I,3]])
 comm=lambda A,B:A*B-B*A
 dC=-I*comm(G,Cb);dO=-I*comm(G,O)
 zero('observable_contact_cancellation',S.trace(dC*O+Cb*dO))
 matrixzero('uniform_clock_does_not_change_prepared_covariance',comm(hb,Cb))
 be,en,eps=S.symbols('be en eps',positive=True)
 ceq=1/(1+S.exp(be*(1+eps)*en))
 zero('instantaneous_KMS_is_different',S.diff(ceq,eps).subs(eps,0)+be*en*S.exp(be*en)/(1+S.exp(be*en))**2)
 check('instantaneous_KMS_extra_term_nonzero',S.diff(ceq,eps).subs({eps:0,be:1,en:1})!=0)
 # Initial generator contact cancels the coordinate covariance at t=0 exactly.
 D0=S.eye(4)/7
 matrixzero('initial_generator_contact_retained',-I*comm(G,Cb)+(D0+I*comm(G,Cb))-D0)

def pulse_functions():
 t=S.symbols('t');p=256*t**4*(1-t)**4
 fs=[S.lambdify(t,S.diff(p,t,n),'numpy') for n in range(3)]
 def pulse(n,t):return fs[n](t) if 0<t<1 else 0.
 return pulse

def matrix_fixtures():
 rng=np.random.default_rng(PROTOCOL['matrix_seed']);energies=np.array(PROTOCOL['energies']);h=np.diag(energies);n=len(energies)
 def herm(scale):
  z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));return scale*(z+z.conj().T)/2
 Q0=herm(.3);Q1=herm(.08);O=herm(.4);initial=herm(.015);pulse=pulse_functions()
 comm=lambda a,b:a@b-b@a;times=np.array([0.,.15,.3,.5,.72,.9,1.,1.2])
 specs=[(beta,kind) for beta in PROTOCOL['inverse_temperatures'] for kind in ['coordinate','independent_state']]+[(2.,'initial_generator')]
 for idx,(beta,kind) in enumerate(specs):
  cov=np.diag(1/(1+np.exp(beta*energies)));offset=.2 if kind=='initial_generator' else 0.
  d0=initial if kind=='independent_state' else np.zeros((n,n),complex)
  gen=lambda t:(offset+pulse(0,t))*Q0+pulse(1,t)*Q1
  gd=lambda t:pulse(1,t)*Q0+pulse(2,t)*Q1
  vertex=lambda t:gd(t)-1j*comm(gen(t),h)
  D=d0+1j*comm(gen(0),cov)
  def expected(t):
   evolution=np.exp(-1j*energies*t)
   return -1j*comm(gen(t),cov)+evolution[:,None]*D*evolution.conj()[None,:]
  errors=[]
  for j,(rtol,atol,step) in enumerate(zip(PROTOCOL['ODE_rtol'],PROTOCOL['ODE_atol'],PROTOCOL['ODE_max_step'])):
   def rhs(t,flat):return (-1j*(comm(h,flat.reshape(n,n))+comm(vertex(t),cov))).ravel()
   sol=solve_ivp(rhs,(0.,1.2),d0.ravel(),method='DOP853',rtol=rtol,atol=atol,max_step=step,dense_output=True)
   check('retarded_matrix_ODE_'+str(idx)+'_'+str(j),sol.success)
   err=max(np.linalg.norm(sol.sol(t).reshape(n,n)-expected(t)) for t in times);errors.append(float(err))
  check('retarded_matrix_coarse_'+str(idx),errors[0]<PROTOCOL['coarse_error_tolerance'],errors[0])
  check('retarded_matrix_fine_'+str(idx),errors[1]<PROTOCOL['fine_error_tolerance'],errors[1])
  check('retarded_matrix_refinement_'+str(idx),errors[1]<PROTOCOL['refinement_ratio']*errors[0]+PROTOCOL['refinement_floor'],errors)
  if kind=='coordinate':
   contact=max(abs(np.trace((-1j*comm(gen(t),cov))@O+cov@(-1j*comm(gen(t),O)))) for t in times)
   check('transported_observable_contact_'+str(idx),contact<1e-12,contact)
   check('coordinate_pulse_has_no_final_state_change_'+str(idx),np.linalg.norm(expected(1.2))<1e-12)
  elif kind=='independent_state':check('independent_initial_state_is_retained_'+str(idx),abs(np.linalg.norm(expected(1.2))-np.linalg.norm(initial))<1e-12)
  else:check('nonzero_initial_generator_contact_needed',np.linalg.norm(1j*comm(gen(0),cov))>1e-3)
  NUMERICS.append({'fixture':'retarded_covariance_ODE','beta':beta,'kind':kind,'errors':errors})
  if idx==0:
   wrong=solve_ivp(lambda t,q:(-1j*(comm(h,q.reshape(n,n))+comm(gd(t),cov))).ravel(),(0.,1.2),d0.ravel(),method='DOP853',rtol=1e-10,atol=1e-12,max_step=.02,dense_output=True)
   discrepancy=max(np.linalg.norm(wrong.sol(t).reshape(n,n)-expected(t)) for t in times)
   check('omitted_shift_commutator_detected',discrepancy>1e-3,discrepancy)
   check('omitted_observable_contact_detected',max(abs(np.trace((-1j*comm(gen(t),cov))@O)) for t in times)>1e-3)
 cov=np.diag(1/(1+np.exp(2*energies)));diff=energies[:,None]-energies[None,:]
 for idx,(gr,gi) in enumerate(PROTOCOL['frequency_fixtures']):
  gamma=complex(gr,gi);G=Q0+gamma*Q1;V=gamma*G-1j*comm(G,h)
  computed=-1j*comm(V,cov)/(gamma+1j*diff);target=-1j*comm(G,cov)
  error=float(np.linalg.norm(computed-target));check('complex_frequency_Ward_'+str(idx),error<PROTOCOL['frequency_tolerance'],error)
  NUMERICS.append({'fixture':'frequency_covariance','gamma':[gr,gi],'error':error})

def pullback_fixtures():
 alpha=np.array([[0,-1j],[1j,0]],complex);beta=np.array([[0,1],[1,0]],complex);pulse=pulse_functions()
 for idx,(t,y,mass) in enumerate(PROTOCOL['pullback_fixtures']):
  ks=[.8,-1.1];amplitudes=[.7+.2j,-.3+.4j];modes=[]
  for j,k in enumerate(ks):
   energies,vectors=np.linalg.eigh(alpha*k+beta*mass);which=1 if j==0 else 0
   modes.append((k,energies[which],vectors[:,which]*amplitudes[j]))
  def wave(t,y):return sum(v*np.exp(1j*k*y-1j*E*t) for k,E,v in modes)
  psi=wave(t,y);py=sum(1j*k*v*np.exp(1j*k*y-1j*E*t) for k,E,v in modes)
  hpsi=sum(E*v*np.exp(1j*k*y-1j*E*t) for k,E,v in modes)
  chi=np.exp(-y*y);chiy=-2*y*chi
  a=pulse(0,t)*chi;ay=pulse(0,t)*chiy
  b=pulse(1,t)*y*chi;by=pulse(1,t)*(chi+y*chiy)
  exact=-1j*a*hpsi+b*py+by*psi/2-ay*(alpha@psi)/2
  def transported(eps):
   W=(1+eps*by)*np.eye(2)-eps*ay*alpha
   eigen,vec=np.linalg.eigh(W)
   if min(eigen)<=0:raise AssertionError('nonspacelike pullback fixture')
   root=(vec*np.sqrt(eigen))@vec.conj().T
   return root@wave(t+eps*a,y+eps*b)
  errors=[]
  for eps in PROTOCOL['pullback_amplitudes']:
   numeric=(transported(eps)-transported(-eps))/(2*eps)
   errors.append(float(np.linalg.norm(numeric-exact)/(1+np.linalg.norm(exact))))
  check('independent_Dirac_pullback_refinement_'+str(idx),all(b<PROTOCOL['pullback_refinement_ratio']*a for a,b in zip(errors,errors[1:])),errors)
  check('independent_Dirac_pullback_limit_'+str(idx),errors[-1]<PROTOCOL['pullback_final_tolerance'],errors[-1])
  NUMERICS.append({'fixture':'exact_wave_coordinate_pullback','t':t,'y':y,'mass':mass,'errors':errors})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 transported clock Dirac response','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();matrix_fixtures();pullback_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 transported clock Dirac response','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'canonical_coordinate_generator_derived':good,'transported_prepared_state_response_derived':good,
  'coordinate_spectral_cancellation_derived':good,'observable_contacts_retained':good,
  'off_shell_density_contacts_derived':good,'actual_projector_flux_retained':good,'initial_state_contacts_retained':good,
  'instantaneous_rethermalization_adopted':False,'counterterm_changed':False,'seed_residual_discarded':False,
  'curvature_free_column_identified_with_pure_gauge':False,'genuine_boundary_shape_response_evaluated':False,
  'physical_clock_closure_solved':False,'global_shape_inverse_proved':False,'physical_coupled_moment_evaluated':False,
  'physical_metric_boundary_selected':False,'actual_metric_causal_inverse_proved':False,'actual_metric_mode_constructed':False,
  'actual_charged_stability_proved':False,'actual_wormhole_instability_proved':False,'retained_fast_factor_removed':False,
  'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
