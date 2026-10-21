"""Check ADM geometry, boundary flux and fixed-domain retarded readout."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1];STEM='adm-dirac-boundary-response'
OWN=['checks/adm_dirac_boundary_response.py','notes/'+STEM+'.md','evidence/adm-dirac-source-review.json']
PREVIOUS='transported-clock-dirac-response';ROWS=[];NUMERICS=[]
PROTOCOL={'scope':'Exact vertices and domain contacts; physical kernel not evaluated',
 'seed':20261021,'energies':[-3.,-1.,-1.,.4,2.2,4.],
 'betas':[.5,2.,5.],'gammas':[[1.,0.],[.4,.8],[2.,.2]],
 'frequency_tolerance':3e-12,'ODE_rtol':[1e-7,1e-10],'ODE_atol':[1e-9,1e-12],
 'ODE_max_step':[.04,.02],'fine_tolerance':3e-9,'refinement_ratio':.15,'refinement_floor':2e-11,
 'matrix_steps':[.04,.02,.01],'matrix_times':[.2,.5,.8],
 'matrix_refinement_ratio':.35,'matrix_refinement_floor':2e-10,'matrix_tolerance':2e-4,
 'physical_kernel_evaluated':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,x):
 x=S.simplify(S.sympify(x).doit());check(name,x==0,x)

def mz(name,x):
 vals=[S.simplify(v) for v in x];check(name,all(v==0 for v in vals),vals)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==90 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('source_scope',len(rev['primary_sources'])==3 and not rev['new_boundary_coupling_selected'] and not rev['genuine_boundary_shape_response_evaluated'])

def geometry():
 t,y=S.symbols('t y',real=True);I=S.I
 N,L,r=[S.Function(z)(t,y) for z in ['N','L','r']]
 shift=S.Function('shift')(t,y);lam=S.symbols('lambda',real=True)
 alpha=S.Matrix([[0,-I],[I,0]]);beta=S.Matrix([[0,1],[1,0]]);ident=S.eye(2)
 g0=beta*alpha;g1=beta
 mz('gamma_time',g0*g0+ident);mz('gamma_space',g1*g1-ident);mz('gamma_anticommutator',g0*g1+g1*g0)
 aa=S.diff(N,y)/(N*L);kk=(S.diff(L,t)-S.diff(L*shift,y))/(N*L)
 # Cartan two-form coefficients in dt wedge dy, evaluated directly.
 omega_t=aa*N+kk*L*shift;omega_y=kk*L
 zero('cartan_time_coframe',-S.diff(N,y)+omega_t*L-omega_y*L*shift)
 zero('cartan_space_coframe',S.diff(L,t)-S.diff(L*shift,y)-omega_y*N)
 phi=S.Matrix([S.Function('f')(t,y),S.Function('g')(t,y)])
 chi=phi/S.sqrt(L);v=N/L;mass=lam*N/r
 # Solve the covariant first-order equation for chi_t using its spin connection.
 angular=-I*alpha*lam/r
 chit=shift*S.diff(chi,y)+g0*g1*(v*S.diff(chi,y)+N*aa*chi/2)-N*kk*chi/2+N*g0*angular*chi
 covariant=g0*(chit-shift*S.diff(chi,y))/N+g1*S.diff(chi,y)/L+kk*g0*chi/2+aa*g1*chi/2+angular*chi
 mz('covariant_equation_solved_for_time',covariant)
 phit=S.sqrt(L)*chit+S.diff(L,t)*phi/(2*L)
 target=-alpha*(v*S.diff(phi,y)+S.diff(v,y)*phi/2)+shift*S.diff(phi,y)+S.diff(shift,y)*phi/2-I*beta*mass*phi
 mz('covariant_to_canonical_ADM_evolution',phit-target)
 # Direct action of the radius connection cancels under chi=r Psi.
 f=S.Function('q')(t,y)
 for z in [t,y]:zero('radius_connection_cancellation_'+str(z),r*(S.diff(f/r,z)+S.diff(r,z)*f/r**2)-S.diff(f,z))
 eps,eta=S.symbols('eps eta');di,dj,wi,wj,si,sj=S.symbols('di dj wi wj si sj',real=True)
 speed=S.exp(eps*di+eta*dj);m=S.exp(eps*wi+eta*wj);sh=eps*si+eta*sj
 zero('kinetic_first_vertex',S.diff(speed,eps).subs({eps:0,eta:0})-di)
 zero('kinetic_mixed_contact',S.diff(speed,eps,eta).subs({eps:0,eta:0})-di*dj)
 zero('mass_mixed_contact',S.diff(m,eps,eta).subs({eps:0,eta:0})-wi*wj)
 zero('affine_shift_second_contact',S.diff(sh,eps,eta))
 ca,cb,cp,z=S.symbols('ca cb cp z',real=True)
 zero('conformal_speed_invariance',S.exp((ca+z)-(cb+z))-S.exp(ca-cb))
 zero('conformal_mass_invariance',S.exp((ca+z)-(cp+z))-S.exp(ca-cp))
 # Direct continuity equation, without assuming a constant lapse or shift.
 psi=S.Matrix([S.Function('u')(y),S.Function('w')(y)])
 K=v*alpha-shift*ident;Q=-I*S.diff(K,y)/2+beta*mass
 # Formal adjoint zeroth coefficient must be Q + i K_y = Q^dagger.
 mz('canonical_Green_identity',Q+I*S.diff(K,y)-Q.xreplace({I:-I}).T)

def boundary():
 I=S.I;alpha=S.Matrix([[0,-I],[I,0]]);beta=S.Matrix([[0,1],[1,0]])
 A=S.kronecker_product(S.diag(1,-1),alpha);chi=S.kronecker_product(S.diag(1,-1),beta)
 ident=S.eye(4);P=(ident+chi)/2;J=A*chi
 mz('paired_projector',P*P-P);mz('paired_zero_flux',P*A*P)
 mz('rotation_generator_antihermitian',J.H+J);mz('rotation_generator_square',J*J+ident)
 th=S.symbols('theta',real=True);u=S.sin(th);c=S.cos(th)
 # The physical branch is -pi/2 < theta < pi/2, hence c>0.
 Pu=(ident+c*chi+u*A)/2
 mz('moving_projector_idempotent',Pu*Pu-Pu);mz('moving_projector_hermitian',Pu.H-Pu)
 zero('moving_projector_half_rank',S.trace(Pu)-2)
 mz('moving_projector_current',Pu*A*Pu-u*Pu)
 R=S.cos(th/2)*ident+S.sin(th/2)*J
 mz('unitary_trivialization',R.H*R-ident);mz('rotation_transports_projector',R*P*R.H-Pu)
 vv=S.symbols('speed',positive=True)
 mz('fixed_boundary_principal',R.H*(vv*A-vv*u*ident)*R-vv*c*A-vv*u*(chi-ident))
 mz('old_projector_flux_defect',P*(vv*A-vv*u*ident)*P+vv*u*P)
 eta=S.symbols('rapidity',real=True);D=S.cosh(eta/2)*ident+S.sinh(eta/2)*A
 mz('rest_frame_boost_flux',D*(S.cosh(eta)*A-S.sinh(eta)*ident)*D-A)
 # Verify range equality using u=tanh eta without square-root branch ambiguities.
 Peta=(ident+chi/S.cosh(eta)+S.tanh(eta)*A)/2
 mz('boosted_range_is_projector_range',(ident-Peta)*D*P)
 e,f=S.symbols('e f');di,dj,si,sj=S.symbols('di dj si sj',real=True)
 theta=S.asin((e*si+f*sj)*S.exp(-e*di-f*dj))
 zero('angle_first',S.diff(theta,e).subs({e:0,f:0})-si)
 zero('angle_mixed',S.diff(theta,e,f).subs({e:0,f:0})+si*dj+sj*di)
 comm=lambda x,z:x*z-z*x
 Xi=si*J/2;Xj=sj*J/2;Xij=-(si*dj+sj*di)*J/2
 Ki=di*A-si*ident;Kj=dj*A-sj*ident
 first=Ki+comm(A,Xi)
 second=di*dj*A+comm(Ki,Xj)+comm(Kj,Xi)+(comm(comm(A,Xi),Xj)+comm(comm(A,Xj),Xi))/2+comm(A,Xij)
 mz('first_boundary_contact',first-di*A-si*(chi-ident))
 mz('mixed_boundary_contact',second-(di*dj-si*sj)*A)
 mz('all_fixed_domain_boundary_vertices_flux_free',P*first*P+P*second*P)
 wrong=second-comm(A,Xij)
 check('missing_second_projector_contact_detected',any(S.simplify(z)!=0 for z in P*wrong*P))
 check('missing_first_projector_contact_detected',any(S.simplify(z)!=0 for z in P*Ki*P))

def fixtures():
 rng=np.random.default_rng(PROTOCOL['seed']);energies=np.array(PROTOCOL['energies']);h=np.diag(energies);n=len(energies)
 comm=lambda a,b:a@b-b@a
 def herm(scale):
  z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));return scale*(z+z.conj().T)/2
 B=[herm(.15),herm(.17)];W=[[herm(.07),herm(.05)],[None,herm(.08)]];W[1][0]=W[0][1]
 J=1j*herm(.2);X=[.7*J,-.4*J];X2=[[.2*J,.3*J],[.3*J,-.1*J]]
 A=[B[i]+comm(h,X[i]) for i in range(2)];D=[-1j*x for x in X]
 T=[[W[i][j]+comm(B[i],X[j])+comm(B[j],X[i])+(comm(comm(h,X[i]),X[j])+comm(comm(h,X[j]),X[i]))/2+comm(h,X2[i][j]) for j in range(2)] for i in range(2)]
 diff=energies[:,None]-energies[None,:]
 maxwrong=0.;maxinitial=0.
 for bi,beta in enumerate(PROTOCOL['betas']):
  C=np.diag(1/(1+np.exp(beta*energies)))
  for gi,(gr,gm) in enumerate(PROTOCOL['gammas']):
   gamma=complex(gr,gm);den=gamma+1j*diff;errors=[]
   for i in range(2):
    for j in range(2):
     dc=-1j*comm(A[j]+gamma*D[j],C)/den;db=-1j*comm(B[j],C)/den
     target=np.trace(C@W[i][j])+np.trace(B[i]@db)
     computed=np.trace(C@T[i][j])+np.trace((A[i]-gamma*D[i])@dc)
     errors.append(abs(computed-target))
     wrong=np.trace(C@T[i][j])+np.trace((A[i]+gamma*D[i])@dc)
     maxwrong=max(maxwrong,abs(wrong-target))
   err=float(max(errors));check('retarded_frequency_invariance_'+str(bi)+'_'+str(gi),err<PROTOCOL['frequency_tolerance'],err)
   NUMERICS.append({'fixture':'full_matrix_retarded_kernel','beta':beta,'gamma':[gr,gm],'error':err})
  # Independent time evolution in the unrotated and fixed representations.
  times=np.array([.12,.25,.4,.55,.7,.85,1.,1.2]);j=1;errors=[];i=0
  def pulse(t):return 256*t**4*(1-t)**4 if 0<t<1 else 0.
  def pd(t):return 1024*t**3*(1-t)**3*(1-2*t) if 0<t<1 else 0.
  for rtol,atol,step in zip(PROTOCOL['ODE_rtol'],PROTOCOL['ODE_atol'],PROTOCOL['ODE_max_step']):
   def rhs(t,z):
    cs=z.reshape(2,n,n);vs=[B[j]*pulse(t),A[j]*pulse(t)+D[j]*pd(t)]
    return np.array([-1j*(comm(h,cs[k])+comm(vs[k],C)) for k in range(2)]).ravel()
   sol=solve_ivp(rhs,(0.,1.2),np.zeros((2,n,n),complex).ravel(),method='DOP853',rtol=rtol,atol=atol,max_step=step,dense_output=True)
   check('matrix_ODE_completed_'+str(bi)+'_'+str(len(errors)),sol.success)
   err=0.
   for t in times:
    db,dc=sol.sol(t).reshape(2,n,n);dct=rhs(t,sol.sol(t)).reshape(2,n,n)[1]
    err=max(err,np.linalg.norm(dc-db+pulse(t)*comm(X[j],C)))
    force=np.trace(A[i]@dc-D[i]@dct)+pulse(t)*np.trace(C@T[i][j])
    bare=np.trace(B[i]@db)+pulse(t)*np.trace(C@W[i][j]);err=max(err,abs(force-bare))
   errors.append(float(err))
  check('matrix_ODE_fine_'+str(bi),errors[-1]<PROTOCOL['fine_tolerance'],errors)
  check('matrix_ODE_refinement_'+str(bi),errors[-1]<PROTOCOL['refinement_ratio']*errors[0]+PROTOCOL['refinement_floor'],errors)
  NUMERICS.append({'fixture':'retarded_time_evolution_and_force','beta':beta,'errors':errors})
  c0=herm(.02);gamma=.7+.4j;dc0=c0/(gamma+1j*diff)
  # Direct evolution derivative is -i[h,delta C]; Laplace includes c0.
  force_direct=np.trace(A[i]@dc0+1j*D[i]@comm(h,dc0))
  force_laplace=np.trace((A[i]-gamma*D[i])@dc0+D[i]@c0)
  check('initial_derivative_contact_'+str(bi),abs(force_direct-force_laplace)<1e-13)
  maxinitial=max(maxinitial,abs(np.trace(D[i]@c0)))
 check('wrong_readout_sign_detected',maxwrong>1e-3,maxwrong)
 check('missing_initial_contact_detected',maxinitial>1e-4,maxinitial)
 # Independent exact exponential conjugation, including source velocities.
 timewrong=0.;secondwrong=0.
 for idx,t in enumerate(PROTOCOL['matrix_times']):
  l1=np.sin(t);l2=np.cos(1.3*t);d1=np.cos(t);d2=-1.3*np.sin(1.3*t)
  expected1=l1*A[0]+d1*D[0]
  expected2=l1*l2*T[0][1]-1j*(d1*l2+l1*d2)*X2[0][1]
  def full(e,f):
   x=e*l1*X[0]+f*l2*X[1]+e*f*l1*l2*X2[0][1]
   xt=e*d1*X[0]+f*d2*X[1]+e*f*(d1*l2+l1*d2)*X2[0][1]
   bare=h+e*l1*B[0]+f*l2*B[1]+e*f*l1*l2*W[0][1]
   rot=expm(x);return rot.conj().T@bare@rot-1j*xt
  err1=[];err2=[]
  for step in PROTOCOL['matrix_steps']:
   fd1=(full(step,0)-full(-step,0))/(2*step)
   fd2=(full(step,step)-full(step,-step)-full(-step,step)+full(-step,-step))/(4*step**2)
   err1.append(float(np.linalg.norm(fd1-expected1)))
   err2.append(float(np.linalg.norm(fd2-expected2)))
  for name,errs in [('first',err1),('mixed',err2)]:
   check('exact_conjugation_'+name+'_refinement_'+str(idx),all(b<PROTOCOL['matrix_refinement_ratio']*a+PROTOCOL['matrix_refinement_floor'] for a,b in zip(errs,errs[1:])),errs)
   check('exact_conjugation_'+name+'_limit_'+str(idx),errs[-1]<PROTOCOL['matrix_tolerance'],errs)
  timewrong=max(timewrong,np.linalg.norm(expected1-l1*A[0]))
  secondwrong=max(secondwrong,np.linalg.norm(comm(h,X2[0][1])))
  NUMERICS.append({'fixture':'exact_exponential_Hamiltonian','t':t,'first_errors':err1,'mixed_errors':err2})
 check('omitted_time_connection_detected',timewrong>1e-3,timewrong)
 check('omitted_mixed_generator_detected',secondwrong>1e-3,secondwrong)

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 ADM Dirac boundary response','input_count':len(ps)}),flush=True)
 try:provenance();geometry();boundary();fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 ADM Dirac boundary response','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'ADM_Hamiltonian_derived':good,'metric_first_second_vertices_derived':good,'geometric_projector_transport_derived':good,
  'fixed_domain_contacts_derived':good,'retarded_readout_derived':good,'initial_state_contacts_retained':good,
  'instantaneous_rethermalization_adopted':False,'counterterm_changed':False,'seed_residual_discarded':False,
  'new_boundary_coupling_selected':False,'genuine_boundary_shape_response_evaluated':False,
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
