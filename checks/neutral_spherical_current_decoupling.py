"""S168 fixed-domain particle-hole symmetry and independent CAR fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
STEM='neutral-spherical-current-decoupling'
OWN=['checks/neutral_spherical_current_decoupling.py','notes/'+STEM+'.md','evidence/neutral-current-source-review.json']
PREVIOUS='spherical-retarded-constraint-recovery'
ROWS=[];NUMERICS=[]
PROTOCOL={
 'scope':'Neutral spherical mixed metric/current response in the selected state and domain',
 'Galerkin_basis':'beta plus: 1,sqrt(2)cos(pi x); beta minus: sqrt(2)sin(pi x), x in [0,1]',
 'quadrature_grids':[64,128],'quadrature_tolerance':1e-12,
 'thermal_beta':0.7,'angular_eigenvalue':1.4,'mass_profile':'2+cos(pi x)/3',
 'metric_probe':'2/5+x/4','current_probe':'1+x/3','quench_amplitude':0.17,
 'times':[0.0,0.3,0.8],'matrix_tolerance':2e-11,
 'electric_bias_control':0.2,'Fock_modes_massive_pair':6,'Fock_modes_LLL':3,
 'fixtures_are_exterior_spectral_data':False,'formal_verification':False,
}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,expr):
 if isinstance(expr,S.MatrixBase):
  values=[S.simplify(e) for e in expr];check(name,all(e==0 for e in values))
 else:
  value=S.simplify(S.sympify(expr));check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():
 return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json',
  'notes/stationary-paired-dirac-phase.md','notes/uniform-dirac-reference.md',
  'notes/coupled-fermion-domain.md','notes/dirac-spherical-source-reduction.md',
  'notes/spherical-anomaly-action.md','notes/thermal-shape-response-and-subtractions.md'})

def provenance():
 p=prior()
 check('predecessor_passed',p['status']=='passed' and len(p['checks'])==89 and all(r['passed'] for r in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 for f in ['coupled-fermion-domain','dirac-spherical-source-reduction','spherical-anomaly-action','thermal-shape-response-and-subtractions']:
  check('required_pinned_'+f,'notes/'+f+'.md' in p['input_sha256_after'])
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('source_scope',len(review['primary_sources'])==1 and not review['electric_sector_stability_proved']
       and not review['physical_metric_boundary_selected'] and not review['order_reduction_adopted'])

def exact_algebra():
 I=S.eye(2);sx=S.Matrix([[0,1],[1,0]]);sy=S.Matrix([[0,-S.I],[S.I,0]]);sz=S.diag(1,-1)
 kp=S.kronecker_product
 lam,m,at,ax=S.symbols('lambda m at ax',real=True)
 D=lam*sz;Gamma=sx
 zero('angular_chirality',Gamma*D+D*Gamma)
 zero('angular_real_basis',D.conjugate()-D)
 zero('angular_real_chirality',Gamma.conjugate()-Gamma)
 C=kp(I,sx,Gamma);alpha=kp(sz,sy,I);mass=kp(I,sx,D)
 chi=kp(sz,sx,I);Pi=(S.eye(8)+chi)/2
 transform=lambda B:C*B.conjugate()*C
 zero('antiunitary_square',C*C.conjugate()-S.eye(8))
 zero('actual_folded_projector',transform(Pi)-Pi)
 zero('folded_boundary_isotropic',Pi*alpha*Pi)
 zero('kinetic_coefficient_sign',transform(-S.I*alpha)-S.I*alpha)
 zero('mass_coefficient_sign',transform(mass)+mass)
 zero('electric_density_parity',transform(S.eye(8))-S.eye(8))
 zero('electric_current_parity',transform(alpha)-alpha)
 zero('metric_vertex_parity',transform(m*mass)+m*mass)
 zero('electric_source_reversal',transform(m*mass-at*S.eye(8)-ax*alpha)
      +(m*mass+at*S.eye(8)+ax*alpha))
 for side,proj in [('left',kp((I+sz)/2,I,I)),('right',kp((I-sz)/2,I,I))]:
  zero('individual_'+side+'_charge_parity',transform(proj)-proj)
  zero('individual_'+side+'_metric_parity',transform(proj*mass)+proj*mass)
 for sign in [-1,1]:
  C0=sign*sx;P0=(I+sign*sx)/2
  zero('LLL_projector_'+str(sign),C0*P0.conjugate()*C0-P0)
  zero('LLL_kinetic_'+str(sign),C0*(-S.I*sy).conjugate()*C0-S.I*sy)
  zero('LLL_current_'+str(sign),C0*sy.conjugate()*C0-sy)
 changed=(I+S.Rational(3,5)*sx+S.Rational(4,5)*sz)/2
 zero('changed_projector_still_a_projector',changed**2-changed)
 zero('changed_projector_still_isotropic',changed*sy*changed)
 check('changed_projector_breaks_symmetry',sx*changed.conjugate()*sx!=changed)
 check('single_mass_sign_breaks_symmetry',sx*(m*sx).conjugate()*sx!=-m*sx)
 betaT,energy=S.symbols('betaT energy',real=True)
 zero('thermal_particle_hole_identity',1/(1+S.exp(-betaT*energy))+1/(1+S.exp(betaT*energy))-1)
 check('zero_energy_occupation',S.Rational(1,2)==1/(1+S.exp(0)))
 f,q,r,c=S.symbols('f q r c',nonzero=True,real=True)
 local=c*(f**2+q**2/(2*r**4))
 zero('retained_local_gauge_evenness',local.subs(f,-f)-local)
 zero('retained_local_metric_current_contact',S.diff(local,f,r).subs(f,0))
 x=S.symbols('x',real=True)
 fs=S.sqrt(2)*S.sin(S.pi*x)
 zero('Galerkin_beta_minus_left_trace',fs.subs(x,0))
 zero('Galerkin_beta_minus_right_trace',fs.subs(x,1))

def galerkin(grid):
 z,w=np.polynomial.legendre.leggauss(grid);x=(z+1)/2;w=w/2
 f=np.array([np.ones_like(x),np.sqrt(2)*np.cos(np.pi*x),np.sqrt(2)*np.sin(np.pi*x)])
 df=np.array([np.zeros_like(x),-np.sqrt(2)*np.pi*np.sin(np.pi*x),np.sqrt(2)*np.pi*np.cos(np.pi*x)])
 labels=np.array([1,1,-1]);kin=np.zeros((3,3));jx=np.zeros((3,3),complex)
 m=2+np.cos(np.pi*x)/3;v=.4+x/4;probe=1+x/3
 def beta_form(profile):
  return np.array([[labels[i]*np.sum(w*f[i]*f[j]*profile) if labels[i]==labels[j] else 0
                    for j in range(3)] for i in range(3)])
 jt=np.array([[np.sum(w*f[i]*f[j]*probe) if labels[i]==labels[j] else 0
               for j in range(3)] for i in range(3)])
 for i in range(3):
  for j in range(3):
   if labels[i]!=labels[j]:
    kin[i,j]=labels[i]*np.sum(w*f[i]*df[j])
    jx[i,j]=1j*labels[i]*np.sum(w*f[i]*f[j]*probe)
 return kin,beta_form(m),beta_form(m*v),jt,jx

def ann_operators(n):
 dim=2**n;ops=[]
 for i in range(n):
  out=np.zeros((dim,dim),complex)
  for state in range(dim):
   if (state>>i)&1:
    lower=state&((1<<i)-1);out[state^(1<<i),state]=(-1)**lower.bit_count()
  ops.append(out)
 return ops

def op(B,ann):
 out=-np.trace(B)*np.eye(ann[0].shape[0],dtype=complex)/2
 for i in range(len(ann)):
  for j in range(len(ann)):
   if abs(B[i,j])>1e-16:out+=B[i,j]*(ann[i].conj().T@ann[j])
 return out

def matrix_function(h,fun):
 e,U=np.linalg.eigh(h);return (U*fun(e))@U.conj().T

def numerical_checks():
 values=[galerkin(n) for n in PROTOCOL['quadrature_grids']]
 err=max(np.max(np.abs(a-b)) for a,b in zip(*values))
 check('Galerkin_grid_agreement',err<PROTOCOL['quadrature_tolerance'],err)
 K,M,O,jt,jx=values[-1];sx=np.array([[0,1],[1,0]],complex);sz=np.diag([1.,-1.])
 Srad=np.diag([1.,1.,-1.]);lam=PROTOCOL['angular_eigenvalue'];tol=PROTOCOL['matrix_tolerance']
 H=np.kron(K,np.eye(2))+np.kron(M,lam*sz);V=np.kron(O,lam*sz)
 C=np.kron(Srad,sx);currents=[np.kron(jt,np.eye(2)),np.kron(jx,np.eye(2))]
 check('Galerkin_Hermitian',np.max(abs(H-H.conj().T))<tol)
 check('Galerkin_particle_hole',np.max(abs(C@H.conj()@C+H))<tol)
 check('Galerkin_metric_vertex',np.max(abs(C@V.conj()@C+V))<tol)
 betaT=PROTOCOL['thermal_beta']
 covariance=matrix_function(H,lambda e:1/(1+np.exp(betaT*e)))
 check('thermal_covariance_symmetry',np.max(abs(C@covariance.conj()@C-(np.eye(6)-covariance)))<tol)
 ann=ann_operators(6);Id=np.eye(64)
 car_error=max(np.max(abs(ann[i]@ann[j].conj().T+ann[j].conj().T@ann[i]-(Id if i==j else 0))) for i in range(6) for j in range(6))
 check('independent_Fock_CAR',car_error<tol)
 HF=op(H,ann);VF=op(V,ann)
 rho=matrix_function(HF,lambda e:np.exp(-betaT*e));rho/=np.trace(rho)
 transformed=[sum(C[i,j]*ann[j].conj().T for j in range(6)) for i in range(6)]
 for label,B in [('Hamiltonian',H),('metric',V),('density',currents[0]),('current',currents[1])]:
  lhs=op(B,transformed);rhs=-op(C@B.conj()@C,ann)
  check('CAR_bilinear_transform_'+label,np.max(abs(lhs-rhs))<tol)
  expected=np.trace((covariance-.5*np.eye(6))@B)
  observed=np.trace(rho@op(B,ann))
  check('Fock_thermal_dictionary_'+label,abs(observed-expected)<tol,abs(observed-expected))
 for j,B in enumerate(currents):
  BF=op(B,ann)
  check('current_parity_'+str(j),np.max(abs(C@B.conj()@C-B))<tol)
  check('neutral_thermal_current_'+str(j),abs(np.trace(rho@BF))<tol,np.trace(rho@BF))
  for t in PROTOCOL['times']:
   U=matrix_function(H,lambda e:np.exp(1j*t*e));Bt=U@B@U.conj().T
   UF=matrix_function(HF,lambda e:np.exp(1j*t*e));BtF=UF@BF@UF.conj().T
   one=np.trace(covariance@(Bt@V-V@Bt));full=np.trace(rho@(BtF@VF-VF@BtF))
   tag=f'{j}_{t}'
   check('mixed_retarded_one_particle_'+tag,abs(one)<tol,one)
   check('mixed_retarded_Fock_'+tag,abs(full)<tol,full)
   check('independent_retarded_agreement_'+tag,abs(one-full)<tol,abs(one-full))
   electric=np.trace(covariance@(Bt@B-B@Bt))
   NUMERICS.append({'current':j,'time':t,'mixed_one_particle_abs':float(abs(one)),
                    'mixed_Fock_abs':float(abs(full)),'electric_commutator_abs':float(abs(electric))})
 check('electric_response_not_discarded',max(x['electric_commutator_abs'] for x in NUMERICS)>1e-5)
 H1=H+PROTOCOL['quench_amplitude']*V;H1F=op(H1,ann)
 for t in PROTOCOL['times'][1:]:
  U=matrix_function(H1,lambda e:np.exp(-1j*t*e));Qt=U@covariance@U.conj().T
  UF=matrix_function(H1F,lambda e:np.exp(-1j*t*e));rhot=UF@rho@UF.conj().T
  check('quench_preserves_covariance_'+str(t),np.max(abs(C@Qt.conj()@C-(np.eye(6)-Qt)))<tol)
  for j,B in enumerate(currents):
   one=np.trace((Qt-.5*np.eye(6))@B);full=np.trace(rhot@op(B,ann))
   check(f'causal_metric_quench_current_{j}_{t}',max(abs(one),abs(full),abs(one-full))<tol,max(abs(one),abs(full)))
 # The zero-mode angular sector is a massless radial Dirac field, not a Weyl half-channel.
 Q0=matrix_function(K,lambda e:1/(1+np.exp(betaT*e)));ann0=ann_operators(3);H0F=op(K,ann0)
 rho0=matrix_function(H0F,lambda e:np.exp(-betaT*e));rho0/=np.trace(rho0)
 check('LLL_thermal_covariance',np.max(abs(Srad@Q0.conj()@Srad-(np.eye(3)-Q0)))<tol)
 for j,B in enumerate([jt,jx]):
  check('LLL_Fock_current_'+str(j),abs(np.trace(rho0@op(B,ann0)))<tol)
 # Controls retain a valid finite Hamiltonian but deliberately violate a theorem hypothesis.
 biased=matrix_function(H,lambda e:1/(1+np.exp(betaT*(e-PROTOCOL['electric_bias_control']))))
 check('electric_bias_breaks_neutral_covariance',np.max(abs(C@biased.conj()@C-(np.eye(6)-biased)))>1e-3)
 check('electric_bias_gives_charge',np.real(np.trace(biased-.5*np.eye(6)))>1e-3)
 single=K+2*lam*Srad
 singleQ=matrix_function(single,lambda e:1/(1+np.exp(betaT*e)))
 check('omitted_mass_sign_has_spurious_charge',abs(np.trace(singleQ-.5*np.eye(3)))>.1,np.trace(singleQ-.5*np.eye(3)))

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
 before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 neutral current decoupling','input_count':len(ps)}),flush=True)
 try:provenance();exact_algebra();numerical_checks()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 neutral current decoupling','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'numpy':np.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'fixed_q_particle_hole_symmetry_derived':good,'actual_zero_source_projector_preserved':good,
  'neutral_HH_metric_current_cross_response_zero':good,'neutral_spherical_linear_sector_consistent':good,
  'electric_sector_response_zero':False,'electric_sector_stability_proved':False,
  'physical_metric_boundary_selected':False,'actual_metric_causal_inverse_proved':False,
  'actual_charged_retarded_transfer_proved':False,'actual_wormhole_instability_proved':False,
  'retained_fast_factor_removed':False,'order_reduction_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
