"""Initial covariance, metric-jet and Ward-source preparation diagnostics."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='causal-preparation-source-map'
PREVIOUS='ward-sources-and-fast-compatibility'
OWN=['checks/causal_preparation_source_map.py','notes/'+STEM+'.md','evidence/causal-preparation-source-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-30','precision_tolerance':'1e-30',
 'Laplace_points':[['2','.7'],['3','0']],'pulse_decay':'1.3','thermal_beta':'.7',
 'energies':['-1','0','2'],'covariance_radius':'.01',
 'pulse_harmonics':[1,2],'pulse_x':['-3','-1.5','-.5'],'kappa':'1',
 'pulse_fixture':'4096*(v+2)^6*(-1-v)^6 on [-2,-1], zero elsewhere; C5 control',
 'theorem_pulse':'any nonnegative nonzero smooth compact function in (-2,-1)',
 'actual_constrained_initial_data_evaluated':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.simplify(S.sympify(value).doit());check(name,value==0,value)

def zeromat(name,value):
 value=value.applyfunc(S.simplify);check(name,value==S.zeros(*value.shape),value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==66 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['initial_preparation_map_derived'] and review['compact_reconstruction_offset_in_range']
       and not review['actual_physical_initial_data_admissibility_proved'])

def initial_identities():
 s,a,Delta,A,D,kappa,c0=S.symbols('s a Delta A D kappa c0',nonzero=True)
 ceff=c0+S.I*kappa*D
 rhs=(ceff-S.I*kappa*(A+s*D)/(s+a))/(s+S.I*Delta)
 direct=c0/(s+S.I*Delta)-S.I*kappa*(A-a*D)/((s+a)*(s+S.I*Delta))
 zero('nonzero_initial_amplitude_covariance',rhs-direct)
 omitted=(c0-S.I*kappa*(A+s*D)/(s+a))/(s+S.I*Delta)
 zero('initial_amplitude_omission_defect',rhs-omitted-S.I*kappa*D/(s+S.I*Delta))
 zero('readout_initial_contact',(A-s*D)*direct+D*c0-(A*direct-D*(s*direct-c0)))
 check('initial_amplitude_defect_nonzero',S.simplify(rhs-omitted)!=0)
 t=S.symbols('t',nonnegative=True);sp,ap=S.symbols('s_pos a_pos',positive=True)
 P0=S.Matrix([[2,1],[1,3]]);P2=S.Matrix([[4,-2],[-2,5]]);P4=S.Matrix([[3,1],[1,-4]])
 X=S.exp(-ap*t)*S.Matrix([1+t+2*t**2+3*t**3,2-t+t**2-2*t**3])
 def lap_component(f):
  poly=S.Poly(S.expand(S.simplify(S.exp(ap*t)*f)),t)
  return sum(coef*S.factorial(power[0])/(sp+ap)**(power[0]+1) for power,coef in poly.terms())
 lap=lambda v:v.applyfunc(lap_component)
 jets=[S.diff(X,t,j).subs(t,0) for j in range(4)]
 J=P2*(sp*jets[0]+jets[1])+P4*(sp**3*jets[0]+sp**2*jets[1]+sp*jets[2]+jets[3])
 lhs=lap(P0*X+P2*S.diff(X,t,2)+P4*S.diff(X,t,4))
 zeromat('full_local_initial_polynomial',lhs-((P0+sp**2*P2+sp**4*P4)*lap(X)-J))
 check('fourth_order_initial_contact_not_zero',P4*jets[3]!=S.zeros(2,1))
 # Noncommuting positive reconstruction matrix and lapse multiplier.
 K=S.Matrix([[4,1,0],[1,5,1],[0,1,6]]);F=S.diag(2,3,7);gamma=S.Rational(5,2)
 q0=S.Matrix([1,2,-1]);q1=S.Matrix([3,1,4]);w0=S.Matrix([2,-1,3]);w1=S.Matrix([-1,2,1])
 Jkin=gamma*q0+q1+(gamma*w0+w1)/3;G=(K+gamma**2*S.eye(3)).inv()
 R=G*F/6;UI=6*F.inv()*Jkin;offset=G*Jkin
 zeromat('affine_kinematic_reconstruction',R*UI-offset)
 zeromat('affine_initial_curvature_contact',6*F.inv()*(K+gamma**2*S.eye(3))*offset
         -6*F.inv()*(gamma*q0+q1)-2*F.inv()*(gamma*w0+w1))
 check('shape_initial_jets_cannot_be_dropped',G*(gamma*w0+w1)/3!=S.zeros(3,1))
 # Rectangular controls verify the full row cancellation, not a physical index theorem.
 L=S.Matrix([[1,2],[0,1],[2,-1]]);DD=S.Matrix([[1,2],[3,-1],[2,1]])
 Kdag=S.Matrix([[2,1],[-1,3]]);epsilon=S.Rational(2,7);L0=L-DD*Kdag/epsilon
 u=S.Matrix([2,-3]);z=L.T.nullspace()[0];fI=L0*u;well=Kdag*u
 zeromat('full_affine_source_in_range',fI+DD*well/epsilon-L*u)
 zero('full_cokernel_offset_cancellation',(z.T*(fI+DD*well/epsilon))[0])
 check('trace_only_offset_does_not_cancel',(z.T*fI)[0]!=0,(z.T*fI)[0])
 FIXTURES.append({'fixture':'initial_jet_and_affine_range_controls','trace_only_defect':str((z.T*fI)[0]),
                  'physical_initial_constraint_admissibility_claimed':False})

def state_fixtures():
 for case,(sr,si) in enumerate(PROTOCOL['Laplace_points']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    s=mp.mpc(sr,si);a=mp.mpf(PROTOCOL['pulse_decay']);beta=mp.mpf(PROTOCOL['thermal_beta'])
    en=list(map(mp.mpf,PROTOCOL['energies']));occ=[1/(1+mp.exp(beta*e)) for e in en]
    C=mp.diag(occ);A=mp.matrix([[1,mp.mpc('.2','.3'),'.4'],[mp.mpc('.2','-.3'),-2,mp.mpc('.1','.2')],['.4',mp.mpc('.1','-.2'),3]])
    D=mp.matrix([['.3',mp.mpc('.1','-.2'),'-.2'],[mp.mpc('.1','.2'),'-.1','.3'],['-.2','.3','.2']])
    c0=mp.matrix([['.1',mp.mpc('.2','.1'),0],[mp.mpc('.2','-.1'),'-.2','.15'],[0,'.15','.1']])
    T=mp.matrix([[1,'.2',0],['.2',2,'.1'],[0,'.1',-1]])
    trace=lambda M:sum(M[i,i] for i in range(3));comm=lambda M:M*C-C*M
    resolvent=lambda M:mp.matrix([[M[i,j]/(s+1j*(en[i]-en[j])) for j in range(3)] for i in range(3)])
    eta=mp.mpf(PROTOCOL['covariance_radius'])
    vals=[e for sign in [-1,1] for e in mp.eighe(C+sign*eta*c0,eigvals_only=True)]
    check(f'finite_covariance_positivity_{case}_{digits}',min(vals)>0 and max(vals)<1)
    check(f'finite_covariance_charge_neutral_{case}_{digits}',abs(trace(c0))<mp.mpf('1e-32'))
    R0=resolvent(c0);bound=mp.norm(c0)/s.real
    check(f'unitary_resolvent_bound_{case}_{digits}',mp.norm(R0)<=bound,mp.nstr(mp.norm(R0)/bound,digits))
    ceff=c0+1j*comm(D)
    K=trace(C*T)-1j*trace((A-s*D)*resolvent(comm(A+s*D)))
    Jc=trace((A-s*D)*resolvent(ceff))+trace(D*c0)
    prediction=K/(s+a)+Jc
    V=A-a*D
    def cov(t,derivative=False):
     entries=[]
     for i in range(3):
      row=[]
      for j in range(3):
       delta=en[i]-en[j];phase=mp.exp(-1j*delta*t);decay=mp.exp(-a*t)
       forcing=-1j*(occ[j]-occ[i])*V[i,j]
       if derivative:val=-1j*delta*c0[i,j]*phase+forcing*(-a*decay+1j*delta*phase)/(1j*delta-a)
       else:val=c0[i,j]*phase+forcing*(decay-phase)/(1j*delta-a)
       row.append(val)
      entries.append(row)
     return mp.matrix(entries)
    check(f'covariance_initial_value_{case}_{digits}',mp.norm(cov(0)-c0)==0)
    direct=mp.quad(lambda t:mp.exp(-s*t)*(trace(A*cov(t))-trace(D*cov(t,True))+trace(C*T)*mp.exp(-a*t)),[0,1,3,8,mp.inf])
    error=abs(direct-prediction)/(1+abs(prediction))
    check(f'initial_preparation_time_integral_{case}_{digits}',error<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(error,digits))
    missing_amp=trace((A-s*D)*resolvent(1j*comm(D)))
    missing_contact=trace(D*c0)
    check(f'initial_amplitude_contact_detected_{case}_{digits}',abs(missing_amp)>mp.mpf('1e-4'),mp.nstr(abs(missing_amp),digits))
    check(f'initial_readout_contact_detected_{case}_{digits}',abs(missing_contact)>mp.mpf('1e-4'),mp.nstr(abs(missing_contact),digits))
    levels.append({'precision':digits,'value_real':mp.nstr(prediction.real,digits),'value_imag':mp.nstr(prediction.imag,digits),
                   'scaled_error':mp.nstr(error,digits),'omitted_initial_amplitude':mp.nstr(abs(missing_amp),digits),
                   'omitted_initial_readout':mp.nstr(abs(missing_contact),digits)})
  with mp.workdps(80):
   va=mp.mpc(levels[0]['value_real'],levels[0]['value_imag']);vb=mp.mpc(levels[1]['value_real'],levels[1]['value_imag'])
   err=abs(va-vb)/(1+abs(vb));check('preparation_precision_'+str(case),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'nonzero_initial_amplitude_time_integral','s':[sr,si],'levels':levels,'charged_continuum_kernel_evaluated':False})

def ward_pulse_fixtures():
 x=S.symbols('x',real=True);gam=S.symbols('gamma',positive=True)
 b=4096*(x+2)**6*(-1-x)**6
 primitive=-S.exp(-gam*x)*sum(S.diff(b,x,k)/gam**(k+1) for k in range(13))
 zero('pulse_primitive_derivative',S.diff(primitive,x)-S.exp(-gam*x)*b)
 Ehat=S.exp(gam*x)*(primitive.subs(x,-1)-primitive)/2
 zero('finite_time_Ward_initial_energy',gam*Ehat-S.diff(Ehat,x)-b/2)
 zero('finite_time_Ward_initial_momentum',gam*Ehat-S.diff(Ehat,x)-b/2)
 E0,J0,W,o0,oh,th,nprime,lprime,Ehat0,Jhat0,Ep,Jp=S.symbols('E0 J0 W o0 oh th nprime lprime Ehat0 Jhat0 Ep Jp')
 F=gam*W*oh+E0-W*o0;G=-W*((nprime-lprime)*oh+lprime*th)+J0
 zero('general_initial_Ward_first',(gam*Ehat0-E0-Jp-W*(gam*oh-o0))-(gam*Ehat0-Jp-F))
 zero('general_initial_Ward_second',(gam*Jhat0-J0-Ep+W*((nprime-lprime)*oh+lprime*th))-(gam*Jhat0-Ep-G))
 fun=S.lambdify((x,gam),primitive,'mpmath')
 for gamma in PROTOCOL['pulse_harmonics']:
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    g=mp.mpf(gamma);bump=lambda v:4096*(v+2)**6*(-1-v)**6 if -2<v<-1 else mp.mpf(0)
    moment=mp.quad(lambda v:mp.exp(-g*v)*bump(v),[-2,mp.mpf('-1.5'),-1])
    analytic=fun(mp.mpf(-1),g)-fun(mp.mpf(-2),g)
    err=abs(moment-analytic)/(1+abs(moment))
    check(f'pulse_horizon_moment_{gamma}_{digits}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
    check(f'nonzero_horizon_moment_{gamma}_{digits}',moment>0,mp.nstr(moment,digits))
    for k,xs in enumerate(PROTOCOL['pulse_x']):
     xx=mp.mpf(xs);lo=max(mp.mpf(0),-2-xx);hi=-1-xx
     val=mp.quad(lambda t:mp.exp(-g*t)*bump(t+xx)/2,[lo,(lo+hi)/2,hi]) if hi>lo else mp.mpf(0)
     lower=max(xx,mp.mpf(-2))
     target=mp.exp(g*xx)*(fun(mp.mpf(-1),g)-fun(lower,g))/2 if lower<-1 else mp.mpf(0)
     er=abs(val-target)/(1+abs(target))
     check(f'pulse_time_transform_{gamma}_{digits}_{k}',er<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(er,digits))
    check(f'first_harmonic_thermal_pole_{gamma}_{digits}',((gamma-2<0)==(gamma==1)))
    levels.append({'precision':digits,'horizon_moment':mp.nstr(moment,digits),'spectral_tensor_V_power':gamma-2})
  with mp.workdps(80):
   err=abs(mp.mpf(levels[0]['horizon_moment'])-mp.mpf(levels[1]['horizon_moment']))/(1+mp.mpf(levels[1]['horizon_moment']))
   check('pulse_precision_'+str(gamma),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'finite_time_regular_ingoing_pulse','gamma':gamma,'levels':levels,'Dirac_state_realizability_claimed':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 causal preparation source map','input_count':len(ps)}),flush=True)
 try:provenance();initial_identities();state_fixtures();ward_pulse_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 causal preparation source map','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed',
  'exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'initial_preparation_map_derived':good,'compact_reconstruction_offset_in_range':good,'finite_time_thermal_identification_refuted':good,
  'actual_physical_initial_data_admissibility_proved':False,'actual_adjoint_profiles_evaluated':False,'actual_charged_two_point_kernel_evaluated':False,
  'physical_fast_branch_excluded':False,'finite_counterterm_changed':False,'new_boundary_prescription_adopted':False,'new_initial_surface_subtraction_adopted':False,
  'order_reduction_adopted':False,'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S075_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]};data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
