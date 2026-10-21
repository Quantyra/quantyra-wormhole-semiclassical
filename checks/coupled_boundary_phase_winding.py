"""Boundary traces, nonnormal extension algebra and a continuous comparator."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]; STEM='coupled-boundary-phase-winding'
PREVIOUS='high-frequency-zero-fast-injectivity'
OWN=['checks/coupled_boundary_phase_winding.py','notes/'+STEM+'.md','evidence/coupled-boundary-phase-source-review.json']
ROWS=[]; FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-29','precision_tolerance':'1e-29',
 'Bessel_cases':[['1.3','2','.7'],['3.7','5','1.2']],'Bessel_columns':['k','gamma','y0'],
 'scaled_cutoffs':['.01','.001','.0001'],'trace_error_coefficient':'20','refinement_ratio_ceiling':'.03',
 'row_model_gamma':[2,4,8],'row_model_remainder_bound':'1',
 'constant_potential_cases':[['1.3','.6','.7','.2'],['3.7','-.4','1.2','-.3']],
 'constant_potential_columns':['k','v0','y0','theta'],'phase_shift_gamma':[4,8,16],'root_indices':[5,6,7],
 'actual_continuous_causal_family_evaluated':False,'actual_benchmark_calibrated':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)}; ROWS.append(row); print(json.dumps(row),flush=True)
 if not ok: raise AssertionError(name)

def zero(name,value):
 value=S.simplify(S.sympify(value).doit()); check(name,value==0,value)

def zeromat(name,value):
 value=value.applyfunc(S.simplify); check(name,value==S.zeros(*value.shape),value)

def prior(): return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths(): return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior(); check('predecessor_passed',p['status']=='passed' and len(p['checks'])==79 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('conditional_scope',review['thermal_maximal_kernel_and_phase_proved'] and review['relative_boundary_trace_gamma_minus_two_bound']
       and review['continuous_resonance_theorem_conditional'] and not review['actual_continuous_family_hypothesis_verified'])

def boundary_algebra():
 y,k,y0,g=S.symbols('y k y0 gamma',positive=True); c,s=S.symbols('A_c A_s',real=True)
 vc=S.sqrt(y)*S.cos(k*S.log(y/y0)); vs=S.sqrt(y)*S.sin(k*S.log(y/y0))
 zero('outer_Green_pairing_sign',vc*S.diff(vs,y)-S.diff(vc,y)*vs-k)
 d=k*S.log(g); rot=S.Matrix([[S.cos(d),S.sin(d)],[-S.sin(d),S.cos(d)]])
 zeromat('trace_dilation_rotation_orthogonal',rot.T*rot-S.eye(2))
 phase=S.symbols('phase',real=True)
 zero('trace_dilation_coefficients',c*S.cos(phase+d)+s*S.sin(phase+d)
      -(rot*S.Matrix([c,s]))[0]*S.cos(phase)-(rot*S.Matrix([c,s]))[1]*S.sin(phase))
 beta,theta=S.symbols('beta theta',real=True)
 zero('reflecting_denominator_phase',S.cos(theta)*S.cos(beta)+S.sin(theta)*S.sin(beta)-S.sin(beta+S.pi/2-theta))

def row_models():
 V=S.Matrix([[S.Rational(1,5),S.Rational(1,10),0,0],[-S.Rational(1,10),S.Rational(3,20),0,0],
             [S.Rational(1,20),S.Rational(1,5),0,0]])
 A0=S.Matrix([[0,1,0,0],[0,0,1,0],[0,0,0,1]]); b=S.Matrix([1,0,0,0]); tau=S.Matrix([[1,0,0,0],[0,1,0,0]])
 norm2=lambda v:S.simplify((v.T*v)[0]); bound=S.simplify(sum(t*t for t in V))
 check('row_remainder_bound',bound<1,bound)
 for gamma in PROTOCOL['row_model_gamma']:
  tag=str(gamma); eta=S.Rational(1,gamma**2); A=gamma**2*A0; Q=A0.T/gamma**2; L=A+V
  h=(S.eye(4)+Q*V).inv()*b; J=(S.eye(4)+Q*V).inv()*Q
  zeromat('maximal_homogeneous_'+tag,L*h)
  zero('maximal_normalization_'+tag,(b.T*h)[0]-1)
  zeromat('maximal_right_inverse_'+tag,L*J-S.eye(3))
  zeromat('maximal_kernel_decomposition_'+tag,S.eye(4)-J*L-h*b.T)
  check('maximal_distance_bound_'+tag,norm2(h-b)<=(eta/(1-eta))**2,norm2(h-b))
  check('maximal_trace_nonzero_'+tag,norm2(tau*h)>0)
  ell=S.Matrix([[2,1]]); den=(ell*tau*h)[0]
  check('generic_phase_denominator_'+tag,den!=0,den)
  Jtheta=J-h*(ell*tau*J)/den
  zeromat('generic_phase_inverse_'+tag,L*Jtheta-S.eye(3))
  zeromat('generic_phase_trace_'+tag,ell*tau*Jtheta)
  ellres=S.Matrix([[-h[1],h[0]]]); basis=h.row_join(S.eye(4)[:,2:4])
  zeromat('resonant_domain_trace_'+tag,ellres*tau*basis)
  check('resonant_range_equals_minimal_'+tag,(L*basis).rank()==2 and (L*basis)[0,:]==S.zeros(1,3))
  check('resonant_kernel_one_dimensional_'+tag,basis.rank()-(L*basis).rank()==1)
  bad=S.Matrix([1,2,-1]); check('resonant_arbitrary_source_obstruction_'+tag,(L*basis).row_join(bad).rank()==3)
  FIXTURES.append({'fixture':'nonnormal_finite_row_model','gamma':gamma,'h':[str(v) for v in h],
   'resonant_row':[str(v) for v in ellres],'distance_squared':str(norm2(h-b)),
   'generic_denominator':str(den),'infinite_dimensional_theorem_inferred_from_rank':False,'charged_kernel_approximation':False})

def trace_at(k,rho,y0,cut):
 y=cut/rho; t=k*mp.log(y/y0)
 u=lambda x:mp.sqrt(x)*mp.besselk(1j*k,rho*x)
 a=u(y)/mp.sqrt(y); d=(y*mp.diff(u,y)-u(y)/2)/(k*mp.sqrt(y))
 return [mp.cos(t)*a-mp.sin(t)*d,mp.sin(t)*a+mp.cos(t)*d]

def analytic_trace(k,rho,y0):
 angle=mp.arg(mp.gamma(1j*k))-k*mp.log(rho*y0/2); amp=abs(mp.gamma(1j*k))
 return [amp*mp.cos(angle),amp*mp.sin(angle)]

def Bessel_traces():
 for case,(ks,gs,ys) in enumerate(PROTOCOL['Bessel_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    k,g,y0=map(mp.mpf,[ks,gs,ys]); target=analytic_trace(k,g,y0); amp=mp.sqrt(sum(v*v for v in target)); records=[]; errors=[]
    for j,cs in enumerate(PROTOCOL['scaled_cutoffs']):
     cut=mp.mpf(cs); val=trace_at(k,g,y0,cut); err=mp.sqrt(sum(abs(a-b)**2 for a,b in zip(val,target)))/amp
     check(f'Bessel_trace_error_{case}_{digits}_{j}',err<mp.mpf(PROTOCOL['trace_error_coefficient'])*cut**2,mp.nstr(err,digits))
     errors.append(err); records.append({'cutoff':cs,'trace_real':[mp.nstr(v.real,digits) for v in val],
      'trace_imag':[mp.nstr(v.imag,digits) for v in val],'relative_error':mp.nstr(err,digits)})
    for j in [1,2]: check(f'Bessel_trace_refinement_{case}_{digits}_{j}',errors[j]/errors[j-1]<mp.mpf(PROTOCOL['refinement_ratio_ceiling']),mp.nstr(errors[j]/errors[j-1],digits))
    levels.append({'precision':digits,'cutoffs':records})
  with mp.workdps(80):
   differences=[]
   for aa,bb in zip(levels[0]['cutoffs'],levels[1]['cutoffs']):
    for ar,ai,br,bi in zip(aa['trace_real'],aa['trace_imag'],bb['trace_real'],bb['trace_imag']):
     av=mp.mpc(ar,ai); bv=mp.mpc(br,bi); differences.append(abs(av-bv)/(1+abs(bv)))
   check('Bessel_trace_precision_'+str(case),max(differences)<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(max(differences),80))
  FIXTURES.append({'fixture':'Bessel_boundary_trace','k':ks,'gamma':gs,'y0':ys,'levels':levels})

def constant_potentials():
 for case,(ks,vs,ys,ts) in enumerate(PROTOCOL['constant_potential_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    k,v,y0,theta=map(mp.mpf,[ks,vs,ys,ts]); phase_records=[]; roots=[]; root_records=[]
    for gamma in PROTOCOL['phase_shift_gamma']:
     rho=mp.sqrt(gamma**2+v); shift=-k*mp.log1p(v/gamma**2)/2
     a=analytic_trace(k,mp.mpf(gamma),y0); b=analytic_trace(k,rho,y0)
     numerical=mp.atan2(a[0]*b[1]-a[1]*b[0],a[0]*b[0]+a[1]*b[1])
     err=abs(numerical-shift)
     check(f'constant_phase_identity_{case}_{digits}_{gamma}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
     check(f'constant_phase_bound_{case}_{digits}_{gamma}',abs(shift)<=k*abs(v)/gamma**2,mp.nstr(abs(shift)*gamma**2/(k*abs(v)),digits))
     phase_records.append({'gamma':gamma,'phase_shift':mp.nstr(shift,digits)})
    for j in PROTOCOL['root_indices']:
     pure=(2/y0)*mp.exp((mp.arg(mp.gamma(1j*k))+mp.pi/2-theta+j*mp.pi)/k)
     root=mp.sqrt(pure*pure-v); rho=mp.sqrt(root*root+v); tr=analytic_trace(k,rho,y0)
     residual=abs(mp.cos(theta)*tr[0]+mp.sin(theta)*tr[1])/abs(mp.gamma(1j*k))
     check(f'constant_root_phase_{case}_{digits}_{j}',residual<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(residual,digits))
     rel=abs(root/pure-1); check(f'constant_root_displacement_{case}_{digits}_{j}',rel<=abs(v)/pure**2,mp.nstr(rel*pure**2/abs(v),digits))
     val=trace_at(k,rho,y0,mp.mpf('.0001')); cutres=abs(mp.cos(theta)*val[0]+mp.sin(theta)*val[1])/abs(mp.gamma(1j*k))
     check(f'constant_root_direct_trace_{case}_{digits}_{j}',cutres<mp.mpf(PROTOCOL['trace_error_coefficient'])*mp.mpf('1e-8'),mp.nstr(cutres,digits))
     roots.append(root); root_records.append({'index':j,'pure':mp.nstr(pure,digits),'root':mp.nstr(root,digits),'direct_trace_residual':mp.nstr(cutres,digits)})
    ratioerr=abs(roots[-1]/roots[-2]/mp.exp(mp.pi/k)-1); previouspure=mp.mpf(root_records[-2]['pure'])
    check(f'constant_root_ratio_{case}_{digits}',ratioerr<=2*abs(v)/previouspure**2,mp.nstr(ratioerr,digits))
    levels.append({'precision':digits,'phase_shifts':phase_records,'roots':root_records})
  with mp.workdps(80):
   errors=[abs(mp.mpf(a['root'])-mp.mpf(b['root']))/(1+abs(mp.mpf(b['root']))) for a,b in zip(levels[0]['roots'],levels[1]['roots'])]
   check('constant_root_precision_'+str(case),max(errors)<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(max(errors),80))
  FIXTURES.append({'fixture':'constant_bounded_potential','k':ks,'v0':vs,'y0':ys,'theta':ts,'levels':levels,'actual_charged_response':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n'); sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths(); hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}; before=hashes(); error=None
 print(json.dumps({'event':'start','story':'S168 coupled boundary phase winding','input_count':len(ps)}),flush=True)
 try: provenance(); boundary_algebra(); row_models(); Bessel_traces(); constant_potentials()
 except Exception: error=traceback.format_exc(); print(error,file=sys.stderr,flush=True)
 after=hashes(); good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 coupled boundary phase winding','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'thermal_maximal_kernel_and_phase_proved':good,'relative_boundary_trace_gamma_minus_two_bound':good,'continuous_resonance_theorem_conditional':good,
  'actual_continuous_family_hypothesis_verified':False,'fixed_phase_Matsubara_resonance_proved':False,'actual_numeric_threshold_certified':False,
  'actual_charged_two_point_kernel_evaluated':False,'physical_fast_branch_excluded':False,'finite_counterterm_changed':False,
  'new_boundary_prescription_adopted':False,'physical_source_projection_adopted':False,'order_reduction_adopted':False,
  'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,'charged_Lorentzian_instability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN): run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json'); data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}; data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n'); print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__': raise SystemExit(main())
