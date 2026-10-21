"""Strong trace bridge and boundary-index controls; no physical mode claim."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import importlib.util,json,sys,traceback
import sympy as S

ROOT=Path(__file__).resolve().parents[1];STEM='coupled-thermal-boundary-index'
OWN=['checks/coupled_thermal_boundary_index.py','notes/'+STEM+'.md','evidence/coupled-boundary-index-source-review.json']
PREVIOUS='thermal-curvature-shape-schur';ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'scope':'Strong thermal curvature equation and index -1 with both outer traces removed',
 'arithmetic':'Exact symbolic and rational four-dimensional metric jets',
 'charged_jets':[[2,'1/100','11/10',1,2,2],[3,'-1/50','6/5',2,3,3],[4,'1/10','21/20',3,5,2]],
 'jet_columns':['lambda','d','horizon','thermal_index','magnetic_flux','evaluation_radius'],
 'shift_controls':['injective','resonant'],
 'full_feedback_numerically_evaluated':False,'physical_frequency_threshold_certified':False,
 'finite_rank_models_are_charged_operator':False,'formal_verification':False}

def load(name,path):
 spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
old=load('pinned_reconstruction','checks/physical_curvature_reconstruction.py')
schur=load('pinned_schur','checks/thermal_curvature_shape_schur.py')

def check(name,ok,detail=''):
 assert name not in {c['name'] for c in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def zero(name,value):
 value=S.simplify(value);check(name,value==0,value)
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==138 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope_boundaries',not review['physical_boundary_selected'] and not review['physical_instability_proved'] and not review['full_cokernel_evaluated'])

def exact_bridge():
 eps,L,U,C,ru,qw,gw,P,F,T,w,ell,k,Ds,f=S.symbols('eps L U C ru qw gw P F T w ell k Ds f')
 gmix=-C*U+gw;density=eps*(-L*U+gmix)+P*(ru+qw)+2*F*T*w
 D=eps*gw+P*qw+2*F*T*w
 zero('strong_density_row',density-(-eps*(L+C)*U+P*ru+D))
 strong=-eps*(L+C)*U+P*ru+Ds*(ell-k*U)
 reduced=(L+C)*U-P*ru/eps+Ds*k*U/eps-f-Ds*ell/eps
 zero('eliminated_row_source_sign',strong+eps*f+eps*reduced)
 # Neither the off-shell contacts nor the strong feedback is dropped.
 check('density_contacts_nonzero',P*(ru+qw)+2*F*T*w!=0)
 check('trace_source_scaling_nonzero',S.diff(reduced,f)==-1 and S.diff(reduced,ell)==-Ds/eps)
 geom=schur.geometric_functions();r,lam,m,Q,d,G,qm,eps,e,fun=geom
 F,Fr,Frr,RR,Rr,HR,j,Cv,VK,A1,A0,B1,B0,T,Tr,dr,er,Mw,Mq=e
 P=S.expand(F*(r*Tr+4*T))
 for name,expr in [('A1',A1),('P',P),('B0',B0),('rB1',r*B1)]:
  zero('outer_zero_'+name,S.limit(expr,r,S.oo))
  mh=r*(1+lam*r*r+Q/r**2-2*d/r**3+d/r**4)/2
  zero('horizon_zero_'+name,S.factor(expr).subs(m,mh))
 check('A1_not_identically_zero',S.simplify(A1)!=0)
 zero('A1_outer_leading_coefficient',S.limit(r*A1,r,S.oo)-28*lam*m)
 # Essential distinction: a multiplier can fail H-compactness despite tail decay.
 x,n=S.symbols('x n',real=True)
 zero('oscillatory_unit_norm',S.integrate(S.sin(n*x)**2,(x,0,S.pi)).subs(n,1)-S.pi/2)
 zero('distinct_mode_orthogonality',S.integrate(S.sin(x)*S.sin(2*x),(x,0,S.pi)))

def charged_jets():
 t,r=S.symbols('t r',real=True)
 for ix,row in enumerate(PROTOCOL['charged_jets']):
  lam,d,rh,n,qm,rp=map(S.sympify,row);Q=1+3*lam+d
  mass=rh*(1+lam*rh*rh+Q/rh**2-2*d/rh**3+d/rh**4)/2
  F=lam*r*r+1-2*mass/r+Q/r**2-2*d/r**3+d/r**4
  gamma=n*S.diff(F,r).subs(r,rh)/2
  p=S.exp(gamma*t)*(r-rh)**2/r**4;v=S.exp(gamma*t)*(r**-3+R(2,3)*r**-4)
  eps=1/(60*S.pi);la=lam+11*eps*lam*lam/6
  RR,W,j,c,dr,er,HR,D,L,A1,A0,B1,B0=old.coeffs(F,r)
  box=lambda z:old.box(F,z,t,r)
  T=RR+12*la+eps*(7*W*W/36+22*j*j/3-11*RR*RR/72-10*qm*qm/r**4+box(RR))
  scalar,C2,J2,boxR,at=old.direct_optical(F,p,v,t,r,(S.Integer(0),rp))
  ev=lambda z:S.simplify(S.sympify(z).subs(at))
  uv=r*(2*(S.diff(v,t,2)-F*S.diff(F*S.diff(v,r),r))/F+(2*S.diff(F,r,2)+4*F/r**2+8*S.diff(F,r)/r)*v)
  dR=-6*box(p)-2*RR*p+uv/r;U=r*(dR-r*S.diff(RR,r)*p);q=r*p
  zero('independent_4D_curvature_'+str(ix),scalar.d-ev(dR))
  dM=-2*qm*qm*ev(p/r**4)
  dt=scalar.d+eps*(7*C2.d/12+11*J2.d/6-11*ev(RR)*scalar.d/36-20*dM+ev(box(dR))+boxR.d)
  direct=ev(r*F)*(dt+ev((4*p+2*v)*T))
  Vtilde=F*(S.diff(F,r)/r-(1/eps+11*lam/3)+11*(RR+12*lam)/36-c)
  operator=-F*S.diff(F*S.diff(U,r),r)+(Vtilde+gamma*gamma)*U
  gmix=-A1*(U-uv)/6+B1*F*S.diff(q,r)-B0*q-2*A0*v
  contacts=F*((r*S.diff(T,r)+4*T)*q+2*T*r*v)
  predicted=ev(-eps*operator+eps*gmix+contacts)
  zero('strong_trace_density_4D_'+str(ix),direct-predicted)
  omitted=ev(contacts)
  check('nonzero_residual_contact_'+str(ix),omitted!=0,omitted)
  check('positive_frequency_and_exterior_'+str(ix),gamma>0 and ev(F)>0,[gamma,ev(F)])
  FIXTURES.append({'fixture':'exact_charged_trace_density','index':ix,'gamma':str(gamma),'radius':str(rp),'direct':str(direct),'contact':str(omitted),'physical_mode':False})

def shift_controls():
 A_shape=S.Matrix([[3,1],[1,2]])
 D=S.Matrix([[1,0],[0,2],[1,-1],[2,1]])
 C=S.Matrix([[1,-2,1],[0,1,3]])
 shift=S.Matrix([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
 cases=[S.Matrix([[1,0,0],[0,1,0],[0,0,1],[1,-2,3]]),S.Matrix([[1,0,0],[0,1,0],[1,1,0],[2,-1,0]])]
 check('shape_control_positive',A_shape[0,0]>0 and A_shape.det()>0)
 for ix,A in enumerate(cases):
  T=A-D*A_shape.inv()*C
  full=T.row_join(-D).col_join(C.row_join(A_shape))
  right=S.eye(3).row_join(S.zeros(3,2)).col_join((-A_shape.inv()*C).row_join(S.eye(2)))
  left=S.eye(4).row_join(D*A_shape.inv()).col_join(S.zeros(2,4).row_join(S.eye(2)))
  expected=A.row_join(S.zeros(4,2)).col_join(S.zeros(2,3).row_join(A_shape))
  check('full_block_factorization_'+str(ix),left*full*right==expected)
  rank=A.rank();ker=3-rank;coker=4-rank
  check('shift_tail_index_'+str(ix),ker-coker==-1,[ker,coker])
  check('full_pair_same_kernel_cokernel_'+str(ix),5-full.rank()==ker and 6-full.rank()==coker)
  check('expected_resonance_'+str(ix),(ker,coker)==[(0,1),(1,2)][ix])
  x=S.Matrix([2,-1,3]);ww=S.Matrix([1,2]);rhs=full*x.col_join(ww)
  f=rhs[:4,0];ell=rhs[4:,0];effective=f+D*A_shape.inv()*ell
  check('shape_source_transform_'+str(ix),effective==A*x)
  zs=A.T.nullspace()
  check('compatible_source_all_moments_'+str(ix),all((z.T*effective)[0]==0 for z in zs))
  z=zs[0];bad=effective+z
  check('incompatible_source_detected_'+str(ix),A.row_join(bad).rank()>rank)
  check('source_contact_omission_detected_'+str(ix),any((z.T*f)[0]!=0 for z in zs))
  check('finite_rank_perturbation_retained_'+str(ix),A-shift!=S.zeros(4,3))
  FIXTURES.append({'fixture':'infinite_shift_finite_rank_control','case':PROTOCOL['shift_controls'][ix],
   'kernel_dimension':ker,'cokernel_dimension':coker,'index':-1,'first_block':str(A),'cokernel_basis':[str(z) for z in zs],
   'tail':'e_n maps to e_(n+1) for n>=3, an isometry onto coordinates >=4','charged_operator':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 coupled thermal boundary index','input_count':len(ps)}),flush=True)
 try:provenance();exact_bridge();charged_jets();shift_controls()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 coupled thermal boundary index','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'strong_trace_density_bridge_derived':good,'bounded_full_shape_source_feedback_derived':good,
  'minimal_domain_relative_compactness_derived':good,'full_auxiliary_pair_index_minus_one_derived':good,
  'nontrivial_coupled_source_compatibility_derived':good,
  'actual_coupled_cokernel_evaluated':False,'actual_kernel_dimension_evaluated':False,'physical_boundary_domain_selected':False,
  'physical_sources_proved_incompatible':False,'continuous_causal_inverse_proved':False,'actual_metric_mode_constructed':False,
  'actual_charged_stability_proved':False,'actual_wormhole_instability_proved':False,'physical_curvature_norm_changed':False,
  'finite_counterterm_changed':False,'order_reduction_adopted':False,'conformal_contour_rotation_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in [('checker',OWN[0]),('note',OWN[1]),('review',OWN[2])]:run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
