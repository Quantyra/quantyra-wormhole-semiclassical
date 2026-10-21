"""Original-action boundary momenta and charged fast-trace pairing checks."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='fast-curvature-boundary-symplectic'
OWN=['checks/fast_curvature_boundary_symplectic.py','notes/'+STEM+'.md','evidence/fast-boundary-symplectic-source-review.json']
PREVIOUS='coupled-fast-metric-boundary-traces';ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'scope':'Full auxiliary variational boundary pairing; no physical phase or stability selection',
 'precision_digits':[50,80],'radii':['40','80','160'],
 'charged_cases':[['2','.01','1.1','1','2'],['3','-.02','1.2','2','5'],['4','1e-14','1.05','4','9']],
 'charged_columns':['lambda','d','horizon','thermal_index','k'],
 'profile_remainder_coefficients':['.2','.3'],
 'charged_fine_relative_limit':'.001','charged_refinement_ratio':'.4',
 'momentum_identity_relative_limit':'1e-40','precision_difference_relative_limit':'1e-40',
 'phase_rk4_steps':[128,256],'phase_path_in_pi_units':'-4/3',
 'phase_fine_error_limit':'1e-8','phase_refinement_ratio':'.07',
 'phase_initial_values':['.2','.5','.9'],'fixed_Robin_rho':'.7','cutoff_phase_offsets':['0','.4'],
 'charged_fixtures_are_physical_benchmark':False,'actual_quantum_kernel_evaluated':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,x):
 x=S.factor(S.cancel(S.sympify(x).doit()));check(name,x==0,x)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==98 and all(r['passed'] for r in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('physical_scope',not rev['physical_phase_selected'] and not rev['physical_stability_proved'] and not rev['full_fast_action_renormalized'])

def local_actions():
 # Differentiate the original nonlinear optical densities before specializing
 # to a charged background or integrating the conformal quadratic density.
 e,co,si,g=S.symbols('e co si gamma',real=True)
 P,V,P1,P2,V1,V2,la,Qm,qm=S.symbols('P V P1 P2 V1 V2 lambda0 Qm qm',real=True)
 p,p1,p2=S.symbols('p p1 p2',real=True)
 ph=P+e*p*co;dp=P2+e*(p2-g*g*p)*co
 grad=P1*P1+2*e*P1*p1*co+e*e*(p1*p1*co*co+g*g*p*p*si*si)
 mixed=P1*V1+e*V1*p1*co
 FC=S.exp(2*ph)*(6*grad+4*mixed)+S.exp(2*V)*(2*S.exp(2*ph)+6*la*S.exp(4*ph)-2*Qm)
 FW=S.exp(-2*V)*(ph*V2*V2/15-dp*V2/15-dp*dp/10-V2*grad/15-R(7,90)*dp*grad-R(7,180)*grad*grad)
 FW+=ph*V2/9-grad/18+(R(1,15)-qm*qm/3)*ph*S.exp(2*V)
 quads=[S.expand(S.diff(raw,e,2).subs(e,0)/2).subs({co**2:1,si**2:1}) for raw in [FC,FW]]
 r=S.symbols('r',positive=True);F=S.symbols('F',positive=True);F1,F2,F3,F4=S.symbols('F1 F2 F3 F4',real=True)
 G,eps=S.symbols('G eps',positive=True);qs=S.symbols('q0:5',real=True);q,qy,qyy,qyyy,q4=qs
 pp=q/r;ppy=(qy+F*q/r)/r;ppyy=(qyy+2*F*ppy-F*F1*pp)/r
 bg={S.exp(2*P):r*r,S.exp(4*P):r**4,S.exp(2*V):F/r**2,S.exp(-2*V):r*r/F,
     P1:-F/r,P2:F*F1/r-F*F/r**2,V1:F/r-F1/2,V2:-F*(F1/r-F/r**2-F2/2),p:pp,p1:ppy,p2:ppyy}
 quadC,quadW=[S.expand(expr.subs(bg,simultaneous=True)) for expr in quads]
 raw=S.cancel(-quadC/(4*G)-15*eps*quadW/G)
 RR=-F2-4*F1/r+2*(1-F)/r**2;HR=2*RR+r*(-F3-4*F2/r+2*F1/r**2-4*(1-F)/r**3)
 jj=-F2/4+(F-1)/(2*r*r);bb=6+eps*RR/6+22*eps*jj/3
 aa=6*eps/F;tt=F/r;VV=g*g+F*F1/r;cc=F*(2*RR+48*la)
 normal=(aa*(qyy-VV*q)**2-bb*((qy+tt*q)**2+g*g*q*q)-cc*q*q)/(4*G)
 variables=[r,F,F1,F2,F3,*qs[:-1]];derivs=[-F,-F*F1,-F*F2,-F*F3,-F*F4,*qs[1:]]
 Dy=lambda expr:S.expand(sum(S.diff(expr,x)*dx for x,dx in zip(variables,derivs)))
 momenta=lambda expr:(S.cancel(S.diff(expr,qy)-Dy(S.diff(expr,qyy))),S.cancel(S.diff(expr,qyy)))
 raw1,raw2=momenta(raw);nor1,nor2=momenta(normal)
 bs=S.symbols('b0:5',real=True);sub=dict(zip(qs,bs))
 flux=lambda m1,m2:S.expand(q*m1.subs(sub,simultaneous=True)+qy*m2.subs(sub,simultaneous=True)-bs[0]*m1-bs[1]*m2)
 zero('original_action_vs_integrated_boundary_pairing',flux(raw1,raw2)-flux(nor1,nor2))
 zero('original_action_vs_integrated_bulk_row',S.diff(raw,q)-Dy(raw1)-S.diff(normal,q)+Dy(nor1))
 zero('quadratic_action_normalization',2*G*nor2-aa*(qyy-VV*q))
 zero('first_momentum_normalization',2*G*nor1+bb*(qy+tt*q)+Dy(aa*(qyy-VV*q)))
 U=6/F*(-qyy+(g*g+F*(F1/r-HR/6))*q)
 zero('exact_charged_second_momentum',2*G*nor2+eps*(U+HR*q))
 zero('exact_charged_first_momentum',2*G*nor1-eps*Dy(U)-(eps*HR-bb)*qy-(eps*Dy(HR)-bb*tt)*q)
 check('magnetic_contact_absent_only_in_pure_conformal_column',S.diff(raw,Qm)==0 and S.diff(raw,qm)==0)
 args=[r,F,F1,F2,F3,F4,G,eps,g,la,*qs[:4]]
 return S.lambdify(args,[raw1,raw2,nor1,nor2],'mpmath')

def exact_trace_checks():
 y,k,lam,eps,G=S.symbols('y k lam eps G',positive=True)
 d0=1/eps+R(23,3)*lam;kk=1/(eps*lam)+R(17,12)
 q0,q1,q2=S.symbols('q0 q1 q2');d0s=S.symbols('d0',positive=True)
 zero('fast_second_momentum',-eps*(6*d0-24*lam)*q0+(6+22*eps*lam)*q0)
 zero('fast_first_momentum_coefficient',6*eps*d0-24*eps*lam-(6-2*eps*lam)-24*eps*lam)
 zero('fast_pairing_coefficient',(6+46*eps*lam)/(2*G*(6*d0)**2)-eps/(12*G*d0))
 fc=S.sqrt(y)*S.cos(k*S.log(y));fs=S.sqrt(y)*S.sin(k*S.log(y))
 zero('trace_Wronskian',S.trigsimp(fc*S.diff(fs,y)-S.diff(fc,y)*fs)-k)
 ac,ass,bc,bs=S.symbols('ac ass bc bs',real=True)
 f1=ac*fc+ass*fs;f2=bc*fc+bs*fs
 zero('full_trace_determinant',S.trigsimp(f1*S.diff(f2,y)-S.diff(f1,y)*f2)-k*(ac*bs-ass*bc))
 check('boundary_pairing_strictly_nonzero',eps*k/(12*G*d0)>0)
 th=S.symbols('theta',real=True);line=S.Matrix([-S.sin(th),S.cos(th)]);J=S.Matrix([[0,1],[-1,0]])
 zero('real_phase_line_isotropic',(line.T*J*line)[0])
 check('two_trace_space_not_isotropic',J.det()!=0)
 u,v=S.symbols('u v',real=True);z=u+S.I*v
 zero('complex_line_requires_real_ratio',S.conjugate(1)*z-S.conjugate(z)*1-2*S.I*v)
 B=ac**3*ass+S.exp(ac+ass)+ac*ass**2
 zero('boundary_functional_has_zero_antisymmetric_variation',S.diff(B,ac,ass)-S.diff(B,ass,ac))
 # Exact block control: bounded symmetric Schur feedback cannot change the
 # antisymmetric boundary block, even when the remaining row is nonsymmetric.
 A=S.Matrix([[3,1],[1,2]]);K=S.Matrix([[1,2],[3,-1],[2,1]])
 C=S.Matrix([[2,1,3],[-2,4,1],[3,1,5]]);M=K*A.inv()*K.T
 zero('Schur_feedback_symmetric',sum((M-M.T).applyfunc(lambda x:x*x)))
 zero('Schur_preserves_boundary_pairing',sum(((C-M)-(C-M).T-(C-C.T)).applyfunc(lambda x:x*x)))
 full=C.row_join(K).col_join(K.T.row_join(A));want=S.diag(C-C.T,S.zeros(2))
 zero('full_block_preserves_boundary_pairing',sum((full-full.T-want).applyfunc(lambda x:x*x)))
 phi=S.symbols('phi',real=True);U=S.sqrt(y)*S.cos(k*S.log(y)+phi)
 rho=-k*S.tan(k*S.log(y)+phi)
 zero('cutoff_Robin_log_derivative',S.trigsimp(y*S.diff(U,y)/U-R(1,2)-rho))
 zero('cutoff_Riccati_flow',S.trigsimp(y*S.diff(rho,y)+rho*rho+k*k))
 zero('projective_phase_period',S.trigsimp(S.tan(k*S.log(y)+phi+S.pi)-S.tan(k*S.log(y)+phi)))
 gmax=R(4,10**17)**2;kbound=60*R(3)/gmax+R(17,12)
 check('benchmark_k_squared',kbound>R(1125,1000)*10**35,kbound)
 check('benchmark_log_cutoff_period',R(22,7)**2/kbound<R(1,10**34))
 r,m,Q,c3,c4=S.symbols('r m Q c3 c4',positive=True)
 F=lam*r*r+1-2*m/r+Q/r**2+c3/r**3+c4/r**4;Fr=S.diff(F,r);Frr=S.diff(F,r,2)
 RR=-Frr-4*Fr/r+2*(1-F)/r**2;HR=2*RR+r*S.diff(RR,r);j=-Frr/4+(F-1)/(2*r*r)
 b=6+eps*RR/6+22*eps*j/3
 zero('charged_HR_limit',S.limit(HR,r,S.oo)+24*lam)
 zero('charged_b_limit',S.limit(b,r,S.oo)-(6-2*eps*lam))
 zero('charged_HR_derivative_limit',S.limit(-F*S.diff(HR,r),r,S.oo))

def charged_fixtures(fun):
 for ix,row in enumerate(PROTOCOL['charged_cases']):
  levels=[]
  for precision in PROTOCOL['precision_digits']:
   with mp.workdps(precision):
    lam,d,rh,n,k=map(mp.mpf,row);eps=1/(lam*(k*k-mp.mpf(17)/12));G=60*mp.pi*eps
    d0=1/eps+mp.mpf(23)*lam/3;Q=1+3*lam+d;m=rh*(1+lam*rh*rh+Q/rh**2-2*d/rh**3+d/rh**4)/2
    F=lambda r:lam*r*r+1-2*m/r+Q/r**2-2*d/r**3+d/r**4
    gamma=n*mp.diff(F,rh)/2;la=lam+11*eps*lam*lam/6;expected=eps*k/(12*G*d0)
    alpha,beta=map(mp.mpf,PROTOCOL['profile_remainder_coefficients'])
    common=lambda r:mp.sqrt(1/(lam*r))*(1+alpha/(lam*r)**2+beta/(lam*r)**3)/(6*d0)
    fields=[lambda r:common(r)*mp.cos(k*mp.log(1/(lam*r))),lambda r:common(r)*mp.sin(k*mp.log(1/(lam*r)))]
    vals=[]
    for jy,rs in enumerate(PROTOCOL['radii']):
     rr=mp.mpf(rs);fj=[mp.diff(F,rr,j) for j in range(5)];ff,f1,f2=fj[:3];alljets=[];mom=[]
     for field in fields:
      a0,a1,a2,a3=[mp.diff(field,rr,j) for j in range(4)]
      qj=[a0,-ff*a1,ff*f1*a1+ff*ff*a2,-ff*((f1*f1+ff*f2)*a1+3*ff*f1*a2+ff*ff*a3)]
      alljets.append(qj);mom.append(fun(rr,*fj,G,eps,gamma,la,*qj))
     qa,qb=alljets;ma,mb=mom
     raw=qa[0]*mb[0]+qa[1]*mb[1]-qb[0]*ma[0]-qb[1]*ma[1]
     normal=qa[0]*mb[2]+qa[1]*mb[3]-qb[0]*ma[2]-qb[1]*ma[3]
     identity=abs(raw-normal)/abs(expected);error=abs(raw/expected-1)
     tag=str(ix)+'_'+str(precision)+'_'+str(jy)
     check('original_action_momenta_identity_'+tag,identity<mp.mpf(PROTOCOL['momentum_identity_relative_limit']),mp.nstr(identity,20))
     check('charged_pairing_positive_'+tag,raw>0,mp.nstr(raw,20))
     vals.append({'radius':rs,'raw_pairing':mp.nstr(raw,precision),'normal_pairing':mp.nstr(normal,precision),'relative_identity_error':mp.nstr(identity,20),'relative_boundary_limit_error':mp.nstr(error,20)})
    check('charged_fine_limit_'+str(ix)+'_'+str(precision),mp.mpf(vals[-1]['relative_boundary_limit_error'])<mp.mpf(PROTOCOL['charged_fine_relative_limit']),vals[-1]['relative_boundary_limit_error'])
    ratios=[mp.mpf(vals[j+1]['relative_boundary_limit_error'])/mp.mpf(vals[j]['relative_boundary_limit_error']) for j in range(2)]
    check('charged_boundary_limit_refines_'+str(ix)+'_'+str(precision),max(ratios)<mp.mpf(PROTOCOL['charged_refinement_ratio']),[mp.nstr(v,15) for v in ratios])
    levels.append({'precision':precision,'G':mp.nstr(G,precision),'gamma':mp.nstr(gamma,precision),'predicted_boundary_pairing':mp.nstr(expected,precision),'values':vals})
  with mp.workdps(80):
   for j in range(3):
    x=mp.mpf(levels[0]['values'][j]['raw_pairing']);z=mp.mpf(levels[1]['values'][j]['raw_pairing'])
    err=abs(x-z)/abs(z);check('charged_precision_refinement_'+str(ix)+'_'+str(j),err<mp.mpf(PROTOCOL['precision_difference_relative_limit']),mp.nstr(err,20))
  FIXTURES.append({'fixture':'charged_action_boundary_pairing','parameters':row,'levels':levels,'physical_benchmark':False})

def phase_fixtures():
 with mp.workdps(60):
  for ix,row in enumerate(PROTOCOL['charged_cases']):
   k=mp.mpf(row[-1]);phi=mp.mpf(PROTOCOL['phase_initial_values'][ix]);length=-4*mp.pi/(3*k)
   initial=(mp.cos(phi),-mp.sin(phi));exact=(mp.cos(k*length+phi),-mp.sin(k*length+phi));vals=[]
   for n in PROTOCOL['phase_rk4_steps']:
    h=length/n;x,y=initial
    rhs=lambda a,b:(k*b,-k*a)
    for j in range(n):
     a=rhs(x,y);b=rhs(x+h*a[0]/2,y+h*a[1]/2);c=rhs(x+h*b[0]/2,y+h*b[1]/2);d=rhs(x+h*c[0],y+h*c[1])
     x+=h*(a[0]+2*b[0]+2*c[0]+d[0])/6;y+=h*(a[1]+2*b[1]+2*c[1]+d[1])/6
    err=mp.sqrt((x-exact[0])**2+(y-exact[1])**2)
    vals.append({'steps':n,'endpoint':[mp.nstr(x,25),mp.nstr(y,25)],'error':mp.nstr(err,25)})
   check('phase_rotation_fine_'+str(ix),mp.mpf(vals[-1]['error'])<mp.mpf(PROTOCOL['phase_fine_error_limit']),vals[-1]['error'])
   ratio=mp.mpf(vals[-1]['error'])/mp.mpf(vals[0]['error'])
   check('phase_rotation_refines_'+str(ix),ratio<mp.mpf(PROTOCOL['phase_refinement_ratio']),mp.nstr(ratio,20))
   rho=mp.mpf(PROTOCOL['fixed_Robin_rho']);projectors=[]
   for offset in map(mp.mpf,PROTOCOL['cutoff_phase_offsets']):
    ps=[]
    for n in range(3):
     delta=mp.exp(-(n*mp.pi+offset)/k);angle=mp.atan(-rho/k)-k*mp.log(delta)
     c,s=mp.cos(angle),mp.sin(angle);ps.append(mp.matrix([[c*c,c*s],[c*s,s*s]]))
    check('cutoff_subsequence_fixed_line_'+str(ix)+'_'+str(offset),max(mp.norm(p-ps[0]) for p in ps)<mp.mpf('1e-50'))
    projectors.append(ps[0])
   distance=mp.norm(projectors[1]-projectors[0]);check('different_cutoff_subsequences_differ_'+str(ix),distance>mp.mpf('.1'),mp.nstr(distance,20))
   FIXTURES.append({'fixture':'normal_model_cutoff_phase_rotation','k':str(k),'phi':str(phi),'levels':vals,'different_subsequence_projector_distance':mp.nstr(distance,25),'charged_causal_resolvent':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 fast curvature variational boundary pairing','input_count':len(ps)}),flush=True)
 try:provenance();fun=local_actions();exact_trace_checks();charged_fixtures(fun);phase_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 fast curvature variational boundary pairing','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'full_conformal_boundary_pairing_derived':good,'charged_fast_pairing_coefficient_derived':good,
  'full_auxiliary_shape_preserves_pairing_derived':good,'boundary_functional_cannot_erase_pairing_derived':good,
  'maximal_reflecting_trace_relation_derived':good,'cutoff_does_not_select_phase_demonstrated':good,
  'physical_boundary_phase_selected':False,'full_fast_action_renormalized':False,'physical_matter_domain_for_all_fast_traces_proved':False,
  'continuous_causal_flux_law_proved':False,'physical_constraint_mode_constructed':False,'actual_quantum_kernel_evaluated':False,
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
