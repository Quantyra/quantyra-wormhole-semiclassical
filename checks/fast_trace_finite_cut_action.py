"""Retained finite-cut action on fast conformal paths; no physical exclusion."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='fast-trace-finite-cut-action'
OWN=['checks/fast_trace_finite_cut_action.py','notes/'+STEM+'.md','evidence/fast-finite-cut-source-review.json']
PREVIOUS='fast-curvature-boundary-symplectic';ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'scope':'Literal finite-cut conformal path test and second-order tail obligation; no physical phase selected',
 'precision_digits':[50,80],'reference_B':'1','base_radius':'100','cycle_indices':[0,1,2],
 'annulus_phase_width_in_pi_units':'1/4','charged_cases':[['2','.01','1.1','2','2'],['3','-.02','1.2','3','5'],['4','1e-14','1.05','5','9']],
 'charged_columns':['lambda','d','horizon','magnetic_flux','k'],
 'fine_relative_limit':'.001','coarse_relative_limit':'.01','refinement_bound':'exp(-pi/k)',
 'precision_relative_limit':'1e-35','contour_nodes':[16,32],'contour_radius':'.03','contour_relative_limit':'1e-28',
 'static_numerical_fixtures':True,'thermal_extension_is_ordinary_power_count':True,
 'fixtures_are_physical_benchmark':False,'physical_boundary_completion_adopted':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.factor(S.cancel(S.sympify(value).doit()));check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==105 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('scope',not rev['physical_fast_branch_excluded'] and not rev['nonlinear_boundary_completion_adopted'] and not rev['physical_stability_proved'])

def boundary_conversion():
 a,ax,axx,s,sx,sxx,sig,Y=S.symbols('a ax axx s sx sxx sig Y')
 R2=-2*(ax+a*a);R0=R2+2*Y;Rx=-2*(axx+2*a*ax)
 Fcal=sig*(-6*axx+R(44,7)*a**3+10*a*Y-9*sxx-9*a*sx-3*ax*s+6*a*a*s+5*Y*s+R(14,3)*a*s*s+R(7,2)*s**3)
 Fcal+=R(18,7)*a*a*s-6*a*sx-R(33,7)*a*s*s-9*s*sx-R(12,7)*s**3
 Jcal=s*(R0-3*sx-3*a*s-2*s*s)-sig*(Rx-3*sxx-3*ax*s-3*a*sx-s*R0+4*a*s*s+3*s**3)
 flux=R(7,540)*(R0*s-sig*Rx)-R(11,180)*(sig*(sxx+ax*s+a*sx)-(sx+a*s)*s+(R2-4*Y)*sig*s/3)
 Ccal=-Fcal-R(7,6)*Jcal+90*flux+3*(sig*Rx-s*R2)
 Gcal=sig*(12*a*ax+R(44,7)*a**3+10*a*Y)-6*ax*s-R(24,7)*a*a*s-6*a*sx-R(96,7)*a*s*s-18*s*sx-R(85,21)*s**3
 zero('literal_heat_path_plus_all_conversion_fluxes',Ccal+Gcal)
 check('omitted_conversion_flux_detected',S.expand(-Fcal+Gcal)!=0)

def exact_actions():
 r,lam,G,F=S.symbols('r lam G F',positive=True);Fr,Frr,qm=S.symbols('Fr Frr qm',real=True)
 p,pr,prr,eta=S.symbols('p pr prr eta',real=True)
 eps=G/(60*S.pi);lam0=lam+R(11,6)*eps*lam*lam
 Rseed=-Frr-4*Fr/r+2*(1-F)/r**2;Rhat=-r*r*Frr+2*r*Fr-2*F
 N=S.sqrt(F)/r;a=(r*Fr/2-F)/S.sqrt(F);ax=-Rhat/2-a*a
 sig=S.log(r)+eta*p;ss=S.sqrt(F)+eta*r*S.sqrt(F)*pr
 sx=r*Fr/2+eta*((r*F+r*r*Fr/2)*pr+r*r*F*prr)
 uu=sx+a*ss;vv=ss*ss
 zero('canonical_base_Box_sigma',uu.subs(eta,0)-(r*Fr-F))
 zero('canonical_conformal_Box_variation',S.diff(uu,eta)-r*r*(F*prr+Fr*pr))
 LW=(sig*((Rhat+2)**2/60-R(11,90)*Rhat-qm*qm/3)+Rhat*uu/30-uu*uu/10+(Rhat/30-R(1,18))*vv-R(7,90)*uu*vv-R(7,180)*vv*vv)/(4*S.pi*r*r)
 box=F*prr+(Fr+2*F/r)*pr
 LEH=r*r*(S.exp(2*eta*p)*(Rseed-6*eta*box-6*eta*eta*F*pr*pr)+6*lam0*S.exp(4*eta*p))/(4*G)
 K=Fr/(2*S.sqrt(F))+2*S.sqrt(F)/r
 GHY=r*r*S.sqrt(F)*S.exp(2*eta*p)*(K+3*eta*S.sqrt(F)*pr)/(2*G)
 EHct=-S.sqrt(F)*(S.sqrt(lam)*r*r*S.exp(3*eta*p)+S.exp(eta*p)/(2*S.sqrt(lam)))/G
 Gcal=sig*(12*a*ax+R(44,7)*a**3+10*a)-6*ax*ss-R(24,7)*a*a*ss-6*a*sx-R(96,7)*a*ss*ss-18*ss*sx-R(85,21)*ss**3
 BW=-N*Gcal/(360*S.pi)
 c0=-221*lam**R(3,2)/(15120*S.pi**2);c1=29*S.sqrt(lam)/(20160*S.pi**2)
 loopct=4*S.pi*S.sqrt(F)*(c0*r*r*S.exp(3*eta*p)+2*c1*S.exp(eta*p))
 bulk=LEH+LW;boundary=GHY+EHct+BW+loopct
 coeff=lambda expr,n:S.expand(S.diff(expr,eta,n).subs(eta,0)/S.factorial(n))
 L1,L2,B1,B2=coeff(bulk,1),coeff(bulk,2),coeff(boundary,1),coeff(boundary,2)
 t=S.symbols('t',positive=True);pz,pzz=S.symbols('pz pzz',real=True)
 normal={F:lam*r*r,Fr:2*lam*r,Frr:2*lam,pr:pz/r,prr:(pzz-pz)/r**2,G:60*S.pi*t/lam}
 leading=[]
 for ix,expr in enumerate([L1,L2,B1,B2]):
  scaled=(expr*(r if ix<2 else 1)*G/(lam*r**3)).subs(normal,simultaneous=True)
  leading.append(S.simplify(S.limit(scaled,r,S.oo)))
 targets=[11*t*p-(9+25*t)*pzz/6-(27+53*t)*pz/6,
  (6+22*t)*p*p-9*p*pz-3*p*pzz-(3+17*t)*pz*pz/2-R(16,3)*t*pz*pzz-R(3,2)*t*pzz*pzz,
  3*pz/2+t*(18*pzz+R(337,7)*pz-R(442,7)*p)/6,
  3*p*pz-R(3,2)*p*p+t*(R(211,7)*pz*pz+18*pz*pzz-R(663,7)*p*p)/6]
 for name,value,target in zip(['linear_bulk','quadratic_bulk','linear_boundary','quadratic_boundary'],leading,targets):zero('full_retained_'+name,value-target)
 args=[r,lam,G,qm,F,Fr,Frr,p,pr,prr]
 return S.lambdify(args,[L1,L2,B1,B2],'mpmath'),S.lambdify([*args,eta],[bulk,boundary],'mpmath'),(p,pz,pzz,t,targets)

def exact_limits(data):
 p,pz,pzz,t,(L1,L2,B1,B2)=data
 f,h=S.symbols('f h',real=True);k2=1/t+R(17,12)
 primitive=-(9+25*t)*pz/6+11*t*p/3
 zero('linear_bulk_primitive',S.diff(primitive,p)*pz+S.diff(primitive,pz)*pzz+3*primitive-L1)
 zero('linear_complete_boundary',primitive+B1-R(3,7)*t*(7*pzz+9*pz-16*p))
 sub={p:f,pz:h-R(3,2)*f,pzz:-k2*f-3*h+R(9,4)*f}
 zero('static_fast_linear_coefficient',(7*pzz+9*pz-16*p).subs(sub)+ (7*k2+R(55,4))*f+12*h)
 check('static_linear_amplitude_map_invertible',(7*k2+R(55,4))**2+144*k2>0)
 l2=S.expand(L2.subs(sub));pp=-(6*t+R(3,2))*f*h+(R(41,18)*t+R(5,12))*f*f
 zero('quadratic_bulk_primitive',S.diff(pp,f)*h-S.diff(pp,h)*k2*f-l2)
 Q=-(3001*t+546)*f*f/504-(71*t+21)*f*h/14-R(167,42)*t*h*h
 zero('full_quadratic_cutoff_term',pp+B2.subs(sub)-Q)
 check('quadratic_double_phase_sine_nonzero',(71*t+21)>0)
 z,k,phi=S.symbols('z k phi',real=True);ff=S.cos(k*z+phi);hh=S.diff(ff,z)
 Aq=-(3001*t+546)/504;Bq=-(71*t+21)/14;Cq=-R(167,42)*t
 mean=(Aq+Cq*k*k)/2;qc=(Aq-Cq*k*k)/2;qs=-Bq*k/2
 zero('quadratic_phase_decomposition',S.trigsimp(Q.subs({f:ff,h:hh})-mean-qc*S.cos(2*(k*z+phi))-qs*S.sin(2*(k*z+phi))))
 H=S.Function('H')(z)
 zero('second_order_tail_boundary_operator',S.exp(3*z)*(7*S.diff(S.exp(-3*z)*H,z,2)+9*S.diff(S.exp(-3*z)*H,z)-16*S.exp(-3*z)*H)-(7*S.diff(H,z,2)-33*S.diff(H,z)+20*H))
 M=S.Matrix([[20-28*k*k,-66*k],[66*k,20-28*k*k]])
 zero('repair_matrix_determinant',M.det()-(20-28*k*k)**2-4356*k*k)
 check('repair_matrix_invertible',S.expand(M.det())==784*k**4+3236*k*k+400)
 hc,hs=S.symbols('hc hs');trial=hc*S.cos(2*k*z)+hs*S.sin(2*k*z)
 out=7*S.diff(trial,z,2)-33*S.diff(trial,z)+20*trial
 zero('repair_matrix_from_derivatives',S.expand(out-(M*S.Matrix([hc,hs]))[0]*S.cos(2*k*z)-(M*S.Matrix([hc,hs]))[1]*S.sin(2*k*z)))
 check('linear_omitted_bulk_integrable',-R(3,2)<-1)
 check('quadratic_omitted_bulk_integrable',-3<-1)
 check('nonlinear_bulk_remainder_integrable',2-3*R(3,2)<-1)
 check('nonlinear_boundary_remainder_vanishes',3-3*R(3,2)<0)
 zero('nonzero_harmonic_first_average',S.integrate(S.cos(z),(z,0,2*S.pi)))
 zero('nonzero_harmonic_second_average',S.integrate(S.cos(z)**2,(z,0,2*S.pi))/(2*S.pi)-R(1,2))

def fixtures(coeff,raw):
 for ix,row in enumerate(PROTOCOL['charged_cases']):
  levels=[]
  for precision in PROTOCOL['precision_digits']:
   with mp.workdps(precision):
    lam,d,rh,qm,k=map(mp.mpf,row);t=1/(k*k-mp.mpf(17)/12);G=60*mp.pi*t/lam
    Qh=1+3*lam+d;m=rh*(1+lam*rh*rh+Qh/rh**2-2*d/rh**3+d/rh**4)/2
    lapse=lambda r:lam*r*r+1-2*m/r+Qh/r**2-2*d/r**3+d/r**4
    radius=mp.mpf(PROTOCOL['base_radius']);phase=-k*mp.log(radius)
    f=lambda z:mp.cos(k*z+phase);df=lambda z:-k*mp.sin(k*z+phase)
    field=lambda r:r**mp.mpf('-1.5')*f(mp.log(r))
    Aq=-(3001*t+546)/504;Bq=-(71*t+21)/14;Cq=-mp.mpf(167)*t/42
    quad=lambda z:Aq*f(z)**2+Bq*f(z)*df(z)+Cq*df(z)**2
    qc=(Aq-Cq*k*k)/2;qs=-Bq*k/2;aa=20-28*k*k;bb=66*k;det=aa*aa+bb*bb
    hc=(-7/(3*t))*(aa*qc+bb*qs)/det;hs=(-7/(3*t))*(-bb*qc+aa*qs)/det
    h=lambda z:hc*mp.cos(2*(k*z+phase))+hs*mp.sin(2*(k*z+phase))
    repair=lambda r:r**-3*h(mp.log(r))
    def args(r,profile):
     return [r,lam,G,qm,lapse(r),mp.diff(lapse,r),mp.diff(lapse,r,2),profile(r),mp.diff(profile,r),mp.diff(profile,r,2)]
    vals=[]
    for level,n in enumerate(PROTOCOL['cycle_indices']):
     left=radius*mp.exp(n*mp.pi/k);right=left*mp.exp(mp.pi/(4*k))
     def values(r):
      c=coeff(*args(r,field));d2=coeff(*args(r,repair));return c[0],c[1],c[1]+d2[0]
     points=[left+(right-left)*j/4 for j in range(5)]
     bulk=[mp.quad(lambda r:values(r)[j],points) for j in range(3)]
     cl,cr=coeff(*args(left,field)),coeff(*args(right,field));dl,dr=coeff(*args(left,repair)),coeff(*args(right,repair))
     actual=[bulk[0]+cr[2]-cl[2],bulk[1]+cr[3]-cl[3],bulk[2]+cr[3]+dr[2]-cl[3]-dl[2]]
     first=lambda r:lam*lam/(140*mp.pi)*r**mp.mpf('1.5')*(-(7*k*k+mp.mpf(55)/4)*f(mp.log(r))-12*df(mp.log(r)))
     expected=[first(right)-first(left),lam/G*(quad(mp.log(right))-quad(mp.log(left)))]
     errs=[abs(actual[j]/expected[j]-1) for j in range(2)]+[abs(actual[2]/expected[1])]
     tag=str(ix)+'_'+str(precision)+'_'+str(level)
     check('annulus_predictions_nonzero_'+tag,abs(expected[0])>mp.mpf('1e-6') and abs(expected[1])>mp.mpf('1e-6'))
     check('annulus_coarse_control_'+tag,max(errs)<mp.mpf(PROTOCOL['coarse_relative_limit']),[mp.nstr(e,18) for e in errs])
     vals.append({'left':mp.nstr(left,precision),'right':mp.nstr(right,precision),'actual':[mp.nstr(a,precision) for a in actual],
      'predicted':[mp.nstr(a,precision) for a in expected],'relative_errors':[mp.nstr(e,25) for e in errs]})
    for j,name in enumerate(['linear','quadratic','second_order_repair']):
     check('fine_annulus_'+name+'_'+str(ix)+'_'+str(precision),mp.mpf(vals[-1]['relative_errors'][j])<mp.mpf(PROTOCOL['fine_relative_limit']))
     ratios=[mp.mpf(vals[n+1]['relative_errors'][j])/mp.mpf(vals[n]['relative_errors'][j]) for n in range(2)]
     check('annulus_refines_'+name+'_'+str(ix)+'_'+str(precision),max(ratios)<mp.exp(-mp.pi/k),[mp.nstr(x,18) for x in ratios])
    # Extract Taylor coefficients from the unexpanded local finite-cut action.
    aargs=args(radius,field);target=coeff(*aargs);contours=[];ar=mp.mpf(PROTOCOL['contour_radius'])
    for nn in PROTOCOL['contour_nodes']:
     sums=[mp.mpc(0) for _ in range(4)]
     for j in range(nn):
      ph=2*mp.pi*j/nn;v=raw(*aargs,ar*mp.exp(1j*ph))
      for hidx,(block,power) in enumerate([(0,1),(0,2),(1,1),(1,2)]):sums[hidx]+=v[block]*mp.exp(-1j*power*ph)/(nn*ar**power)
     err=max(abs(sums[j]-target[j])/(1+abs(target[j])) for j in range(4))
     check('unexpanded_action_Taylor_'+str(ix)+'_'+str(precision)+'_'+str(nn),err<mp.mpf(PROTOCOL['contour_relative_limit']),mp.nstr(err,20))
     contours.append({'nodes':nn,'scaled_error':mp.nstr(err,25)})
    levels.append({'precision':precision,'values':vals,'Taylor_contours':contours,'repair_coefficients':[mp.nstr(hc,25),mp.nstr(hs,25)]})
  with mp.workdps(80):
   for j in range(3):
    a=levels[0]['values'][j]['actual'];b=levels[1]['values'][j]['actual']
    err=max(abs(mp.mpf(x)-mp.mpf(y))/(1+abs(mp.mpf(y))) for x,y in zip(a,b))
    check('precision_refinement_'+str(ix)+'_'+str(j),err<mp.mpf(PROTOCOL['precision_relative_limit']),mp.nstr(err,20))
  FIXTURES.append({'fixture':'charged_literal_cutoff_action','parameters':row,'levels':levels,'second_order_tail_is_physical_solution':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 fast traces in the retained finite-cut action','input_count':len(ps)}),flush=True)
 try:provenance();boundary_conversion();coeff,raw,data=exact_actions();exact_limits(data);fixtures(coeff,raw)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 fast traces in the retained finite-cut action','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'literal_static_linear_boundary_failure_derived':good,'literal_thermal_quadratic_cutoff_oscillation_derived':good,
  'fixed_phase_alone_insufficient_for_action_completion_derived':good,'nonlinear_second_order_tail_obligation_derived':good,
  'physical_fast_branch_excluded':False,'nonlinear_boundary_completion_adopted':False,'full_coupled_boundary_action_solved':False,
  'regulator_independence_for_moving_cuts_proved':False,'actual_shape_response_at_second_order_evaluated':False,
  'continuous_causal_stability_proved':False,'physical_constraint_mode_constructed':False,'actual_wormhole_instability_proved':False,
  'finite_counterterm_changed':False,'physical_curvature_norm_changed':False,'order_reduction_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in [('checker',OWN[0]),('note',OWN[1]),('review',OWN[2])]:run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
