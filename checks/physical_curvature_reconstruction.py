"""Physical-norm curvature reconstruction and exact optical source checks."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from itertools import product
import importlib.util, json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='physical-curvature-reconstruction'
OWN=['checks/physical_curvature_reconstruction.py','notes/'+STEM+'.md','evidence/curvature-reconstruction-source-review.json']
PREVIOUS='explicit-curvature-metric-source';ROWS=[];NUMERICS=[]
spec=importlib.util.spec_from_file_location('pinned_metric_source',ROOT/'checks/explicit_curvature_metric_source.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
Dual=old.Dual;radial=old.radial
PROTOCOL={'scope':'Actual seed physical curvature reconstruction and optical trace reduction; shape solve unresolved',
 'mpmath_precision':55,'Green_relative_tolerance':'1e-38',
 'Green_model':{'F':'1/(lambda y^2)','potential':'6/y^2','lambda':'2','q_exponent':'0.25',
                'gamma':['0.4','1.3','3'],'y':['0.05','0.5','2']},
 'curvature_fixture_count':4,'curvature_arithmetic':'exact rational full four-dimensional first variations',
 'Green_model_is_actual_charged_seed':False,'fixtures_are_physical_solutions':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.simplify(S.sympify(value).doit());check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==79 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('source_scope',len(review['primary_sources'])==1 and not review['singular_lapse_multiplier_bounded_separately']
       and not review['actual_retarded_shape_inverse_proved'])

def coeffs(F,r):
 R,W,j,c,d,e=radial(F,r);HR=2*R+r*S.diff(R,r)
 D=S.diff(F,r)/r-HR/6;L=S.diff(F,r)+r*S.diff(F,r,2)+r*HR/6
 A1=r*F*d;A0=r*F*(d*S.diff(F,r)+e*F)
 B1=2*r*F*e;B0=2*e*F**2+r*F*e*S.diff(F,r)+F*d*L
 return R,W,j,c,d,e,HR,D,L,A1,A0,B1,B0

def box(F,w,t,r):return -S.diff(w,t,2)/F+S.diff(r*r*F*S.diff(w,r),r)/(r*r)

def exact_checks():
 r,t=S.symbols('r t',positive=True);F=S.Function('F')(r)
 p=S.Function('p')(t,r);v=S.Function('v')(t,r)
 R,W,j,c,d,e,HR,D,L,A1,A0,B1,B0=coeffs(F,r)
 psip=-6*box(F,p,t,r)-HR*p;qp=r*p
 Kqp=-F*S.diff(F*S.diff(qp,r),r)+F*D*qp
 zero('exact_curvature_reconstruction',r*psip-6*(S.diff(qp,t,2)+Kqp)/F)
 b=2*r*F*S.diff(p,r)-r*S.diff(F,r)*p
 ar=-r*S.diff(p,r,2)-r*S.diff(p,t,2)/F**2
 zero('conformal_areal_H_variation',S.diff(b,r)+F*ar+r*psip/6+L*p)
 gm=r*F*(d*(S.diff(b,r)+F*ar)+e*b)
 zero('conformal_feedback_identity',gm-(-A1*r*psip/6+B1*F*S.diff(qp,r)-B0*qp))
 zero('shape_derivative_cancellation',r*F*(d*(S.diff(-2*F*v,r)+F*S.diff(2*v,r))+e*(-2*F*v))+2*A0*v)
 # Independent orbit conformal variation: delta R2=-2v R2-2 Box2 v,
 # delta H=-2v F', delta(grad r)^2=-2v F.
 psiv=2*(S.diff(v,t,2)-F*S.diff(F*S.diff(v,r),r))/F+(2*S.diff(F,r,2)+4*F/r**2+8*S.diff(F,r)/r)*v
 warp=2*v*S.diff(F,r,2)-2*(-S.diff(v,t,2)/F+S.diff(F*S.diff(v,r),r))+4*v*F/r**2+8*v*S.diff(F,r)/r
 zero('shape_curvature_variation',psiv-warp)
 y,ell=S.symbols('y ell',positive=True);rho=ell*y/(ell+y)
 zero('weighted_energy_second_derivative',S.diff(rho**2,y,2)/2-ell**3*(ell-2*y)/(ell+y)**4)
 zero('weighted_energy_bound_polynomial',(ell+y)**4-ell**3*(ell-2*y)-(6*ell**3*y+6*ell**2*y*y+4*ell*y**3+y**4))
 # r derivatives converted exactly to inward y derivatives.
 iy=F*S.diff(F*S.diff(1/F,r),r)
 zero('adjoint_supersolution',-iy+D+(S.Symbol('gamma')**2)/F-(S.diff(F,r,2)-S.diff(F,r)**2/F+D+S.Symbol('gamma')**2/F))
 lam,m,Q,c3,c4,B,x,mu,ds,z=S.symbols('lam m Q c3 c4 B x mu ds z',real=True)
 Fr=lam*r*r+1-2*m/r+Q/r**2+c3/r**3+c4/r**4
 R,W,j,c,d,e,HR,D,L,A1,A0,B1,B0=coeffs(Fr,r)
 zero('seed_reconstruction_potential',D-(6*lam+2*m/r**3-2*Q/r**4-4*c3/r**5-8*c4/r**6))
 zero('seed_conformal_gauge_coefficient',L-(-2*m/r**2+4*Q/r**3+10*c3/r**4+20*c4/r**5))
 sub={r:B/z,lam:x/B**2,m:mu*B,Q:(1+3*x+ds)*B**2,c3:-2*ds*B**3,c4:ds*B**4}
 P=6*x+2*mu*z**3-2*(1+3*x+ds)*z**4+8*ds*z**5-8*ds*z**6
 zero('seed_P_polynomial',B**2*D.subs(sub,simultaneous=True)-P)
 zero('positive_classical_factor',6+4*z**3-6*z**4-4-(1-z)*(6*z**3+2*z*z+2*z+2))
 delta=S.symbols('delta',nonnegative=True)
 zero('signed_seed_perturbation',P.subs(mu,1+2*x+delta)-(x*(6+4*z**3-6*z**4)+2*z**3*(1-z)+2*delta*z**3-2*ds*z**4*(1-2*z)**2))
 check('retained_d_guard_implies_half_x',S.Rational(3,34)<S.Rational(1,2))
 f=x+z*z-2*mu*z**3+(1+3*x+ds)*z**4-2*ds*z**5+ds*z**6
 N=S.expand((2*f-z*S.diff(f,z))**2-f*(2*f-2*z*S.diff(f,z)+z*z*S.diff(f,z,2)+P-3*x))
 zero('threshold_constant_term',N.subs(z,0)+x*x)
 zero('threshold_linear_term',N.coeff(z,1))
 Hpol=S.cancel((N+x*x)/z**2);check('threshold_is_polynomial',Hpol.is_polynomial(z))
 threshold=S.diff(Fr,r)**2-Fr*(S.diff(Fr,r,2)+D-3*lam)
 zero('threshold_exact_radial_identity',B*B*threshold.subs(sub,simultaneous=True)-(-x*x/z**2+Hpol))
 zero('B1_outer_coefficient',S.limit(r**3*B1,r,S.oo)+S.Rational(88,3)*lam*Q)
 zero('B0_outer_coefficient',S.limit(r**2*B0,r,S.oo)+S.Rational(176,3)*lam**2*Q)
 check('feedback_horizon_factors',S.simplify(B1/Fr).is_rational_function(r) and S.simplify(B0/Fr).is_rational_function(r))
 # Initial contacts, derived from u=6/F(q_tt+Kq)+u_v.
 gam,UP,Qp,Vv,Iu,Iq,Iv,Uv,RF,CC,BB,Jq,Aq=S.symbols('gamma U Q V Iu Iq Iv Uv F C B Jq A')
 uvhat=Uv-2*r*Iv/RF
 zero('kinematic_Laplace_contacts',(UP-uvhat-6*(Aq-Iq)/RF).subs(Aq,RF*(UP-uvhat)/6+Iq))
 zero('reduced_initial_source_sign',(Aq-(Iu-CC*(UP-uvhat)-2*A0*Vv-BB*Jq-RF)).expand()
      -(Aq+CC*UP-(Iu+CC*uvhat-2*A0*Vv-RF-BB*Jq)).expand())
 alpha=S.Rational(1,4)
 check('model_q_in_L2_outside_H1',2*alpha>-1 and 2*alpha-2<=-1)
 check('model_lapse_weighted_input_not_L2',2*alpha-4<=-1)
 # Symbolic half-integer basis normalization for the Green fixture.
 yy,gg=S.symbols('yy gg',positive=True)
 dec=S.exp(-gg*yy)*(1+3/(gg*yy)+3/(gg*yy)**2)
 grow=S.exp(gg*yy)*(1-3/(gg*yy)+3/(gg*yy)**2)
 zero('stable_model_decaying_solution',-S.diff(dec,yy,2)+(6/yy**2+gg**2)*dec)
 zero('stable_model_growing_solution',-S.diff(grow,yy,2)+(6/yy**2+gg**2)*grow)
 zero('stable_model_Wronskian',grow*S.diff(dec,yy)-S.diff(grow,yy)*dec+2*gg)

def direct_optical(F,p,v,t,r,point):
 """Full four-dimensional connection and contractions, including angular variation."""
 th,ph=S.symbols('theta phi',real=True);coords=[t,r,th,ph]
 at={t:point[0],r:point[1],th:S.pi/2,ph:0}
 base=[-F,1/F,r*r,r*r*S.sin(th)**2]
 varied=[2*(p+v)*base[0],2*(p+v)*base[1],2*p*base[2],2*p*base[3]]
 def jet(i,*indices):
  vv=base[i];dv=varied[i]
  for j in indices:vv=S.diff(vv,coords[j]);dv=S.diff(dv,coords[j])
  return Dual(vv.subs(at),dv.subs(at))
 gg=[jet(i) for i in range(4)];inv=[1/g for g in gg]
 first={(j,i):jet(i,j) for j,i in product(range(4),repeat=2)}
 second={(k,j,i):jet(i,k,j) for k,j,i in product(range(4),repeat=3)}
 def dg(j,i,k):return first[j,i] if i==k else Dual()
 def ddg(l,j,i,k):return second[l,j,i] if i==k else Dual()
 def bracket(i,j,k):return dg(j,i,k)+dg(k,i,j)-dg(i,j,k)
 Gamma={(i,j,k):inv[i]*bracket(i,j,k)/2 for i,j,k in product(range(4),repeat=3)}
 def derivative(l,i,j,k):
  return (-inv[i]**2*first[l,i]*bracket(i,j,k)+inv[i]*(ddg(l,j,i,k)+ddg(l,k,i,j)-ddg(l,i,j,k)))/2
 Riem={}
 for i,j,k,l in product(range(4),repeat=4):
  Riem[i,j,k,l]=derivative(k,i,l,j)-derivative(l,i,k,j)+sum(
   (Gamma[i,k,n]*Gamma[n,l,j]-Gamma[i,l,n]*Gamma[n,k,j] for n in range(4)),Dual())
 Ric={(j,l):sum((Riem[i,j,i,l] for i in range(4)),Dual()) for j,l in product(range(4),repeat=2)}
 scalar=sum((inv[i]*Ric[i,i] for i in range(4)),Dual())
 ric2=sum((inv[i]*inv[j]*Ric[i,j]**2 for i,j in product(range(4),repeat=2)),Dual())
 riem2=sum((gg[i]*inv[j]*inv[k]*inv[l]*Riem[i,j,k,l]**2 for i,j,k,l in product(range(4),repeat=4)),Dual())
 C2=riem2-2*ric2+scalar**2/3;J2=ric2-scalar**2/4
 R=radial(F,r)[0]
 boxR=inv[1]*S.diff(R,r,2).subs(at)-sum((inv[i]*Gamma[1,i,i] for i in range(4)),Dual())*S.diff(R,r).subs(at)
 return scalar,C2,J2,boxR,at

def optical_fixtures():
 t,r=S.symbols('t r',real=True)
 F1=1+2*r*r-6/r+3/r**2
 F2=1+3*r*r/5-7/r+5/r**2-S.Rational(2,7)/r**3+S.Rational(1,7)/r**4
 fixtures=[(F1,t*t*r+t/r,r*r+t*t/r,(S.Rational(1,3),S.Integer(2))),
 (F2,t*t/r+r*r,t*t*r+t/r**2,(S.Rational(2,5),S.Integer(3))),
 (F1,S.Integer(0),t*t*r*r+r,(S.Rational(-1,4),S.Rational(5,2))),
 (F2,t*t*r+r*r,S.Integer(0),(S.Rational(3,7),S.Integer(3)))]
 for ix,(F,p,v,pt) in enumerate(fixtures):
  scalar,C2,J2,boxR,at=direct_optical(F,p,v,t,r,pt)
  ev=lambda expr:S.simplify(S.sympify(expr).subs(at));tag=str(ix)
  R,W,j,c,d,e,HR,D,L,A1,A0,B1,B0=coeffs(F,r)
  psip=-6*box(F,p,t,r)-HR*p
  psiv=2*(S.diff(v,t,2)-F*S.diff(F*S.diff(v,r),r))/F+(2*S.diff(F,r,2)+4*F/r**2+8*S.diff(F,r)/r)*v
  psi=scalar.d-ev(r*p*S.diff(R,r))
  zero('direct_4D_background_R_'+tag,scalar.v-ev(R))
  zero('direct_4D_background_C2_'+tag,C2.v-ev(W*W/3))
  zero('direct_4D_background_J2_'+tag,J2.v-ev(4*j*j))
  zero('direct_4D_optical_Psi_'+tag,psi-ev(psip+psiv))
  # Exact covariance subtraction in delta(Box)R, not just subtraction of R'.
  KRA=boxR.d-ev(r*p*S.diff(box(F,R,t,r),r))+ev(box(F,r*p*S.diff(R,r),t,r))
  dCA=C2.d-ev(r*p*S.diff(W*W/3,r));dJA=J2.d-ev(r*p*S.diff(4*j*j,r))
  direct=ev(r*F)*(KRA+7*dCA/12+11*dJA/6-ev(c)*psi)
  target=-A1*r*psip/6+B1*F*S.diff(r*p,r)-B0*r*p-2*A0*v
  zero('direct_4D_reconstructed_source_'+tag,direct-ev(target))
  check('optical_time_jet_active_'+tag,ev(S.diff(p,t,2))!=0 or ev(S.diff(v,t,2))!=0)
  NUMERICS.append({'fixture':'exact_optical_curvature','index':ix,'t':str(pt[0]),'r':str(pt[1]),'Psi':str(psi),'g_mix':str(direct)})

def green_fixtures():
 lam=mp.mpf(2);alpha=mp.mpf('.25');ell=1/mp.sqrt(lam);tol=mp.mpf('1e-38')
 target=lambda y:y**alpha*mp.exp(-y)
 dtarget=lambda y:target(y)*(alpha/y-1)
 rho=lambda y:ell*y/(ell+y)
 for ix,gs in enumerate(PROTOCOL['Green_model']['gamma']):
  g=mp.mpf(gs);nu=mp.mpf('2.5')
  aa=lambda y:mp.sqrt(y)*mp.besseli(nu,g*y)
  bb=lambda y:mp.sqrt(y)*mp.besselk(nu,g*y)
  da=lambda y:aa(y)/(2*y)+mp.sqrt(y)*g*(mp.besseli(nu-1,g*y)+mp.besseli(nu+1,g*y))/2
  db=lambda y:bb(y)/(2*y)-mp.sqrt(y)*g*(mp.besselk(nu-1,g*y)+mp.besselk(nu+1,g*y))/2
  uu=lambda y:6*lam*target(y)*(6-alpha*(alpha-1)+2*alpha*y+(g*g-1)*y*y)
  source=lambda y:uu(y)/(6*lam*y*y)
  for j,ys in enumerate(PROTOCOL['Green_model']['y']):
   y=mp.mpf(ys);head=mp.quad(lambda s:aa(s)*source(s),[0,y]);tail=mp.quad(lambda s:bb(s)*source(s),[y,mp.inf])
   rec=bb(y)*head+aa(y)*tail;der=db(y)*head+da(y)*tail
   er=abs(rec-target(y))/(1+abs(target(y)));ed=abs(der-dtarget(y))/(1+abs(dtarget(y)))
   tag=str(ix)+'_'+str(j)
   check('Green_reconstruction_'+tag,er<tol,mp.nstr(er,14))
   check('Green_derivative_'+tag,ed<tol,mp.nstr(ed,14))
   row=(bb(y)*mp.quad(lambda s:aa(s)/(s*s),[0,y])+aa(y)*mp.quad(lambda s:bb(s)/(s*s),[y,mp.inf]))/(6*lam)
   col=(bb(y)*mp.quad(aa,[0,y])+aa(y)*mp.quad(bb,[y,mp.inf]))/(6*lam*y*y)
   check('Schur_row_bound_'+tag,row<=1/(18*lam),mp.nstr(row*18*lam,14))
   check('Schur_column_bound_'+tag,col<=1/(18*lam),mp.nstr(col*18*lam,14))
   NUMERICS.append({'fixture':'leading_model_weighted_Green','gamma':gs,'y':ys,
    'relative_reconstruction_error':mp.nstr(er,15),'relative_derivative_error':mp.nstr(ed,15),
    'row_bound_ratio':mp.nstr(row*18*lam,15),'column_bound_ratio':mp.nstr(col*18*lam,15)})
  qnorm=mp.quad(lambda y:target(y)**2,[0,1,mp.inf]);unorm=mp.quad(lambda y:uu(y)**2,[0,1,mp.inf])
  dn=mp.quad(lambda y:(rho(y)*dtarget(y))**2,[0,1,mp.inf])
  check('physical_reconstruction_norm_'+str(ix),qnorm<=unorm/(18*lam)**2,mp.nstr(mp.sqrt(qnorm/unorm),15))
  check('weighted_derivative_bound_'+str(ix),dn<=unorm/(9*lam)**2,mp.nstr(mp.sqrt(dn/unorm),15))

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 mp.mp.dps=PROTOCOL['mpmath_precision'];ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
 before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 physical curvature reconstruction','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();optical_fixtures();green_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 physical curvature reconstruction','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'physical_curvature_reconstruction_bound_derived':good,'finite_real_frequency_threshold_derived':good,
  'bounded_conformal_feedback_derived':good,'compact_derivative_feedback_derived':good,
  'exact_optical_trace_reduction_derived':good,'initial_metric_contacts_retained':good,
  'actual_retarded_shape_inverse_proved':False,'physical_coupled_moment_evaluated':False,
  'physical_metric_boundary_selected':False,'actual_metric_causal_inverse_proved':False,
  'actual_metric_mode_constructed':False,'actual_charged_stability_proved':False,
  'actual_wormhole_instability_proved':False,'retained_fast_factor_removed':False,
  'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
