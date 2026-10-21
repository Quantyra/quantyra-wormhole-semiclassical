"""Global thermal shape form: exact identities and independent fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import numpy as np
import scipy
import mpmath as mp
from scipy.integrate import solve_ivp, solve_bvp
from numpy.polynomial.legendre import leggauss
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM='global-thermal-curvature-free-form'
OWN=['checks/global_thermal_curvature_free_form.py','notes/'+STEM+'.md',
     'evidence/global-thermal-shape-source-review.json']
PREVIOUS='global-classical-curvature-free-form';ROWS=[];NUMERICS=[]
R=S.Rational
PROTOCOL={
 'scope':'Global full thermal curvature-free core coercivity; no full inverse or causal stability',
 'seed_fixtures':[[2.,.01,1.1,1,.08,.18,2.],
                  [3.,-.02,1.2,2,.35,.75,3.],
                  [4.,1e-14,1.05,4,1.,1.8,5.]],
 'fixture_columns':['lambda','d','horizon','Matsubara_index','support_left','support_right','magnetic_degeneracy_for_identity'],
 'left_endpoints':[1e-4,5e-5], 'BVP_tolerances':[1e-6,1e-8],
 'quadrature_nodes':[48,96], 'max_nodes':50000,
 'geometry_rtol':2e-13,'geometry_atol':2e-15,
 'right_tail_exponent':12.,'time_quadrature_points':32,
 'action_steps':[2e-3,1e-3], 'normalized_identity_limit':2e-5,
 'nonlinear_difference_evaluation':'exact polynomial times exponential, cancellation-free sinh/cosh symmetric difference; 60-digit direct-density point controls',
 'normalized_refinement_limit':2e-7,
 'IMS_modes':[1,3,8], 'IMS_quadratures':[64,128], 'IMS_identity_limit':2e-10,
 'benchmark_threshold_numerically_enclosed':False,
 'fixtures_are_physical_modes':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,x):
 x=S.simplify(S.sympify(x).doit());check(name,x==0,x)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==73 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('source_scope',len(rev['primary_sources'])==2 and rev['time_contact_generalization_explicit']
       and not rev['full_global_shape_inverse_proved'] and not rev['continuous_causal_bound_proved'])

def exact_checks():
 e,om,rr,uu,vv,yb,kk,pp,tt,de,ds,ww,jj=S.symbols('e omega R U V yB k p T DeltaT DeltaS w j',real=True)
 rv=S.exp(-2*e*tt)*(rr-2*e*de)
 uv=S.exp(-2*e*tt)*(uu+e*ds)
 vvar=S.exp(-2*e*tt)*(vv+2*e*ww+e*e*jj)
 raw=S.exp(2*e*tt)*((om+e*pp)*(rv**2/60-yb*rv/18+kk)+rv*(uv+vvar)/30
      -yb*vvar/18-uv**2/10-R(7,90)*uv*vvar-R(7,180)*vvar**2)
 ff=ds+de/3
 DD=pp*(-rr/15+yb/9)+tt*((6*om+1)*rr/45+11*vv/135)-11*ww/135
 LL=tt*(-rr/15+2*uu/5+7*vv/45)-7*ww/45
 JJ=rr/30-7*(uu+vv)/90-yb/18;WC=-2*rr/15+14*(uu+vv)/45
 ZZ=tt**2*(om*rr**2/30+rr*(uu+vv)/15-uu**2/5-7*uu*vv/45-7*vv**2/90+2*kk*om)+pp*tt*(-rr**2/30+2*kk)
 completed=-ff**2/10+(om+R(1,6))*de**2/15+de*DD+ff*LL+JJ*jj+WC*tt*ww-7*ww**2/45+ZZ
 zero('unrestricted_Weyl_Hessian',S.diff(raw,e,2).subs(e,0)/2-completed)
 zero('Young_shift',R(1,6)+R(3,20)-R(19,60))
 zero('Young_D_square',de**2/100+25*DD**2-de*DD-(de/10-5*DD)**2)
 zero('Young_L_square',ff**2/10-ff*LL-(ff**2/20-5*LL**2)-(ff-10*LL)**2/20)
 # Independent local jets for a full nonzero Fourier mode on a warped cylinder.
 th,tp,tpp,hp,W,Wp,T,Ty,Tyy,g,ainv=S.symbols('theta tp tpp hp W Wp T Ty Tyy g ainv',real=True)
 lap=Tyy+hp*Ty-g*g*ainv*T
 chis=[(S.cos(th),-S.sin(th)*tp,-S.cos(th)*tp**2-S.sin(th)*tpp),
       (S.sin(th),S.cos(th)*tp,-S.sin(th)*tp**2+S.cos(th)*tpp)]
 point=S.expand(sum(W*(c*lap+2*cp*Ty+(cpp+hp*cp)*T)**2 for c,cp,cpp in chis)-W*lap**2)
 expected=6*W*tp**2*Ty**2+2*W*tp**2*ainv*g**2*T**2
 expected+=(8*W*tp*tpp+4*W*hp*tp**2+2*Wp*tp**2)*T*Ty
 expected+=W*(tp**4+(tpp+hp*tp)**2)*T*T
 flux=-2*((Wp*tp**2+2*W*tp*tpp+hp*W*tp**2)*T*Ty+W*tp**2*(Ty**2+T*Tyy))
 zero('full_weighted_IMS_with_flux',S.trigsimp(S.expand(point-expected-flux)))
 zero('time_term_is_required',S.diff(expected,g,2)-4*W*tp**2*ainv*T*T)
 check('static_formula_not_valid_for_nonzero_modes',S.diff(expected,g)!=0)
 s=S.symbols('s',real=True);beta=462*s**6-1980*s**7+3465*s**8-3080*s**9+1386*s**10-252*s**11
 zero('partition_left_value',beta.subs(s,0));zero('partition_right_value',beta.subs(s,1)-1)
 check('partition_C5',all(S.diff(beta,s,j).subs(s,t)==0 for j in range(1,6) for t in [0,1]))
 # Rational coefficient bounds, including their integrated constants.
 checks={
  'HR_negative':24-36*R(1,10**12)>0,
  'D_S_coefficient':R(13,15)+R(1,9)<1,
  'D_T_coefficient':R(91,45)+R(55,135)<3,
  'L_T_coefficient':R(13,15)+R(44,5)+R(7,9)<11,
  'J_coefficient':R(13,30)+R(7*27,90)+R(1,18)<3,
  'W_coefficient':R(26,15)+R(14*27,45)<11,
  'D_integral_budget':3*R(10,27)+27+R(121,6075)*R(50,9)<30,
  'L_integral_budget':242+R(98,2025)*R(50,9)<250,
  'J_integral_budget':3*R(10,9)<4,
  'mixed_gradient_budget':R(50,9)<9,
  'Z_integral_budget':R(169,30)+R(351,15)+R(484,5)+R(770,45)+R(175,90)+R(169,30)<160,
  'loworder_ledger':25*30+5*250+4+33+160==2197,
  'lambda_cost':R(2197,6)<400,
  'magnetic_cost':R(2,9)<1,
  'benchmark_Weyl_cost':400*R(4,10**17)**2+R(4,10**27)<R(1,10**24),
  'canonical_margin':R(33,600)-R(3,1000)*R(22,7)>R(1,100),
  'interior_margin':R(139,10)-6-R(19,60)>R(15,2),
  'full_classical_margin':R(1,32)-R(1,10**24)>R(1001,1000)/40}
 for name,ok in checks.items():check(name,ok)
 A,B,z=S.symbols('A B z',positive=True)
 zero('thermal_absorption',R(1,100)-R(1,400)-R(1,400)-R(1,200))
 # A constant metric rescaling s=C/B^2 has dvol_bar=s*dvol_hat,
 # Delta_bar=Delta_hat/s. Check both derivative and optical-volume factors.
 C,BB,FF,r=S.symbols('C BB FF r',positive=True);scale=C/BB**2
 zero('J2_physical_rescaling',C*scale/scale**2-BB**2)
 zero('J0_optical_rescaling',scale*BB**2*FF/r**2-C*FF/r**2)
 zero('physical_frequency_warp',(S.sqrt(C*FF)/r)**2-C*FF/r**2)

def ims_fixtures():
 x=S.symbols('x');theta=S.pi*(462*x**6-1980*x**7+3465*x**8-3080*x**9+1386*x**10-252*x**11)/2
 thetaf=S.lambdify(x,[theta,S.diff(theta,x),S.diff(theta,x,2)],'numpy')
 for mode in PROTOCOL['IMS_modes']:
  vals=[]
  for count in PROTOCOL['IMS_quadratures']:
   zz,wt=leggauss(count);y=(zz+1)/2;weights=wt/2
   th,tp,tpp=thetaf(y);h=.2*np.sin(2*y);hp=.4*np.cos(2*y);a=np.exp(h)
   W=1.3+.2*y+.1*y*y;Wp=.2+.2*y
   T=1+.3*np.sin(3*y)+.2*y;Ty=.9*np.cos(3*y)+.2;Tyy=-2.7*np.sin(3*y)
   lap=Tyy+hp*Ty-mode**2/a**2*T
   direct=-W*lap*lap
   for c,cp,cpp in [(np.cos(th),-np.sin(th)*tp,-np.cos(th)*tp**2-np.sin(th)*tpp),
                    (np.sin(th),np.cos(th)*tp,-np.sin(th)*tp**2+np.cos(th)*tpp)]:
    direct+=W*(c*lap+2*cp*Ty+(cpp+hp*cp)*T)**2
   time=2*W*tp**2/a**2*mode**2*T*T
   rhs=6*W*tp**2*Ty**2+time+(8*W*tp*tpp+4*W*hp*tp**2+2*Wp*tp**2)*T*Ty+W*(tp**4+(tpp+hp*tp)**2)*T*T
   lhs=float(np.dot(weights*a,direct));right=float(np.dot(weights*a,rhs));missing=float(np.dot(weights*a,time))
   err=abs(lhs-right)/(1+abs(lhs));tag=str(mode)+'_'+str(count)
   check('warped_IMS_direct_'+tag,err<PROTOCOL['IMS_identity_limit'],err)
   check('warped_IMS_missing_time_detected_'+tag,missing/(1+abs(lhs))>1e-3,missing)
   vals.append({'quadrature':count,'direct':lhs,'integrated':right,'time_term':missing,'normalized_error':err})
  diff=abs(vals[0]['direct']-vals[1]['direct'])/(1+abs(vals[1]['direct']))
  check('warped_IMS_refinement_'+str(mode),diff<PROTOCOL['IMS_identity_limit'],diff)
  NUMERICS.append({'fixture':'warped_weighted_IMS','mode':mode,'refinements':vals,'normalized_refinement':diff})

def seed_fixture(index,row):
 lam,d,rh,n,lo,hi,qmag=row;qh=1+3*lam+d
 mass=rh*(1+lam*rh**2+qh/rh**2-2*d/rh**3+d/rh**4)/2
 frh=2*lam*rh+2*mass/rh**2-2*qh/rh**3+6*d/rh**4-4*d/rh**5
 gamma=n*frh/2;L=hi+PROTOCOL['right_tail_exponent']/gamma
 check('positive_thermal_frequency_'+str(index),gamma>0,gamma)
 geo=solve_ivp(lambda y,u:lam+u*u-2*mass*u**3+qh*u**4-2*d*u**5+d*u**6,
       [0,L],[0.],rtol=PROTOCOL['geometry_rtol'],atol=PROTOCOL['geometry_atol'],dense_output=True,max_step=.01)
 check('geometry_solve_'+str(index),geo.success,geo.message)
 def fields(y):
  r=1/geo.sol(y)[0]
  F=lam*r*r+1-2*mass/r+qh/r**2-2*d/r**3+d/r**4
  Fr=2*lam*r+2*mass/r**2-2*qh/r**3+6*d/r**4-4*d/r**5
  Frr=2*lam-4*mass/r**3+6*qh/r**4-24*d/r**5+20*d/r**6
  RR=-12*lam+4*d/r**5-6*d/r**6;Rr=-20*d/r**6+36*d/r**7
  HR=2*RR+r*Rr;Cv=24*lam+4/r**2-8*d/r**5+12*d/r**6
  return r,F,Fr,Frr,HR,Cv,F*(Fr/r-HR/6)
 def shape(y,order=0):
  x=(2*y-lo-hi)/(hi-lo);z=np.clip(x,-1,1);u=1-z*z
  val=[u**6,-12*z*u**5,-12*u**4*(1-11*z*z)][order]
  return np.where(abs(x)<1,val*(2/(hi-lo))**order,0.)
 def src(y):
  r,F,Fr,Frr,HR,Cv,VK=fields(y)
  return r*(-shape(y,2)+(gamma**2+F*Cv/2)*shape(y))
 vals=[]
 for level,(left,tol,nquad) in enumerate(zip(PROTOCOL['left_endpoints'],PROTOCOL['BVP_tolerances'],PROTOCOL['quadrature_nodes'])):
  grid=np.unique(np.concatenate([np.geomspace(left,lo,75),np.linspace(lo,hi,180),np.linspace(hi,L,350)]))
  def ode(y,state):return np.vstack([state[1],(fields(y)[-1]+gamma**2)*state[0]+src(y)/3])
  slope=3/left+(gamma**2+2*lam)*left/7
  sol=solve_bvp(ode,lambda a,b:np.array([a[1]-slope*a[0],b[0]]),grid,np.zeros((2,len(grid))),tol=tol,max_nodes=PROTOCOL['max_nodes'])
  tag=str(index)+'_'+str(level);check('reconstruction_BVP_'+tag,sol.success,sol.message)
  zz,ww=leggauss(nquad)
  edges=np.unique(np.concatenate([np.geomspace(left,lo,12),np.linspace(lo,hi,16),np.linspace(hi,L,40)]))
  y=np.concatenate([(a+b)/2+(b-a)*zz/2 for a,b in zip(edges[:-1],edges[1:])])
  weights=np.concatenate([(b-a)*ww/2 for a,b in zip(edges[:-1],edges[1:])])
  r,F,Fr,Frr,HR,Cv,VK=fields(y);v=shape(y);vy=shape(y,1);vyy=shape(y,2)
  qc,qcy=sol.sol(y);p=qc/r;py=(qcy+F*qc/r)/r
  pyy=2*F*py/r+(gamma**2-F*HR/6)*p+src(y)/(3*r)
  # B=1: physical canonical curvature and radius jets from the charged metric.
  fac=F/r**2;om=np.log(r);rhhat=-r*r*Frr+2*r*Fr-2*F
  U=r*Fr-F;V=F;wy=-F/r;k=1/15-qmag*qmag/3
  nt=PROTOCOL['time_quadrature_points'];angle=2*np.pi*np.arange(nt)/nt
  co=np.cos(angle)[None,:];si=np.sin(angle)[None,:]
  col=lambda x:x[:,None]
  pp=col(p)*co;tt=col(v)*co
  ds=col(2*r*py-r*r*HR*p/6+src(y)*r/(3*F))*co
  de=col((vyy-gamma**2*v)/fac)*co
  wg=col(wy*py/fac)*co
  jj=(col(py)*co)**2+(gamma*col(p)*si)**2;jj/=col(fac)
  om,rr,uu,vv=map(col,[om,rhhat,U,V])
  ff=ds+de/3
  DD=pp*(-rr/15+1/9)+tt*((6*om+1)*rr/45+11*vv/135)-11*wg/135
  LL=tt*(-rr/15+2*uu/5+7*vv/45)-7*wg/45
  JJ=rr/30-7*(uu+vv)/90-1/18;WC=-2*rr/15+14*(uu+vv)/45
  ZZ=tt**2*(om*rr**2/30+rr*(uu+vv)/15-uu**2/5-7*uu*vv/45-7*vv**2/90+2*k*om)+pp*tt*(-rr**2/30+2*k)
  completed=-ff**2/10+(om+1/6)*de**2/15+de*DD+ff*LL+JJ*jj+WC*tt*wg-7*wg**2/45+ZZ
  integrate=lambda arr:float(np.dot(weights*fac,2*np.mean(arr,axis=1)))
  target=-integrate(completed)/(2*np.pi)
  # Direct finite differences of the finite, nonlinear Weyl density with
  # its changing volume and conformal metric, without the Hessian formula.
  # Windows longdouble has only float64 precision. Factoring the EXACT
  # finite action avoids subtracting its large background outside v's support.
  def mul(a,b):
   out=[np.zeros_like(pp) for _ in range(len(a)+len(b)-1)]
   for i,x in enumerate(a):
    for j,z in enumerate(b):out[i+j]+=x*z
   return out
  def add(*terms):
   out=[np.zeros_like(pp) for _ in range(max(len(a) for a in terms))]
   for a in terms:
    for i,x in enumerate(a):out[i]+=x
   return out
  scaled=lambda a,s:[s*x for x in a]
  rp=[rr,-2*de];up=[uu,ds];vp=[vv,2*wg,jj];op=[om,pp]
  poly=add(scaled(mul(op,mul(rp,rp)),1/60),scaled(mul(rp,add(up,vp)),1/30),
           scaled(mul(up,up),-1/10),scaled(mul(up,vp),-7/90),scaled(mul(vp,vp),-7/180))
  qtwo=pp*de/9-jj/18
  def symmetric_difference(eps):
   xx=2*eps*tt;sh=np.sinh(xx);ch=np.cosh(xx);half=np.sinh(xx/2)**2
   return (4*poly[0]*half-2*eps*poly[1]*sh+2*eps**2*poly[2]*ch
           -2*eps**3*poly[3]*sh+2*eps**4*poly[4]*ch+2*eps**2*qtwo
           +4*k*om*half+2*k*eps*pp*sh)/eps**2
  # Three independent direct high-precision nonlinear point evaluations,
  # including both conformal tails; no completed Hessian is used here.
  point_errors=[]
  with mp.workdps(60):
   step=PROTOCOL['action_steps'][1];stable=symmetric_difference(step)
   for yy in [lo/2,(lo+hi)/2,hi+1/gamma]:
    iy=int(np.argmin(abs(y-yy)));it=3
    oval,rval,uval,vval=[mp.mpf(float(a[iy,0])) for a in [om,rr,uu,vv]]
    pval,tval,sval,dval,wval,jval=[mp.mpf(float(a[iy,it])) for a in [pp,tt,ds,de,wg,jj]]
    kval=mp.mpf(k)
    def rawdensity(eps):
     ex=mp.exp(-2*eps*tval);rnew=ex*(rval-2*eps*dval)
     unew=ex*(uval+eps*sval);vnew=ex*(vval+2*eps*wval+eps*eps*jval)
     return mp.exp(2*eps*tval)*((oval+eps*pval)*(rnew*rnew/60-rnew/18+kval)
       +rnew*(unew+vnew)/30-vnew/18-unew*unew/10-mp.mpf(7)*unew*vnew/90-mp.mpf(7)*vnew*vnew/180)
    ee=mp.mpf(step);ref=(rawdensity(ee)+rawdensity(-ee)-2*rawdensity(0))/ee**2
    point_errors.append(float(abs(mp.mpf(float(stable[iy,it]))-ref)/(1+abs(ref))))
  check('direct_60digit_nonlinear_points_'+tag,max(point_errors)<2e-10,point_errors)
  diffs=[]
  for step in PROTOCOL['action_steps']:
   hessian=symmetric_difference(step)
   direct=-integrate(hessian)/(4*np.pi);err=abs(direct-target)/(1+abs(target))
   check('nonlinear_Weyl_difference_'+tag+'_'+str(step),err<PROTOCOL['normalized_identity_limit'],err)
   diffs.append({'step':step,'value':direct,'normalized_error':err})
  check('nonlinear_step_refinement_'+tag,diffs[1]['normalized_error']<.4*diffs[0]['normalized_error']+1e-8,diffs)
  gp=float(np.dot(weights,r*r*(py*py+gamma**2*p*p)))
  sourceenergy=float(-3*np.dot(weights,src(y)*qc))
  gpvia=sourceenergy/9+float(np.dot(weights,F*HR*qc*qc))/6
  gp_error=abs(gp-gpvia)/(1+abs(gp))
  check('global_conformal_gradient_identity_'+tag,gp_error<PROTOCOL['normalized_refinement_limit'],gp_error)
  check('global_conformal_gradient_upper_'+tag,gp<=sourceenergy/9+1e-7, [gp,sourceenergy/9])
  tail=float(sol.sol(hi+1/gamma)[0]);check('conformal_tail_retained_'+tag,abs(tail)>1e-10,tail)
  vals.append({'left_endpoint':left,'BVP_tolerance':tol,'BVP_nodes':len(sol.x),'quadrature':nquad,
    'completed_Weyl_Hessian':target,'nonlinear_differences':diffs,'direct_60digit_point_errors':point_errors,'conformal_gradient_norm':gp,
    'conformal_gradient_identity_error':gp_error,'conformal_tail_sample':tail})
 ref=abs(vals[1]['completed_Weyl_Hessian']-vals[0]['completed_Weyl_Hessian'])/(1+abs(vals[1]['completed_Weyl_Hessian']))
 check('charged_Weyl_refinement_'+str(index),ref<PROTOCOL['normalized_refinement_limit'],ref)
 NUMERICS.append({'fixture':'charged_thermal_Weyl_action','parameters':row,'gamma':gamma,
  'mass':mass,'right_endpoint':L,'refinements':vals,'normalized_refinement':ref,
  'threshold_certified_for_fixture':False,'physical_kernel_evaluated':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 global thermal curvature-free form','input_count':len(ps)}),flush=True)
 try:
  provenance();exact_checks();ims_fixtures()
  for index,row in enumerate(PROTOCOL['seed_fixtures']):seed_fixture(index,row)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 global thermal curvature-free form','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
  'numpy':np.__version__,'scipy':scipy.__version__,'protocol':PROTOCOL,'checks':ROWS,'numerics':NUMERICS,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'global_full_thermal_core_coercivity_derived':good,'time_contacts_retained':good,
  'global_conformal_tails_retained':good,'full_form_continuous_on_stronger_completion_derived':good,
  'benchmark_frequency_threshold_numerically_enclosed':False,'full_global_shape_inverse_proved':False,
  'full_form_closed_in_own_energy_norm_proved':False,'continuous_causal_bound_proved':False,
  'full_quantum_kernel_evaluated':False,'physical_metric_boundary_domain_selected':False,
  'actual_metric_causal_inverse_proved':False,'actual_metric_mode_constructed':False,
  'actual_charged_stability_proved':False,'actual_wormhole_instability_proved':False,
  'counterterm_changed':False,'massive_tower_discarded':False,'retained_fast_factor_removed':False,
  'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in [('checker',OWN[0]),('note',OWN[1]),('review',OWN[2])]:run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
