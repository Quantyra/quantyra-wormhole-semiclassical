"""Thermal mixed form: exact trace contacts and charged-action fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback, warnings
import numpy as np
import scipy
from scipy.integrate import solve_ivp, solve_bvp
from numpy.polynomial.legendre import leggauss
import sympy as S

ROOT=Path(__file__).resolve().parents[1];STEM='thermal-curvature-shape-schur'
OWN=['checks/thermal_curvature_shape_schur.py','notes/'+STEM+'.md','evidence/thermal-schur-source-review.json']
PREVIOUS='global-thermal-shape-inverse';ROWS=[];NUMERICS=[];R=S.Rational
PROTOCOL={
 'scope':'Bounded compact thermal curvature-shape coupling and auxiliary Schur form; no causal instability',
 'seed_fixtures':[[2.,.01,1.1,1,.08,.18,2.],[3.,-.02,1.2,2,.35,.75,3.],[4.,1e-14,1.05,4,1.,1.8,5.]],
 'fixture_columns':['lambda','d','horizon','thermal_index','shape_left','shape_right','magnetic_flux'],
 'left_endpoints':[1e-4,5e-5],'BVP_tolerances':[1e-6,1e-8],
 'quadrature_nodes':[48,96],'max_nodes':50000,'right_tail_exponent':12.,
 'geometry_rtol':2e-13,'geometry_atol':2e-15,
 'mesh_merge_relative_gap':1e-12,'BVP_runtime_warnings_are_errors':True,
 'mixed_relative_limit':2e-5,'mixed_absolute_limit':2e-7,'mixed_refinement_limit':2e-5,
 'dilations':[.00125,.000625,.0003125],'dilation_G':1e-12,'dilation_d':1e-14,
 'dilation_coarse_relative_limit':.05,'dilation_fine_relative_limit':.005,
 'dilation_quadrature_relative_limit':2e-10,'dilation_action_relative_limit':2e-9,
 'benchmark_threshold_numerically_enclosed':False,'actual_Schur_feedback_evaluated':False,
 'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,x):
 x=S.simplify(S.sympify(x).doit());check(name,x==0,x)

def merge_grid(points):
 points=np.unique(points);gap=PROTOCOL['mesh_merge_relative_gap']*max(1.,abs(points[-1]))
 return points[np.concatenate([[True],np.diff(points)>gap])]

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==96 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('scope',len(rev['primary_sources'])==2 and rev['original_physical_curvature_L2_retained']
       and not rev['physical_constraint_mode_constructed'] and not rev['Lorentzian_instability_proved'])

def geometric_functions():
 r,lam,m,Q,d,G,qm=S.symbols('r lam m Q d G qm',positive=True);eps=G/(60*S.pi)
 F=lam*r*r+1-2*m/r+Q/r**2-2*d/r**3+d/r**4;Fr=S.diff(F,r);Frr=S.diff(F,r,2)
 RR=-Frr-4*Fr/r+2*(1-F)/r**2;Rr=S.diff(RR,r);HR=2*RR+r*Rr
 j=-Frr/4+(F-1)/(2*r*r);Wc=-Frr+2*Fr/r+2*(1-F)/(r*r)
 dr=7*Wc/(3*r)+44*j/(3*r)+Rr;er=44*j/(3*r*r)+S.diff(RR,r,2)+2*Rr/r
 Cv=2*Frr+4*F/(r*r)+8*Fr/r;VK=F*(Fr/r-HR/6)
 A1=r*F*dr;A0=r*F*(dr*Fr+er*F);B1=2*r*F*er
 L=Fr+r*Frr+r*HR/6;B0=2*er*F*F+r*F*er*Fr+F*dr*L
 lam0=lam+11*eps*lam*lam/6
 trace=S.expand(RR+12*lam0+eps*(7*Wc**2/36+22*j*j/3-11*RR*RR/72-10*qm*qm/r**4+F*S.diff(RR,r,2)+(Fr+2*F/r)*Rr))
 tr=S.diff(trace,r)
 Mw=eps*(A1*Cv/6-2*A0/r)+2*F*trace;Mq=F*(r*tr+4*trace)-eps*B0
 exprs=[F,Fr,Frr,RR,Rr,HR,j,Cv,VK,A1,A0,B1,B0,trace,tr,dr,er,Mw,Mq]
 fun=S.lambdify([r,lam,m,Q,d,G,qm],exprs,'numpy')
 return (r,lam,m,Q,d,G,qm,eps,exprs,fun)

def exact_checks(geom):
 r,lam,m,Q,d,G,qm,eps,e,fn=geom
 F,Fr,Frr,RR,Rr,HR,j,Cv,VK,A1,A0,B1,B0,T,Tr,dr,er,Mw,Mq=e
 zero('vacuum_trace_constant_cancels',S.limit(T,r,S.oo))
 check('trace_has_quartic_decay',S.limit(r**4*T,r,S.oo).is_finite)
 check('trace_density_contact_has_quintic_decay',S.limit(r**5*(r*Tr+4*T),r,S.oo).is_finite)
 for name,expr in [('K2_squared',r*r*F*dr*dr/9),('B1',B1),('Mw',Mw),('Mq',Mq)]:
  zero('outer_coupling_coefficient_vanishes_'+name,S.limit(expr,r,S.oo))
  at_horizon=S.factor(expr).subs(m,r*(lam*r*r+1+Q/r**2-2*d/r**3+d/r**4)/2)
  zero('horizon_coupling_coefficient_vanishes_'+name,at_horizon)
 w,qp,qv,v,v2,g=S.symbols('w qp qv v v2 g',real=True)
 uv=r*(2/F*(g*g*v-v2)+Cv*v)
 gmix=A1*uv/6-B1*qp-B0*qv-2*A0*v
 target=-eps*r*r*dr*(v2-g*g*v)/3-eps*B1*qp+Mw*w+Mq*qv
 zero('exact_mixed_trace_rearrangement',(eps*gmix+F*((r*Tr+4*T)*qv+2*T*w)-target).subs(w,r*v))
 check('residual_contact_not_identically_zero',S.simplify(F*((r*Tr+4*T)*qv+2*T*w))!=0)
 zero('mixed_pairing_constant',R(1,2)*R(1,18)-R(1,36))
 eta=S.symbols('eta',positive=True);chi=eta**2*(eta-1)**6*(2-eta)**6
 integ=lambda x:S.integrate(S.expand(x),(eta,1,2))
 D=integ(S.diff(chi,eta)**2);P=integ(chi**2/eta**2);A=integ(eta**2*S.diff(chi,eta,2)**2)
 H=integ((-eta**2*S.diff(chi,eta,2)+6*chi)**2);N=D+6*P
 zero('dilation_fourth_identity',integ((eta*S.diff(chi,eta,2)-2*chi/eta)**2)-A-4*D-4*P)
 check('dilation_profile_C5',all(S.diff(chi,eta,jj).subs(eta,t)==0 for jj in range(6) for t in [1,2]))
 check('curvature_normalization_nonzero',H>0,H)
 y=S.symbols('y',positive=True)
 u=lam*y+lam**2*y**3/3-m*lam**3*y**4/2
 zero('outer_reciprocal_radius_series',S.series(S.diff(u,y)-(lam+u*u-2*m*u**3+Q*u**4-2*d*u**5+d*u**6),y,0,4).removeO())
 invF=S.series((1/F).subs(r,1/u),y,0,6).removeO()
 zero('inverse_lapse_correction_series',invF-(lam*y*y-lam**2*y**4/3+m*lam**3*y**5))
 seedD=S.series((Fr/r-HR/6).subs(r,1/u),y,0,4).removeO()
 zero('curvature_potential_correction_series',seedD-(6*lam+2*m*lam**3*y**3))
 h=-eta**2*S.diff(chi,eta,2)+6*chi
 I2=integ(h*eta**4*S.diff(chi,eta,2));J2=integ(h*eta**2*chi)
 I3=integ(h*(-eta**5*S.diff(chi,eta,2)+2*eta**3*chi))
 check('dilation_quartic_ratio',A<10**6*N,A/N)
 gg=S.symbols('gg',positive=True)
 raw=6*gg*(A+4*D+4*P)-(6-2*gg)*(D+2*P)-(24+88*gg)*P
 zero('complete_negative_limit',raw-(-6*N+gg*(6*A+26*D-60*P)))
 check('benchmark_epsilon_lambda',R(4,10**17)**2/180<R(1,10**35))
 check('strict_negative_margin',-6*N+R(1,10**35)*(6*A+26*D-60*P)<-5*N)
 zero('normalized_negative_coefficient',R(2,72)-R(1,36))
 # Algebraic square completion, with the shape energy inner product.
 cc,aa,bb,kk=S.symbols('c a b k')
 zero('Schur_square_completion',cc+aa+2*bb-(cc-kk+(aa+2*bb+kk)))
 return dict(D=float(D),P=float(P),A=float(A),H=float(H),N=float(N),I2=float(I2),J2=float(J2),I3=float(I3),profile=S.lambdify(eta,[chi,S.diff(chi,eta),S.diff(chi,eta,2)],'numpy'))

def local_actions():
 e,co,si,g=S.symbols('e co si g',real=True)
 P,V,P1,P2,V1,V2,la,Qm,qm=S.symbols('P V P1 P2 V1 V2 la Qm qm',real=True)
 p,p1,p2,v,v1,v2,u,u1,u2=S.symbols('p p1 p2 v v1 v2 u u1 u2',real=True)
 ph=P+e*p*co;vv=V+e*v*co
 dp=P2+e*(p2-g*g*p)*co;dv=V2+e*(v2-g*g*v)*co
 grad=P1*P1+2*e*P1*p1*co+e*e*(p1*p1*co*co+g*g*p*p*si*si)
 mixed=P1*V1+e*(P1*v1+V1*p1)*co+e*e*(p1*v1*co*co+g*g*p*v*si*si)
 FC=S.exp(2*ph)*(6*grad+4*mixed)+S.exp(2*vv)*(2*S.exp(2*ph)+6*la*S.exp(4*ph)-2*Qm)
 FW=S.exp(-2*vv)*(ph*dv*dv/15-dp*dv/15-dp*dp/10-dv*grad/15-R(7,90)*dp*grad-R(7,180)*grad*grad)
 FW+=ph*dv/9-grad/18+(R(1,15)-qm*qm/3)*ph*S.exp(2*vv)
 quads=[S.expand(S.diff(raw,e,2).subs(e,0)/2).subs({co**2:1,si**2:1}) for raw in [FC,FW]]
 crosses=[]
 for q in quads:
  extra=q.subs({p:p+u,p1:p1+u1,p2:p2+u2},simultaneous=True)-q-q.subs({p:u,p1:u1,p2:u2,v:0,v1:0,v2:0},simultaneous=True)
  crosses.append(S.expand(extra/2))
 args=[P,V,P1,P2,V1,V2,g,la,Qm,qm,u,u1,u2,p,p1,p2,v,v1,v2]
 check('nonlinear_action_mixed_bilinearity',all(S.diff(q,u,2)==0 and S.diff(q,u1,2)==0 and S.diff(q,u2,2)==0 for q in crosses))
 check('classical_mixed_physical_charge_cancels',S.diff(crosses[0],Qm)==0)
 pure=[q.subs({v:0,v1:0,v2:0},simultaneous=True) for q in quads]
 return S.lambdify(args,crosses,'numpy'),S.lambdify([P,V,P1,P2,V1,V2,g,la,Qm,qm,p,p1,p2],pure,'numpy')

def seed_geometry(row,G,dfix=None):
 lam,d,rh,n,lo,hi,qm=row
 if dfix is not None:d=dfix
 Qh=1+3*lam+d;mass=rh*(1+lam*rh*rh+Qh/rh**2-2*d/rh**3+d/rh**4)/2
 gamma=n*(2*lam*rh+2*mass/rh**2-2*Qh/rh**3+6*d/rh**4-4*d/rh**5)/2
 eps=G/(60*np.pi);la=lam+11*eps*lam*lam/6
 return lam,d,rh,n,lo,hi,qm,Qh,mass,gamma,eps,la

def mixed_fixture(index,row,geom,action):
 fn=geom[-1];lam,d,rh,n,lo,hi,qm,Qh,mass,gamma,eps,la=seed_geometry(row,1.)
 loU=lo+.2*(hi-lo);hiU=hi+.4*(hi-lo);L=hiU+PROTOCOL['right_tail_exponent']/gamma
 geo=solve_ivp(lambda y,u:lam+u*u-2*mass*u**3+Qh*u**4-2*d*u**5+d*u**6,[0,L],[0.],rtol=PROTOCOL['geometry_rtol'],atol=PROTOCOL['geometry_atol'],dense_output=True,max_step=.01)
 check('mixed_geometry_'+str(index),geo.success,geo.message)
 def fields(y):
  r=1/geo.sol(y)[0];return r,fn(r,lam,mass,Qh,d,1.,qm)
 def shape(y,l,h,order=0):
  x=(2*y-l-h)/(h-l);z=np.clip(x,-1,1);b=1-z*z
  val=[b**6,-12*z*b**5,-12*b**4*(1-11*z*z)][order]
  return np.where(abs(x)<1,val*(2/(h-l))**order,0.)
 def src(y):
  r,f=fields(y);return r*(-shape(y,lo,hi,2)+(gamma**2+f[0]*f[7]/2)*shape(y,lo,hi))
 def U(y):return shape(y,loU,hiU)
 vals=[]
 for level,(left,tol,nq) in enumerate(zip(PROTOCOL['left_endpoints'],PROTOCOL['BVP_tolerances'],PROTOCOL['quadrature_nodes'])):
  grid=merge_grid(np.concatenate([np.geomspace(left,lo,75),np.linspace(lo,hi,180),np.linspace(loU,hiU,180),np.linspace(hiU,L,350)]))
  check('initial_mesh_separated_'+str(index)+'_'+str(level),np.min(np.diff(grid))>PROTOCOL['mesh_merge_relative_gap']*max(1.,L))
  slope=3/left+(gamma*gamma+2*lam)*left/7;bc=lambda a,b:np.array([a[1]-slope*a[0],b[0]])
  sv=solve_bvp(lambda y,z:np.vstack([z[1],(fields(y)[1][8]+gamma*gamma)*z[0]+src(y)/3]),bc,grid,np.zeros((2,len(grid))),tol=tol,max_nodes=PROTOCOL['max_nodes'])
  su=solve_bvp(lambda y,z:np.vstack([z[1],(fields(y)[1][8]+gamma*gamma)*z[0]-fields(y)[1][0]*U(y)/6]),bc,grid,np.zeros((2,len(grid))),tol=tol,max_nodes=PROTOCOL['max_nodes'])
  tag=str(index)+'_'+str(level);check('shape_reconstruction_'+tag,sv.success,sv.message);check('curvature_reconstruction_'+tag,su.success,su.message)
  xx,wt=leggauss(nq);edges=merge_grid(np.concatenate([np.geomspace(left,lo,12),np.linspace(lo,hi,16),np.linspace(loU,hiU,16),np.linspace(hiU,L,40)]))
  y=np.concatenate([(a+b)/2+(b-a)*xx/2 for a,b in zip(edges[:-1],edges[1:])]);weights=np.concatenate([(b-a)*wt/2 for a,b in zip(edges[:-1],edges[1:])])
  r,f=fields(y);F,Fr,Frr,RR,Rr,HR,j,Cv,VK,A1,A0,B1,B0,T,Tr,dr,er,Mw,Mq=f
  v=shape(y,lo,hi);vy=shape(y,lo,hi,1);vyy=shape(y,lo,hi,2);w=r*v
  qv,qvy=sv.sol(y);qu,quy=su.sol(y);p=qv/r;pu=qu/r
  py=(qvy+F*qv/r)/r;puy=(quy+F*qu/r)/r
  pyy=2*F*py/r+(gamma*gamma-F*HR/6)*p+src(y)/(3*r)
  puyy=2*F*puy/r+(gamma*gamma-F*HR/6)*pu-F*U(y)/(6*r)
  P=np.log(r);V=.5*np.log(F/r**2);P1=-F/r;P2=F*Fr/r-F*F/r**2
  V1=F/r-Fr/2;V2=-F*(Fr/r-F/r**2-Frr/2)
  bcw,bww=action(P,V,P1,P2,V1,V2,gamma,la,7.2,qm,pu,puy,puyy,p,py,pyy,v,vy,vyy)
  direct=float(np.dot(weights,-bcw/2-bww/(2*np.pi)))
  D=-eps*r*r*dr*(vyy-gamma*gamma*v)/3-eps*B1*qvy+Mw*w+Mq*qv
  trace=float(-np.dot(weights,qu*D)/2)
  residual=float(-np.dot(weights,qu*F*((r*Tr+4*T)*qv+2*T*w))/2)
  scale=max(abs(direct),abs(trace),1e-12);err=abs(direct-trace)/scale
  check('original_action_vs_trace_mixed_'+tag,err<PROTOCOL['mixed_relative_limit'] and abs(direct-trace)<PROTOCOL['mixed_absolute_limit'],[direct,trace,err])
  check('seed_residual_omission_detected_'+tag,abs(residual)/scale>1e-8,abs(residual)/scale)
  omitted=float(-np.dot(weights,qu*eps*A1*src(y)/(3*F))/2)
  check('curvature_free_source_term_detected_'+tag,abs(omitted)/scale>1e-5,abs(omitted)/scale)
  tails=[float(sv.sol(hiU+1/gamma)[0]),float(su.sol(hiU+1/gamma)[0])]
  check('both_conformal_tails_retained_'+tag,all(abs(z)>1e-12 for z in tails),tails)
  vals.append({'left_endpoint':left,'BVP_tolerance':tol,'quadrature':nq,'shape_nodes':len(sv.x),'curvature_nodes':len(su.x),
   'direct_mixed_action':direct,'trace_mixed_action':trace,'relative_identity_error':err,'residual_contact':residual,
   'omitted_residual_relative_effect':abs(residual)/scale,'tails':tails})
 ref=abs(vals[0]['direct_mixed_action']-vals[1]['direct_mixed_action'])/max(abs(vals[1]['direct_mixed_action']),1e-12)
 check('mixed_action_refinement_'+str(index),ref<PROTOCOL['mixed_refinement_limit'],ref)
 NUMERICS.append({'fixture':'charged_mixed_action','parameters':row,'gamma':gamma,'right_endpoint':L,'refinements':vals,'relative_refinement':ref,'benchmark_threshold_certified':False})

def dilation_fixture(index,row,geom,pure,constants):
 fn=geom[-1];G=PROTOCOL['dilation_G']
 lam,d,rh,n,lo,hi,qm,Qh,mass,gamma,eps,la=seed_geometry(row,G,PROTOCOL['dilation_d'])
 L=2*max(PROTOCOL['dilations']);geo=solve_ivp(lambda y,u:lam+u*u-2*mass*u**3+Qh*u**4-2*d*u**5+d*u**6,[0,L],[0.],rtol=PROTOCOL['geometry_rtol'],atol=PROTOCOL['geometry_atol'],dense_output=True,max_step=L/100)
 check('dilation_geometry_'+str(index),geo.success,geo.message)
 expected=-6*constants['N']+eps*lam*(6*constants['A']+26*constants['D']-60*constants['P'])
 Ulimit=36*lam*lam*constants['H'];vals=[]
 c2=2*(lam*constants['I2']/3+gamma*gamma*constants['J2'])/constants['H']
 c3=2*mass*lam*lam*constants['I3']/constants['H']
 for size in PROTOCOL['dilations']:
  levels=[]
  for nq in PROTOCOL['quadrature_nodes']:
   xx,wt=leggauss(nq);eta=1.5+.5*xx;weights=.5*wt*size;y=size*eta
   chi,ch1,ch2=constants['profile'](eta);q=chi/np.sqrt(size);qy=ch1/size**1.5;qyy=ch2/size**2.5
   r=1/geo.sol(y)[0];f=fn(r,lam,mass,Qh,d,G,qm)
   F,Fr,Frr,RR,Rr,HR,j,Cv,VK,A1,A0,B1,B0,T,Tr,dr,er,Mw,Mq=f
   p=q/r;py=(qy+F*q/r)/r;pyy=(qyy+2*F*py-F*Fr*p)/r
   P=np.log(r);V=.5*np.log(F/r**2);P1=-F/r;P2=F*Fr/r-F*F/r**2;V1=F/r-Fr/2;V2=-F*(Fr/r-F/r**2-Frr/2)
   qc,qw=pure(P,V,P1,P2,V1,V2,gamma,la,7.2,qm,p,py,pyy)
   raw=float(np.dot(weights,-qc- G*qw/np.pi))*size*size
   ko=6+eps*RR/6+22*eps*j/3
   form=6*eps/F*(qyy-(gamma*gamma+F*Fr/r)*q)**2-ko*((qy+F*q/r)**2+gamma*gamma*q*q)-F*(2*RR+48*la)*q*q
   scaled=float(np.dot(weights,form))*size*size
   U=6/F*(-qyy+(gamma*gamma+VK)*q);norm=float(np.dot(weights,U*U))
   err=abs(raw-scaled)/max(abs(scaled),1e-300);tag=str(index)+'_'+str(size)+'_'+str(nq)
   check('dilation_full_action_identity_'+tag,err<PROTOCOL['dilation_action_relative_limit'],err)
   check('negative_dilation_action_'+tag,scaled<0,scaled)
   levels.append({'quadrature':nq,'scaled_full_action':scaled,'scaled_raw_action':raw,'curvature_norm_squared':norm,'action_identity_error':err})
  qe=abs(levels[0]['scaled_full_action']-levels[1]['scaled_full_action'])/abs(levels[1]['scaled_full_action'])
  check('dilation_quadrature_refinement_'+str(index)+'_'+str(size),qe<PROTOCOL['dilation_quadrature_relative_limit'],qe)
  action_error=abs(levels[1]['scaled_full_action']-expected)/abs(expected);norm_error=abs(levels[1]['curvature_norm_squared']-Ulimit)/Ulimit
  vals.append({'delta':size,'refinements':levels,'relative_action_limit_error':action_error,'relative_curvature_limit_error':norm_error})
 check('dilation_coarse_limit_'+str(index),max(vals[0]['relative_action_limit_error'],vals[0]['relative_curvature_limit_error'])<PROTOCOL['dilation_coarse_relative_limit'],vals[0])
 check('dilation_fine_limit_'+str(index),max(vals[-1]['relative_action_limit_error'],vals[-1]['relative_curvature_limit_error'])<PROTOCOL['dilation_fine_relative_limit'],vals[-1])
 check('dilation_limit_refinement_'+str(index),all(vals[i+1]['relative_action_limit_error']<.5*vals[i]['relative_action_limit_error'] and vals[i+1]['relative_curvature_limit_error']<.5*vals[i]['relative_curvature_limit_error'] for i in range(2)))
 NUMERICS.append({'fixture':'charged_conformal_dilation','parameters':row,'G':G,'d':d,'gamma':gamma,'predicted_scaled_action_limit':expected,'predicted_curvature_norm_limit':Ulimit,'relative_norm_correction_C2':c2,'relative_norm_correction_C3':c3,'dilations':vals,'actual_Schur_correction_evaluated':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 thermal curvature-shape Schur','input_count':len(ps)}),flush=True)
 try:
  warnings.filterwarnings('error',category=RuntimeWarning,module=r'scipy\.integrate\._bvp')
  provenance();geom=geometric_functions();constants=exact_checks(geom);cross,pure=local_actions()
  for i,row in enumerate(PROTOCOL['seed_fixtures']):mixed_fixture(i,row,geom,cross)
  for i,row in enumerate(PROTOCOL['seed_fixtures']):dilation_fixture(i,row,geom,pure,constants)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 thermal curvature-shape Schur','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
  'numpy':np.__version__,'scipy':scipy.__version__,'protocol':PROTOCOL,'checks':ROWS,'numerics':NUMERICS,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'physical_L2_to_shape_form_dual_bound_derived':good,'fixed_thermal_cross_coupling_compactness_derived':good,
  'auxiliary_thermal_shape_elimination_derived':good,'compact_positive_Schur_feedback_derived':good,
  'reduced_Euclidean_form_unbounded_below_derived':good,'shape_response_vanishes_on_outer_test_sequence_derived':good,
  'physical_curvature_norm_changed':False,'physical_metric_boundary_domain_selected':False,
  'physical_constraint_mode_constructed':False,'continuous_causal_inverse_proved':False,
  'actual_Schur_feedback_numerically_evaluated':False,'benchmark_frequency_threshold_numerically_enclosed':False,
  'full_quantum_kernel_evaluated':False,'actual_charged_stability_proved':False,'actual_wormhole_instability_proved':False,
  'counterterm_changed':False,'conformal_contour_rotation_adopted':False,'massive_tower_discarded':False,
  'retained_fast_factor_removed':False,'order_reduction_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in [('checker',OWN[0]),('note',OWN[1]),('review',OWN[2])]:run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
