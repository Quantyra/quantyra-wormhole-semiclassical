"""Global classical curvature-free form: identities, comparison and fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import numpy as np
import scipy
from scipy.integrate import solve_ivp, solve_bvp
from numpy.polynomial.legendre import leggauss
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM='global-classical-curvature-free-form'
OWN=['checks/global_classical_curvature_free_form.py','notes/'+STEM+'.md',
     'evidence/global-classical-shape-source-review.json']
PREVIOUS='first-thermal-horizon-flux';ROWS=[];NUMERICS=[]
PROTOCOL={
 'scope':'Global large-real-frequency classical component, not the full quantum shape inverse or physical stability',
 'seed_fixtures':[[2.,.01,1.1,.02,7.2,4.,.08,.18],
                  [3.,-.02,1.2,.01,10.3,12.,.35,.75],
                  [4.,1e-14,1.05,1e-10,13.,32.,1.,1.8]],
 'fixture_columns':['lambda','d','horizon','lambda0_minus_lambda','physical_magnetic_length_squared','gamma','support_left','support_right'],
 'left_endpoints':[1e-4,5e-5], 'BVP_tolerances':[1e-6,1e-8],
 'quadrature_nodes':[48,96], 'max_nodes':50000,
 'geometry_rtol':2e-13,'geometry_atol':2e-15,
 'right_tail_exponent':32., 'normalized_identity_limit':2e-7,
 'normalized_refinement_limit':2e-7,
 'shape':'(1-x^2)^6 for |x|<1, with factored derivatives; piecewise C5',
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
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==65 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('source_and_scope',len(rev['primary_sources'])==1 and rev['physical_charge_kept_separate'] and not rev['full_global_shape_inverse_proved'])
 check('development_failures_retained',len(rev['development_runs'])==2 and all(r['checker_source_utf8'] for r in rev['development_runs']))

def exact_checks():
 r,lam,m,Qh,c3,c4,Q=S.symbols('r lam m Qh c3 c4 Q',positive=True)
 F=lam*r*r+1-2*m/r+Qh/r**2+c3/r**3+c4/r**4
 R=-S.diff(F,r,2)-4*S.diff(F,r)/r+2*(1-F)/r**2
 HR=2*R+r*S.diff(R,r)
 Cv=2*S.diff(F,r,2)+4*F/r**2+8*S.diff(F,r)/r
 VK=F*(S.diff(F,r)/r-HR/6);a=F/r
 c=F*(S.diff(F,r,2)+4*S.diff(F,r)/r+HR/6)
 zero('seed_scalar_curvature',R+12*lam+2*c3/r**5+6*c4/r**6)
 zero('shape_curvature_coefficient',Cv-24*lam-4/r**2-4*c3/r**5-12*c4/r**6)
 w=S.Function('w')(r);g=S.symbols('g',positive=True);dy=lambda e:-F*S.diff(e,r)
 v=w/r
 Bv=-dy(dy(v))+g*g*v+F*Cv*v/2
 Bt=-dy(dy(w))-2*a*dy(w)+(g*g+VK+c)*w
 zero('unitary_shape_source',r*Bv-Bt)
 Sp=-dy(dy(v))+2*a*dy(v)+(g*g-F*HR/6)*v
 zero('unitary_reconstruction',r*Sp+dy(dy(w))-(g*g+VK)*w)
 zero('weighted_gradient_potential',a*a-dy(a)-F*S.diff(F,r)/r)
 z=S.symbols('z',positive=True)
 ys=z/lam-z**3/(3*lam**2)+m*z**4/(2*lam**2)
 for name,expr,normal,target in [
  ('VK',VK,6/ys**2,2*lam),('c',c,6/ys**2,2*lam),
  ('P',F*(1/r**2+3*lam-Q/r**4),3/ys**2,2*lam),
  ('Wg',F*S.diff(F,r)/r,2/ys**2,2*lam/3),
  ('b',a,1/ys,0),('b_derivative',dy(a),-1/ys**2,2*lam/3)]:
  zero('outer_bounded_coefficient_'+name,S.limit(expr.subs(r,1/z)-normal,z,0)-target)
 RR,Rr,F0,rr,la0,p,vv,kin=S.symbols('RR Rr F0 rr la0 p vv kin')
 Tc=RR+12*la0
 qpp=-(6*kin+F0*rr**2*(2*RR+48*la0)*p*p)/2
 contact=-F0*rr**2*((rr*Rr+4*Tc)*p*p+2*Tc*p*vv)/2
 transformed=-qpp+2*contact
 target=3*kin-F0*rr**2*(RR+rr*Rr/2)*p*p-F0*rr**2*((2*Tc+rr*Rr/2)*p*p+2*Tc*p*vv)
 zero('classical_residual_congruence',transformed-target)
 check('dropping_residual_changes_form',S.simplify(contact)!=0)
 d,B=S.symbols('d B',positive=True)
 deltaR=S.simplify((R+12*lam).subs({c3:-2*d*B**3,c4:d*B**4,r:B/z}))
 zero('seed_small_curvature',deltaR-d*(4*z**5-6*z**6)/B**2)
 rrder=S.simplify((r*S.diff(R,r)).subs({c3:-2*d*B**3,c4:d*B**4,r:B/z}))
 zero('seed_small_radial_curvature',rrder-d*(-20*z**5+36*z**6)/B**2)
 e=S.Rational(1,1000);sig=S.Rational(1,10**9)
 margin=(S.Rational(1,9)-16*e)/3-sig*(S.Rational(10,27)+8)
 check('coercivity_margin_above_one_over_32',margin>S.Rational(1,32),margin)
 check('physical_gradient_coercivity_one_over_40',margin/(1+e)>S.Rational(1,40))
 check('energy_upper_below_ten',(3+e)**2/(1-e)<10)
 check('lapse_weight_bound',(S.Rational(8,3)+e)/3<1)
 gmax=(S.Rational(4,10**17))**2
 shift=11*gmax/(360*3)
 check('cosmological_shift_budget',shift<sig)
 check('classical_trace_residual_budget',12*shift+10*S.Rational(1,10**12)<sig)
 check('classical_density_residual_budget',24*shift+48*S.Rational(1,10**12)<sig)
 # Verify the direct optical action's second derivative independently.
 eps,P,W,pp,vp,Py,Wy,ppy,vpy,freq,lp,Qp=S.symbols('eps P W pp vp Py Wy ppy vpy freq lp Qp')
 raw=S.exp(2*(P+eps*pp))*(6*(Py+eps*ppy)**2+4*(Py+eps*ppy)*(Wy+eps*vpy)
     +eps**2*freq**2*(6*pp**2+4*pp*vp))
 raw+=S.exp(2*(W+eps*vp))*(2*S.exp(2*(P+eps*pp))+6*lp*S.exp(4*(P+eps*pp))-2*Qp)
 explicit=S.exp(2*P)*(4*pp**2*(6*Py**2+4*Py*Wy)+4*pp*(12*Py*ppy+4*Py*vpy+4*ppy*Wy)
     +12*ppy**2+8*ppy*vpy+12*freq**2*pp**2+8*freq**2*pp*vp)
 explicit+=S.exp(2*W)*(2*S.exp(2*P)*(2*pp+2*vp)**2+6*lp*S.exp(4*P)*(4*pp+2*vp)**2-2*Qp*(2*vp)**2)
 zero('direct_action_second_derivative',S.diff(raw,eps,2).subs(eps,0)-explicit)
 zero('physical_charge_second_derivative',S.diff(explicit,Qp)+8*S.exp(2*W)*vp**2)

def matrix_fixtures():
 base=np.array([[2.,.2+.1j,0],[.2-.1j,3.,.3],[0,.3,4.]])
 perturb=np.array([[.1,.05j,.03],[-.05j,-.08,.04],[.03,.04,.06]])
 f=np.array([1.+.1j,-.2+.3j,.6]);df=np.array([.03,-.04j,.02])
 for i,gamma in enumerate([1.,3.,10.]):
  K0=base+gamma**2*np.eye(3);K=K0+perturb
  eps=np.linalg.norm(perturb,2)/gamma**2
  E0=float(np.vdot(f,np.linalg.solve(K0,f)).real)
  dn=float(np.sqrt(np.vdot(df,np.linalg.solve(K0,df)).real))
  actual=float(np.vdot(f+df,np.linalg.solve(K,f+df)).real)
  lower=(E0-2*np.sqrt(E0)*dn)/(1+eps);upper=(np.sqrt(E0)+dn)**2/(1-eps)
  check('noncommuting_variational_lower_'+str(i),actual>=lower-1e-13,[lower,actual])
  check('noncommuting_variational_upper_'+str(i),actual<=upper+1e-13,[actual,upper])
  check('noncommuting_control_'+str(i),np.linalg.norm(K0@perturb-perturb@K0)>.01)
  NUMERICS.append({'fixture':'noncommuting_form_comparison','gamma':gamma,'lower':lower,'actual':actual,'upper':upper})

def seed_fixture(index,row):
 lam,d,rh,dl,Q,gamma,lo,hi=row
 qh=1+3*lam+d
 mass=rh*(1+lam*rh**2+qh/rh**2-2*d/rh**3+d/rh**4)/2
 L=hi+PROTOCOL['right_tail_exponent']/gamma
 fun=lambda y,u:lam+u*u-2*mass*u**3+qh*u**4-2*d*u**5+d*u**6
 geo=solve_ivp(fun,[0,L],[0.],rtol=PROTOCOL['geometry_rtol'],atol=PROTOCOL['geometry_atol'],dense_output=True,max_step=.01)
 check('geometry_solve_'+str(index),geo.success,geo.message)
 def fields(y):
  u=geo.sol(y)[0];r=1/u
  F=lam*r*r+1-2*mass/r+qh/r**2-2*d/r**3+d/r**4
  Fr=2*lam*r+2*mass/r**2-2*qh/r**3+6*d/r**4-4*d/r**5
  R=-12*lam+4*d/r**5-6*d/r**6;Rr=-20*d/r**6+36*d/r**7
  HR=2*R+r*Rr;Cv=24*lam+4/r**2-8*d/r**5+12*d/r**6
  return r,F,Fr,R,Rr,Cv,F*(Fr/r-HR/6)
 def shape(y,order=0):
  x=(2*y-lo-hi)/(hi-lo);z=np.clip(x,-1,1);u=1-z*z
  val=[u**6,-12*z*u**5,-12*u**4*(1-11*z*z)][order]
  return np.where(abs(x)<1,val*(2/(hi-lo))**order,0.)
 def src(y):
  r,F,Fr,R,Rr,Cv,VK=fields(y)
  return r*(-shape(y,2)+(gamma**2+F*Cv/2)*shape(y))
 values=[]
 for level,(left,tol,nquad) in enumerate(zip(PROTOCOL['left_endpoints'],PROTOCOL['BVP_tolerances'],PROTOCOL['quadrature_nodes'])):
  grid=np.unique(np.concatenate([np.geomspace(left,lo,75),np.linspace(lo,hi,180),np.linspace(hi,L,350)]))
  def ode(y,state):return np.vstack([state[1],(fields(y)[-1]+gamma**2)*state[0]+src(y)/3])
  slope=3/left+(gamma**2+2*lam)*left/7
  sol=solve_bvp(ode,lambda a,b:np.array([a[1]-slope*a[0],b[0]]),grid,np.zeros((2,len(grid))),tol=tol,max_nodes=PROTOCOL['max_nodes'])
  tag=str(index)+'_'+str(level)
  check('reconstruction_BVP_'+tag,sol.success,sol.message)
  zz,ww=leggauss(nquad)
  edges=np.unique(np.concatenate([np.geomspace(left,lo,12),np.linspace(lo,hi,16),np.linspace(hi,L,40)]))
  y=np.concatenate([(a+b)/2+(b-a)*zz/2 for a,b in zip(edges[:-1],edges[1:])])
  weights=np.concatenate([(b-a)*ww/2 for a,b in zip(edges[:-1],edges[1:])])
  r,F,Fr,R,Rr,Cv,VK=fields(y);v=shape(y);vy=shape(y,1)
  q,qp=sol.sol(y);p=q/r;py=(qp+F*q/r)/r;w=r*v
  Py=-F/r;Wy=F/r-Fr/2
  kinetic=r*r*(4*p*p*(6*Py*Py+4*Py*Wy)+4*p*(12*Py*py+4*Py*vy+4*py*Wy)+12*py*py+8*py*vy+12*gamma**2*p*p+8*gamma**2*p*v)
  potential=F/(r*r)*(2*r*r*(2*p+2*v)**2+6*(lam+dl)*r**4*(4*p+2*v)**2-2*Q*(2*v)**2)
  direct=-np.dot(weights,kinetic+potential)/4
  E=-3*np.dot(weights,src(y)*q)
  Egradient=9*np.dot(weights,qp*qp+(VK+gamma**2)*q*q)
  Tc=R+12*(lam+dl)
  contact=np.dot(weights,F*((2*Tc+r*Rr/2)*q*q+2*Tc*q*w))
  via=E/3-2*np.dot(weights,F*(1/r**2+3*(lam+dl)-Q/r**4)*w*w)-contact
  norm=np.dot(weights,r*r*(vy*vy+gamma**2*v*v))
  identity=abs(direct-via)/(1+abs(direct));energy=abs(E-Egradient)/(1+abs(E))
  check('direct_action_vs_congruence_'+tag,identity<PROTOCOL['normalized_identity_limit'],identity)
  check('reconstruction_energy_identity_'+tag,energy<PROTOCOL['normalized_identity_limit'],energy)
  if index<2:check('residual_omission_detected_'+tag,abs(contact)/(1+abs(direct))>1e-5,contact)
  qtail=float(sol.sol(hi+1/gamma)[0])
  check('global_conformal_tail_retained_'+tag,abs(qtail)>1e-9,qtail)
  tail_estimate=float(np.dot(weights,abs(src(y)))*np.exp(-32)/(6*gamma*np.sqrt(2*gamma)))
  values.append({'left_endpoint':left,'tolerance':tol,'quadrature_nodes':nquad,'BVP_nodes':len(sol.x),
   'direct_action':float(direct),'congruence_action':float(via),'source_energy':float(E),'gradient_energy':float(Egradient),
   'residual_contact':float(contact),'weighted_gradient_norm':float(norm),'action_to_norm_ratio':float(direct/norm),
   'normalized_action_error':float(identity),'normalized_energy_error':float(energy),
   'conformal_tail_sample':qtail,'quadrature_estimate_of_analytic_tail_bound':tail_estimate})
 refinement=abs(values[1]['direct_action']-values[0]['direct_action'])/(1+abs(values[1]['direct_action']))
 check('action_refinement_'+str(index),refinement<PROTOCOL['normalized_refinement_limit'],refinement)
 NUMERICS.append({'fixture':'charged_classical_action','parameters':row,'comparison_charge_coefficient':qh,
  'mass':mass,'right_endpoint':L,'refinements':values,'normalized_refinement':refinement,
  'threshold_certified_for_fixture':False,'physical_kernel_evaluated':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 global classical curvature-free form','input_count':len(ps)}),flush=True)
 try:
  provenance();exact_checks();matrix_fixtures()
  for index,row in enumerate(PROTOCOL['seed_fixtures']):seed_fixture(index,row)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 global classical curvature-free form','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
  'numpy':np.__version__,'scipy':scipy.__version__,'protocol':PROTOCOL,'checks':ROWS,'numerics':NUMERICS,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'global_classical_curvature_free_coercivity_derived':good,'global_classical_variational_inverse_derived':good,
  'classical_off_shell_contacts_retained':good,'physical_magnetic_charge_distinction_retained':good,
  'benchmark_frequency_threshold_numerically_enclosed':False,'full_global_shape_inverse_proved':False,
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
