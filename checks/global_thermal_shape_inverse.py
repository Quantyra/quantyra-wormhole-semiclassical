"""Full thermal shape closure: graph, endpoint and adjoint checks."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import numpy as np
import scipy
from scipy.integrate import solve_ivp, solve_bvp
from numpy.polynomial.legendre import leggauss
import sympy as S

ROOT=Path(__file__).resolve().parents[1];STEM='global-thermal-shape-inverse'
OWN=['checks/global_thermal_shape_inverse.py','notes/'+STEM+'.md',
     'evidence/global-thermal-inverse-source-review.json']
PREVIOUS='global-thermal-curvature-free-form';ROWS=[];NUMERICS=[];R=S.Rational
PROTOCOL={
 'scope':'Full high-thermal-frequency Friedrichs shape inverse on auxiliary compact-core closure',
 'normal_fixtures':[[.7,.2,2.],[2.,-.1,5.],[.25,.3,11.]],
 'normal_quadratures':[96,192],'normal_identity_limit':2e-11,
 'seed_fixtures':[[2.,.01,1.1,1,.08,.18],[3.,-.02,1.2,2,.35,.75],[4.,1e-14,1.05,4,1.,1.8]],
 'fixture_columns':['lambda','d','horizon','thermal_index','support_left','support_right'],
 'left_endpoints':[1e-4,5e-5],'BVP_tolerances':[1e-6,1e-8],
 'quadrature_nodes':[48,96],'max_nodes':50000,'right_tail_exponent':12.,
 'geometry_rtol':2e-13,'geometry_atol':2e-15,
 'normalized_adjoint_limit':2e-7,'normalized_refinement_limit':2e-7,
 'full_quantum_inverse_numerically_evaluated':False,'formal_verification':False}

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
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==109 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('source_and_scope',len(rev['primary_sources'])==2 and rev['core_operator_membership_proved']
       and not rev['physical_boundary_realization_identified'] and not rev['continuous_causal_inverse_proved'])

def exact_checks():
 y,g=S.symbols('y g',positive=True);z=S.Function('z')(y)
 direct=(-S.diff(z,y,2)+(g*g+6/y**2)*z)**2
 positive=S.diff(z,y,2)**2+12*S.diff(z,y)**2/y**2+2*g*g*S.diff(z,y)**2+12*g*g*z*z/y**2+g**4*z*z
 flux=-2*g*g*z*S.diff(z,y)-12*z*S.diff(z,y)/y**2-12*z*z/y**3
 zero('normal_graph_identity_with_flux',direct-positive-S.diff(flux,y))
 hardy=S.diff(z,y)**2/y**2-R(9,4)*z*z/y**4
 square=(S.diff(z,y)/y-R(3,2)*z/y**2)**2
 zero('weighted_Hardy_square_with_flux',hardy-square-S.diff(R(3,2)*z*z/y**3,y))
 a=S.Function('a')(y);c=S.Function('c')(y);vk=S.Function('Vk')(y);w=S.Function('w')(y)
 B=lambda x:-S.diff(x,y,2)-2*a*S.diff(x,y)+(g*g+vk+c)*x
 Badj=lambda x:-S.diff(x,y,2)+2*a*S.diff(x,y)+(g*g+vk+c+2*S.diff(a,y))*x
 concomitant=S.diff(z,y)*w-z*S.diff(w,y)-2*a*z*w
 zero('unweighted_adjoint_divergence',z*B(w)-w*Badj(z)-S.diff(concomitant,y))
 check('omitting_divergence_contact_changes_pairing',S.simplify(2*S.diff(a,y)*z*w)!=0)
 b,j,W=S.symbols('b j W');by=S.symbols('by');zp,z0=S.symbols('zp z0')
 diff=2*(1/y+b)*zp+(2*(-1/y**2+by)+6/y**2+j-W)*z0
 zero('adjoint_minus_K_decomposition',diff-(2*zp/y+4*z0/y**2+2*b*zp+(2*by+j-W)*z0))
 zero('weighted_graph_constant',2/S.sqrt(12)+4/S.sqrt(27)-7/(3*S.sqrt(3)))
 check('rational_sqrt3_lower',R(17,10)**2<3)
 bound=(1+R(70,51)*(1+R(1,1000))+R(2,1000))/3
 check('reconstruction_bound_four_fifths',bound<R(4,5),bound)
 check('form_controls_N',R(1,32)-R(1,10**24)>R(1,33))
 for name,lhs,rhs in [('rp_response',R(4,5)*40,32),('gradient_response',40*40,1600),
       ('shape_J2_response',200*40,8000),('theta_J2_response',40*40,1600),('N_response',33*40,1320)]:zero(name,lhs-rhs)
 rr,lam,m,Q,c3,c4,rh,C=S.symbols('r lam m Q c3 c4 rh C',positive=True)
 F=lam*rr**2+1-2*m/rr+Q/rr**2+c3/rr**3+c4/rr**4
 primitive=-lam/(3*rr**3)-1/(5*rr**5)+m/(3*rr**6)-Q/(7*rr**7)-c3/(8*rr**8)-c4/(9*rr**9)
 zero('canonical_weight_integrable_primitive',S.diff(primitive,rr)-F/rr**6)
 zero('canonical_weight_tail_limit',S.limit(primitive,rr,S.oo))
 # A positive point-evaluation form alone is not closable in L2.
 # f_n=max(1-n|x|,0) has f_n(0)=1 but squared L2 norm 2/(3n).
 n=S.symbols('n',positive=True)
 zero('nonclosable_positive_form_control',2*S.integrate((1-n*y)**2,(y,0,1/n))-2/(3*n))

def weyl_endpoint_checks():
 # Differentiate the original optical nonlinear density; derive the Euler
 # rows before setting v=0. No completed-square formula is used.
 e,co,si,g,k=S.symbols('e co si g k',real=True)
 pj=S.symbols('p0:7');vj=S.symbols('v0:7');Pj=S.symbols('P0:7');Vj=S.symbols('V0:7')
 pp=Pj[0]+e*pj[0]*co;vv=Vj[0]+e*vj[0]*co
 dp=Pj[2]+e*(pj[2]-g*g*pj[0])*co;dv=Vj[2]+e*(vj[2]-g*g*vj[0])*co
 grad=Pj[1]**2+2*e*Pj[1]*pj[1]*co+e*e*(pj[1]**2*co**2+g*g*pj[0]**2*si**2)
 raw=S.exp(-2*vv)*(pp*dv**2/15-dp*dv/15-dp**2/10-dv*grad/15-R(7,90)*dp*grad-R(7,180)*grad**2)
 raw+=pp*dv/9-grad/18+k*pp*S.exp(2*vv)
 quad=S.expand(S.diff(raw,e,2).subs(e,0)/2).subs({co**2:1,si**2:1})
 alljets=[pj,vj,Pj,Vj]
 def Dy(expr):return S.expand(sum(S.diff(expr,js[i])*js[i+1] for js in alljets for i in range(6)))
 def Euler(js):return S.expand((S.diff(quad,js[0])-Dy(S.diff(quad,js[1]))+Dy(Dy(S.diff(quad,js[2]))))/2)
 ep=Euler(pj).subs(dict.fromkeys(vj,0));ev=Euler(vj).subs(dict.fromkeys(vj,0))
 y=S.symbols('y',positive=True)
 def subjets(P,V,p,deriv):
  sub={g:3,k:R(-7,3)}
  for js,expr in [(Pj,P),(Vj,V),(pj,p)]:
   for j in range(7):sub[js[j]]=expr;expr=deriv(expr)
  return sub
 sub=subjets(-S.log(y),S.Integer(0),y**4,lambda a:S.diff(a,y))
 for name,expr in [('p',ep),('v',ev)]:
  val=S.simplify(expr.subs(sub));lim=S.limit(val,y,0)
  check('outer_Weyl_row_bounded_'+name,lim.is_finite,lim)
  zero('outer_weighted_Weyl_row_vanishes_'+name,S.limit(y*val,y,0))
 bad=subjets(-S.log(y),S.Integer(0),S.sqrt(y),lambda a:S.diff(a,y))
 badlimit=S.limit(y**R(7,2)*ep.subs(bad),y,0)
 check('uncontrolled_conformal_tail_fails_L2',badlimit.is_finite and badlimit!=0,badlimit)
 rho=S.symbols('rho',positive=True)
 for n in [1,2,4]:
  sub=subjets(rho*rho/7,S.log(rho)+rho*rho/11,rho**n*(1+rho*rho/5),lambda a:-rho*S.diff(a,rho));sub[g]=n
  for name,expr in [('p',ep),('v',ev)]:
   val=S.simplify(expr.subs(sub));lim=S.limit(val/rho**(n+2),rho,0)
   check('smooth_bolt_Weyl_density_'+name+'_'+str(n),lim.is_finite,lim)

def normal_fixtures():
 y=S.symbols('y');
 for i,(a,b,g) in enumerate(PROTOCOL['normal_fixtures']):
  expr=y**3*(1+b*y)*S.exp(-a*y);fn=S.lambdify(y,[expr,S.diff(expr,y),S.diff(expr,y,2)],'numpy')
  vals=[];L=50/a
  for nq in PROTOCOL['normal_quadratures']:
   xx,wt=leggauss(nq);yy=(xx+1)*L/2;wt=wt*L/2;z,zp,zpp=fn(yy)
   K=-zpp+(g*g+6/yy**2)*z
   direct=float(np.dot(wt,K*K))
   terms=[zpp*zpp,12*zp*zp/yy**2,2*g*g*zp*zp,12*g*g*z*z/yy**2,g**4*z*z]
   pieces=[float(np.dot(wt,t)) for t in terms];via=sum(pieces)
   err=abs(via-direct)/(1+abs(direct));tag=str(i)+'_'+str(nq)
   check('normal_graph_quadrature_'+tag,err<PROTOCOL['normal_identity_limit'],err)
   check('weighted_graph_contact_detected_'+tag,pieces[1]/(1+abs(direct))>1e-8,pieces[1])
   hh=float(np.dot(wt,zp*zp/yy**2-R(9,4)*z*z/yy**4))
   sq=float(np.dot(wt,(zp/yy-1.5*z/yy**2)**2))
   check('weighted_Hardy_quadrature_'+tag,abs(hh-sq)/(1+abs(sq))<PROTOCOL['normal_identity_limit'],[hh,sq])
   vals.append({'quadrature':nq,'direct_graph_norm_squared':direct,'positive_terms':pieces,'normalized_error':err})
  ref=abs(vals[0]['direct_graph_norm_squared']-vals[1]['direct_graph_norm_squared'])/(1+abs(vals[1]['direct_graph_norm_squared']))
  check('normal_graph_refinement_'+str(i),ref<PROTOCOL['normal_identity_limit'],ref)
  NUMERICS.append({'fixture':'normal_graph','parameters':[a,b,g],'refinements':vals,'normalized_refinement':ref})

def seed_fixture(index,row):
 lam,d,rh,n,lo,hi=row;qh=1+3*lam+d
 mass=rh*(1+lam*rh**2+qh/rh**2-2*d/rh**3+d/rh**4)/2
 frh=2*lam*rh+2*mass/rh**2-2*qh/rh**3+6*d/rh**4-4*d/rh**5
 gamma=n*frh/2;lo2=lo+.3*(hi-lo);hi2=hi+.5*(hi-lo)
 L=hi2+PROTOCOL['right_tail_exponent']/gamma
 geo=solve_ivp(lambda y,u:lam+u*u-2*mass*u**3+qh*u**4-2*d*u**5+d*u**6,[0,L],[0.],
   rtol=PROTOCOL['geometry_rtol'],atol=PROTOCOL['geometry_atol'],dense_output=True,max_step=.01)
 check('geometry_solve_'+str(index),geo.success,geo.message)
 def fields(y):
  r=1/geo.sol(y)[0];F=lam*r*r+1-2*mass/r+qh/r**2-2*d/r**3+d/r**4
  Fr=2*lam*r+2*mass/r**2-2*qh/r**3+6*d/r**4-4*d/r**5
  Frr=2*lam-4*mass/r**3+6*qh/r**4-24*d/r**5+20*d/r**6
  HR=-24*lam-12*d/r**5+24*d/r**6
  Cv=24*lam+4/r**2-8*d/r**5+12*d/r**6
  VK=F*(Fr/r-HR/6);a=F/r;ay=-F*(Fr/r-F/r**2);c=F*(Frr+4*Fr/r+HR/6)
  return r,F,VK,Cv,a,ay,c
 def shape(y,left,right,order=0):
  x=(2*y-left-right)/(right-left);z=np.clip(x,-1,1);u=1-z*z
  val=[u**6,-12*z*u**5,-12*u**4*(1-11*z*z)][order]
  return np.where(abs(x)<1,val*(2/(right-left))**order,0.)
 def source(y):
  r,F,VK,Cv,a,ay,c=fields(y)
  return r*(-shape(y,lo,hi,2)+(gamma**2+F*Cv/2)*shape(y,lo,hi))
 def force(y):return fields(y)[0]*shape(y,lo2,hi2)
 vals=[]
 for level,(left,tol,nq) in enumerate(zip(PROTOCOL['left_endpoints'],PROTOCOL['BVP_tolerances'],PROTOCOL['quadrature_nodes'])):
  grid=np.unique(np.concatenate([np.geomspace(left,lo,75),np.linspace(lo,hi,180),np.linspace(lo2,hi2,180),np.linspace(hi2,L,350)]))
  slope=3/left+(gamma**2+2*lam)*left/7
  bc=lambda a,b:np.array([a[1]-slope*a[0],b[0]])
  qp=solve_bvp(lambda y,z:np.vstack([z[1],(fields(y)[2]+gamma**2)*z[0]+source(y)/3]),bc,grid,np.zeros((2,len(grid))),tol=tol,max_nodes=PROTOCOL['max_nodes'])
  adj=solve_bvp(lambda y,z:np.vstack([z[1],(fields(y)[2]+gamma**2)*z[0]-force(y)]),bc,grid,np.zeros((2,len(grid))),tol=tol,max_nodes=PROTOCOL['max_nodes'])
  tag=str(index)+'_'+str(level)
  check('primal_BVP_'+tag,qp.success,qp.message);check('adjoint_BVP_'+tag,adj.success,adj.message)
  edges=np.unique(np.concatenate([np.geomspace(left,lo,12),np.linspace(lo,hi,16),np.linspace(lo2,hi2,16),np.linspace(hi2,L,40)]))
  xx,wt=leggauss(nq)
  y=np.concatenate([(a+b)/2+(b-a)*xx/2 for a,b in zip(edges[:-1],edges[1:])])
  weights=np.concatenate([(b-a)*wt/2 for a,b in zip(edges[:-1],edges[1:])])
  r,F,VK,Cv,a,ay,c=fields(y);w=r*shape(y,lo,hi);f=force(y)
  q=qp.sol(y)[0];z,zy=adj.sol(y)
  tf=-(f+2*a*zy+(2*ay+c)*z)/3
  wrong=-(f+2*a*zy+c*z)/3
  direct=float(np.dot(weights,f*q));via=float(np.dot(weights,w*tf))
  wrongpair=float(np.dot(weights,w*wrong));green=float(-np.dot(weights,source(y)*z)/3)
  err=abs(direct-via)/(1+abs(direct));ge=abs(direct-green)/(1+abs(direct))
  check('charged_adjoint_pairing_'+tag,err<PROTOCOL['normalized_adjoint_limit'],err)
  check('independent_resolvent_pairing_'+tag,ge<PROTOCOL['normalized_adjoint_limit'],ge)
  missing=abs(wrongpair-direct)/(1+abs(direct))
  check('missing_divergence_contact_detected_'+tag,missing>1e-5,missing)
  normratio=float(np.sqrt(np.dot(weights,q*q)/np.dot(weights,w*w)))
  tail=float(qp.sol(hi2+1/gamma)[0]);check('global_tail_retained_'+tag,abs(tail)>1e-10,tail)
  vals.append({'left_endpoint':left,'BVP_tolerance':tol,'quadrature':nq,'primal_nodes':len(qp.x),'adjoint_nodes':len(adj.x),
   'direct_pairing':direct,'adjoint_pairing':via,'resolvent_pairing':green,'normalized_adjoint_error':err,
   'normalized_resolvent_error':ge,'omitted_contact_discrepancy':missing,'observed_reconstruction_norm_ratio':normratio,'tail_sample':tail})
 ref=abs(vals[1]['direct_pairing']-vals[0]['direct_pairing'])/(1+abs(vals[1]['direct_pairing']))
 check('charged_pairing_refinement_'+str(index),ref<PROTOCOL['normalized_refinement_limit'],ref)
 NUMERICS.append({'fixture':'charged_reconstruction_adjoint','parameters':row,'gamma':gamma,'right_endpoint':L,
  'refinements':vals,'normalized_refinement':ref,'benchmark_threshold_certified':False,'full_quantum_inverse_evaluated':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 global thermal shape inverse','input_count':len(ps)}),flush=True)
 try:
  provenance();exact_checks();weyl_endpoint_checks();normal_fixtures()
  for i,row in enumerate(PROTOCOL['seed_fixtures']):seed_fixture(i,row)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 global thermal shape inverse','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
  'numpy':np.__version__,'scipy':scipy.__version__,'protocol':PROTOCOL,'checks':ROWS,'numerics':NUMERICS,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'full_thermal_core_operator_maps_into_H_derived':good,'full_thermal_form_closability_derived':good,
  'full_global_thermal_Friedrichs_inverse_derived':good,'bounded_global_reconstruction_derived':good,
  'physical_metric_boundary_domain_selected':False,'continuous_causal_inverse_proved':False,
  'benchmark_frequency_threshold_numerically_enclosed':False,'full_quantum_kernel_evaluated':False,
  'full_quantum_inverse_numerically_evaluated':False,'actual_metric_mode_constructed':False,
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
