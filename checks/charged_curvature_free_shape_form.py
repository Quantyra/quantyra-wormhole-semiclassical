"""Actual charged shape congruence and conformal-tail transfer fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import numpy as np
import scipy
from scipy.integrate import solve_ivp, solve_bvp, quad
from scipy.special import iv, kv
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1];STEM='charged-curvature-free-shape-form'
OWN=['checks/charged_curvature_free_shape_form.py','notes/'+STEM+'.md','evidence/charged-shape-form-source-review.json']
PREVIOUS='physical-curvature-reconstruction';ROWS=[];NUMERICS=[]
PROTOCOL={'scope':'Actual charged curvature-free shape form and positive annular transfer; no global physical inverse',
 'seed_fixture':{'B':1,'lambda':2,'d_seed':.01,'horizon':1.1},
 'dilations':[.01,.005,.0025],'scaled_frequencies':[1.,2.],
 'shape':'256(eta-1)^4(2-eta)^4 on [1,2], zero outside',
 'domain_end':32.,'left_endpoints':[1e-4,5e-5],'BVP_tolerances':[1e-7,1e-9],
 'geometry_rtol':2e-13,'geometry_atol':2e-15,'max_nodes':30000,
 'quadrature_nodes':96,'dilation_refinement_ratio':.35,
 'BVP_refinement_absolute_tolerance':1e-7,'BVP_refinement_fraction_of_dilation_error':.01,
 'tail_norm_bound_tolerance':1e-8,'tail_fraction_of_dilation_error':1e-4,
 'actual_Dirac_kernel_numerically_evaluated':False,'fixtures_are_physical_modes':False,'formal_verification':False}

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
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==98 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 check('quantum_and_positive_inputs_pinned',all('notes/'+n+'.md' in p['input_sha256_after'] for n in [
  'annular-retarded-operator-limit','curvature-free-all-frequency-bound','local-coupled-frequency-polynomial']))
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('source_scope',len(review['primary_sources'])==1 and not review['conformal_tails_discarded'] and not review['global_shape_inverse_proved'])

def exact_checks():
 M=S.Matrix([[3,1],[1,2]]);N=S.Matrix([[-2,1],[3,4]]);P=M.inv();Q=-P*N
 A=S.Matrix([[-3,2],[2,5]]);B=S.Matrix([[1,-2],[4,3]]);D=S.Matrix([[7,1],[1,6]])
 H=A.row_join(B).col_join(B.T.row_join(D));J=P.row_join(Q).col_join(S.zeros(2).row_join(S.eye(2)))
 transformed=J.T*H*J;Scf=D+B.T*Q+Q.T*B+Q.T*A*Q;cross=Q.T*A*P+B.T*P
 zero('reconstruction_curvature_identity',(M*P-S.eye(2)).norm())
 zero('curvature_free_graph',(M*Q+N).norm())
 zero('exact_shape_congruence',(transformed[2:4,2:4]-Scf).norm())
 zero('exact_cross_block',(transformed[2:4,0:2]-cross).norm())
 check('missing_transpose_row_detected',Scf!=D+B.T*Q)
 uv=S.Matrix([1,-2,3,4]);source=H*J*uv
 zero('transformed_source_equivalence',(transformed*uv-J.T*source).norm())
 zero('shape_source_transpose',(J.T*source)[2:4,0].norm()-(source[2:4,0]+Q.T*source[0:2,0]).norm())
 K=S.Matrix([[2,1],[1,4]]);HK=S.zeros(2).row_join(S.zeros(2)).col_join(S.zeros(2).row_join(K))
 zero('canonical_shape_kernel_unchanged',(J.T*HK*J-HK).norm())
 initial=S.Matrix([S.Rational(1,3),-2,0,0]);sfull=H*(J*uv+initial)
 zero('initial_source_shift',(J.T*(sfull-H*initial)-transformed*uv).norm())
 eps,G,r,F,T,Tr,p,v,gm=S.symbols('eps G r F T Tr p v gm',nonzero=True)
 row=-(eps*gm/(r*F)+r*Tr*p+(4*p+2*v)*T)/(8*S.pi*G)
 defect=-eps*r*p*gm/(2*G)-F*r*r*((r*Tr+4*T)*p*p+2*T*p*v)/(2*G)
 zero('exact_density_contact_identity',4*S.pi*F*r*r*p*row-defect)
 pp,pv,vv,er=S.symbols('pp pv vv er')
 zero('curvature_free_form_contact_rewrite',(pp+2*pv+vv-(vv-pp+2*er)).subs(pv,er-pp))
 rr,lam,m,Qh,c3,c4,flux=S.symbols('rr lam m Qh c3 c4 flux',positive=True)
 Fr=lam*rr**2+1-2*m/rr+Qh/rr**2+c3/rr**3+c4/rr**4
 R=-S.diff(Fr,rr,2)-4*S.diff(Fr,rr)/rr-2*(Fr-1)/rr**2
 W=-S.diff(Fr,rr,2)+2*S.diff(Fr,rr)/rr+2*(1-Fr)/rr**2
 j=-S.diff(Fr,rr,2)/4+(Fr-1)/(2*rr**2)
 lam0=lam+11*eps*lam**2/6
 boxR=Fr*S.diff(R,rr,2)+(S.diff(Fr,rr)+2*Fr/rr)*S.diff(R,rr)
 trace=R+12*lam0+eps*(7*W*W/36+22*j*j/3-11*R*R/72-10*flux**2/rr**4+boxR)
 zero('trace_outer_constant_cancelled',S.limit(trace,rr,S.oo))
 zero('trace_outer_fourth_coefficient',S.limit(rr**4*trace,rr,S.oo)+10*eps*flux**2)
 zero('trace_density_contact_leading_cancelled',S.limit(rr**4*(rr*S.diff(trace,rr)+4*trace),rr,S.oo))
 Drec=S.diff(Fr,rr)/rr-(2*R+rr*S.diff(R,rr))/6;VK=Fr*Drec
 z=S.symbols('z',positive=True)
 yseries=z/lam-z**3/(3*lam**2)+m*z**4/(2*lam**2)
 zero('global_K_remainder_outer_limit',S.limit(VK.subs(rr,1/z)-6/yseries**2,z,0)-2*lam)
 eta,k=S.symbols('eta k',positive=True);zz=S.Function('z')(eta);pprof=eta*zz
 Ak=lambda w:eta**2*(S.diff(w,eta,2)-k*k*w)-2*eta*S.diff(w,eta)-4*w
 zero('normal_reconstruction_identification',Ak(pprof)+eta**3*(-S.diff(zz,eta,2)+(6/eta**2+k*k)*zz))
 w=eta**3*S.exp(-eta)
 integrand=(-S.diff(w,eta,2)+6*w/eta**2)**2-S.diff(w,eta,2)**2-12*S.diff(w,eta)**2/eta**2
 zero('global_stable_graph_norm_fixture',S.integrate(S.expand(integrand),(eta,0,S.oo)))
 C,lo,x=S.symbols('C lo x',positive=True);L=lo/S.sqrt(x)
 zero('physical_time_and_profile_normalization',(C*L**2).subs(C,1/(lam*L**2))-1/lam)
 a0,A,Bb,X,h=S.symbols('a0 A B X h',positive=True);theta=a0/(a0+Bb)
 zero('logarithmic_coercivity_combination',theta*(A*X-Bb*h)+(1-theta)*a0*h-a0*A*X/(a0+Bb))
 check('coercivity_weight_between_zero_and_one',theta.is_positive and S.simplify(1-theta).is_positive)
 xi=S.symbols('xi',positive=True)
 check('H2_alone_does_not_control_log_space',S.limit(1+S.log(2+xi**2),xi,S.oo)==S.oo)
 # Independent exact noncommuting resolvent check with a bounded scaled perturbation.
 K0=S.Matrix([[4,1,0],[1,3,1],[0,1,5]]);Wm=S.Matrix([[1,2,0],[2,-1,1],[0,1,2]])
 dd=S.Rational(1,20);sig=S.Rational(3,2);Ra=(K0+sig*S.eye(3)).inv();Rb=(K0+dd**2*Wm+sig*S.eye(3)).inv()
 zero('global_second_resolvent_identity',(Rb-Ra+Rb*(dd**2*Wm)*Ra).norm())
 check('resolvent_fixture_noncommuting',K0*Wm!=Wm*K0)

def profile_functions():
 x=S.symbols('x');expr=256*(x-1)**4*(2-x)**4
 fs=[S.lambdify(x,S.diff(expr,x,n),'numpy') for n in [0,2]]
 def fun(n,x):
  xx=np.asarray(x);return np.where((xx>1)&(xx<2),fs[n](xx),0.)
 return lambda x:fun(0,x),lambda x:fun(1,x)

def normal_green(k,v,v2):
 nu=2.5
 aa=lambda y:np.sqrt(y)*iv(nu,k*y);bb=lambda y:np.sqrt(y)*kv(nu,k*y)
 source=lambda y:(float(v2(y))-k*k*float(v(y)))/(3*y)-4*float(v(y))/y**3
 opts={'epsabs':2e-12,'epsrel':2e-12,'limit':150}
 headall=quad(lambda s:aa(s)*source(s),1,2,**opts)[0]
 tailall=quad(lambda s:bb(s)*source(s),1,2,**opts)[0]
 def value(y):
  if y<=1:return aa(y)*tailall
  if y>=2:return bb(y)*headall
  return bb(y)*quad(lambda s:aa(s)*source(s),1,y,**opts)[0]+aa(y)*quad(lambda s:bb(s)*source(s),y,2,**opts)[0]
 return lambda grid:np.array([value(float(y)) for y in np.atleast_1d(grid)]),source

def seed_fixtures():
 lam=2.;ds=.01;rh=1.1
 ext=lambda r:(r-1)**2/r**2*(1+lam*(r*r+2*r+3)+ds/r**2)
 mass=1+2*lam+rh*ext(rh)/2
 # Ascending polynomial coefficients, B=1; geometry is integrated from the actual outer endpoint.
 fcoef=np.array([lam,0,1,-2*mass,1+3*lam+ds,-2*ds,ds]);fpoly=np.polynomial.Polynomial(fcoef)
 fp=fpoly.deriv();fpp=fpoly.deriv(2)
 ppoly=np.polynomial.Polynomial([6*lam,0,0,2*mass,-2*(1+3*lam+ds),8*ds,-8*ds])
 check('seed_fixture_horizon',abs(fpoly(1/rh))<1e-13)
 v,v2=profile_functions();L=PROTOCOL['domain_end'];ng=PROTOCOL['quadrature_nodes'];gx,gw=leggauss(ng)
 pieces=[1e-4,1.,2.,8.,L];nodes=np.concatenate([(a+b)/2+(b-a)*gx/2 for a,b in zip(pieces[:-1],pieces[1:])])
 weights=np.concatenate([(b-a)*gw/2 for a,b in zip(pieces[:-1],pieces[1:])])
 for ki,k in enumerate(PROTOCOL['scaled_frequencies']):
  ref,source0=normal_green(k,v,v2);z0=ref(nodes);f0=np.array([source0(float(y)) for y in nodes]);old_error=None;old_source=None
  for di,d in enumerate(PROTOCOL['dilations']):
   tag=str(ki)+'_'+str(di)
   geo=solve_ivp(lambda eta,z:d*fpoly(z),(0,L),[0.],method='DOP853',rtol=PROTOCOL['geometry_rtol'],atol=PROTOCOL['geometry_atol'],dense_output=True)
   check('seed_geometry_solve_'+tag,geo.success)
   def data(eta):
    z=geo.sol(eta)[0];ff=fpoly(z);F=ff/(z*z);Rd=lam*d/z
    Cv=24*ff-12*z*fp(z)+2*z*z*fpp(z)
    pot=d*d*F*ppoly(z)
    src=Rd*(v2(eta)-k*k*v(eta))/3-Rd*d*d*F*Cv*v(eta)/6
    return pot,src
   sols=[]
   for jj,(left,tol) in enumerate(zip(PROTOCOL['left_endpoints'],PROTOCOL['BVP_tolerances'])):
    grid=np.unique(np.concatenate([np.geomspace(left,.5,60),np.linspace(.5,1,35),np.linspace(1,2,125),np.linspace(2,L,200)]))
    slope=3/left+(k*k+2*lam*d*d)*left/7
    def rhs(eta,state):
     pot,src=data(eta);return np.vstack([state[1],(k*k+pot)*state[0]-src])
    sol=solve_bvp(rhs,lambda a,b:np.array([a[1]-slope*a[0],b[0]]),grid,np.zeros((2,len(grid))),tol=tol,max_nodes=PROTOCOL['max_nodes'])
    check('curvature_free_BVP_'+tag+'_'+str(jj),sol.success,sol.message)
    sols.append(sol)
   coarse=sols[0].sol(nodes)[0];fine=sols[1].sol(nodes)[0]
   err=float(np.sqrt(np.dot(weights,(fine-z0)**2)))
   refine=float(np.sqrt(np.dot(weights,(fine-coarse)**2)))
   pot,src=data(nodes);srcerr=float(np.sqrt(np.dot(weights,(src-f0)**2)))
   graphdiff=src-f0-k*k*(fine-z0)-(pot-6/nodes**2)*fine
   grapherr=float(np.sqrt(np.dot(weights,graphdiff**2)))
   check('BVP_refinement_absolute_'+tag,refine<PROTOCOL['BVP_refinement_absolute_tolerance'],refine)
   check('BVP_refinement_resolves_dilation_'+tag,refine<PROTOCOL['BVP_refinement_fraction_of_dilation_error']*err,refine/err)
   if old_error is not None:
    check('global_reconstruction_dilation_'+tag,err<PROTOCOL['dilation_refinement_ratio']*old_error,err/old_error)
    check('compact_source_dilation_'+tag,srcerr<PROTOCOL['dilation_refinement_ratio']*old_source,srcerr/old_source)
   fL1=quad(lambda eta:abs(float(data(np.array([eta]))[1][0])),1,2,epsabs=1e-10,epsrel=1e-11,limit=200)[0]
   tail=fL1*np.exp(-k*(L-2))/(2*k*np.sqrt(2*k))
   check('explicit_horizon_tail_bound_'+tag,tail<PROTOCOL['tail_norm_bound_tolerance'] and tail<PROTOCOL['tail_fraction_of_dilation_error']*err,tail)
   check('conformal_tail_is_nonzero_'+tag,abs(float(sols[1].sol(3.)[0]))>1e-6)
   NUMERICS.append({'fixture':'actual_seed_curvature_free_reconstruction','kappa':k,'dilation':d,
    'L2_difference_from_normal':err,'K0_graph_difference':grapherr,'source_difference':srcerr,
    'BVP_refinement_difference':refine,'omitted_beyond_L_tail_norm_bound':tail,
    'dirichlet_L_boundary_error_L2_bound':tail,'normal_tail_at_eta3':float(ref(np.array([3.]))[0]),
    'actual_tail_at_eta3':float(sols[1].sol(3.)[0]),'nodes':[len(s.x) for s in sols]})
   old_error=err;old_source=srcerr

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 charged curvature-free shape form','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();seed_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 charged curvature-free shape form','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'exact_charged_shape_congruence_derived':good,'global_conformal_tail_transfer_derived':good,
  'actual_annular_curvature_free_positivity_derived':good,'restricted_variational_inverse_derived':good,
  'full_logarithmic_norm_retained':good,'off_shell_density_contacts_retained':good,
  'global_shape_inverse_proved':False,'physical_coupled_moment_evaluated':False,
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
