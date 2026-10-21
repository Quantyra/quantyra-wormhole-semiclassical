"""Boundary-clock reconstruction: exact weighted bounds and independent BVP fixtures."""
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

ROOT=Path(__file__).resolve().parents[1];STEM='boundary-clock-curvature-reconstruction'
OWN=['checks/boundary_clock_curvature_reconstruction.py','notes/'+STEM+'.md','evidence/boundary-clock-reconstruction-source-review.json']
PREVIOUS='global-optical-gauge-and-boundary-clock';ROWS=[];NUMERICS=[]
PROTOCOL={'scope':'Auxiliary weighted reconstruction with clock boundary trace and physical-L2 curvature feedback',
 'seed_fixture':{'B':1,'lambda':2,'d_seed':.01,'horizon':1.1},'real_frequencies':[64.,128.,256.],
 'profile':'one below eta=1; 1-35s^4+84s^5-70s^6+20s^7 for s=eta-1 in [0,1]; zero above eta=2',
 'subtraction_cutoff_scales':[1.,1.5],'left_endpoints':[1e-4,5e-5],
 'right_endpoints':[24.,32.],'BVP_tolerances':[1e-7,1e-10],'max_nodes':30000,
 'geometry_rtol':2e-13,'geometry_atol':2e-15,'quadrature_nodes':96,
 'BVP_refinement_tolerance':2e-7,'cutoff_independence_tolerance':2e-7,
 'normal_Green_tolerance':2e-7,'boundary_trace_tolerance':1e-5,
 'finite_domain_fixtures_are_physical_boundaries':False,'physical_curvature_norm_changed':False,
 'actual_boundary_shape_response_evaluated':False,'actual_Dirac_kernel_evaluated':False,'formal_verification':False}

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
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==101 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 check('reconstruction_inputs_pinned',all('notes/'+n+'.md' in p['input_sha256_after'] for n in ['physical-curvature-reconstruction','explicit-curvature-metric-source','charged-curvature-free-shape-form']))
 review=json.loads((ROOT/OWN[2]).read_bytes());check('source_scope',len(review['primary_sources'])==1 and not review['physical_curvature_norm_changed'] and not review['actual_boundary_shape_response_evaluated'])

def exact_checks():
 r,B,lam,g=S.symbols('r B lam g',positive=True);F=S.Function('F')(r);z=B/r
 dr=lambda a:S.diff(a,r);dy=lambda a:-F*dr(a)
 h=z**(-S.Rational(3,2));b=z**S.Rational(1,2)/F
 hterm=S.Rational(3,2)*F*dr(F)/r+S.Rational(3,4)*F*F/r**2
 bterm=dr(F)**2-F*dr(dr(F))+F*dr(F)/(2*r)+3*F*F/(4*r*r)
 zero('first_weighted_supersolution_derivative',dy(dy(h))/h-hterm)
 zero('second_weighted_supersolution_derivative',dy(dy(b))/b-bterm)
 zero('weighted_Schur_dual_relation',F*b-z*z*h)
 zz,x,m,d=S.symbols('z x m d',positive=True)
 f=x+zz**2-2*m*zz**3+(1+3*x+d)*zz**4-2*d*zz**5+d*zz**6
 P=6*x+2*m*zz**3-2*(1+3*x+d)*zz**4+8*d*zz**5-8*d*zz**6
 fz=S.diff(f,zz);fzz=S.diff(f,zz,2)
 N1=f*(S.Rational(15,4)*f-S.Rational(3,2)*zz*fz+x-P)
 N2=(2*f-zz*fz)**2-f*(2*f-2*zz*fz+zz**2*fzz)+f*(2*f-zz*fz)/2+3*f*f/4+x*f-f*P
 Fr=f/zz**2;Fprime=(2*f-zz*fz)/(B*zz);Fsecond=(2*f-2*zz*fz+zz**2*fzz)/B**2
 M1=3*Fr*Fprime/(2*(B/zz))+3*Fr**2/(4*(B/zz)**2)+(x-P)*Fr/B**2
 M2=Fprime**2-Fr*Fsecond+Fr*Fprime/(2*(B/zz))+3*Fr**2/(4*(B/zz)**2)+(x-P)*Fr/B**2
 zero('first_threshold_polynomial_identity',B*B*zz*zz*M1-N1)
 zero('second_threshold_polynomial_identity',B*B*zz*zz*M2-N2)
 for i,N in enumerate([N1,N2]):
  zero('threshold_outer_constant_'+str(i),N.subs(zz,0)+5*x*x/4)
  zero('threshold_no_linear_term_'+str(i),S.diff(N,zz).subs(zz,0))
  H=S.cancel((N+5*x*x/4)/zz**2);check('threshold_polynomial_'+str(i),H.is_polynomial(zz))
 zero('weighted_energy_derivative',dy(dy(z**4))/(2*z*z)-(6*dy(z)**2+2*z*dy(dy(z))))
 mf,mw,un,qn=S.symbols('mf mw un qn',positive=True)
 zero('derivative_bound_constant',(mf*un*qn/6+mw*qn**2).subs(qn,un/(6*lam))-un**2*(mf/lam+mw/lam**2)/36)
 chi=S.Function('chi')(r)
 R=-dr(dr(F))-4*dr(F)/r-2*(F-1)/r**2;HR=2*R+r*dr(R);D=dr(F)/r-HR/6
 Cv=2*dr(dr(F))+4*F/r**2+8*dr(F)/r
 Lop=lambda a:-dy(dy(a))+(F*D+g*g)*a
 Uv=lambda a:r*(2*(g*g*a-dy(dy(a)))/F+Cv*a)
 uchi=4*r*(dy(dy(chi))-g*g*chi)/F-12*dy(chi)+r*(Cv+HR)*chi
 zero('exact_boundary_lift_curvature',6*Lop(-r*chi)/F+Uv(chi)-uchi)
 delta=S.Function('delta')(r)
 ud=4*r*(dy(dy(delta))-g*g*delta)/F-12*dy(delta)+r*(Cv+HR)*delta
 zero('cutoff_change_compensated',Uv(delta)-ud-6*Lop(r*delta)/F)
 Q,c3,c4,mm=S.symbols('Q c3 c4 mm');seed=lam*r*r+1-2*mm/r+Q/r**2+c3/r**3+c4/r**4
 replace=lambda a:S.simplify(a.subs(F,seed).doit())
 zero('charged_lift_cancellation',replace(Cv+HR)-(4/r**2+10*c3/r**5+36*c4/r**6))
 uone=replace(uchi.subs(chi,1).doit())
 zero('boundary_lift_outer_coefficient',S.limit(r*uone,r,S.oo)-4*(lam-g*g)/lam)
 # Here y~1/(lambda r), so r*u_chi ->4(lambda-gamma^2)/lambda.
 y=S.symbols('y',positive=True)
 check('clock_q_in_weighted_space',S.integrate(y*y*(1/y)**2,(y,0,1))==1)
 check('clock_q_outside_old_space',S.integrate((1/y)**2,(y,0,1))==S.oo)
 check('homogeneous_singular_branch_excluded',S.integrate(y*y*(1/y**2)**2,(y,0,1))==S.oo)
 eta=S.symbols('eta',positive=True);cc=S.Function('chi')(eta)
 e0source=2*(cc-S.diff(cc,eta,2))/(3*eta)+2*S.diff(cc,eta)/eta**2
 normal_total=(S.diff(cc,eta,2)-cc)/(3*eta)-4*cc/eta**3
 normal_lift=-cc/eta
 zero('stable_normal_subtraction',normal_total-(-S.diff(normal_lift,eta,2)+(6/eta**2+1)*normal_lift)-e0source)
 alpha=S.symbols('alpha');zero('normal_particular_boundary_slope',(6*alpha/eta-2/(3*eta)).subs(alpha,S.Rational(1,9)))
 U,uv,iq,iv,rr,ff=S.symbols('U uv iq iv rr ff')
 zero('all_initial_contacts_in_single_source',ff*(U-uv+(6*iq+2*rr*iv)/ff)/6-(ff*(U-uv)/6+iq+rr*iv/3))
 cg,cn,fn=S.symbols('cg cn fn',positive=True)
 zero('clock_composition_uniform_bound',g*(fn/S.sqrt(2*g))*(cg/S.sqrt(g))-cg*fn/S.sqrt(2))

def cutoff_functions():
 s=S.symbols('s');expr=1-35*s**4+84*s**5-70*s**6+20*s**7
 fs=[S.lambdify(s,S.diff(expr,s,n),'numpy') for n in range(3)]
 def get(n,eta,scale=1.):
  eta=np.asarray(eta);ss=eta/scale-1
  values=fs[n](ss)/scale**n
  out=np.where(eta<=scale,1. if n==0 else 0.,np.where(eta<2*scale,values,0.))
  return out.item() if out.ndim==0 else out
 return get

def numerical_fixtures():
 cf=cutoff_functions();lam=2.;ds=.01;rh=1.1
 ext=lambda r:(r-1)**2/r**2*(1+lam*(r*r+2*r+3)+ds/r**2)
 mass=1+2*lam+rh*ext(rh)/2;Qh=1+3*lam+ds
 fp=np.polynomial.Polynomial([lam,0,1,-2*mass,Qh,-2*ds,ds]);fz=fp.deriv();fzz=fp.deriv(2)
 pp=np.polynomial.Polynomial([6*lam,0,0,2*mass,-2*Qh,8*ds,-8*ds])
 zp=np.polynomial.Polynomial([0,1]);N1=fp*(3.75*fp-1.5*zp*fz+lam-pp)
 N2=(2*fp-zp*fz)**2-fp*(2*fp-2*zp*fz+zp*zp*fzz)+.5*fp*(2*fp-zp*fz)+.75*fp*fp+lam*fp-fp*pp
 hnorms=[float(sum(abs((n+1.25*lam*lam).coef[2:]))) for n in [N1,N2]]
 old=(2*fp-zp*fz)**2-fp*(2*fp-2*zp*fz+zp*zp*fzz+pp-3*lam)
 threshold=max(np.sqrt(1+max(hnorms)),np.sqrt(1+sum(abs((old+lam*lam).coef[2:]))))
 check('fixture_frequencies_above_both_thresholds',min(PROTOCOL['real_frequencies'])>threshold,threshold)
 points=np.array([.1,.3,.7,1.,1.4,1.8,2.,2.5,3.,5.,8.,16.])
 gx,gw=leggauss(PROTOCOL['quadrature_nodes']);nodes=[];weights=[]
 for left,right in [(0.,1.),(1.,2.),(2.,3.),(3.,8.),(8.,24.)]:
  nodes.extend(left+(gx+1)*(right-left)/2);weights.extend(gw*(right-left)/2)
 nodes=np.array(nodes);weights=np.array(weights)
 for gamma in [None]+PROTOCOL['real_frequencies']:
  normal=gamma is None;tag='normal' if normal else str(gamma);d=0. if normal else 1/gamma
  if not normal:
   geo=solve_ivp(lambda eta,z:d*fp(z),(0.,32.),[0.],method='DOP853',rtol=PROTOCOL['geometry_rtol'],atol=PROTOCOL['geometry_atol'],dense_output=True)
   check('charged_geometry_'+tag,geo.success)
  def data(eta):
   eta=np.asarray(eta)
   if normal:return 1/eta,1/eta**2,6/eta**2,-24/eta**2,np.zeros_like(eta),None
   z=geo.sol(eta)[0];f=fp(z);F=f/z**2;Rd=lam*d/z
   Cv=24*f-12*z*fz(z)+2*z*z*fzz(z)
   HR=-24*lam-12*ds*z**5+24*ds*z**6
   HS=4*z*z-20*ds*z**5+36*ds*z**6
   return Rd,lam*d*d*F,d*d*F*pp(z),d*d*F*HR,d*d*F*HS,(z,F,Cv)
  sols=[]
  for scale in PROTOCOL['subtraction_cutoff_scales']:
   refined=[]
   for j,(left,end,tol) in enumerate(zip(PROTOCOL['left_endpoints'],PROTOCOL['right_endpoints'],PROTOCOL['BVP_tolerances'])):
    mesh=np.unique(np.concatenate([np.geomspace(left,.5,60),np.linspace(.5,1,30),np.linspace(1,3,180),np.linspace(3,end,160)]))
    alpha=(1-lam*d*d)/9
    def rhs(eta,state):
     Rd,J,pot,Hterm,HS,_=data(eta);V=cf(0,eta);Vp=cf(1,eta);Vpp=cf(2,eta)
     delta=cf(0,eta,scale)-V;dp=cf(1,eta,scale)-Vp;dpp=cf(2,eta,scale)-Vpp
     forcing=2*Rd*(V-Vpp)/3+2*J*Vp-Rd*HS*V/6+Rd*(delta-dpp)+2*J*dp-Rd*Hterm*delta/6
     return np.vstack([state[1],(pot+1)*state[0]-forcing])
    sol=solve_bvp(rhs,lambda a,b:np.array([a[1]-3*a[0]/left+2*alpha,b[0]]),mesh,np.zeros((2,len(mesh))),tol=tol,max_nodes=PROTOCOL['max_nodes'])
    check('boundary_BVP_'+tag+'_'+str(scale)+'_'+str(j),sol.success,sol.message)
    refined.append(sol)
   Rd,_,_,_,_,_=data(points)
   discrepancy=float(np.max(np.abs(refined[1].sol(points)[0]-refined[0].sol(points)[0])))
   check('boundary_BVP_refinement_'+tag+'_'+str(scale),discrepancy<PROTOCOL['BVP_refinement_tolerance'],discrepancy)
   sols.append(refined[-1])
  Rd,J,pot,Hterm,HS,extra=data(points)
  w1=sols[0].sol(points)[0]-Rd*cf(0,points,1.)
  w2=sols[1].sol(points)[0]-Rd*cf(0,points,1.5)
  difference=float(np.max(np.abs(w1-w2)))
  check('boundary_lift_cutoff_independence_'+tag,difference<PROTOCOL['cutoff_independence_tolerance'],difference)
  if normal:
   aa=lambda eta:np.sqrt(eta)*iv(2.5,eta);bb=lambda eta:np.sqrt(eta)*kv(2.5,eta)
   forcing=lambda eta:2*(cf(0,eta)-cf(2,eta))/(3*eta)+2*cf(1,eta)/eta**2
   green=lambda eta:bb(eta)*quad(lambda s:aa(s)*forcing(s),0,min(eta,2.),epsabs=1e-12,epsrel=1e-12,points=[1.] if min(eta,2.)>1 else None)[0]+(aa(eta)*quad(lambda s:bb(s)*forcing(s),eta,2.,epsabs=1e-12,epsrel=1e-12,points=[1.] if eta<1 else None)[0] if eta<2 else 0.)
   ref=np.array([green(eta) for eta in points]);err=float(np.max(np.abs(ref-sols[0].sol(points)[0])))
   check('independent_normal_Bessel_Green',err<PROTOCOL['normal_Green_tolerance'],err)
   NUMERICS.append({'fixture':'normal_boundary_reconstruction','Green_error':err,'cutoff_difference':difference})
   continue
  Rd,J,pot,Hterm,HS,extra=data(nodes);z,F,Cv=extra;V=cf(0,nodes);Vp=cf(1,nodes);Vpp=cf(2,nodes)
  e,ep=sols[0].sol(nodes);w=e-Rd*V;wp=ep+J*V-Rd*Vp
  q=gamma*w/lam;qy=gamma*gamma*wp/lam
  uv=(2*gamma*gamma*(V-Vpp)/F+Cv*V)/z
  qnorm=float(np.sqrt(np.dot(weights,(z*q)**2)/gamma));unorm=float(np.sqrt(np.dot(weights,(z*uv)**2)/gamma))
  check('weighted_reconstruction_bound_'+tag,qnorm<unorm/(6*lam),qnorm/unorm)
  left=PROTOCOL['left_endpoints'][-1];Rl,*_=data(np.array([left]));pl=sols[0].sol(left)[0]/Rl[0]-1
  check('nonzero_boundary_trace_'+tag,abs(pl+1)<PROTOCOL['boundary_trace_tolerance'],pl)
  Fr=(2*fp(z)-z*fz(z))/z
  dR=z**4*(28*mass-128*Qh*z/3+440*ds*z*z/3-100*ds*z**3)
  eR=z**6*(-44*Qh/3+460*ds*z/3-246*ds*z*z)
  L=-2*mass*z*z+4*Qh*z**3-20*ds*z**4+20*ds*z**5
  A1=F*dR/z;A0=F*(dR*Fr+eR*F)/z
  B1=2*F*eR/z;B0=2*eR*F*F+F*eR*Fr/z+F*dR*L
  feedback=A1*uv/6-B1*qy-B0*q-2*A0*V
  gnorm=float(np.sqrt(np.dot(weights,feedback**2)/gamma))
  check('boundary_feedback_finite_'+tag,np.isfinite(gnorm) and gnorm>0,gnorm)
  NUMERICS.append({'fixture':'charged_boundary_reconstruction','gamma':gamma,'common_frequency_threshold':threshold,
   'cutoff_difference':difference,'weighted_q_norm':qnorm,'weighted_source_norm':unorm,
   'boundary_p':float(pl),'bulk_curvature_feedback_norm':gnorm,'sqrt_gamma_feedback_norm':float(np.sqrt(gamma)*gnorm),
   'finite_domain_norms_are_interval_certified':False,'physical_boundary_response_evaluated':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 boundary clock curvature reconstruction','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();numerical_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 boundary clock curvature reconstruction','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'weighted_reconstruction_bound_derived':good,'exact_boundary_lift_derived':good,'cutoff_independence_derived':good,
  'extended_curvature_feedback_bounded':good,'clock_composition_bound_derived':good,'initial_contacts_retained':good,
  'physical_curvature_norm_changed':False,'actual_boundary_shape_response_evaluated':False,
  'physical_clock_closure_solved':False,'global_shape_inverse_proved':False,'physical_coupled_moment_evaluated':False,
  'physical_metric_boundary_selected':False,'actual_metric_causal_inverse_proved':False,'actual_metric_mode_constructed':False,
  'actual_charged_stability_proved':False,'actual_wormhole_instability_proved':False,'retained_fast_factor_removed':False,
  'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
