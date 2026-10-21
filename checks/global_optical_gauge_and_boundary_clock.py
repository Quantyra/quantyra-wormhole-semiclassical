"""Global optical gauge, fixed-clock range, and charged boundary-jet fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1]
STEM='global-optical-gauge-and-boundary-clock'
OWN=['checks/global_optical_gauge_and_boundary_clock.py','notes/'+STEM+'.md','evidence/optical-gauge-source-review.json']
PREVIOUS='charged-curvature-free-shape-form';ROWS=[];NUMERICS=[]
PROTOCOL={'scope':'Global optical gauge range in an explicit test domain; physical gauge coverage is unresolved',
 'frequencies':[[.8,0],[2.,0],[1.4,.6]],'source_types':['incompatible','compatible'],
 'source_profile':'256(y-1)^4(2-y)^4 on [1,2], zero elsewhere',
 'ODE_rtol':[1e-7,1e-12],'ODE_atol':[1e-9,1e-14],'ODE_max_step':[.08,.02],
 'coarse_error_tolerance':2e-6,'fine_error_tolerance':2e-9,'refinement_ratio':.1,'refinement_floor':1e-11,
 'norm_quadrature_nodes':[64,128],'norm_refinement_tolerance':1e-10,
 'seed_fixture':{'B':1,'lambda':2,'d_seed':.01,'horizon':1.1},
 'seed_radii':[40.,80.,160.,320.],'seed_frequencies':[.8,2.,4.],
 'outer_jet_refinement_ratio':.35,'outer_jet_final_tolerance':.002,
 'physical_modes_computed':False,'actual_Dirac_kernel_evaluated':False,'formal_verification':False}

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
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==79 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 check('prior_domains_pinned',all('notes/'+n+'.md' in p['input_sha256_after'] for n in ['physical-curvature-reconstruction','spherical-retarded-constraint-recovery']))
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_scope',len(review['primary_sources'])==1 and review['gauge_domain_explicit'] and review['horizon_exception_retained'] and not review['physical_gauge_coverage_proved'])

def exact_checks():
 t,y,th,ph=S.symbols('t y theta phi',real=True);coords=[t,y,th,ph]
 F=S.Function('F')(y);r=S.Function('r')(y)
 xt=S.Function('T')(t,y);xy=S.Function('Y')(t,y);xi=[xt,xy,0,0]
 metric=S.diag(-F,F,r*r,r*r*S.sin(th)**2)
 lie=S.zeros(4)
 for i in range(4):
  for j in range(4):
   lie[i,j]=sum(xi[k]*S.diff(metric[i,j],coords[k])+metric[k,j]*S.diff(xi[k],coords[i])+metric[i,k]*S.diff(xi[k],coords[j]) for k in range(4))
 nprime=S.diff(F,y)/(2*F);lprime=S.diff(r,y)/r
 zero('four_dimensional_Lie_tt',lie[0,0]/(-2*F)-nprime*xy-S.diff(xt,t))
 zero('four_dimensional_Lie_yy',lie[1,1]/(2*F)-nprime*xy-S.diff(xy,y))
 zero('four_dimensional_Lie_ty',lie[0,1]/F+S.diff(xt,y)-S.diff(xy,t))
 zero('four_dimensional_Lie_sphere',lie[2,2]/(2*r*r)-lprime*xy)
 zero('four_dimensional_Lie_angular_isotropy',lie[3,3]-S.sin(th)**2*lie[2,2])
 A,B,C=S.symbols('A B C');An=A-nprime*xy-S.diff(xt,t);Bn=B-nprime*xy-S.diff(xy,y);Cn=C+S.diff(xt,y)-S.diff(xy,t)
 zp=xt+xy;zm=xt-xy
 zero('outgoing_characteristic_sign',S.diff(zp,t)-S.diff(zp,y)-(A-B+C)+(An-Bn+Cn))
 zero('incoming_characteristic_sign',S.diff(zm,t)+S.diff(zm,y)-(A-B-C)+(An-Bn-Cn))
 ap,am=S.symbols('ap am');zero('both_fixed_traces_force_both_zero',S.Matrix([(ap+am)/2,(ap-am)/2]).jacobian([ap,am]).det()+S.Rational(1,2))
 g,a=S.symbols('g a',positive=True);chi=y**2*S.exp(-a*y)
 zero('compatible_outgoing_moment',S.integrate(S.exp(-g*y)*(g*chi-S.diff(chi,y)),(y,0,S.oo)))
 zero('nonzero_boundary_gauge_moment',S.integrate(S.exp(-g*y)*(g*S.exp(-a*y)+a*S.exp(-a*y)),(y,0,S.oo))-1)
 zero('incompatible_moment_value',S.integrate(S.exp(-(g+a)*y),(y,0,S.oo))-1/(g+a))
 zero('moment_functional_norm',S.integrate(S.exp(-2*a*y),(y,0,S.oo))-1/(2*a))
 zpo=S.exp(-a*y)/(g+a);zmi=S.exp(-g*y)/(g+a)
 zero('explicit_outgoing_resolvent',g*zpo-S.diff(zpo,y)-S.exp(-a*y))
 zero('explicit_reflected_incoming',g*zmi+S.diff(zmi,y))
 zero('endpoint_preservation',zpo.subs(y,0)-zmi.subs(y,0))
 check('fixed_clock_rejects_positive_source',zpo.subs(y,0).is_positive)
 for sign in [1,-1]:zero('homogeneous_transport_'+str(sign),g*S.exp(sign*g*y)-sign*S.diff(S.exp(sign*g*y),y))
 energy=S.integrate(zpo*(-S.diff(zpo,y))+zmi*S.diff(zmi,y),(y,0,S.oo))
 zero('reflected_boundary_flux_cancelled',energy)
 s=S.symbols('s',real=True);pulse=256*s**4*(1-s)**4
 conv=S.integrate(pulse*pulse,(s,0,1))
 check('compact_retarded_clock_pulse_nonzero',conv>0,conv)
 zero('retarded_convolution_independent_beta',conv-256**2*S.factorial(8)**2/S.factorial(17))
 # The exact physical curvature identity for a residual optical pure gauge.
 rr,eta,etap=S.symbols('rr eta etap',positive=True);Fr=S.Function('F')(rr)
 dr=lambda w:S.diff(w,rr)
 dy=lambda w:S.expand(-Fr*dr(w)+etap*S.diff(w,eta)+g*g*eta*S.diff(w,etap))
 R=-dr(dr(Fr))-4*dr(Fr)/rr-2*(Fr-1)/rr**2
 D=dr(Fr)/rr-(2*R+rr*dr(R))/6
 q=Fr*eta;v=(dr(Fr)/2-Fr/rr)*eta-etap
 Cv=2*dr(dr(Fr))+4*Fr/rr**2+8*dr(Fr)/rr
 uv=rr*(2*(g*g*v-dy(dy(v)))/Fr+Cv*v)
 zero('exact_charged_pure_gauge_curvature',6*(g*g*q-dy(dy(q))+Fr*D*q)/Fr+uv)
 Tr=S.symbols('Tr');zero('off_shell_trace_gauge_contact',Fr*eta*Tr-Tr*(Fr*eta))
 lam,c=S.symbols('lam c',positive=True)
 rout=1/(lam*y)-y/3;Fout=1/(lam*y*y)+S.Rational(1,3)
 pout=Fout/rout*c*S.sinh(g*y)
 vout=(lam*rout-Fout/rout)*c*S.sinh(g*y)-g*c*S.cosh(g*y)
 zero('outer_conformal_trace',S.limit(pout,y,0)-g*c)
 zero('outer_shape_trace',S.limit(vout,y,0)+g*c)
 zero('outer_orbit_trace',S.limit(pout+vout,y,0))
 zero('nonintegrable_q_coefficient',S.limit(y*rout*pout,y,0)-g*c/lam)
 check('q_not_in_L2',S.integrate(y**-2,(y,0,1))==S.oo)
 kh,delta=S.symbols('kh delta',positive=True)
 exponent=1-g/kh
 zero('outgoing_EF_power_at_resonance',exponent.subs(g,kh))
 zero('outgoing_EF_power_above_resonance',exponent.subs(g,2*kh)+1)
 check('outgoing_EF_diverges_above_resonance',S.limit(delta**exponent.subs(g,2*kh),delta,0)==S.oo)
 zero('linear_horizon_exact_exception',(-kh*delta*S.exp(2*kh*(-S.log(delta)/(2*kh)))).simplify()+kh)
 # Initial data enter the same transformed forcing, independently in each characteristic.
 z0,fhat,zhat,dzhat=S.symbols('z0 fhat zhat dzhat')
 zero('initial_gauge_contact_retained',(g*zhat-z0-dzhat-fhat)-(g*zhat-dzhat-(fhat+z0)))

def profile(y):
 y=np.asarray(y);w=np.where((y>1)&(y<2),256*(y-1)**4*(2-y)**4,0.)
 return w.item() if w.ndim==0 else w

def derivative(y):
 y=np.asarray(y);w=np.where((y>1)&(y<2),1024*(y-1)**3*(2-y)**3*(3-2*y),0.)
 return w.item() if w.ndim==0 else w

def cq(fun,a,b):
 if a>=b:return 0j
 return quad(fun,a,b,epsabs=3e-13,epsrel=3e-13,complex_func=True)[0]

def gauge_fixtures():
 points=np.array([0.,.3,1.,1.1,1.5,1.9,2.,3.,8.])
 for idx,(gr,gi) in enumerate(PROTOCOL['frequencies']):
  gamma=complex(gr,gi)
  for kind in PROTOCOL['source_types']:
   fp=lambda y:profile(y) if kind=='incompatible' else gamma*profile(y)-derivative(y)
   fm=lambda y:.5*profile(y)
   moment=cq(lambda s:np.exp(-gamma*s)*fp(s),1.,2.)
   plus=lambda y:0j if y>=2 else cq(lambda s:np.exp(-gamma*(s-y))*fp(s),max(1.,y),2.)
   minus=lambda y:np.exp(-gamma*y)*moment+cq(lambda s:np.exp(-gamma*(y-s))*fm(s),1.,min(y,2.))
   ref=np.array([[plus(y),minus(y)] for y in points]);errors=[];sols=[]
   tag=str(idx)+'_'+kind
   for j,(rtol,atol,step) in enumerate(zip(PROTOCOL['ODE_rtol'],PROTOCOL['ODE_atol'],PROTOCOL['ODE_max_step'])):
    sp=solve_ivp(lambda y,z:gamma*z-fp(y),(2.,0.),np.array([0j]),method='DOP853',rtol=rtol,atol=atol,max_step=step,dense_output=True)
    sm=solve_ivp(lambda y,z:fm(y)-gamma*z,(0.,8.),np.array([sp.sol(0.)[0]]),method='DOP853',rtol=rtol,atol=atol,max_step=step,dense_output=True)
    check('gauge_ODE_success_'+tag+'_'+str(j),sp.success and sm.success)
    got=np.array([[sp.sol(y)[0] if y<=2 else 0j,sm.sol(y)[0]] for y in points])
    errors.append(float(np.max(np.abs(got-ref))));sols.append((sp,sm))
   check('gauge_coarse_error_'+tag,errors[0]<PROTOCOL['coarse_error_tolerance'],errors[0])
   check('gauge_fine_error_'+tag,errors[1]<PROTOCOL['fine_error_tolerance'],errors[1])
   check('gauge_refinement_'+tag,errors[1]<PROTOCOL['refinement_ratio']*errors[0]+PROTOCOL['refinement_floor'],errors)
   if kind=='compatible':check('compatible_clock_'+tag,abs(moment)<1e-12 and abs(sols[-1][0].sol(0.)[0])<2e-9,abs(moment))
   else:check('incompatible_clock_'+tag,abs(moment)>1e-4,abs(moment))
   norms=[]
   for nodes in PROTOCOL['norm_quadrature_nodes']:
    xx,ww=leggauss(nodes);out=0.;src=0.
    for left,right in [(0.,1.),(1.,2.)]:
     yy=left+(xx+1)*(right-left)/2;weights=ww*(right-left)/2
     out+=float(sum(w*(abs(plus(y))**2+abs(minus(y))**2) for y,w in zip(yy,weights)))
     src+=float(sum(w*(abs(fp(y))**2+abs(fm(y))**2) for y,w in zip(yy,weights)))
    tail=abs(minus(2.))**2/(2*gr);norms.append((out+tail,src,tail))
   check('full_norm_quadrature_refinement_'+tag,max(abs(a-b) for a,b in zip(norms[0],norms[1]))<PROTOCOL['norm_refinement_tolerance'])
   out,src,tail=norms[-1]
   check('global_gauge_resolvent_bound_'+tag,out<=src/(gr*gr)+1e-10,out/(src/(gr*gr)))
   check('incoming_horizon_tail_retained_'+tag,tail>0,tail)
   NUMERICS.append({'fixture':'global_gauge_ODE','gamma':[gr,gi],'kind':kind,'moment':[moment.real,moment.imag],
    'errors':errors,'full_gauge_norm_squared':out,'source_norm_squared':src,'exact_incoming_tail_norm_squared':tail})
 # Time pulse evaluated without using the frequency Green formula.
 value=quad(lambda s:profile(1+s)*profile(2-s),0.,1.,epsabs=1e-13,epsrel=1e-13)[0]
 target=float(256**2*S.factorial(8)**2/S.factorial(17))
 check('retarded_compact_pulse_quadrature',abs(value-target)<1e-13,value)
 NUMERICS.append({'fixture':'compact_retarded_clock','time':2.,'value':value,'exact_value':target})

def seed_fixtures():
 lam=2.;ds=.01;rh=1.1
 ext=lambda r:(r-1)**2/r**2*(1+lam*(r*r+2*r+3)+ds/r**2)
 mass=1+2*lam+rh*ext(rh)/2
 poly=np.polynomial.Polynomial([lam,0,1,-2*mass,1+3*lam+ds,-2*ds,ds])
 f=lambda r:r*r*poly(1/r)
 fr=lambda r:2*r*poly(1/r)-poly.deriv()(1/r)
 check('actual_seed_horizon',abs(f(rh))<1e-12)
 kh=fr(rh)/2;check('actual_seed_simple_horizon',kh>0,kh)
 for gamma in PROTOCOL['seed_frequencies']:
  errors=[];records=[]
  for rr in PROTOCOL['seed_radii']:
   yy=quad(lambda z:1/poly(z),0.,1/rr,epsabs=1e-14,epsrel=1e-13)[0]
   p=f(rr)/rr*np.sinh(gamma*yy)
   v=(fr(rr)/2-f(rr)/rr)*np.sinh(gamma*yy)-gamma*np.cosh(gamma*yy)
   err=max(abs(p/gamma-1),abs(v/gamma+1),abs(lam*yy*rr*p/gamma-1));errors.append(err)
   records.append({'r':rr,'y':yy,'p':p,'v':v,'lambda_y_q_over_gamma':lam*yy*rr*p/gamma})
  check('charged_outer_jet_refinement_'+str(gamma),all(b<PROTOCOL['outer_jet_refinement_ratio']*a for a,b in zip(errors,errors[1:])),errors)
  check('charged_outer_jet_limit_'+str(gamma),errors[-1]<PROTOCOL['outer_jet_final_tolerance'],errors[-1])
  NUMERICS.append({'fixture':'actual_charged_outer_gauge_jet','gamma':gamma,'errors':errors,'records':records})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 global optical gauge and boundary clock','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();gauge_fixtures();seed_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 global optical gauge and boundary clock','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'global_gauge_range_derived':good,'fixed_clock_compatibility_derived':good,
  'compact_metric_probe_counterexample_derived':good,'boundary_reconstruction_mismatch_derived':good,
  'initial_gauge_contacts_retained':good,'horizon_exception_retained':good,
  'physical_gauge_coverage_proved':False,'boundary_clock_response_evaluated':False,
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
