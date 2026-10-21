"""Verify clock subtraction, charged outer jets and logarithmic core density."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import numpy as np
import scipy
from scipy.integrate import solve_ivp, solve_bvp

ROOT=Path(__file__).resolve().parents[1];STEM='boundary-clock-response-extension'
OWN=['checks/boundary_clock_response_extension.py','notes/'+STEM+'.md','evidence/boundary-response-extension-source-review.json']
PREVIOUS='adm-dirac-boundary-response';ROWS=[];NUMERICS=[]
PROTOCOL={'scope':'Outer input extension and compact bulk readout, not a physical boundary-force evaluation',
 'gammas':[64.,128.],'eta_points':[.04,.02,.01,.005],
 'BVP_left':[1e-4,5e-5],'BVP_end':[24.,32.],'BVP_tolerance':[1e-7,1e-10],
 'BVP_max_nodes':25000,'BVP_refinement_tolerance':1e-6,
 'jet_refinement_ratio':.45,'jet_refinement_floor':1e-6,'jet_final_tolerance':2e-4,
 'cutoff_scales':[.08,.04,.02,.01],'FFT_sizes':[8192,16384],'FFT_box_length':8.,
 'FFT_relative_tolerance':2e-6,'tail_refinement_ratio':.75,
 'physical_kernel_evaluated':False,'formal_verification':False}

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
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==81 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('continuity_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_continuity_inputs']))
 check('source_and_causality_scope',len(review['primary_sources'])==2 and not review['cosh_lift_claimed_causal'] and not review['boundary_force_evaluated'])

def exact_checks():
 t,y=S.symbols('t y',positive=True);g,lam=S.symbols('g lam',positive=True)
 c=S.Function('c')(t);rho=S.Function('rho')(y)
 a=rho*(c+S.diff(c,t,2)*y*y/2);b=rho*(S.diff(c,t)*y+S.diff(c,t,3)*y**3/6)
 zero('causal_lift_kinetic',S.diff(a,t)-S.diff(b,y)+S.diff(rho,y)*(S.diff(c,t)*y+S.diff(c,t,3)*y**3/6))
 zero('causal_lift_shift',S.diff(b,t)-S.diff(a,y)-rho*S.diff(c,t,4)*y**3/6+S.diff(rho,y)*(c+S.diff(c,t,2)*y*y/2))
 zero('causal_boundary_clock',a.subs(y,0).subs(rho.subs(y,0),1)-c)
 zero('causal_boundary_location',b.subs(y,0))
 ae=rho*S.cosh(g*y)/g;be=rho*S.sinh(g*y)/g
 zero('exact_optical_lift_kinetic',g*ae-S.diff(be,y)+S.diff(rho,y)*S.sinh(g*y)/g)
 zero('exact_optical_lift_shift',g*be-S.diff(ae,y)+S.diff(rho,y)*S.cosh(g*y)/g)
 zero('polynomial_vs_cosh_time_order',S.limit((S.cosh(g*y)/g-1/g-g*y*y/2)/y**4,y,0)-g**3/24)
 zero('polynomial_vs_sinh_space_order',S.limit((S.sinh(g*y)/g-y-g*g*y**3/6)/y**5,y,0)-g**4/120)
 z,m,Q,c3,c4=S.symbols('z m Q c3 c4',real=True)
 F=lam/z**2+1-2*m*z+Q*z*z+c3*z**3+c4*z**4
 dr=lambda f:-z*z*S.diff(f,z)
 R=-dr(dr(F))-4*z*dr(F)+2*z*z*(1-F);HR=2*R+dr(R)/z
 D=z*dr(F)-HR/6;pot=S.expand(F*D)
 zser=lam*y+lam**2*y**3/3-m*lam**3*y**4/2+(2*lam**3/15+Q*lam**4/5)*y**5
 # Laurent coefficients are extracted before substitution to avoid an enormous nested rational series.
 zero('inverse_radius_IVP',S.series(S.diff(zser,y)-(lam+zser**2-2*m*zser**3+Q*zser**4),y,0,5).removeO())
 invz=1/(lam*y)-y/3+m*lam*y*y/2
 zero('inverse_radius_series',S.series(zser*invz-1,y,0,4).removeO())
 zero('reconstruction_potential_pole',pot.coeff(z,-2)-6*lam*lam)
 zero('reconstruction_potential_z_constant',pot.coeff(z,0)-6*lam)
 zero('reconstruction_potential_y_constant',-4*lam+pot.coeff(z,0)-2*lam)
 nyminusly=S.expand(F*z-dr(F)/2)
 zero('geometry_difference_exact',nyminusly-(z-3*m*z*z+2*Q*z**3+S.Rational(5,2)*c3*z**4+3*c4*z**5))
 ndiff=lam*y-3*m*lam*lam*y*y
 ly=-1/y-2*lam*y/3+3*m*lam*lam*y*y/2
 vp=S.series(S.cosh(g*y)+ndiff*S.sinh(g*y)/g,y,0,4).removeO()
 pp=S.series(ly*S.sinh(g*y)/g,y,0,4).removeO()
 zero('coordinate_shape_quadratic',vp.coeff(y,2)-g*g/2-lam)
 zero('coordinate_shape_cubic',vp.coeff(y,3)+3*m*lam*lam)
 zero('coordinate_conformal_quadratic',pp.coeff(y,2)+g*g/6+2*lam/3)
 zero('coordinate_conformal_cubic',pp.coeff(y,3)-3*m*lam*lam/2)
 Fu=S.expand(-2*g*g/(3*z)+F*(4*z+10*c3*z**4+36*c4*z**5)/6)
 leading=Fu.coeff(z,-1)*invz+Fu.coeff(z,1)*lam*y+Fu.coeff(z,0)
 sourcepole=2*(lam-g*g)/(3*lam);sourceone=(2*g*g+4*lam)/9
 zero('source_pole',S.expand(leading).coeff(y,-1)-sourcepole)
 zero('source_linear',S.expand(leading).coeff(y,1)-sourceone)
 a1=(lam-g*g)/(9*lam);alog=-(g*g+lam)*(g*g+2*lam)/(45*lam)
 K0=lambda q:-S.diff(q,y,2)+6*q/y**2
 zero('regular_particular_power',K0(y)-6/y)
 zero('resonant_particular_log',K0(y**3*S.log(y))+5*y)
 zero('particular_leading_coefficient',6*a1-sourcepole)
 zero('particular_log_coefficient',-5*alog+(2*lam+g*g)*a1-sourceone)
 pb2=-lam*a1;pr2=S.simplify(pb2-pp.coeff(y,2));vr2=-vp.coeff(y,2)
 zero('regular_conformal_coefficient',pr2-5*(g*g+2*lam)/18)
 zero('regular_shape_coefficient',vr2+(g*g+2*lam)/2)
 zero('curvature_free_quadratic_relation',36*pr2+20*vr2)
 check('genuine_remainder_not_zero',pr2.is_positive is True)
 # Resolvent identity at the spectral level, including the singular input map.
 E,G=S.symbols('E G',positive=True)
 zero('reference_frequency_resolvent_identity',(1+(G*G-g*g)/(E+g*g))/(E+G*G)-1/(E+g*g))
 d,k=S.symbols('d k',positive=True)
 zero('Fourier_fourth_order_scaling',d**5*(k/d)**4-d*k**4)
 zero('raw_boundary_value_wrong_scaling',d*(k/d)**4-k**4/d**3)
 # Direct local integrals for the leading regular profiles; cutoff derivatives have the same scaling.
 zero('weighted_conformal_tail',S.integrate((2*y/y)**2+(y*y/y**2)**2+4,(y,0,d))-9*d)
 zero('logarithmic_tail_integral',S.integrate(-S.log(y),(y,0,d))-d*(1-S.log(d)))
 c0,v0,ch=S.symbols('c0 v0 ch')
 zero('initial_clock_contact',g*((v0+c0)/g)-c0-v0)
 check('omitted_clock_initial_value_detected',S.simplify(g*(v0/g)-c0-v0)!=0)

def cutoff(n,x):
 x=np.asarray(x);u=x-1.
 coeff=np.polynomial.Polynomial([1.,0.,0.,0.,-35.,84.,-70.,20.]).deriv(n)
 return np.where(x<=1,1. if n==0 else 0.,np.where(x<2,coeff(u),0.))

def reconstruction_fixtures():
 lam=2.;ds=.01;rh=1.1
 ext=lambda r:(r-1)**2/r**2*(1+lam*(r*r+2*r+3)+ds/r**2)
 mass=1+2*lam+rh*ext(rh)/2;Q=1+3*lam+ds
 fp=np.polynomial.Polynomial([lam,0,1,-2*mass,Q,-2*ds,ds]);fz=fp.deriv()
 potpoly=np.polynomial.Polynomial([6*lam,0,0,2*mass,-2*Q,8*ds,-8*ds])
 points=np.array(PROTOCOL['eta_points'])
 for gamma in PROTOCOL['gammas']:
  d=1/gamma
  geo=solve_ivp(lambda eta,z:d*fp(z),(0.,32.),[0.],method='DOP853',rtol=1e-12,atol=1e-14,dense_output=True)
  check('charged_geometry_'+str(gamma),geo.success)
  def data(eta):
   z=geo.sol(eta)[0];f=fp(z);F=f/z**2
   return z,F,lam*d/z,lam*d*d*F,d*d*F*potpoly(z),d*d*F*(4*z*z-20*ds*z**5+36*ds*z**6)
  sols=[];alpha=(1-lam*d*d)/9
  for idx,(left,end,tol) in enumerate(zip(PROTOCOL['BVP_left'],PROTOCOL['BVP_end'],PROTOCOL['BVP_tolerance'])):
   mesh=np.unique(np.concatenate([np.geomspace(left,.5,90),np.linspace(.5,3,230),np.linspace(3,end,160)]))
   def rhs(eta,q):
    z,F,Rd,J,pot,HS=data(eta)
    forcing=2*Rd*(cutoff(0,eta)-cutoff(2,eta))/3+2*J*cutoff(1,eta)-Rd*HS*cutoff(0,eta)/6
    return np.vstack([q[1],(pot+1)*q[0]-forcing])
   sol=solve_bvp(rhs,lambda a,b:np.array([a[1]-3*a[0]/left+2*alpha,b[0]]),mesh,np.zeros((2,len(mesh))),tol=tol,max_nodes=PROTOCOL['BVP_max_nodes'])
   check('charged_BVP_'+str(gamma)+'_'+str(idx),sol.success,sol.message);sols.append(sol)
  z,F,Rd,J,pot,HS=data(points);Fr=(2*fp(z)-z*fz(z))/z
  ly=-F*z;ndiff=F*z-Fr/2
  pcoord=ly*np.sinh(points)/gamma
  vcoord=np.cosh(points)+ndiff*np.sinh(points)/gamma
  regular=[-1+sol.sol(points)[0]/Rd-pcoord for sol in sols]
  error=float(np.max(np.abs(regular[1]-regular[0])/points**2))
  check('charged_BVP_coefficient_refinement_'+str(gamma),error<PROTOCOL['BVP_refinement_tolerance'],error)
  pc=5*(1+2*lam*d*d)/18;vc=-(1+2*lam*d*d)/2
  pe=np.abs(regular[-1]/points**2-pc);ve=np.abs((1-vcoord)/points**2-vc)
  for name,errs in [('conformal',pe),('shape',ve)]:
   check('boundary_'+name+'_jet_refinement_'+str(gamma),all(b<PROTOCOL['jet_refinement_ratio']*a+PROTOCOL['jet_refinement_floor'] for a,b in zip(errs,errs[1:])),errs)
   check('boundary_'+name+'_jet_limit_'+str(gamma),errs[-1]<PROTOCOL['jet_final_tolerance'],errs[-1])
  NUMERICS.append({'fixture':'charged_boundary_subtraction','gamma':gamma,'eta':points.tolist(),'p_errors':pe.tolist(),'v_errors':ve.tolist(),'BVP_difference':error})

def Fourier_fixtures():
 records=[];raw=[]
 for delta in PROTOCOL['cutoff_scales']:
  values=[];rawvalues=[]
  for n in PROTOCOL['FFT_sizes']:
   length=PROTOCOL['FFT_box_length'];dx=length/n;x=(np.arange(n)-n//2)*dx
   psi=x*x*np.exp(-(delta*x)**2)*cutoff(0,np.abs(x))
   freq=2*np.pi*np.fft.fftfreq(n,d=dx);dk=2*np.pi/length
   ft=np.fft.fft(psi)*dx/np.sqrt(2*np.pi)
   ft_raw=np.fft.fft(cutoff(0,np.abs(x)))*dx/np.sqrt(2*np.pi)
   symbol=(1+(freq/delta)**2)**2*(1+np.log(2+(freq/delta)**2))
   values.append(float(delta**5*np.sum(symbol*abs(ft)**2)*dk/2))
   rawvalues.append(float(delta*np.sum(symbol*abs(ft_raw)**2)*dk/2))
  rel=abs(values[1]-values[0])/values[1]
  check('Fourier_grid_refinement_'+str(delta),rel<PROTOCOL['FFT_relative_tolerance'],rel)
  ratio=values[-1]/(delta*(1+abs(np.log(delta))))
  check('regular_tail_finite_scaled_norm_'+str(delta),np.isfinite(ratio) and ratio<1000,ratio)
  records.append(values[-1]);raw.append(rawvalues[-1])
  NUMERICS.append({'fixture':'logarithmic_cutoff_norm','delta':delta,'norm_squared':values,'scaled_ratio':ratio,'raw_boundary_norm_squared':rawvalues})
 check('regular_tail_norm_decreases',all(b<PROTOCOL['tail_refinement_ratio']*a for a,b in zip(records,records[1:])),records)
 check('raw_boundary_cutoff_diverges',all(b>4*a for a,b in zip(raw,raw[1:])),raw)

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 boundary response extension','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();reconstruction_fixtures();Fourier_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 boundary response extension','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'causal_clock_jet_lift_derived':good,'regular_boundary_remainder_derived':good,'compact_bulk_response_extension_derived':good,
  'complex_old_L2_reconstruction_derived':good,'initial_clock_contacts_retained':good,
  'cosh_lift_claimed_causal':False,'genuine_boundary_column_discarded':False,
  'counterterm_changed':False,'seed_residual_discarded':False,'boundary_force_evaluated':False,
  'physical_kernel_numerically_evaluated':False,'physical_clock_closure_solved':False,'global_shape_inverse_proved':False,
  'actual_metric_mode_constructed':False,'actual_charged_stability_proved':False,'actual_wormhole_instability_proved':False,
  'retained_fast_factor_removed':False,'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
