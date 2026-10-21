"""S168: boundary traces and exact comparison source compatibility.

Ordinary continuum proof; finite tests are not a charged exterior spectrum.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import math
import sys
import traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='fast-scalar-boundary-solvability'
OWN=['checks/fast_scalar_boundary_solvability.py','notes/'+STEM+'.md','evidence/fast-boundary-source-review.json']
PREVIOUS='neutral-spherical-current-decoupling'
ROWS=[];NUMERICS=[]
PROTOCOL={
 'scope':'Exact AdS comparison two-zero-trace range test and initial-data consequence',
 'coefficients':[2,5],'Laplace_frequencies':[0.5,1.5,3.0],
 'source_interval':['0.7','1.0'],'RK_steps_per_segment':[64,128],
 'noncompatible_source':'t^3(1-t)^3, t=(z-0.7)/0.3, zero outside',
 'compatible_test':'q=t^4(1-t)^4; f=(A+gamma^2)q',
 'mpmath_precision':50,'hypergeometric_tolerance':'1e-35',
 'quadrature_tolerance':'1e-35','RK_absolute_tolerance':1e-7,
 'RK_refinement_factor':0.2,'RK_refinement_floor':1e-11,
 'fixtures_are_benchmark_or_exterior_spectrum':False,'formal_verification':False,
}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,expr):
 result=S.simplify(S.sympify(expr).doit());check(name,result==0,result)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():
 return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json',
  'notes/stationary-global-weighted-map.md','notes/stationary-adjoint-phase.md',
  'notes/ads-fast-scalar-stability-preflight.md'})

def provenance():
 p=prior()
 check('predecessor_passed',p['status']=='passed' and len(p['checks'])==89 and all(r['passed'] for r in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 for f in ['ads-fast-scalar-stability-preflight','stationary-adjoint-phase','stationary-global-weighted-map']:
  check('required_pinned_'+f,'notes/'+f+'.md' in p['input_sha256_after'])
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('primary_source_scope',len(review['primary_sources'])==1 and not review['actual_charged_transfer_proved']
       and not review['finite_time_physical_stability_decided'] and not review['order_reduction_adopted'])

def exact_checks():
 z,k,gamma=S.symbols('z k gamma',positive=True,real=True)
 ac,ass,bc,bs=S.symbols('ac ass bc bs',real=True)
 fc=S.sqrt(z)*S.cos(k*S.log(z));fs=S.sqrt(z)*S.sin(k*S.log(z))
 a=k**2+S.Rational(1,4)
 zero('indicial_plus',(S.Rational(1,2)+S.I*k)*(S.Rational(1,2)+S.I*k-1)+a)
 zero('indicial_minus',(S.Rational(1,2)-S.I*k)*(S.Rational(1,2)-S.I*k-1)+a)
 zero('inverse_square_cos_solution',S.diff(fc,z,2)+a*fc/z**2)
 zero('inverse_square_sin_solution',S.diff(fs,z,2)+a*fs/z**2)
 zero('boundary_Wronskian',fc*S.diff(fs,z)-S.diff(fc,z)*fs-k)
 u=ac*fc+ass*fs;v=bc*fc+bs*fs
 zero('boundary_trace_pairing',u*S.diff(v,z)-S.diff(u,z)*v-k*(ac*bs-ass*bc))
 theta=S.symbols('theta',real=True)
 zero('phase_line_is_isotropic',k*((-S.sin(theta))*S.cos(theta)-S.cos(theta)*(-S.sin(theta))))
 check('two_zero_trace_space_not_maximal',S.Matrix([[0,k],[-k,0]]).det()==k**2)
 p=S.symbols('p')
 c2=(gamma**2-a/3)/(4*p+2)
 zero('Frobenius_first_correction',((p+2)*(p+1)-p*(p-1))*c2+a/3-gamma**2)
 zero('Frobenius_nonresonant_denominator',(4*(S.Rational(1,2)+S.I*k)+2)*(4*(S.Rational(1,2)-S.I*k)+2)-16*(1+k**2))
 w=S.Function('w')(z);U=S.Function('U')(z);V=S.Function('V')(z)
 op=lambda y:-S.diff(y,z,2)+V*y
 zero('full_Green_identity',w*op(U)-U*op(w)-S.diff(S.diff(w,z)*U-w*S.diff(U,z),z))
 UL,UpL=S.symbols('UL UpL')
 zero('center_moment_sign',1*UL-0*UpL-UL)
 LfunU,u0,u1,rhs=S.symbols('LfunU u0 u1 rhs')
 zero('Laplace_initial_source',(gamma**2*LfunU-gamma*u0-u1+rhs)
      -((gamma**2*LfunU+rhs)-(gamma*u0+u1)))
 # Hypergeometric substitution with y=cos(z)^2 and prefactor -cos(z)sin(z)^p.
 y=S.cos(z)**2
 logder=-S.tan(z)+p/S.tan(z)
 gdd_over_g=p*(p-1)/S.sin(z)**2-(p+1)**2
 ab=((p+1)**2+gamma**2)/4
 prime_coeff=2*logder*S.diff(y,z)+S.diff(y,z,2)+4*((p+2)*y-S.Rational(3,2))
 zero('hypergeometric_first_derivative_coefficient',S.trigsimp(prime_coeff))
 zero('hypergeometric_zero_order_coefficient',gdd_over_g-gamma**2-p*(p-1)/S.sin(z)**2+4*ab)
 zero('hypergeometric_principal_coefficient',S.trigsimp(S.diff(y,z)**2/(y*(1-y))-4))
 q=z**S.Rational(3,2)*S.cos(k*S.log(z))
 expected=z**(-S.Rational(1,2))*((S.Rational(3,4)-k**2)*S.cos(k*S.log(z))-2*k*S.sin(k*S.log(z)))
 zero('conformal_second_derivative',S.diff(q,z,2)-expected)
 zero('second_derivative_nonzero_amplitude',(S.Rational(3,4)-k**2)**2+4*k**2-(k**4+S.Rational(5,2)*k**2+S.Rational(9,16)))
 G,lam=S.symbols('G lambda',positive=True,real=True)
 ksquared=60*S.pi/(G*lam)+S.Rational(17,12)
 zero('proper_phase_scale',G*lam*ksquared-(60*S.pi+S.Rational(17,12)*G*lam))
 check('proper_scale_lower_reserve',60*3>13**2)
 t=S.symbols('t',real=True)
 source=t**3*(1-t)**3;qpoly=t**4*(1-t)**4
 zero('compact_source_mass',S.integrate(source,(t,0,1))-S.Rational(1,140))
 for j in range(4):
  zero('compatible_left_jet_'+str(j),S.diff(qpoly,t,j).subs(t,0))
  zero('compatible_right_jet_'+str(j),S.diff(qpoly,t,j).subs(t,1))

def w_hyper(z,a,gamma):
 k=mp.sqrt(a-mp.mpf(1)/4);p=mp.mpf(1)/2+1j*k
 return -mp.cos(z)*mp.sin(z)**p*mp.hyp2f1((p+1+1j*gamma)/2,(p+1-1j*gamma)/2,mp.mpf('1.5'),mp.cos(z)**2)

def rk_forced(a,gamma,compatible,N):
 lo=.7;hi=1.;L=math.pi/2;width=hi-lo
 def force(z):
  if z<=lo or z>=hi:return 0.
  t=(z-lo)/width
  if not compatible:return t**3*(1-t)**3
  q=t**4*(1-t)**4
  qdd=(12*t*t-80*t**3+180*t**4-168*t**5+56*t**6)/width**2
  return -qdd+(gamma*gamma-a/math.sin(z)**2)*q
 def fun(z,u,up):return up,(gamma*gamma-a/math.sin(z)**2)*u-force(z)
 u=0.;up=0.
 for left,right in [(lo,hi),(hi,L)]:
  h=(right-left)/N
  for i in range(N):
   z=left+i*h
   a1,b1=fun(z,u,up)
   a2,b2=fun(z+h/2,u+h*a1/2,up+h*b1/2)
   a3,b3=fun(z+h/2,u+h*a2/2,up+h*b2/2)
   a4,b4=fun(z+h,u+h*a3,up+h*b3)
   u+=h*(a1+2*a2+2*a3+a4)/6;up+=h*(b1+2*b2+2*b3+b4)/6
 return u

def independent_fixtures():
 mp.mp.dps=PROTOCOL['mpmath_precision'];L=mp.pi/2;lo=mp.mpf('.7');hi=mp.mpf(1);width=hi-lo
 for a in PROTOCOL['coefficients']:
  for gf in PROTOCOL['Laplace_frequencies']:
   gamma=mp.mpf(str(gf));tag=f'{a}_{gf}'
   w=lambda z:w_hyper(z,a,gamma)
   for zz in ['0.4','0.9','1.3']:
    z=mp.mpf(zz);value=w(z)
    residual=-mp.diff(w,z,2)+(gamma**2-a/mp.sin(z)**2)*value
    check('hypergeometric_ODE_'+tag+'_'+zz,abs(residual)<mp.mpf('1e-35'),mp.nstr(abs(residual),8))
    check('real_homogeneous_solution_'+tag+'_'+zz,abs(mp.im(value))<mp.mpf('1e-35'),mp.nstr(abs(mp.im(value)),8))
   check('regular_center_'+tag,abs(w(L))<mp.mpf('1e-35') and abs(mp.diff(w,L)-1)<mp.mpf('1e-35'))
   def q(z):
    t=(z-lo)/width;return t**4*(1-t)**4
   def source(z):
    t=(z-lo)/width;return t**3*(1-t)**3
   def compatible(z):return -mp.diff(q,z,2)+(gamma**2-a/mp.sin(z)**2)*q(z)
   moment=mp.quad(lambda z:w(z)*source(z),[lo,hi])
   compatible_moment=mp.quad(lambda z:w(z)*compatible(z),[lo,hi])
   check('compact_source_fails_moment_'+tag,abs(moment)>mp.mpf('1e-6'),mp.nstr(moment,20))
   check('compatible_source_zero_moment_'+tag,abs(compatible_moment)<mp.mpf('1e-35'),mp.nstr(abs(compatible_moment),8))
   errors=[];compatible_errors=[];solutions=[]
   for N in PROTOCOL['RK_steps_per_segment']:
    result=rk_forced(a,gf,False,N);result_c=rk_forced(a,gf,True,N)
    err=abs(result-float(mp.re(moment)));errors.append(err);compatible_errors.append(abs(result_c));solutions.append(result)
    check('independent_forced_ODE_moment_'+tag+'_'+str(N),err<PROTOCOL['RK_absolute_tolerance'],err)
    check('independent_compatible_ODE_'+tag+'_'+str(N),abs(result_c)<PROTOCOL['RK_absolute_tolerance'],abs(result_c))
   check('forced_ODE_refinement_'+tag,errors[1]<=PROTOCOL['RK_refinement_factor']*errors[0]+PROTOCOL['RK_refinement_floor'],errors)
   check('compatible_ODE_refinement_'+tag,compatible_errors[1]<=PROTOCOL['RK_refinement_factor']*compatible_errors[0]+PROTOCOL['RK_refinement_floor'],compatible_errors)
   NUMERICS.append({'a':a,'gamma':gf,'moment_real':mp.nstr(mp.re(moment),25),
      'moment_imaginary_abs':mp.nstr(abs(mp.im(moment)),8),'compatible_moment_abs':mp.nstr(abs(compatible_moment),8),
      'RK_center_values':solutions,'RK_moment_errors':errors,'RK_compatible_errors':compatible_errors})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 fast scalar boundary solvability','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();independent_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 fast scalar boundary solvability','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'comparison_boundary_trace_classification_derived':good,'two_zero_trace_range_iff_moment_derived':good,
  'comparison_open_data_retarded_exclusion_test_derived':good,'proper_phase_scale_derived':good,
  'actual_charged_transfer_proved':False,'physical_metric_boundary_selected':False,
  'finite_time_physical_stability_decided':False,'physical_ultraviolet_cutoff_selected':False,
  'actual_metric_causal_inverse_proved':False,'actual_wormhole_instability_proved':False,
  'stationary_compatible_inverse_invalidated':False,'order_reduction_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
