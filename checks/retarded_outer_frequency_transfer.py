"""S168: continuous positive-imaginary-frequency outer core-form transfer.

Ordinary continuum proof with independent spectral fixtures, not a mode search.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='retarded-outer-frequency-transfer'
OWN=['checks/retarded_outer_frequency_transfer.py','notes/'+STEM+'.md','evidence/retarded-outer-source-review.json']
PREVIOUS='fast-scalar-boundary-solvability'
ROWS=[];NUMERICS=[]
PROTOCOL={
 'scope':'Four actual nonzero thermal samples transfer the outer core-form limit to continuous positive imaginary frequencies',
 'matrix_measure':'exp(-t) [[1,t],[t,t^2]] dt on [0,infinity)',
 'matrix_thermal_indices':[2,7,31],'fractions':['0.2','0.5','0.9'],'thermal_spacing':'0.7',
 'mpmath_precision':60,'quadrature_identity_tolerance':'1e-40','matrix_order_tolerance':'1e-40',
 'scaled_d':['1/20','1/40','1/80','1/160'],'scaled_kappa':'1.3','scaled_offset':'0.37',
 'scaled_refinement_factor':'0.65','scaled_refinement_floor':'1e-40',
 'fixtures_are_actual_exterior_spectrum':False,'formal_verification':False,
}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,expr):
 result=S.simplify(S.sympify(expr).doit());check(name,result==0,result)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior()
 check('predecessor_passed',p['status']=='passed' and len(p['checks'])==126 and all(r['passed'] for r in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 required=['thermal-shape-response-and-subtractions','local-coupled-frequency-polynomial',
  'uniform-logarithmic-frequency-comparison','outer-coupled-scalar-separation',
  'curvature-free-all-frequency-bound','neutral-spherical-current-decoupling','spherical-retarded-constraint-recovery']
 check('required_inputs_pinned',all('notes/'+f+'.md' in p['input_sha256_after'] for f in required))
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('primary_source_scope',len(review['primary_sources'])==1 and not review['actual_metric_causal_inverse_proved']
  and not review['all_positive_frequencies_uniformly_controlled'] and not review['off_axis_complex_frequency_transfer_proved']
  and not review['order_reduction_adopted'])

def weights(x,nodes):
 ans=[]
 for j,a in enumerate(nodes):
  numerator=1;denominator=1
  for i,b in enumerate(nodes):
   if i!=j:numerator*=x-b;denominator*=a-b
  ans.append(numerator/denominator)
 return ans

def exact_checks():
 x,t,a,b,c=S.symbols('x t a b c',positive=True)
 nodes=[a,b,c];ws=weights(x,nodes)
 for j in range(3):zero('quadratic_contact_reproduction_'+str(j),sum(w*z**j for w,z in zip(ws,nodes))-x**j)
 f=lambda z:-1/(t+z)
 product=(x-a)*(x-b)*(x-c)
 zero('exact_four_denominator_identity',f(x)-sum(w*f(z) for w,z in zip(ws,nodes))-product/((t+x)*(t+a)*(t+b)*(t+c)))
 zero('third_resolvent_derivative',S.diff(f(x),x,3)-6/(t+x)**4)
 n,h,theta=S.symbols('n h theta',positive=True,real=True)
 bn=lambda j:h*h*j*j
 x=bn(n)+(bn(n+1)-bn(n))*theta
 lower=[bn(n),bn(n+1),bn(n+2)];upper=[bn(n-1),bn(n),bn(n+1)]
 lw=[S.factor(w) for w in weights(x,lower)];uw=[S.factor(w) for w in weights(x,upper)]
 # These factorizations prove the signs for n>=2, 0<=theta<=1.
 expected_lower=[(1-theta)*(4*n+4-(2*n+1)*theta)/(4*n+4),
  theta*(4*n+4-(2*n+1)*theta)/(2*n+3),
  -(2*n+1)**2*theta*(1-theta)/((4*n+4)*(2*n+3))]
 expected_upper=[-(2*n+1)**2*theta*(1-theta)/(4*n*(2*n-1)),
  ((2*n-1)+(2*n+1)*theta)*(1-theta)/(2*n-1),
  ((2*n-1)+(2*n+1)*theta)*theta/(4*n)]
 for name,actual,expected in [('lower',lw,expected_lower),('upper',uw,expected_upper)]:
  for j in range(3):zero(name+'_signed_weight_'+str(j),actual[j]-expected[j])
  zero(name+'_weight_sum',sum(actual)-1)
  check(name+'_left_endpoint_weights',[S.simplify(w.subs(theta,0)) for w in actual]==([1,0,0] if name=='lower' else [0,1,0]))
  check(name+'_right_endpoint_weights',[S.simplify(w.subs(theta,1)) for w in actual]==([0,1,0] if name=='lower' else [0,0,1]))
 zero('theta_product_maximum',S.Rational(1,4)-theta*(1-theta)-(theta-S.Rational(1,2))**2)
 lowermax=(2*n+1)**2/(4*(4*n+4)*(2*n+3))
 uppermax=(2*n+1)**2/(16*n*(2*n-1))
 zero('lower_negative_weight_slack',S.Rational(1,8)-lowermax-(12*n+10)/(8*(4*n+4)*(2*n+3)))
 zero('upper_negative_weight_slack',S.Rational(25,96)-uppermax-(n-2)*(26*n+3)/(96*n*(2*n-1)))
 zero('lower_absolute_sum',lw[0]+lw[1]-lw[2]-(1-2*lw[2]))
 zero('upper_absolute_sum',-uw[0]+uw[1]+uw[2]-(1-2*uw[0]))
 zero('lower_absolute_sum_bound',1+2*S.Rational(1,8)-S.Rational(5,4))
 zero('upper_absolute_sum_bound',1+2*S.Rational(25,96)-S.Rational(73,48))
 zero('positive_lower_remainder_product',S.prod([x-z for z in lower])
  -h**6*(2*n+1)**2*theta*(1-theta)*(4*n+4-(2*n+1)*theta))
 zero('negative_upper_remainder_product',S.prod([x-z for z in upper])
  +h**6*(2*n+1)**2*theta*(1-theta)*(2*n-1+(2*n+1)*theta))
 check('all_anchors_positive_for_n_ge_2',bn(n-1).subs(n,2)>0)
 check('n_one_would_require_forbidden_static_sample',bn(n-1).subs(n,1)==0)
 wrong_lower=(-lw[2]).subs({n:2,theta:S.Rational(1,2)})
 check('negative_weight_detects_lower_anchor_only_error',wrong_lower>0,wrong_lower)
 # Actual D is zero at all anchors, but the valid anchor lower bounds [0,0,-1]
 # would falsely claim the above positive lower bound if signed weights were ignored.
 r=S.symbols('r',nonnegative=True)
 offsets=[-1-r,-r,1-r,2-r]
 check('four_anchor_scaled_spacing_bound',all(max(abs(o.subs(r,0)),abs(o.subs(r,1)))<=2 for o in offsets))
 u=S.Matrix([1,t]);M=u*u.T
 zero('matrix_measure_determinant',M.det())
 check('matrix_measure_positive_trace',S.Poly(S.trace(M),t).all_coeffs()==[1,0,1])
 comm=M.subs(t,1)*M.subs(t,2)-M.subs(t,2)*M.subs(t,1)
 check('matrix_measure_not_simultaneously_diagonal',comm!=S.zeros(2),comm)
 d,ss=S.symbols('d ss',positive=True)
 zero('scaled_atomic_measure_weight',-d**(-6)/(d**(-2)+ss)+d**(-4)/(1+d*d*ss))

def maxabs(A):return max(abs(A[i,j]) for i in range(A.rows) for j in range(A.cols))
def eigs(A):return mp.eigsy(A,eigvals_only=True)
def opnorm(A):return max(abs(v) for v in eigs(A))
def interpolate(x,nodes,values):
 out=mp.zeros(2)
 for w,M in zip(weights(x,nodes),values):out+=w*M
 return out

def continuous_matrix_fixtures():
 mp.mp.dps=PROTOCOL['mpmath_precision'];tol=mp.mpf(PROTOCOL['quadrature_identity_tolerance'])
 ordertol=mp.mpf(PROTOCOL['matrix_order_tolerance']);h=mp.mpf(PROTOCOL['thermal_spacing'])
 A0=mp.matrix([[1,2],[2,-1]]);A1=mp.matrix([[2,-1],[-1,-3]]);A2=mp.matrix([[-1,mp.mpf('.3')],[mp.mpf('.3'),2]])
 poly=lambda x:A0+x*A1+x*x*A2
 def integral(x,nodes=()):
  entries=[]
  for power in range(3):
   def integrand(t):
    den=t+x
    for b in nodes:den*=t+b
    return mp.exp(-t)*t**power/den
   entries.append(mp.quad(integrand,[0,1,4,mp.inf]))
  return mp.matrix([[entries[0],entries[1]],[entries[1],entries[2]]])
 for n in PROTOCOL['matrix_thermal_indices']:
  anchors=[h*h*j*j for j in range(n-1,n+3)]
  values=[poly(a)-integral(a) for a in anchors]
  for sf in PROTOCOL['fractions']:
   theta=mp.mpf(sf);x=anchors[1]+theta*(anchors[2]-anchors[1]);D=poly(x)-integral(x)
   low=interpolate(x,anchors[1:],values[1:]);up=interpolate(x,anchors[:3],values[:3]);tag=f'{n}_{sf}'
   rlow=mp.fprod(x-a for a in anchors[1:])*integral(x,anchors[1:])
   rup=mp.fprod(x-a for a in anchors[:3])*integral(x,anchors[:3])
   errlo=maxabs(D-low-rlow);errup=maxabs(D-up-rup)
   check('independent_continuous_lower_remainder_'+tag,errlo<tol,mp.nstr(errlo,8))
   check('independent_continuous_upper_remainder_'+tag,errup<tol,mp.nstr(errup,8))
   minlo=min(eigs(D-low));minup=min(eigs(up-D))
   check('matrix_Loewner_enclosure_'+tag,minlo>=-ordertol and minup>=-ordertol,[mp.nstr(minlo,12),mp.nstr(minup,12)])
   target=mp.matrix([[mp.mpf('.4'),mp.mpf('.2')],[mp.mpf('.2'),mp.mpf('-.7')]])
   E=max(opnorm(v-target) for v in values);error=opnorm(D-target)
   check('matrix_error_transfer_'+tag,error<=mp.mpf(73)/48*E+ordertol,mp.nstr(error/E,12))
   NUMERICS.append({'fixture':'continuous_matrix','n':n,'fraction':sf,'lower_identity_error':mp.nstr(errlo,12),
    'upper_identity_error':mp.nstr(errup,12),'lower_gap_min_eigenvalue':mp.nstr(minlo,18),
    'upper_gap_min_eigenvalue':mp.nstr(minup,18),'target_error_over_max_anchor_error':mp.nstr(error/E,18)})

def scaled_fixtures():
 h=mp.mpf(PROTOCOL['thermal_spacing']);kappa=mp.mpf(PROTOCOL['scaled_kappa']);offset=mp.mpf(PROTOCOL['scaled_offset'])
 tol=mp.mpf(PROTOCOL['matrix_order_tolerance'])
 H0=mp.matrix([[1,mp.mpf('.3')],[mp.mpf('.3'),-2]])
 H1=mp.matrix([[mp.mpf('.4'),mp.mpf('-.2')],[mp.mpf('-.2'),mp.mpf('.7')]])
 H2=mp.matrix([[mp.mpf('.1'),mp.mpf('.05')],[mp.mpf('.05'),mp.mpf('-.03')]])
 P=mp.matrix([[1,mp.mpf('.4')],[mp.mpf('.4'),2]]);R=mp.matrix([[mp.mpf('.2'),mp.mpf('-.1')],[mp.mpf('-.1'),mp.mpf('.3')]])
 L=lambda y:H0+y*H1+y*y*H2-P/(1+y)
 target=L(kappa*kappa);previous=None;previous_E=None
 check('scaled_positive_atomic_matrix',min(eigs(P))>0)
 for ds in PROTOCOL['scaled_d']:
  numerator,denominator=ds.split('/');d=mp.mpf(numerator)/mp.mpf(denominator)
  gamma=kappa/d+offset*h;n=int(mp.floor(gamma/h));x=gamma*gamma
  anchors=[h*h*j*j for j in range(n-1,n+3)]
  # Direct polynomial plus a positive atomic spectral measure at t=d^-2.
  def direct(sigma):return d**(-4)*H0+d**(-2)*sigma*H1+sigma*sigma*H2+d**(-2)*R-d**(-6)*P/(d**(-2)+sigma)
  scaled=d**4*direct(x);values=[d**4*direct(a) for a in anchors]
  low=interpolate(x,anchors[1:],values[1:]);up=interpolate(x,anchors[:3],values[:3])
  error=opnorm(scaled-target);E=max(opnorm(v-target) for v in values)
  check('scaled_resolvent_equals_profile_'+ds,maxabs(scaled-L(d*d*x)-d*d*R)<tol)
  check('scaled_four_anchors_near_target_'+ds,max(abs(d*mp.sqrt(a)-d*gamma) for a in anchors)<=2*d*h+tol)
  check('scaled_Loewner_enclosure_'+ds,min(eigs(scaled-low))>=-tol and min(eigs(up-scaled))>=-tol)
  check('scaled_error_transfer_'+ds,error<=mp.mpf(73)/48*E+tol,mp.nstr(error/E,16))
  if previous is not None:
   factor=mp.mpf(PROTOCOL['scaled_refinement_factor']);floor=mp.mpf(PROTOCOL['scaled_refinement_floor'])
   check('scaled_actual_error_refines_'+ds,error<factor*previous+floor,mp.nstr(error/previous,12))
   check('scaled_anchor_error_refines_'+ds,E<factor*previous_E+floor,mp.nstr(E/previous_E,12))
  NUMERICS.append({'fixture':'scaled_atomic_matrix','d':ds,'n':n,'actual_error':mp.nstr(error,20),
   'max_anchor_error':mp.nstr(E,20),'error_ratio':mp.nstr(error/E,18)})
  previous=error;previous_E=E

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 retarded outer frequency transfer','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();continuous_matrix_fixtures();scaled_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 retarded outer frequency transfer','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'four_sample_causal_enclosure_derived':good,'uniform_signed_weight_bound_derived':good,
  'actual_seed_retarded_outer_core_limit_derived':good,'fixed_finite_profile_matrix_limit_derived':good,
  'all_local_contacts_retained':good,'actual_metric_causal_inverse_proved':False,
  'actual_metric_mode_constructed':False,'physical_metric_boundary_selected':False,
  'all_positive_frequencies_uniformly_controlled':False,'off_axis_complex_frequency_transfer_proved':False,
  'actual_charged_stability_proved':False,'actual_wormhole_instability_proved':False,
  'retained_fast_factor_removed':False,'order_reduction_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
