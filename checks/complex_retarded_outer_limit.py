"""S168 complex outer core limit: exact identities and spectral fixtures.

The analytic convergence and weighted weak-limit proofs are in the note.
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
STEM='complex-retarded-outer-limit'
OWN=['checks/complex_retarded_outer_limit.py','notes/'+STEM+'.md','evidence/complex-outer-source-review.json']
PREVIOUS='retarded-outer-frequency-transfer'
ROWS=[];NUMERICS=[]
PROTOCOL={
 'scope':'Complex core-form convergence and weighted scaled spectral measure for fixed profiles',
 'mpmath_precision':65,'identity_tolerance':'1e-42',
 'complex_points':[['0.4','0.9'],['-0.5','0.2'],['-2','0.4'],['4','-1.2']],
 'continuous_measure':'exp(-t) [[1,t],[t,t^2]] dt',
 'positive_anchors':[1,2,3,4],'log_spectral_atoms':['0.25','1','4'],
 'contour_center':['0.6','0.8'],'contour_radius':'0.15','contour_nodes':[32,64],
 'contour_absolute_tolerances':['1e-22','1e-45'],'contour_refinement_factor':'1e-8',
 'scaled_d':['1/10','1/20','1/40','1/80'],'scaled_refinement_factor':'0.65',
 'fixtures_are_exterior_spectrum':False,'formal_verification':False,
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
 check('predecessor_passed',p['status']=='passed' and len(p['checks'])==98 and all(r['passed'] for r in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 required=['thermal-shape-response-and-subtractions','local-coupled-frequency-polynomial',
  'uniform-logarithmic-frequency-comparison','outer-coupled-scalar-separation',
  'curvature-free-all-frequency-bound','neutral-spherical-current-decoupling','spherical-retarded-constraint-recovery']
 check('required_inputs_pinned',all('notes/'+f+'.md' in p['input_sha256_after'] for f in required))
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('source_scope',len(review['primary_sources'])==1 and not review['actual_metric_causal_inverse_proved']
  and not review['physical_metric_boundary_selected'] and not review['uniform_infinite_profile_space_control_proved']
  and not review['spectral_cut_boundary_values_controlled'])

def weights(z):return [(z-2)*(z-3)/2,-(z-1)*(z-3),(z-1)*(z-2)/2]
def P3(z):return (z-1)*(z-2)*(z-3)
def W(t):return (t+1)*(t+2)*(t+3)*(t+4)

def exact_checks():
 z=S.symbols('z');t,d,A,r=S.symbols('t d A r',positive=True)
 for k in range(3):zero('all_quadratic_contacts_'+str(k),sum(w*j**k for w,j in zip(weights(z),[1,2,3]))-z**k)
 f=lambda x:-1/(t+x)
 zero('complex_resolvent_remainder',f(z)-sum(w*f(j) for w,j in zip(weights(z),[1,2,3]))-P3(z)/((t+z)*(t+1)*(t+2)*(t+3)))
 zero('four_anchor_divided_difference',(f(4)-3*f(3)+3*f(2)-f(1))/6-1/W(t))
 zero('scaled_measure_mass_factor',-d**4/(t+z/d**2)+d**6/(d**2*t+z))
 zero('weighted_transform_identity',1/((t+z)*(t+1)*(t+2)*(t+3))-(t+4)/((t+z)*W(t)))
 c,k=S.symbols('c k',real=True)
 zero('sector_denominator_slack',t*t+r*r+2*t*r*k-c*c*(t+r)**2
  -((1-c*c)*(t-r)**2+2*t*r*(k-(2*c*c-1))))
 u=S.symbols('u',nonnegative=True)
 zero('compactification_transform',(t+4)/(t+z)-1/(1+(z-4)/(t+4)))
 for j in range(5):
  zero('compact_moment_derivative_'+str(j),S.diff((t+4)/(t+z),z,j).subs(z,4)-(-1)**j*S.factorial(j)/(t+4)**j)
 logf=(r+z)**2*S.log(r+z)/2
 zero('log_multiplier_third_derivative',S.diff(logf,z,3)-1/(r+z))
 zero('density_integral_normalization',S.integrate(u**2/(u+A)**4,(u,0,S.oo))-1/(3*A))
 zero('spectral_density_coefficient',6/(60*S.pi)/(3*A)-1/(30*S.pi*A))
 x,y=S.symbols('x y',real=True)
 bad=(x*x+y*y)*S.log(x*x+y*y)/4
 dbar=(S.diff(bad,x)+S.I*S.diff(bad,y))/2
 check('modulus_square_is_not_holomorphic',S.simplify(dbar.subs({x:1,y:1}))!=0)
 M=S.Matrix([[z-1,1],[1,0]])
 zero('full_matrix_remains_invertible',M.det()+1)
 check('compressed_zero_is_not_full_pole',M[0,0].subs(z,1)==0 and M.subs(z,1).det()==-1)

def fro(M):return mp.sqrt(sum(abs(x)**2 for x in M))
def maxabs(M):return max(abs(x) for x in M)
def matquad(fn,breaks):return mp.matrix([[mp.quad(lambda t:fn(t)[i,j],breaks) for j in range(2)] for i in range(2)])
def interp(z,values):
 result=mp.zeros(2)
 for w,v in zip(weights(z),values[:3]):result+=w*v
 return result
def complex_points():return [mp.mpc(a,b) for a,b in PROTOCOL['complex_points']]

def continuous_fixtures():
 mp.mp.dps=PROTOCOL['mpmath_precision'];tol=mp.mpf(PROTOCOL['identity_tolerance'])
 A0=mp.matrix([[1,2],[2,-1]]);A1=mp.matrix([[2,-1],[-1,-3]]);A2=mp.matrix([[-1,mp.mpf('.3')],[mp.mpf('.3'),2]])
 poly=lambda z:A0+z*A1+z*z*A2
 measure=lambda t:mp.exp(-t)*mp.matrix([[1,t],[t,t*t]])
 breaks=[0,mp.mpf('.5'),1,2,4,8,mp.inf]
 def closed(z):
  q=mp.exp(z)*mp.e1(z)
  return poly(z)-mp.matrix([[q,1-z*q],[1-z*q,1-z+z*z*q]])
 anchors=[closed(mp.mpf(j)) for j in [1,2,3,4]]
 G=(anchors[3]-3*anchors[2]+3*anchors[1]-anchors[0])/6
 Gi=matquad(lambda t:measure(t)/W(t),breaks)
 check('continuous_positive_four_sample_matrix',min(mp.eigsy(G,eigvals_only=True))>0)
 check('continuous_sample_mass_equals_integral',maxabs(G-Gi)<tol,mp.nstr(maxabs(G-Gi),10))
 for ix,z in enumerate(complex_points()):
  direct=poly(z)-matquad(lambda t:measure(t)/(t+z),breaks)
  remainder=P3(z)*matquad(lambda t:measure(t)/((t+z)*(t+1)*(t+2)*(t+3)),breaks)
  err=maxabs(direct-interp(z,anchors)-remainder)
  check('complex_continuous_closed_form_'+str(ix),maxabs(direct-closed(z))<tol)
  check('complex_continuous_remainder_'+str(ix),err<tol,mp.nstr(err,10))
  third=mp.matrix([[mp.diff(lambda zz:closed(zz)[i,j],z,3) for j in range(2)] for i in range(2)])
  thi=6*matquad(lambda t:measure(t)/(t+z)**4,breaks)
  derr=maxabs(third-thi)
  check('complex_continuous_third_derivative_'+str(ix),derr<tol,mp.nstr(derr,10))
  distance=abs(z) if z.real>=0 else abs(z.imag)
  boundM=1+abs(4-z)/distance
  bound=sum(abs(w)*fro(v) for w,v in zip(weights(z),anchors[:3]))+abs(P3(z))*boundM*(G[0,0]+G[1,1])
  check('complex_matrix_local_bound_'+str(ix),fro(direct)<=bound,mp.nstr(fro(direct)/bound,12))
  theta=mp.pi-abs(mp.arg(z));sectorM=max(1,4/abs(z))/mp.sin(theta/2)
  probe_ts=[mp.mpf(0),mp.mpf('.1'),abs(z),mp.mpf(4),mp.mpf(100)]
  check('sector_ratio_bound_'+str(ix),all((t+4)/abs(t+z)<=sectorM+tol for t in probe_ts))
  NUMERICS.append({'fixture':'continuous_matrix','z':[str(z.real),str(z.imag)],'remainder_error':mp.nstr(err,12),
   'third_derivative_error':mp.nstr(derr,12),'bound_ratio':mp.nstr(fro(direct)/bound,18)})

def logarithmic_fixtures():
 tol=mp.mpf(PROTOCOL['identity_tolerance']);c=1/(30*mp.pi)
 rs=[mp.mpf(r) for r in PROTOCOL['log_spectral_atoms']]
 Ms=[mp.matrix([[1,mp.mpf('.2')],[mp.mpf('.2'),mp.mpf('.5')]]),
     mp.matrix([[mp.mpf('.3'),mp.mpf('-.1')],[mp.mpf('-.1'),1]]),
     mp.matrix([[mp.mpf('.6'),mp.mpf('.3')],[mp.mpf('.3'),mp.mpf('.8')]])]
 A0=mp.matrix([[1,mp.mpf('.4')],[mp.mpf('.4'),-2]]);A1=mp.matrix([[mp.mpf('.2'),-1],[-1,mp.mpf('.3')]])
 A2=mp.matrix([[mp.mpf('-.1'),mp.mpf('.2')],[mp.mpf('.2'),mp.mpf('.4')]])
 def full(z):
  out=A0+z*A1+z*z*A2
  for r,M in zip(rs,Ms):out+=c*(r+z)**2*mp.log(r+z)*M/2
  return out
 def derivative(z,m):
  out=A1+2*z*A2 if m==1 else 2*A2 if m==2 else mp.zeros(2)
  for r,M in zip(rs,Ms):
   a=r+z;v=a*mp.log(a)+a/2 if m==1 else mp.log(a)+mp.mpf('1.5') if m==2 else 1/a
   out+=c*v*M
  return out
 def integrate_density(kernel):
  out=mp.zeros(2)
  for r,M in zip(rs,Ms):
   # Integrate the independent density (t-r)^2/(60pi), t>=r.
   val=mp.quad(lambda u:u*u*kernel(r+u),[0,1,4,mp.inf])/(60*mp.pi)
   out+=val*M
  return out
 check('normal_log_matrix_weights_positive',all(min(mp.eigsy(M,eigvals_only=True))>0 for M in Ms))
 check('normal_log_matrix_weights_noncommute',maxabs(Ms[0]*Ms[1]-Ms[1]*Ms[0])>mp.mpf('.01'))
 anchors=[full(mp.mpf(j)) for j in [1,2,3,4]]
 G=(anchors[3]-3*anchors[2]+3*anchors[1]-anchors[0])/6
 check('normal_density_weighted_mass',maxabs(G-integrate_density(lambda t:1/W(t)))<tol)
 for ix,z in enumerate(complex_points()):
  rem=P3(z)*integrate_density(lambda t:1/((t+z)*(t+1)*(t+2)*(t+3)))
  err=maxabs(full(z)-interp(z,anchors)-rem)
  derr=maxabs(derivative(z,3)-6*integrate_density(lambda t:1/(t+z)**4))
  check('normal_density_complex_remainder_'+str(ix),err<tol,mp.nstr(err,10))
  check('normal_density_complex_derivative_'+str(ix),derr<tol,mp.nstr(derr,10))
  NUMERICS.append({'fixture':'normal_log_density','z':[str(z.real),str(z.imag)],'remainder_error':mp.nstr(err,12),
    'third_derivative_error':mp.nstr(derr,12)})
 center=mp.mpc(*PROTOCOL['contour_center']);radius=mp.mpf(PROTOCOL['contour_radius']);previous=None
 for N,tols in zip(PROTOCOL['contour_nodes'],PROTOCOL['contour_absolute_tolerances']):
  vals=[full(center+radius*mp.exp(2j*mp.pi*j/N)) for j in range(N)];errors=[]
  for m in [1,2,3]:
   approx=mp.zeros(2)
   for j,val in enumerate(vals):approx+=val*mp.exp(-2j*mp.pi*j*m/N)
   approx*=mp.factorial(m)/(N*radius**m)
   error=maxabs(approx-derivative(center,m));errors.append(error)
   check('independent_Cauchy_derivative_'+str(N)+'_'+str(m),error<mp.mpf(tols),mp.nstr(error,12))
  if previous is not None:check('Cauchy_contour_refinement',max(errors)<mp.mpf(PROTOCOL['contour_refinement_factor'])*max(previous))
  NUMERICS.append({'fixture':'Cauchy_contour','nodes':N,'derivative_errors':[mp.nstr(e,16) for e in errors]});previous=errors

def scaled_fixtures_and_controls():
 tol=mp.mpf(PROTOCOL['identity_tolerance']);P=mp.matrix([[1,mp.mpf('.3')],[mp.mpf('.3'),2]])
 Q=mp.matrix([[mp.mpf('.5'),mp.mpf('-.2')],[mp.mpf('-.2'),1]])
 poly=lambda z:mp.matrix([[z*z-1,z+mp.mpf('.2')],[z+mp.mpf('.2'),2-z*z/3]])
 limit=lambda z:poly(z)-P/(2+z)
 previous=None;previous_moment=None
 for ds in PROTOCOL['scaled_d']:
  a,b=ds.split('/');d=mp.mpf(a)/mp.mpf(b);T=d**(-2)
  def full(z):return poly(z)-(1+d)*P/(2+d+z)+d*P3(z)*(T+4)/(T+z)*Q
  errors=[fro(full(z)-limit(z)) for z in complex_points()]
  moment_errors=[]
  for j in range(4):
   moment=(1+d)*P/(W(2+d)*(6+d)**j)+d*Q/(T+4)**j
   target=P/(W(2)*6**j);moment_errors.append(fro(moment-target))
  vals=[full(mp.mpf(j)) for j in [1,2,3,4]]
  G=(vals[3]-3*vals[2]+3*vals[1]-vals[0])/6
  check('scaled_moving_measure_mass_'+ds,maxabs(G-((1+d)*P/W(2+d)+d*Q))<tol)
  for ix,z in enumerate(complex_points()):
   integral=(1+d)*P/((2+d+z)*(3+d)*(4+d)*(5+d))+d*(T+4)/(T+z)*Q
   check('scaled_positive_measure_remainder_'+ds+'_'+str(ix),maxabs(full(z)-interp(z,vals)-P3(z)*integral)<tol)
  if previous is not None:
   factor=mp.mpf(PROTOCOL['scaled_refinement_factor'])
   check('scaled_complex_error_refinement_'+ds,max(errors)<factor*previous)
   check('scaled_weak_moment_refinement_'+ds,max(moment_errors)<factor*previous_moment)
  NUMERICS.append({'fixture':'scaled_positive_measure','d':ds,'max_complex_error':mp.nstr(max(errors),20),
   'weighted_moment_errors':[mp.nstr(e,16) for e in moment_errors],'weighted_far_mass_scalar':str(d),
   'unweighted_far_mass_scalar':mp.nstr(d*W(T),20)})
  previous=max(errors);previous_moment=max(moment_errors)
 n=S.symbols('n',positive=True,integer=True);z=S.symbols('z')
 zero('analytic_only_control_negative_third_derivative',S.diff(S.exp(-n*z),z,3)+n**3*S.exp(-n*z))
 check('positive_axis_only_control_fails_complex_bound',abs(mp.exp(-8*mp.mpc(-1,1)))>mp.exp(7)
  and mp.exp(-8)<mp.exp(-4))

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 complex retarded outer limit','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();continuous_fixtures();logarithmic_fixtures();scaled_fixtures_and_controls()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 complex retarded outer limit','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'complex_compact_response_bound_derived':good,'actual_complex_outer_core_limit_derived':good,
  'frequency_derivative_convergence_derived':good,'weighted_scaled_spectral_limit_derived':good,
  'all_local_contacts_retained':good,'actual_metric_causal_inverse_proved':False,
  'actual_metric_mode_constructed':False,'physical_metric_boundary_selected':False,
  'uniform_infinite_profile_space_control_proved':False,'spectral_cut_boundary_values_controlled':False,
  'quantitative_annular_threshold_proved':False,'actual_charged_stability_proved':False,
  'actual_wormhole_instability_proved':False,'retained_fast_factor_removed':False,'order_reduction_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
