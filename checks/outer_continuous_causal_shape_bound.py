"""Four thermal anchors, angular derivative bounds and continuous controls."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]; STEM='outer-continuous-causal-shape-bound'
PREVIOUS='coupled-boundary-phase-winding'
OWN=['checks/outer_continuous_causal_shape_bound.py','notes/'+STEM+'.md','evidence/outer-continuous-causal-source-review.json']
ROWS=[]; FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-29','precision_tolerance':'1e-29',
 'measure_cases':[[20,'.07','.2'],[37,'.11','.65']],'measure_columns':['n','omega_T','fraction_between_squared_nodes'],
 'positive_measure':'t^2 exp(-t) (1,t)(1,t)^T dt on [0,infinity)',
 'log_cases':[['2.3','1.7',['1','-.4']],['5.2','3.1',['.7','1.2']]],
 'log_A':[['2','.3'],['.3','1']],'log_B':[['1.1','.2'],['.2','.8']],
 'log_steps':['.002','.001'],'finite_difference_tolerance':'1e-6','refinement_ratio_ceiling':'.3',
 'tower_q':[3,7],'tower_C':5,'tower_cutoffs':[16,64],
 'horizon_exponents':['2.5','3.25'],'horizon_c':'.3','horizon_cutoffs':['.2','.1'],
 'actual_charged_kernel_evaluated':False,'global_causal_shape_bound_verified':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)}; ROWS.append(row); print(json.dumps(row),flush=True)
 if not ok: raise AssertionError(name)

def zero(name,value):
 value=S.simplify(S.sympify(value).doit()); check(name,value==0,value)

def prior(): return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths(): return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior(); check('predecessor_passed',p['status']=='passed' and len(p['checks'])==137 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['actual_outer_continuous_compensation_proved'] and review['outer_supported_full_causal_core_bound_proved']
       and not review['global_continuous_shape_coercivity_proved'] and not review['actual_continuous_Schur_family_verified'])

def lagrange(nodes,s,f):
 return sum(f(x)*S.prod((s-y)/(x-y) for y in nodes if y!=x) for x in nodes)

def exact_algebra():
 a,b,c,d,s,t=S.symbols('a b c d s t',positive=True); f=lambda x:-1/(t+x)
 up=lagrange([a,b,c],s,f); low=lagrange([b,c,d],s,f)
 zero('upper_resolvent_remainder',f(s)-up-(s-a)*(s-b)*(s-c)/((t+s)*(t+a)*(t+b)*(t+c)))
 zero('lower_resolvent_remainder',f(s)-low-(s-b)*(s-c)*(s-d)/((t+s)*(t+b)*(t+c)*(t+d)))
 zero('bracket_width_identity',up-low-(d-a)*(s-b)*(c-s)/((t+a)*(t+b)*(t+c)*(t+d)))
 for power in [0,1,2]:
  for tag,nodes in [('upper',[a,b,c]),('lower',[b,c,d])]: zero('local_polynomial_'+tag+'_'+str(power),lagrange(nodes,s,lambda x:x**power)-s**power)
 z=S.symbols('z',nonnegative=True); n=z+20
 budgets={'upper_Lebesgue':S.Rational(1,2)*8*n*(2*n-1)-(2*n+1)**2,
  'lower_Lebesgue':S.Rational(1,2)*2*(4*n+4)*(2*n+3)-(2*n+1)**2,
  'node_relative_displacement':n*n/4-4*n-4,
  'scalar_width':48*n**3-S.Rational(6,5)*(6*n+3)*(2*n+1)**2}
 for name,p in budgets.items():
  co=S.Poly(S.expand(p),z).all_coeffs(); check('uniform_'+name,all(v>=0 for v in co) and co[-1]>0,co)
 check('warp_ratio_bound',1/(1-S.Rational(16,100))<S.Rational(6,5))
 zero('anchor_fourth_budget',S.Rational(1,500)*S.Rational(3,2)*S.Rational(13,10)**2-S.Rational(507,100000))
 zero('combined_fourth_budget',S.Rational(507,100000)+S.Rational(1,1000)-S.Rational(607,100000))
 margin=S.Rational(27,600)/S.Rational(22,7)-S.Rational(607,100000)
 check('actual_outer_margin',margin>0,margin)
 zero('scalar_alpha_cubic_moment',S.integrate(t**3*(1-t)**3,(t,0,1))-S.Rational(1,140))
 alpha=S.symbols('alpha',real=True)
 logconstant=alpha*S.log(S.pi)+2*S.loggamma((1-alpha)/4)-2*S.loggamma((1+alpha)/4)
 zero('Pitt_logarithmic_Fourier_conversion',-S.diff(logconstant,alpha).subs(alpha,0)+S.log(2*S.pi)-S.polygamma(0,S.Rational(1,4))-S.log(2))
 check('retained_logarithmic_constant_lower_bound',-1-S.Rational(7,5)-S.Rational(11,7)>-4)
 rho,nu,cc=S.symbols('rho nu c',positive=True); p=rho**nu*(1+cc*rho**2)
 zero('noninteger_horizon_Laplacian',S.diff(p,rho,2)+S.diff(p,rho)/rho-nu**2*p/rho**2-4*cc*(nu+1)*rho**nu)

def tower_bounds():
 m,q,C=S.symbols('m q C',positive=True)
 zero('angular_derivative_telescoping',2*(q+2*m)*C**2/(m*m*(m+q)**2)-(2*C**2/q)*(1/m**2-1/(m+q)**2))
 for q in PROTOCOL['tower_q']:
  C=PROTOCOL['tower_C']; full=S.Rational(2*C*C,q)*sum(S.Rational(1,j*j) for j in range(1,q+1))
  for N in PROTOCOL['tower_cutoffs']:
   partial=sum(S.Rational(2*(q+2*j)*C*C,j*j*(j+q)**2) for j in range(1,N+1))
   tail=S.Rational(2*C*C,q)*sum(S.Rational(1,j*j) for j in range(N+1,N+q+1))
   zero(f'angular_tail_identity_{q}_{N}',partial+tail-full)
   bound=S.Rational(2*C*C,N*(N+q)); check(f'angular_tail_bound_{q}_{N}',0<tail<bound,str(tail/bound))
   check(f'complete_angular_bound_{q}_{N}',full<S.Rational(4*C*C,q))
   FIXTURES.append({'fixture':'complete_angular_tail','q':q,'C':C,'N':N,'tail':str(tail),'bound':str(bound),'ratio':str(S.N(tail/bound,40))})

def mp_interp(nodes,s,values):
 out=mp.zeros(2)
 for x,fx in zip(nodes,values): out+=fx*mp.fprod((s-y)/(x-y) for y in nodes if y!=x)
 return out

def eigmin2(M): return (M[0,0]+M[1,1]-mp.sqrt((M[0,0]-M[1,1])**2+4*M[0,1]**2))/2
def serial(M,digits): return [[mp.nstr(M[i,j],digits) for j in range(M.cols)] for i in range(M.rows)]
def parse(M): return mp.matrix([[mp.mpf(x) for x in row] for row in M])

def measure_brackets():
 for case,(n,ws,rs) in enumerate(PROTOCOL['measure_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    omega=mp.mpf(ws); frac=mp.mpf(rs); nodes=[omega**2*j*j for j in [n-1,n,n+1,n+2]]; a,b,c,d=nodes; s=b+frac*(c-b)
    P=lambda x:mp.matrix([[1-x/10+x*x/100,mp.mpf('.2')+x/10+x*x/50],[mp.mpf('.2')+x/10+x*x/50,2+3*x/10+7*x*x/100]])
    def spec(x): return mp.matrix([[mp.quad(lambda t:t**(2+i+j)*mp.exp(-t)/(t+x),[0,1,mp.inf]) for j in range(2)] for i in range(2)])
    def lap(x): return mp.matrix([[mp.quad(lambda u:mp.exp(-x*u)*mp.factorial(2+i+j)/(1+u)**(3+i+j),[0,1,mp.inf]) for j in range(2)] for i in range(2)])
    vals=[P(x)-spec(x) for x in nodes]; target=P(s)-lap(s); direct=P(s)-spec(s)
    err=mp.norm(target-direct)/(1+mp.norm(direct)); check(f'measure_independent_Laplace_{case}_{digits}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
    upper=mp_interp(nodes[:3],s,vals[:3]); lower=mp_interp(nodes[1:],s,vals[1:])
    lo=eigmin2(target-lower); hi=eigmin2(upper-target)
    check(f'actual_lower_bracket_{case}_{digits}',lo>0,mp.nstr(lo,digits)); check(f'actual_upper_bracket_{case}_{digits}',hi>0,mp.nstr(hi,digits))
    dd=mp.matrix([[mp.quad(lambda t:t**(2+i+j)*mp.exp(-t)/mp.fprod(t+x for x in nodes),[0,1,mp.inf]) for j in range(2)] for i in range(2)])
    width=(d-a)*(s-b)*(c-s)*dd; err=mp.norm(upper-lower-width)/(1+mp.norm(width))
    check(f'independent_bracket_width_{case}_{digits}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
    wrongvals=[P(x)+spec(x) for x in nodes]; wrong=P(s)+lap(s); wrongupper=mp_interp(nodes[:3],s,wrongvals[:3])
    check(f'wrong_bubble_sign_detected_{case}_{digits}',eigmin2(wrongupper-wrong)<0)
    levels.append({'precision':digits,'s':mp.nstr(s,digits),'actual':serial(target,digits),'lower_gap_min':mp.nstr(lo,digits),'upper_gap_min':mp.nstr(hi,digits)})
  with mp.workdps(80):
   aa=parse(levels[0]['actual']); bb=parse(levels[1]['actual']); err=mp.norm(aa-bb)/(1+mp.norm(bb))
   check('measure_precision_'+str(case),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'positive_matrix_measure','n':n,'omega_T':ws,'fraction':rs,'levels':levels,'charged_spectral_measure':False})

def log_derivatives():
 for case,(ss,ms,vs) in enumerate(PROTOCOL['log_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    A=parse(PROTOCOL['log_A']); B=parse(PROTOCOL['log_B']); vec=mp.matrix([mp.mpf(x) for x in vs]); s=mp.mpf(ss); mass2=mp.mpf(ms)
    check(f'noncommuting_log_fixture_{case}_{digits}',mp.norm(A*B-B*A)>mp.mpf('.01'))
    def phi(x):
     ev,U=mp.eigsy(A+x*B+mass2*mp.eye(2)); matrix=U*mp.diag([mp.log(v) for v in ev])*U.T
     return (vec.T*matrix*vec)[0]
    def integrand(t):
     R=(A+s*B+(mass2+t)*mp.eye(2))**-1
     return 6*(vec.T*R*B*R*B*R*B*R*vec)[0]
    actual=mp.quad(integrand,[0,1,mp.inf]); bmax=max(mp.eigsy(B,eigvals_only=True)); bound=2*bmax**3*mp.norm(vec)**2/mass2**3
    check(f'positive_log_third_derivative_{case}_{digits}',0<actual<=bound,mp.nstr(actual/bound,digits))
    records=[]; errors=[]
    for j,hs in enumerate(PROTOCOL['log_steps']):
     step=mp.mpf(hs); fd=(phi(s+2*step)-2*phi(s+step)+2*phi(s-step)-phi(s-2*step))/(2*step**3)
     err=abs(fd-actual)/(1+abs(actual)); errors.append(err)
     check(f'log_derivative_finite_difference_{case}_{digits}_{j}',err<mp.mpf(PROTOCOL['finite_difference_tolerance']),mp.nstr(err,digits))
     records.append({'step':hs,'finite_difference':mp.nstr(fd,digits),'scaled_error':mp.nstr(err,digits)})
    ratio=errors[1]/errors[0]; check(f'log_derivative_refinement_{case}_{digits}',ratio<mp.mpf(PROTOCOL['refinement_ratio_ceiling']),mp.nstr(ratio,digits))
    check(f'log_factor_six_required_{case}_{digits}',abs(actual-actual/6)/(1+abs(actual))>100*errors[-1])
    levels.append({'precision':digits,'integral':mp.nstr(actual,digits),'bound':mp.nstr(bound,digits),'finite_differences':records})
  with mp.workdps(80):
   aa=mp.mpf(levels[0]['integral']); bb=mp.mpf(levels[1]['integral']); err=abs(aa-bb)/(1+abs(bb))
   check('log_precision_'+str(case),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'noncommuting_log_derivative','s':ss,'mass_squared':ms,'vector':vs,'levels':levels,'charged_operator':False})

def horizon_controls():
 rho=S.symbols('rho',positive=True)
 for case,ns in enumerate(PROTOCOL['horizon_exponents']):
  nu=S.Rational(ns); cc=S.Rational(PROTOCOL['horizon_c']); p=rho**nu*(1+cc*rho**2); density=rho*S.diff(p,rho,2)**2
  primitive=S.integrate(density,rho); exact=S.lambdify(rho,primitive,'mpmath'); levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    nn=mp.mpf(ns); c=mp.mpf(PROTOCOL['horizon_c']); pfun=lambda r:r**nn*(1+c*r*r); vals=[]
    for j,ds in enumerate(PROTOCOL['horizon_cutoffs']):
     cut=mp.mpf(ds); value=mp.quad(lambda r:r*mp.diff(pfun,r,2)**2,[0,cut]); target=exact(cut)
     err=abs(value-target)/(1+abs(target)); check(f'noninteger_horizon_density_{case}_{digits}_{j}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
     vals.append(mp.nstr(value,digits))
    check(f'noninteger_horizon_integrability_{case}_{digits}',mp.mpf(vals[1])/mp.mpf(vals[0])<mp.mpf('.13'))
    levels.append({'precision':digits,'integrals':vals})
  with mp.workdps(80):
   errors=[abs(mp.mpf(a)-mp.mpf(b))/(1+abs(mp.mpf(b))) for a,b in zip(levels[0]['integrals'],levels[1]['integrals'])]
   check('horizon_precision_'+str(case),max(errors)<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(max(errors),80))
  FIXTURES.append({'fixture':'noninteger_horizon_density','nu':ns,'levels':levels,'smooth_periodic_Euclidean_harmonic_claimed':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n'); sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths(); hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}; before=hashes(); error=None
 print(json.dumps({'event':'start','story':'S168 outer continuous causal shape bound','input_count':len(ps)}),flush=True)
 try: provenance(); exact_algebra(); tower_bounds(); measure_brackets(); log_derivatives(); horizon_controls()
 except Exception: error=traceback.format_exc(); print(error,file=sys.stderr,flush=True)
 after=hashes(); good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 outer continuous causal shape bound','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'actual_outer_continuous_compensation_proved':good,'outer_supported_full_causal_core_bound_proved':good,
  'four_anchor_interpolation_controls_negative_weights':good,'Beckner_source_locator_corrected_constant_unchanged':good,
  'global_continuous_shape_coercivity_proved':False,'actual_continuous_Schur_family_verified':False,'actual_numeric_threshold_certified':False,
  'actual_charged_two_point_kernel_evaluated':False,'physical_fast_branch_excluded':False,'finite_counterterm_changed':False,
  'new_boundary_prescription_adopted':False,'physical_source_projection_adopted':False,'order_reduction_adopted':False,
  'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,'charged_Lorentzian_instability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN): run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json'); data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}; data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n'); print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__': raise SystemExit(main())
