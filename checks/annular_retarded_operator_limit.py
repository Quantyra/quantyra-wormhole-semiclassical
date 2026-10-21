"""S168 annular logarithmic form-norm transfer; ordinary proof plus fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='annular-retarded-operator-limit'
OWN=['checks/annular_retarded_operator_limit.py','notes/'+STEM+'.md','evidence/annular-operator-source-review.json']
PREVIOUS='complex-retarded-outer-limit'
ROWS=[];NUMERICS=[]
PROTOCOL={
 'scope':'Annular X=H2 x H2(log)^(1/2) operator-norm limit, not a global exterior inverse',
 'mpmath_precision':60,'identity_tolerance':'1e-38',
 'matrix_log_perturbations':['1','0.5','0.25'],'matrix_refinement_factor':'0.6',
 'step_a':'1.3','step_height':'0.7','support_b':'0.5','step_locations':['1.5','2.5','3.5'],
 'step_resolvent_shifts':['0','0.5','4'],'observation_source_pairs':[['-0.3','0.2'],['0.3','-0.1']],
 'measure_dimensions':[3,7,11],'complex_points':[['0.4','0.9'],['-2','0.4'],['4','-1.2']],
 'interpolation_frequencies':[1,2,17,257],
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
 check('predecessor_passed',p['status']=='passed' and len(p['checks'])==93 and all(r['passed'] for r in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 required=['uniform-logarithmic-frequency-comparison','outer-coupled-scalar-separation',
  'retarded-outer-frequency-transfer','canonical-frequency-core-bound','local-coupled-frequency-polynomial']
 check('required_uniform_inputs_pinned',all('notes/'+f+'.md' in p['input_sha256_after'] for f in required))
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('source_scope',len(review['primary_sources'])==2 and not review['global_exterior_operator_transfer_proved']
  and not review['physical_metric_boundary_selected'] and not review['quantitative_physical_collar_radius_proved'])

def exact_checks():
 t,A,k,c,L,M,R,b,x,y,q=S.symbols('t A k c L M R b x y q',positive=True)
 u0=t*t*A/(1+t*t*A);u1=1/(1+t*t*A)
 zero('K_decomposition',u0+u1-1)
 zero('K_minimizer_value',u0*u0+t*t*A*u1*u1-t*t*A/(1+t*t*A))
 z=S.symbols('z',real=True)
 zero('K_minimizer_stationarity',S.diff(z*z+t*t*A*(1-z)**2,z).subs(z,u0))
 zero('critical_log_interpolation',2*S.integrate(t*A/(1+t*t*A),(t,0,1))-S.log(1+A))
 zero('tail_resolvent_time_integral',S.integrate((c+t)**(-S.Rational(7,4)),(t,0,S.oo))-4/(3*c**S.Rational(3,4)))
 zero('tail_coefficient',M*S.sqrt(L)/2*4/(3*c**S.Rational(3,4))-2*M*S.sqrt(L)/(3*c**S.Rational(3,4)))
 zero('free_kernel_ODE',-S.diff(S.exp(-k*(x-y))/(2*k),x,2)+k*k*S.exp(-k*(x-y))/(2*k))
 zero('free_kernel_derivative_jump',-S.Rational(1,2)-S.Rational(1,2)+1)
 hs_one=S.integrate(S.exp(-2*k*(x-y))/(4*k*k),(x,R,S.oo))
 zero('free_kernel_one_tail_integral',hs_one-S.exp(-2*k*(R-y))/(8*k**3))
 zero('step_reflection_coefficient',(k-q)/(k+q)+(q-k)/(k+q))
 reflection=(k-q)/(k+q)
 zero('step_value_match',1+reflection-2*k/(k+q))
 zero('step_derivative_match',k*(reflection-1)+q*(1+reflection))
 zero('V_bound',S.Rational(42,2)+S.Rational(12,10)**2/4-S.Rational(534,25))
 check('V_bound_below_22',S.Rational(534,25)<22)
 check('V_prime_bound_below_2000',S.Rational(3270,2)+S.Rational(12,10)*42/2<2000)
 d=S.symbols('d',positive=True)
 zero('logarithmic_cutoff_tail',S.exp(-S.sqrt(c)*(2*S.log(1/d)/S.sqrt(c)))-d*d)
 zero('toy_arbitrary_epsilon_loss',d**S.Rational(2,3)+d*d/(d**S.Rational(2,3))**2+d**4/(d**S.Rational(2,3))**3
   -(2*d**S.Rational(2,3)+d*d))
 check('fixed_loss_does_not_vanish',S.limit(S.Rational(1,1000)+d*d+d**4,d,0)==S.Rational(1,1000))
 n=S.symbols('n',positive=True)
 check('H2_alone_misses_log_order',S.limit(S.log(1+n*n),n,S.oo)==S.oo)
 check('log_form_weight_controls_multiplier',S.limit(S.log(1+n*n)/(1+S.log(2+n*n)),n,S.oo)==1)
 for dim in [8,16,32]:
  P=S.zeros(dim);P[dim-1,dim-1]=1
  check('fixed_compression_misses_unit_norm_'+str(dim),P[:4,:4]==S.zeros(4) and S.trace(P*P)==1)

def maxabs(A):return max(abs(x) for x in A)
def opnorm(A):return max(mp.svd(A,compute_uv=False))
def logmat(A):
 values,Q=mp.eigsy(A);return Q*mp.diag([mp.log(v) for v in values])*Q.T
def matrix_log_fixtures():
 mp.mp.dps=PROTOCOL['mpmath_precision'];tol=mp.mpf(PROTOCOL['identity_tolerance'])
 B=mp.diag([2,5,9]);W=mp.matrix([['.3','.4','0'],['.4','-.2','.1'],['0','.1','.5']]);I=mp.eye(3)
 logB=logmat(B);previous=None;previous_form=None
 Ninv=mp.diag([1/(v*mp.sqrt(1+mp.log(v+1))) for v in [2,5,9]])
 check('matrix_fixture_noncommuting',maxabs(B*W-W*B)>mp.mpf('.1'))
 for delta_s in PROTOCOL['matrix_log_perturbations']:
  delta=mp.mpf(delta_s);A=B+delta*W;logA=logmat(A)
  check('matrix_common_gap_'+delta_s,min(mp.eigsy(A,eigvals_only=True))>mp.mpf('1.5'))
  def integrand(t):return (B+t*I)**-1-(A+t*I)**-1
  integral=mp.matrix([[mp.quad(lambda t:integrand(t)[i,j],[0,1,4,mp.inf]) for j in range(3)] for i in range(3)])
  err=maxabs(logA-logB-integral)
  check('independent_matrix_log_integral_'+delta_s,err<tol,mp.nstr(err,10))
  shift=mp.mpf('1.3');RA=(A+shift*I)**-1;RB=(B+shift*I)**-1
  check('noncommutative_resolvent_order_'+delta_s,maxabs(RB-RA-RA*(delta*W)*RB)<tol)
  error=opnorm(logA-logB);bound=delta*opnorm(W)/mp.mpf('1.5')
  check('matrix_log_norm_bound_'+delta_s,error<=bound,mp.nstr(error/bound,12))
  formerror=opnorm(Ninv*(A*logA*A-B*logB*B)*Ninv)
  if previous is not None:
   factor=mp.mpf(PROTOCOL['matrix_refinement_factor'])
   check('matrix_log_refinement_'+delta_s,error<factor*previous)
   check('critical_form_refinement_'+delta_s,formerror<factor*previous_form)
  NUMERICS.append({'fixture':'matrix_log','delta':delta_s,'integral_error':mp.nstr(err,16),
   'log_norm_error':mp.nstr(error,20),'critical_form_norm_error':mp.nstr(formerror,20),'bound_ratio':mp.nstr(error/bound,18)})
  previous=error;previous_form=formerror

def step_green_matched(x,y,R,a,height,t):
 k=mp.sqrt(a+t);q=mp.sqrt(a+height+t);ey=mp.exp(k*y);iy=1/ey;eR=mp.exp(k*R);iR=1/eR
 coeff=mp.matrix([[ey,-ey,-iy,0],[-k*ey,k*ey,-k*iy,0],[0,eR,iR,-1],[0,k*eR,-k*iR,q]])
 sol=mp.lu_solve(coeff,mp.matrix([0,-1,0,0]))
 return sol[0]*mp.exp(k*x) if x<y else sol[1]*mp.exp(k*x)+sol[2]*mp.exp(-k*x)

def step_fixtures():
 tol=mp.mpf(PROTOCOL['identity_tolerance']);a=mp.mpf(PROTOCOL['step_a']);height=mp.mpf(PROTOCOL['step_height'])
 b=mp.mpf(PROTOCOL['support_b']);length=2*b;previous=None
 for R_s in PROTOCOL['step_locations']:
  R=mp.mpf(R_s);maxerr=0
  for ts in PROTOCOL['step_resolvent_shifts']:
   t=mp.mpf(ts);k=mp.sqrt(a+t);q=mp.sqrt(a+height+t)
   for ix,(xs,ys) in enumerate(PROTOCOL['observation_source_pairs']):
    x=mp.mpf(xs);y=mp.mpf(ys)
    explicit=(mp.exp(-k*abs(x-y))+(k-q)/(k+q)*mp.exp(-k*(2*R-x-y)))/(2*k)
    error=abs(step_green_matched(x,y,R,a,height,t)-explicit);maxerr=max(maxerr,error)
    check('independent_step_ODE_'+R_s+'_'+ts+'_'+str(ix),error<tol,mp.nstr(error,10))
  def refl(k):return (mp.sqrt(k*k+height)-k)/(mp.sqrt(k*k+height)+k)
  lower=mp.sqrt(a)
  trace=mp.quad(lambda k:refl(k)*(mp.exp(-2*k*(R-b))-mp.exp(-2*k*(R+b)))/(2*k),[lower,lower+1,mp.inf])
  box=mp.quad(lambda k:refl(k)*(mp.exp(-k*(R-b))-mp.exp(-k*(R+b)))**2/(k*k*length),[lower,lower+1,mp.inf])
  bound=2*height*mp.sqrt(length)*mp.exp(-mp.sqrt(a)*(R-b))/(3*a**mp.mpf('.75'))
  check('localized_log_positive_trace_bound_'+R_s,0<box<=trace<=bound,[mp.nstr(v,16) for v in [box,trace,bound]])
  if previous is not None:check('distant_step_log_tail_decreases_'+R_s,trace<previous/5)
  NUMERICS.append({'fixture':'distant_step','R':R_s,'matched_ODE_error':mp.nstr(maxerr,16),
   'unit_box_log_form':mp.nstr(box,20),'localized_log_trace':mp.nstr(trace,20),'proved_tail_bound':mp.nstr(bound,20)})
  previous=trace

def interpolation_and_measure_fixtures():
 tol=mp.mpf(PROTOCOL['identity_tolerance'])
 As=[mp.mpf(v) for v in PROTOCOL['interpolation_frequencies']];coeff=[mp.mpf(1),mp.mpf('-.3'),mp.mpf('.02'),mp.mpf('.001')]
 base=sum(A*A*abs(v)**2 for A,v in zip(As,coeff))
 quad=base+2*mp.quad(lambda t:sum(A*A*(t*t*A/(1+t*t*A))*abs(v)**2 for A,v in zip(As,coeff))/t,[0,mp.mpf('.01'),mp.mpf('.1'),1])
 direct=sum(A*A*(1+mp.log(1+A))*abs(v)**2 for A,v in zip(As,coeff))
 check('independent_log_interpolation_integral',abs(quad-direct)<tol,mp.nstr(abs(quad-direct),12))
 NUMERICS.append({'fixture':'critical_interpolation','quadrature_error':mp.nstr(abs(quad-direct),16),'norm_squared':mp.nstr(direct,20)})
 ts=[mp.mpf('.2'),mp.mpf(2),mp.mpf(13)]
 for dim in PROTOCOL['measure_dimensions']:
  vectors=[mp.matrix([1 for i in range(dim)]),mp.matrix([i+1 for i in range(dim)]),mp.matrix([(-1)**i for i in range(dim)])]
  vectors=[v/mp.norm(v) for v in vectors]
  matrices=[scale*v*v.T for scale,v in zip([1,mp.mpf('.7'),mp.mpf('.3')],vectors)]
  G=sum(matrices,mp.zeros(dim));gnorm=opnorm(G);trace=sum(G[i,i] for i in range(dim))
  check('operator_bound_uses_no_trace_factor_'+str(dim),gnorm<trace and maxabs(matrices[0]*matrices[1]-matrices[1]*matrices[0])>mp.mpf('.001'))
  for ix,(re,im) in enumerate(PROTOCOL['complex_points']):
   z=mp.mpc(re,im);rs=[(t+4)/(t+z) for t in ts];maximum=max(abs(r) for r in rs)
   T=sum((r*M for r,M in zip(rs,matrices)),mp.zeros(dim));norm=opnorm(T)
   check('dimension_free_positive_measure_'+str(dim)+'_'+str(ix),norm<=maximum*gnorm+tol,mp.nstr(norm/(maximum*gnorm),12))
   u=vectors[0]+1j*vectors[2];w=vectors[1]-1j*vectors[0]
   value=abs((u.H*T*w)[0]);cs=maximum*mp.sqrt(mp.re((u.H*G*u)[0])*mp.re((w.H*G*w)[0]))
   check('positive_measure_cross_Cauchy_Schwarz_'+str(dim)+'_'+str(ix),value<=cs+tol)
   NUMERICS.append({'fixture':'positive_operator_measure','dimension':dim,'z':[re,im],
    'operator_bound_ratio':mp.nstr(norm/(maximum*gnorm),18),'total_operator_norm':mp.nstr(gnorm,18),'total_trace':mp.nstr(trace,18)})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 annular retarded operator limit','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();matrix_log_fixtures();step_fixtures();interpolation_and_measure_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 annular retarded operator limit','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'localized_logarithm_tail_bound_derived':good,'critical_logarithmic_form_norm_retained':good,
  'actual_annular_operator_norm_limit_derived':good,'complex_operator_derivative_limit_derived':good,
  'all_local_contacts_retained':good,'global_exterior_operator_transfer_proved':False,
  'physical_metric_boundary_selected':False,'actual_metric_causal_inverse_proved':False,
  'actual_metric_mode_constructed':False,'quantitative_physical_collar_radius_proved':False,
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
