"""Global radial graph, logarithmic resolvents, localization and form budgets."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]; STEM='global-continuous-causal-shape-inverse'
PREVIOUS='outer-continuous-causal-shape-bound'
OWN=['checks/global_continuous_causal_shape_inverse.py','notes/'+STEM+'.md','evidence/global-continuous-causal-source-review.json']
ROWS=[]; FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-28','precision_tolerance':'1e-28',
 'bolt_cases':[['10','3.25'],['17.3','4.5']],'bolt_warp':'cos(Y), -pi/2<Y<pi/2',
 'bolt_profile':'cos(Y)^nu (1+sin(Y)^2/5)',
 'bolt_quadrature_coordinate':'u=sin(Y) in [0,1], a dY=du; exact real endpoint',
 'matrix_A':[['2','.7'],['.7','1']],'matrix_B':[['1','0'],['0','1000']],
 'matrix_s':['10','17.3'],'matrix_vector':['.8','-.3'],
 'resolvent_t':['.03','.4','3'],'resolvent_weights':['.7','1.3','.2'],
 'finite_difference_steps':['.002','.001'],'finite_difference_tolerance':'1e-6','refinement_ratio_ceiling':'.3',
 'jensen_s':['3.7','25'],'jensen_A':[['2','-1','0'],['-1','2','-1'],['0','-1','2']],
 'jensen_B':[['1','0','0'],['0','7','0'],['0','0','100']],
 'jensen_angles':['.1','.7','1.3'],
 'actual_charged_kernel_evaluated':False,'parameter_continuity_proved':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.simplify(S.sympify(value).doit());check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==101 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['global_continuous_shape_coercivity_proved'] and review['closed_auxiliary_causal_shape_inverse_proved']
       and review['uniform_continuous_ray_strong_feedback_bound_proved'] and not review['parameter_continuity_proved']
       and not review['physical_boundary_state_admissibility_proved'] and not review['goal_complete'])

def exact_algebra():
 y,s=S.symbols('y s',real=True);a=S.Function('a')(y);v=S.Function('v')(y);u=S.Function('u')(y)
 B=a**-2;A=lambda f:-S.diff(a*S.diff(f,y),y)/a; H=lambda f:A(f)+s*B*f
 K=4*S.diff(a,y)**2-2*a*S.diff(a,y,2)
 zero('radial_multiplier_Laplacian',-A(B)-K*B**2)
 lhs=a*H(v)**2;rhs=a*A(v)**2+2*s*a*B*S.diff(v,y)**2+a*(s*s-s*K)*B**2*v*v
 flux=s*a*S.diff(B,y)*v*v-2*s*a*B*v*S.diff(v,y)
 zero('exact_graph_density_and_flux',lhs-rhs-S.diff(flux,y))
 potential=(s+a*S.diff(a,y,2)/2-S.diff(a,y)**2/4)/a**2
 zero('half_density_potential',S.sqrt(a)*H(u/S.sqrt(a))+S.diff(u,y,2)-potential*u)
 theta=S.Function('theta')(y);chis=[S.cos(theta),S.sin(theta)]
 zero('exact_continuum_IMS',sum(c*H(c*v) for c in chis)-H(v)-S.diff(theta,y)**2*v)
 z=S.symbols('z',positive=True);t=S.symbols('t',nonnegative=True)
 zero('factored_log_resolvent_identity',S.integrate(t*(1-1/(1+z*t)),(t,0,1))-(S.log(1+z)-z+z*z/2)/z**2)
 n,q=S.symbols('n q',integer=True,positive=True);m=S.symbols('m',integer=True,positive=True)
 zero('angular_count_exact',S.summation(2*(q+2*m),(m,1,n))-2*n*(n+q)-2*n)
 zero('angular_count_margin',3*n*(n+q)-(2*n*(n+q)+2*n)-n*(n+q-2))
 zero('counting_resolvent_integral',S.integrate(6*t/(t+z)**3,(t,0,S.oo))-3/z)
 zero('alpha_moment',S.integrate(t*t*(1-t)**2,(t,0,1))-S.Rational(1,30))
 zero('log_derivative_budget',S.Rational(3,2)*S.Rational(1,30)*S.Rational(1,2)-S.Rational(1,40))
 zero('product_third_derivative_budget',6*(1+2)**2/S.Integer(40)-S.Rational(27,20))
 x=S.symbols('x',nonnegative=True);nn=x+20
 for name,p in {'relative_nodes':nn*nn/4-4*nn-4,
                'remainder_width':48*nn**3-(6*nn+3)*(2*nn+1)**2}.items():
  co=S.Poly(S.expand(p),x).all_coeffs();check('uniform_'+name,all(c>=0 for c in co) and co[-1]>0,co)
 zero('anchor_fourth_budget',S.Rational(3,2)**3/500-S.Rational(27,4000))
 zero('intermediate_graph_budget',S.Rational(9,4)/S.Rational(3,4)**3-S.Rational(16,3))
 zero('interpolation_budget',S.Rational(16,3)*S.Rational(27,20)*2-S.Rational(72,5))
 check('interpolation_cost_below_one_thousandth',S.Rational(72,5*3*20**3)<S.Rational(1,1000))
 zero('continuous_fourth_budget',S.Rational(27,4000)+S.Rational(1,1000)-S.Rational(31,4000))
 zero('lower_order_absorption',S.Rational(31,4000)+S.Rational(1,1000)-S.Rational(7,800))
 zero('Jensen_scalar_margin',S.Rational(33,20)-S.Rational(1,20)-S.Rational(8,5))
 zero('scalar_compensation_coefficient',S.Rational(8,5)/30-S.Rational(4,75))
 check('KLMN_strict_relative_bound',S.Rational(525,3200)*S.Rational(22,7)<S.Rational(3,5))
 check('actual_fourth_reserve',S.Rational(4,75)-S.Rational(7,800)*S.Rational(22,7)>S.Rational(1,200))
 check('classical_residual_reserve',S.Rational(1,32)-S.Rational(1,10**24)>S.Rational(1,33))
 check('physical_gradient_reserve',S.Rational(1,33)>S.Rational(1001,1000)/40)
 ab=S.cos(y);kb=S.trigsimp(4*S.diff(ab,y)**2-2*ab*S.diff(ab,y,2));kv=S.trigsimp(S.diff(ab,y)**2/4-ab*S.diff(ab,y,2)/2)
 zero('smooth_bolt_KB',kb-(2+2*S.sin(y)**2))
 zero('smooth_bolt_KV',kv-(S.Rational(1,2)-S.sin(y)**2/4))
 check('smooth_bolt_unbounded_multiplier',S.limit(1/S.cos(y)**2,y,S.pi/2,dir='-')==S.oo)
 uu,pp,cc=S.symbols('u p c',positive=True);dd=1-uu**2;ff=1+cc*uu**2;vv=dd**pp*ff
 directA=-S.diff(dd*S.diff(vv,uu),uu)
 explicitA=dd**(pp-1)*((2*pp*dd-4*pp**2*uu**2)*ff+(4+8*pp)*cc*uu**2*dd-2*cc*dd**2)
 positive_d=S.symbols('positive_d',positive=True)
 reduced=S.simplify((directA-explicitA).subs(dd,positive_d)/positive_d**(pp-1))
 zero('real_bolt_coordinate_operator',S.factor(reduced.subs(positive_d,dd)))

def parse(M):return mp.matrix([[mp.mpf(x) for x in row] for row in M])
def serial(M,digits):return [[mp.nstr(M[i,j],digits) for j in range(M.cols)] for i in range(M.rows)]
def eigmin(M):return min(mp.eigsy((M+M.T)/2,eigvals_only=True))
def normop(M):return max(abs(x) for x in mp.eigsy((M+M.T)/2,eigvals_only=True))
def matrixlog(M):
 ev,U=mp.eigsy(M);return U*mp.diag([mp.log(x) for x in ev])*U.T

def bolt_quadrature():
 for case,(ss,ns) in enumerate(PROTOCOL['bolt_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    s=mp.mpf(ss);nu=mp.mpf(ns);p=nu/2;c=mp.mpf(1)/5
    d=lambda u:1-u*u;f=lambda u:1+c*u*u
    Av=lambda u:d(u)**(p-1)*((2*p*d(u)-4*p*p*u*u)*f(u)+(4+8*p)*c*u*u*d(u)-2*c*d(u)**2)
    Hv=lambda u:d(u)**(p-1)*((s+2*p*d(u)-4*p*p*u*u)*f(u)+(4+8*p)*c*u*u*d(u)-2*c*d(u)**2)
    vu=lambda u:u*d(u)**(p-1)*(2*c*d(u)-2*p*f(u))
    bnfun=lambda u:d(u)**(2*p-2)*f(u)**2
    grid=[0,mp.mpf('.5'),1]
    direct=2*mp.quad(lambda u:Hv(u)**2,grid)
    rebuilt=2*mp.quad(lambda u:Av(u)**2+2*s*vu(u)**2+(s*s-s*(2+2*u*u))*bnfun(u),grid)
    bn=2*mp.quad(bnfun,grid)
    err=abs(direct-rebuilt)/(1+abs(direct));check(f'bolt_independent_graph_integral_{case}_{digits}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
    ratio=s*(s-4)*bn/direct;check(f'bolt_graph_lower_bound_{case}_{digits}',0<ratio<=1,mp.nstr(ratio,digits))
    wrong=2*mp.quad(lambda u:Av(u)**2+2*s*vu(u)**2+s*s*bnfun(u),grid)
    check(f'bolt_curvature_term_required_{case}_{digits}',abs(wrong-direct)/direct>mp.mpf('.01'))
    levels.append({'precision':digits,'direct':mp.nstr(direct,digits),'rebuilt':mp.nstr(rebuilt,digits),'lower_bound_ratio':mp.nstr(ratio,digits)})
  with mp.workdps(80):
   a0=mp.mpf(levels[0]['direct']);a1=mp.mpf(levels[1]['direct']);err=abs(a0-a1)/(1+abs(a1))
   check('bolt_precision_'+str(case),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'smooth_bolt_graph','s':ss,'nu':ns,'levels':levels,'charged_geometry':False})

def factored_derivatives():
 for case,ss in enumerate(PROTOCOL['matrix_s']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    A=parse(PROTOCOL['matrix_A']);B=parse(PROTOCOL['matrix_B']);s=mp.mpf(ss);v=mp.matrix(list(map(mp.mpf,PROTOCOL['matrix_vector'])))
    ts=list(map(mp.mpf,PROTOCOL['resolvent_t']));weights=list(map(mp.mpf,PROTOCOL['resolvent_weights']));I=mp.eye(2);M=sum(weights)/4
    check(f'noncommuting_factor_fixture_{case}_{digits}',mp.norm(A*B-B*A)>1)
    def gmat(x):return sum((w*(I-(I+t*(A+x*B))**-1) for t,w in zip(ts,weights)),mp.zeros(2))
    H=A+s*B;h=H*v;b=B*v;derivs=[]
    for j in [1,2,3]:
     positive=sum((w*mp.factorial(j)*t**j*((I+t*H)**-1*B)**j*(I+t*H)**-1 for t,w in zip(ts,weights)),mp.zeros(2))
     ratio=normop(positive)/(mp.factorial(j)*M/s**j)
     check(f'factor_derivative_order_{case}_{digits}_{j}',eigmin(positive)>=0 and ratio<=1,mp.nstr(ratio,digits))
     derivs.append((-1)**(j+1)*positive)
    direct=(h.T*derivs[2]*h+6*b.T*derivs[1]*h+6*b.T*derivs[0]*b)[0]
    phi=lambda x:(((A+x*B)*v).T*gmat(x)*((A+x*B)*v))[0]
    errors=[];fdrows=[]
    for hs in PROTOCOL['finite_difference_steps']:
     step=mp.mpf(hs);fd=(phi(s+2*step)-2*phi(s+step)+2*phi(s-step)-phi(s-2*step))/(2*step**3)
     err=abs(fd-direct)/(1+abs(direct));errors.append(err)
     check(f'factor_finite_difference_{case}_{digits}_{hs}',err<mp.mpf(PROTOCOL['finite_difference_tolerance']),mp.nstr(err,digits))
     fdrows.append({'step':hs,'value':mp.nstr(fd,digits),'scaled_error':mp.nstr(err,digits)})
    ratio=errors[1]/errors[0];check(f'factor_refinement_{case}_{digits}',ratio<mp.mpf(PROTOCOL['refinement_ratio_ceiling']),mp.nstr(ratio,digits))
    omitted=(h.T*derivs[2]*h)[0];check(f'factor_product_terms_required_{case}_{digits}',abs(omitted-direct)/(1+abs(direct))>100*errors[-1])
    levels.append({'precision':digits,'third_derivative':mp.nstr(direct,digits),'finite_differences':fdrows})
  with mp.workdps(80):
   aa=mp.mpf(levels[0]['third_derivative']);bb=mp.mpf(levels[1]['third_derivative']);err=abs(aa-bb)/(1+abs(bb))
   check('factor_precision_'+str(case),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'noncommuting_factored_derivative','s':ss,'levels':levels,'charged_operator':False})

def jensen_localization():
 for case,ss in enumerate(PROTOCOL['jensen_s']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    A=parse(PROTOCOL['jensen_A']);B=parse(PROTOCOL['jensen_B']);s=mp.mpf(ss);H=A+s*B;I=mp.eye(3)
    angles=list(map(mp.mpf,PROTOCOL['jensen_angles']));chis=[mp.diag([fun(x) for x in angles]) for fun in [mp.cos,mp.sin]]
    comp=sum((c*H*c for c in chis),mp.zeros(3));W=comp-H;wn=normop(W);gap=s*min(mp.eigsy(B,eigvals_only=True))
    loc=sum((c*matrixlog(H)*c for c in chis),mp.zeros(3))
    first=eigmin(matrixlog(comp)-loc);second=eigmin(matrixlog(H)+wn/gap*I-loc)
    check(f'jensen_operator_concavity_{case}_{digits}',first>0,mp.nstr(first,digits))
    check(f'jensen_localization_budget_{case}_{digits}',second>0,mp.nstr(second,digits))
    check(f'jensen_noncommuting_partition_{case}_{digits}',mp.norm(H*chis[0]-chis[0]*H)>mp.mpf('.1'))
    check(f'jensen_wrong_direction_detected_{case}_{digits}',eigmin(loc-matrixlog(comp))<0)
    levels.append({'precision':digits,'log_H':serial(matrixlog(H),digits),'Jensen_gap_min':mp.nstr(first,digits),'budget_gap_min':mp.nstr(second,digits)})
  with mp.workdps(80):
   aa=parse(levels[0]['log_H']);bb=parse(levels[1]['log_H']);err=mp.norm(aa-bb)/(1+mp.norm(bb))
   check('jensen_precision_'+str(case),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'noncommuting_Jensen_localization','s':ss,'levels':levels,'discrete_W_is_continuum_multiplier':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 global continuous causal shape inverse','input_count':len(ps)}),flush=True)
 try:provenance();exact_algebra();bolt_quadrature();factored_derivatives();jensen_localization()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 global continuous causal shape inverse','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'global_continuous_shape_coercivity_proved':good,'closed_auxiliary_causal_shape_inverse_proved':good,
  'uniform_continuous_ray_strong_feedback_bound_proved':good,
  'parameter_continuity_proved':False,'actual_continuous_Schur_family_verified':False,'actual_numeric_threshold_certified':False,
  'actual_charged_two_point_kernel_evaluated':False,'physical_fast_branch_excluded':False,'finite_counterterm_changed':False,
  'new_boundary_prescription_adopted':False,'physical_source_projection_adopted':False,'order_reduction_adopted':False,
  'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,'charged_Lorentzian_instability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]};data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
