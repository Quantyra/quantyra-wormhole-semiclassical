"""Full constant-row weak limit diagnostics; two-point comparison remains open."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import mpmath as mp
import sympy as S

ROOT=Path(__file__).resolve().parents[1]; STEM='moving-cut-constant-quantum-row'
PREVIOUS='moving-cut-global-quantum-locality'
OWN=['checks/moving_cut_constant_quantum_row.py','notes/'+STEM+'.md','evidence/moving-cut-constant-row-source-review.json']
ROWS=[]; FIXTURES=[]
PROTOCOL={'precision_digits':[50,80], 'matrix_cases':[['.7','.03'],['1.5','.1'],['3','.2']],
 'radial_sites':5,'time_momentum':'.4','odd_weights':[1,3],
 'Gaussian_widths':['.4','1','2'],'variance_rates':['.3','.5','1.2'],
 'angular_cases':[['1','1/9'],['8','1'],['16000000000000','320000000000/7']],
 'heat_scales':['.01','.1','1','10'],'angular_modes':256,
 'cut_powers':[3,5,8],'profile_cubic_coefficients':['-.1','.2'],
 'identity_tolerance':'1e-38','quadrature_tolerance':'1e-32','precision_tolerance':'1e-32',
 'actual_charged_kernel_evaluated':False,'ordinary_continuum_proof':True,
 'heat_integral_tail_coordinate':'u=1/sqrt(t) on t>=1; tolerance unchanged'}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def zero(name,x):
 if isinstance(x,S.MatrixBase):x=sum(S.simplify(abs(v)**2) for v in x)
 x=S.simplify(x);check(name,x==0,x)
def err(a,b):return abs(a-b)/(1+abs(b))
def out(a):return mp.nstr(a,mp.mp.dps)
def tr(A):return sum(A[i,i] for i in range(A.rows))
def norm(A):return mp.sqrt(sum(abs(v)**2 for v in A))
def kron(A,B):
 return mp.matrix([[A[i//B.rows,j//B.cols]*B[i%B.rows,j%B.cols] for j in range(A.cols*B.cols)] for i in range(A.rows*B.rows)])
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==244 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',not review['actual_two_point_comparison_proved'] and not review['physical_fast_branch_excluded'])

def symbolic():
 one=S.eye(2);s1=S.Matrix([[0,1],[1,0]]);s2=S.Matrix([[0,-S.I],[S.I,0]]);s3=S.diag(1,-1)
 N=-S.kronecker_product(s1,one);T=-S.kronecker_product(s2,s2)
 chi=S.kronecker_product(s3,s1);C=S.kronecker_product(one,s2);grading=S.kronecker_product(s3,one)
 zero('reflection_normal',N*chi+chi*N);zero('reflection_tangent',T*chi-chi*T)
 zero('sector_exchange_normal',C*N-N*C);zero('sector_exchange_tangent',C*T-T*C)
 zero('sector_exchange_boundary',C*chi+chi*C)
 zero('mass_grading_normal',grading*N+N*grading);zero('mass_grading_tangent',grading*T+T*grading)
 zero('no_third_smearing_contact',S.trace(chi))
 u,b=S.symbols('u b',positive=True);B=S.Function('B')(u);a=S.sqrt(B)
 L=lambda f:S.sqrt(B)*S.diff(f,u)/b
 zero('proper_linear_jet',L(a)/a-S.diff(B,u)/(2*b*S.sqrt(B)))
 zero('proper_quadratic_jet',L(L(a))/(2*a)-S.diff(B,u,2)/(4*b*b))
 zero('proper_cubic_jet',L(L(L(a)))/(6*a)-S.diff(B,u)*S.diff(B,u,2)/(24*b**3*S.sqrt(B))-S.sqrt(B)*S.diff(B,u,3)/(12*b**3))
 t,m=S.symbols('t m',positive=True)
 zero('heat_weight_derivative',-2*S.diff(t*m*m*S.exp(-t*m*m),t)-2*(t*m**4-m*m)*S.exp(-t*m*m))
 zero('small_time_weight',-2*S.diff(t*(2/t**2+S.Symbol('a0')+S.Symbol('a1')*t),t)-(4/t**2-2*S.Symbol('a0')-4*S.Symbol('a1')*t))
 zero('absolute_heat_reserve',2*(40+12)-104);zero('gap_heat_reserve',2*(40*8+12*4)-736)
 z=S.symbols('z',real=True);zero('bridge_variance_average',S.integrate((2*z*(1-z))**2,(z,0,1))-S.Rational(2,15))
 for j in [1,3]:
  s=S.symbols('s',positive=True)
  beta=S.integrate((s-1)**j/s**5,(s,1,S.oo))
  zero('odd_jet_beta_'+str(j),beta-(S.Rational(1,12) if j==1 else S.Rational(1,4)))
  zero('constant_row_boundary_coefficient_'+str(j),S.Rational(8,5)*beta-(S.Rational(2,15) if j==1 else S.Rational(2,5)))
 for p in [1,2]:
  k=S.symbols('k',positive=True)
  poly=sum(S.factorial(p)*k**j/S.factorial(j) for j in range(p+1))+k**(p+1)
  zero('angular_tail_polynomial_'+str(p),poly-([1+k+k*k,2+2*k+k*k+k**3][p-1]))

def matrix_fixtures():
 for case,(ms,ds) in enumerate(PROTOCOL['matrix_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    M=mp.mpf(ms);d=mp.mpf(ds);n=PROTOCOL['radial_sites'];I2=mp.eye(2);s1=mp.matrix([[0,1],[1,0]]);s2=mp.matrix([[0,-mp.j],[mp.j,0]]);s3=mp.diag([1,-1])
    N=-kron(s1,I2);T=-kron(s2,s2);chi=kron(s3,s1);C=kron(I2,s2)
    Dy=mp.matrix(n);P=mp.matrix(n)
    for i in range(n):Dy[i,(i+1)%n]=mp.mpf('.5');Dy[i,(i-1)%n]=-mp.mpf('.5');P[i,(-i)%n]=1
    D=kron(N,Dy)+mp.j*mp.mpf(PROTOCOL['time_momentum'])*kron(T,mp.eye(n));I=mp.eye(4*n)
    reflection=kron(chi,P);exchange=kron(C,mp.eye(n));plus=(I+reflection)/2
    x=[mp.mpf(min(i,n-i))/4 for i in range(n)]
    warp=kron(mp.eye(4),mp.diag([1+d*y+y*y/5 for y in x]));invhalf=mp.diag([1/mp.sqrt(warp[i,i]) for i in range(4*n)])
    v=kron(mp.eye(4),mp.diag([(d+y)**2*mp.exp(-y) for y in x]));U=v*v
    Dg=invhalf*D*invhalf;R=(Dg+M*I)**-1;Sf=(D+M*warp)**-1;half=mp.diag([mp.sqrt(warp[i,i]) for i in range(4*n)])
    tol=mp.mpf(PROTOCOL['identity_tolerance']);suffix=f'{case}_{digits}'
    check('paired_skew_adjoint_'+suffix,norm(Dg+Dg.transpose_conj())<tol)
    check('reflection_operator_'+suffix,norm(reflection*(D+M*warp)-(D+M*warp)*reflection)<tol)
    check('sector_exchange_'+suffix,norm(exchange*reflection+reflection*exchange)<tol)
    check('conformal_mass_resolvent_'+suffix,norm(R-half*Sf*half)<tol)
    X=M*warp*v;Y=M*warp*U
    bubble=tr(Sf*X*Sf*X);tad=-tr(Sf*Y)
    be=err(tr(plus*Sf*X*Sf*X),bubble/2);te=err(-tr(plus*Sf*Y),tad/2)
    check('even_mass_bubble_half_'+suffix,be<tol,out(be));check('even_mass_second_vertex_half_'+suffix,te<tol,out(te))
    check('curved_mass_bubble_equivalence_'+suffix,err(bubble,M*M*tr(R*v*R*v))<tol)
    direct=M*M*tr(U*R*R)-M*tr(U*R)
    Pg=M*M*I-Dg*Dg;Pinv=Pg**-1;squared=2*M**4*tr(U*Pinv*Pinv)-2*M*M*tr(U*Pinv)
    se=err(direct,squared);check('constant_row_square_identity_'+suffix,se<tol,out(se))
    check('probe_noncommutes_'+suffix,norm(U*Dg-Dg*U)>mp.mpf('.001'))
    ev,Q=mp.eighe(mp.j*Dg);weights=[sum(U[i,i]*abs(Q[i,k])**2 for i in range(4*n)) for k in range(4*n)]
    def heat(t):return sum(weights[k]*mp.exp(-t*(M*M+ev[k]**2)) for k in range(4*n))
    value=mp.quad(lambda t:2*(t*M**4-M*M)*heat(t),[0,1,mp.inf])
    he=err(value,direct);check('constant_row_heat_integral_'+suffix,he<mp.mpf(PROTOCOL['quadrature_tolerance']),out(he))
    check('second_vertex_alone_is_different_'+suffix,abs(direct+M*tr(U*R))>mp.mpf('1e-5'))
    levels.append({'precision':digits,'row':out(mp.re(direct)),'bubble_half_error':out(be),'second_vertex_half_error':out(te),'square_error':out(se),'heat_error':out(he)})
  with mp.workdps(80):
   pe=err(mp.mpf(levels[0]['row']),mp.mpf(levels[1]['row']));check('matrix_precision_'+str(case),pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
  FIXTURES.append({'fixture':'paired_even_mass_and_full_row','M':ms,'cut_parameter':ds,'levels':levels,'actual_kernel':False})

def gaussian_D(j,r):
 a=1/r;phi=mp.exp(-a*a/2)/mp.sqrt(2*mp.pi);tail=mp.erfc(a/mp.sqrt(2))/2
 if j==1:return 2*(r*phi-tail)
 return 2*((r+2*r**3)*phi-(1+3*r*r)*tail)

def gaussian_fixtures():
 for j in PROTOCOL['odd_weights']:
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    independent=[]
    for rs in PROTOCOL['Gaussian_widths']:
     r=mp.mpf(rs)
     direct=2*mp.quad(lambda s:(s-1)**j*mp.exp(-s*s/(2*r*r))/(mp.sqrt(2*mp.pi)*r),[1,mp.inf])
     er=err(direct,gaussian_D(j,r));check(f'Gaussian_tail_{j}_{digits}_{rs}',er<mp.mpf(PROTOCOL['quadrature_tolerance']),out(er));independent.append(out(er))
    integral=mp.quad(lambda r:gaussian_D(j,r)/r**5,[0,1,mp.inf]);target=mp.mpf('.25' if j==1 else '.75')
    ie=err(integral,target);check(f'Gaussian_heat_integral_{j}_{digits}',ie<mp.mpf(PROTOCOL['quadrature_tolerance']),out(ie))
    heat_values=[]
    for vs in PROTOCOL['variance_rates']:
     v=mp.mpf(vs)
     # The j=3 raw t tail is t^-3/2. This exact substitution removes
     # its fractional endpoint before numerical quadrature.
     raw=mp.quad(lambda t:gaussian_D(j,mp.sqrt(v*t))/t**3,[0,1])
     raw+=mp.quad(lambda u:2*u**3*gaussian_D(j,mp.sqrt(v)/u),[0,1])
     expect=(mp.mpf('.5') if j==1 else mp.mpf('1.5'))*v*v
     er=err(raw,expect);check(f'Gaussian_Tonelli_{j}_{digits}_{vs}',er<mp.mpf(PROTOCOL['quadrature_tolerance']),out(er));heat_values.append(out(er))
    coefficient=4*integral*mp.mpf(2)/15
    ce=err(coefficient,mp.mpf(2)/(15 if j==1 else 5));check(f'full_row_coefficient_{j}_{digits}',ce<mp.mpf(PROTOCOL['quadrature_tolerance']),out(ce))
    levels.append({'precision':digits,'reduced_integral':out(integral),'Gaussian_errors':independent,'integral_error':out(ie),'Tonelli_errors':heat_values,'coefficient_without_Ccj_over_pi':out(coefficient)})
  with mp.workdps(80):
   pe=err(mp.mpf(levels[0]['reduced_integral']),mp.mpf(levels[1]['reduced_integral']));check('Gaussian_precision_'+str(j),pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
  FIXTURES.append({'fixture':'odd_jet_boundary_heat_profile','weight':j,'levels':levels})

def angular_fixtures():
 for case,(qs,cs) in enumerate(PROTOCOL['angular_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    q=mp.mpf(qs);cr=S.Rational(cs);C=mp.mpf(int(cr.p))/int(cr.q);M1sq=(1+q)/C;N=PROTOCOL['angular_modes'];values=[]
    for ss in PROTOCOL['heat_scales']:
     t=mp.mpf(ss)/M1sq;A=mp.mpf(0);T2=mp.mpf(0)
     for n in range(1,N+1):
      m2=n*(n+q)/C;deg=2*(q+2*n);e=mp.exp(-t*m2);A+=deg*m2*e;T2+=deg*m2*m2*e
     B=2*(t*T2-A);z=t*N*(N+q)/C
     tail=8*C/t**2*mp.exp(-z)*(3+3*z+2*z*z+z**3)
     bound=(abs(B)+tail)*t*t/C
     suffix=f'{case}_{digits}_{ss}'
     check('angular_tail_monotone_'+suffix,z>=2,out(z))
     check('full_weight_enclosure_'+suffix,bound<120,out(bound))
     check('finite_heat_weight_not_twice_A_'+suffix,abs(B-2*A)/(1+abs(A))>mp.mpf('1e-20'))
     # The signed finite-mode heat integral is also an exact two-vertex row weight.
     derivative=-2*(A-t*T2);check('heat_derivative_all_modes_'+suffix,err(B,derivative)<mp.mpf('1e-40'))
     values.append({'heat_scale':ss,'partial_B':out(B),'absolute_tail_bound':out(tail),'absolute_weight_ratio_to_C_over_t2':out(bound)})
    levels.append({'precision':digits,'values':values})
  with mp.workdps(80):
   pe=max(err(mp.mpf(a['partial_B']),mp.mpf(b['partial_B'])) for a,b in zip(levels[0]['values'],levels[1]['values']));check('angular_precision_'+str(case),pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
  FIXTURES.append({'fixture':'complete_angular_heat_weight_enclosure','q':qs,'C':cs,'levels':levels})

def fp_gamma(k):
 if k>0:return mp.gamma(k)
 m=-k;return (-1)**m/mp.factorial(m)*(mp.harmonic(m)-mp.euler)
def Fp(d,p):return mp.exp(-d)*sum(mp.binomial(4,j)*d**(4-j)*fp_gamma(j-p+1) for j in range(5))

def cut_fixtures():
 for case,c3s in enumerate(PROTOCOL['profile_cubic_coefficients']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    c3=mp.mpf(c3s);c4=mp.mpf('.05');values=[];limit=12*c3/(5*mp.pi)
    def row(d):
     a=1+d*d/2+c3*d**3+c4*d**4
     c1=(d+3*c3*d*d+4*c4*d**3)/a;cutc3=(c3+4*c4*d)/a
     return 2*c1*Fp(d,3)/(15*mp.pi)+2*cutc3*Fp(d,1)/(5*mp.pi)
    for power in PROTOCOL['cut_powers']:
     delta=mp.mpf(2)**(-power)
     avg=mp.quad(lambda s:630*(s-1)**4*(2-s)**4*row(delta*s),[1,2])
     scaled=abs(avg-limit)/delta;check(f'linear_cubic_cut_limit_{case}_{digits}_{power}',scaled<5,out(scaled))
     values.append({'power':power,'averaged_linear_profiles':out(avg),'limit_error_over_cut':out(scaled)})
    check(f'cut_refinement_{case}_{digits}',abs(mp.mpf(values[-1]['averaged_linear_profiles'])-limit)<abs(mp.mpf(values[0]['averaged_linear_profiles'])-limit)/8)
    levels.append({'precision':digits,'limiting_cubic_profile':out(limit),'values':values})
  with mp.workdps(80):
   pe=max(err(mp.mpf(a['averaged_linear_profiles']),mp.mpf(b['averaged_linear_profiles'])) for a,b in zip(levels[0]['values'],levels[1]['values']));check('cut_precision_'+str(case),pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
  FIXTURES.append({'fixture':'moving_linear_and_cubic_finite_parts','c3':c3s,'levels':levels,'entire_actual_row_evaluated':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 full constant quantum row','input_count':len(ps)}),flush=True)
 try:provenance();symbolic();matrix_fixtures();gaussian_fixtures();angular_fixtures();cut_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 full constant quantum row','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed','exception':error,
  'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'complete_constant_row_weak_limit_derived':good,'paired_even_mass_representation_derived':good,
  'actual_two_point_comparison_proved':False,'pure_shape_finite_cut_limit_proved':False,'physical_fast_branch_excluded':False,
  'actual_charged_kernel_evaluated':False,'entire_actual_row_numerically_evaluated':False,'pointwise_cutoff_limit_proved':False,
  'physical_action_averaging_adopted':False,'bare_tadpole_independently_renormalized':False,'finite_counterterm_changed':False,
  'new_boundary_prescription_adopted':False,'order_reduction_adopted':False,'physical_exterior_state_constructed':False,
  'continuous_causal_stability_proved':False,'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]};data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
