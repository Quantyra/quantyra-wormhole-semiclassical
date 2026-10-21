"""Diagnostics for the moving-cut canonical normal form, not physical stability."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import mpmath as mp
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM='smooth-local-quantum-comparison';PREVIOUS='reflected-boundary-kernel-limit'
OWN=['checks/smooth_local_quantum_comparison.py','notes/'+STEM+'.md','evidence/smooth-local-quantum-source-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-32','precision_tolerance':'1e-32',
 'exponential_pairs':[['1','2'],['.4','1.7'],['2.3','.8']],
 'Riesz_parameters':['.75','1.25','1.75'],
 'finite_part_cases':[['.03',3,'.7'],['.2',4,'1.3'],['.6',6,'2.3']],
 'cut_powers':[3,5,7],'density_epsilon':'.1','principal_constant':'.7',
 'Taylor_split':'0.125','Taylor_terms':40,'fixed_subtraction_scale':1,
 'ordinary_continuum_proof':True,'actual_charged_two_point_kernel_evaluated':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def zero(name,x):x=S.simplify(x);check(name,x==0,x)
def out(x):return mp.nstr(x,mp.mp.dps)
def err(a,b):return abs(a-b)/(1+abs(b))
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==76 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['canonical_normal_form_derived'] and not review['full_pure_shape_limit_proved'] and not review['physical_fast_branch_excluded'])

def symbolic():
 mx,my,vx,vy=S.symbols('mx my vx vy')
 zero('density_difference_square',mx*my*(vx-vy)**2-(mx*vx-my*vy)**2-mx*(my-mx)*vx**2-my*(mx-my)*vy**2)
 a,b,z=S.symbols('a b z',positive=True)
 R=a*b*(a**(3-2*z)-b**(3-2*z))/((b*b-a*a)*S.cos(S.pi*z))
 J=a*a*b*b/(a+b)-a*b*b-b*a*a
 L=a*b*(a**3*S.log(a)-b**3*S.log(b))/(b*b-a*a)
 zero('bilinear_J_contact',R.subs(z,0)-J)
 zero('bilinear_L_continuation',-S.diff(R,z).subs(z,0)/2-L)
 bulk=a**3*b**3*(S.log(a)/a-S.log(b)/b)/(b*b-a*a)
 zero('bilinear_L_boundary_terms',bulk-a*b*b*S.log(b)-b*a*a*S.log(a)-L)
 f2,fp=S.symbols('f2 fp')
 zero('Mellin_second_derivative_contact',S.expand((z-1)*(z-2)*(f2/(2*z)+fp)).coeff(z,0)-(2*fp-S.Rational(3,2)*f2))
 x=S.symbols('x');m=S.Function('m')(x);v=S.Function('v')(x)
 zero('density_principal_boundary',S.diff(m*v,x)*S.diff(m*v,x,2)-m*m*S.diff(v,x)*S.diff(v,x,2)
      -(m*S.diff(m,x)*v*S.diff(v,x,2)+2*m*S.diff(m,x)*S.diff(v,x)**2
      +(m*S.diff(m,x,2)+2*S.diff(m,x)**2)*v*S.diff(v,x)+S.diff(m,x)*S.diff(m,x,2)*v*v))
 zero('second_shape_contact',S.diff(v*v,x,2)-2*v*S.diff(v,x,2)-2*S.diff(v,x)**2)
 h0,h1,h2,h3,s=S.symbols('h0 h1 h2 h3 s')
 h=lambda t:h0+h1*t+h2*t*t/2+h3*t**3/6
 zero('midpoint_density_linear',S.expand((h(s/2)+h(-s/2))/2-h0).coeff(s,1))
 zero('midpoint_density_quadratic',S.expand((h(s/2)+h(-s/2))/2-h0)-h2*s*s/8)
 zero('H2_extension_value',3-2-1);zero('H2_extension_derivative',-3+4-1)

def compare_levels(name,levels,keys):
 with mp.workdps(80):
  pe=max(err(mp.mpf(levels[0][k]),mp.mpf(levels[1][k])) for k in keys)
  check(name+'_precision',pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))

def exponential_fixtures():
 for i,(aa,bb) in enumerate(PROTOCOL['exponential_pairs']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    a=mp.mpf(aa);b=mp.mpf(bb);values=[]
    for ss in PROTOCOL['Riesz_parameters']:
     s=mp.mpf(ss)
     low=mp.quad(lambda p:p**(4-2*s)/((p*p+a*a)*(p*p+b*b)),[0,1])
     # p=w^-2 makes every selected tail exponent a nonnegative integer.
     high=mp.quad(lambda w:2*w**(4*s-3)/((1+a*a*w**4)*(1+b*b*w**4)),[0,1])
     direct=2*a*b/mp.pi*(low+high)
     exact=a*b*(a**(3-2*s)-b**(3-2*s))/((b*b-a*a)*mp.cos(mp.pi*s))
     er=err(direct,exact);check(f'Riesz_integral_{i}_{digits}_{ss}',er<mp.mpf(PROTOCOL['identity_tolerance']),out(er));values.append(out(direct))
    low=mp.quad(lambda p:mp.log(p)/((p*p+a*a)*(p*p+b*b)),[0,1])
    high=mp.quad(lambda t:-t*t*mp.log(t)/((1+a*a*t*t)*(1+b*b*t*t)),[0,1])
    bulk=2*a**3*b**3/mp.pi*(low+high)
    from_boundary=bulk-a*b*b*mp.log(b)-b*a*a*mp.log(a)
    exact=a*b*(a**3*mp.log(a)-b**3*mp.log(b))/(b*b-a*a)
    er=err(from_boundary,exact);check(f'Fourier_bilinear_contact_{i}_{digits}',er<mp.mpf(PROTOCOL['identity_tolerance']),out(er))
    levels.append({'precision':digits,'L':out(from_boundary),'contact_error':out(er),'Riesz_values':values})
  compare_levels('exponential_'+str(i),levels,['L'])
  FIXTURES.append({'fixture':'bilinear_Riesz_and_Fourier_contacts','a':aa,'b':bb,'levels':levels})

def poly_value(c,x):return sum(a*x**j for j,a in enumerate(c))
def derivative(c,a):return [(j+1)*(c[j+1] if j+1<len(c) else 0)-a*c[j] for j in range(len(c))]
def taylor(c,a,n):return sum(c[k]*(-a)**(n-k)/mp.factorial(n-k) for k in range(min(n,len(c)-1)+1))
def moment(c,a):return sum(c[j]*mp.factorial(j)/a**(j+1) for j in range(len(c)))
def conv(c,e):return [sum(c[j]*e[k-j] for j in range(max(0,k-len(e)+1),min(k,len(c)-1)+1)) for k in range(len(c)+len(e)-1)]
def T_exact(c,a):return -c[0]*mp.log(a)+sum(c[j]*mp.factorial(j-1)/a**j for j in range(1,len(c)))
def T_direct(c,a):
 h=mp.mpf(PROTOCOL['Taylor_split'])
 near=sum((taylor(c,a,n)-c[0]*(-1)**n/mp.factorial(n))*h**n/n for n in range(1,PROTOCOL['Taylor_terms']+1))
 return near+mp.quad(lambda x:(poly_value(c,x)*mp.exp(-a*x)-c[0]*mp.exp(-x))/x,[h,1,mp.inf])
def FP3_direct(c,a):
 h=mp.mpf(PROTOCOL['Taylor_split']);f0=c[0];f1=taylor(c,a,1);f2=taylor(c,a,2)
 near=sum(taylor(c,a,n)*h**(n-2)/(n-2) for n in range(3,PROTOCOL['Taylor_terms']+1))
 mid=mp.quad(lambda x:(poly_value(c,x)*mp.exp(-a*x)-f0-f1*x-f2*x*x)/x**3,[h,1])
 far=mp.quad(lambda x:poly_value(c,x)*mp.exp(-a*x)/x**3,[1,mp.inf])
 return near+mid+far-f0/2-f1

def finite_part_fixtures():
 for i,(ds,k,aa) in enumerate(PROTOCOL['finite_part_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    d=mp.mpf(ds);a=mp.mpf(aa);c=[mp.binomial(k,j)*d**(k-j) for j in range(k+1)];c2=derivative(derivative(c,a),a)
    direct=T_direct(c2,a);fp=FP3_direct(c,a);rhs=2*fp+(mp.euler-mp.mpf('1.5'))*c2[0]
    er=err(direct,rhs);check(f'FP_second_derivative_{i}_{digits}',er<mp.mpf(PROTOCOL['identity_tolerance']),out(er))
    er2=err(direct,T_exact(c2,a));check(f'T_integral_{i}_{digits}',er2<mp.mpf(PROTOCOL['identity_tolerance']),out(er2))
    levels.append({'precision':digits,'T':out(direct),'FP3':out(fp),'identity_error':out(er),'T_error':out(er2)})
  compare_levels('finite_part_'+str(i),levels,['T','FP3'])
  FIXTURES.append({'fixture':'Mellin_derivative_contact_direct_subtraction','cut':ds,'degree':k,'rate':aa,'levels':levels})

def shifted(k,a,d,scale=1):return [scale*mp.exp(-a*d)*mp.binomial(k,j)*d**(k-j) for j in range(k+1)]
def cosine(c,a,p):return sum(c[k]*mp.factorial(k)*mp.re((a+1j*p)**(-k-1)) for k in range(len(c)))
def cosine_inverse(c,a,t):return sum(c[k]*mp.factorial(k)*mp.re((t/(a*t+1j))**(k+1)) for k in range(len(c)))
def LJexp(a,b):
 return (a*b*(a**3*mp.log(a)-b**3*mp.log(b))/(b*b-a*a),-a*b*(a*a+a*b+b*b)/(a+b))

def density_fixtures():
 levels=[]
 for digits in PROTOCOL['precision_digits']:
  with mp.workdps(digits):
   eps=mp.mpf(PROTOCOL['density_epsilon']);ell=mp.mpf(PROTOCOL['principal_constant']);values=[]
   for power in [None]+PROTOCOL['cut_powers']:
    d=mp.mpf(0) if power is None else mp.mpf(2)**(-power)
    f=shifted(2,mp.mpf(1),d,eps);f1=derivative(f,1);f2=derivative(f1,1)
    g=[(shifted(4,mp.mpf(2),d),mp.mpf(2)),(shifted(6,mp.mpf(3),d,eps),mp.mpf(3))]
    g1=[(derivative(c,a),a) for c,a in g];g2=[(derivative(c,a),a) for c,a in g1]
    low=mp.quad(lambda p:mp.log(p)*cosine(f2,1,p)*sum(cosine(c,a,p) for c,a in g2),[0,1])
    high=mp.quad(lambda t:-mp.log(t)*cosine_inverse(f2,1,t)*sum(cosine_inverse(c,a,t) for c,a in g2)/t**2,[0,1])
    L=2/mp.pi*(low+high)-f1[0]*sum(T_exact(c,a) for c,a in g2)-sum(c[0] for c,a in g1)*T_exact(f2,1)
    J=sum(moment(conv(f2,c),1+a) for c,a in g2)+f1[0]*sum(c[0] for c,a in g2)+sum(c[0] for c,a in g1)*f2[0]
    refs=[]
    for index in (0,1):
     ref=eps*mp.diff(lambda a,b:mp.exp(-(a+b)*d)*LJexp(a,b)[index],(mp.mpf(1),mp.mpf(2)),(2,4))
     ref+=eps*eps*mp.diff(lambda a,b:mp.exp(-(a+b)*d)*LJexp(a,b)[index],(mp.mpf(1),mp.mpf(3)),(2,6));refs.append(ref)
    erL=err(L,refs[0]);erJ=err(J,refs[1]);label='zero' if power is None else str(power)
    check(f'density_Fourier_vs_parameter_L_{digits}_{label}',erL<mp.mpf(PROTOCOL['identity_tolerance']),out(erL))
    check(f'density_bulk_vs_parameter_J_{digits}_{label}',erJ<mp.mpf(PROTOCOL['identity_tolerance']),out(erJ))
    values.append({'cut':out(d),'L':out(L),'J':out(J),'row':out(L-ell*J),'L_error':out(erL),'J_error':out(erJ)})
   limit=mp.mpf(values[0]['row']);first=mp.mpf(values[1]['row']);last=mp.mpf(values[-1]['row'])
   check('density_nonzero_limit_'+str(digits),abs(limit)>mp.mpf('.001'),out(limit))
   ratio=abs(last-limit)/abs(first-limit);check('density_cut_refinement_'+str(digits),ratio<mp.mpf('.125'),out(ratio))
   levels.append({'precision':digits,'values':values,'refinement_ratio':out(ratio)})
 with mp.workdps(80):
  pe=max(err(mp.mpf(a['row']),mp.mpf(b['row'])) for a,b in zip(levels[0]['values'],levels[1]['values']));check('density_precision',pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
 FIXTURES.append({'fixture':'nonzero_density_row_moving_cut','m':'1+epsilon Y^2 exp(-Y)','v':'Y^2 exp(-Y)','levels':levels,'actual_shape_solution':False})

def anchor_fixtures():
 levels=[]
 for digits in PROTOCOL['precision_digits']:
  with mp.workdps(digits):
   values=[]
   for power in [None]+PROTOCOL['cut_powers']:
    d=mp.mpf(0) if power is None else mp.mpf(2)**(-power);c=shifted(2,mp.mpf(1),d);c1=derivative(c,1);c2=derivative(c1,1);square=conv(c1,c1)
    direct=mp.quad(lambda x:mp.log(x)*poly_value(square,x)*mp.exp(-2*x),[0,1,mp.inf])
    exact=sum(square[k]*mp.factorial(k)/2**(k+1)*(mp.digamma(k+1)-mp.log(2)) for k in range(len(square)))
    er=err(direct,exact);label='zero' if power is None else str(power)
    check(f'log_anchor_integral_{digits}_{label}',er<mp.mpf(PROTOCOL['identity_tolerance']),out(er))
    collar=mp.quad(lambda x:-mp.log(x)*poly_value(square,x)*mp.exp(-2*x),[0,1])
    bound=moment(square,2)+moment(conv(c2,c2),2);ratio=collar/bound
    check(f'log_anchor_H2_bound_{digits}_{label}',0<ratio<1,out(ratio))
    values.append({'cut':out(d),'value':out(direct),'identity_error':out(er),'collar_bound_ratio':out(ratio)})
   levels.append({'precision':digits,'values':values})
 with mp.workdps(80):
  pe=max(err(mp.mpf(a['value']),mp.mpf(b['value'])) for a,b in zip(levels[0]['values'],levels[1]['values']));check('anchor_precision',pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
 FIXTURES.append({'fixture':'half_line_logarithmic_second_anchor','levels':levels})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 smooth local quantum comparison','input_count':len(ps)}),flush=True)
 try:provenance();symbolic();exponential_fixtures();finite_part_fixtures();density_fixtures();anchor_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 smooth local quantum comparison','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed','exception':error,
  'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'canonical_normal_form_derived':good,'density_row_retained':True,'fixed_derivative_contacts_retained':True,
  'full_pure_shape_limit_proved':False,'physical_fast_branch_excluded':False,'actual_charged_two_point_kernel_evaluated':False,
  'rough_finite_amplitude_family_differentiability_proved':False,'physical_action_averaging_adopted':False,
  'finite_counterterm_changed':False,'new_boundary_prescription_adopted':False,'order_reduction_adopted':False,
  'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]};data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
