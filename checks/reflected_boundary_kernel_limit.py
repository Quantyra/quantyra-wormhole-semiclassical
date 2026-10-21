"""Reflected corner correction diagnostics; smooth local comparison remains open."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import mpmath as mp
import sympy as S

ROOT=Path(__file__).resolve().parents[1];STEM='reflected-boundary-kernel-limit';PREVIOUS='moving-cut-constant-quantum-row'
OWN=['checks/reflected_boundary_kernel_limit.py','notes/'+STEM+'.md','evidence/reflected-boundary-source-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'Bessel_radii':['.01','.1','1','10','100'],
 'convolution_cases':[['.2','1'],['.7','2']],'tag_orders':[1,3],
 'finite_part_cases':[['.04','.2'],['.2','.35'],['.6','.8']],
 'cut_powers':[3,6,9],'spike_powers':[2,4,6],
 'Taylor_split':'1/8','Taylor_terms':40,'identity_tolerance':'1e-32',
 'precision_tolerance':'1e-32','actual_charged_two_point_kernel_evaluated':False,
 'ordinary_continuum_proof':True,'fixed_subtraction_scale':1,
 'angular_convolution_coordinate':'cos(theta)=1-w^2, 0<=w<=1'}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def err(a,b):return abs(a-b)/(1+abs(b))
def out(x):return mp.nstr(x,mp.mp.dps)
def zero(name,x):x=S.simplify(x);check(name,x==0,x)
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==247 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',not review['smooth_local_to_normal_comparison_proved'] and not review['physical_fast_branch_excluded'])

def symbolic():
 r,mu=S.symbols('r mu',positive=True)
 zero('linear_tag_radial_reserve',2*S.pi*S.integrate(S.exp(-mu*r/2),(r,0,S.oo))-4*S.pi/mu)
 zero('cubic_tag_radial_reserve',2*S.pi*S.integrate(r*r*S.exp(-mu*r/2),(r,0,S.oo))-32*S.pi/mu**3)
 zero('coincident_boundary_linear_moment',2*S.integrate(S.exp(-2*mu*r),(r,0,S.oo))-1/mu)
 zero('coincident_boundary_cubic_moment',S.Rational(4,3)*S.integrate(r*r*S.exp(-2*mu*r),(r,0,S.oo))-1/(3*mu**3))
 u=S.symbols('u',positive=True)
 zero('corner_angular_mass',S.integrate(1-S.log(u),(u,0,1))-2)
 zero('zero_cut_corner_value',S.integrate((1-S.log(u))*(3*u*u-u**4)/8,(u,0,1))-S.Rational(41,300))
 a,b,h=S.symbols('a b h');p=S.Function('p')
 F=(p(a)-p(b))**2
 zero('corner_value',F.subs({a:0,b:0}))
 zero('corner_first_jet',S.diff(F,a).subs({a:0,b:0}))
 zero('corner_second_jet',S.diff(F,a,2).subs({a:0,b:0})-2*S.Subs(S.Derivative(p(a),a),a,0)**2)
 zero('corner_mixed_jet',S.diff(F,a,b).subs({a:0,b:0})+2*S.Subs(S.Derivative(p(a),a),a,0)**2)
 x=S.symbols('x');bump=x*x*(1-x*x)**6/2
 zero('spike_curvature',S.diff(bump,x,2).subs(x,0)-1)
 zero('spike_H2_exponent',2*13-3*8-2)
 return S.lambdify(x,[bump,S.diff(bump,x),S.diff(bump,x,2)],'mpmath'),S.lambdify((),S.integrate(S.diff(bump,x,2)**2,(x,-1,1)),'mpmath')

def majorant_fixtures():
 levels=[]
 for digits in PROTOCOL['precision_digits']:
  with mp.workdps(digits):
   L=2*mp.sqrt(2)/(1/mp.sqrt(2)-mp.mpf('.25'));epsilon=mp.mpf(1)/28
   check('weighted_L1_reserve_'+str(digits),L<7,out(L))
   series=sum((n+1)*(epsilon*L)**n for n in range(256));closed=1/(1-epsilon*L)**2
   check('Born_series_'+str(digits),err(series,closed)<mp.mpf(PROTOCOL['identity_tolerance']))
   factor=mp.sqrt(2)/mp.pi*closed;check('full_kernel_majorant_'+str(digits),factor<1,out(factor))
   reserve=64+512/mp.e**2;check('crossing_cubic_reserve_'+str(digits),reserve<256,out(reserve))
   ratios=[]
   for rs in PROTOCOL['Bessel_radii']:
    r=mp.mpf(rs);actual=(mp.besselk(0,r)+mp.besselk(1,r))/(2*mp.pi)
    bound=mp.sqrt(2)/mp.pi*mp.exp(-r/mp.sqrt(2))/r;ratio=actual/bound
    check('flat_Dirac_envelope_'+str(digits)+'_'+rs,ratio<=1,out(ratio));ratios.append(out(ratio))
   levels.append({'precision':digits,'weighted_L1':out(L),'Born_factor':out(factor),'Bessel_ratios':ratios,'crossing_cubic_reserve':out(reserve)})
 FIXTURES.append({'fixture':'Dirac_Born_majorant','levels':levels})
 for case,(us,ms) in enumerate(PROTOCOL['convolution_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    u=mp.mpf(us);mu=mp.mpf(ms);values=[]
    for j in PROTOCOL['tag_orders']:
     def radial(c):
      if c==0:return mp.mpf(0)
      assert c>0
      return sum(mp.binomial(j,k)*c**k*(-u)**(j-k)*mp.gammainc(k,2*mu*u/c,mp.inf)/(2*mu)**k for k in range(j+1))
     value=4*mp.quad(lambda w:radial(1-w*w)/mp.sqrt(2-w*w),[0,1])
     # The radial integral was evaluated analytically; retain the complete
     # angular integral and compare with the independent envelope (4).
     bound=(4*mp.pi/mu if j==1 else 32*mp.pi/mu**3)*mp.exp(-mu*u)
     ratio=value/bound;check(f'tagged_convolution_{case}_{digits}_{j}',0<ratio<1,out(ratio))
     values.append({'j':j,'value':out(value),'bound_ratio':out(ratio)})
    levels.append({'precision':digits,'values':values})
  with mp.workdps(80):
   pe=max(err(mp.mpf(a['value']),mp.mpf(b['value'])) for a,b in zip(levels[0]['values'],levels[1]['values']));check('convolution_precision_'+str(case),pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
  FIXTURES.append({'fixture':'tagged_half_plane_convolution','distance':us,'mu':ms,'levels':levels})

def coefficients(d,t):
 a=t;b=1-t
 ca=[mp.binomial(4,k)*d**(4-k)*a**k for k in range(5)]
 cb=[mp.binomial(4,k)*d**(4-k)*b**k for k in range(5)]
 cc=[-2*sum(mp.binomial(2,i)*mp.binomial(2,k-i)*d**(4-k)*a**i*b**(k-i) for i in range(max(0,k-2),min(2,k)+1)) for k in range(5)]
 return ca,cb,cc
def gamma_fp(k,a):
 if k>=3:return mp.gamma(k-2)/a**(k-2)
 m=2-k;return (-1)**m/mp.factorial(m)*a**m*(mp.harmonic(m)-mp.euler-mp.log(a))
def analytic_fp(d,t):
 if t==mp.mpf('.5'):return mp.mpf(0)
 ca,cb,cc=coefficients(d,t)
 if t==0:return mp.exp(-2*d)*sum(cb[k]*gamma_fp(k,2)+cc[k]*gamma_fp(k,1) for k in range(5))
 if t==1:return mp.exp(-2*d)*sum(ca[k]*gamma_fp(k,2)+cc[k]*gamma_fp(k,1) for k in range(5))
 return mp.exp(-2*d)*sum(ca[k]*gamma_fp(k,2*t)+cb[k]*gamma_fp(k,2-2*t)+cc[k]*gamma_fp(k,1) for k in range(5))
def direct_fp(d,t):
 ca,cb,cc=coefficients(d,t);h=mp.mpf(1)/8
 series=mp.mpf(0)
 for n in range(3,PROTOCOL['Taylor_terms']+1):
  fn=mp.exp(-2*d)*sum((ca[k]*(-2*t)**(n-k)+cb[k]*(-2+2*t)**(n-k)+cc[k]*(-1)**(n-k))/mp.factorial(n-k) for k in range(min(4,n)+1))
  series+=fn*h**(n-2)/(n-2)
 f2=(2*t-1)**2*mp.exp(-2*d)*(2*d-d*d)**2
 def F(s):return ((d+t*s)**2*mp.exp(-d-t*s)-(d+(1-t)*s)**2*mp.exp(-d-(1-t)*s))**2
 return series+mp.quad(lambda s:(F(s)-f2*s*s)/s**3,[h,1])+mp.quad(lambda s:F(s)/s**3,[1,mp.inf])

def finite_part_fixtures():
 for case,(ds,ts) in enumerate(PROTOCOL['finite_part_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    d=mp.mpf(ds);t=mp.mpf(ts);direct=direct_fp(d,t);analytic=analytic_fp(d,t);er=err(direct,analytic)
    check(f'corner_finite_part_{case}_{digits}',er<mp.mpf(PROTOCOL['identity_tolerance']),out(er))
    levels.append({'precision':digits,'value':out(analytic),'scaled_error':out(er)})
  with mp.workdps(80):
   pe=err(mp.mpf(levels[0]['value']),mp.mpf(levels[1]['value']));check('finite_part_precision_'+str(case),pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
  FIXTURES.append({'fixture':'corner_Fourier_free_finite_part','cut':ds,'angular_ratio':ts,'levels':levels})
 levels=[]
 for digits in PROTOCOL['precision_digits']:
  with mp.workdps(digits):
   values=[];limit=mp.mpf(41)/300
   for power in PROTOCOL['cut_powers']:
    d=mp.mpf(2)**(-power)
    # Symmetry removes the angular logarithm at t=1/2 without sampling it.
    value=mp.quad(lambda u:(1-mp.log(u))*analytic_fp(d,(1+u)/2) if u else mp.mpf(0),[0,mp.mpf('.5'),1])
    scaled=abs(value-limit)/d;check(f'ordinary_corner_cut_{digits}_{power}',scaled<5,out(scaled))
    values.append({'power':power,'value':out(value),'limit_error_over_cut':out(scaled)})
   check('corner_cut_refinement_'+str(digits),abs(mp.mpf(values[-1]['value'])-limit)<abs(mp.mpf(values[0]['value'])-limit)/16)
   levels.append({'precision':digits,'limit':out(limit),'values':values})
 with mp.workdps(80):
  pe=max(err(mp.mpf(a['value']),mp.mpf(b['value'])) for a,b in zip(levels[0]['values'],levels[1]['values']));check('corner_cut_precision',pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
 FIXTURES.append({'fixture':'ordinary_logarithmic_corner_cut_limit','levels':levels})

def spike_fixtures(bump,bumpnorm):
 levels=[]
 for digits in PROTOCOL['precision_digits']:
  with mp.workdps(digits):
   values=[]
   for power in PROTOCOL['spike_powers']:
    dn=mp.mpf(2)**(-power);ell=dn**8;amp=dn**13;z0=mp.mpf('-.5');d=dn+z0*ell
    alpha=mp.mpf('.8');beta=mp.mpf('.2');rt=alpha-beta;t=mp.mpf('.7');s=ell*t
    f=2*d+amp/ell*bump(z0)[1]
    e=rt*s*s+amp*(bump(z0+alpha*t)[0]-bump(z0+beta*t)[0]-rt*t*bump(z0)[1])
    rem=abs(2*rt*s*f*e+e*e)
    en2=ell*mp.quad(lambda z:(2+amp/ell**2*bump(z)[2])**2,[z0,z0+t]);en=mp.sqrt(en2)
    bound=mp.mpf(4)/3*s**mp.mpf('2.5')*abs(f)*en+mp.mpf(4)/9*s**3*en2
    ratio=rem/bound;check(f'concentrated_H2_remainder_{digits}_{power}',ratio<=1,out(ratio))
    cost=amp*amp/ell**3*bumpnorm();check(f'concentrated_energy_cost_{digits}_{power}',err(cost,dn*dn*bumpnorm())<mp.mpf('1e-40'))
    peak=2+dn**(-3);values.append({'power':power,'second_derivative_peak':out(peak),'H2_cost':out(cost),'remainder_ratio':out(ratio)})
   levels.append({'precision':digits,'values':values})
 with mp.workdps(80):
  pe=max(err(mp.mpf(a['remainder_ratio']),mp.mpf(b['remainder_ratio'])) for a,b in zip(levels[0]['values'],levels[1]['values']));check('spike_precision',pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
 FIXTURES.append({'fixture':'concentrated_H2_corner_remainder','levels':levels,'actual_shape_solution':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 reflected boundary kernel limit','input_count':len(ps)}),flush=True)
 try:provenance();bump,bumpnorm=symbolic();majorant_fixtures();finite_part_fixtures();spike_fixtures(bump,bumpnorm)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 reflected boundary kernel limit','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed','exception':error,
  'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'reflected_boundary_correction_limit_derived':good,'corner_contacts_retained':True,
  'smooth_local_to_normal_comparison_proved':False,'pure_shape_finite_cut_limit_proved':False,'physical_fast_branch_excluded':False,
  'actual_charged_two_point_kernel_evaluated':False,'rough_finite_amplitude_family_differentiability_proved':False,
  'physical_action_averaging_adopted':False,'finite_counterterm_changed':False,'new_boundary_prescription_adopted':False,
  'order_reduction_adopted':False,'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]};data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
