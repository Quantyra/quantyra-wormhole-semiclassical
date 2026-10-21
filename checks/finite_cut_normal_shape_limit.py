"""Moving-cut normal shape model; actual charged finite-cut comparison remains open."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='finite-cut-normal-shape-limit';PREVIOUS='energy-shape-nonlinear-trace-forcing'
OWN=['checks/finite_cut_normal_shape_limit.py','notes/'+STEM+'.md','evidence/finite-cut-normal-source-review.json']
ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'precision_digits':[50,80],'frequencies':['0','.5','2','5'],
 'exponential_rates':['1','2','5'],'cut_powers':[8,12,16],
 'angular_q':['1','8','16000000000000'],'angular_t_over_q2':['1e-8','.01','.25','1','4','100'],
 'angular_tail_N':256,'contour_radius':'.01','contour_nodes':[24,48],
 'quadrature_scaled_limit':'1e-30','precision_scaled_limit':'1e-32',
 'moving_cut_sqrt_bound':'20','moving_cut_fine_error_limit':'.01',
 'model_only':True,'actual_charged_hessian_evaluated':False,
 'normal_limit_is_ordinary_argument':True,'physical_action_averaging_adopted':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)}
 ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.factor(S.cancel(S.sympify(value).doit()));check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==230 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('scope',not rev['actual_finite_cut_comparison_proved'] and not rev['actual_energy_domain_transfer_proved'] and not rev['physical_fast_branch_excluded'])

def symbolic():
 eta,sigma=S.symbols('eta sigma');T=S.symbols('T:4')
 def D(f):return sum(S.diff(f,T[j])*T[j+1] for j in range(3))
 def coefficient(f):return S.expand(S.diff(f.subs({t:eta*t for t in T},simultaneous=True),eta,2).subs(eta,0)/2)
 scalar=-2*S.exp(-2*T[0])*T[2];K=S.exp(-T[0])*T[1]
 rn=-S.exp(-T[0])*D(scalar);heat_section=-scalar/2
 boundary=S.exp(T[0])*(-3*rn+5*scalar*K+22*heat_section*K-R(40,7)*K**3)
 box_flux=-3*D(scalar)
 zero('heat_boundary_quadratic',coefficient(boundary)-12*T[0]*T[3]-24*T[1]*T[2])
 zero('heat_bulk_flux_quadratic',coefficient(box_flux)+12*T[0]*T[3]+12*T[1]*T[2])
 zero('heat_net_boundary',coefficient(boundary+box_flux)-12*T[1]*T[2])
 zero('heat_bulk_fourth',coefficient(-R(3,2)*S.exp(2*T[0])*scalar**2)+6*T[2]**2)
 proper_a=K;proper_ap=S.exp(-T[0])*D(proper_a)
 wz=S.exp(T[0])*sigma*(12*proper_a*proper_ap+R(44,7)*proper_a**3)
 zero('independent_Weyl_boundary_log',coefficient(wz)-12*sigma*T[1]*T[2])
 z=S.symbols('z');V=S.Function('V')(z);A=R(3,2)*V-S.diff(V,z)
 g=S.diff(V,z,2)-2*S.diff(V,z)+R(3,4)*V
 zero('critical_profile_Hardy_identity',S.diff(A,z)-A/2+g)
 u,v=S.symbols('u v',positive=True)
 phi=2*u/(u*u+v);psi=2*u/(u*u+v)**2
 zero('phi_second_derivative',S.diff(phi,u,2)-4*u*(u*u-3*v)/(u*u+v)**3)
 zero('psi_second_derivative',S.diff(psi,u,2)-24*u*(u*u-v)/(u*u+v)**4)
 zero('phi_derivative_max',S.diff(phi,u).subs(u,0)-2/v)
 zero('phi_derivative_min',S.diff(phi,u).subs(u,S.sqrt(3*v))+1/(4*v))
 zero('psi_derivative_max',S.diff(psi,u).subs(u,0)-2/v**2)
 zero('psi_derivative_min',S.diff(psi,u).subs(u,S.sqrt(v))+1/(2*v**2))
 zero('trapezoid_phi_constant',R(1,4)*R(5,2)*R(4,3)-R(5,6))
 zero('trapezoid_psi_constant',R(1,4)*3*R(4,3)**2-R(4,3))
 alpha=S.symbols('alpha');zero('angular_integral_constant',8*S.integrate(alpha*(1-alpha),(alpha,0,1))-R(4,3))
 p,k,d=S.symbols('p k d',positive=True)
 fourier=d*d*k/(k*k+p*p)+2*d*(k*k-p*p)/(k*k+p*p)**2+2*(k**3-3*k*p*p)/(k*k+p*p)**3
 ap=2*d-k*d*d
 reg=S.factor(S.cancel(p**4*fourier**2-ap**2));hf=S.factor(S.cancel(-p*p*fourier-ap))
 zero('Fourier_cusp_tail',S.limit(p*p*fourier,p,S.oo)+ap)
 zero('subtracted_tail_integrable',S.limit(reg,p,S.oo))
 return S.lambdify((p,k,d),reg,'mpmath'),S.lambdify((p,k,d),hf,'mpmath')

def err(a,b):return abs(a-b)/(1+abs(b))
def string(x):return mp.nstr(x,mp.mp.dps)
def kernel_fixtures():
 for xs in PROTOCOL['frequencies']:
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    x=mp.mpf(xs);splits=[0,1,2,4,8,16,32,64,128,mp.inf]
    kern=lambda r:1/(2*mp.sinh(r/2))
    direct=mp.digamma(mp.mpf('.5'))+mp.quad(lambda r:kern(r)*(-mp.expm1(-1j*x*r)),splits)
    target=mp.digamma(mp.mpf('.5')+1j*x)
    car=2*mp.quad(lambda r:mp.cos(x*r)/(2*mp.cosh(r/2)),splits)
    car_target=mp.pi/mp.cosh(mp.pi*x)
    e=max(err(direct,target),err(car,car_target))
    check('kernel_'+xs+'_'+str(digits),e<mp.mpf(PROTOCOL['quadrature_scaled_limit']),mp.nstr(e,25))
    check('boundary_multiplier_bound_'+xs+'_'+str(digits),abs(target)<=2+5*abs(x))
    def riesz(s):
     return (mp.gamma(mp.mpf('.5')-1j*x)/mp.gamma(mp.mpf('.5')-1j*x+2*s)+mp.gamma(mp.mpf('.5')+1j*x-2*s)/mp.gamma(mp.mpf('.5')+1j*x))/(2*mp.cos(mp.pi*s))
    contours=[]
    for n in PROTOCOL['contour_nodes']:
     rad=mp.mpf(PROTOCOL['contour_radius'])
     der=mp.fsum(riesz(rad*mp.exp(2j*mp.pi*j/n))*mp.exp(-2j*mp.pi*j/n) for j in range(n))/(n*rad)
     ce=err(-der/2,mp.re(target));check('Riesz_contour_'+xs+'_'+str(digits)+'_'+str(n),ce<mp.mpf('1e-25'),mp.nstr(ce,25))
     contours.append({'nodes':n,'scaled_error':mp.nstr(ce,25)})
    levels.append({'precision':digits,'real_symbol':string(mp.re(direct)),'imaginary_symbol':string(mp.im(direct)),'Carleman':string(car),'scaled_error':mp.nstr(e,25),'contours':contours})
  with mp.workdps(80):
   pe=err(mp.mpf(levels[0]['real_symbol']),mp.mpf(levels[1]['real_symbol']))
   check('kernel_precision_'+xs,pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,25))
  FIXTURES.append({'fixture':'Mellin_and_boundary_kernels','frequency':xs,'levels':levels,'actual_quantum_response':False})
 with mp.workdps(80):
  norm=mp.quad(lambda r:r/(2*mp.sinh(r/2)),[0,1,4,16,mp.inf])
  check('crossing_kernel_L1_mass',err(norm,mp.pi**2/2)<mp.mpf('1e-60'),mp.nstr(norm,25))

def cut_values(k,d,reg,hf):
 fac=2/mp.pi*mp.exp(-2*k*d);splits=[0,1,k+1,mp.inf]
 J=fac*mp.quad(lambda p:reg(p,k,d),splits)
 L=fac*mp.quad(lambda p:mp.log(p)*reg(p,k,d),splits)
 H=fac*mp.quad(lambda p:mp.log(p)*hf(p,k,d)**2,splits)
 a=mp.exp(-k*d)*(2*d-k*d*d);bp=k*k*d*d-4*k*d+2;b=mp.exp(-k*d)*bp
 t=mp.exp(-k*d)*mp.quad(lambda x:bp*mp.exp(-x)*mp.expm1((1-k)*x)/x+mp.exp(-k*x)*(2*k*k*d-4*k+k*k*x),[0,1,mp.inf])
 h2=lambda y:mp.exp(-2*k*y)*(k*k*y*y-4*k*y+2)**2
 mass=mp.quad(h2,[d,1,mp.inf]);bulk=mp.quad(lambda y:mp.log(y)*h2(y),[d,1,mp.inf])
 contact=2*mp.log(d)*a*b if d else mp.mpf(0)
 return {'J':J,'L':L,'decomposition_error':err(L,H-2*a*t),'local_contact_error':err(J,mass+2*a*b),'P':L+bulk+contact,'boundary':-2*a*t+contact,'mass':mass}

def cut_fixtures(reg,hf):
 for ks in PROTOCOL['exponential_rates']:
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    k=mp.mpf(ks)
    density=lambda x:(x*x+mp.mpf('.25'))*(x*x+mp.mpf('2.25'))/(2*k*mp.cosh(mp.pi*x))
    symbol=lambda x:mp.re(mp.digamma(mp.mpf('.5')+1j*x))-mp.pi/(2*mp.cosh(mp.pi*x))
    limit=2*mp.quad(lambda x:symbol(x)*density(x),[0,1,3,10,mp.inf])
    origin=cut_values(k,mp.mpf(0),reg,hf)
    e0=err(origin['P'],limit);check('Mellin_limit_'+ks+'_'+str(digits),e0<mp.mpf(PROTOCOL['quadrature_scaled_limit']),mp.nstr(e0,25))
    values=[]
    for n in PROTOCOL['cut_powers']:
     d=mp.mpf(2)**(-n);v=cut_values(k,d,reg,hf)
     e=max(v['decomposition_error'],v['local_contact_error'])
     check('finite_part_'+ks+'_'+str(digits)+'_'+str(n),e<mp.mpf(PROTOCOL['quadrature_scaled_limit']),mp.nstr(e,25))
     ce=abs(v['P']-limit);ratio=ce/mp.sqrt(d)
     check('cut_diagnostic_'+ks+'_'+str(digits)+'_'+str(n),ratio<mp.mpf(PROTOCOL['moving_cut_sqrt_bound']),mp.nstr(ratio,25))
     if n==PROTOCOL['cut_powers'][-1]:check('fine_cut_'+ks+'_'+str(digits),ce<mp.mpf(PROTOCOL['moving_cut_fine_error_limit']),mp.nstr(ce,25))
     values.append({'power':n,'P':string(v['P']),'boundary':string(v['boundary']),'absolute_limit_error':mp.nstr(ce,25),'sqrt_scaled_error':mp.nstr(ratio,25),'finite_part_scaled_error':mp.nstr(e,25)})
    scale_cut=mp.mpf(2)**(-PROTOCOL['cut_powers'][0])
    scaled=cut_values(mp.mpf(1),k*scale_cut,reg,hf)['P']
    se=err(k*mp.mpf(values[0]['P']),scaled)
    check('dilation_'+ks+'_'+str(digits),se<mp.mpf(PROTOCOL['quadrature_scaled_limit']),mp.nstr(se,25))
    levels.append({'precision':digits,'limit':string(limit),'limit_scaled_error':mp.nstr(e0,25),'dilation_scaled_error':mp.nstr(se,25),'values':values})
  with mp.workdps(80):
   pe=max(err(mp.mpf(a['P']),mp.mpf(b['P'])) for a,b in zip(levels[0]['values'],levels[1]['values']))
   check('cut_precision_'+ks,pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,25))
  FIXTURES.append({'fixture':'moving_cut_polynomial_exponential','rate':ks,'profile':'v(y)=y^2 exp(-k y)','levels':levels,'actual_shape_solution':False})

def angular_values(q,t):
 w=mp.sqrt(q*q/4-t);xp=1+q/2+w;xm=1+q/2-w
 F=-2*(mp.digamma(xp)+mp.digamma(xm))
 Fp=(mp.polygamma(1,xp)-mp.polygamma(1,xm))/w if w else 2*mp.polygamma(2,1+q/2)
 return mp.re(F),mp.re(Fp)

def angular_fixtures():
 for qs in PROTOCOL['angular_q']:
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    q=mp.mpf(qs);values=[]
    for rs in PROTOCOL['angular_t_over_q2']:
     t=q*q*mp.mpf(rs);F,Fp=angular_values(q,t)
     ratio=max(abs(F+2*mp.log(t)),abs(t*Fp+2))*t/(8*q)
     check('angular_bound_'+qs+'_'+str(digits)+'_'+rs,ratio<=1,mp.nstr(ratio,25))
     tailratio=None
     if q<=8:
      N=PROTOCOL['angular_tail_N'];u=mp.mpf(N)+q/2;vv=t-q*q/4
      ph=lambda n:(2*n+q)/(n*(n+q)+t)
      ps=lambda n:(2*n+q)/(n*(n+q)+t)**2
      approx=2*mp.fsum(ph(n) for n in range(1,N+1))-2*mp.log(N*(N+q)+t)-ph(N)
      approxd=-2*mp.fsum(ps(n) for n in range(1,N+1))-2/(N*(N+q)+t)+ps(N)
      eb=-2*(vv-u*u)/(u*u+vv)**2/4
      ed=-2*(vv-3*u*u)/(u*u+vv)**3/4
      tailratio=max(abs(F-approx)/eb,abs(Fp-approxd)/ed)
      check('angular_independent_sum_'+qs+'_'+str(digits)+'_'+rs,tailratio<=1,mp.nstr(tailratio,25))
     values.append({'t_over_q2':rs,'F':string(F),'Fprime':string(Fp),'bound_ratio':mp.nstr(ratio,25),'tail_enclosure_ratio':mp.nstr(tailratio,25) if tailratio is not None else None})
    levels.append({'precision':digits,'values':values})
  with mp.workdps(80):
   pe=max(max(err(mp.mpf(a['F']),mp.mpf(b['F'])),err(mp.mpf(a['Fprime']),mp.mpf(b['Fprime']))) for a,b in zip(levels[0]['values'],levels[1]['values']))
   check('angular_precision_'+qs,pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,25))
  FIXTURES.append({'fixture':'full_angular_sum_remainder','q':qs,'levels':levels,'curved_quantum_comparison':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 finite-cut normal shape limit','input_count':len(ps)}),flush=True)
 try:
  provenance();reg,hf=symbolic();kernel_fixtures();cut_fixtures(reg,hf);angular_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 finite-cut normal shape limit','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'principal_boundary_contact_derived':good,'normal_model_weak_limit_derived':good,'full_flat_tower_remainder_improved':good,
  'actual_finite_cut_comparison_proved':False,'actual_energy_domain_transfer_proved':False,
  'pure_shape_finite_cut_limit_proved':False,'full_boundary_action_solved':False,'physical_action_averaging_adopted':False,
  'physical_fast_branch_excluded':False,'new_boundary_prescription_adopted':False,'actual_second_vertex_evaluated':False,
  'rough_finite_amplitude_action_differentiability_proved':False,'full_shape_equation_solved':False,
  'static_source_compatibility_proved':False,'regulator_independence_for_moving_cuts_proved':False,
  'continuous_causal_stability_proved':False,'physical_constraint_mode_constructed':False,'actual_wormhole_instability_proved':False,
  'finite_counterterm_changed':False,'order_reduction_adopted':False,'formal_verification':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in [('checker',OWN[0]),('note',OWN[1]),('review',OWN[2])]:run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
