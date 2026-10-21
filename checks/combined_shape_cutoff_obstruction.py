"""Complete local shape weights and conditional cutoff obstruction diagnostics."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='combined-shape-cutoff-obstruction';PREVIOUS='smooth-local-quantum-comparison'
OWN=['checks/combined_shape_cutoff_obstruction.py','notes/'+STEM+'.md','evidence/combined-shape-cutoff-source-review.json']
ROWS=[];FIXTURES=[];AUDITS=[];R=S.Rational
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-32','precision_tolerance':'1e-32',
 'jet_cases':[['.12','.2','-.1','.03','-.04','.05','.3','-.2','1.1','-.9','1.3','.1','2'],
              ['-.2','.1','.3','.07','.12','-.08','-.1','.4','2.2','-1.4','2.1','.2','3']],
 'jet_columns':['s0','s1','s2','v0','v1','v2','beta','beta_prime','omega','omega_prime','omega_second','sphere_curvature','magnetic_flux'],
 'contour_nodes':[24,32],'contour_radius':'.03','contour_tolerance':'1e-30',
 'profiles':[['.3','.4','2','3'],['.6','.5','3','1.5']],
 'profile_columns':['shape_s_extra_decay','shape_v_extra_decay','shape_s_log_frequency','shape_v_log_frequency'],
 'cut_Z':[2,5,8],'spike_powers':[3,6,9], 'trace_t':['.05','.2','1'],
 'fixed_frequency_transfer':'exact bulk frequency polynomial and full boundary derivative weights',
 'actual_shape_solution_evaluated':False,'ordinary_continuum_proof':True,
 'physical_fast_branch_excluded':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def zero(name,value):value=S.factor(S.cancel(S.sympify(value).doit()));check(name,value==0,value)
def out(x):return mp.nstr(x,mp.mp.dps)
def err(a,b):return abs(a-b)/(1+abs(b))
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==96 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes());check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['full_pure_shape_weak_limit_derived'] and review['conditional_fast_expansion_obstructed'] and not review['physical_fast_branch_excluded'])

def weight_audit(label,expr,shapes,background,weights,dimension,log_symbol=None):
 variables=list(shapes)+list(background);poly=S.Poly(S.expand(expr),*variables);records=[]
 for powers,coef in poly.terms():
  sp=powers[:len(shapes)];bp=powers[len(shapes):]
  degree=sum(sp);jet=sum((j%3)*n for j,n in enumerate(sp));weight=sum(weights[x]*n for x,n in zip(background,bp))
  log_power=0 if log_symbol is None else bp[list(background).index(log_symbol)]
  extra=weight-jet+dimension
  records.append({'monomial':str(S.Mul(coef,*(x**n for x,n in zip(variables,powers) if n))),
   'shape_degree':degree,'shape_derivative_weight':jet,'coefficient_power':weight,'extra_power':extra,'log_power':log_power})
 check(label+'_quadratic',all(r['shape_degree']==2 for r in records),len(records))
 check(label+'_weights',all(r['extra_power']>=(2 if r['log_power'] else 0) for r in records),min(r['extra_power'] for r in records))
 check(label+'_log_degree',all(r['log_power']<=1 for r in records))
 AUDITS.append({'part':label,'dimension':dimension,'monomials':records})

def derive():
 eta=S.symbols('eta');s=S.symbols('s:3');v=S.symbols('v:3');b,b1,W,W1,W2,ks,q=S.symbols('b b1 W W1 W2 ks q')
 shapes=[*s,*v];bg=[b,b1,W,W1,W2,ks,q];weights={b:1,b1:0,W:0,W1:-1,W2:-2,ks:0,q:0}
 e=S.exp(-eta*v[0]);A=-e*(b+eta*v[1]);Ax=e**2*(b1+eta*v[2]-eta*v[1]*(b+eta*v[1]))
 T=-e*(W1+eta*s[1]);Tx=e**2*(W2+eta*s[2]-eta*v[1]*(W1+eta*s[1]))
 Rh=-2*e**2*(b1+b*b+eta*(v[2]+b*v[1]));U=e**2*(W2+b*W1+eta*(s[2]+b*s[1]));V=e**2*(W1+eta*s[1])**2
 L=(W+eta*s[0])*((Rh+2*ks)**2/60-R(11,90)*ks*Rh-q*q*ks*ks/3)+Rh*U/30-U*U/10+(Rh/30-ks/18)*V-R(7,90)*U*V-R(7,180)*V*V
 GG=(W+eta*s[0])*(12*A*Ax+R(44,7)*A**3+10*ks*A)-6*Ax*T-R(24,7)*A*A*T-6*A*Tx-R(96,7)*A*T*T-18*T*Tx-R(85,21)*T**3
 rawbulk=-R(15,2)*S.exp(2*eta*v[0])*L;rawboundary=S.exp(eta*v[0])*GG/12
 hess=lambda f:S.expand(S.diff(f,eta,2).subs(eta,0))
 qb=hess(rawbulk);qt=hess(rawboundary)
 rembulk=S.expand(qb+W*(v[2]+b*v[1])**2);remboundary=S.expand(qt+2*W*v[1]*v[2])
 zero('bulk_highest_log_coefficient',S.diff(qb,W,v[2],v[2])+2)
 zero('boundary_highest_log_coefficient',S.diff(qt,W,v[1],v[2])+2)
 weight_audit('Weyl_bulk_remainder',rembulk,shapes,bg,weights,4,W)
 weight_audit('Weyl_boundary_remainder',remboundary,shapes,bg,weights,3,W)
 E1,E2,E3,E4,lam0,Q,c0,c1=S.symbols('E1 E2 E3 E4 lam0 Q c0 c1')
 rawC=E2*S.exp(2*eta*s[0])*(-2*(b1+b*b)-2*eta*(v[2]+b*v[1])-6*(W2+b*W1+eta*(s[2]+b*s[1]))-6*(W1+eta*s[1])**2)
 rawC+=2*ks*E2*S.exp(2*eta*(s[0]+v[0]))+6*lam0*E4*S.exp(eta*(4*s[0]+2*v[0]))-2*Q*S.exp(2*eta*v[0])
 rawBC=-E2*S.exp(2*eta*s[0])*(b+3*W1+eta*(v[1]+3*s[1]))+c0*E3*S.exp(eta*(3*s[0]+v[0]))+c1*E1*S.exp(eta*(s[0]+v[0]))
 more=[E1,E2,E3,E4,lam0,Q,c0,c1];cw={**weights,E1:-1,E2:-2,E3:-3,E4:-4,lam0:0,Q:0,c0:0,c1:0}
 weight_audit('classical_bulk',hess(rawC),shapes,bg+more,cw,4,W)
 weight_audit('GHY_and_intrinsic_counterterms',hess(rawBC),shapes,bg+more,cw,3,W)
 # Audit the full retained optical frequency polynomial, without using a
 # static boundary polynomial as a dynamic formula.
 wb2=S.symbols('wb2');J1=-24*W*wb2*v[0]**2+12*W*v[0]*v[2]+W1**2*(7*s[0]**2+14*s[0]*v[0]+12*v[0]**2)
 J1+=-14*W1*s[0]*s[1]-12*W1*s[1]*v[0]+W2*(7*s[0]**2+36*s[0]*v[0]+12*v[0]**2)+wb2*(6*s[0]**2+24*s[0]*v[0])-18*s[0]*s[2]-6*s[0]*v[2]-6*s[2]*v[0]
 fw={**weights,wb2:0};weight_audit('fixed_gamma2_Weyl',J1+ks*(s[0]**2+2*s[0]*v[0]),shapes,bg+[wb2],fw,4,W)
 weight_audit('fixed_gamma4_Weyl',3*s[0]**2+2*s[0]*v[0]-2*W*v[0]**2,shapes,bg,weights,4,W)
 # Independently retained normal z-polynomials, with t=1. Convert radial
 # jets at Y=1: h_Y=-h_z, h_YY=h_zz+h_z.
 z=S.symbols('z');p=S.symbols('p:3');w=S.symbols('w:3');ep=S.symbols('ep')
 un=1+ep*(p[1]+p[2]);vn=(1+ep*p[1])**2;bn=ep*(w[1]+w[2])
 normalL=15*S.exp(-2*ep*w[0])*((z+ep*p[0])*bn*bn/15-un*bn/15-un*un/10-bn*vn/15-R(7,90)*un*vn-R(7,180)*vn*vn)
 an=ep*w[1];axn=ep*(w[1]+w[2])-ep**2*w[1]**2;sn=1+ep*p[1];sxn=(1-ep*w[1])*(1+ep*p[1])+ep*p[2]
 Gn=(z+ep*p[0])*(12*an*axn+R(44,7)*an**3)-6*axn*sn-R(24,7)*an*an*sn-6*an*sxn-R(96,7)*an*sn*sn-18*sn*sxn-R(85,21)*sn**3
 normalB=-S.exp(-2*ep*w[0])*Gn/6
 nc=lambda f:S.expand(S.diff(f,ep,2).subs(ep,0)/2)
 lb=nc(normalL);bt=nc(normalB)
 rep={b:0,b1:0,W:z,W1:-1,W2:1,ks:0,**{s[j]:[p[0],-p[1],p[2]+p[1]][j] for j in range(3)},**{v[j]:[w[0],-w[1],w[2]+w[1]][j] for j in range(3)}}
 zero('independent_complete_normal_bulk',qb.subs(rep,simultaneous=True)+lb)
 zero('independent_complete_normal_boundary',qt.subs(rep,simultaneous=True)+bt)
 zero('normal_bulk_log',S.diff(lb,z)-(w[1]+w[2])**2)
 zero('normal_boundary_log',S.diff(bt,z)+2*w[1]*(w[1]+w[2]))
 rb=S.expand(lb-z*(w[1]+w[2])**2);rt=S.expand(bt+2*z*w[1]*(w[1]+w[2]))
 check('normal_remainders_no_log',not rb.has(z) and not rt.has(z))
 args=[*s,*v,*bg]
 raw=S.lambdify([*args,eta],[rawbulk,rawboundary],'mpmath');coeff=S.lambdify(args,[qb,qt],'mpmath')
 return raw,coeff,rb,rt,[*p,*w]

def density_and_bound_identities():
 Y=S.symbols('Y',positive=True);m=S.Function('m')(Y);v=S.Function('v')(Y);om=S.symbols('om');freq=S.symbols('freq');a=S.Function('a')(Y)
 u=m*v;W=S.diff(m,Y,2)/m+freq**2/a**2
 zero('exact_density_Laplacian',m*(S.diff(v,Y,2)+2*S.diff(m,Y)*S.diff(v,Y)/m-freq**2*v/a**2)-(S.diff(u,Y,2)-W*u))
 A,B,omega,logY=S.symbols('A B omega logY')
 zero('bulk_compensator_difference',-omega*(A-B)**2-logY*A*A-(-(omega+logY)*A*A+2*omega*A*B-omega*B*B))
 zero('boundary_density_product',S.diff(u,Y)*S.diff(u,Y,2)-m*m*S.diff(v,Y)*S.diff(v,Y,2)
      -(m*S.diff(m,Y)*v*S.diff(v,Y,2)+2*m*S.diff(m,Y)*S.diff(v,Y)**2+(m*S.diff(m,Y,2)+2*S.diff(m,Y)**2)*v*S.diff(v,Y)+S.diff(m,Y)*S.diff(m,Y,2)*v*v))
 t=S.symbols('t',positive=True);D=168*(5*t+3)*(11*t+3);N=78815*t**3+74391*t**2+16137*t+81
 zero('strict_nonzero_trace_margin',N/D-R(3,56)-t*(78815*t*t+73896*t+15705)/D)
 check('trace_margin_coefficients_positive',all(c>0 for c in S.Poly(78815*t*t+73896*t+15705,t).all_coeffs()))
 i,j=S.symbols('i j',integer=True)
 zero('critical_bulk_log_coordinate_power',(i+j-4)+(R(3,2)-i)+(R(3,2)-j)+1)
 zero('critical_boundary_log_coordinate_power',(i+j-3)+(R(3,2)-i)+(R(3,2)-j))

def contour_fixtures(raw,coeff):
 for case,row in enumerate(PROTOCOL['jet_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    args=list(map(mp.mpf,row));exact=coeff(*args);records=[]
    for n in PROTOCOL['contour_nodes']:
     radius=mp.mpf(PROTOCOL['contour_radius']);vals=[0j,0j]
     for j in range(n):
      eta=radius*mp.exp(2j*mp.pi*j/n);value=raw(*args,eta)
      for ch in range(2):vals[ch]+=2*value[ch]/(eta*eta*n)
     errors=[err(vals[ch],exact[ch]) for ch in range(2)]
     check(f'unexpanded_Weyl_Hessian_{case}_{digits}_{n}',max(errors)<mp.mpf(PROTOCOL['contour_tolerance']),out(max(errors)))
     records.append({'nodes':n,'scaled_errors':[out(e) for e in errors]})
    levels.append({'precision':digits,'hessian':[out(x) for x in exact],'contours':records})
  with mp.workdps(80):
   pe=max(err(mp.mpf(a),mp.mpf(b)) for a,b in zip(levels[0]['hessian'],levels[1]['hessian']));check('contour_precision_'+str(case),pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
  FIXTURES.append({'fixture':'unexpanded_local_Weyl_Hessian','jets':row,'levels':levels})

def quadratic_evaluator(poly,variables):
 terms=[]
 for powers,c in S.Poly(poly,*variables).terms():
  indices=[i for i,n in enumerate(powers) for _ in range(n)];assert len(indices)==2
  terms.append((indices[0],indices[1],int(S.numer(c)),int(S.denom(c))))
 def evaluate(values):return sum(mp.mpf(n)/d*values[i]*values[j] for i,j,n,d in terms)
 def integral(modes,left,right=None):
  total=0j
  for i,j,n,d in terms:
   for ci,li in modes[i]:
    for cj,lj in modes[j]:
     lam=li+lj;fac=(-mp.exp(lam*left) if right is None else mp.exp(lam*right)-mp.exp(lam*left))/lam
     total+=mp.mpf(n)/d*ci*cj*fac
  return mp.re(total)
 def envelope(modes,z,integrated):
  total=mp.mpf(0)
  for i,j,n,d in terms:
   for ci,li in modes[i]:
    for cj,lj in modes[j]:
     lam=li+lj;total+=abs(mp.mpf(n)/d*ci*cj)*mp.exp(mp.re(lam)*z)/(abs(lam) if integrated else 1)
  return total
 return evaluate,integral,envelope

def local_limit_fixtures(rb,rt,variables):
 bulk,integral,bulkbound=quadratic_evaluator(rb,variables);boundary,_,boundarybound=quadratic_evaluator(rt,variables)
 for case,row in enumerate(PROTOCOL['profiles']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    al,be,nu,mu=map(mp.mpf,row);modes=[]
    for rate,freq,amplitude in [(al,nu,mp.mpf('.3')),(be,mu,mp.mpc('.2','.1'))]:
     lam=-rate+1j*freq
     for j in range(3):
      c=amplitude*(lam-mp.mpf('1.5'))**j/2;modes.append([(c,lam),(mp.conj(c),mp.conj(lam))])
    def values(z):return [mp.re(sum(c*mp.exp(lam*z) for c,lam in mode)) for mode in modes]
    limit=integral(modes,0);records=[]
    for cut in PROTOCOL['cut_Z']:
     z=mp.mpf(cut);direct=mp.quad(lambda x:bulk(values(x)),list(range(cut+1)));exact=integral(modes,0,z)
     er=err(direct,exact);check(f'local_bulk_integral_{case}_{digits}_{cut}',er<mp.mpf(PROTOCOL['identity_tolerance']),out(er))
     value=direct+boundary(values(z));tail=-integral(modes,z)+boundary(values(z))
     er2=err(value-limit,tail);check(f'local_remainder_tail_{case}_{digits}_{cut}',er2<mp.mpf(PROTOCOL['identity_tolerance']),out(er2))
     enclosure=bulkbound(modes,z,True)+boundarybound(modes,z,False)
     check(f'local_tail_enclosure_{case}_{digits}_{cut}',abs(value-limit)<=enclosure,out(abs(value-limit)/enclosure))
     records.append({'Z':cut,'remainder_action':out(value),'limit_error':out(abs(value-limit)),'limit_error_enclosure':out(enclosure),'identity_error':out(er)})
    ratio=mp.mpf(records[-1]['limit_error_enclosure'])/mp.mpf(records[0]['limit_error_enclosure']);check('local_cut_refinement_'+str(case)+'_'+str(digits),ratio<mp.mpf('.1'),out(ratio))
    levels.append({'precision':digits,'limit':out(limit),'values':records,'refinement_ratio':out(ratio)})
  with mp.workdps(80):
   pe=err(mp.mpf(levels[0]['limit']),mp.mpf(levels[1]['limit']));check('local_limit_precision_'+str(case),pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
  FIXTURES.append({'fixture':'complete_normal_Weyl_remainder_cutoff','profile':row,'levels':levels,'actual_shape_solution':False})

def weak_boundary_fixtures():
 x=S.symbols('x');psi=x*x*(1-x*x)**6/2
 zero('spike_second_jet',S.diff(psi,x,2).subs(x,0)-1)
 zero('spike_H2_scaling',2*13-3*8-2)
 fun=S.lambdify(x,[S.diff(psi,x),S.diff(psi,x,2)],'mpmath');levels=[]
 for digits in PROTOCOL['precision_digits']:
  with mp.workdps(digits):
   records=[]
   for n in PROTOCOL['spike_powers']:
    d=mp.mpf(2)**(-n);ell=d**8;amp=d**13
    def jets(x):
     y=d+ell*x;p1,p2=fun(x);return y,2*y+amp/ell*p1,2+amp/ell**2*p2
    # Subdivide the polynomial support. The Cauchy--Schwarz bound below
    # controls the absolute logarithmic pairing without sampling its sign.
    intervals=[mp.mpf(-1)+mp.mpf(j)/32 for j in range(65)]
    first=ell*mp.quad(lambda x:(jets(x)[1]/jets(x)[0])**2,intervals)
    second=ell*mp.quad(lambda x:jets(x)[2]**2,intervals)
    signed=ell*mp.quad(lambda x:jets(x)[1]*jets(x)[2]/jets(x)[0],intervals)
    bound=mp.sqrt(first*second);peak=2*d*(2+d**(-3))
    check(f'critical_boundary_pairing_bound_{digits}_{n}',abs(signed)<=bound,out(abs(signed)/bound))
    check(f'critical_boundary_small_log_mass_{digits}_{n}',bound<100*d**5,out(bound/d**5))
    records.append({'power':n,'pointwise_peak':out(peak),'log_L1_upper_bound':out(bound),'signed_pairing':out(signed)})
   check('boundary_peak_growth_'+str(digits),mp.mpf(records[-1]['pointwise_peak'])>100*mp.mpf(records[0]['pointwise_peak']))
   levels.append({'precision':digits,'values':records})
 with mp.workdps(80):
  pe=max(err(mp.mpf(a['log_L1_upper_bound']),mp.mpf(b['log_L1_upper_bound'])) for a,b in zip(levels[0]['values'],levels[1]['values']));check('weak_boundary_precision',pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
 FIXTURES.append({'fixture':'concentrated_critical_boundary_contact','levels':levels,'ordinary_pointwise_limit_claimed':False})

def obstruction_fixtures():
 for case,ts in enumerate(PROTOCOL['trace_t']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    t=mp.mpf(ts);k=mp.sqrt(1/t+mp.mpf(17)/12)
    rc=-(78815*t**3+74391*t*t+16137*t+81)/(168*(5*t+3)*(11*t+3))
    rs=-k*(6655*t**3+111*t*t-3015*t-567)/(84*(5*t+3)*(11*t+3));amplitude=mp.sqrt(rc*rc+rs*rs)
    check(f'trace_amplitude_margin_{case}_{digits}',amplitude>mp.mpf(3)/56,out(amplitude))
    bump=lambda x:mp.exp(-1/(1-x*x)) if abs(x)<1 else mp.mpf(0)
    norm=mp.quad(bump,[-1,0,1]);response=mp.quad(lambda x:bump(x)*mp.cos(mp.pi*x/4),[-1,0,1])/norm
    check(f'compact_test_response_{case}_{digits}',mp.sqrt(mp.mpf('.5'))<response<1,out(response))
    phase=mp.arg(mp.mpc(rc,rs));h=mp.pi/(8*k);vals=[]
    for opposite in (0,1):
     center=(phase+20*mp.pi+opposite*mp.pi)/(2*k)
     val=mp.quad(lambda x:bump(x)*(rc*mp.cos(2*k*(center+h*x))+rs*mp.sin(2*k*(center+h*x))),[-1,0,1])/norm
     target=(-1)**opposite*amplitude*response;er=err(val,target)
     check(f'translated_obstruction_test_{case}_{digits}_{opposite}',er<mp.mpf(PROTOCOL['identity_tolerance']),out(er));vals.append(out(val))
    levels.append({'precision':digits,'amplitude':out(amplitude),'test_response':out(response),'opposite_phase_tests':vals})
  with mp.workdps(80):
   pe=err(mp.mpf(levels[0]['amplitude']),mp.mpf(levels[1]['amplitude']));check('obstruction_precision_'+str(case),pe<mp.mpf(PROTOCOL['precision_tolerance']),out(pe))
  FIXTURES.append({'fixture':'nonzero_translated_trace_obstruction','t':ts,'levels':levels,'physical_instability_claimed':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 combined shape cutoff obstruction','input_count':len(ps)}),flush=True)
 try:
  provenance();raw,coeff,rb,rt,variables=derive();density_and_bound_identities();contour_fixtures(raw,coeff);local_limit_fixtures(rb,rt,variables);weak_boundary_fixtures();obstruction_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 combined shape cutoff obstruction','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed','exception':error,
  'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'monomial_weight_audits':AUDITS,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'full_pure_shape_weak_limit_derived':good,'original_local_action_retained':True,'conditional_fast_expansion_obstructed':good,
  'physical_fast_branch_excluded':False,'all_nonlinear_families_excluded':False,'actual_charged_two_point_kernel_evaluated':False,
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
