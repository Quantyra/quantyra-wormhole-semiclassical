"""Full energy-class shape forcing in the nonlinear trace; pure shape action open."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];STEM='energy-shape-nonlinear-trace-forcing'
OWN=['checks/energy_shape_nonlinear_trace_forcing.py','notes/'+STEM+'.md','evidence/energy-shape-trace-source-review.json']
PREVIOUS='first-order-mixed-boundary-response';ROWS=[];FIXTURES=[];R=S.Rational;POLYNOMIALS=[]
PROTOCOL={'precision_digits':[50,80],'charged_radii':['80','160','320'],
 'charged_cases':[['2','.01','1.1','2','2'],['3','-.02','1.2','3','5'],['4','1e-14','1.05','5','9']],
 'charged_columns':['lambda','d','horizon','magnetic_flux','k'],
 'contour_nodes':[16,32],'contour_radius':'.03','contour_scaled_limit':'1e-30',
 'charged_coarse_limit':'.01','charged_fine_limit':'.001','charged_r2_error_bound':'100',
 'precision_scaled_limit':'1e-35','quadrature_scaled_limit':'1e-35',
 'shape_A':['.3','.25','2'],'shape_V':['.2','.2','3'],
 'shape_profile_columns':['amplitude','extra_decay','log_frequency'],
 'spike_indices':[4,8,12],'green_points':['-.4','0','.4','2','4','8'],
 'static_numerical_fixtures':True,'thermal_transfer_is_ordinary_coefficient_bound':True,
 'actual_shape_solution_evaluated':False,'canonical_determinant_evaluated':False,
 'physical_benchmark_evaluated':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.factor(S.cancel(S.sympify(value).doit()));check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==123 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('scope',not rev['physical_fast_branch_excluded'] and not rev['pure_shape_finite_cut_limit_proved'] and not rev['new_boundary_prescription_adopted'] and not rev['physical_stability_proved'])

def quadratic_derivation():
 t,z,eta=S.symbols('t z eta');w=S.symbols('w:5');v=S.symbols('v:5');P=S.symbols('P:5');V=S.symbols('V:5')
 u=1+w[1]+w[2];q=(1+w[1])**2;b=v[1]+v[2]
 LC=(S.exp(2*w[0])*(6*q+4*(1+w[1])*v[1])+6*(1+R(11,6)*t)*S.exp(4*w[0]+2*v[0]))/4
 LW=15*t*S.exp(-2*v[0])*((z+w[0])*b*b/15-u*b/15-u*u/10-b*q/15-7*u*q/90-7*q*q/180)
 sub={**{w[i]:eta*w[i] for i in range(5)},**{v[i]:eta*v[i] for i in range(5)}}
 jets={**{w[i]:sum(S.binomial(i,j)*(-R(3,2))**(i-j)*P[j] for j in range(i+1)) for i in range(5)},**{v[i]:sum(S.binomial(i,j)*(-R(3,2))**(i-j)*V[j] for j in range(i+1)) for i in range(5)}}
 F=[S.expand(S.diff(S.diff(LC+LW,w[i]).subs(sub,simultaneous=True),eta,2).subs(eta,0)/2).subs(jets,simultaneous=True).expand() for i in range(3)]
 def D(q):return S.diff(q,z)+sum(S.diff(q,w[i])*w[i+1]+S.diff(q,v[i])*v[i+1] for i in range(4))
 n=S.exp(-w[0]-v[0]);K=[n*(1+w[1]+v[1]),n*(1+w[1]),n*(1+w[1])]
 DK=[n*D(x) for x in K]
 Ric=[-sum(DK[i]+K[i]**2 for i in range(3))]+[-DK[i]-K[i]*sum(K) for i in range(3)]
 scalar=S.simplify(sum(Ric));ric2=S.expand(sum(x*x for x in Ric))
 riem2=4*sum((DK[i]+K[i]**2)**2 for i in range(3))+4*sum((K[i]*K[j])**2 for i in range(3) for j in range(i+1,3))
 weyl=S.simplify(riem2-2*ric2+scalar*scalar/3);jsq=S.simplify(ric2-scalar*scalar/4)
 box=S.exp(-2*w[0]-2*v[0])*(D(D(scalar))+(3+2*w[1])*D(scalar))
 trace=scalar+12+22*t+t*(R(7,12)*weyl+R(11,6)*jsq-R(11,72)*scalar*scalar+box)
 raw=S.exp(4*w[0]+2*v[0])*trace/2
 second=S.expand(S.diff(raw.subs(sub,simultaneous=True),eta,2).subs(eta,0)/2).subs(jets,simultaneous=True).expand()
 def DP(ff):return S.diff(ff,z)+sum(S.diff(ff,a[i])*a[i+1] for a in [P,V] for i in range(4))
 force=S.expand(F[0]-DP(F[1])+DP(DP(F[2])))

 zero('independent_full_coupled_quadratic_trace',second-force)
 for i,ff in enumerate(F):
  zero('no_explicit_radial_log_'+str(i),S.diff(ff,z))
  check('quadratic_second_jet_polynomial_'+str(i),S.Poly(ff,*P,*V).total_degree()==2 and not any(ff.has(j) for j in [*P[3:],*V[3:]]))
  POLYNOMIALS.append(str(ff))
 f,g=S.symbols('f g');k2=1/t+R(17,12)
 fast={**{P[i]:[f,g,-k2*f,-k2*g,k2*k2*f][i] for i in range(5)},**{x:0 for x in V}}
 Nt=(R(471,4)*t+R(289,4)+6/t)*f*f+(R(35,3)*t-5)*f*g-(21*t+3)*g*g
 zero('unchanged_persistent_fast_forcing',force.subs(fast,simultaneous=True)-Nt)
 target=t*(-21*P[0]**2/8+7*P[0]*P[1]/2-5*P[0]*V[0]/2-7*P[1]**2/6-22*P[1]*V[0]/3+6*P[2]*V[0]-41*V[0]**2/6-4*V[0]*V[1]+2*V[0]*V[2])
 zero('explicit_second_flux_polynomial',F[2]-target)
 zero('L1_second_shape_jet_square',S.Poly(F[0],*P,*V).coeff_monomial(V[2]**2)-t)
 A=S.symbols('A:5');shape={P[i]:P[i]+A[i] for i in range(5)}
 for i,ff in enumerate(F):
  error=S.Poly(S.expand(ff.subs(shape,simultaneous=True)-ff.subs({x:0 for x in V})),*P,*A,*V)
  check('error_contains_shape_in_every_monomial_'+str(i),all(sum(m[5:])>=1 for m,c in error.terms()))
 check('nonzero_mixed_forcing_retained',any(S.diff(force,P[i],V[j])!=0 for i in range(5) for j in range(5)))
 # Normal optical/logarithmic norm identities.
 y=S.symbols('y',positive=True);H=S.Function('H')(-S.log(y));h=y**R(3,2)*H
 zero('weighted_H1_log_identity',S.diff(h/y,y)*S.sqrt(y)-(H/2-S.Subs(S.Derivative(S.Function('H')(z),z),z,-S.log(y))))
 zero('H2_log_identity',S.diff(h,y,2)*S.sqrt(y)-(S.Subs(S.Derivative(S.Function('H')(z),z,2),z,-S.log(y))-2*S.Subs(S.Derivative(S.Function('H')(z),z),z,-S.log(y))+3*H/4))
 args=[t,*P,*V]
 return S.lambdify(args,F,'mpmath'),S.lambdify(args,force,'mpmath')


def green_derivation():
 x,t,k=S.symbols('x t k',real=True);c=1/(3+23*t)
 plus=c*S.exp(-x);minus=c*S.exp(3*x/2)*(S.cos(k*x)-5*S.sin(k*x)/(2*k))
 replacement={t:1/(k*k-R(17,12))}
 for label,h in [('right',plus),('left',minus)]:
  ode=3*t*S.diff(h,x,3)-6*t*S.diff(h,x,2)+(2*t+3)*S.diff(h,x)+(11*t+3)*h
  zero('Green_homogeneous_'+label,S.expand(ode).subs(replacement))
 for j in [0,1]:zero('Green_continuity_'+str(j),(S.diff(plus,x,j)-S.diff(minus,x,j)).subs(x,0))
 zero('Green_second_jump',((S.diff(plus,x,2)-S.diff(minus,x,2)).subs(x,0)-1/(3*t)).subs(replacement))
 zero('response_derivative_contact',3/(3*t)-1/t)
 D=S.symbols('D');C=3*t*D*D-9*t*D+11*t+3
 zero('tail_operator_factorization',C-3*t*((D-R(3,2))**2+1/t+R(17,12)))
 hh,jj=S.symbols('hh jj');X=3*(D-4)*hh+(D-6)*jj
 zero('shared_tail_row',-(D-4)*(D+1)*C*hh-(D-6)*(D+1)*C*jj/3+(D+1)*C*X/3)
 zero('shared_tail_action',3*t*(7*D-5)*(D-4)*hh/7+t*(7*D-5)*(D-6)*jj/7-t*(7*D-5)*X/7)
 rc=-(78815*t**3+74391*t*t+16137*t+81)/(168*(5*t+3)*(11*t+3))
 zero('persistent_gap_polynomial',(-rc-R(3,56))*168*(5*t+3)*(11*t+3)-(78815*t**3+73896*t*t+15705*t))
 zero('spike_H2_energy_exponent',2*R(7,4)-3-R(1,2))
 zero('spike_quadratic_fourth_jet_exponent',2*R(7,4)-4+R(1,2))
 return S.lambdify([x,t,k],[[S.diff(plus,x,j) for j in range(4)],[S.diff(minus,x,j) for j in range(4)]],'mpmath')

def charged_trace_function():
 r,lam,t,qm=S.symbols('r lam t qm',positive=True);ff=S.symbols('F:5');w=S.symbols('w:5');v=S.symbols('v:5')
 F,Fr=ff[:2]
 def D(a):return S.diff(a,r)+sum(S.diff(a,x[i])*x[i+1] for x in [ff,w,v] for i in range(4))
 n=S.sqrt(F)*S.exp(-w[0]-v[0]);Kt=n*(w[1]+v[1]+Fr/(2*F));Ks=n*(w[1]+1/r)
 at=S.expand(-n*D(Kt)-Kt*Kt);ass=S.expand(-n*D(Ks)-Ks*Ks);bs=-Kt*Ks;cs=S.exp(-2*w[0])/r**2-Ks*Ks
 eig=[at+2*ass,at+2*bs,ass+bs+cs,ass+bs+cs]
 scalar=S.simplify(sum(eig));ric2=sum(a*a for a in eig);riem2=4*(at*at+2*ass*ass+2*bs*bs+cs*cs)
 weyl=S.simplify(riem2-2*ric2+scalar*scalar/3);jsq=ric2-scalar*scalar/4
 box=S.exp(-2*w[0]-2*v[0])*(F*D(D(scalar))+(Fr+2*F/r+2*F*w[1])*D(scalar))
 trace=scalar+12*lam*(1+11*t/6)+(t/lam)*(7*weyl/12+11*jsq/6-11*scalar*scalar/72-10*qm*qm*S.exp(-4*w[0])/r**4+box)
 vacuum={**{x:0 for x in w},**{x:0 for x in v}}
 zero('charged_seed_scalar',scalar.subs(vacuum)-(-ff[2]-4*Fr/r+2*(1-F)/r**2))
 zero('charged_seed_Weyl',weyl.subs(vacuum)-(-ff[2]+2*Fr/r+2*(1-F)/r**2)**2/3)
 raw=S.exp(4*w[0]+2*v[0])*trace/(2*lam)
 return S.lambdify([r,lam,t,qm,*ff,*w,*v],raw,'mpmath')

def fall(x,n):
 out=1
 for j in range(n):out*=x-j
 return out

def charged_fixtures(raw,force):
 for ix,row in enumerate(PROTOCOL['charged_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    lam,d,rh,qm,k=map(mp.mpf,row);t=1/(k*k-mp.mpf(17)/12)
    Q=1+3*lam+d;m=rh*(1+lam*rh**2+Q/rh**2-2*d/rh**3+d/rh**4)/2
    aa,da,oa=map(mp.mpf,PROTOCOL['shape_A']);av,dv,ov=map(mp.mpf,PROTOCOL['shape_V'])
    scale=(471*t/4+mp.mpf(289)/4+6/t)+abs(35*t/3-5)*k+(21*t+3)*k*k
    records=[]
    for ri,rstring in enumerate(PROTOCOL['charged_radii']):
     r=mp.mpf(rstring);z=mp.log(r);nu=-mp.mpf('1.5')+1j*k;nua=-mp.mpf('1.5')-da+1j*oa;nuv=-mp.mpf('1.5')-dv+1j*ov
     fj=[sum(c*fall(n,j)*r**(n-j) for c,n in [(lam,2),(1,0),(-2*m,-1),(Q,-2),(-2*d,-3),(d,-4)]) for j in range(5)]
     w=[mp.re(fall(nu,j)*r**(nu-j)+aa*fall(nua,j)*r**(nua-j)) for j in range(5)]
     v=[mp.re(av*fall(nuv,j)*r**(nuv-j)) for j in range(5)]
     P=[mp.re((1j*k)**j*mp.exp(1j*k*z)+aa*(-da+1j*oa)**j*mp.exp((-da+1j*oa)*z)) for j in range(5)]
     V=[mp.re(av*(-dv+1j*ov)**j*mp.exp((-dv+1j*ov)*z)) for j in range(5)]
     fn=lambda eta:raw(r,lam,t,qm,*fj,*[eta*x for x in w],*[eta*x for x in v])
     actual=r**3*mp.diff(fn,0,2)/2;expected=force(t,*P,*V);error=abs(actual-expected)/scale
     tag=str(ix)+'_'+str(digits)+'_'+str(ri)
     check('charged_normal_limit_'+tag,error<mp.mpf(PROTOCOL['charged_coarse_limit']),mp.nstr(error,25))
     contours=[]
     for nodes in PROTOCOL['contour_nodes']:
      rad=mp.mpf(PROTOCOL['contour_radius']);coef=sum(fn(rad*mp.exp(2j*mp.pi*j/nodes))*mp.exp(-4j*mp.pi*j/nodes)/(nodes*rad*rad) for j in range(nodes))
      ce=abs(r**3*coef-actual)/(1+abs(actual))
      check('charged_unexpanded_'+tag+'_'+str(nodes),ce<mp.mpf(PROTOCOL['contour_scaled_limit']),mp.nstr(ce,25))
      contours.append({'nodes':nodes,'scaled_error':mp.nstr(ce,25)})
     records.append({'radius':rstring,'actual':mp.nstr(actual,digits),'normal':mp.nstr(expected,digits),'scaled_error':mp.nstr(error,25),'contours':contours})
    check('charged_fine_'+str(ix)+'_'+str(digits),mp.mpf(records[-1]['scaled_error'])<mp.mpf(PROTOCOL['charged_fine_limit']))
    for j,record in enumerate(records):
     weighted=mp.mpf(record['scaled_error'])*mp.mpf(record['radius'])**2
     check('charged_r2_error_'+str(ix)+'_'+str(digits)+'_'+str(j),weighted<mp.mpf(PROTOCOL['charged_r2_error_bound']),mp.nstr(weighted,25))
    levels.append({'precision':digits,'values':records})
  with mp.workdps(80):
   pe=max(abs(mp.mpf(a['actual'])-mp.mpf(b['actual']))/(1+abs(mp.mpf(b['actual']))) for a,b in zip(levels[0]['values'],levels[1]['values']))
   check('charged_precision_'+str(ix),pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,25))
  FIXTURES.append({'fixture':'charged_full_coupled_quadratic_trace','case':row,'levels':levels,'physical_solution':False})

def bump_functions():
 x=S.symbols('x');b=(1-x*x)**6
 return S.lambdify(x,[S.diff(b,x,j) for j in range(5)],'mpmath')

def distribution_fixtures(Fs,force,bump):
 for ix,row in enumerate(PROTOCOL['charged_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    k=mp.mpf(row[-1]);t=1/(k*k-mp.mpf(17)/12);records=[]
    for n in PROTOCOL['spike_indices']:
     ell=mp.mpf(2)**(-n);amp=ell**mp.mpf('1.75')
     def val(x):
      z=ell*x;P=[mp.re((1j*k)**j*mp.exp(1j*k*z)) for j in range(5)];V=[amp*jj/ell**j for j,jj in enumerate(bump(x))]
      ff=Fs(t,*P,*V);base=Fs(t,*P,*[0]*5);err=[ff[j]-base[j] for j in range(3)]
      test=mp.exp(-z*z);testj=[test,-2*z*test,(4*z*z-2)*test]
      weak=sum(err[j]*testj[j] for j in range(3));direct=(force(t,*P,*V)-force(t,*P,*[0]*5))*test
      return weak,direct,force(t,*P,*V)-force(t,*P,*[0]*5)
     parts=[-1,-mp.mpf('.5'),0,mp.mpf('.5'),1]
     weak=ell*mp.quad(lambda x:val(x)[0],parts);direct=ell*mp.quad(lambda x:val(x)[1],parts)
     err=abs(weak-direct)/(1+abs(weak));center=val(0)[2]
     energy=amp*amp/ell**3*mp.quad(lambda x:bump(x)[2]**2,parts)
     check('weak_pairing_identity_'+str(ix)+'_'+str(digits)+'_'+str(n),err<mp.mpf(PROTOCOL['quadrature_scaled_limit']),mp.nstr(err,25))
     records.append({'index':n,'weak_forcing':mp.nstr(weak,digits),'direct_pairing':mp.nstr(direct,digits),'scaled_pairing_error':mp.nstr(err,25),'pointwise_center':mp.nstr(center,digits),'second_derivative_energy':mp.nstr(energy,digits)})
    for j in range(2):
     ratio=abs(mp.mpf(records[j+1]['weak_forcing'])/mp.mpf(records[j]['weak_forcing']))
     check('weak_forcing_decay_'+str(ix)+'_'+str(digits)+'_'+str(j),ratio<mp.mpf('.7'),mp.nstr(ratio,25))
    growth=abs(mp.mpf(records[-1]['pointwise_center'])/mp.mpf(records[0]['pointwise_center']))
    check('pointwise_forcing_not_used_'+str(ix)+'_'+str(digits),growth>2,mp.nstr(growth,25))
    levels.append({'precision':digits,'values':records})
  with mp.workdps(80):
   pe=max(abs(mp.mpf(a['weak_forcing'])-mp.mpf(b['weak_forcing']))/(1+abs(mp.mpf(b['weak_forcing']))) for a,b in zip(levels[0]['values'],levels[1]['values']))
   check('weak_pairing_precision_'+str(ix),pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,25))
  FIXTURES.append({'fixture':'concentrated_H2_shape_trace','k':row[-1],'levels':levels,'actual_shape_solution':False})

def green_fixtures(kernels,bump):
 for ix,row in enumerate(PROTOCOL['charged_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    k=mp.mpf(row[-1]);t=1/(k*k-mp.mpf(17)/12);records=[]
    def ker(x):return kernels(x,t,k)[0 if x>=0 else 1]
    for xs in PROTOCOL['green_points']:
     x=mp.mpf(xs);splits=sorted(set([mp.mpf(-1),mp.mpf('-.5'),mp.mpf(0),mp.mpf('.5'),mp.mpf(1)]+([x] if -1<x<1 else [])))
     def values(s):
      b=bump(s);g=ker(x-s);rr=[mp.mpf('.3')*b[0],mp.mpf('.2')*b[1],b[0]];src=mp.mpf('.3')*b[0]+mp.mpf('.8')*b[2]
      return 3*g[0]*src,3*(g[0]*rr[0]-g[1]*rr[1]+g[2]*rr[2]),3*g[1]*src,3*(g[1]*rr[0]-g[2]*rr[1]+g[3]*rr[2])
     vals=[mp.quad(lambda s:values(s)[j],splits) for j in range(4)];contact=bump(x)[0]/t if abs(x)<1 else 0
     err=max(abs(vals[0]-vals[1])/(1+abs(vals[0])),abs(vals[2]-vals[3]-contact)/(1+abs(vals[2])))
     check('Green_convolution_'+str(ix)+'_'+str(digits)+'_'+xs,err<mp.mpf(PROTOCOL['quadrature_scaled_limit']),mp.nstr(err,25))
     if abs(x)<1:check('Green_contact_required_'+str(ix)+'_'+str(digits)+'_'+xs,abs(vals[2]-vals[3])>mp.mpf('.01'))
     if x>1:check('Green_right_stable_tail_'+str(ix)+'_'+str(digits)+'_'+xs,abs(vals[2]+vals[0])<mp.mpf('1e-35'),mp.nstr(abs(vals[2]+vals[0]),25))
     records.append({'x':xs,'response':mp.nstr(vals[0],digits),'derivative':mp.nstr(vals[2],digits),'derivative_contact':mp.nstr(contact,digits),'scaled_error':mp.nstr(err,25)})
    levels.append({'precision':digits,'values':records})
  with mp.workdps(80):
   pe=max(abs(mp.mpf(a['response'])-mp.mpf(b['response']))/(1+abs(mp.mpf(b['response']))) for a,b in zip(levels[0]['values'],levels[1]['values']))
   check('Green_precision_'+str(ix),pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,25))
  FIXTURES.append({'fixture':'tail_Green_with_distribution_contact','k':row[-1],'levels':levels,'physical_solution':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 energy-shape nonlinear trace forcing','input_count':len(ps)}),flush=True)
 try:
  provenance();Fs,force=quadratic_derivation();kernels=green_derivation();raw=charged_trace_function();bump=bump_functions()
  charged_fixtures(raw,force);distribution_fixtures(Fs,force,bump);green_fixtures(kernels,bump)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 energy-shape nonlinear trace forcing','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'complete_coupled_trace_divergence_form_derived':good,'energy_shape_forcing_class_derived':good,
  'bounded_tail_persistent_combination_preserved':good,'Green_derivative_contact_verified':good,
  'quadratic_trace_flux_polynomials':POLYNOMIALS,
  'pure_shape_finite_cut_limit_proved':False,'full_boundary_action_solved':False,
  'pointwise_forcing_decay_assumed':False,'physical_action_averaging_adopted':False,
  'physical_fast_branch_excluded':False,'full_first_order_shape_response_included':False,'new_boundary_prescription_adopted':False,
  'rough_finite_amplitude_action_differentiability_proved':False,'full_shape_equation_solved':False,
  'static_source_compatibility_proved':False,'regulator_independence_for_moving_cuts_proved':False,
  'continuous_causal_stability_proved':False,'physical_constraint_mode_constructed':False,'actual_wormhole_instability_proved':False,
  'finite_counterterm_changed':False,'physical_curvature_norm_changed':False,'order_reduction_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in [('checker',OWN[0]),('note',OWN[1]),('review',OWN[2])]:run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
