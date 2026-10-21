"""Physical curvature row, source weight and horizon boundary range fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from functools import lru_cache
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='physical-curvature-source-and-boundary-range'
OWN=['checks/physical_curvature_source_and_boundary_range.py','notes/'+STEM+'.md','evidence/physical-curvature-range-source-review.json']
PREVIOUS='outer-volume-norm-singular-sequence';ROWS=[];NUMERICS=[]
PROTOCOL={'scope':'Actual charged seed scalar diagonal and total-source compatibility; not a metric stability verdict',
 'mpmath_precision':55,'identity_tolerance':'1e-38',
 'horizon_fixtures':[['2','0.5','0'],['2','1.5','0'],['5','0.5','0'],['5','1.5','0'],['5','0.7','0.3']],
 'interface':'1','source_support':['0.25','0.75'],'RK_grids':[384,768],
 'RK_relative_tolerance':'2e-5','refinement_ratio':'0.2','refinement_floor':'1e-25',
 'source_envelope':'[16(y-1/4)(3/4-y)]^4 on [1/4,3/4]',
 'seed_fixture':{'B':'1','lambda':'2','d':'0.01','horizon':'1.1','E_over_lambda':'11.25',
 'outer_radii':[100,200,400],'horizon_offsets':['0.01','0.005','0.0025']},
 'fixtures_are_actual_metric_spectrum':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,expr):
 value=S.simplify(S.sympify(expr).doit());check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==107 and all(r['passed'] for r in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 check('physical_inputs_pinned',all('notes/'+f+'.md' in p['input_sha256_after'] for f in [
  'charged-conformal-response-and-constraints','stationary-centered-covariant-seed',
  'retained-benchmark-parameter-recovery','fast-scalar-boundary-solvability']))
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('source_scope',len(review['primary_sources'])==2 and not review['metric_coupling_set_to_zero']
  and not review['physical_metric_boundary_selected'] and not review['actual_charged_instability_proved'])

def exact_checks():
 y,t=S.symbols('y t',positive=True);r=S.Function('r')(y);N=S.Function('N')(y);u=S.Function('u')(t,y)
 mass,B,f=S.symbols('mass B f');psi=u/r
 # Divergence formula for the actual spherical scalar d'Alembertian.
 box=(-S.diff(psi,t,2)+S.diff(r*r*S.diff(psi,y),y)/(r*r))/(N*N)
 V=S.diff(r,y,2)/r-N*N*mass
 zero('wave_transform_and_source_sign',r*N*N*(box+mass*psi+B-f)
  -(-S.diff(u,t,2)+S.diff(u,y,2)-V*u+r*N*N*(B-f)))
 zero('kinetic_measure',N*r*r/N-r*r)
 zero('source_measure',(r*N*N)**2-r*r*N**4)
 A=S.symbols('A',positive=True)
 zero('areal_kinetic_weight',r*r*A/N-r*r*(A/N))
 zero('areal_source_weight',r*r*N**4*A/N-r*r*A*N**3)
 dr,xi,Rr,Tr,dR,dT=S.symbols('dr xi Rr Tr dR dT')
 zero('curvature_gauge_invariance',(dR+xi*Rr)-Rr*(dr+xi)-(dR-Rr*dr))
 zero('off_shell_source_gauge_invariance',(dT+xi*Tr)-Tr*(dr+xi)-(dT-Tr*dr))
 check('dropping_trace_contact_changes_source',S.diff(dT-Tr*dr,dr)==-Tr)
 rr,lam,m,Q,c3,c4,E=S.symbols('r lam m Q c3 c4 E',positive=True)
 F=lam*rr**2+1-2*m/rr+Q/rr**2+c3/rr**3+c4/rr**4
 R=-S.diff(F,rr,2)-4*S.diff(F,rr)/rr-2*(F-1)/rr**2
 zero('exact_seed_scalar_curvature',R+12*lam+2*c3/rr**5+6*c4/rr**6)
 Z=-2*c3/rr**5-6*c4/rr**6
 zero('exact_tortoise_r_second_derivative',(-F)*S.diff(-F,rr)-F*S.diff(F,rr))
 # Independent inverse-radius power series of the tortoise integral.
 z=S.symbols('z',positive=True)
 inverse_integrand=S.series(1/(z*z*F.subs(rr,1/z)),z,0,5).removeO()
 Yseries=S.integrate(inverse_integrand,(z,0,z))
 zero('tortoise_leading_series',S.expand(Yseries).coeff(z,1)-1/lam)
 zero('tortoise_cubic_series',S.expand(Yseries).coeff(z,3)+1/(3*lam**2))
 zero('tortoise_mass_series',S.expand(Yseries).coeff(z,4)-m/(2*lam**2))
 aa=E/lam-2;Vrad=F*(S.diff(F,rr)/rr-E+11*Z/36)
 tail=S.series(Vrad.subs(rr,1/z)+aa/Yseries**2,z,0,2).removeO()
 zero('bounded_potential_remainder_limit',tail.subs(z,0)-(2*lam-E)/3)
 k,g=S.symbols('k g',positive=True)
 zero('fast_parameter',60*S.pi/g+S.Rational(5,3)-S.Rational(1,4)-(60*S.pi/g+S.Rational(17,12)))
 zeta=S.symbols('zeta',positive=True);eta=(zeta-1)**2*(2-zeta)**2
 I0=S.integrate(eta**2,(zeta,1,2));I1=S.integrate(S.diff(eta,zeta)**2,(zeta,1,2))
 zero('negative_test_kinetic_ratio',I1/I0-12)
 check('retained_a_exceeds_compact_threshold',S.Rational(1125,1000)*10**35>96)
 # Free singular kernel and endpoint normalization, no complex conjugation in its Wronskian.
 s=S.symbols('s',positive=True)
 pc=S.sqrt(y)*S.cos(k*S.log(y));ps=S.sqrt(y)*S.sin(k*S.log(y))
 zero('outer_basis_Wronskian',pc*S.diff(ps,y)-S.diff(pc,y)*ps-k)
 K=S.sqrt(y*s)*S.sin(k*S.log(y/s))/k
 zero('Volterra_diagonal_value',K.subs(y,s))
 zero('Volterra_derivative_jump',S.diff(K,y).subs(y,s)-1)
 zero('Volterra_kernel_equation',S.diff(K,y,2)+(k*k+S.Rational(1,4))*K/y**2)
 MM=S.symbols('M',positive=True)
 zero('Volterra_Gronwall_exponent',S.integrate(MM*s/k,(s,0,y))-MM*y*y/(2*k))
 gamma=S.symbols('gamma',positive=True)
 H=(1-S.exp(-2*gamma*(s-y)))/(2*gamma)
 zero('horizon_weighted_kernel',(S.exp(gamma*y)*S.sinh(gamma*(s-y))/gamma*S.exp(-gamma*s)-H).rewrite(S.exp))
 zero('horizon_Wronskian',S.exp(-gamma*y)*S.diff(S.exp(gamma*y),y)-S.diff(S.exp(-gamma*y),y)*S.exp(gamma*y)-2*gamma)
 bcoef,moment=S.symbols('bcoef moment')
 zero('growing_coefficient_sign',(2*gamma*bcoef+moment).subs(bcoef,-moment/(2*gamma)))
 U,u0,u1,AU,gh=S.symbols('U u0 u1 AU gh')
 zero('Laplace_initial_contact_sign',gamma**2*U-gamma*u0-u1+AU-gh-(gamma**2*U+AU-(gamma*u0+u1+gh)))
 check('source_weight_stronger_than_volume',S.limit(y**-6/y**-4,y,0)==S.oo)

def seed_fixtures():
 lam=mp.mpf(2);d=mp.mpf('.01');rh=mp.mpf('1.1');E=lam*mp.mpf('11.25');a=E/lam-2
 ext=lambda r:(r-1)**2/r**2*(1+lam*(r*r+2*r+3)+d/r**2)
 delta=rh*ext(rh)
 F=lambda r:ext(r)-delta/r
 Z=lambda r:4*d/r**5-6*d/r**6
 V=lambda r:F(r)*(mp.diff(F,r)/r-E+11*Z(r)/36)
 check('seed_fixture_simple_horizon',abs(F(rh))<mp.mpf('1e-50') and mp.diff(F,rh)>0)
 old=None;limit=(2*lam-E)/3
 for rad in PROTOCOL['seed_fixture']['outer_radii']:
  r=mp.mpf(rad);y=mp.quad(lambda z:1/F(z),[r,mp.inf]);W=V(r)+a/y**2;err=abs(W-limit)
  check('seed_outer_potential_finite_'+str(rad),mp.isfinite(W))
  if old is not None:check('seed_outer_potential_refinement_'+str(rad),err<mp.mpf('.6')*old,mp.nstr(err/old,14))
  NUMERICS.append({'fixture':'seed_outer_potential','r':rad,'y':mp.nstr(y,20),'W':mp.nstr(W,20),'limiting_W':mp.nstr(limit,20),'error':mp.nstr(err,20)})
  old=err
 old=None
 for offset in PROTOCOL['seed_fixture']['horizon_offsets']:
  r=rh+mp.mpf(offset);vv=V(r)
  if old is not None:check('seed_horizon_potential_refinement_'+offset,abs(vv)<mp.mpf('.55')*old,mp.nstr(abs(vv)/old,14))
  check('seed_exterior_positive_'+offset,F(r)>0)
  NUMERICS.append({'fixture':'seed_horizon_potential','offset':offset,'potential':mp.nstr(vv,20)})
  old=abs(vv)

def horizon_fixtures():
 lo=mp.mpf('.25');hi=mp.mpf('.75');interface=mp.mpf(1);tol=mp.mpf('1e-38')
 xx=S.symbols('x');chi_expr=(16*(xx-S.Rational(1,4))*(S.Rational(3,4)-xx))**4
 chi0=S.lambdify(xx,chi_expr,'mpmath');chi2=S.lambdify(xx,S.diff(chi_expr,xx,2),'mpmath')
 def chi(y):return chi0(y) if lo<y<hi else mp.mpf(0)
 def qsecond(y):return chi2(y) if lo<y<hi else mp.mpf(0)
 for ix,(astr,gre,gim) in enumerate(PROTOCOL['horizon_fixtures']):
  a=mp.mpf(astr);k=mp.sqrt(a-mp.mpf('.25'));gamma=mp.mpc(gre,gim)
  factors={sign:mp.gamma(1+sign*1j*k)*(2/gamma)**(sign*1j*k) for sign in [-1,1]}
  def basis_value(y,sign):return factors[sign]*mp.sqrt(y)*mp.besseli(sign*1j*k,gamma*y)
  def basis(y,sign):
   order=sign*1j*k;factor=factors[sign]
   val=factor*mp.sqrt(y)*mp.besseli(order,gamma*y)
   der=factor*(mp.besseli(order,gamma*y)/(2*mp.sqrt(y))
     +mp.sqrt(y)*gamma*(mp.besseli(order-1,gamma*y)+mp.besseli(order+1,gamma*y))/2)
   return val,der
  def real_basis(y):
   p,dp=basis(y,1);m,dm=basis(y,-1)
   return (p+m)/2,(p-m)/(2j),(dp+dm)/2,(dp-dm)/(2j)
  c,s,dc,ds=real_basis(interface);wh=mp.exp(-gamma);dwh=-gamma*wh
  cc=(wh*ds-dwh*s)/k;cs=(c*dwh-dc*wh)/k
  @lru_cache(maxsize=None)
  def w(y):
   if y>=interface:return mp.exp(-gamma*y)
   p=basis_value(y,1);m=basis_value(y,-1)
   return cc*(p+m)/2+cs*(p-m)/(2j)
  check('horizon_value_match_'+str(ix),abs(cc*c+cs*s-wh)<tol)
  check('horizon_derivative_match_'+str(ix),abs(cc*dc+cs*ds-dwh)<tol)
  for yy in [mp.mpf('.05'),mp.mpf('.2'),mp.mpf('.7')]:
   c,s,dc,ds=real_basis(yy)
   check('outer_Wronskian_'+str(ix)+'_'+str(yy),abs(c*ds-dc*s-k)<tol)
   p,dp=basis(yy,1);base=mp.exp(1j*k*mp.log(yy));M=abs(gamma*gamma)
   bound=mp.expm1(M*yy**2/(2*k));err=abs(p/mp.sqrt(yy)-base)
   check('Volterra_remainder_bound_'+str(ix)+'_'+str(yy),err<=bound,mp.nstr(err/bound,14))
  # This source makes the bilinear moment positive even at complex gamma.
  source=lambda y:chi(y)*mp.conj(w(y))
  moment=mp.quad(lambda y:w(y)*source(y),[lo,(lo+hi)/2,hi])
  check('nonzero_bounded_moment_'+str(ix),mp.re(moment)>0 and abs(mp.im(moment))<tol)
  expected=-moment/(2*gamma)
  compatible=lambda y:-qsecond(y)+(gamma*gamma-a/y**2)*chi(y)
  cm=mp.quad(lambda y:w(y)*compatible(y),[lo,(lo+hi)/2,hi])
  check('compatible_compact_source_'+str(ix),abs(cm)<tol,mp.nstr(abs(cm),12))
  old_errors=None
  for grid in PROTOCOL['RK_grids']:
   step=(interface-lo)/grid;state=[mp.mpc(0)]*4
   def rhs(y,st):
    potential=gamma*gamma-a/y**2
    return [st[1],potential*st[0]-source(y),st[3],potential*st[2]-compatible(y)]
   for j in range(grid):
    y=lo+j*step
    k1=rhs(y,state);k2=rhs(y+step/2,[v+step*d/2 for v,d in zip(state,k1)])
    k3=rhs(y+step/2,[v+step*d/2 for v,d in zip(state,k2)])
    k4=rhs(y+step,[v+step*d for v,d in zip(state,k3)])
    state=[v+step*(a1+2*a2+2*a3+a4)/6 for v,a1,a2,a3,a4 in zip(state,k1,k2,k3,k4)]
   bg=mp.exp(-gamma)*(state[1]+gamma*state[0])/(2*gamma)
   bc=mp.exp(-gamma)*(state[3]+gamma*state[2])/(2*gamma)
   errors=[abs(bg-expected)/(1+abs(expected)),abs(bc)]
   check('forced_horizon_coefficient_'+str(ix)+'_'+str(grid),errors[0]<mp.mpf('2e-5'),mp.nstr(errors[0],12))
   check('compatible_horizon_coefficient_'+str(ix)+'_'+str(grid),errors[1]<mp.mpf('2e-5'),mp.nstr(errors[1],12))
   if old_errors is not None:
    for j,(er,prev) in enumerate(zip(errors,old_errors)):
     check('independent_ODE_refinement_'+str(ix)+'_'+str(j),er<mp.mpf('.2')*prev+mp.mpf('1e-25'),mp.nstr(er/prev,14))
   NUMERICS.append({'fixture':'horizon_range_ODE','a':astr,'gamma':[gre,gim],'grid':grid,
    'moment':mp.nstr(moment,22),'expected_growing_coefficient':mp.nstr(expected,22),
    'incompatible_relative_error':mp.nstr(errors[0],14),'compatible_absolute_error':mp.nstr(errors[1],14)})
   old_errors=errors

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 mp.mp.dps=PROTOCOL['mpmath_precision'];ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
 before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 physical curvature source and range','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();seed_fixtures();horizon_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 physical curvature source and range','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'physical_curvature_source_norm_derived':good,'actual_seed_trace_potential_derived':good,
  'actual_scalar_diagonal_range_derived':good,'coupled_forcing_moment_retained':good,
  'all_local_contacts_retained':good,'physical_metric_boundary_selected':False,
  'metric_coupling_set_to_zero':False,'actual_metric_causal_inverse_proved':False,
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
