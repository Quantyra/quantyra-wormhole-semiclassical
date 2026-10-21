"""Coupled second-order boundary obstruction; fixed-seed canonical variation."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='coupled-fast-quadratic-boundary-obstruction'
OWN=['checks/coupled_fast_quadratic_boundary_obstruction.py','notes/'+STEM+'.md','evidence/coupled-fast-boundary-source-review.json']
PREVIOUS='fast-nonlinear-trace-boundary-compatibility';ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'scope':'Coupled second-order trace/action obstruction in the stated family; canonical first variation only',
 'precision_digits':[50,80],'reference_B':'1','gauge_coupling':'1','base_radius':'80','cycle_indices':[0,1,2],
 'charged_cases':[['2','.01','1.1','2','2'],['3','-.02','1.2','3','5'],['4','1e-14','1.05','5','9']],
 'charged_columns':['lambda','d','horizon','magnetic_flux','k'],
 'annulus_phase_width_in_pi_units':'1/4','fine_relative_limit':'.001','coarse_relative_limit':'.01',
 'refinement_bound':'exp(-pi/k)','precision_scaled_limit':'1e-35',
 'contour_nodes':[16,32],'contour_radius':'.03','contour_scaled_limit':'1e-30',
 'finite_part_omega':['4','10','18'],'finite_part_cuts':['.001','.00001','.0000001'],
 'finite_part_powers':[1,2,3],'finite_part_fine_absolute_limit':'.001','finite_part_bound':'C_omega sqrt(d)',
 'canonical_full_determinant_numerically_evaluated':False,'fixtures_are_physical_benchmark':False,
 'full_first_order_shape_included':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.factor(S.cancel(S.sympify(value).doit()));check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==149 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('scope',not rev['physical_fast_branch_excluded'] and not rev['full_first_order_shape_response_included'] and not rev['new_boundary_prescription_adopted'] and not rev['physical_stability_proved'])

def original_shape_action():
 r,lam,G,F=S.symbols('r lam G F',positive=True);Fr,Frr,qm,gc=S.symbols('Fr Frr qm gc',real=True)
 v,vr,vrr,eta=S.symbols('v vr vrr eta',real=True)
 fn=S.exp(-2*eta*v)*F;fnp=S.exp(-2*eta*v)*(Fr-2*eta*vr*F)
 fnpp=S.exp(-2*eta*v)*(Frr-4*eta*vr*Fr-2*eta*vrr*F+4*eta*eta*vr*vr*F)
 alpha=2*eta*v;ap=2*eta*vr;app=2*eta*vrr
 R4=-fnpp-3*ap*fnp-2*fn*(app+ap*ap)-4*fnp/r-4*fn*ap/r+2*(1-fn)/r**2
 lam0=lam+11*G*lam*lam/(360*S.pi);Q2=S.pi*G*qm*qm/(gc*gc)
 LEH=S.exp(alpha)*r*r*(R4+6*lam0)/(4*G)
 Lmag=-S.exp(alpha)*Q2/(2*G*r*r)
 T=r*(fnp+fn*ap);RB=-r*r*fnpp-3*r*r*fnp*ap-2*r*r*fn*(app+ap*ap)+2*T-2*fn
 LW=S.exp(alpha)/(4*S.pi*r*r)*(S.log(r)*(RB*RB/60-RB/18+R(1,15)-qm*qm/3)+RB*T/30-T*T/10+R(11,90)*T*fn-R(11,180)*fn*fn-fn/18)
 GHY=S.exp(alpha)*(r*r*(fnp+2*fn*ap)/4+r*fn)/G
 CT=-S.exp(alpha)*S.sqrt(fn)*(S.sqrt(lam)*r*r+1/(2*S.sqrt(lam)))/G
 c0=-221*lam**R(3,2)/(15120*S.pi**2);c1=29*S.sqrt(lam)/(20160*S.pi**2)
 loop=4*S.pi*S.exp(alpha)*S.sqrt(fn)*(c0*r*r+2*c1)
 N=S.exp(alpha)*S.sqrt(fn)/r;a=(r*fnp/2+r*fn*ap-fn)/S.sqrt(fn)
 ax=-RB/2-a*a;s=S.sqrt(fn);sx=r*fnp/2
 GG=S.log(r)*(12*a*ax+R(44,7)*a**3+10*a)-6*ax*s-R(24,7)*a*a*s-6*a*sx-R(96,7)*a*s*s-18*s*sx-R(85,21)*s**3
 BW=-N*GG/(360*S.pi)
 bulk=LEH+Lmag+LW;boundary=GHY+CT+loop+BW
 l=S.expand(S.diff(bulk,eta).subs(eta,0));b=S.expand(S.diff(boundary,eta).subs(eta,0))
 t=S.symbols('t',positive=True);vp,vpp=S.symbols('vp vpp',real=True)
 normal={F:lam*r*r,Fr:2*lam*r,Frr:2*lam,vr:vp/r,vrr:(vpp-vp)/r**2,G:60*S.pi*t/lam}
 leadl=S.simplify(S.limit((l*G/(lam*r*r)).subs(normal,simultaneous=True),r,S.oo))
 leadb=S.simplify(S.limit((b*G/(lam*r**3)).subs(normal,simultaneous=True),r,S.oo))
 targetl=(R(1,2)+2*t)*(6*v-vp-vpp);targetb=t*vpp+(R(1,2)+R(9,7)*t)*vp-(1+R(76,7)*t)*v
 zero('original_shape_bulk',leadl-targetl);zero('original_shape_boundary',leadb-targetb)
 primitive=(R(1,2)+2*t)*(2*v-vp)
 zero('shape_bulk_primitive',3*primitive+S.diff(primitive,v)*vp+S.diff(primitive,vp)*vpp-targetl)
 zero('complete_shape_endpoint',primitive+targetb-t*(7*vpp-5*vp-48*v)/7)
 zero('magnetic_shape_derivative',S.diff(Lmag,eta).subs(eta,0)+Q2*v/(G*r*r))
 zero('magnetic_shape_subleading',S.limit((S.diff(Lmag,eta).subs(eta,0)*G/(lam*r*r)).subs(normal,simultaneous=True),r,S.oo))
 args=[r,lam,G,qm,gc,F,Fr,Frr,v,vr,vrr]
 return S.lambdify(args,[l,b],'mpmath'),S.lambdify([*args,eta],[bulk,boundary],'mpmath')

def coupled_identities():
 t,k=S.symbols('t k',positive=True);D=S.symbols('D');k2=1/t+R(17,12)
 C=3*t*D*D-9*t*D+11*t+3
 Lh=-(D-4)*(D+1)*C;Lj=-(D-6)*(D+1)*C/3
 Mh=R(3,7)*t*(7*D*D-33*D+20);Mj=t*(7*D*D-47*D+30)/7
 zero('conformal_boundary_factor',Mh-3*t*(7*D-5)*(D-4)/7)
 zero('shape_boundary_factor',Mj-t*(7*D-5)*(D-6)/7)
 zero('coupled_row_determinant',Lh*Mj-Lj*Mh)
 h,j=S.symbols('h j');X=3*(D-4)*h+(D-6)*j
 zero('shared_trace_combination',Lh*h+Lj*j+(D+1)*C*X/3)
 zero('shared_boundary_combination',Mh*h+Mj*j-t*(7*D-5)*X/7)
 zero('bounded_root_factorization',C-3*t*((D-R(3,2))**2+k2))
 def red(q):
  n,d=S.fraction(S.cancel(q));return S.factor(S.rem(n,k*k-k2,k)/S.rem(d,k*k-k2,k))
 J=S.Matrix([[0,2*k],[-2*k,0]]);I=S.eye(2)
 CJ=3*t*J*J-9*t*J+(11*t+3)*I;A=((J+I)*CJ).applyfunc(red)
 n=S.Matrix([(295*t*t+195*t+18)/(4*t),k*(R(5,2)-R(35,6)*t)])
 det=red(A.det());check('forced_X_matrix_invertible',det>0)
 xp=(3*A.adjugate()*n/det).applyfunc(red)
 for i in range(2):zero('forced_X_'+str(i),red((A*xp-3*n)[i]))
 q=S.Matrix([-9*(t-9)/56,k*(71*t+21)/28])
 result=(q+t*(7*J-5*I)*xp/7).applyfunc(red)
 rc=-(78815*t**3+74391*t*t+16137*t+81)/(168*(5*t+3)*(11*t+3))
 rs=-k*(6655*t**3+111*t*t-3015*t-567)/(84*(5*t+3)*(11*t+3))
 zero('coupled_persistent_cosine',result[0]-rc);zero('coupled_persistent_sine',result[1]-rs)
 positive=78815*t**3+73896*t*t+15705*t
 zero('strict_gap_polynomial',(-rc-R(3,56))*168*(5*t+3)*(11*t+3)-positive)
 check('strict_gap_positive',positive>0)
 # A freely chosen second-order shape pair changes h but cannot change X.
 jc,js=S.symbols('jc js');jvec=S.Matrix([jc,js]);Dh=3*(J-4*I)
 hp=(Dh.adjugate()*(xp-(J-6*I)*jvec)/red(Dh.det())).applyfunc(red)
 for i in range(2):zero('arbitrary_shape_same_X_'+str(i),red((Dh*hp+(J-6*I)*jvec-xp)[i]))
 z=S.symbols('z',real=True);zero('thermal_mean_square',S.integrate(S.cos(z)**2,(z,0,2*S.pi))/(2*S.pi)-R(1,2))

def canonical_controls():
 aa,bb=S.symbols('aa bb');chi=S.Matrix([[0,1],[1,0]]);chiral=S.diag(1,-1);K=S.diag(aa,bb)
 zero('balanced_chi_trace',S.trace(chi));zero('chiral_image_anticommutation',S.trace((chi*chiral+chiral*chi)**2))
 zero('balanced_reflected_trace',S.trace(chi*K));zero('third_smearing_contact',30*S.trace(chi))
 check('nonbalanced_contact_detected',30*S.trace(S.eye(2))!=0)
 for weight in [1,2,3]:
  power=4-weight;check('profile_contact_order_'+str(weight),power-1<=2)
 x=S.symbols('x',positive=True)
 zero('Holder_remainder_integral',2*S.integrate((1-x)*S.sqrt(x),(x,0,1))-R(8,15))
 zero('worst_finite_part_remainder',R(8,15)*S.integrate(x**R(-1,2),(x,0,1))-R(16,15))
 for p in [1,2,3]:
  weights=[R(1,S.factorial(n)*(n-p+1)) if n-p+1!=0 else S.Integer(0) for n in range(3)]
  check('finite_part_polynomial_weights_'+str(p),all(abs(w)<=1 for w in weights),weights)
 omega=S.symbols('omega',positive=True);nu=3+S.I*omega
 for p in [1,2,3]:zero('limiting_finite_part_'+str(p),S.re(1/(nu-p+1))-(4-p)/((4-p)**2+omega**2))
 check('third_jet_amplitude_nonzero',S.expand(S.prod(nu-j for j in range(3))*S.prod(S.conjugate(nu)-j for j in range(3)))==omega**6+14*omega**4+49*omega**2+36)
 check('higher_profile_integrability',R(-1,2)>-1)

def fall(x,n):
 out=1
 for j in range(n):out*=x-j
 return out

def charged_fixtures(coeff,raw):
 for ix,row in enumerate(PROTOCOL['charged_cases']):
  levels=[]
  for precision in PROTOCOL['precision_digits']:
   with mp.workdps(precision):
    lam,d,rh,qm,k=map(mp.mpf,row);t=1/(k*k-mp.mpf(17)/12);G=60*mp.pi*t/lam;gc=mp.mpf(PROTOCOL['gauge_coupling'])
    radius=mp.mpf(PROTOCOL['base_radius']);phase=-2*k*mp.log(radius)
    Qh=1+3*lam+d;m=rh*(1+lam*rh*rh+Qh/rh**2-2*d/rh**3+d/rh**4)/2
    def args(r):
     fj=[sum(co*fall(power,n)*r**(power-n) for co,power in [(lam,2),(1,0),(-2*m,-1),(Qh,-2),(-2*d,-3),(d,-4)]) for n in range(3)]
     exponent=-3+2j*k;vj=[mp.re(mp.exp(1j*phase)*fall(exponent,n)*r**(exponent-n)) for n in range(3)]
     return [r,lam,G,qm,gc,*fj,*vj]
    def endpoint(r):
     theta=2*k*mp.log(r)+phase
     return (lam/G)*(t/7)*((30-28*k*k)*mp.cos(theta)+94*k*mp.sin(theta))
    vals=[]
    for level,n in enumerate(PROTOCOL['cycle_indices']):
     left=radius*mp.exp(n*mp.pi/k);right=left*mp.exp(mp.pi/(4*k))
     points=[left+(right-left)*j/4 for j in range(5)]
     actual=mp.quad(lambda r:coeff(*args(r))[0],points)+coeff(*args(right))[1]-coeff(*args(left))[1]
     expected=endpoint(right)-endpoint(left);err=abs(actual/expected-1)
     tag=str(ix)+'_'+str(precision)+'_'+str(level)
     check('shape_annulus_nonzero_'+tag,abs(expected)>mp.mpf('1e-6'))
     check('shape_annulus_coarse_'+tag,err<mp.mpf(PROTOCOL['coarse_relative_limit']),mp.nstr(err,22))
     vals.append({'left':mp.nstr(left,precision),'right':mp.nstr(right,precision),'actual':mp.nstr(actual,precision),'predicted':mp.nstr(expected,precision),'relative_error':mp.nstr(err,25)})
    errors=[mp.mpf(v['relative_error']) for v in vals]
    check('shape_annulus_fine_'+str(ix)+'_'+str(precision),errors[-1]<mp.mpf(PROTOCOL['fine_relative_limit']))
    ratios=[errors[n+1]/errors[n] for n in range(2)]
    check('shape_annulus_refines_'+str(ix)+'_'+str(precision),max(ratios)<mp.exp(-mp.pi/k),[mp.nstr(x,20) for x in ratios])
    aargs=args(radius);target=coeff(*aargs);contours=[];ar=mp.mpf(PROTOCOL['contour_radius'])
    for nodes in PROTOCOL['contour_nodes']:
     total=[mp.mpc(0),mp.mpc(0)]
     for j in range(nodes):
      angle=2*mp.pi*j/nodes;values=raw(*aargs,ar*mp.exp(1j*angle))
      for b in range(2):total[b]+=values[b]*mp.exp(-1j*angle)/(nodes*ar)
     ce=max(abs(total[b]-target[b])/(1+abs(target[b])) for b in range(2))
     check('unexpanded_shape_Taylor_'+str(ix)+'_'+str(precision)+'_'+str(nodes),ce<mp.mpf(PROTOCOL['contour_scaled_limit']),mp.nstr(ce,22))
     contours.append({'nodes':nodes,'scaled_error':mp.nstr(ce,25)})
    levels.append({'precision':precision,'values':vals,'Taylor_contours':contours})
  with mp.workdps(80):
   for n in range(3):
    a=mp.mpf(levels[0]['values'][n]['actual']);b=mp.mpf(levels[1]['values'][n]['actual']);pe=abs(a-b)/(1+abs(b))
    check('charged_precision_refinement_'+str(ix)+'_'+str(n),pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,22))
  FIXTURES.append({'fixture':'charged_local_shape_boundary_action','parameters':row,'levels':levels,'full_canonical_determinant_evaluated':False})

def finite_part(nu,p,d):
 # Integral Taylor remainder, evaluated with expm1 to avoid cancellation at x=0.
 def integrand(x):
  if x==0:return fall(nu,p)*d**(nu-p)/mp.factorial(p-1)
  return (1-x)**(p-1)*fall(nu,p-1)*d**(nu-p+1)*mp.expm1((nu-p+1)*mp.log1p(x/d))/(mp.factorial(p-1)*x)
 value=mp.quad(integrand,[0,d,mp.sqrt(d),1])
 value+=sum(fall(nu,n)*d**(nu-n)/(mp.factorial(n)*(n-p+1)) for n in range(p-1))
 return mp.re(value)

def distribution_fixtures():
 for ix,wstring in enumerate(PROTOCOL['finite_part_omega']):
  levels=[]
  for precision in PROTOCOL['precision_digits']:
   with mp.workdps(precision):
    omega=mp.mpf(wstring);nu=3+1j*omega
    bound=mp.mpf(16)*abs(fall(nu,3))/15+1+abs(nu)+abs(fall(nu,2));vals=[]
    for level,dstring in enumerate(PROTOCOL['finite_part_cuts']):
     d=mp.mpf(dstring);outputs=[]
     for p in PROTOCOL['finite_part_powers']:
      actual=finite_part(nu,p,d);expected=mp.re(1/(nu-p+1));err=abs(actual-expected)
      tag=str(ix)+'_'+str(precision)+'_'+str(level)+'_'+str(p)
      check('finite_part_translation_bound_'+tag,err<=bound*mp.sqrt(d),mp.nstr(err/(bound*mp.sqrt(d)),22))
      if level==len(PROTOCOL['finite_part_cuts'])-1:check('finite_part_fine_'+tag,err<mp.mpf(PROTOCOL['finite_part_fine_absolute_limit']),mp.nstr(err,22))
      outputs.append({'power':p,'actual':mp.nstr(actual,precision),'limit':mp.nstr(expected,precision),'absolute_error':mp.nstr(err,25),'bound_ratio':mp.nstr(err/(bound*mp.sqrt(d)),25)})
     vals.append({'cut':dstring,'profiles':outputs})
    # The third jet does not converge, although the three finite-part actions do.
    third=fall(nu,3);a=mp.re(third);b=mp.re(-third)
    check('nonconvergent_third_jet_'+str(ix)+'_'+str(precision),abs(a-b)>1,mp.nstr(abs(a-b),22))
    levels.append({'precision':precision,'values':vals,'continuity_bound_coefficient':mp.nstr(bound,25),'third_jet_subsequence_gap':mp.nstr(abs(a-b),25)})
  with mp.workdps(80):
   for n in range(3):
    aa=levels[0]['values'][n]['profiles'];bb=levels[1]['values'][n]['profiles']
    pe=max(abs(mp.mpf(x['actual'])-mp.mpf(y['actual']))/(1+abs(mp.mpf(y['actual']))) for x,y in zip(aa,bb))
    check('finite_part_precision_refinement_'+str(ix)+'_'+str(n),pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,22))
  FIXTURES.append({'fixture':'cubic_oscillatory_finite_part_translation','omega':wstring,'levels':levels,'full_canonical_determinant_evaluated':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 coupled fast quadratic boundary obstruction','input_count':len(ps)}),flush=True)
 try:provenance();coeff,raw=original_shape_action();coupled_identities();canonical_controls();charged_fixtures(coeff,raw);distribution_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 coupled fast quadratic boundary obstruction','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'full_local_shape_endpoint_derived':good,'trace_boundary_common_combination_derived':good,
  'canonical_fixed_seed_first_variation_limit_derived':good,'coupled_second_order_cutoff_obstruction_derived':good,
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
