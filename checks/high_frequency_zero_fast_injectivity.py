"""Scaled compatible Green inverse and nonnormal perturbation diagnostics."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='high-frequency-zero-fast-injectivity'
PREVIOUS='causal-preparation-source-map'
OWN=['checks/high_frequency_zero_fast_injectivity.py','notes/'+STEM+'.md','evidence/high-frequency-zero-fast-source-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-29','precision_tolerance':'1e-29',
 'Bessel_cases':[['1.3','1'],['3.7','2']],'Bessel_columns':['k','gamma'],
 'Green_readouts':['.5','1.4','2.5'],'Wronskian_readout':'1.1',
 'shift_gamma':[2,4,8],'Neumann_orders':[1,3,6],'shift_remainder_norm_bound':'1',
 'compact_profile':'4096*(y-1)^6*(2-y)^6 on [1,2]; C5 control in D_min',
 'actual_C_a_and_M_star_evaluated':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.simplify(S.sympify(value).doit());check(name,value==0,value)

def zeromat(name,value):
 value=value.applyfunc(S.simplify);check(name,value==S.zeros(*value.shape),value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==78 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['high_frequency_auxiliary_injectivity_proved'] and review['single_compatibility_condition_above_threshold']
       and not review['actual_numeric_threshold_certified'])

def scaling_identities():
 y,z,s,gamma,k,a=S.symbols('y z s gamma k a',positive=True);f=S.Function('f')
 ug=S.sqrt(gamma)*f(gamma*y)
 rhs=gamma**2*S.sqrt(gamma)*(-S.diff(f(z),z,2)-a*f(z)/z**2+f(z)).subs(z,gamma*y)
 zero('exact_minimal_operator_dilation',-S.diff(ug,y,2)-a*ug/y**2+gamma**2*ug-rhs)
 zero('unitary_dilation_Jacobian',gamma/gamma-1)
 zero('scaled_Green_kernel_prefactor',S.sqrt(gamma*y)*S.sqrt(gamma*s)/gamma-S.sqrt(y*s))
 b=S.sqrt(y)*f(gamma*y)
 model=-S.diff(b,y,2)-(k*k+S.Rational(1,4))*b/y**2+gamma**2*b
 bessel=-gamma**2*S.sqrt(y)*(S.diff(f(z),z,2)+S.diff(f(z),z)/z+(k*k/z**2-1)*f(z)).subs(z,gamma*y)
 zero('imaginary_order_Bessel_reduction',model-bessel)

def infinite_shift_controls():
 I=S.I;V=S.Matrix([[S.Rational(1,5)+I/10,-S.Rational(1,7)],[S.Rational(1,9),I/8],[-I/11,S.Rational(1,6)]])
 shift=S.Matrix([[0,0],[1,0],[0,1]]);b=S.Matrix([1,0,0])
 norm2=lambda v:S.simplify((v.conjugate().T*v)[0])
 Vbound=S.simplify(sum(S.conjugate(t)*t for t in V))
 check('finite_rank_remainder_Frobenius_bound',Vbound<1,Vbound)
 for gamma in PROTOCOL['shift_gamma']:
  tag=str(gamma);eta=S.Rational(1,gamma**2);A=gamma**2*shift;T=shift.T/gamma**2
  L=A+V;z=(S.eye(3)+T.conjugate().T*V.conjugate().T).inv()*b
  J=(S.eye(2)+T*V).inv()*T;Pi=S.eye(3)-b*z.conjugate().T
  check('perturbation_strict_threshold_'+tag,eta<1,eta)
  zeromat('compatible_left_inverse_'+tag,J*L-S.eye(2))
  zeromat('true_adjoint_null_'+tag,L.conjugate().T*z)
  zero('adjoint_normalization_'+tag,(b.conjugate().T*z)[0]-1)
  zeromat('unprojected_source_residual_'+tag,L*J-Pi)
  zeromat('oblique_projection_idempotence_'+tag,Pi*Pi-Pi)
  check('adjoint_norm_bound_'+tag,norm2(z)<=1/(1-eta)**2)
  check('adjoint_distance_bound_'+tag,norm2(z-b)<=(eta/(1-eta))**2)
  gram=(L.conjugate().T*L-(gamma**2-1)**2*S.eye(2)).applyfunc(S.simplify)
  check('compatible_inverse_norm_bound_'+tag,gram[0,0]>0 and gram.det()>0)
  N=T.conjugate().T*V.conjugate().T;records=[]
  for order in PROTOCOL['Neumann_orders']:
   approx=S.zeros(3,1);term=b
   for n in range(order+1):approx+=term;term=-N*term
   error2=norm2(z-approx);bound2=(eta**(order+1)/(1-eta))**2
   ratio=S.simplify(error2/bound2)
   check('Neumann_tail_bound_'+tag+'_'+str(order),error2<=bound2,S.N(ratio,30))
   records.append({'order':order,'squared_error':str(error2),'squared_bound':str(bound2),'bound_ratio':str(S.N(ratio,30))})
  f=S.Matrix([1+I,2,-1+2*I]);residual=(f-L*J*f).applyfunc(S.simplify)
  check('incompatible_source_not_solved_'+tag,norm2(residual)>0)
  FIXTURES.append({'fixture':'finite_rank_infinite_shift','gamma':gamma,'eta_bound':str(eta),
   'adjoint_distance_squared':str(norm2(z-b)),'residual_coefficient':str(S.simplify((z.conjugate().T*f)[0])),
   'Neumann_tails':records,'tail':'gamma^2 e_n maps to e_(n+1) for every n>=2',
   'physical_charged_kernel_approximation':False})

def Bessel_fixtures():
 y,a,g=S.symbols('y a gamma',positive=True);u=4096*(y-1)**6*(2-y)**6
 forcing=S.lambdify((y,a,g),-S.diff(u,y,2)-a*u/y**2+g*g*u,'mpmath')
 for case,(ks,gs) in enumerate(PROTOCOL['Bessel_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    k=mp.mpf(ks);gamma=mp.mpf(gs);aa=k*k+mp.mpf('.25');nu=1j*k
    h=lambda x:mp.sqrt(x)*mp.besseli(nu,gamma*x)
    b=lambda x:mp.sqrt(x)*mp.besselk(nu,gamma*x)
    f=lambda x:forcing(x,aa,gamma)
    yr=mp.mpf(PROTOCOL['Wronskian_readout'])
    wr=h(yr)*mp.diff(b,yr)-mp.diff(h,yr)*b(yr)
    error=abs(wr+1)
    check(f'Bessel_Wronskian_{case}_{digits}',error<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(error,digits))
    residual=-mp.diff(b,yr,2)-aa*b(yr)/yr**2+gamma**2*b(yr)
    error=abs(residual)/(1+abs(b(yr)))
    check(f'Bessel_homogeneous_equation_{case}_{digits}',error<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(error,digits))
    moment=mp.quad(lambda x:b(x)*f(x),[1,mp.mpf('1.5'),2])
    check(f'compatible_horizon_moment_{case}_{digits}',abs(moment)<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(abs(moment),digits))
    yy=mp.mpf('.7');ss=mp.mpf('1.8')
    kernel=h(yy)*b(ss)
    scaled=mp.sqrt(gamma*yy)*mp.besseli(nu,gamma*yy)*mp.sqrt(gamma*ss)*mp.besselk(nu,gamma*ss)/gamma
    error=abs(kernel-scaled)/(1+abs(kernel))
    check(f'Green_kernel_scaling_{case}_{digits}',error<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(error,digits))
    values=[]
    for j,ys in enumerate(PROTOCOL['Green_readouts']):
     yy=mp.mpf(ys);left=min(yy,mp.mpf(2));right=max(yy,mp.mpf(1))
     intleft=mp.quad(lambda x:h(x)*f(x),[1,(1+left)/2,left]) if left>1 else mp.mpf(0)
     intright=mp.quad(lambda x:b(x)*f(x),[right,(right+2)/2,2]) if right<2 else mp.mpf(0)
     value=b(yy)*intleft+h(yy)*intright
     target=4096*(yy-1)**6*(2-yy)**6 if 1<yy<2 else mp.mpf(0)
     error=abs(value-target)/(1+abs(target))
     check(f'compatible_Green_reconstruction_{case}_{digits}_{j}',error<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(error,digits))
     values.append({'y':ys,'real':mp.nstr(value.real,digits),'imag':mp.nstr(value.imag,digits)})
    levels.append({'precision':digits,'reconstructed':values,'moment_real':mp.nstr(moment.real,digits),'moment_imag':mp.nstr(moment.imag,digits)})
  with mp.workdps(80):
   errors=[]
   for va,vb in zip(levels[0]['reconstructed'],levels[1]['reconstructed']):
    aa=mp.mpc(va['real'],va['imag']);bb=mp.mpc(vb['real'],vb['imag']);errors.append(abs(aa-bb)/(1+abs(bb)))
   check('Bessel_precision_'+str(case),max(errors)<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(max(errors),80))
  FIXTURES.append({'fixture':'inverse_square_compatible_Green','k':ks,'gamma':gs,'levels':levels,
                   'benchmark_large_order_threshold_evaluated':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 high frequency zero-fast injectivity','input_count':len(ps)}),flush=True)
 try:provenance();scaling_identities();infinite_shift_controls();Bessel_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 high frequency zero-fast injectivity','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'high_frequency_auxiliary_injectivity_proved':good,'single_compatibility_condition_above_threshold':good,'compatible_inverse_gamma_minus_two_bound':good,
  'actual_numeric_threshold_certified':False,'actual_adjoint_profiles_numerically_evaluated':False,'actual_charged_two_point_kernel_evaluated':False,
  'physical_fast_branch_excluded':False,'finite_counterterm_changed':False,'new_boundary_prescription_adopted':False,'physical_source_projection_adopted':False,
  'order_reduction_adopted':False,'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]};data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
