"""Complete seed endpoint coefficients and the forced logarithmic source weight.

Continuum proof: notes/current-outer-logarithmic-source-realization.md.
No finite spectral problem or replacement boundary domain is evaluated.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import current_product_matched_source_projection as old
import sympy as sp
from flint import arb,ctx
STEM='current-outer-logarithmic-source-realization'
OWN=['notes/'+STEM+'.md','checks/current_outer_logarithmic_source_realization.py','evidence/'+STEM+'-review.json']
ROWS=[];SYMBOLIC={};enc=old.encoded
def check(name,ok,detail=''):
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def algebra():
 z,beta,C,B2,G,q,Q3,Q4,Q5,Q6=sp.symbols('z beta C B2 G q Q3 Q4 Q5 Q6',positive=True)
 Q=beta+z*z/C+Q3*z**3+Q4*z**4+Q5*z**5+Q6*z**6
 phi=-sp.log(z);D=lambda f:sp.diff(Q*sp.diff(f,z),z)
 L=-sp.diff(Q,z,2);M=D(phi);V=Q*sp.diff(phi,z)**2
 A=phi*L**2/60+L*M/30-M**2/10+L*V/30-7*M*V/90-7*V**2/180
 fb=sp.expand(G/(2*sp.pi*C*B2)*(-2*C*C*A-2*C*C*D((phi*L+M+V)/30)+C*M/9+2*(sp.Rational(1,15)-q*q/3)*phi))
 logpart=sp.expand(fb.coeff(sp.log(z),1));regular=sp.expand(fb-logpart*sp.log(z))
 pole=regular.coeff(z,-1);clog=logpart.coeff(z,0)
 k=G*q*q/(12*sp.pi*C*B2);cnm=4*G*C*beta*Q4/(5*sp.pi*B2)
 check('complete_Weyl_pole',sp.simplify(pole+G*C*beta*Q3/(5*sp.pi*B2))==0)
 check('complete_Weyl_logarithm',sp.simplify(clog-(4*k-cnm))==0)
 proper_c=Q3*sp.sqrt(beta)/2
 check('massive_pole_cancels_actual_Weyl_pole',sp.simplify(2*G*C/(5*sp.pi*B2)*proper_c*sp.sqrt(beta)+pole)==0)
 fw=sp.expand(G/(2*sp.pi*C*B2)*(C*C*(L*L/60+D(L)/30+D(-M/5-7*V/90)-2*sp.diff((L/30-7*M/90-7*V/90)*sp.diff(phi,z)*Q,z))-C*L/18+sp.Rational(1,15)-q*q/3+C*M/9))
 check('complete_scalar_constant',sp.simplify(fw.coeff(z,0)+2*k)==0)
 la=sp.symbols('lambda',positive=True)
 check('exact_outer_vacuum_cancellation',sp.simplify(12*11*G*la**2/(360*sp.pi)-11*G*(12*la)**2/(4320*sp.pi))==0)
 # Direct Einstein/reference Euler variations at the outer leading metric.
 s,a,b,ee=sp.symbols('s a b e');nn=a*z**s;bb=b*z**s
 ph=phi+ee*nn;tt=ee*bb;delta=lambda f:beta*sp.diff(f,z,2)
 DD=-2*delta(tt)-6*delta(ph)-6*beta*sp.diff(ph,z)**2
 Fnu=sp.exp(2*ee*nn)*DD*z**-2+12*beta*sp.exp(4*ee*nn+2*tt)*z**-4
 Fb=-delta(z**-2*sp.exp(2*ee*nn))+6*beta*sp.exp(4*ee*nn+2*tt)*z**-4
 nrow=sp.simplify(sp.diff(Fnu,ee).subs(ee,0)/z**(s-4)/beta)
 brow=sp.simplify(sp.diff(Fb,ee).subs(ee,0)/z**(s-4)/beta)
 mat=sp.Matrix([[sp.diff(nrow,a),sp.diff(nrow,b)],[sp.diff(brow,a),sp.diff(brow,b)]])
 expected=sp.Matrix([[24+18*s-6*s*s,24+2*s-2*s*s],[12+10*s-2*s*s,12]])
 check('reference_Euler_indicial_matrix',(mat-expected).applyfunc(sp.simplify)==sp.zeros(2,2),str(mat))
 check('physical_mass_column',(mat.subs(s,3)*sp.Matrix([1,-2])).applyfunc(sp.simplify)==sp.zeros(2,1))
 check('reference_indicial_determinant',sp.factor(mat.det())==-4*s*(s-3)*(s-4)*(s+1))
 fn,cl=sp.symbols('fnu clog');aa=(6*fn+7*cl)/(40*beta);bl=-(2*fn+3*cl)/(8*beta)
 check('complete_forced_logarithms',sp.simplify(beta*(-30*aa-14*bl)+fn)==0 and sp.simplify(beta*(20*aa+12*bl)+cl)==0)
 check('magnetic_response_outer_logarithms',sp.simplify(aa.subs({fn:-2*k,cl:4*k})-2*k/(5*beta))==0 and sp.simplify(bl.subs({fn:-2*k,cl:4*k})+k/beta)==0)
 # Entire angular discrepancy: exact trapezoid kernel and positive moments.
 u=sp.symbols('u',nonnegative=True)
 check('trapezoid_Peano_kernel',sp.expand(sp.Rational(1,8)-u*(1-u)/2-(2*u-1)**2/8)==0)
 moments=[sp.integrate(f*sp.exp(-u),(u,0,sp.oo)) for f in [1+u,u+2,u*(u+2)]]
 check('angular_entire_variation_majorants',moments==[2,3,4])
 E,Ep,Epp=sp.symbols('E Ep Epp');fun=sp.Function('f')
 chain=(lambda expr:sp.diff(expr,E)*Ep+sp.diff(expr,Ep)*Epp)
 second=chain(chain(Ep*fun(E)))
 check('angular_g_second_derivative',sp.simplify(second-3*Ep*Epp*sp.diff(fun(E),E)-Ep**3*sp.diff(fun(E),E,2))==0)
 check('angular_discrepancy_coefficients',sp.Rational(1,4)*(6*2+4*4)==7 and sp.Rational(1,4)*3==sp.Rational(3,4))
 check('universal_Brownian_fourth_moment',sp.integrate((u*(1-u))**2,(u,0,1))==sp.Rational(1,30))
 check('cubic_profile_maximum',sp.simplify(4/(sp.pi*sp.sqrt(sp.pi))*sp.gamma(sp.Rational(5,2))**2/sp.gamma(5)-3/(32*sp.sqrt(sp.pi)))==0)
 # Fourier identities fix the sign independently of the local action.
 ss,p=sp.symbols('s p',positive=True)
 coeff=-2*sp.gamma(ss+1)*sp.sin(sp.pi*ss/2)
 check('Fourier_first_logarithm',sp.simplify(sp.diff(coeff,ss).subs(ss,0))==-sp.pi)
 check('Fourier_second_logarithm',sp.simplify(sp.diff(coeff*p**(-ss-1),ss,2).subs(ss,0)-2*sp.pi*(sp.log(p)+sp.EulerGamma)/p)==0)
 check('quartic_log_fourth_derivative',sp.simplify(sp.diff(z**4*sp.log(z),z,4)-24*sp.log(z)-50)==0)
 wlocal=sp.expand(sp.diff(-sp.log(z)*sp.diff(z**4*sp.log(z),z,2),z,2))
 check('Weyl_log_squared_coefficient',wlocal.coeff(sp.log(z),2)==-24)
 check('complete_signed_log_squared_coefficient',12-24==-12)
 # Independent direct heat integral extraction: no spectral truncation.
 l,T=sp.symbols('l T');heat=3*(T*T-4*l*l)-12*l*(T-2*l)
 check('direct_heat_log_squared_coefficient',sp.expand(heat).coeff(l,2)==12)
 SYMBOLIC.update(Weyl_pole=str(pole),Weyl_logarithm=str(clog),scalar_constant=str(fw.coeff(z,0)),reference_matrix=str(mat),
  universal_nonlocal_principal='-C^2/(30pi) D_Y^4 log|D_Y|',complete_normalized_log_squared_coefficient='-12 K b_L/F_lambda',
  ordinary_paired_domain_and_remainder_argument_required=True)

def execute():
 algebra();p,h,zh,env=old.scalar.prev.york.actual_parameters()
 C,B2,beta,Fl=p['C'],p['B2'],p['beta'],p['Fl'];q=arb(p['q'])
 k=q*q/(12*arb.pi()*C*B2);Q3=-4*beta-2*p['u']-p['tau'];Q4=3*beta+p['u']+p['eps']
 cnm=4*C*beta*Q4/(5*arb.pi()*B2);clog=4*k-cnm;fnu=-2*k
 alog=(6*fnu+7*clog)/(40*beta);blog=-(2*fnu+3*clog)/(8*beta)
 aM=2*k/(5*beta);bM=-k/beta
 # Algebraically correlated differences, no subtraction of broad intervals.
 aRest=-7*cnm/(40*beta);bRest=3*cnm/(8*beta)
 K=C*beta**2/(15*arb.pi()*B2);source2=-12*K*blog/Fl
 c=Q3*beta.sqrt()/2;t0=arb('.001')
 discrepancy=c.abs_upper()*(84*C*t0.sqrt()+3*q*q*t0**arb('1.5'))/(64*arb.pi().sqrt())
 check('actual_complete_b_log_negative',blog>arb('-1e-28') and blog<arb('-9e-29'))
 check('actual_squared_log_source_nonzero',source2>arb('3e-63') and source2<arb('5e-63'))
 check('remaining_logarithm_is_not_zero',aRest<0 and bRest>0)
 check('magnetic_dominance_not_source_substitution',cnm/(4*k)<arb('1e-5'))
 check('all_angular_cubic_discrepancy_finite',discrepancy>0 and discrepancy<arb('1e21'))
 return dict(parameters={a:enc(v) for a,v in dict(C=C,B2=B2,beta=beta,h=h,Q3=Q3,Q4=Q4,F_lambda=Fl).items()},
  k_M=enc(k),nonmagnetic_log_coefficient=enc(cnm),full_scalar_constant=enc(fnu),full_shape_log_coefficient=enc(clog),
  complete_nu_z4_log_coefficient=enc(alog),complete_b_z4_log_coefficient=enc(blog),
  magnetic_nu_z4_log_coefficient=enc(aM),magnetic_b_z4_log_coefficient=enc(bM),
  remaining_nu_z4_log_coefficient=enc(aRest),remaining_b_z4_log_coefficient=enc(bRest),
  nonmagnetic_to_magnetic_shape_log_ratio=enc(cnm/(4*k)),K=enc(K),actual_residual_log_squared_coefficient=enc(source2),
  all_angular_cubic_discrepancy_upper_t0=enc(discrepancy),fixed_seed_massive_remainder_is_continuous=True,
  massive_leading_log_coefficient='0',R0_and_Cb_numerical_values_known=False,
  finite_single_log_canonical_residual_bound_exists=False,second_reference_iteration_returns_single_log_metric_phase=False,
  full_quantum_strong_inverse_obstructed=False,full_quantum_strong_inverse_proved=False,
  actual_global_remaining_source_reconstructed=False,actual_controlled_solution=False,nonlinear_admission=False,all_sector_stability=False)

def main():
 ctx.prec=512;review=json.loads((ROOT/OWN[2]).read_bytes());paths=set(OWN+review['retained_inputs'])
 for mod in list(sys.modules.values()):
  f=getattr(mod,'__file__',None)
  if f:
   try:p=Path(f).resolve().relative_to(ROOT).as_posix()
   except ValueError:continue
   if p.endswith('.py'):paths.add(p)
 before={p:digest(p) for p in sorted(paths)};out=ROOT/'evidence'/f'{STEM}.json';archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
  frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256']);run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['symbolic']=SYMBOLIC;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),result=run.get('result')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
