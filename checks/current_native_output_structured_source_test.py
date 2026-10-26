"""Native physical output and full structured-source identities, not a solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import factorial,comb
import json,sys,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
STEM='current-native-output-structured-source-test'
OWN=['notes/'+STEM+'.md','checks/current_native_output_structured_source_test.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/stationary-centered-first-order-geometry.md','notes/stationary-completed-metric-response-norm.md',
 'notes/comparison-path-interior-response.md','notes/divided-jet-dirac-response.md',
 'notes/finite-cold-admission-and-response.md','notes/current-controlled-solution-audit-closeout.md',
 'notes/current-finite-local-action-assembly.md','notes/current-schur-residual-phase-domain.md',
 'notes/current-global-mixed-schur-preconditioner.md','notes/current-direct-triple-return-attempt.md',
 'notes/current-interior-translation-jet-bound.md','notes/current-physical-residual-regularity-test.md',
 'evidence/current-physical-residual-regularity-test.json','evidence/current-localized-bolt-adjoint-sensitivity.json']
CHECKS=[]
def sha(p):return sha256(p.read_bytes()).hexdigest()
def check(n,ok):
 CHECKS.append({'name':n,'passed':bool(ok)});print(n+(': passed' if ok else ': FAILED'),flush=True)
 if not ok:raise AssertionError(n)
def enc(x):
 m,r,e=x.mid_rad_10exp(50)
 return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
 if isinstance(x,dict):
  d=x['decimal_enclosure'];return arb(d['mid']+' +/- '+d['rad'])*arb(10)**int(d['exponent'])
 return arb(x.strip('[]'))
def rat(x):return arb(int(x.p))/int(x.q)
def peval(p,x):
 ans=arb(0)
 for c in sp.Poly(p).all_coeffs():ans=ans*x+rat(c)
 return ans
def execute():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('ordinary_argument_frozen',all(sha(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
 pr=json.loads((ROOT/'evidence/current-physical-residual-regularity-test.json').read_bytes())['runs'][-1]
 bp=json.loads((ROOT/'evidence/current-localized-bolt-adjoint-sensitivity.json').read_bytes())['runs'][-1]
 check('previous_physical_test_accepted',pr['status']=='passed')
 check('actual_parameter_receipt_accepted',bp['status']=='passed')
 for run in [pr,bp]:
  for p,h in run['input_sha256_before'].items():
   if p in INPUTS and p.startswith('notes/'):check('unchanged_'+p,sha(ROOT/p)==h)
 pars=bp['result']['parameters'];h,zh,gamma,C,B2,beta=[dec(pars[k]) for k in ['h','zh','gamma','C','B2','beta']]
 s0=h/(1+h);R=arb('.1');Jv=R**2+2*R/100+arb(2).sqrt()/10000
 ts=gamma*Jv/(2*s0);Ms=1/(s0*(1-ts));Ss=s0+gamma*Jv/2
 check('whole_disk_inverse_multiplier_neumann',ts<arb('.007'))
 check('whole_disk_radius_multiplier_target',Ms<arb('2e33'))
 def nonlinear(rho):
  ex=(Ss*rho).exp();tt=rho*ex
  return {'remainder':Ss*rho*rho*ex/2+tt*tt/(2*(1-tt)),
   'first':ex/(1-tt),'second':(1+Ss)*ex/(1-tt)**2,'value':rho+Ss*rho*rho*ex/2+tt*tt/(2*(1-tt))}
 nl=nonlinear(arb('1e-8'))
 check('nonlinear_native_radius_remainder',nl['remainder']<arb('5.000001e-17'))
 check('nonlinear_native_radius_first',nl['first']<arb('1.00000002'))
 check('nonlinear_native_radius_second',nl['second']<arb('1.0000001'))
 qb=arb('1.667e-9');boltdebit=s0*qb*qb*(s0*qb).exp()/2
 check('actual_bolt_nonlinear_conversion_debit',boltdebit<arb('8e-52'))
 # Exact nonlinear physical-radius identity and both variations.
 ss,ww,zz,nn=sp.symbols('s w z nu',positive=True)
 ee=(sp.exp(ss*ww)-1)/ss;vv=sp.log(1+ee)
 check('native_radius_linearization',sp.simplify(sp.diff(vv,ww).subs(ww,0))==1)
 check('native_radius_second_derivative',sp.simplify(sp.diff(vv,ww,2)+(1-ss)*sp.exp(ss*ww)/(1+ee)**2)==0)
 check('native_radius_quadratic_coefficient',sp.simplify(sp.diff(vv,ww,2).subs(ww,0)-(ss-1))==0)
 check('actual_radius_completion',sp.simplify((1-zz)*sp.exp(sp.log((sp.exp(nn)-zz)/(1-zz)))+zz-sp.exp(nn))==0)
 # All regular packet jets, bounded on the WHOLE annulus by Bernstein coefficients.
 v,x,t=sp.symbols('v x t');p=sp.Poly((v-1)**16*(4-v)**16,v)
 kj=[];bern=[]
 for j in range(6):
  sq=sp.Integer(0)
  for b in range(j//2+1):
   term=sp.Rational(factorial(2*b),factorial(b))*sp.diff(p.diff((0,b)).as_expr().subs(v,x*x),x,j-2*b)
   sq+=comb(j,2*b)*term**2
  xx=sp.Poly(sp.expand(sq),x)
  check('radial_packet_tensor_even_'+str(j),all(k[0]%2==0 for k,c in xx.terms()))
  pv=sp.Poly(sum(c*v**(k[0]//2) for k,c in xx.terms()),v)
  pt=sp.Poly(pv.as_expr().subs(v,1+3*t),t);n=pt.degree()
  bb=[sum(pt.nth(jj)*sp.Rational(comb(k,jj),comb(n,jj)) for jj in range(k+1)) for k in range(n+1)]
  bound=max(bb);check('positive_whole_annulus_tensor_bound_'+str(j),bound>0)
  kj.append(rat(bound).sqrt());bern.append(str(bound))
 width=arb('1e-100')
 jnu=sum((width**(6-j)*kj[j]/(factorial(j)*100**j) for j in range(6)),arb(0))
 jw=Ms*jnu;packetnl=nonlinear(jw)
 check('finite_packet_inside_native_radius_ball',jw<arb('1e-40'))
 # Every term of the physical areal derivative, without asymptotic substitution.
 z=sp.Symbol('z');cp=[sp.Integer(1),z]
 for _ in range(5):
  nxt=[sp.Integer(0)]*(len(cp)+1)
  for k,c in enumerate(cp):nxt[k]+=-z*z*sp.diff(c,z);nxt[k+1]+=-z*z*c
  cp=[sp.expand(c) for c in nxt]
 check('original_areal_conversion_top',cp[6]==-z**11)
 zpoint=zh-gamma*width**2*arb('2.5')/2
 jets=[rat(p.diff((0,j)).eval(sp.Rational(5,2)))*(-2/(gamma*width**2))**j*width**6 for j in range(7)]
 cv=[peval(c,zpoint) if c!=0 else arb(0) for c in cp]
 areal=(abs(cv[6]*jets[6])-sum((abs(cv[j]*jets[j]) for j in range(6)),arb(0)))/B2**(arb(5)/2)
 check('finite_packet_physical_areal_derivative_large',areal>arb('1e600'))
 # Free ordered operator algebra; these letters do not represent selected modes.
 A,Ai,CC,Ca,Cai,DD,TT=sp.symbols('A Ai C Ca Cai D T',commutative=False)
 pairs={(A,Ai),(Ai,A),(Ca,Cai),(Cai,Ca)}
 def reduce(expr):
  answer=sp.Integer(0)
  for term in sp.Add.make_args(sp.expand(expr)):
   coeff,rest=term.as_coeff_Mul();word=list(sp.Mul.make_args(rest));changed=True
   while changed:
    changed=False
    for j in range(len(word)-1):
     if (word[j],word[j+1]) in pairs:
      word=word[:j]+word[j+2:];changed=True;break
   answer+=coeff*sp.Mul(*word)
  return sp.expand(answer)
 source=(DD+Ca*Ai*CC)*TT;nu=Cai*DD*TT
 check('complete_signed_scalar_reconstruction',reduce(nu-(-Ai*CC*TT+Cai*source))==0)
 check('complete_original_second_source_row_zero',reduce(-Ca*nu+DD*TT)==0)
 check('complete_original_first_source_row',reduce(-A*nu-CC*TT+A*Cai*source)==0)
 # Exact physical normal factors and the continuous Fourier multiplier bound.
 k,be,cv0,rd=sp.symbols('k beta Cv rD')
 jc=-6*be*(cv0-rd*k*(k-3))*(k-4)*(k+1)
 ca0=2*be*(cv0-rd*k*(k-3))*(6-k)*(k+1)
 check('complete_normal_composition_ratio',sp.cancel(jc/ca0-3*(k-4)/(k-6))==0)
 check('structured_quartic_source_cancels',jc.subs(k,4)==0)
 sig=sp.Rational('4.005');tt=sp.Symbol('t',real=True)
 gap=sp.factor(9-9*((sig-4)**2+tt**2)/((sig-6)**2+tt**2))
 check('normal_H3_norm_at_most_three',sp.simplify(gap-9*((sig-6)**2-(sig-4)**2)/((sig-6)**2+tt**2))==0 and (sig-6)**2>(sig-4)**2)
 check('normal_zero_frequency_value',3*(sig-4)/(6-sig)==sp.Rational(1,133))
 lam=C*beta/B2;eps=1/(60*arb.pi());Cv=1+11*eps*lam/3;rD=(96*arb.pi()*arb('1e22')-eps)*lam
 mixedquartic=20*beta*(Cv-4*rD)
 check('actual_mixed_quartic_nonresonant',mixedquartic>3 and mixedquartic<4)
 # Full P32 prefactor, with original q and degree, no change of norm.
 q=arb('.338');xx=-1-2/q;ch0=arb(1);ch1=xx
 for _ in range(2,30):ch0,ch1=ch1,2*xx*ch1-ch0
 paff=arb('1.1e-10')*(1+q)**3/abs(ch1)
 ceiling=arb('5e-9')/(paff*6*Ms)
 check('conditional_native_slice_ceiling_finite',ceiling>0 and ceiling.is_finite())
 return {'precision_bits':512,
  'whole_bolt_disk_radius_conversion':{k:enc(vv) for k,vv in {'s0':s0,'J5_radius_squared_bound':Jv,
   'relative_multiplier_series_bound':ts,'inverse_radius_J5_multiplier':Ms,'s_J5_bound':Ss,
   'P5_to_P5_radius_multiplier':6*Ms,'nonlinear_remainder_at_1e_8':nl['remainder'],
   'nonlinear_first_bound':nl['first'],'nonlinear_second_bound':nl['second'],
   'actual_bolt_gap_nonlinear_debit_at_qb':boltdebit}.items()},
  'finite_actual_domain_packet':{'width':enc(width),'exact_tensor_Bernstein_upper_squares':bern,
   'canonical_J5_upper':enc(jnu),'linear_native_radius_J5_upper':enc(jw),
   'nonlinear_native_radius_J5_upper':enc(packetnl['value']),'physical_areal_fifth_lower':enc(areal)},
  'complete_structured_source':{'second_row':0,'quartic_source_finitepart':0,'normal_H3_multiplier_bound':3,
   'normal_zero_frequency_absolute_value':'1/133','mixed_quartic_normal_value':enc(mixedquartic),
   'finite_geometry_H3_constant':None,'original_Schur_Y_comparison':None},
  'budget':{'P32_affine_prefactor':enc(paff),'conditional_Katlas_Mcanonical_ceiling_5e_9':enc(ceiling),
   'exact_structured_c4_bolt_debit':0,'actual_numerical_c4_error':None},
  'scope':{'original_geometric_admission_scope_recovered':True,'areal_counterexamples_retracted':False,
   'complete_X_defined':False,'native_radius_multiplier_certified':True,'structured_source_c4_zero':True,
   'full_structured_source_norm_computed':False,'complete_coupled_estimate_certified':False,
   'original_defect_certified':False,'actual_mass_target_certified':False,'quantum_inverse_applied':False,
   'higher_degree_adopted':False,'normalization_improvement_claimed':False,'new_physical_gate_closed':False,
   'construction_obstruction':False}}
if __name__=='__main__':
 mode='acceptance' if '--acceptance' in sys.argv else 'development'
 out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
 hashes={p:sha(ROOT/p) for p in INPUTS};run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
  'input_sha256_before':hashes,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
  'runtime':{'python':sys.version,'python_sha256':sha(Path(sys.executable)),'flint_version':flint.__version__,'sympy_version':sp.__version__}}
 try:
  run['result']=execute();run['input_sha256_after']={p:sha(ROOT/p) for p in INPUTS}
  check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
 except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
 run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
 ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};ledger['runs'].append(run)
 out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
 summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
 with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(summary.encode())
 print(summary,flush=True);sys.exit(0 if run['status']=='passed' else 1)
