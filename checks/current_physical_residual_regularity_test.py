"""Original physical-output and full quantum source regularity tests; not a solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import factorial,comb
import json,sys,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
import sympy as sp
from sympy.functions.combinatorial.numbers import stirling
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from operators.signed_local_core import SignedLocalCore
STEM='current-physical-residual-regularity-test'
OWN=['notes/'+STEM+'.md','checks/current_physical_residual_regularity_test.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['operators/signed_local_core.py','notes/current-finite-local-action-assembly.md',
 'notes/current-bolt-source-domain-test.md','evidence/current-bolt-source-domain-test.json',
 'evidence/current-global-mass-conversion.json','evidence/current-localized-bolt-adjoint-sensitivity.json',
 'notes/current-bolt-resummed-dirac-parametrix.md','notes/current-controlled-solution-audit-closeout.md',
 'notes/current-complete-residual-admission-test.md','notes/current-schur-physical-output-test.md',
 'notes/current-interior-translation-jet-bound.md','notes/current-common-regulator-strong-limit-test.md']
CHECKS=[]
def sh(p):return sha256(p.read_bytes()).hexdigest()
def check(n,ok):
 CHECKS.append({'name':n,'passed':bool(ok)});print(n+(': passed' if ok else ': FAILED'),flush=True)
 if not ok:raise AssertionError(n)
def enc(x):
 m,r,e=x.mid_rad_10exp(50)
 return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
 if isinstance(x,dict):
  e=x['decimal_enclosure'];return arb(e['mid']+' +/- '+e['rad'])*arb(10)**int(e['exponent'])
 return arb(x.strip('[]'))
def prior(stem):
 r=json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]
 check('accepted_'+stem,r['status']=='passed')
 for p,h in r['input_sha256_before'].items():
  if p in INPUTS and p.startswith(('notes/','operators/')):check('unchanged_'+p,sh(ROOT/p)==h)
 return r['result']
def evaluate(expr,env):
 if expr in env:return env[expr]
 if expr==sp.pi:return arb.pi()
 if expr.is_Rational:return arb(int(expr.p))/int(expr.q)
 if expr.is_Add:return sum((evaluate(e,env) for e in expr.args),arb(0))
 if expr.is_Mul:
  ans=arb(1)
  for e in expr.args:ans*=evaluate(e,env)
  return ans
 if expr.is_Pow:
  base,power=expr.args
  return evaluate(base,env)**(int(power) if power.is_Integer else evaluate(power,env))
 if expr.func==sp.log:return evaluate(expr.args[0],env).log()
 raise ValueError(str(expr))
def execute():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('ordinary_argument_and_checker_frozen',all(sh(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
 residual=prior('current-bolt-source-domain-test')
 pars=prior('current-global-mass-conversion')['parameters']
 bp=prior('current-localized-bolt-adjoint-sensitivity')['parameters']
 beta,u,eps,tau=[dec(pars[k]) for k in ['beta','u','eps','tau']]
 C,B2=[dec(bp[k]) for k in ['C','B2']]
 core=SignedLocalCore();z=core.z
 lam=core.C*core.beta/core.B2
 vacuum={core.lam0:lam+11*lam**2/(360*sp.pi)}
 rows=[[sp.Integer(0) for _ in range(5)] for _ in range(2)]
 for allocation in ['Einstein','Maxwell','R2','Weyl']:
  op=core.original_A_coefficients(allocation)
  for a in range(2):
   for k in range(5):rows[a][k]+=op[a][0][k]
 coeff=[[sp.factor(z**k*sum(rows[a][j]*stirling(j,k,kind=2) for j in range(k,5))).subs(vacuum)
         for k in range(5)] for a in range(2)]
 dd=(96*sp.pi*core.cR-1/(60*sp.pi))*core.C/core.B2
 for a,c in enumerate([6,2]):
  check('complete_original_fourth_coefficient_row_'+str(a),sp.simplify(coeff[a][4]-c*dd*z**4*core.Q[0]**2)==0)
 # Independent recursion for the proper fifth derivative, using Q jet symbols.
 q=core.Q
 deriv=lambda f:sp.diff(f,z)+sum(q[j+1]*sp.diff(f,q[j]) for j in range(7))
 lp=[sp.Integer(1)]
 for n in range(5):
  new=[sp.Integer(0)]*(n+2)
  for k,c in enumerate(lp):
   new[k]+=sp.sqrt(q[0])*deriv(c);new[k+1]+=sp.sqrt(q[0])*c
  lp=[sp.simplify(v) for v in new]
 expected=[sp.Integer(0),sp.sqrt(q[0])*(q[2]**2/4+3*q[1]*q[3]/4+q[0]*q[4]/2),
  sp.sqrt(q[0])*(15*q[1]*q[2]/4+5*q[0]*q[3]/2),
  sp.sqrt(q[0])*(15*q[1]**2/4+5*q[0]*q[2]),5*q[0]**sp.Rational(3,2)*q[1],q[0]**sp.Rational(5,2)]
 for k in range(6):check('whole_proper_fifth_derivative_coefficient_'+str(k),sp.simplify(lp[k]-expected[k])==0)
 zsym=sp.Symbol('t');qpoly=core.beta+zsym**2/core.C-(4*core.beta+2/core.C+core.tau)*zsym**3+(3*core.beta+1/core.C+core.eps)*zsym**4-2*core.eps*zsym**5+core.eps*zsym**6
 env0={core.beta:beta,core.C:C,core.B2:B2,core.eps:eps,core.tau:tau,core.cR:arb('1e22'),
       core.qflux:arb(2000)**4,core.gR:arb(2000)**(-4)}
 def env_at(at):
  env={**env0,z:at,zsym:at}
  env.update({q[j]:evaluate(sp.diff(qpoly,zsym,j),env) for j in range(9)})
  return env
 band=arb('.5 +/- .01');envband=env_at(band)
 U=[[evaluate(c,envband).abs_upper() for c in row] for row in coeff]
 check('all_original_coefficients_finite',all(v.is_finite() for row in U for v in row))
 check('actual_band_uniformly_elliptic',envband[q[0]]>0)
 e=arb('1e-100');amplitude=arb('1e-75');ell=arb('.005')
 cutoff=[arb(1)]+[arb(factorial(21)//factorial(21-j))*2**(j-1)/ell**j for j in range(1,5)]
 G=[];gam=arb.pi().sqrt()
 for j in range(5):
  if j:gam*=arb(2*j-1)/2
  G.append(gam.sqrt())
 row_upper=[]
 for a in range(2):
  row_upper.append(amplitude*sum((U[a][k]*e**(4-k)*sum((comb(k,j)*cutoff[j]*e**j*G[k-j] for j in range(k+1)),arb(0)) for k in range(5)),arb(0)))
 Ri=(sum((r*r for r in row_upper),arb(0))).sqrt()
 envpoint=env_at(arb('.5')+e)
 cpoint=[evaluate(c,envpoint) for c in expected]
 lower=amplitude*arb('-.5').exp()*e**(-arb(3)/2)*(6*cpoint[5]-2*e*abs(cpoint[4])-2*e*e*abs(cpoint[3])-e**4*abs(cpoint[1]))
 check('test_point_inside_cutoff_plateau',arb('.5')+e<arb('.505'))
 check('complete_both_row_residual_passes_L2_bolt_target',Ri<arb('8.800093264e-85'))
 boltcoef=dec(residual['complete_mixed_chart_residual_test']['whole_interior_L2_pair_to_bolt'])
 boltdebit=boltcoef*Ri
 check('complete_bolt_upper_estimator_passes_half_margin',boltdebit<arb('8.335e-10'))
 check('actual_proper_fifth_derivative_exceeds_001',lower>arb('.01'))
 check('source_not_assigned_zero',all(r>0 for r in row_upper))
 envcenter=env_at(arb('.5'));Qc=envcenter[q[0]]
 dphys=evaluate(dd,envcenter)
 A=arb(40).sqrt()*dphys*arb('.5')**4*Qc**2*(105*arb.pi().sqrt()/16).sqrt()
 B=6*arb('-.5').exp()*Qc**(arb(5)/2)
 check('strict_asymptotic_coefficient',A>0 and B>0)

 # Original areal lapse output, including all lower derivative terms.
 cp=[sp.Integer(1),z]
 for _ in range(5):
  nxt=[sp.Integer(0)]*(len(cp)+1)
  for k,c in enumerate(cp):
   nxt[k]+=-z**2*sp.diff(c,z);nxt[k+1]+=-z**2*c
  cp=[sp.expand(c) for c in nxt]
 check('physical_areal_top_coefficient',cp[6]==-z**11)
 # Independent direct action on a dummy function.
 fun=sp.Function('nu')(z);direct=fun+z*sp.diff(fun,z)
 for _ in range(5):direct=-z**2*sp.diff(direct,z)
 check('physical_areal_all_coefficients',sp.expand(direct-sum(cp[j]*sp.diff(fun,z,j) for j in range(7)))==0)
 hj=[1,-1,0,2,-2,-6,16]
 zpt=arb('.5')+e
 cval=[evaluate(c,{z:zpt}) for c in cp]
 pscaled=amplitude*arb('-.5').exp()*e**(-arb(5)/2)*(16*zpt**11-sum((abs(cval[j])*abs(hj[j])*e**(6-j) for j in range(1,6)),arb(0)))
 pphysical=pscaled/B2**(arb(5)/2)
 check('actual_physical_areal_fifth_derivative_exceeds_one',pphysical>1)
 Bmid=16*arb('-.5').exp()*arb('.5')**11/B2**(arb(5)/2)
 # Whole two-row regular-bolt limit; no principal coefficient is assumed.
 lam,gam,tip,v=sp.symbols('lambda gamma tip v',positive=True)
 disk_sub={z:tip-gam*lam*v/2,q[0]:gam**2*lam*v,q[1]:-2*gam}
 boltc=[]
 for aa,fac in enumerate([6,2]):
  row=[]
  for k in [2,3,4]:
   poly=sp.Poly(sp.expand(coeff[aa][k].subs(disk_sub)),lam)
   check('bolt_no_lower_scaling_power_row_'+str(aa)+'_order_'+str(k),all(poly.nth(j)==0 for j in range(max(0,k-2))))
   val=sp.factor(poly.nth(k-2)*(-2/gam)**k)
   expected_limit=16*fac*dd*tip**4*{2:2,3:4*v,4:v**2}[k]
   check('complete_bolt_scaled_row_'+str(aa)+'_order_'+str(k),sp.simplify(val-expected_limit)==0)
   row.append(str(val))
  boltc.append(row)
 # Rational integral in v gives the exact full ordered H3 leading norm.
 p=sp.Poly((v-1)**16*(4-v)**16,v)
 h=sp.Poly(v**2*p.diff((0,4)).as_expr()+4*v*p.diff((0,3)).as_expr()+2*p.diff((0,2)).as_expr(),v)
 h2=h.diff((0,2));h3=h.diff((0,3))
 integrand=sp.Poly(192*v*h2.as_expr()**2+192*v**2*h2.as_expr()*h3.as_expr()+64*v**3*h3.as_expr()**2,v)
 ant=integrand.integrate();Jrat=ant.eval(4)-ant.eval(1)
 check('radial_integral_positive_exact',Jrat>0)
 p6=p.diff((0,6)).eval(sp.Rational(5,2))
 check('radial_sixth_jet_exact',p6==-sp.factorial(6)*sp.binomial(16,3)*sp.Rational(9,4)**13)
 # Tensor identity checked directly at (t,0), including all eight ordered entries.
 xx,yy,t=sp.symbols('x y t',real=True);hh=sp.Function('h');rad=hh(xx**2+yy**2)
 tensor=sum(sp.binomial(3,j)*sp.diff(rad,xx,j,yy,3-j)**2 for j in range(4))
 tensoraxis=sp.simplify(tensor.subs({xx:t,yy:0}))
 target=192*t**2*sp.Subs(sp.Derivative(hh(v),(v,2)),v,t**2)**2+192*t**4*sp.Subs(sp.Derivative(hh(v),(v,2)),v,t**2)*sp.Subs(sp.Derivative(hh(v),(v,3)),v,t**2)+64*t**6*sp.Subs(sp.Derivative(hh(v),(v,3)),v,t**2)**2
 check('radial_ordered_third_tensor',sp.simplify(tensoraxis-target)==0)
 gamma,zh=[dec(bp[k]) for k in ['gamma','zh']]
 J=arb.pi()*evaluate(Jrat,{})
 A3=16*dphys*zh**4*(40*J).sqrt()
 Bbolt=64*zh**11*abs(evaluate(p6,{}))/(B2**(arb(5)/2)*gamma**6)
 check('regular_bolt_positive_finite_coefficients',A3>0 and Bbolt>0 and A3.is_finite() and Bbolt.is_finite())
 # Exact UV coefficient on the actual middle and the physical forward factor.
 uv=C*Qc**2/(15*arb.pi()*B2)
 rem_log10=4000+(2*C/B2).log()/arb(10).log()
 check('complete_quantum_uv_log_nonzero',uv>0)
 check('forward_remainder_log_finite',rem_log10.is_finite())
 physical_result={k:enc(val) for k,val in {'areal_fifth_derivative_lower':pphysical,
  'B5_times_areal_fifth_lower':pscaled,'areal_fifth_to_Ri_ratio_lower':pphysical/Ri,
  'middle_physical_asymptotic_coefficient':Bmid,'middle_physical_ratio_coefficient':Bmid/A}.items()}
 bolt_result={'exact_source_limit_coefficients':boltc,'radial_integral_over_pi':str(Jrat),
  'radial_sixth_jet':str(p6),**{k:enc(val) for k,val in {'A3_source_limit':A3,
  'B_bolt_output_limit':Bbolt,'output_source_ratio_coefficient':Bbolt/A3}.items()},
  'ratio_width_exponent':-6,'source_indices_ruled_out':'s<9','endpoint_s9_decided':False}
 quantum_result={k:enc(val) for k,val in {'actual_middle_log_coefficient_magnitude':uv,
  'log10_normalized_forward_remainder_coefficient':rem_log10}.items()}

 return {'physical_areal_output':physical_result,'regular_bolt_H3_test':bolt_result,'quantum_stronger_source_test':quantum_result,'precision_bits':512,'complete_original_column_coefficients':[[str(c) for c in row] for row in coeff],
  'coefficient_upper_bounds':[[enc(v) for v in row] for row in U],
  'finite_actual_domain_test':{k:enc(v) for k,v in {'width':e,'amplitude':amplitude,
   'scalar_row_L2_upper':row_upper[0],'shape_row_L2_upper':row_upper[1],'complete_interior_Ri_upper':Ri,
   'whole_bolt_error_upper':boltdebit,'actual_proper_fifth_derivative_lower':lower,
   'certified_C5_to_residual_ratio_lower':lower/Ri}.items()},
  'asymptotic_comparison_failure':{k:enc(v) for k,v in {'source_coefficient_A':A,'output_coefficient_B':B,
   'ratio_growth_coefficient_B_over_A':B/A}.items()},
  'exact_observations':{'outer_residual_Ro':0,'actual_finitepart_c4':0,'actual_bolt_response':0},
  'scope':{'no_finite_C5_comparison_from_L2_residual_tuple':True,'no_finite_areal_C5_comparison_from_regular_H3':True,'H14_to_quantum_H10_unbounded':True,'complete_original_X_defined':False,'full_coupled_operator_used':True,
   'finite_example_passes_complete_bolt_error_slice':True,'actual_retained_quantum_center_evaluated':False,
   'original_IY_proved_infinite':False,'F3_output_maps_disproved':False,'original_defect_certified':False,
   'actual_mass_target_certified':False,'complete_coupled_estimate_certified':False,'higher_degree_adopted':False,
   'global_H3_to_C5_inverse_bound_certified':False,'normalization_changed':False,
   'construction_obstruction':False,'new_physical_gate_closed':False}}
if __name__=='__main__':
 mode='acceptance' if '--acceptance' in sys.argv else 'development'
 out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
 hashes={p:sh(ROOT/p) for p in INPUTS};run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
  'input_sha256_before':hashes,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
  'runtime':{'python':sys.version,'python_sha256':sh(Path(sys.executable)),'flint_version':flint.__version__,'sympy_version':sp.__version__}}
 try:
  run['result']=execute();run['input_sha256_after']={p:sh(ROOT/p) for p in INPUTS}
  check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
 except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
 run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
 ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};ledger['runs'].append(run)
 out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
 summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
 with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(summary.encode())
 print(summary,flush=True);sys.exit(0 if run['status']=='passed' else 1)
