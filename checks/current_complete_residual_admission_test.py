"""Complete coupled residual/admission counterexample; not a physical solve."""
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
STEM='current-complete-residual-admission-test'
OWN=['notes/'+STEM+'.md','checks/current_complete_residual_admission_test.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['operators/signed_local_core.py','notes/current-finite-local-action-assembly.md',
 'notes/current-bolt-source-domain-test.md','evidence/current-bolt-source-domain-test.json',
 'evidence/current-global-mass-conversion.json','evidence/current-localized-bolt-adjoint-sensitivity.json',
 'notes/current-bolt-resummed-dirac-parametrix.md','notes/current-controlled-solution-audit-closeout.md']
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
 e=arb('1e-50');amplitude=arb('1e-75');ell=arb('.005')
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
 return {'precision_bits':512,'complete_original_column_coefficients':[[str(c) for c in row] for row in coeff],
  'coefficient_upper_bounds':[[enc(v) for v in row] for row in U],
  'finite_actual_domain_test':{k:enc(v) for k,v in {'width':e,'amplitude':amplitude,
   'scalar_row_L2_upper':row_upper[0],'shape_row_L2_upper':row_upper[1],'complete_interior_Ri_upper':Ri,
   'whole_bolt_error_upper':boltdebit,'actual_proper_fifth_derivative_lower':lower,
   'certified_C5_to_residual_ratio_lower':lower/Ri}.items()},
  'asymptotic_comparison_failure':{k:enc(v) for k,v in {'source_coefficient_A':A,'output_coefficient_B':B,
   'ratio_growth_coefficient_B_over_A':B/A}.items()},
  'exact_observations':{'outer_residual_Ro':0,'actual_finitepart_c4':0,'actual_bolt_response':0},
  'scope':{'no_finite_C5_comparison_from_L2_residual_tuple':True,'full_coupled_operator_used':True,
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
