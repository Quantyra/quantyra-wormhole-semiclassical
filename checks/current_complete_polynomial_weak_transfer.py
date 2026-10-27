"""Complete actual-domain P32 approximation transfer; retained signed center is unknown."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='current-complete-polynomial-weak-transfer'
OWN=['notes/'+STEM+'.md','checks/current_complete_polynomial_weak_transfer.py',
 'evidence/'+STEM+'-sources.json','evidence/'+STEM+'-review.json']
NOTES=['current-joint-radius-composition-test','current-coupled-truncation-energy-error',
 'current-interpolated-complete-bolt-test','current-coupled-adjoint-residual-audit',
 'current-adjoint-uniform-estimator-test','current-bolt-adjoint-domain',
 'current-angular-resummed-quantum-operator','current-schur-residual-phase-domain',
 'current-controlled-solution-audit-closeout']
RECEIPTS=['current-joint-radius-composition-test','current-interpolated-complete-bolt-test',
 'current-coupled-truncation-energy-error','current-admission-bottleneck-estimator-test',
 'current-complete-bolt-energy-bound','current-localized-bolt-adjoint-sensitivity']
INPUTS=OWN+['notes/'+n+'.md' for n in NOTES]+['evidence/'+n+'.json' for n in RECEIPTS]+[
 'evidence/current-coupled-mass-adjoint-runtime.json']
CHECKS=[]
def sha(p):return sha256(p.read_bytes()).hexdigest()
def check(n,c):
 ok=bool(c);CHECKS.append({'name':n,'passed':ok});print(n+(': passed' if ok else ': FAILED'),flush=True)
 if not ok:raise AssertionError(n)
def enc(x):
 m,r,e=x.mid_rad_10exp(50)
 return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
 if not isinstance(x,dict):return arb(x.strip('[]'))
 d=x['decimal_enclosure'];return arb(d['mid']+' +/- '+d['rad'])*arb(10)**int(d['exponent'])
def prior(stem):
 r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
 check('accepted_'+stem,r['status']=='passed');return r['result']
def execute():
 ctx.prec=512
 rev=json.loads((ROOT/OWN[-1]).read_bytes())
 check('frozen_ordinary_argument_and_arithmetic',all(sha(ROOT/p)==h for p,h in rev['scientific_sha256'].items()))
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes())
 check('pinned_python',sha(Path(sys.executable))==rt['python_sha256'])
 check('pinned_packages',all(sha(Path(pkg['root'])/p)==h for pkg in rt['packages'].values() for p,h in pkg['files'].items()))
 joint=prior('current-joint-radius-composition-test')['constants']
 point=prior('current-interpolated-complete-bolt-test')
 energy=prior('current-coupled-truncation-energy-error')
 priorrows=prior('current-admission-bottleneck-estimator-test')['results']
 mixed=prior('current-complete-bolt-energy-bound')['constants']
 local=prior('current-localized-bolt-adjoint-sensitivity')['parameters']
 beta=dec(energy['parameters']['beta']);Cnew=dec(joint['Cnew'])
 h=dec(local['h']);ellb=(1+h)/h
 q=arb(169)/500;f=q/(1+q);U=arb('1.1e-10');budget=arb('1.667e-9')
 e1=arb('9.204e-79')*arb('.3961')**128
 e2=arb('9.204e-77')*arb('.3961')**128
 eps=e2*Cnew**2/beta;ef=e1*Cnew/beta
 check('same_actual_form_tail_and_source',eps.overlaps(dec(point['complete_retained_inverse_transport_test']['epsilon_new'])) and ef.overlaps(dec(point['complete_retained_inverse_transport_test']['source_error_energy_dual'])))
 check('same_domain_Neumann_inverse_admissible',eps>0 and eps<arb('.5'))
 inv=1/(1-eps);M=q*inv;dT=eps*q*inv
 def dn(n):return eps*f**n*sum((inv**j for j in range(1,n+1)),arb(0))
 j0=ellb*(dec(joint['K0_composed_point'])*Cnew+dec(joint['Krem_complete_residual_point']))
 j1=dec(priorrows['Kret_original_mixed_point_conversion'])*dec(mixed['KV'])*dec(mixed['Kb'])*q.sqrt()
 CJ=j0*f**3+j1*f**2
 check('both_original_rows_reproduce_accepted_complete_bound',CJ.overlaps(dec(joint['attained_joint_bolt_operator_upper'])))
 DG0=j0*dn(3);DG1=j1*dn(2);DG=DG0+DG1
 # Rational polynomial throughout; no sampling of physical modes.
 t=S.Symbol('t');qr=S.Rational(169,500)
 Q=S.Poly(((1+t)/(1+qr))**3*S.chebyshevt(29,2*t/qr-1),t)
 check('exact_complete_polynomial_degree',Q.degree()==32)
 Dexact=abs(S.chebyshevt(29,-1-2/qr))
 def aa(x):return arb(str(S.numer(x)))/arb(str(S.denom(x)))
 D29=aa(Dexact);h3=(1+q)**3/D29
 check('same_P32_physical_prefactor',h3.overlaps(dec(joint['h3'])))
 check('exact_P32_residual_identity',S.cancel(((1+qr)**3/Dexact)*(t/(1+t))**3*Q.as_expr()-t**3*S.chebyshevt(29,2*t/qr-1)/Dexact)==0)
 coeff=[Q.nth(k) for k in range(33)]
 DQ=sum((k*aa(abs(coeff[k]))*M**(k-1) for k in range(1,33)),arb(0))
 dQ=DQ*dT;dw=(ef+eps*U)*inv
 terms={'joint_map_change':h3*DG*U,
  'complete_polynomial_change':h3*(CJ+DG)*dQ*U,
  'full_affine_response_change':h3*(CJ+DG)*(1+dQ)*dw}
 E32=sum(terms.values(),arb(0));Cstar=budget/(h3*U)
 check('complete_approximation_error_below_original_bolt_margin',E32<budget)
 check('joint_operator_approximation_error_below_conditional_norm_target',DG<Cstar)
 check('signed_retained_center_still_required',h3*U*CJ>budget)
 check('no_false_certification_by_error_transfer',h3*U*(CJ+DG)+E32>budget)
 check('strict_nonzero_complete_errors',all(x>0 for x in [DG0,DG1,dQ,dw,*terms.values()]))
 return {'precision_bits':512,
  'constants':{k:enc(v) for k,v in {'epsilon':eps,'source_dual_error':ef,'q':q,'f':f,'M':M,
   'dT':dT,'dF2':dn(2),'dF3':dn(3),'j0':j0,'j1':j1,'CJ':CJ,
   'DG_scalar_row':DG0,'DG_mixed_row':DG1,'DG_complete':DG,'DQ':DQ,'dQ':dQ,
   'dw':dw,'h3':h3,'E32_complete':E32,'original_bolt_margin':budget,
   'E32_fraction_of_margin':E32/budget,'conditional_center_margin_after_E32':budget-E32,
   'conditional_joint_norm_target':Cstar,'DG_fraction_of_conditional_target':DG/Cstar,
   'positive_estimator_remaining_overrun':h3*U*CJ/budget}.items()},
  'complete_error_debits':{k:enc(v) for k,v in terms.items()},
  'exact_Q32_coefficients_ascending':[str(v) for v in coeff],
  'unchanged_approximation_depth':256,'degree':32,
  'unknowns':{k:None for k in ['retained_signed_joint_norm','retained_signed_affine_correction',
   'retained_quantum_form_application','numerical_solve_and_quadrature_errors',
   'full_original_A_comparison','original_Y_extension','original_defect','actual_kappa_m',
   'all_other_coupled_entries','nonlinear_neighborhood']},
  'scope':{'complete_P32_bolt_transport_difference_certified':True,
   'weak_extension_to_actual_adjoint_uses_density':True,'rough_nonlinear_metric_family_claimed':False,
   'signed_quantum_center_executed':False,'bolt_allowance_closed':False,
   'full_geometric_admission_closed':False,'complete_numerical_application_certified':False,
   'complete_coupled_bound_certified':False,'construction_obstruction':False,'full_goal_complete':False},
  'decision':'Keep actual weak signed method and P32. Existing transport accuracy suffices for complete polynomial bolt transfer. Stop tail work; evaluate the retained signed response.'}
if __name__=='__main__':
 mode='acceptance' if '--acceptance' in sys.argv else 'development'
 output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
 hashes={p:sha(ROOT/p) for p in INPUTS}
 run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
 'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
 'runtime':{'python':sys.version,'executable_sha256':sha(Path(sys.executable))}}
 try:
  run['result']=execute();run['input_sha256_after']={p:sha(ROOT/p) for p in INPUTS}
  check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
 except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
 run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
 ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
 ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
 message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
 with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(message.encode())
 print(message,flush=True)
 if run['status']!='passed':raise SystemExit(1)
