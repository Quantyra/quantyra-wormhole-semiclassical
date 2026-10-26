"""Complete residual-to-bolt comparison; no quantum inverse application."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json,sys,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1];STEM='current-coupled-application-residual-bound'
OWN=['notes/'+STEM+'.md','checks/current_coupled_application_residual_bound.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-complete-bolt-energy-bound.md','evidence/current-complete-bolt-energy-bound.json',
 'notes/current-coupled-truncation-energy-error.md','notes/current-strong-quotient-reconstruction-test.md',
 'notes/current-coupled-schur-graph-attempt.md','notes/current-physical-york-dual-assembly.md',
 'notes/current-global-boundary-moment-return.md','notes/current-implicit-phase-graph-domain.md',
 'notes/current-finite-local-action-assembly.md','evidence/current-localized-bolt-adjoint-sensitivity.json']
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
 check('accepted_'+stem,r['status']=='passed');return r
def execute():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('ordinary_argument_checker_frozen',all(sh(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
 prev=prior('current-complete-bolt-energy-bound');pars=prior('current-localized-bolt-adjoint-sensitivity')['result']['parameters']
 check('accepted_endpoint_proof_unchanged',sh(ROOT/'notes/current-complete-bolt-energy-bound.md')==prev['input_sha256_before']['notes/current-complete-bolt-energy-bound.md'])
 A,C,D,nu,T,rn,rt=sp.symbols('A C D nu T rn rt',nonzero=True)
 check('both_forcing_rows_retained',sp.expand((-A*nu-C*T)-rn).subs(nu,(-C*T-rn)/A)==0 and
       sp.expand((-C*nu+D*T)-rt).subs(nu,(D*T-rt)/C)==0)
 z,sigma,c4=sp.symbols('z sigma c4',positive=True);f=sp.Function('f')(z);k=sp.Function('k')(z);g=sp.Function('g')(z)
 rnu=z**sigma*f;rT=c4*z**4+z**sigma*(f/3+k)
 der=lambda a:z*sp.diff(a,z)
 G=sp.expand(der(rT)-4*rT-(der(rnu)-6*rnu)/3)
 check('exact_global_undifferentiated_row_recovery',sp.simplify(G.subs(sp.diff(k,z),(g-2*f/3-(sigma-4)*k)/z)-z**sigma*g)==0)
 check('quartic_kernel_retained',sp.expand(der(z**4)-4*z**4)==0)
 delta=sp.Rational(1,200);w=sp.symbols('w',nonnegative=True)
 ratio=((delta-2)**2+w)/(9*(delta**2+w))
 check('whole_frequency_scalar_recovery_bound',sp.simplify(sp.diff(ratio,w))<0 and ratio.subs(w,0)==133**2)
 h,gamma,beta,Cgeom,B2=[dec(pars[x]) for x in ['h','gamma','beta','C','B2']]
 Fl=1-384*arb.pi()*arb('1e22')*Cgeom*beta/B2
 check('actual_scalar_and_TT_positive_margins',Fl>arb('.99999998') and beta>arb('.16'))
 old=prev['result']['constants'];Cs=dec(old['scalar_point_Cs']);Kxi=dec(old['Kxi']);Cinf=arb('7.477e50');ellb=(1+h)/h
 Kapp=ellb*(Cs**2+(Kxi*Cinf/2)**2).sqrt();budget=arb('1.667e-9')/2
 Ecap=budget/Kapp
 tr=arb(1+133**2+200**2);det=arb(200**2)
 Krow=((tr+(tr**2-4*det).sqrt())/2).sqrt()
 check('complete_fast_source_recovery_guard',Krow<241)
 MR=((16/arb('23.99')+1225/(arb('.99')*Fl))/(4*(3-arb(5).sqrt()))).sqrt()
 cr=MR/beta;cf=Kapp*cr*Krow;cm=Kapp*cr/arb(5).sqrt()
 check('source_to_energy_physical_factor_retained',cr>100 and cr<130)
 fastcap=budget/cf;momentcap=budget/cm
 eH=arb('1.020e-25');ef=arb('1.365e-78');energy=arb('1.1e-10')
 err=(ef+eH*energy)/(1-eH)
 source_debit=Kapp*ef/(1-eH);H_debit=Kapp*eH*energy/(1-eH);debit=Kapp*err
 check('complete_affine_error_decomposition',debit.overlaps(source_debit+H_debit))
 check('unconditional_application_upper_estimator_excessive',debit>budget)
 check('source_only_upper_estimator_excessive',source_debit>budget)
 check('hessian_only_upper_estimator_excessive',H_debit>budget)
 efcap=Ecap # Hessian approximation exact, source only.
 eHcap=Ecap/(energy+Ecap) # source exact, solve eps*E/(1-eps)<=Ecap.
 check('hessian_only_threshold_identity',(Kapp*eHcap*energy/(1-eHcap)).overlaps(budget))
 check('finitepart_source_has_positive_norm',(1+h)**(-arb(5)/2)/arb(5).sqrt()>arb('.44'))
 p76=prev['result']['polynomial_test']['rows']['76'];tail76=dec(p76['bolt_tail_upper'])
 check('degree76_tail_pass_not_application_pass',tail76<budget and debit>budget)
 return {'precision_bits':512,'source_and_output_maps':{key:enc(val) for key,val in
  {'Kapp_unfiltered_complete_coupled_error':Kapp,'Krow_global_row_recovery':Krow,
   'canonical_return_MR':MR,'physical_source_to_energy_MR_over_beta':cr,
   'fast_residual_to_bolt_coefficient':cf,'finitepart_residual_to_bolt_coefficient':cm}.items()},
  'actual_retained_operator_application_estimator':{key:enc(val) for key,val in
   {'energy_error_upper':err,'whole_bolt_error_upper':debit,'source_only_debit':source_debit,
    'hessian_only_debit':H_debit,'whole_improvement_required':debit/budget,
    'source_improvement_required':source_debit/budget,'hessian_improvement_required':H_debit/budget}.items()},
  'conditional_residual_targets':{key:enc(val) for key,val in
   {'half_bolt_margin':budget,'complete_energy_dual_residual_ceiling':Ecap,
    'global_fast_pair_ceiling_moment_zero':fastcap,'finitepart_ceiling_fast_zero':momentcap,
    'source_error_ceiling_hessian_exact':efcap,'hessian_error_ceiling_source_exact':eHcap}.items()},
  'scope':{'full_coupled_unfiltered_residual_to_bolt_bounded':True,
   'original_two_rows_and_finitepart_retained':True,'full_application_upper_estimator_fails':True,
   'original_Y_identified_with_global_residual_test':False,'original_IY_comparison_certified':False,
   'actual_strong_residual_evaluated':False,'actual_retained_inverse_applied':False,
   'original_defect_certified':False,'graph_physical_application_maps_certified':False,
   'actual_mass_column_certified':False,'complete_coupled_estimate_certified':False,
   'higher_degree_adopted':False,'construction_obstruction':False,'new_physical_gate_closed':False}}
if __name__=='__main__':
 mode='acceptance' if '--acceptance' in sys.argv else 'development';out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
 hashes={p:sh(ROOT/p) for p in INPUTS};run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
 'input_sha256_before':hashes,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
 'runtime':{'python':sys.version,'python_sha256':sh(Path(sys.executable)),'flint_version':flint.__version__,'sympy_version':sp.__version__}}
 try:
  run['result']=execute();run['input_sha256_after']={p:sh(ROOT/p) for p in INPUTS};check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
 except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
 run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
 ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};ledger['runs'].append(run)
 out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode());summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
 with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(summary.encode())
 print(summary,flush=True);sys.exit(0 if run['status']=='passed' else 1)
