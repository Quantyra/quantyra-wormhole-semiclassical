"""Actual-domain evaluator test; not a signed quantum or coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-bolt-resummed-dirac-parametrix'
OWN=['notes/'+STEM+'.md','checks/current_bolt_resummed_dirac_parametrix.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-global-dirac-transport-series.md',
 'notes/current-adaptive-dirac-trace-remainder.md','notes/current-interpolated-complete-bolt-test.md',
 'evidence/current-interpolated-complete-bolt-test.json','evidence/current-global-mass-conversion.json',
 'evidence/current-localized-bolt-adjoint-sensitivity.json','notes/current-bolt-source-domain-test.md',
 'notes/current-third-metric-conversion.md','notes/current-coincident-dirac-action-assembly.md']
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
  if p in INPUTS and p.startswith('notes/'):check('unchanged_'+p,sh(ROOT/p)==h)
 return r['result']
def execute():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('ordinary_argument_and_checker_frozen',all(sh(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
 latest=prior('current-interpolated-complete-bolt-test')
 p=prior('current-global-mass-conversion')['parameters']
 bp=prior('current-localized-bolt-adjoint-sensitivity')['parameters']
 beta,u,eps,tau=[dec(p[k]) for k in ['beta','u','eps','tau']]
 gamma,C,B2=[dec(bp[k]) for k in ['gamma','C','B2']]
 check('actual_positive_temperature_small',gamma>0 and gamma<arb('1e-33'))
 q=arb(2000)**4
 qp=16*beta/9+u/2+4*abs(eps)+3*tau
 qpp=12*beta+2*u+82*abs(eps)+6*tau
 check('whole_polynomial_Qprime_absolute',qp<arb('.3'))
 check('whole_polynomial_Qsecond_absolute',qpp<arb('2.1'))
 check('whole_polynomial_Q_upper',beta+u+abs(eps)<arb('.17'))
 check('central_Q_lower',beta/2-abs(eps)-tau>arb('.04'))
 check('deep_Qprime_strict_sign_by_dominance',12*beta>4*abs(eps) and u>0 and tau>0)
 n0=1
 while arb(n0)*(n0+q)/C<10000:n0+=1
 check('first_entire_angular_subsequence',arb(n0)*(n0+q)/C>10000 and arb(n0-1)*(n0-1+q)/C<10000)
 Rlo=1-arb('.00042');Rhi=1+arb('.00042')
 klower=1/Rhi.sqrt();kupper=(1+1/(4*Rlo)).sqrt()/Rlo.sqrt()
 alower=(2-arb('.00168'))/(8*Rhi.sqrt()*(16*Rhi+arb('.25')))
 check('entire_bolt_rectangle_k_bounds',klower>arb('.5') and kupper<2)
 check('entire_bolt_rectangle_alpha_lower',alower>arb('.01'))
 positive=2*arb(-12).exp()/10000
 central=16*gamma/20
 cross=128*(-200/arb('.17').sqrt()).exp()
 trace_lower=positive-central-cross
 norm_lower=arb(-6).exp()/200
 check('full_signed_trace_lower_not_only_patch',trace_lower>arb('1e-9'))
 check('old_reference_has_mass_independent_norm_lower',norm_lower>arb('1.2e-5'))
 # Volterra comparison is uniform by positivity of the sinh series.
 rho=arb('.2');kk=arb('1.05').sqrt();ratio=(kk*rho).sinh()/(kk*rho)
 flat_rel=(ratio-1)/(2-ratio)/rho**2
 check('actual_proper_distance_a_comparison',ratio<arb('1.01') and 2-ratio>arb('.99'))
 check('actual_proper_distance_flat_error',flat_rel<arb('.22'))
 d=arb(100);flat_coef=arb('.88')*arb(3).sqrt()
 outer_coef=arb(2).sqrt()/(4*arb('.99'))+arb('8.4').sqrt()/(8*d.sqrt())
 check('flat_high_thermal_local_defect',flat_coef<arb('1.53'))
 check('flat_critical_local_defect',arb('.22')*d.sqrt()/(d-gamma/2)<arb('1.53'))
 check('outer_actual_transport_local_defect',outer_coef<arb('.4'))
 check('flat_local_defect_fits_common_envelope',arb('1.53')/d.sqrt()<arb('.4'))
 check('large_thermal_localized_inverse',2*arb(3).sqrt()*gamma/d.sqrt()<1)
 check('patches_disjoint_in_actual_geometry',4/d.sqrt()<2*arb('.99')/arb('.17').sqrt())
 defect=(arb.pi()+arb('.4'))/d.sqrt()
 check('complete_matter_parametrix_defect',defect<arb('.4'))
 # Trace-ideal factors include the critical |omega|=1/2 branch.
 check('critical_Bessel_small_argument_bound',arb(1).exp()*(arb(2).log()+1)<5)
 check('critical_Bessel_large_argument_bound',arb.pi()/4<1)
 hscoef=arb(2).sqrt()*395/d.sqrt()+arb(24)/18
 check('entire_reference_Hilbert_Schmidt_coefficient',hscoef<58)
 outerV=1/(4*arb('.99'))+arb('4.2').sqrt()/(8*d.sqrt())
 check('outer_zero_order_coefficient',outerV<arb('.3'))
 check('full_defect_Hilbert_Schmidt_row',arb('.44')+gamma*(arb('.3')+arb.pi()/2)<1)
 F5=C/3+4*C*gamma/81
 Kq=58*arb(4)**10*F5/(arb.pi()*B2*gamma)
 Cnew=dec(latest['constants']['Cnew']);Knew=dec(latest['constants']['complete_unfiltered_Knew'])
 betaout=dec(bp['beta']);energy=arb('1.1e-10');margin=arb('8.335e-10')
 old1=arb('9.204e-79')*arb('.3961')**128
 old2=arb('9.204e-77')*arb('.3961')**128
 def evaluate(N):
  e1=Kq*arb('.4')**(N-11)/arb('.6')
  e2=Kq*arb('.4')**(N-12)*(arb('.4')/arb('.6')+2*(arb(N)/arb('.6')+arb('.4')/arb('.6')**2))
  ef=(old1+e1)*Cnew/betaout;eh=(old2+e2)*Cnew*Cnew/betaout
  result={'high_source_tail':e1,'high_Hessian_tail':e2,'total_energy_source_error':ef,'total_form_defect':eh}
  if eh<1:
   result.update({'source_bolt_debit':Knew*ef/(1-eh),'Hessian_bolt_debit':Knew*eh*energy/(1-eh),
    'complete_bolt_upper':Knew*(ef+eh*energy)/(1-eh)})
  return result
 selected=None;previous=None
 for N in range(12,1001):
  r=evaluate(N)
  if 'complete_bolt_upper' in r and r['complete_bolt_upper']<margin:
   selected=(N,r);break
  previous=(N,r)
 check('finite_fixed_word_count_preserves_complete_error_slice',selected is not None)
 N,r=selected
 check('previous_estimator_certifiably_excessive',previous[0]==N-1 and previous[1]['complete_bolt_upper']>margin)
 check('source_and_contact_Hessian_both_retained',r['complete_bolt_upper'].overlaps(r['source_bolt_debit']+r['Hessian_bolt_debit']))
 check('old_transport_certificate_not_replaced_by_zero',old1>0 and old2>0 and r['total_form_defect']>old2*Cnew*Cnew/betaout)
 return {'precision_bits':512,'actual_domain_evaluator_failure':{'first_angular_n_with_M_gt100':n0,
  **{k:enc(v) for k,v in {'full_Tr_S0_squared_lower':trace_lower,'positive_patch_lower':positive,
   'entire_central_negative_upper':central,'entire_cross_bolt_negative_upper':cross,
   'mass_independent_S0_norm_lower':norm_lower}.items()},'separately_unregulated_absolute_word_sum_diverges':True},
  'actual_seed_matter_parametrix':{k:enc(v) for k,v in {'proper_distance_flat_error_coefficient':flat_rel,
   'd100_complete_defect_upper':defect,'output_L2_coefficient_d_times_normP':arb(2),
   'HS_sum_coefficient_upper':hscoef,'complete_spectral_F5_upper':F5,'normalized_Kq':Kq}.items()},
  'complete_transport_slice_test':{'first_passing_fixed_high_sector_N':N,'not_gravitational_polynomial_degree':True,
   'passing':{k:enc(v) for k,v in r.items()},'previous_N':N-1,
   'previous':{k:enc(v) for k,v in previous[1].items()},'full_slice_margin':enc(margin),
   'previous_required_upper_estimator_improvement':enc(previous[1]['complete_bolt_upper']/margin),
   'low_sector_source_error_upper':enc(old1),'low_sector_Hessian_error_upper':enc(old2)},
  'scope':{'specified_retained_word_interchange_fails':True,'actual_seed_matter_parametrix_bounded':True,
   'complete_matter_truncation_to_bolt_slice_passes':True,'actual_retained_quantum_center_evaluated':False,
   'actual_trial_residual_evaluated':False,'numerical_inverse_applied':False,'original_defect_certified':False,
   'graph_physical_output_maps_certified':False,'actual_mass_target_certified':False,
   'every_coupled_entry_certified':False,'complete_coupled_estimate_certified':False,
   'higher_gravitational_degree_adopted':False,'uniform_nonlinear_neighborhood_certified':False,
   'normalization_changed':False,'construction_obstruction':False,'new_physical_gate_closed':False}}
if __name__=='__main__':
 mode='acceptance' if '--acceptance' in sys.argv else 'development'
 out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
 hashes={p:sh(ROOT/p) for p in INPUTS};run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
  'input_sha256_before':hashes,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
  'runtime':{'python':sys.version,'python_sha256':sh(Path(sys.executable)),'flint_version':flint.__version__}}
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
