"""Original-bolt contribution of every term; no inverse application."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-admission-bottleneck-estimator-test'
OWN=['notes/'+STEM+'.md','checks/current_admission_bottleneck_estimator_test.py',
     'evidence/'+STEM+'-sources.json','evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/'+n+'.md' for n in ['current-interpolated-complete-bolt-test',
 'current-complete-bolt-energy-bound','current-structured-york-moment-bound',
 'current-original-admission-map-budget-test','current-native-output-structured-source-test',
 'current-adjoint-uniform-estimator-test','current-controlled-solution-audit-closeout',
 'finite-cold-admission-and-response','finite-outer-admission-and-response']]+[
 'evidence/'+n+'.json' for n in ['current-interpolated-complete-bolt-test',
 'current-complete-bolt-energy-bound','current-structured-york-moment-bound',
 'current-localized-bolt-adjoint-sensitivity','current-coupled-mass-adjoint-runtime']]
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
 check('frozen_ordinary_sources',all(sha(ROOT/p)==h for p,h in rev['scientific_sha256'].items()))
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes())
 check('pinned_python',sha(Path(sys.executable))==rt['python_sha256'])
 check('pinned_packages',all(sha(Path(pkg['root'])/p)==h for pkg in rt['packages'].values() for p,h in pkg['files'].items()))
 old=prior('current-complete-bolt-energy-bound')['constants']
 latest=prior('current-interpolated-complete-bolt-test')['constants']
 new=prior('current-structured-york-moment-bound')['constants']
 params=prior('current-localized-bolt-adjoint-sensitivity')['parameters']
 h,gamma=[dec(params[k]) for k in ['h','gamma']]
 q=arb('.338');f=q/(1+q);ellb=(1+h)/h
 Kxi=dec(old['Kxi']);KV=dec(old['KV']);Kb=dec(old['Kb'])
 Cnew=dec(latest['Cnew']);Cs=dec(latest['scalar_point_Cs_unchanged'])
 Kret=ellb*(arb('2.5')*h)**(arb(7)/5)*3000*h**(-arb(7)/5)/gamma.sqrt()
 x=-1-2/q;t0,t1=arb(1),x
 for _ in range(2,30):t0,t1=t1,2*x*t1-t0
 D29=abs(t1);h2=(1+q)**2/D29;h3=(1+q)**3/D29;paff=arb('1.1e-10')*h3
 check('same_P32_and_original_normalization',paff.overlaps(dec(new['p_aff']))
       and (Kret*KV*q.sqrt()).overlaps(dec(latest['J1_unchanged'])))
 AT=ellb*Kxi*Cnew*f**3/2;As=ellb*Cs*q.sqrt()*f**3
 AYpublished=Kret*KV*q.sqrt()*f**2
 AYold=AYpublished*Kb
 AYnew=Kret*dec(new['complete_bar_S2_weighted_bound'])
 check('full_published_three_term_bound_reproduced',(AT+As+AYpublished).overlaps(dec(latest['complete_Cb3_upper'])))
 check('available_Kb_substitution_is_valid_gain',Kb>0 and Kb<1 and AYold<AYpublished)
 check('weighted_outer_source_bound_does_not_improve_bolt_term',AYold<AYnew)
 AY=AYold;budget=arb('1.667e-9')
 debits={'TT_point':paff*AT,'scalar_trace':paff*As,'complete_mixed_York_return':paff*AY}
 total=sum(debits.values(),arb(0));without={}
 for name,value in debits.items():
  check(name+'_positive_estimator_exceeds_whole_bolt_budget',value>budget)
  # Sum the other terms directly, avoiding subtraction of nearly equal balls.
  without[name]=sum((v for k,v in debits.items() if k!=name),arb(0))
  check('perfect_removal_of_'+name+'_still_fails',without[name]>budget)
 joint_target=budget/paff
 application_target=budget/dec(latest['complete_unfiltered_Knew'])
 check('finite_positive_joint_and_application_targets',joint_target>0 and application_target>0)
 row={k:{'physical_debit':enc(v),'ratio_to_whole_bolt_budget':enc(v/budget),
          'fraction_of_complete_positive_estimator':enc(v/total),
          'complete_debit_after_perfect_removal':enc(without[k]),
          'remaining_ratio_after_perfect_removal':enc(without[k]/budget)} for k,v in debits.items()}
 return {'precision_bits':512,'original_output':'L_b(nu,T)=(1+h)nu(L)/h',
  'input':'actual fixed-seed E with ell/N normalization; complete R_S V F^3',
  'complete_positive_estimator_terms':row,
  'results':{k:enc(v) for k,v in {'q':q,'f':f,'h2':h2,'h3':h3,'p_aff':paff,
   'Kret_original_mixed_point_conversion':Kret,'complete_bolt_upper':total,
   'complete_estimator_required_gain':total/budget,
   'joint_Cb3_target_if_whole_margin_free':joint_target,
   'new_weighted_source_bolt_debit':paff*AYnew,
   'weighted_source_bound_ratio_to_existing_bolt_term':AYnew/AY,
   'application_Etot_star_target_if_whole_margin_free':application_target,
   'zero_error_original_defect_product_ceiling':arb('.5')/h3,
   'full_margin_graph_ceiling':arb('5.34235e-12')/paff,
   'limited_two_fast_screen':arb('.028913')+arb('1e-11')*arb('1.0706e10')}.items()},
  'application_source_space':'E_tot* dual to [a(s,s)+ell(xi,xi)]/N; distinct from the TT input E*',
  'unknowns':{k:None for k in ['joint_original_geometric_output_norm','complete_structured_H3',
   'source_test_to_original_Y','Cdef3','I_Y','Cg3','actual_kappa_m','independent_phase_columns',
   'signed_coupled_centers','application_errors','nonlinear_neighborhood']},
  'method_decision':'Stop separated three-term P32 admission estimator and York-only optimization; retain construction/P32 and target joint original radius map.',
  'scope':{'specified_estimator_failure_verified':True,'York_complete_budget_contribution_executed':True,
   'observation_space_and_bolt_domain_specified':True,'complete_nonlinear_XY_pair_constructed':False,
   'actual_operator_norm_lower_bound':False,'actual_signed_correction_computed':False,
   'complete_admission_certified':False,'new_physical_gate_closed':False,
   'construction_obstruction':False,'QCD_theorem_imported':False,'full_goal_complete':False}}
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
