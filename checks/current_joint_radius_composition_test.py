"""Complete structured bolt map after scalar/mixed composition; no signed solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='current-joint-radius-composition-test'
OWN=['notes/'+STEM+'.md','checks/current_joint_radius_composition_test.py',
 'evidence/'+STEM+'-sources.json','evidence/'+STEM+'-review.json']
NOTES=['current-admission-bottleneck-estimator-test','current-schur-residual-phase-domain',
 'current-original-admission-map-budget-test','current-native-output-structured-source-test',
 'current-global-mixed-schur-preconditioner','product-centered-scalar-cold-response',
 'current-complete-source-square-factorization','current-interpolated-complete-bolt-test',
 'current-complete-bolt-energy-bound','current-coupled-truncation-energy-error',
 'current-adjoint-uniform-estimator-test','current-controlled-solution-audit-closeout']
RECEIPTS=['current-admission-bottleneck-estimator-test','current-interpolated-complete-bolt-test',
 'current-coupled-truncation-energy-error','current-localized-bolt-adjoint-sensitivity',
 'current-complete-bolt-energy-bound','current-complete-source-square-factorization']
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
 check('frozen_ordinary_sources',all(sha(ROOT/p)==h for p,h in rev['scientific_sha256'].items()))
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes())
 check('pinned_python',sha(Path(sys.executable))==rt['python_sha256'])
 check('pinned_packages',all(sha(Path(pkg['root'])/p)==h for pkg in rt['packages'].values() for p,h in pkg['files'].items()))
 old=prior('current-admission-bottleneck-estimator-test')
 point=prior('current-interpolated-complete-bolt-test')['constants']
 energy=prior('current-coupled-truncation-energy-error')
 local=prior('current-localized-bolt-adjoint-sensitivity')['parameters']
 mix=prior('current-complete-bolt-energy-bound')['constants']
 beta,C,B2,gamma=[dec(energy['parameters'][k]) for k in ['beta','C','B2','gamma']]
 h=dec(local['h']);ellb=(1+h)/h
 q=arb('.338');f=q/(1+q);epsbar=C/(60*arb.pi()*B2)
 d=96*arb.pi()*arb('1e22')*C/B2-epsbar;cvdelta=11*epsbar*beta/3
 CD=dec(energy['global_shape_Laplacian_constant']);Cnew=dec(point['Cnew'])
 check('actual_positive_geometry_and_scalar_gap',h>0 and gamma>0 and beta>arb('.1666') and d>0)
 T0=arb('17.6');gT=arb('.17').sqrt()*T0+(arb('.17')*T0**2+T0*CD).sqrt()
 DT=CD+2*arb('.17').sqrt()*gT;HT=(DT**2+13*beta*gT**2).sqrt()
 PT=DT+2*gT+arb('2.001')*T0;m4=4*beta
 U0=PT/(3*m4);U1=PT/(6*m4.sqrt());U2=PT/3
 w0=arb('1.1e-27')*beta;t1=arb('1.4e-26')*beta**(arb(3)/2);t2=arb('8.6e-25')*beta**2
 t0=w0+22*epsbar*beta**2;tc0=arb('3.01')*beta;tc1=4*tc0
 k0=(2*d+11*epsbar/6)*w0+22*epsbar/3*arb('1e-26')
 k1=(2*d+11*epsbar/6)*t1+22*epsbar/3*arb('1e-25')
 EQ=2*d*(t0*DT+2*t1*gT+t2*T0)+epsbar*arb('1e12')*(HT+gT+T0)+20/(60*arb.pi())*arb('4e-27')*beta*T0
 QQ=(22*epsbar*tc0/3+k0)*U2+(22*epsbar*tc1/3+k1)*U1+2*w0*U0
 R0=EQ+2*cvdelta*PT+QQ
 S0=1/arb('23.99').sqrt();Ds=(beta/(6*d)).sqrt()+2*arb('.17').sqrt()*(beta/arb('5.99')).sqrt()
 rstar=S0*gamma/(4*Ds);Cs0=4*(S0*Ds/gamma).sqrt()
 check('scalar_point_averaging_interval_is_admissible',rstar>0 and rstar<arb('.99'))
 K0=(1+arb('13.5')+(arb('2.001')+4*beta)/(4*beta))/3
 Krem=Cs0*S0*R0/beta
 x=-1-2/q;t_prev,t_now=arb(1),x
 for _ in range(2,30):t_prev,t_now=t_now,2*x*t_now-t_prev
 D29=abs(t_now);h2=(1+q)**2/D29;h3=(1+q)**3/D29;paff=arb('1.1e-10')*h3
 check('unchanged_P32_prefactor',paff.overlaps(dec(old['results']['p_aff'])))
 Kret=dec(old['results']['Kret_original_mixed_point_conversion'])
 AY=Kret*dec(mix['KV'])*dec(mix['Kb'])*q.sqrt()*f**2
 debits={'composed_L4_inverse_P_return':paff*ellb*K0*Cnew*f**3,
         'complete_composition_residual_return':paff*ellb*Krem*f**3,
         'complete_correlated_mixed_return':paff*AY}
 check('same_entire_mixed_return',debits['complete_correlated_mixed_return'].overlaps(dec(old['complete_positive_estimator_terms']['complete_mixed_York_return']['physical_debit'])))
 total=sum(debits.values(),arb(0));CJ=ellb*(K0*Cnew+Krem)*f**3+AY
 oldtotal=dec(old['results']['complete_bolt_upper']);budget=arb('1.667e-9')
 check('complete_physical_composition_identity',(paff*CJ).overlaps(total))
 check('consequential_complete_estimator_gain',oldtotal/total>arb('1e8'))
 check('new_complete_upper_estimator_still_exceeds_target',total>budget)
 check('perfect_removal_of_composed_scalar_return_still_fails',debits['complete_correlated_mixed_return']>budget)
 check('perfect_removal_of_mixed_return_still_fails',debits['composed_L4_inverse_P_return']+debits['complete_composition_residual_return']>budget)
 check('shared_residual_is_paid_and_subdominant',Krem>0 and Krem<K0*Cnew)
 # Preserve factor order; L4^-1 P is a symbol, never P L4^-1.
 D0,P,Q,E,c=S.symbols('D0 P Q E c',commutative=False)
 # Apply A to u0 using L u0=P/3, without manipulating unrelated factors.
 A_u0=2*(D0+c)*P+Q*S.Symbol('u0',commutative=False)
 C_T=2*D0*P+E
 r=E-2*c*P-Q*S.Symbol('u0',commutative=False)
 check('exact_noncommuting_complete_row_cancellation',S.expand(C_T-A_u0-r)==0)
 check('unchanged_limited_quantum_fast_screen',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
 return {'precision_bits':512,'input_space':'actual TT E, ell/N; T3=VF^3u; full paired-domain F',
  'output':'original bolt row L_b=(1+h)nu(L)/h of joint original radius map',
  'constants':{k:enc(v) for k,v in {'beta':beta,'d':d,'epsbar':epsbar,'CD':CD,'Cnew':Cnew,
   'gT':gT,'DT':DT,'HT':HT,'PT':PT,'U0':U0,'U1':U1,'U2':U2,'EQ':EQ,'QQ':QQ,
   'R0_complete_residual':R0,'Cs0_scalar_point':Cs0,'scalar_averaging_length':rstar,
   'K0_composed_point':K0,'Krem_complete_residual_point':Krem,'h2':h2,'h3':h3,'p_aff':paff,
   'attained_joint_bolt_operator_upper':CJ,'attained_joint_physical_upper':total,
   'complete_estimator_improvement':oldtotal/total,'remaining_overrun_factor':total/budget,
   'conditional_target_if_whole_margin_free':budget/paff,'margin_minus_upper':budget-total,
   'after_perfect_scalar_removal_ratio':debits['complete_correlated_mixed_return']/budget,
   'after_perfect_mixed_removal_ratio':(debits['composed_L4_inverse_P_return']+debits['complete_composition_residual_return'])/budget}.items()},
  'complete_debits':{k:enc(v) for k,v in debits.items()},
  'unknowns':{k:None for k in ['actual_signed_joint_map','higher_original_admission_comparisons','original_Y_source_comparison',
   'Cdef3_IY','graph_comparison','actual_kappa_m','remaining_coupled_columns','application_errors','nonlinear_neighborhood']},
  'decision':'Retain exact joint composition gain; reject its positive two-term P32 bolt certificate; require signed transferred residual, not isolated component improvement.',
  'scope':{'complete_structured_bolt_map_bounded':True,'quantitative_joint_gain_executed':True,
   'bolt_allowance_closed':False,'complete_admission_closed':False,'actual_quantum_inverse_numerically_applied':False,
   'actual_signed_correction_computed':False,'full_A_membership_uses_accepted_fixed_seed_graph_theorem':True,
   'numerical_full_A_norm_certified':False,'construction_obstruction':False,'full_goal_complete':False}}
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
