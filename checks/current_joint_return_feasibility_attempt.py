"""Support the near-return argument; no complete strong certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import sys,json,traceback
DEVELOPMENT='--development' in sys.argv
if DEVELOPMENT:sys.argv.remove('--development')
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-joint-return-feasibility-attempt'
OWN=['notes/'+STEM+'.md','checks/current_joint_return_feasibility_attempt.py','evidence/'+STEM+'-review.json']
CHECKS=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(n,b):
 CHECKS.append(dict(name=n,passed=bool(b)))
 if not b:raise AssertionError(n)
def calc():
 ctx.prec=512
 pars,_,_,_=old.scalar.prev.york.actual_parameters()
 sig=arb('4.005');H=arb('.005006');C=pars['C'];beta=pars['beta'];B2=pars['B2'];g=pars['g']
 q=[arb(1),arb(0),pars['u']/beta,(-4*beta-2*pars['u']-pars['tau'])/beta,(3*beta+pars['u']+pars['eps'])/beta,-2*pars['eps']/beta,pars['eps']/beta]
 K=C*beta**2/(15*arb.pi()*B2);L=10*arb.pi()*(pars['Fl']+23*g/(180*arb.pi()))/g
 Cb=[2/(K*L*arb('.005')*sig*(sig-3)),2/(K*L*sig*(sig-3)),2/(K*L*sig),2/(K*L)]
 Cb4=max(Cb[0],8/(K*L),1/(K*(L*arb(6).log()).sqrt()))
 Cq=arb('2.328');W=Cb+[Cb[1],Cb[2],Cb[3],Cb4,Cq/(K*arb(2).log())]
 WL=[max(W[j]*(2+arb('1e18')).log(),Cq/K*arb('1e18')**(j-8)) for j in range(8)]+[Cq/K]
 d=arb(5).sqrt()/2;cm=arb('.99')*arb(2).sqrt()*d;cp=arb('.99')*arb(10).sqrt()*d
 pref=8*C/(arb.pi()*B2)*(1+pars['gamma']*cp*H)
 def Z(r,j):return sum((comb(j,l)*factorial(r+j-l+2)*factorial(l)/cm**(r+j-l+2)/cp**(l+1) for l in range(j+1)),arb(0))
 Kj=[pref*2**21*factorial(j)*400**j*Z(2,j) for j in range(5)]
 Sj=[pref*2**11*factorial(j)*400**j*Z(1,j) for j in range(5)]
 PN=[{(0,0):arb(1)}]
 for n in range(4):
  out={}
  for (p,m),v in PN[-1].items():
   for j,qj in enumerate(q):
    if qj==0:continue
    out[p+j,m]=out.get((p+j,m),arb(0))+qj*(sig-n+p)*v
    out[p+j,m+1]=out.get((p+j,m+1),arb(0))+qj*v
  PN.append(out)
 V=[]
 for n in range(5):V.append(sum((abs(v)*arb('.011')**p*arb('1.001')**(abs(sig-n+p)+1)*max(arb(1),2**(sig-n))*W[m] for (p,m),v in PN[n].items()),arb(0))*arb('2.02').sqrt())
 vi=arb('2.5').sqrt()*arb('3.01').sqrt()*arb(3).sqrt()*arb('3.004')**sig*W[0]
 parts=[]
 for n in range(5):
  row=[Kj[n]*vi,arb(0),arb(0),arb(0)]
  for r in range(1,n):
   for ell in range(n-r):
    mm=n-r-1-ell
    row[1]+=comb(n-r-1,ell)*arb(2*(r+ell+1))/(r+1)*Kj[r+ell]*V[mm]
  for ell in range(n+1):row[2]+=comb(n,ell)*Sj[ell]*V[n-ell]
  if n==4:row[3]=arb('.004')*H**3*(Kj[0]*vi+Sj[0]*V[0])
  parts.append(row)
 cs=[[arb(1)]]
 for ell in range(1,4):cs.append([(5+j-sig)*(cs[-1][j] if j<ell else 0)+(cs[-1][j-1] if j else 0) for j in range(ell+1)])
 by=[arb(0)]*4
 for ell in range(4):
  for r in range(ell+2):
   co=arb('2.02') if r==0 else arb('4.04')*cs[ell][r-1]
   outer=100*arb(comb(3,ell)).sqrt()*co*max(arb(1),2**(4+r-sig))
   for k in range(4):by[k]+=outer*parts[r][k]
 eta=sum(by,arb(0));du=arb('.001')/arb('.499')
 rate=2*(1-(arb('.002').exp()-arb('-.002').exp()))-2*arb('.002').exp()*(arb('.001').exp()-1)-8*du/arb(5).sqrt()
 check('complex_rate_whole_family',rate>arb('1.98'))
 check('two_interface_germ_displacement',du+arb('.09')<arb('.1'))
 check('germ_difference_and_two_arm_inflation',arb('.09').exp()<arb('1.1') and arb('.18').exp()<arb('1.2'))
 check('full_mass_germ_value_and_lipschitz',16*arb('2.21')<36 and 8*arb('.004').exp()*arb('1.1')<36 and 32*arb('1.1')<36)
 check('whole_two_interface_bubble_guard',2*arb('.008').exp()*16*4*((400*5+2*2*10*5)+(400*20+2*2*40*20))<2**21)
 check('whole_two_interface_seagull_guard',2*arb('.004').exp()*4*4*(9+36)<2**11)
 check('enlarged_germ_arms',32/(arb('1.98')*arb(5).sqrt())+arb('.75')<10 and 160/(arb('1.98')*arb(5).sqrt())+arb('1.8')<40)
 check('entire_cylinder_seed_mass_bound',sum((abs(q[k])*arb('.05')**k for k in range(2,7)),arb(0))<arb('.001'))
 check('fixed_seed_germ_inside_one_sided_polynomial',2*H*arb('1.001')+2/(400*arb(18))<arb('.021'))
 check('fixed_seed_inverse_Rouche_guard',sum((abs(q[k])*arb('.021')**k for k in range(2,7)),arb(0))<arb('.001') and 2*(1-arb('.001'))>1)
 check('mass_Taylor_germ_guard',arb('.001')+arb('.002')*sum((arb(1)/arb(7200)**j/factorial(j) for j in range(4)),arb(0))<arb('.004'))
 check('local_input_overlap_three',3/arb('.5')<8 and 3/arb('.5')>4)
 check('fixed_output_input_coordinate_guards',2*H*arb('1.001')<arb('.011') and 3*arb('1.001')<arb('3.004') and 3*arb('1.001')/arb('.999')<arb('3.01'))
 check('full_annular_near_return_allowance',eta<arb('.01755'))
 check('fourth_mass_vertex_retained',by[3]>0 and by[3]<arb('1e-20'))
 check('near_metric_variation_fast_allowance',arb('.01755')*arb('5.681e9')<99710000)
 check('point_evaluation_estimator_fails_to_sum',(sig-8+(8-sig)).contains(0) and 2**0==1)
 budget=arb('.1478262320966')-arb('.01755')-arb('1e-11')*99710000
 check('remaining_original_budget',budget.overlaps(arb('.1292791320966')))
 check('remaining_example_background_budget',(budget-arb('8.123e-37')*arb('1e33')).overlaps(arb('.1284668320966')))
 check('original_seed_and_third_ledger',(arb('.0032800008424')+arb('.01755')).overlaps(arb('.0208300008424')) and 881016306+99710000==980726306)
 return dict(outcome='unresolved_complete_strong_estimate',actual_complete_strong_bound_executed=False,
  near_return_allowance=old.encoded(eta),near_return_parts=[old.encoded(x) for x in by],scaled_bubble_guards=[old.encoded(x) for x in Kj],scaled_seagull_guards=[old.encoded(x) for x in Sj],
  rate=old.encoded(rate),near_return_upper='.01755',near_fast_metric_variation_upper='99710000',
  near_allowance_is_direct_sum_of_jets_not_discontinuous_global_H3_operator=True,
  local_L2_overlap_replaces_divergent_point_evaluation_estimator=True,
  all_angular_thermal_modes_jumps_seagull_retained=True,relative_common_regulator_identified_on_fixed_annuli=True,
  uniform_on_stated_admitted_family=True,original_X_family_admission=False,
  remaining_budget=old.encoded(budget),remaining_with_example_B_cap=old.encoded(budget-arb('8.123e-4')),
  equal_example_eta=old.encoded((budget-arb('8.123e-4'))/2),equal_example_C=old.encoded((budget-arb('8.123e-4'))/arb('2e-11')),
  optimistic_RX_alpha_ceiling=old.encoded(arb('.25')/arb('.0208300008424')),
  optimistic_RX_alpha_squared_ceiling=old.encoded(arb('5e7')/980726306),
  example_B_cap_is_target_not_measurement=True,exploratory_coordinate_or_far_values_credited=False,
  full_eta_dd_known=False,full_C_dd_known=False,B_cap_known=False,original_X_graph_comparisons_known=False,
  original_target_certified=False,complete_upper_bound_known_to_fail=False,construction_obstruction=False,
  actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,
  physical_domain_changed=False,normalization_changed=False,projection_reevaluated=False,audit_status='closed',goal_status='active',
  next_complete_logarithmic_return_calculation_executed=False)

def main():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 prev=json.loads((ROOT/'evidence/current-restored-frozen-response.json').read_bytes())['runs'][-1]
 paths=set(OWN+review['retained_inputs'])|set(prev['input_sha256_after'])
 before={p:digest(p) for p in sorted(paths)}
 out=ROOT/'evidence'/f'{STEM}{"-development" if DEVELOPMENT else ""}.json'
 archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
  frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
  check('previous_scientific_sources_preserved',all(digest(p)==h for p,h in prev['input_sha256_after'].items()))
  run['result']=calc();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=CHECKS;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 if not DEVELOPMENT:
  with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 result=run.get('result',{})
 print(json.dumps(dict(status=run['status'],checks=len(CHECKS),inputs=len(before),exception=run.get('exception'),near_return=result.get('near_return_upper'),near_metric=result.get('near_fast_metric_variation_upper'),remaining=result.get('remaining_budget')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
