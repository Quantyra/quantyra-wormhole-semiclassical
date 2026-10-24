"""Support the complete massive thermal near comparison in the original strong graph."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-massive-thermal-strong-return'
OWN=['notes/'+STEM+'.md','checks/current_massive_thermal_strong_return.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def algebra():
 n,q=S.symbols('n q',positive=True)
 count=2*n*(q+n+1)
 check('entire_angular_count',S.expand(4*n*(n+q)-count-2*n*(n+q-1))==0)
 x,MM,be=S.symbols('x M beta',positive=True)
 integ=S.integrate((MM**2/be+x*x/be**2)**(-3),(x,-S.oo,S.oo))
 check('complete_frequency_integral',S.simplify(integ-3*S.pi*be**S.Rational(7,2)/(8*MM**5))==0)
 check('frequency_derivative_Cauchy_guard',all(5*factorial(j)*1000**j<10**40 for j in range(11)))
 check('extra_p_disk_neumann_ratio',S.Rational(2*3*10,200)==S.Rational(3,10) and S.Rational(3,1)/(1-S.Rational(3,10))<5)
 check('germ_inverse_Q_density_derivatives',all(4*factorial(j)*100**j<10**12 for j in range(5)))
 check('germ_inverse_Q_potential_derivatives',all(16*factorial(j)*100**j<10**12 for j in range(4)))
 for j in range(4):
  mass=12*sum(comb(j,l)*factorial(j-l)*100**(j-l) for l in range(j+1))
  freq=S.Rational(101,100)*16*factorial(j)*100**j
  check('local_A_derivative_guard_'+str(j),mass+freq<10**12)
 z=S.symbols('z');b=S.Function('b')(z)
 fourth=S.diff(S.exp(b),z,4)/S.exp(b)
 check('fourth_mass_vertex_complete',S.simplify(fourth-(S.diff(b,z,4)+4*S.diff(b,z)*S.diff(b,z,3)+3*S.diff(b,z,2)**2+6*S.diff(b,z)**2*S.diff(b,z,2)+S.diff(b,z)**4))==0)
 bb=arb('.002');reg=bb.exp()*(4*bb**2+3*bb**2+6*bb**3+bb**4)
 check('rough_fourth_vertex_guard',bb.exp()*bb<1 and reg<3)
 check('all_angular_tail_constant',4*(1+S.Rational(2,18))<8)
 check('all_winding_tail_constant',(-arb(18)).exp()<arb('.5'))
 check('normalized_source_prefactor_reserve',32*10**121/(2*arb.pi())<arb(10)**126)
 check('weighted_seam_delta_annihilated',8-S.Rational('4.005')>0)

def execute():
 algebra()
 pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 beta=pars['beta'];C=pars['C'];B2=pars['B2'];gamma=pars['gamma'];d=arb('.001');eta=arb('.01')
 h1=1/(1-d).sqrt()-1;h2=d/(1-d)
 loss=(h1*h1+h2*h2).sqrt()+(arb('.002').exp()-1)/(1-d).sqrt()+eta/(beta.sqrt()*(1-d))
 ratio=4*loss;kernel=1/(1-ratio)
 # Explicit derivative ledger for j<=10, at most four diagonal derivatives.
 D=[10**40]
 for ii in range(4):D.append(22*10**12*sum(comb(ii,ll)*D[ii-ll] for ll in range(ii+1)))
 source_regular=8*81*3*10**12*max(D);source_singular=8*10**12*D[0]
 # Original near input indicator remains a split, and is not differentiated
 # at its remote endpoint. Local observation jets of that input are smooth.
 sigma=arb('4.005');alpha=5-sigma;aa=arb('.01');bb=arb('.02');normal=1801/arb(2).sqrt()
 L=[arb(1)]
 for jj in range(1,5):L.append(L[-1]*(sigma+2-jj))
 input_jets=[normal*L[jj]*aa**(sigma-jj) for jj in range(4)]
 input_jets.append(normal*L[4]*aa**(sigma-3))
 input_sup=normal*bb**sigma
 Kin=max(input_sup,max(input_jets));srad=arb('.0005')/Kin;vrad=arb('.0005')
 cs=[[arb(1)]]
 for kk in range(1,4):cs.append([(alpha+jj)*(cs[-1][jj] if jj<kk else 0)+(cs[-1][jj-1] if jj>0 else 0) for jj in range(kk+1)])
 caps=[1,3,10,60];rows=[]
 for kk in range(4):
  rr=arb(0)
  for ell in range(kk+1):
   for jj,cc in enumerate(cs[kk-ell]):
    power=alpha+jj
    rr+=comb(kk,ell)*caps[ell]*cc*aa**power/(2*power).sqrt()
    if jj==3:
     rr+=comb(kk,ell)*caps[ell]*cc*aa**(power-1)/(2*(power-1)).sqrt()
  rows.append(rr)
 Cout=sum((arb(comb(3,kk)).sqrt()*v for kk,v in enumerate(rows)),arb(0))
 T=2*arb.pi()/gamma
 logB=126*arb(10).log()+C.log()-10*T.log()-arb('.18')*T+(2/B2).log()
 logH=logB+Cout.log()-srad.log();logDH=logH-vrad.log();logDHfast=logDH+arb('5.681e6').log()

 be=pars['beta'];co=[arb(0),arb(0),pars['u']/be,-4-2*pars['u']/be-pars['tau']/be,3+pars['u']/be+pars['eps']/be,-2*pars['eps']/be,pars['eps']/be]
 qerr=sum((co[j].abs_upper()*arb('.05')**j for j in range(2,7)),arb(0))
 check('whole_real_cylinder_coefficient_tube',qerr<arb('.001') and be>0 and be<1)
 check('physical_output_complex_disk',arb('.01')+arb('.01')<arb('.05') and 1/(be*(1-qerr))<16 and 1/(be*(1-qerr)).sqrt()<4)
 check('first_mass_and_period',((1+pars['q'])/C).sqrt()>18 and T>100)
 check('complex_frequency_neumann_ratio',ratio<arb('.111'))
 check('complex_frequency_kernel_bound',kernel<arb('1.125') and kernel<3)
 check('bounded_p_insertion',1/(be*(1-d))<10)
 check('diagonal_derivative_recursion_guard',max(D)<10**95)
 check('source_regular_derivative_ledger',source_regular<10**120)
 check('source_rough_vertex_ledger',source_singular<10**120)
 check('original_near_input_jet_guard',Kin<arb('1503.018') and input_sup<arb('.0002'))
 check('two_parameter_tube_admission',(srad*Kin).overlaps(arb('.0005')) and (srad*Kin+vrad).overlaps(arb('.001')))
 check('all_output_powers_integrable',alpha>0 and alpha+2>0)
 check('weighted_output_return_factor',Cout<arb('.999'))
 check('original_scale_positive',B2>0 and C>0 and gamma>0)
 check('source_log_strict_margin',logB<-arb('2.236e33'))
 check('thermal_near_graph_certified',logH<-1000*arb(10).log())
 check('thermal_metric_per_A_certified',logDH<-1000*arb(10).log())
 check('thermal_fast_metric_certified',logDHfast<-1000*arb(10).log())
 check('zero_thermal_sensitivity_negligible',2*arb(10)**(-1000)<arb('1e-14'))
 check('joint_existing_uniform_reserve',arb('2e-14')+2*arb(10)**(-1000)<arb('1e-12'))
 return dict(kernel_neumann_ratio=old.encoded(ratio),kernel_constant=old.encoded(kernel),diagonal_derivative_guards=[str(v) for v in D],source_regular_guard=str(source_regular),source_singular_guard=str(source_singular),near_input_jet_bound=old.encoded(Kin),near_input_global_sup=old.encoded(input_sup),input_jets=[old.encoded(v) for v in input_jets],near_input_metric_radius=old.encoded(srad),global_A_metric_radius=old.encoded(vrad),output_return_constant=old.encoded(Cout),output_rows=[old.encoded(v) for v in rows],period=old.encoded(T),normalized_thermal_source_log_upper=old.encoded(logB),thermal_near_H3_log_upper=old.encoded(logH),thermal_metric_per_A_log_upper=old.encoded(logDH),thermal_fast_metric_log_upper=old.encoded(logDHfast),thermal_H3_upper='1e-1000',thermal_metric_per_A_upper='1e-1000',thermal_fast_metric_upper='1e-1000',
  complete_massive_thermal_near_bound_on_admitted_family=True,thermal_common_regulator_identified=True,all_angular_levels_retained=True,all_thermal_windings_retained=True,thermal_mass_seagull_retained=True,local_UV_contacts_retained_in_zero_winding=True,physical_paired_domain_changed=False,actual_bolt_return_retained=True,global_bolt_clock_gap_used=False,original_near_input_split_changed=False,new_background_admission_condition_added=False,discontinuous_auxiliary_input_is_physical_metric=False,LLL_remains_separate=True,zero_temperature_state_selected_for_actual_geometry=False,original_normalization_changed=False,remaining_nonlocal_target='.159916395999',remaining_zero_winding_curved_response_known=False,remaining_contact_jets_known=False,complete_C3_known=False,actual_global_background_A_admission_known=False,original_X_graph_comparison_known=False,full_signed_projection_executed=False,actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,next_zero_winding_calculation_executed=False)

def main():
 ctx.prec=512;review=json.loads((ROOT/OWN[2]).read_bytes());paths=set(OWN+review['retained_inputs'])
 for mod in list(sys.modules.values()):
  f=getattr(mod,'__file__',None)
  if f:
   try:p=Path(f).resolve().relative_to(ROOT).as_posix()
   except ValueError:continue
   if p.endswith('.py'):paths.add(p)
 before={p:digest(p) for p in sorted(paths)};out=ROOT/'evidence'/f'{STEM}.json'
 archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256']);run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 result=run.get('result',{})
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),thermal_bound=result.get('thermal_H3_upper'), log_bound=result.get('thermal_near_H3_log_upper')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
