"""Support the ordinary localized full-response proof; no global inverse solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json,sys,traceback,math
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-interior-translation-jet-bound'
OWN=['notes/'+STEM+'.md','checks/current_interior_translation_jet_bound.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-schur-polynomial-return-test.md','evidence/current-schur-polynomial-return-test.json',
 'notes/current-seed-full-hessian-realization.md','notes/strong-interior-quantum-remainder.md',
 'notes/compact-massive-logarithmic-graph.md','notes/general-interior-dirac-response.md',
 'notes/divided-jet-dirac-response.md','notes/finite-middle-source-four.md',
 'notes/stationary-curved-quantum-tt.md','notes/current-coupled-schur-graph-attempt.md',
 'evidence/current-global-mass-conversion.json','evidence/current-bolt-frequency-domain-test.json']
CHECKS=[]
def check(name,ok):
 CHECKS.append({'name':name,'passed':bool(ok)});print(name+(': passed' if ok else ': FAILED'),flush=True)
 if not ok:raise AssertionError(name)
def enc(x):
 m,r,e=x.mid_rad_10exp(50)
 return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
 d=x['decimal_enclosure'];return arb(d['mid'],d['rad'])*arb(10)**int(d['exponent'])
def execute():
 ctx.prec=512
 rv=json.loads((ROOT/OWN[2]).read_bytes())
 for p,h in rv['scientific_sha256'].items():check('frozen_'+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)
 check('prior_full_polynomial_evidence',json.loads((ROOT/INPUTS[4]).read_bytes())['runs'][-1]['status']=='passed')
 oldp=json.loads((ROOT/INPUTS[-2]).read_bytes())['runs'][-1]
 oldg=json.loads((ROOT/INPUTS[-1]).read_bytes())['runs'][-1]
 check('retained_actual_seed_parameters_passed',oldp['status']==oldg['status']=='passed')
 check('actual_seed_gap_and_conversion_guards',dec(oldp['result']['parameters']['Fl'])>arb('.99') and dec(oldg['result']['numerics']['gamma'])>arb('5e-34') and dec(oldg['result']['numerics']['M1'])>18)
 fac=math.factorial
 nrm=sp.factorial(65)/sp.factorial(32)**2
 check('cutoff_normalization_guard',nrm<=65*2**64)
 check('far_gluing_cutoff_gradient',nrm/4**32<7 and 7*2048<10**5)
 bound=arb(str(nrm))*arb('1.001')**33*arb('2.001')**32/33+1
 check('complex_step_polynomial_guard',bound<arb('1e32'))
 check('finite_cutoff_regularities',32>=17 and 65==2*33-1)
 check('cap_seventeen_jet_guard',all(arb(fac(j))*arb('1e32')*arb('1e6')**j<arb('1e150') for j in range(18)))
 check('localization_seventeen_jet_guard',all(arb(fac(j))*arb('1e32')*arb('1e9')**j<arb('1e201') for j in range(18)))
 rho=arb('1e-250')
 tube=sum(arb('1e201')*rho**j/arb(fac(j)) for j in range(1,13))
 check('finite_parameter_C5_tube',tube<arb('1e-40'))
 budgets=[12*128,12*128**2+900*128,12*128**3+900*128**2+2*10**6*128,12*128**4+900*128**3+2*10**6*128**2+10**8*128]
 check('all_outer_Taylor_divisions',all(a<b for a,b in zip(budgets,[10**4,2*10**6,10**10,10**14])))
 check('complete_Born_inflation',100**4<100**5 and 5**5<100**5)
 check('raw_complex_source_reserve',(arb(2)/3+1)*arb('1e17')*arb(100)**5<arb('1e28'))
 cap_comm=arb(fac(12))*rho**(-12)*arb('1e150')
 check('operator_Cauchy_twelfth_bound',cap_comm<arb('1e3160'))
 check('actual_gluing_separation',sp.Rational(1,64)-sp.Rational(1,128)-sp.Rational(1,4096)>sp.Rational(1,176))
 coerc=arb('-.054').exp()-arb(2)/18
 check('endpoint_coercivity',coerc>arb('.8'))
 worst=arb(0)
 for j in range(1,13):
  mm=j+3; LL=8192*mm; BB=20000*mm
  for k in range(1,j+2):
   ss=sum(arb(math.comb(k,b))*arb(2 if b==1 else 2*fac(b)*100**b)/arb(BB)**b for b in range(1,k+1))
   ratio=arb('1.3')*(arb(LL)/BB+ss)
   if ratio>worst:worst=ratio
  check('endpoint_recurrence_order_'+str(j),all(arb('1.3')*(arb(LL)/BB+sum(arb(math.comb(k,b))*arb(2 if b==1 else 2*fac(b)*100**b)/arb(BB)**b for b in range(1,k+1)))<1 for k in range(1,j+2)))
 check('endpoint_row_constants',all(arb('1e14')*(j+3)*arb(20000*(j+3))**j<arb('1e90') for j in range(1,13)))
 ratio=arb('.0612')*4*arb.pi()
 check('all_Born_kernel_ratio',ratio<arb('.8'))
 env=arb('.6')*(1+arb('.8'))/(1-arb('.8'))**3
 check('whole_kernel_convolution_envelope',env<200)
 check('complete_return_angular_sum_guard',8*200**2<10**6)
 check('full_return_twelfth_budget',3320+10+300+6+50<4000)
 check('full_cap_and_return_fit_twelfth_budget',arb('1e3160')+arb('1e3686')<arb('1e4000'))
 check('remote_input_proper_separation',arb('-.027').exp()/1024>arb(1)/2048)
 check('far_row_rescaling_guard',arb('9.8304e8')*256**4<arb('1e20') and 530*256<10**6)
 check('far_derivative_row_guard',all(arb('1e35')*(j+3)*arb(10**6*(j+3))**j<arb('1e120') for j in range(11)))
 check('far_nested_recurrence_guard',all(arb('1.3')*(arb(131072*(j+3))/(10**6*(j+3))+sum(arb(math.comb(k,b))*arb(2 if b==1 else 2*fac(b)*100**b)/arb(10**6*(j+3))**b for b in range(1,k+1)))<1 for j in range(1,11) for k in range(1,j+2)))
 check('remote_angular_and_H10_budget',8*(10**6//2)**2==2*10**12 and 240+50+13+20<500)
 m0=4*(3-arb(5).sqrt())
 value_embed=(2450*2*arb.pi()/(m0*arb('.99')**2*arb('5e-34'))).sqrt()
 check('original_energy_to_global_value_norm',m0>3 and value_embed<arb('1e19'))
 mat=sp.Matrix([[sp.Rational(k)**j/sp.factorial(j) for j in range(12)] for k in range(12)])
 check('finite_interpolation_determinant',mat.det()==1)
 inv=mat.inv(); rowmax=max(sum(abs(v) for v in inv.row(j)) for j in range(12))
 check('finite_interpolation_inverse_guard',rowmax<10**140)
 check('finite_Taylor_remainder_guard',sp.Rational(11)**12/sp.factorial(12)<10**5)
 check('finite_interpolation_constant',12*10**140*(1+10**5)<10**150)
 gap=sp.Rational(1)-sp.Rational(8-sp.Rational(13,4),12)
 check('strict_derivative_gain',gap==sp.Rational(29,48))
 ds=1/(1-arb(2)**(-arb(29)/24))
 check('infinite_dyadic_square_sum',ds<2)
 check('all_output_orders_fit_declared_constant',150+300+sp.Rational(3700*10,12)+50<4000)
 check('last_fourth_order_bootstrap_reserve',14-sp.Rational(53,4)==sp.Rational(3,4))
 logK=arb(4000).log()/arb(10).log()+4000
 logL=arb(4)/3*(logK+arb(2).log()/arb(10).log())
 loglow=logK+arb(53)/4*logL
 check('absorption_half_exact',sp.Rational(4,3)*sp.Rational(3,4)==1)
 return {'precision_bits':512,'localized_operator':'psi (H_mass-P_sigma) chi; actual paired finite-domain kernel',
  'base_H13over4_to_L2_bound':'10^300 C','twelfth_commutator_H8_to_L2_bound':'10^4000 C',
  'strong_remainder_bound':'||A u||H^s <= 10^4000 C ||u||H^(s+13/4), integer 0<=s<=10',
  'derivative_loss':'13/4','fourth_order_reserve':'3/4','finite_parameter_radius':'1e-250',
  'maximum_endpoint_recurrence_ratio':enc(worst),'kernel_Neumann_ratio_upper':enc(ratio),
  'remote_input_H10_bound':'10^500 C times global gbar L2 norm',
  'original_energy_to_global_gbar_L2_upper':enc(value_embed),
  'remote_input_H10_from_original_energy_bound':'10^519 C E(xi), input T=V xi',
  'finite_interpolation_inverse_row_norm_exact':str(rowmax),'dyadic_square_sum_upper':enc(ds),
  'sensitivity_with_principal_factor_4000_over_C':{'log10_K':enc(logK),'log10_frequency_split_L':enc(logL),'log10_low_frequency_coefficient_s10':enc(loglow),'is_complete_coupled_output_bound':False},
  'scope':{'actual_paired_domain_retained':True,'common_regulator_full_angular_thermal_sums':True,
   'density_contacts_retained':True,'localized_input_finite_domain_returns_included':True,
   'interior_differentiated_massive_response_bounded':True,'global_remote_input_map_bounded':True,
   'global_coupled_output_constants_bounded':False,'original_strong_delta_certified':False,
   'physical_inverse_numerically_applied':False,'actual_kappa_m_computed':False,
   'nonlinear_neighborhood_transfer':False,'nonlinear_admission_contraction':False,
   'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
   'all_sector_stability':False,'actual_controlled_solution':False,'construction_obstruction':False,
   'full_goal_complete':False}}
if __name__=='__main__':
 mode='acceptance' if '--acceptance' in sys.argv else 'development'
 output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
 hashes={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
 run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
  'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
  'runtime':{'python':sys.version,'python_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),'flint_version':flint.__version__,'sympy_version':sp.__version__}}
 try:
  run['result']=execute();run['input_sha256_after']={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
  check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
 except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
 run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
 ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]};ledger['runs'].append(run)
 output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
 summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
 with (ROOT/'evidence'/(STEM+'.log')).open('ab') as out:out.write(summary.encode())
 print(summary,flush=True)
 if run['status']!='passed':raise SystemExit(1)
