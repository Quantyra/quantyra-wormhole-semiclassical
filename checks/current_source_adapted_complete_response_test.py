"""Complete fixed-trial comparison; explicit global kernel rows, no signed center."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback,argparse
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
ROOT=Path(__file__).resolve().parents[1]
STEM='current-source-adapted-complete-response-test'
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS=['notes/'+STEM+'.md','checks/current_source_adapted_complete_response_test.py',
 'checks/current_quartic_ritz_rank_test.py',
 'notes/stationary-curved-quantum-tt.md','notes/current-geodesic-two-sided-remainder.md',
 'notes/current-source-adapted-fractional-return.md','notes/current-source-adapted-local-return.md',
 'notes/current-positive-time-quantum-operator.md','notes/current-first-quantum-solve-attempt.md',
 'notes/current-coupled-truncation-energy-error.md',
 'evidence/current-source-adapted-fractional-return.json',
 'evidence/current-source-adapted-fractional-return-verification.json',
 'evidence/current-source-adapted-local-return-attempt4.json',
 'evidence/current-source-adapted-geodesic-first-attempt.json',
 'evidence/current-positive-time-quantum-operator.json',
 'evidence/current-complete-polynomial-weak-transfer.json',
 'evidence/current-global-green-trial-enclosure.json',
 'evidence/current-global-mass-conversion.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']
def digest(p):return sha256(p.read_bytes()).hexdigest()
def result(name):
 d=json.loads((ROOT/'evidence'/(name+'.json')).read_bytes())
 return d.get('runs',[d])[-1]['result']
def cs(name):return {k:dec(v) for k,v in result(name)['constants'].items()}
def up(v):return arb(abs(v).upper())
def execute():
 ctx.prec=512
 pi=arb.pi();rt2=arb(2).sqrt();ee=arb(1).exp();q=arb(2000)**4
 p={k:dec(v) for k,v in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 C,B2,beta=[p[k] for k in ['C','B2','beta']]
 mass={k:dec(v) for k,v in result('current-global-mass-conversion')['parameters'].items()}
 assert beta+abs(mass['u'])<arb(1)/5 and p['L']<=1
 assert mass['eps']<0 and mass['tau']>0
 assert ((q+1)/C).sqrt()>18
 shape=result('current-source-adapted-fractional-return');old=cs('current-source-adapted-fractional-return')
 S=[dec(v) for v in shape['derivative_accounting']['reconstructed_reference_derivative_constants']]
 D0=dec(shape['connection_tensor_derivative_bounds'][0])
 Z=cs('current-source-adapted-geodesic-first-attempt')['original_model_fractional_scale']
 U0=up(S[0]*Z);U1=up(arb(10).sqrt()*S[1]*Z);U2=up(10*(S[2]+D0*S[1])*Z)
 CV=cs('current-source-adapted-local-return-attempt4')['CV']
 epsilon=cs('current-complete-polynomial-weak-transfer')['epsilon']
 energy_factor=(1-epsilon)*arb('.99')*(1-arb('1.179e-8'))
 CG=(arb('1e5')*beta/energy_factor).sqrt()
 CD=cs('current-positive-time-quantum-operator')['CD_t']
 flow_test2=(3*CD**2/8+CG**2/4).sqrt()
 Lreal=2*rt2/(1/rt2-arb(1)/4);theta=arb('.026')*Lreal
 Kreal=rt2/pi/(1-theta)**2
 real_row_guard=2*pi*(arb('1.026')/arb('.974'))*Kreal**2
 assert theta<1 and real_row_guard<8
 Lcomplex=2*rt2*arb(5)/3/(arb('.75')/rt2-arb(1)/8)
 theta_c=arb('1.25')*arb('.026')*Lcomplex
 Kcomplex=5*rt2/(3*pi)/(1-theta_c)**2
 complex_row_guard=2*pi*(arb('1.026')/arb('.974'))*Kcomplex**2
 assert theta_c<1 and complex_row_guard<32
 sl=arb(260);sf=arb(512);cl=arb(1)/396;cf=arb(1)/789;cc=arb(1)/792
 def angular_e1(c):return 4/c**2 # coefficient of C
 def angular_e1_over_m2(c):
  assert 18*c<1
  l=(1/(18*c)).log();return 4*(l*l+(1+2/ee)*l+3/ee)
 sum_R2=8*sl**2*angular_e1(cl)
 sum_R2_over_m2=8*sl**2*angular_e1_over_m2(cl)
 local_cross=40*sum_R2;local_square=2*sl**2*sum_R2_over_m2
 far=16*sf**2*angular_e1(cf)
 geodesic_far=2*pi*18*arb(256)**4
 G0=2*(local_cross+local_square+far+geodesic_far)
 Pglue=(arb(40)/3)*32*sl**2*angular_e1(cc)
 factor=2/(B2*beta)
 split=arb('1e-12');K=arb('1e13');R0=arb(1)/400;b=(split/(4*K)).sqrt()
 assert b<R0
 alpha=split/3;bg=K*(arb(1)/4+2*(R0/b).log())
 debits={
  'local_gain_two':factor*2*pi*C*(alpha*U2**2+bg*U1**2).sqrt()*(alpha*flow_test2**2+bg*CG**2).sqrt(),
  'finite_second_contact':factor*arb('1e13')*C*U1*CG,
  'global_row_subtracted_gluing':factor*G0*C*U0*CV,
  'actual_constant_contact_cap':factor*(arb('1e18')*q*q/C)*(1+1/arb(5).sqrt())*U0*CV,
  'actual_constant_contact_gluing':factor*Pglue*C*U0*CV,
  'known_curvature_polynomial':factor*100*C*(U1*CG/(36*pi)+816*U0*CV/(960*pi)),
  'shared_vacuum_uncertainty':arb('1e-21')/energy_factor*old['actual_trial_Et_upper'],
  'retained_energy_comparison':epsilon/(1-epsilon)*old['actual_trial_Et_upper']}
 total=sum(debits.values(),arb(0));screen=arb('1e33');account=arb('4.308597e32')
 reference=old['complete_exact_analytic_first_reference_residual']
 rest=sum((v for k,v in debits.items() if k!='actual_constant_contact_cap'),arb(0))
 half=account/2
 thresholds={
  'cap_contact_multiplier_for_HALF_residual_account':half/(debits['actual_constant_contact_cap']/arb('1e18')),
  'cap_contact_report_reduction_for_HALF_account':debits['actual_constant_contact_cap']/half,
  'all_other_debits_report_reduction_for_HALF_account':rest/half,
  'remaining_without_cap_over_entire_residual_account':rest/account,
  'unweighted_zero_order_density_ceiling_if_IT_ALONE_spent_account':account/(factor*U0*CV),
  'gradient_density_ceiling_if_IT_ALONE_spent_account':account/(factor*U1*CG)}
 vals=dict(source_L2_upper=U0,source_gradient_upper=U1,source_Hessian_upper=U2,
  test_weighted_L2_constant=CV,test_gradient_constant=CG,test_Laplacian_constant=CD,
  real_Born_ratio=theta,real_squared_row_prefactor=real_row_guard,
  complex_Born_ratio=theta_c,complex_squared_row_prefactor=complex_row_guard,
  local_cross_coefficient=local_cross,local_square_coefficient=local_square,far_coefficient=far,
  geodesic_far_coefficient=geodesic_far,complete_global_gluing_coefficient=G0,
  constant_contact_gluing_coefficient=Pglue,
  global_gluing_coefficient_improvement=arb('6e24')/G0,
  constant_contact_gluing_coefficient_improvement=arb('3.3e24')/Pglue,
  complete_actual_minus_reference_upper=total,
  complete_actual_first_residual_upper=total+reference,
  complete_residual_upper_over_screen=(total+reference)/screen,
  complete_difference_upper_over_remaining_account=total/account,
  old_to_new_complete_upper_ratio=old['known_coarse_actual_minus_reference_upper']/total,
  other_debits_without_cap_upper=rest,
  dominant_cap_fraction=debits['actual_constant_contact_cap']/total)
 assert all(v.is_finite() and v>=0 for v in [*debits.values(),*vals.values(),*thresholds.values()])
 return dict(constants={k:enc(v) for k,v in vals.items()},
  debits={k:enc(v) for k,v in debits.items()},
  conditional_sensitivity_NOT_bounds_or_physical_margins={k:enc(v) for k,v in thresholds.items()},
  actual_signed_nonzero_center=None,numerical_application_error=None,signed_Gt=None,
  scope=dict(complete_fixed_trial_upper_estimator=True,actual_domain_form_extension=True,
   full_angular_thermal_contacts_and_returns_retained=True,complete_actual_first_inverse=False,
   original_admission=False,nonlinear_existence=False,compatibility=False,physical_matching=False,
   validity_noise=False,all_sector_stability=False,construction_obstruction=False,goal_complete=False)),total+reference<screen
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
 out=ROOT/'evidence'/(args.receipt+'.json');assert not out.exists(),'Preserve immutable receipts'
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes());pkg=rt['packages']['flint']
 assert digest(Path(sys.executable))==rt['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=rt['python_sha256'])
 try:
  rec['result'],passed=execute()
  rec['status']='complete_fixed_trial_estimator_passes' if passed else 'complete_fixed_trial_estimator_excessive_signed_evaluation_open'
  after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(input_sha256_after=after,protected_paths_unchanged=len(protected),pinned_flint_files_unchanged=len(pkg['files']))
 except Exception:rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()},
  debits={k:v['ball'] for k,v in rec.get('result',{}).get('debits',{}).items()}),indent=2))
 if rec['status']=='failed':raise SystemExit(1)
