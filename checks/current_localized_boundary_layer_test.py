"""Actual fifth-jet density estimator test; full quantum column unresolved."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
ROOT=Path(__file__).resolve().parents[1]
STEM='current-localized-boundary-layer-test'
INPUTS=['notes/'+STEM+'.md','checks/current_localized_boundary_layer_test.py',
 'checks/current_quartic_ritz_rank_test.py','notes/current-localized-short-time-trial.md',
 'notes/current-rough-paired-variation.md','notes/current-local-boundary-heat-remainder.md',
 'notes/current-curved-paired-heat-bridge.md','notes/current-disk-heat-propagation.md',
 'notes/current-smeared-heat-projection-domain.md',
 'evidence/current-global-green-trial-enclosure.json',
 'evidence/current-localized-short-time-trial.json',
 'evidence/current-positive-time-quantum-operator.json',
 'evidence/current-first-quantum-solve-attempt.json',
 'evidence/current-first-inverse-source-screen.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
def digest(p):return sha256(p.read_bytes()).hexdigest()
def result(name):
 d=json.loads((ROOT/'evidence'/(name+'.json')).read_bytes())
 r=d.get('runs',[d])[-1];assert r['status']!='failed'
 return r['result']
def execute():
 ctx.prec=512
 p={k:dec(v) for k,v in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 L,gamma,beta,C,B2=[p[k] for k in ['L','gamma','beta','C','B2']]
 B0=gamma+L
 a4=-arb(15)/2*beta*gamma**3/B0**4
 a5=-36*beta*gamma**3/B0**5
 c5=a5*beta**(arb(5)/2)
 assert c5<0
 pi=arb.pi();phi0=-15/(32*pi.sqrt())
 window_min=-phi0-1/(4*pi)
 lower=(arb(1)/8).sqrt()*window_min
 upper=(32*arb(7).sqrt()/(147*pi**3)).sqrt()
 assert 0<lower<upper
 asymptotic_lower=beta**(arb(1)/4)*abs(c5)*lower
 B5=abs(c5)*upper
 old=result('current-localized-short-time-trial')['constants']
 plus=result('current-positive-time-quantum-operator')['constants']
 CLt=dec(plus['CL_t']);k=dec(plus['total_operator_difference_upper'])
 Elower=dec(old['trial_Et_lower']);vnorm=dec(old['original_v_L2'])
 t0=arb('.001');pref=4*C*CLt/(B2*beta)
 d5=pref*arb(4)/3*t0**(arb(3)/4)*B5
 tau5=d5/Elower
 tau_req=dec(old['conditional_relative_column_screen'])
 tau6=dec(old['attained_leading_column_upper'])
 Ghat=dec(result('current-first-quantum-solve-attempt')['constants']['proved_source_dual_upper'])
 Kfirst=dec(result('current-first-inverse-source-screen')['constants']['first_residual_to_original_bolt_multiplier'])
 prefrel=pref*CLt
 oldpower=prefrel*arb(4)/7*t0**(arb(7)/4)
 newpower=prefrel*arb(4)/3*t0**(arb(3)/4)
 newceil=(tau_req-tau6)/newpower
 splitceil=(tau_req-tau6-tau5)/oldpower
 assert splitceil>0
 constants=dict(actual_z4_coefficient=a4,actual_z5_coefficient=a5,actual_proper_y5_coefficient=c5,
  complete_profile_value_at_zero=phi0,profile_abs_lower_on_zero_to_one_eighth=window_min,
  profile_L2_window_lower=lower,profile_L2_halfline_upper=upper,
  actual_remainder_tminus7over4_liminf_lower=asymptotic_lower,
  explicit_layer_original_L2_t7over4_coefficient_upper=B5,
  explicit_layer_full_angular_Et_dual_upper=d5,
  explicit_layer_relative_column_upper=tau5,
  explicit_layer_fraction_of_conditional_column_screen=tau5/tau_req,
  explicit_layer_original_first_error_envelope=Kfirst*Ghat*tau5/(1-k),
  power_change_cost_for_same_boundary_coefficient=newpower/oldpower,
  conditional_unsplit_t7over4_boundary_relative_ceiling_OTHER_ZERO=newceil,
  conditional_NEW_remainder_t11over4_boundary_relative_ceiling_OTHER_ZERO=splitceil,
  explicit_layer_relative_L2_coefficient=B5/vnorm)
 assert all(v.is_finite() for v in constants.values())
 return dict(precision_bits=ctx.prec,constants={k:enc(v) for k,v in constants.items()},
  original_L2_density_ansatz_no_finite_constants=True,
  argument_basis='actual nonzero fifth jet and complete paired parabolic scaling limit; arithmetic only verifies coefficients and budget',
  weak_Et_dual_t11over4_bound_disproved=False,actual_inverse_or_construction_obstructed=False,
  actual_minus_layer_remainder_constant=None,actual_Ritz_coefficient=None,complete_signed_first_residual=None,
  scope=dict(actual_trial_and_domain_retained=True,complete_profile_Fourier_integrals=True,
   explicit_layer_all_angular_regulator_removal=True,full_actual_quantum_remainder_enclosed=False,
   complete_quantum_certificate=False,new_full_admission_margin=False,
   nonlinear_existence_matching_validity_stability=False,goal_complete=False))
if __name__=='__main__':
 out=ROOT/'evidence'/(STEM+'.json');assert not out.exists()
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes());pkg=rt['packages']['flint']
 assert digest(Path(sys.executable))==rt['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=rt['python_sha256'])
 try:
  rec['result']=execute()
  after={p:digest(ROOT/p) for p in INPUTS};assert before==after
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(status='unsplit_density_estimator_rejected_explicit_layer_bounded_full_response_unresolved',
   input_sha256_after=after,protected_paths_unchanged=len(protected))
 except Exception:
  import traceback
  rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 keys=['actual_proper_y5_coefficient','actual_remainder_tminus7over4_liminf_lower',
  'explicit_layer_relative_column_upper','explicit_layer_original_first_error_envelope',
  'power_change_cost_for_same_boundary_coefficient',
  'conditional_unsplit_t7over4_boundary_relative_ceiling_OTHER_ZERO']
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:rec['result']['constants'][k]['ball'] for k in keys if 'result' in rec}),indent=2))
 if rec['status']=='failed':raise SystemExit(1)
