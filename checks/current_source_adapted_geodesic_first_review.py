"""Author review: source bound, strong fractional coefficient, unresolved full inverse."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
ROOT=Path(__file__).resolve().parents[1]
STEM='current-source-adapted-geodesic-first-attempt'
def digest(p):return sha256(p.read_bytes()).hexdigest()
def up(v):return arb(v.upper())
def lo(v):return arb(v.lower())

ctx.prec=512
path=ROOT/'evidence'/(STEM+'-verification.json')
assert not path.exists()
d=json.loads((ROOT/'evidence'/(STEM+'.json')).read_bytes())
assert d['status']=='source_adapted_full_source_enclosed_actual_reference_inverse_unresolved'
assert d['input_sha256_before']==d['input_sha256_after']
for p,h in d['input_sha256_before'].items():assert digest(ROOT/p)==h,p
assert d['ordinary_argument_utf8']==(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8')
assert d['arithmetic_source_utf8']==(ROOT/'checks/current_source_adapted_geodesic_first_attempt.py').read_text(encoding='utf-8')
assert digest(Path(sys.executable))==d['python_sha256']
baseline=json.loads(Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json').read_bytes())
assert all(digest(Path(p))==h for p,h in baseline.items())
r=d['result'];c={k:dec(v) for k,v in r['constants'].items()}
assert c['old_full_mixed_remainder_debit']>arb('2.55e33')
assert c['complete_actual_source_debit']<arb('5.640405e32')
assert c['source_debit_over_first_stage_screen']<arb('.566591')
assert c['norm_gain']>arb('4.5242')
assert c['actual_cutoff_scalar_Es_upper']<arb('3.488669e48')
assert c['classical_A']>8 and c['classical_B']>17
assert c['K0']>arb('.4444444')
assert c['actual_bolt_rho4_log_coefficient']<arb('-2.2622e63')
assert c['model_fractional_profile_upper']<arb('1.224745')
assert c['actual_reference_fractional_output_coefficient']<arb('1.554891e-32')
assert c['full_first_residual_cost_per_Cflow']<arb('1.395154e17')
known=up(c['complete_actual_source_debit'])+up(c['complete_model_density_remainder'])
assert known<arb('5.642e32')
remaining=lo(arb('9.955e32')-known)
assert remaining>arb('4.313e32')
ceiling=lo(remaining/up(c['full_first_residual_cost_per_Cflow']))
norm_ceiling=lo(remaining/up(c['actual_reference_fractional_output_coefficient']))
assert ceiling>arb('3.09e15')
# At exponent six, the same integral recovers the prior accepted coefficient.
mu=arb(8)
def bm(nu):return arb(2)**(mu-3)*(mu/2).gamma()**2*(mu/2+nu).gamma()*(mu/2-nu).gamma()/mu.gamma()
six=(bm(1)-bm(0))/(720*2*arb.pi()**2)
assert six.overlaps(1/(525*arb.pi()**2))
# Algebraic high-mass cancellation, including Cv-1, at three exact symbols.
# The identity holds for all p by lambda_n=1/(d L^2); these check the arithmetic.
for p in [arb(0),arb(1),arb('1e9')]:
 lam=c['large_mass_lambda'];cv=c['Cv']
 assert ((cv+p/lam)/(1+p/lam)-(1+(cv-1)/(1+p/lam))).contains(0)
unknowns=['actual_flow_norm','actual_flow_comparison_constant','complete_finite_return_error',
 'actual_numerical_application_error','complete_first_reference_residual',
 'actual_complete_quantum_correction','signed_Gt']
assert all(r[k] is None for k in unknowns)
review=dict(
 source='Both full scalar and mixed identities retained. Exact Gram moments are used only for y_hat norm; original cutoff, geometry and full mixed remainder are all charged in Et*.',
 comparison_domain='Actual trial is defined by exact finite-Q York reconstruction; model resolvents define profiles only. Zero mass/log data, fixed phase and regular-bolt form closure retained.',
 compression_correction='A0 includes local rows, contacts, York pullback and Et Riesz identification. The prior shorthand equality with the nonlocal compression is corrected; prior numerical results do not depend on it.',
 regularity='Nonzero actual rho^4 log(rho) coefficient excludes H5/C5-source estimates. Dyadic local bounds establish fractional domain membership but not its numerical norm; final physical C5 admission remains unproved.',
 spectral_measure='Boundary measure is recovered from the L2 resolvent vector, Q Gram identity and digamma reflection. Its complete density is <=1/2. This describes only the auxiliary comparison, not the physical geodesic operator.',
 fractional_operator='Spectral theorem on actual Liouville flow; Bessel moments have mu=13/2>2; complete angular majorant converges. Exact UV logarithm and all flow returns retained; local contacts separate.',
 remaining_gap='The actual fractional source norm and complete local/finite York return remain unknown. A conditional ceiling cannot replace either estimate. Later actual-minus-reference quantum applications remain unexecuted.',
 decision='Retain source-adapted trial and fractional strong evaluation; reject an H5 source proof. Next bound the joint actual-domain return and fractional column in the complete residual inequality.',
 claims='No first inverse, signed Gt, original admission, existence, matching, validity/noise or stability is certified. No construction obstruction.',
 verification='Author ordinary-argument and outward-arithmetic/provenance review; no independent or Lean review.')
out=dict(status='complete_source_and_fractional_operator_reviewed_actual_inverse_unresolved',
 finished_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(ROOT/'evidence'/(STEM+'.json')),
 review_source_sha256=digest(Path(__file__)),review_source_utf8=Path(__file__).read_text(encoding='utf-8'),
 checked_inputs=len(d['input_sha256_before']),protected_paths_unchanged=len(baseline),review=review,
 constants={k:enc(v) for k,v in dict(known_first_stage_upper=known,
  remainder_before_UNKNOWN_entries=remaining,
  conditional_Cflow_ceiling_RETURN_AND_NUM_ZERO=ceiling,
  conditional_actual_flow_norm_ceiling_RETURN_AND_NUM_ZERO=norm_ceiling,
  recovered_sixth_coefficient=six).items()},
 first_reference_inverse_certified=False,complete_original_admission=False,goal_complete=False)
path.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf-8'))
print(json.dumps(dict(status=out['status'],checked_inputs=out['checked_inputs'],
 constants={k:v['ball'] for k,v in out['constants'].items()}),indent=2))
