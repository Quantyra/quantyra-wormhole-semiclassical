"""Author verification of the ordinary range test and frozen outward evidence."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc

ROOT=Path(__file__).resolve().parents[1]
STEM='current-localized-first-residual-range-test'
def digest(p):return sha256(p.read_bytes()).hexdigest()
def read(n):
 d=json.loads((ROOT/'evidence'/(n+'.json')).read_bytes())
 return d.get('runs',[d])[-1]
def up(x):return arb(x.upper())
def lo(x):return arb(x.lower())
def mag(x):return up(abs(x))

ctx.prec=512
path=ROOT/'evidence'/(STEM+'-verification.json')
assert not path.exists(),'Preserve immutable verification.'
d=read(STEM);assert d['status']=='actual_source_one_column_residual_screen_rejected'
assert d['input_sha256_before']==d['input_sha256_after']
for p,h in d['input_sha256_before'].items():assert digest(ROOT/p)==h,p
assert digest(Path(sys.executable))==d['python_sha256']
baseline=json.loads(Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json').read_bytes())
assert all(digest(Path(p))==h for p,h in baseline.items())
assert d['ordinary_argument_utf8']==(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8')
assert d['arithmetic_source_utf8']==(ROOT/'checks/current_localized_first_residual_range_test.py').read_text(encoding='utf-8')
for i in [1,2]:assert read(STEM+'-attempt'+str(i))['status']=='failed'
r=d['result'];cs={k:dec(v) for k,v in r['constants'].items()}
ts=[{k:dec(v) for k,v in t.items() if isinstance(v,dict) and 'decimal_enclosure' in v} for t in r['tests']]
G1,D1,E1=[ts[0][k] for k in ['source_signed','response_signed','test_energy_upper']]
G2,D2,E2=[ts[1][k] for k in ['source_signed','response_signed','test_energy_upper']]
At,ga=cs['source_amplitude'],cs['physical_measure_gamma']
det=D2*G1-D1*G2
assert det<arb('-4.8') and det>arb('-12.504')
lower=lo(lo(abs(det))*lo(At)*lo(ga)/up(mag(D2)*up(E1)+mag(D1)*up(E2)))
assert lower>arb('1.95e42')
assert cs['complete_first_residual_lower'].overlaps(lower) or cs['complete_first_residual_lower']>lower
allowances=[]
for G,D,E in [(G1,D1,E1),(G2,D2,E2)]:
 assert D>0
 tol=arb('1e33')*up(E)/(At*ga)
 allowed=(G+arb(0,mag(tol)))/D
 allowances.append(allowed)
assert up(allowances[0])<lo(allowances[1])
fs={k:dec(v) for k,v in read('current-first-inverse-source-screen')['result']['constants'].items()}
first={k:dec(v) for k,v in read('current-first-quantum-solve-attempt')['result']['constants'].items()}
K=fs['first_residual_to_original_bolt_multiplier']
assert K*lower>arb('.0624')
assert K*lower/arb('1.667e-9')>arb('3.74e7')
lam=cs['explicit_trial_coefficient_over_At']
explicit_second=G2-lam*D2
assert explicit_second>0
explicit_lower=lo(lo(At)*lo(ga)*lo(explicit_second)/up(E2))
explicit_upper=up(first['proved_source_dual_upper']+mag(lam*At)*up(cs['psi_full_energy_upper']))
next_residual=lo((1-first['complete_geodesic_form_error'])*first['conditional_common_stage_output_error_ceiling'])
assert first['exact_application_first_residual_error']<arb('2.431e19')
assert next_residual>arb('9.955e32')
# Coarse comparison ledger: ordinary termwise product estimates, not a
# substitute for the actual-domain argument or a new quantum mode cutoff.
ledger={'Einstein':arb('1e4'),'magnetic':arb('1e3'),'R2':arb('1e6'),
        'Weyl':arb('2e7'),'canonical_local':arb('1e5')}
assert sum(ledger.values(),arb(0))<arb('1e8')
review=dict(
 first_domain_step='Both York identities and endpoint data follow from exact I_b derivative; cutoff fields are C6 and bolt-analytic. Scalar y_hat is used only in its proved Es domain.',
 energy_comparison='Original dz and z weights retained. Full geodesic logarithm uses Dom(X^3), complete angular sum, and the accepted .0022 reference comparison; no physical C5 inference.',
 complete_source='Exact polynomial moments on the actual compact test support; full mixed remainder bounded by the accepted adjoint. Near-integer coefficient and all large-degree endpoint errors retained.',
 complete_quantum='Accepted fixed-psi short-time column plus positive-time operator, shared vacuum, a2, LLL, contact reallocation and transport error all charged. No independent phase or neighborhood extension.',
 cancellation='Exact response cancels in h=D2*xi1-D1*xi2; full interval determinant stays negative. This is an actual-source lower test valid for every scalar c.',
 coefficient_scope='Explicit c=At times the recorded exact binary-rational center. No claim that this is the Ritz coefficient. Lower bound covers all coefficient errors staying in span{psi}.',
 mathematical_decision='Change the one-dimensional correction range. Retain the accepted localized quantum evaluator and use a complete source-adapted principal inverse; its applications remain unexecuted.',
 claims_boundary='No lower bound on actual bolt error, no complete Gt or admission margin, no inverse/construction obstruction, no existence or stability result.',
 verification_scope='Author ordinary-argument review and frozen outward-arithmetic/provenance verification; no independent or Lean verification.')
out=dict(status='actual_source_fixed_range_failure_verified_full_inverse_open',
 finished_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(ROOT/'evidence'/(STEM+'.json')),
 review_source_sha256=digest(Path(__file__)),review_source_utf8=Path(__file__).read_text(encoding='utf-8'),
 checked_input_count=len(d['input_sha256_before']),protected_paths_unchanged=len(baseline),
 preserved_failed_attempts=2,review=review,constants={k:enc(v) for k,v in dict(
  first_allowed_coefficient_interval=allowances[0],second_allowed_coefficient_interval=allowances[1],
  coefficient_gap_lower=lo(allowances[1])-up(allowances[0]),
  verified_universal_first_residual_lower=lower,
  explicit_trial_first_residual_lower=explicit_lower,explicit_trial_complete_residual_upper=explicit_upper,
  verified_original_first_estimator_floor=K*lower,
  next_geodesic_application_residual_target=next_residual,
  alternative_exact_application_analytic_error=first['exact_application_first_residual_error']).items()},
 complete_first_inverse_certified=False,complete_Gt_certified=False,original_admission_closed=False,
 nonlinear_existence_compatibility_matching_validity_noise_stability_closed=False,goal_complete=False)
path.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf-8'))
print(json.dumps(dict(status=out['status'],checked_inputs=out['checked_input_count'],
 constants={k:v['ball'] for k,v in out['constants'].items()}),indent=2))
