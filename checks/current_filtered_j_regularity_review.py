"""Author review: structured-range estimator failure, with no Dj promotion."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
import argparse
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_filtered_j_regularity_test as m
arb, ctx, dec, enc = m.arb, m.ctx, m.dec, m.enc
ROOT = m.ROOT
digest = lambda p: sha256(p.read_bytes()).hexdigest()

ap = argparse.ArgumentParser()
ap.add_argument('--receipt', default=m.STEM+'-attempt2')
args = ap.parse_args()
source = ROOT/'evidence'/(args.receipt+'.json')
raw = json.loads(source.read_bytes())
assert raw['status'] == 'sequential_strong_estimator_rejected_actual_Dj_unresolved'
assert raw['input_sha256_before'] == raw['input_sha256_after']
assert all(digest(ROOT/p) == h for p, h in raw['input_sha256_after'].items())
assert raw['ordinary_argument_utf8'] == (ROOT/'notes'/(m.STEM+'.md')).read_text(encoding='utf-8')
assert raw['arithmetic_source_utf8'] == (ROOT/'checks/current_filtered_j_regularity_test.py').read_text(encoding='utf-8')
protected = json.loads(m.BASELINE.read_bytes())
assert all(digest(Path(p)) == h for p, h in protected.items())
rt = json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
pkg = rt['packages']['flint']
assert digest(Path(sys.executable)) == rt['python_sha256']
assert all(digest(Path(pkg['root'])/p) == h for p, h in pkg['files'].items())

ctx.prec = 768
r = raw['result']
c = {k: dec(x) for k, x in r['constants'].items()}
prior = m.result('current-correlated-signed-trial-test-attempt3')
assert r['exact_nominal_parameters'] == prior['exact_nominal_parameters']
v = {k: m.qdec(x) for k, x in r['exact_nominal_parameters'].items()}
assert c['retained_test_embedding_CLt'] > 1
assert c['complete_strong_column_coefficient'].overlaps(c['Kgeo']*c['Cflow'])
assert c['band_end_P'].overlaps(c['original_gamma']/c['averaging_radius'])
assert c['band_end_P']/2 > 100
assert c['band_density_lower_guard'] > arb(1)/4
assert c['D_over_p_squared_guard'] < 6

profiles = r['whole_spectrum_profile_bounds']
assert len(profiles) == 4
for k, profile in enumerate(profiles):
    assert profile['stage'] == k
    N = [m.qdec(x) for x in profile['numerator']]
    D = [m.qdec(x) for x in profile['denominator']]
    assert len(D)-len(N) == k+1 and all(x > 0 for x in D)
    g, M = m.qdec(profile['Bernstein_guard']), m.qdec(profile['coefficient_sum'])
    assert g > 0 and M == sum(abs(x) for x in N)
    assert all(x >= g*comb(len(D)-1, j) for j, x in enumerate(D))
    assert all(abs(x) <= M*comb(len(N)-1, j) for j, x in enumerate(N))
    # Exact full-spectrum integral constants, independently of the formula
    # used in execute(): half int_0^infty (1+p)^(-2(k+1)) dp.
    half_L2_integral = Fraction(1, 4*k+2)
    expected = (c['original_gamma']*m.qb(half_L2_integral)).sqrt()*m.qb(M/g)
    assert expected.overlaps(dec(profile['physical_profile_L2_upper']))
    if k >= 2:
        half_X_integral = Fraction(1, 4*k-7)
        expected = (c['original_gamma']*m.qb(half_X_integral)).sqrt()*m.qb(M/g)
        assert expected.overlaps(dec(profile['full_spectrum_X_upper']))
    if k:
        previous = profiles[k-1]
        assert N == m.mul([m.qdec(x) for x in previous['numerator']],
                          [v['kappa']*v['Cv'], v['kappa']*v['delta']])
        assert D == m.mul([m.qdec(x) for x in previous['denominator']],
                          [v['d0'], v['d1'], Fraction(1)])

# The final analytic reference is EXACTLY the previously tested averaged
# source profile, including every pole and the original amplitude.
av = prior['averaged_reference_comparison']
assert [m.qdec(x) for x in profiles[3]['numerator']] == [v['At']*m.qdec(x) for x in av['numerator']]
assert profiles[3]['denominator'] == av['denominator']

lower = r['strong_estimator_lower_tests']
assert dec(lower[0]['X_lower']) > arb('1.993e181')
assert dec(lower[1]['X_lower']) > arb('6.635e95')
assert dec(lower[0]['complete_estimator_geodesic_debit_lower']) > arb('2.683e157')
assert dec(lower[1]['complete_estimator_geodesic_debit_lower']) > arb('3.536e72')
assert dec(lower[0]['debit_lower_over_Dj_target']) > arb('1.801e123')
assert dec(lower[1]['debit_lower_over_Dj_target']) > arb('2.373e38')
assert c['direct_test_coefficient_Atest'] < arb('1.542e20')
assert c['remaining_geodesic_source_charge_Btest'] < arb('2.971e26')
assert c['Btest_over_Dj_target'] < arb('1.994e-8')
ceiling = c['optimistic_Ktest_ceiling_ALL_other_errors_zero']
assert arb('9.667e13') < ceiling < arb('9.668e13')
assert (c['direct_test_coefficient_Atest']*ceiling+c['remaining_geodesic_source_charge_Btest']).overlaps(c['prescribed_Dj_target'])

for key in ['actual_Dj', 'Ktest', 'Erest', 'actual_signed_Zr']:
    assert r[key] is None
for key in ['actual_quantum_response_or_Dj_bounded', 'direct_filtered_test_regularity_proved',
            'complete_alternative_bound_proved', 'available_physical_margin_established',
            'actual_bolt_screen_rejected', 'original_admission', 'nonlinear_existence',
            'exact_compatibility', 'physical_matching', 'validity_noise', 'all_sector_stability',
            'construction_obstruction', 'goal_complete']:
    assert r['scope'][key] is False

history = []
for path in sorted((ROOT/'evidence').glob(m.STEM+'*.json')):
    if path == source or 'verification' in path.name:
        continue
    prior_raw = json.loads(path.read_bytes())
    history.append(dict(path=str(path.relative_to(ROOT)), sha256=digest(path), status=prior_raw['status']))
out = ROOT/'evidence'/(args.receipt+'-verification.json')
assert not out.exists()
record = dict(status='structured_sequential_estimator_failure_reviewed_Dj_open',
    reviewed_utc=datetime.now(timezone.utc).isoformat(), receipt_sha256=digest(source),
    review_source_sha256=digest(Path(__file__)), frozen_inputs_checked=len(raw['input_sha256_after']),
    protected_paths_checked=len(protected),
    review_type='author ordinary domain, four-row identity, exact averaged profile, spectral integral and budget review',
    independent_specialist_review=False, lean_verification=False,
    exact_previous_reference_identity_checked=True, physical_averaging_radius_unchanged=True,
    sequential_estimator_rejected=True, actual_Dj_bounded=False,
    filtered_phase_regularity_unresolved=True, construction_obstruction=False,
    goal_complete=False, preserved_receipts=history)
out.write_bytes((json.dumps(record, indent=2)+'\n').encode('utf-8'))
print(json.dumps(record, indent=2))
