"""Full sequential Dj estimator test; no actual signed quantum certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
import argparse
import json
import sys
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb, ctx, dec, enc

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-filtered-j-regularity-test'
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
PRIOR = 'evidence/current-correlated-signed-trial-test-attempt3.json'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def result(stem):
    raw = json.loads((ROOT / 'evidence' / (stem + '.json')).read_bytes())
    return raw.get('runs', [raw])[-1]['result']


def constants(stem):
    return {k: dec(v) for k, v in result(stem)['constants'].items()}


def qdec(v):
    return Fraction(int(v['numerator']), int(v['denominator']))


def qb(v):
    return arb(v.numerator) / arb(v.denominator)


def qe(v):
    return dict(numerator=str(v.numerator), denominator=str(v.denominator))


def mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def poles(values):
    out = [Fraction(1)]
    for v in values:
        out = mul(out, [v, Fraction(1)])
    return out


def execute():
    ctx.prec = 768
    prior = result('current-correlated-signed-trial-test-attempt3')
    v = {k: qdec(x) for k, x in prior['exact_nominal_parameters'].items()}
    p = {k: dec(x) for k, x in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
    green = constants('current-filtered-mixed-green-residual')
    geo = constants('current-source-adapted-geodesic-first-attempt')
    frac = constants('current-source-adapted-fractional-return')
    detail = result('current-source-adapted-fractional-return')['derivative_accounting']
    S0 = dec(detail['reconstructed_reference_derivative_constants'][0])
    Kgeo = geo['actual_reference_fractional_output_coefficient']
    CLt = constants('current-positive-time-quantum-operator')['CL_t']
    assert CLt > 1  # Kgeo is a conservative transposed bilinear coefficient.
    Cflow = frac['actual_reconstruction_fractional_comparison']
    coefficient = Kgeo*Cflow
    qt = arb('.338')/(1-green['epsilon'])
    ft = qt/(1+qt)
    gamma = p['gamma']
    radius = arb('1e-110')
    P = gamma/radius
    target = arb('1.49e34')
    A0 = 3*v['alpha']/(4*v['gamma']*v['A'])
    assert A0 > 0 and v['kappa'] > 0 and v['delta'] > 0
    assert P/2 > 100
    assert all(0 < v[k] < 100 for k in ['twice_a', 'muS', 'muT', 'muY', 'scalar_c'])
    density_at_band_start = (arb.pi()*(P/(2*p['a'])-arb(1)/4).sqrt()).tanh()/2
    assert density_at_band_start > arb(1)/4
    D_ratio = 1 + qb(v['d1'])/(P/2) + qb(v['d0'])/(P/2)**2
    assert D_ratio < 6
    lower_tests = []
    for k in [0, 1]:
        exponent = arb(7)/2-2*k
        integral = (P**exponent-(P/2)**exponent)/exponent
        Xlower = (gamma*qb(A0)**2/1024*(qb(v['kappa']*v['delta'])/6)**(2*k)*integral).sqrt()
        debit_lower = ft**(3-k)*coefficient*Xlower
        assert debit_lower > target
        lower_tests.append(dict(stage=k, X_lower=enc(Xlower),
            complete_estimator_geodesic_debit_lower=enc(debit_lower),
            debit_lower_over_Dj_target=enc(debit_lower/target),
            equivalent_original_bolt_debit_lower=enc(green['pref_t']*debit_lower),
            spectral_power_integral=enc(integral)))

    N = [-A0*x for x in poles([v['muY'], v['muS']+2*v['scalar_c']])]
    D = poles([v['twice_a'], v['muS'], v['muT']])
    lin = [v['kappa']*v['Cv'], v['kappa']*v['delta']]
    quad = [v['d0'], v['d1'], Fraction(1)]
    profiles = []
    norm_values = []
    for k in range(4):
        degree = 2*k+3
        assert len(N) == k+3 and len(D) == degree+1
        assert all(x > 0 for x in D)
        g = min(x/Fraction(comb(degree, j)) for j, x in enumerate(D))
        M = sum(abs(x) for x in N)
        Lbound = (gamma/(2*(2*k+1))).sqrt()*qb(M/g)
        item = dict(stage=k, numerator=[qe(x) for x in N], denominator=[qe(x) for x in D],
            Bernstein_guard=qe(g), coefficient_sum=qe(M), physical_profile_L2_upper=enc(Lbound))
        Xbound = None
        if k >= 2:
            # Half the integral of (1+p)^(9/2-2(k+1)) over p>=0.
            Xbound = (gamma/(4*k-7)).sqrt()*qb(M/g)
            item['full_spectrum_X_upper'] = enc(Xbound)
        profiles.append(item)
        norm_values.append((Lbound, Xbound))
        N, D = mul(N, lin), mul(D, quad)

    L0, L1 = norm_values[0][0], norm_values[1][0]
    X2, X3 = norm_values[2][1], norm_values[3][1]
    Atest = Kgeo*S0*(qt*L0+L1)
    Btest = coefficient*(ft*X2+X3)
    optimistic = (target-Btest)/Atest
    vals = dict(original_gamma=gamma, averaging_radius=radius, band_end_P=P,
        band_density_lower_guard=density_at_band_start, D_over_p_squared_guard=D_ratio,
        original_qt=qt, original_ft=ft, Kgeo=Kgeo, Cflow=Cflow,
        retained_test_embedding_CLt=CLt,
        complete_strong_column_coefficient=coefficient,
        actual_zero_derivative_reconstruction_S0=S0,
        direct_test_coefficient_Atest=Atest,
        remaining_geodesic_source_charge_Btest=Btest,
        Btest_over_Dj_target=Btest/target,
        optimistic_Ktest_ceiling_ALL_other_errors_zero=optimistic,
        prescribed_Dj_target=target,
        source_norm_reduction_X0_lower_over_X3_upper=dec(lower_tests[0]['X_lower'])/X3,
        source_norm_reduction_X1_lower_over_X3_upper=dec(lower_tests[1]['X_lower'])/X3)
    assert all(x.is_finite() and x > 0 for x in vals.values())
    return dict(constants={k: enc(x) for k, x in vals.items()},
        exact_nominal_parameters=prior['exact_nominal_parameters'],
        strong_estimator_lower_tests=lower_tests, whole_spectrum_profile_bounds=profiles,
        actual_Dj=None, Ktest=None, Erest=None, actual_signed_Zr=None,
        scope=dict(actual_averaged_structured_trials_used=True,
            original_physical_norms_and_radius_cost_retained=True,
            complete_sequential_additive_estimator_rejected=True,
            full_angular_thermal_action_retained_in_Kgeo=True,
            exact_all_four_residual_identity=True,
            actual_quantum_response_or_Dj_bounded=False,
            direct_filtered_test_regularity_proved=False,
            complete_alternative_bound_proved=False,
            available_physical_margin_established=False,
            actual_bolt_screen_rejected=False, original_admission=False,
            nonlinear_existence=False, exact_compatibility=False, physical_matching=False,
            validity_noise=False, all_sector_stability=False,
            construction_obstruction=False, goal_complete=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--receipt', default=STEM)
    args = ap.parse_args()
    out = ROOT/'evidence'/(args.receipt+'.json')
    assert not out.exists(), 'Preserve immutable receipts'
    old = json.loads((ROOT/PRIOR).read_bytes())
    inputs = list(dict.fromkeys([
        'notes/'+STEM+'.md', 'checks/current_filtered_j_regularity_test.py', PRIOR,
        'evidence/current-correlated-signed-trial-test-attempt3-verification.json',
        'notes/current-source-adapted-fractional-return.md',
        'notes/current-source-adapted-geodesic-first-attempt.md',
        'notes/current-seed-full-hessian-realization.md',
        'notes/current-geodesic-two-sided-remainder.md',
        *old['input_sha256_after']]))
    before = {p: digest(ROOT/p) for p in inputs}
    protected = json.loads(BASELINE.read_bytes())
    rec = dict(started_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256_before=before,
        ordinary_argument_utf8=(ROOT/inputs[0]).read_text(encoding='utf-8'),
        arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'))
    try:
        assert all(digest(ROOT/p) == h for p, h in old['input_sha256_after'].items())
        assert all(digest(Path(p)) == h for p, h in protected.items())
        rt = json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
        pkg = rt['packages']['flint']
        assert digest(Path(sys.executable)) == rt['python_sha256']
        assert all(digest(Path(pkg['root'])/p) == h for p, h in pkg['files'].items())
        rec['result'] = execute()
        after = {p: digest(ROOT/p) for p in inputs}
        assert after == before
        assert all(digest(Path(p)) == h for p, h in protected.items())
        rec.update(status='sequential_strong_estimator_rejected_actual_Dj_unresolved',
            input_sha256_after=after, protected_paths_unchanged=len(protected),
            python_sha256=rt['python_sha256'], pinned_flint_files_unchanged=len(pkg['files']))
    except Exception:
        rec.update(status='failed', traceback=traceback.format_exc())
    rec['finished_utc'] = datetime.now(timezone.utc).isoformat()
    out.write_bytes((json.dumps(rec, indent=2)+'\n').encode('utf-8'))
    print(json.dumps(dict(status=rec['status'], error=rec.get('traceback'),
        constants={k: x['ball'] for k, x in rec.get('result', {}).get('constants', {}).items()},
        lower_tests=rec.get('result', {}).get('strong_estimator_lower_tests')), indent=2))
    if rec['status'] == 'failed':
        raise SystemExit(1)
