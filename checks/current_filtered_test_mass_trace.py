"""Source-specific full-response bridge and exact mixed mass-trace coefficient."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import argparse
import json
import sys
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_filtered_j_regularity_test as previous
arb, ctx, dec, enc = previous.arb, previous.ctx, previous.dec, previous.enc
ROOT = previous.ROOT
STEM = 'current-filtered-test-mass-trace'
BASELINE = previous.BASELINE
runtime = json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['sympy']['root']).parent))
import sympy as sp


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def result(stem):
    raw = json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())
    return raw.get('runs', [raw])[-1]['result']


def constants(stem):
    return {k: dec(v) for k, v in result(stem)['constants'].items()}


def evaluate(expr, values):
    if expr == sp.pi:
        return arb.pi()
    if expr.is_Rational:
        return arb(int(expr.p))/arb(int(expr.q))
    if expr.is_Symbol:
        return values[str(expr)]
    if expr.is_Add:
        return sum((evaluate(v, values) for v in expr.args), arb(0))
    if expr.is_Mul:
        out = arb(1)
        for v in expr.args:
            out *= evaluate(v, values)
        return out
    if expr.is_Pow and expr.exp.is_Integer:
        return evaluate(expr.base, values)**int(expr.exp)
    raise AssertionError('Unexpected coefficient expression: '+str(expr))


def execute():
    ctx.prec = 768
    p = {k: dec(v) for k, v in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
    mp = {k: dec(v) for k, v in result('current-global-mass-conversion')['parameters'].items()}
    green = constants('current-filtered-mixed-green-residual')
    epsilon = constants('current-complete-polynomial-weak-transfer')['epsilon']
    CU = constants('current-source-adapted-local-return-attempt4')['CU']
    alpha = green['alpha_original_radius_cost']
    q = arb('.338')
    f = q/(1+q)
    radius = arb('1e-110')
    Wmax = 2*(arb(4)/3).exp()
    density = p['L']**2*(Wmax/radius).sqrt()
    J = alpha*density*(CU*(1+epsilon).sqrt()+(q/arb('23.99')).sqrt())
    inverse_factor = epsilon*(4-6*epsilon+4*epsilon**2-epsilon**3)/(1-epsilon)**4
    bridge = (1+epsilon).sqrt()*f**3*inverse_factor*J
    target = arb('1.49e34')
    assert 0 < epsilon < arb('1e-70')
    assert bridge < arb('1e14') and bridge/target < arb('1e-20')

    phase_raw = json.loads((ROOT/'evidence/current-schur-residual-phase-domain.json').read_bytes())
    phase_run = phase_raw.get('runs', [phase_raw])[-1]
    assert phase_run['status'] == 'passed'
    phase = phase_run['result']
    symbols = {k: sp.Symbol(k) for k in ['z', 'C', 'B2', 'beta', 'eps', 'cR', 'tau', 'qflux']}
    symbols['pi'] = sp.pi
    z, C, B2, beta, eps, cR, qflux = [symbols[k] for k in ['z', 'C', 'B2', 'beta', 'eps', 'cR', 'qflux']]
    rows = {k: sp.sympify(v, locals=symbols) for k, v in phase['complete_mixed_coefficients'].items()}
    valuations = {}
    for key, expr in rows.items():
        numerator, denominator = sp.fraction(sp.cancel(expr))
        assert not denominator.has(z)
        pol = sp.Poly(numerator, z)
        valuations[key] = None if pol.is_zero else min(k[0] for k, _ in pol.terms())
    assert valuations == phase['whole_mixed_outer_valuations']
    assert valuations['xi0'] == 4 and valuations['J_Y'] == 4
    assert all(valuations[k] >= 3 for k in ['xi1', 'xi2', 'b_radial'])
    c04 = sp.cancel(sp.expand(rows['xi0']).coeff(z, 4))
    lam = C*beta/B2
    Fstar = 1-384*sp.pi*cR*lam+23*lam/(180*sp.pi)
    expected = 4*(3*beta+eps+1/C)*Fstar-qflux**2/(2*sp.pi*B2*C)
    assert sp.cancel(c04-expected) == 0
    values = dict(C=p['C'], B2=p['B2'], beta=p['beta'], eps=mp['eps'],
                  cR=arb('1e22'), qflux=arb(2000)**4)
    coefficient = evaluate(expected, values)
    direct = evaluate(c04, values)
    assert coefficient.overlaps(direct)
    assert arb('1.99') < coefficient < arb('2.01')
    physical_coefficient = coefficient/(8*arb.pi()*p['B2']/p['C'])
    assert physical_coefficient > 0
    vals = dict(retained_full_form_error=epsilon, prescribed_radius=radius,
        normalized_density_supremum_upper=Wmax, physical_Hmu_density_upper=density,
        complete_source_Edual_upper=J, four_inverse_perturbation_factor=inverse_factor,
        full_to_retained_j_energy_error=bridge,
        bridge_over_complete_Dj_target=bridge/target,
        bridge_original_bolt_charge=green['pref_t']*bridge,
        conditional_full_response_error_ceiling=target-bridge,
        full_E_to_retained_Et_norm_comparison=(1+epsilon).sqrt(),
        conditional_full_E_error_ceiling=(target-bridge)/(1+epsilon).sqrt(),
        complete_mixed_c04=coefficient,
        physical_mixed_adjoint_leading_coefficient=physical_coefficient)
    assert all(v.is_finite() and v > 0 for v in vals.values())
    return dict(constants={k: enc(v) for k, v in vals.items()},
        complete_mixed_valuations=valuations,
        extracted_c04=str(c04), factored_c04=str(expected),
        exact_coefficient_identity=True,
        actual_Dj=None, actual_signed_Zr=None, numerical_mass_trace_norm=None,
        complete_signed_mass_pairing=None,
        scope=dict(complete_averaged_j_response_transfer_certified=True,
            original_density_and_Riesz_map_retained=True,
            transport_depth_or_P32_degree_changed=False,
            ordinary_full_filtered_range_mass_nonzero=True,
            full_unweighted_flow_comparison_has_no_finite_constant=True,
            retained_filtered_boundary_trace_proved=False,
            actual_Dj_or_Zr_bounded=False, complete_mass_phase_budget_closed=False,
            original_admission=False, nonlinear_existence=False,
            exact_compatibility=False, physical_matching=False,
            validity_noise=False, all_sector_stability=False,
            construction_obstruction=False, goal_complete=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--receipt', default=STEM)
    args = ap.parse_args()
    out = ROOT/'evidence'/(args.receipt+'.json')
    assert not out.exists(), 'Preserve immutable receipts'
    old = json.loads((ROOT/'evidence/current-filtered-j-regularity-test-attempt2.json').read_bytes())
    inputs = list(dict.fromkeys([
        'notes/'+STEM+'.md', 'checks/current_filtered_test_mass_trace.py',
        'evidence/current-filtered-j-regularity-test-attempt2.json',
        'evidence/current-filtered-j-regularity-test-attempt2-verification.json',
        'notes/current-schur-residual-phase-domain.md',
        'notes/recentered-reference-inverse.md',
        'notes/current-global-boundary-moment-return.md',
        'notes/current-seed-full-hessian-realization.md',
        'notes/current-canonical-interior-schur-bound.md',
        'evidence/current-schur-residual-phase-domain.json',
        'evidence/current-schur-residual-phase-domain-review.json',
        'operators/york_mixed_return.py', 'operators/signed_local_core.py',
        *old['input_sha256_after']]))
    before = {p: digest(ROOT/p) for p in inputs}
    protected = json.loads(BASELINE.read_bytes())
    rec = dict(started_utc=datetime.now(timezone.utc).isoformat(), input_sha256_before=before,
        ordinary_argument_utf8=(ROOT/inputs[0]).read_text(encoding='utf-8'),
        arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'))
    try:
        assert all(digest(ROOT/p) == h for p, h in old['input_sha256_after'].items())
        assert all(digest(Path(p)) == h for p, h in protected.items())
        assert digest(Path(sys.executable)) == runtime['python_sha256']
        for pkg in runtime['packages'].values():
            assert all(digest(Path(pkg['root'])/p) == h for p, h in pkg['files'].items())
        phase_review = json.loads((ROOT/'evidence/current-schur-residual-phase-domain-review.json').read_bytes())
        assert all(digest(ROOT/p) == h for p, h in phase_review['scientific_sha256'].items())
        rec['result'] = execute()
        after = {p: digest(ROOT/p) for p in inputs}
        assert after == before
        assert all(digest(Path(p)) == h for p, h in protected.items())
        rec.update(status='full_response_bridge_passes_unweighted_full_range_estimator_rejected',
            input_sha256_after=after, protected_paths_unchanged=len(protected),
            python_sha256=runtime['python_sha256'],
            pinned_packages={k: len(v['files']) for k, v in runtime['packages'].items()})
    except Exception:
        rec.update(status='failed', traceback=traceback.format_exc())
    rec['finished_utc'] = datetime.now(timezone.utc).isoformat()
    out.write_bytes((json.dumps(rec, indent=2)+'\n').encode('utf-8'))
    print(json.dumps(dict(status=rec['status'], error=rec.get('traceback'),
        constants={k: v['ball'] for k, v in rec.get('result', {}).get('constants', {}).items()},
        factored_coefficient=rec.get('result', {}).get('factored_c04')), indent=2))
    if rec['status'] == 'failed':
        raise SystemExit(1)
