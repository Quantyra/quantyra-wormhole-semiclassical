"""Actual orbit-invariant lower bound for the direct flow inverse candidate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-flow-inverse-compression-test'
RUNTIME = ROOT/'evidence/current-coupled-mass-adjoint-runtime.json'
runtime = json.loads(RUNTIME.read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/'+STEM+'.md', 'checks/current_flow_inverse_compression_test.py',
    'notes/current-geodesic-two-sided-remainder.md',
    'notes/current-global-geodesic-return.md',
    'notes/current-geodesic-flow-multiplier.md',
    'notes/current-schur-form-preconditioner-test.md',
    'notes/current-local-reference-range-test.md',
    'notes/current-weak-complement-inverse-attempt.md',
    'evidence/current-geodesic-two-sided-remainder.json',
    'evidence/current-structured-source-graph-test.json',
    'evidence/current-quartic-ritz-rank-test-final.json',
    'evidence/current-outer-input-source-test-final.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(value):
    d = value['decimal_enclosure']
    return arb(int(d['mid']), int(d['rad'])) * arb(10)**int(d['exponent'])


def enc(value):
    m, r, e = value.mid_rad_10exp(60)
    return {'ball': value.str(40), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def polyeval(p, x):
    return sum((v*x**k for k, v in p.items()), arb(0))


def execute():
    ctx.prec = 512
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for pkg in runtime['packages'].values():
        for name, want in pkg['files'].items():
            assert digest(Path(pkg['root'])/name) == want
    geo = json.loads((ROOT/INPUTS[8]).read_bytes())['runs'][-1]
    source = json.loads((ROOT/INPUTS[9]).read_bytes())
    rank = json.loads((ROOT/INPUTS[10]).read_bytes())
    outer = json.loads((ROOT/INPUTS[11]).read_bytes())
    assert geo['status'] == rank['status'] == 'passed'
    assert source['status'] == 'source_graph_bound_checked_complete_application_unresolved'
    ep = next(x for x in geo['numerics'] if x['kind'] == 'actual_seed_correlated_endpoint')
    gp = next(x for x in geo['numerics'] if x['kind'] == 'actual_finite_seed_two_sided_parameters')
    Q = {k: dec(v) for k, v in zip([0, 2, 3, 4, 5, 6], ep['Q_coefficients'])}
    Qp = {k-1: k*v for k, v in Q.items() if k}
    beta, C, B2 = Q[0], dec(gp['C']), dec(gp['B2'])
    eps = dec(rank['result']['constants']['epsilon'])
    eq = dec(source['constants']['actual_reference_comparison'])
    fl = dec(gp['F_lambda'])
    assert dec(ep['h']) > 0 and ep['exact_Q_bolt'] == '0'
    assert dec(gp['F_eff']) < 0

    # Outward bounds over the entire fixed interval, never point samples.
    qmin, qmax, q1min, q1max = map(arb, ['.112', '.118', '.244', '.256'])
    left, right = arb('.49'), arb('.51')
    cells = 256
    for i in range(cells):
        lo = left+(right-left)*i/cells
        hi = left+(right-left)*(i+1)/cells
        cell = lo.union(hi)
        qv, qp = polyeval(Q, cell), polyeval(Qp, cell)
        assert qv > qmin and qv < qmax
        assert -qp > q1min and -qp < q1max
    assert qmax < beta
    assert -Q[3] > 0 and Q[4] > 0
    dmin = (-Q[3]-arb(4)/3*Q[4]*right)*left**2
    dmin -= 2*Q[2].abs_upper()*right
    dmin -= sum((k*Q[k].abs_upper()*right**(k-1) for k in [5, 6]), arb(0))
    assert dmin > arb('.078')
    tmax = 4*(right/dmin).sqrt()

    phi2_exact = sum((Fraction(2*(-1)**j*comb(16, j), 2*j+1)
                      for j in range(17)), Fraction(0))
    assert phi2_exact > 0
    phi2 = arb(phi2_exact.numerator)/phi2_exact.denominator
    lphi = 16*(arb(7)/16)**6*(arb(22)/5)
    aphi = arb(3)/20*lphi/(arb(7)/20).sqrt()
    l0 = aphi**2*q1min/(5*arb.pi()*q1max*qmax.sqrt()*tmax)
    k2 = (1-eps)*C/(50*arb.pi()*B2*beta)
    g1 = 100000*beta/((1-eps)*arb('.99')*fl)
    d2 = (1/k2.sqrt()+q1max*(g1/qmin).sqrt())/qmin
    dual_coefficient = phi2*d2**2
    h = arb('1e-16')
    assert arb('.5')-h > left and arb('.5')+h < right
    quadratic_lower = l0*h
    dual_squared_upper = dual_coefficient*h**4
    inverse_norm_lower = quadratic_lower/dual_squared_upper
    full_defect_lower = (1-eq)*inverse_norm_lower-1
    desired_complement_defect = arb('.4')
    inverse_norm_ceiling = (1+desired_complement_defect)/(1-eq)
    required_response_suppression = inverse_norm_lower/inverse_norm_ceiling
    assert l0 > 0 and inverse_norm_lower > 1 and full_defect_lower > desired_complement_defect

    # Current output/source allowances are inputs, not improved errors.
    out_constants = outer['result']['constants']
    assert isinstance(out_constants, dict)
    outer_error = dec(out_constants['attained_outer_S_omission_error'])
    outer_debit = dec(out_constants['attained_outer_FULL_residual_debit'])
    inner_ceiling = dec(out_constants['remaining_inner_application_allowance_NOT_attained'])
    assert (outer_error+inner_ceiling).overlaps(arb('5e-32'))

    z = arb(1)/2
    qhalf = polyeval(Q, z)
    leading_coefficient = C*qhalf**2/(15*arb.pi()*B2*beta)
    finite_part = (z*(C*qhalf).sqrt()).log()-dec(gp['Lscale'])-arb(6)/5
    assert leading_coefficient > 0 and finite_part > 0
    principal_quadratic_upper = h**4*phi2/(leading_coefficient*finite_part)
    witness_response_gain = quadratic_lower/principal_quadratic_upper
    assert witness_response_gain > required_response_suppression
    vals = dict(
        original_beta=beta, qmin=qmin, qmax=qmax, q1min=q1min, q1max=q1max,
        period_difference_lower=dmin, full_crossing_period_upper=tmax,
        Phi2=phi2, phi_second_derivative_lower=lphi, Abel_integral_lower=aphi,
        orbit_invariant_quadratic_coefficient_lower=l0,
        original_energy_fourth_reserve=k2, original_energy_gradient_bound=g1,
        original_energy_to_interior_second_derivative=d2,
        source_dual_squared_coefficient_upper=dual_coefficient,
        h=h, candidate_quadratic_lower=quadratic_lower,
        actual_source_dual_squared_upper=dual_squared_upper,
        candidate_norm_lower=inverse_norm_lower,
        complete_reference_defect_lower=full_defect_lower,
        desired_complement_defect=desired_complement_defect,
        necessary_inverse_norm_ceiling=inverse_norm_ceiling,
        required_response_suppression=required_response_suppression,
        actual_reference_comparison=eq,
        compress_first_interior_log_coefficient=leading_coefficient,
        compress_first_interior_finite_part=finite_part,
        frozen_principal_quadratic_upper=principal_quadratic_upper,
        frozen_principal_witness_response_gain_lower=witness_response_gain,
        preserved_outer_S_omission_error=outer_error,
        preserved_outer_full_residual_debit=outer_debit,
        remaining_S_inner_ceiling_NOT_attained=inner_ceiling,
    )
    return dict(constants={k: enc(v) for k, v in vals.items()},
                exact_Phi2=str(phi2_exact), whole_interval_cells=cells,
                full_flow_inverse_has_no_bounded_energy_dual_extension=True,
                bounded_coarse_or_endpoint_repairs_cannot_cure_divergence=True,
                physical_or_structured_source_counterexample=False,
                complete_structured_application_error_computed=False,
                actual_full_residual_computed=False,
                physical_margin_closed=False,
                conditional_ceiling_promoted_to_constant=False,
                construction_obstruction=False,
                inverse_after_compression_certified=False,
                goal_complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/'+STEM+'.json')
    args = parser.parse_args()
    out = ROOT/args.output
    if out.exists():
        raise SystemExit('Refusing to replace an immutable receipt.')
    before = {p: digest(ROOT/p) for p in INPUTS}
    protected = json.loads(BASELINE.read_bytes())
    run = dict(started_utc=datetime.now(timezone.utc).isoformat(),
               input_sha256_before=before, precision_bits=512,
               frozen_sources_utf8={p: (ROOT/p).read_text(encoding='utf-8') for p in INPUTS[:2]},
               previous_goal_turn='no progress: requested email review supplied no bound')
    try:
        for p, want in protected.items():
            assert digest(Path(p)) == want, p
        run['result'] = execute()
        run['input_sha256_after'] = {p: digest(ROOT/p) for p in INPUTS}
        assert before == run['input_sha256_after']
        assert all(digest(Path(p)) == want for p, want in protected.items())
        run['protected_paths_unchanged'] = len(protected)
        run['status'] = 'candidate_rejected_no_bounded_energy_dual_extension'
    except Exception as exc:
        run['status'] = 'failed'
        run['exception'] = type(exc).__name__+': '+str(exc)
    run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    out.write_text(json.dumps(run, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': run['status'], 'receipt': str(out),
                      'exception': run.get('exception'),
                      'constants': run.get('result', {}).get('constants', {})}, indent=2))
    if run['status'] == 'failed':
        raise SystemExit(1)
