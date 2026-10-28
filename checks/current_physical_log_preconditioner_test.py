"""Outward arithmetic supporting the physical-log whole-defect rejection."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import argparse
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_quartic_ritz_rank_test as h

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-physical-log-preconditioner-test'
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS = [
    'notes/'+STEM+'.md', 'checks/current_physical_log_preconditioner_test.py',
    'checks/current_quartic_ritz_rank_test.py',
    'notes/recentered-reference-inverse.md',
    'notes/current-seed-full-hessian-realization.md',
    'notes/stationary-weighted-tt-inverse.md',
    'notes/current-quartic-ritz-rank-test.md',
    'notes/current-global-scalar-log-comparison.md',
    'notes/current-structured-source-graph-test.md',
    'evidence/current-quartic-ritz-rank-test-final.json',
    'evidence/current-finite-local-action-assembly.json',
    'evidence/current-global-scalar-log-comparison-accepted.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def execute():
    h.ctx.prec = 512
    a = h.arb
    assert digest(Path(sys.executable)) == h.runtime['python_sha256']
    for pkg in h.runtime['packages'].values():
        for name, expected in pkg['files'].items():
            assert digest(Path(pkg['root'])/name) == expected
    accepted = json.loads((ROOT/'evidence/current-quartic-ritz-rank-test-final.json').read_bytes())
    assert accepted['status'] == 'passed'
    for path, expected in accepted['input_sha256_before'].items():
        assert digest(ROOT/path) == expected
    assert accepted['input_sha256_before'] == accepted['input_sha256_after']
    scalar_log = json.loads((ROOT/'evidence/current-global-scalar-log-comparison-accepted.json').read_bytes())
    for path, expected in scalar_log['input_sha256_before'].items():
        assert digest(ROOT/path) == expected

    local = h.receipt('current-finite-local-action-assembly')
    ep = h.pick(local, 'actual_seed_correlated_endpoint')
    assert ep['exact_Q_bolt'] == '0'
    end = h.dec(ep['z_h'])
    Q = [a(0)]*7
    for k, v in zip([0, 2, 3, 4, 5, 6], ep['Q_coefficients']):
        Q[k] = h.dec(v)
    beta = Q[0]
    psi = [a(0)]*4+[a(1)]
    # Exact York source vanishes for psi=z^4, J=0.
    assert all(c.is_zero() for c in h.add(h.scale(psi, 4), h.scale(h.shift(h.diff(psi), 1), -1)))
    Lpsi = h.scale(h.add(
        h.shift(h.mul(Q, h.diff(h.diff(psi))), 2),
        h.shift(h.mul(h.diff(Q), h.diff(psi)), 2),
        h.scale(h.shift(h.mul(Q, h.diff(psi)), 1), -2)), -1/beta)
    m1 = h.primitive_at(h.shift(h.mul(Q, h.mul(h.diff(psi), h.diff(psi))), -2), end)/beta
    m3 = h.primitive_at(h.shift(h.mul(Q, h.mul(h.diff(Lpsi), h.diff(Lpsi))), -2), end)/beta
    # Independent moment identity <psi,Lpsi>, with exact zero bolt flux.
    m1_pairing = h.primitive_at(h.shift(h.mul(psi, Lpsi), -4), end)
    assert m1.overlaps(m1_pairing) and m1 > 0 and m3 > 0
    r = accepted['result']
    known = h.dec(r['full_quantum_energy_upper_rows'][1]['known_energy_diagonal'])
    b = h.dec(r['variational_B_lower_matrix'][1][1])
    full_lower = (1-a('.0022'))*(1-a('1.227e-74'))/(1+a('.0022'))*known+b
    quantum_upper = 3*a('1e-33')/(20*a.pi())*(400*m1+14*m3)
    ref_upper = a('1.5')*(1+a('1e-30'))*m1+quantum_upper
    ratio = full_lower/ref_upper
    defect_lower = ratio-1
    correction_lower = full_lower*a(2)/3-ref_upper
    required_factor_lower = full_lower*a(2)/(3*ref_upper)
    gap = a('1.5')*(1-a('1.179e-8'))-60*a('1e-33')/a.pi()
    H_inverse_upper = 8/gap
    alternative_transfer = a('.49')+a('.000654')*a('1.49')/a('.997216')
    assert H_inverse_upper < a('5.34')
    assert full_lower > a('2.516441223600354')
    assert ref_upper < a('0.800000000011634')
    assert defect_lower > a('2.14555152945470')
    assert correction_lower > a('.87762748238860')
    assert required_factor_lower > a('2.09703435296980')
    assert alternative_transfer < a('.491')
    values = dict(m1=m1, m3=m3, known_complete_local_energy=known,
        full_scalar_variational_lower=b, complete_A_Q_energy_lower=full_lower,
        reference_quantum_energy_upper=quantum_upper, reference_energy_upper=ref_upper,
        generalized_Rayleigh_lower=ratio, original_dual_defect_lower=defect_lower,
        necessary_reference_energy_addition=correction_lower,
        necessary_reference_energy_factor=required_factor_lower,
        reference_H_inverse_upper=H_inverse_upper,
        conditional_alternative_transferred_defect=alternative_transfer)
    return {'precision_bits': h.ctx.prec, 'constants': {k:h.enc(v) for k,v in values.items()},
        'scope': {'uniform_exact_reference_rejected': True, 'both_original_equations': True,
                  'complete_quantum_comparison': True, 'finite_mode_replacement': False,
                  'implemented_inverse': False, 'actual_structured_source_error_computed': False,
                  'physical_margin_closed': False, 'construction_obstruction': False,
                  'goal_complete': False},
        'unknown': {'implemented_P0_error': None, 'S_application_error': None,
                    'K_application_error': None, 'signed_G_t': None,
                    'physical_and_graph_output_comparisons': None}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/'+STEM+'.json')
    args = parser.parse_args()
    output = ROOT/args.output
    assert not output.exists(), 'Receipts are immutable; choose a fresh path.'
    protected = json.loads(BASELINE.read_bytes())
    verify = lambda: all(digest(Path(p)) == v for p,v in protected.items())
    assert verify()
    before = {p:digest(ROOT/p) for p in INPUTS}
    data = {'started_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256_before': before, 'ordinary_argument_utf8': (ROOT/INPUTS[0]).read_text(encoding='utf-8'),
        'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8')}
    try:
        data['result'] = execute()
        data['status'] = 'specified_uniform_estimator_rejected_full_goal_unresolved'
    except Exception as exc:
        data['status'] = 'failed'
        data['error'] = repr(exc)
        raise
    finally:
        data['input_sha256_after'] = {p:digest(ROOT/p) for p in INPUTS}
        data['protected_paths_unchanged'] = verify()
        data['protected_path_count'] = len(protected)
        data['finished_utc'] = datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
    assert before == data['input_sha256_after'] and data['protected_paths_unchanged']
    for k,v in data['result']['constants'].items():
        print(k+': '+v['ball'])
    print('Receipt: '+str(output))


if __name__ == '__main__':
    main()
