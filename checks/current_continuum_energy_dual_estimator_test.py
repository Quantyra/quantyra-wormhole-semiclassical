"""Only the error-transfer arithmetic; no quantum response is evaluated."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-continuum-energy-dual-estimator-test'
INPUTS = [
    'notes/' + STEM + '.md',
    'checks/current_continuum_energy_dual_estimator_test.py',
    'notes/current-rough-paired-variation.md',
    'notes/current-curved-paired-heat-bridge.md',
    'notes/current-psi-signed-action-factorization-test.md',
    'notes/current-signed-observation-regularization-test.md',
    'notes/current-adjoint-uniform-estimator-test.md',
    'evidence/current-complete-polynomial-weak-transfer.json',
    'evidence/current-signed-observation-regularization-test.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(value):
    d = value['decimal_enclosure']
    return arb(d['mid'] + ' +/- ' + d['rad']) * arb(10) ** int(d['exponent'])


def enc(value):
    m, r, e = value.mid_rad_10exp(60)
    return {'ball': value.str(40), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def execute():
    ctx.prec = 512
    runtime = json.loads((ROOT / INPUTS[-1]).read_bytes())
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for package in runtime['packages'].values():
        for name, expected in package['files'].items():
            assert digest(Path(package['root']) / name) == expected
    original = json.loads((ROOT / INPUTS[-3]).read_bytes())['runs'][-1]
    assert original['status'] == 'passed'
    epsilon = dec(original['result']['constants']['epsilon'])
    observation = json.loads((ROOT / INPUTS[-2]).read_bytes())
    assert observation['status'] == 'unresolved'
    pref = dec(observation['result']['constants']['pref_t'])
    q = arb('.338')
    qt = q / (1 - epsilon)
    ft = qt / (1 + qt)
    coefficient = pref * ft ** 3
    assert 0 < epsilon < arb('1.227e-74')
    assert coefficient > 0
    # Independent exact rearrangement; interval overlap is arithmetic support.
    alternative = pref * (q / (1 + q - epsilon)) ** 3
    assert coefficient.overlaps(alternative)
    assert coefficient < arb('8.1e-45')
    return {
        'physical_third_functional_error_coefficient': enc(coefficient),
        'f_t': enc(ft),
        'physical_prefactor': enc(pref),
        'formula': 'e_bolt,L <= pref_t * f_t^3 * e_L',
        'analytic_result': 'direct smooth-norm-to-energy-dual estimator has infinite supremum',
        'analytic_evidence': 'ordinary argument in notes/' + STEM + '.md, section 3',
        'B': None, 'epsilon_heat': None, 'signed_transfer_Zr': None,
        'scope': {
            'quantum_response_evaluated': False,
            'complete_heat_remainder_bounded': False,
            'physical_budget_improved': False,
            'existence_or_stability_proved': False,
            'actual_response_unbounded_proved': False,
            'construction_obstruction': False,
        },
    }


if __name__ == '__main__':
    path = ROOT / 'evidence' / (STEM + '.json')
    assert not path.exists(), 'Preserve the existing receipt.'
    before = {p: digest(ROOT / p) for p in INPUTS}
    receipt = {
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'unresolved_complete_response',
        'input_sha256_before': before,
        'ordinary_argument_utf8': (ROOT / INPUTS[0]).read_text(encoding='utf-8'),
        'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8'),
        'result': execute(),
    }
    after = {p: digest(ROOT / p) for p in INPUTS}
    assert before == after
    receipt['input_sha256_after'] = after
    receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
    path.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': receipt['status'],
                     'coefficient': receipt['result']['physical_third_functional_error_coefficient']['ball'],
                     'B': None, 'epsilon_heat': None}, indent=2))
