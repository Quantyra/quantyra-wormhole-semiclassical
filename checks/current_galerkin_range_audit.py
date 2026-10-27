"""Enclose the finite-rank estimator floor, not a signed quantum evaluation."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-galerkin-range-audit'
runtime = json.loads((ROOT / 'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/' + STEM + '.md', 'checks/current_galerkin_range_audit.py',
    'notes/current-spatial-inverse-output-test.md',
    'evidence/current-spatial-inverse-output-test.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def dec(value):
    d = value['decimal_enclosure']
    return arb(int(d['mid']), int(d['rad'])) * arb(10) ** int(d['exponent'])

def enc(value):
    m, r, e = value.mid_rad_10exp(60)
    return {'ball': value.str(35), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}

def execute():
    ctx.prec = 512
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for pkg in runtime['packages'].values():
        for name, expected in pkg['files'].items():
            assert digest(Path(pkg['root']) / name) == expected
    prior = json.loads((ROOT / INPUTS[3]).read_bytes())
    assert prior['status'] == 'error_estimator_test_passed_signed_action_unresolved'
    c = prior['result']['constants']
    eps = dec(c['epsilon'])
    qt = arb('.338') / (1 - eps)
    ft = qt / (1 + qt)
    e0 = (arb('.0022') + eps) / (1 - eps)
    inverse_floor = 1 / (1 + qt)
    preconditioner_floor = 1 / (1 + qt + e0)
    product_floor = qt * inverse_floor
    W0 = dec(c['W0_complete_scalar_physical_upper'])
    W1 = dec(c['W1_complete_mixed_physical_upper'])
    estimator_floor = 7 * W0 + 3 * W1
    budget = arb('1.667e-9')
    target = arb('1e-30')
    factor = (1 + qt).sqrt()
    assert inverse_floor > arb('.747')
    assert preconditioner_floor > arb('.746')
    assert product_floor.overlaps(ft)
    assert product_floor / target > arb('2.5e29')
    assert estimator_floor / budget > arb('1e26')
    assert factor > 1
    return {'constants': {k: enc(v) for k, v in {
        'retained_inverse_finite_rank_error_lower': inverse_floor,
        'geodesic_inverse_finite_rank_error_lower': preconditioner_floor,
        'separated_F_error_estimator_floor_NOT_actual_error': product_floor,
        'separated_F_estimator_floor_over_proposed_target': product_floor / target,
        'physical_upper_estimator_floor_NOT_actual_error': estimator_floor,
        'physical_upper_estimator_floor_over_whole_allowance': estimator_floor / budget,
        'Ritz_range_comparison_PROVED': factor,
        'required_range_tail_NOT_attained': target / factor,
    }.items()}, 'signed_G_evaluated': False, 'range_tail_evaluated': False,
        'physical_requirements_closed': [],
        'claim': 'specified separated estimator fails; complete signed certificate unresolved'}

if __name__ == '__main__':
    output = ROOT / 'evidence' / (STEM + '.json')
    if output.exists():
        raise SystemExit('Immutable receipt exists; refusing overwrite.')
    before = {p: digest(ROOT / p) for p in INPUTS}
    receipt = {'started_utc': datetime.now(timezone.utc).isoformat(),
               'input_sha256_before': before,
               'ordinary_argument_utf8': (ROOT / INPUTS[0]).read_text(encoding='utf-8'),
               'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8')}
    try:
        receipt['result'] = execute()
        receipt['status'] = 'estimator_failure_enclosed_signed_evaluation_unresolved'
    except BaseException as exc:
        receipt['status'] = 'failed'
        receipt['failure'] = type(exc).__name__
        raise
    finally:
        receipt['input_sha256_after'] = {p: digest(ROOT / p) for p in INPUTS}
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        output.write_bytes((json.dumps(receipt, indent=2) + '\n').encode('utf-8'))
    assert before == receipt['input_sha256_after']
    print(json.dumps({'status': receipt['status'], 'constants': {
        k: v['ball'] for k, v in receipt['result']['constants'].items()}}, indent=2))
