"""Arithmetic support for the ordinary weak-response attempt, not its UV proof."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-structured-weak-response-attempt'
RUNTIME = ROOT / 'evidence/current-coupled-mass-adjoint-runtime.json'
runtime = json.loads(RUNTIME.read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/' + STEM + '.md',
    'checks/current_structured_weak_response_attempt.py',
    'notes/current-rough-paired-variation.md',
    'notes/current-canonical-heat-restoration-assembly.md',
    'notes/current-full-domain-contact-assembly.md',
    'notes/current-angular-resummed-quantum-operator.md',
    'notes/current-coupled-truncation-energy-error.md',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-global-mass-conversion.json',
    'evidence/current-complete-polynomial-weak-transfer.json',
    'evidence/current-continuum-energy-dual-estimator-test.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(value):
    d = value['decimal_enclosure']
    return arb(int(d['mid']), int(d['rad'])) * arb(10) ** int(d['exponent'])


def enc(value):
    m, r, e = value.mid_rad_10exp(60)
    return {'ball': value.str(40), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def prior(stem):
    run = json.loads((ROOT / 'evidence' / (stem + '.json')).read_bytes())['runs'][-1]
    assert run['status'] == 'passed'
    return run['result']


def execute():
    ctx.prec = 512
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for package in runtime['packages'].values():
        for name, expected in package['files'].items():
            assert digest(Path(package['root']) / name) == expected
    a = prior('current-coupled-truncation-energy-error')
    beta, C, B2 = [dec(a['parameters'][k]) for k in ['beta', 'C', 'B2']]
    CL, CD, CI = [dec(a[k]) for k in ['global_shape_L2_constant',
        'global_shape_Laplacian_constant', 'global_shape_supremum_constant']]
    epsilon = dec(prior('current-complete-polynomial-weak-transfer')['constants']['epsilon'])
    assert 0 < epsilon < arb('1.227e-74')
    p = prior('current-global-mass-conversion')['parameters']
    u, seed_eps, seed_tau = [dec(p[k]) for k in ['u', 'eps', 'tau']]
    coeff = [beta, arb(0), u, -4*beta-2*u-seed_tau,
             3*beta+u+seed_eps, -2*seed_eps, seed_eps]
    curvature_upper = sum(j*(j-1)*abs(c) for j, c in enumerate(coeff) if j >= 2)
    assert curvature_upper < 20
    assert 12/arb(1).exp()+192/arb(2).exp() < 32
    assert 12/arb(1).exp() < 32
    q = arb(2000)**4
    assert (q+1)/C > 324
    assert CI < arb('7.477e50')
    tau0 = arb('.001')
    cheap_heat = 64*CI**2/(1-epsilon)
    c0 = 2*CL**2/(arb.pi()*(1-epsilon))
    c2 = (4*20**2*CL**2+16*20*CL*CD+8*CD**2)/(480*arb.pi()*(1-epsilon))
    common = 1/(B2*beta)*(4*C/tau0)/324
    raw = common*1521*cheap_heat/tau0
    area = common*c0/tau0**2
    fourth = common*c2
    positive = raw+area+fourth
    old = json.loads((ROOT / 'evidence/current-continuum-energy-dual-estimator-test.json').read_bytes())
    K3 = dec(old['result']['physical_third_functional_error_coefficient'])
    debit = K3*positive
    allowance = arb('1.667e-9')
    D1 = cheap_heat*2*C/(arb.pi()*B2*beta)
    assert cheap_heat < arb('3.578e103')
    assert positive < arb('2.362e79')
    assert debit < arb('1.910e35')
    assert (area+fourth)/raw < arb('1e-70')
    values = {'curvature_coefficient_sum_upper_expression': curvature_upper,
        'C_heat': cheap_heat, 'C0': c0, 'C2': c2,
        'positive_time_raw_debit': raw, 'positive_time_area_debit': area,
        'positive_time_fourth_debit': fourth, 'B_positive': positive,
        'physical_third_functional_error_coefficient': K3,
        'physical_error_per_unit_phi': debit,
        'unit_phi_upper_expression_over_entire_bolt_allowance': debit/allowance,
        'conditional_phi_norm_threshold': allowance/debit,
        'UV_divergent_majorant_coefficient_L1': D1}
    return {'constants': {k: enc(v) for k,v in values.items()},
        'phi_norm': None, 'A_phi': None, 'alpha_heat': None, 'signed_Zr': None,
        'status': 'positive_time_extension; absolute_UV_estimator_fails; full_response_unresolved',
        'scope': {'positive_time_full_spectral_bound_proved_by_note': True,
            'full_regulator_removal_certified': False, 'signed_quantum_response_evaluated': False,
            'physical_budget_improved': False, 'actual_response_unbounded_proved': False,
            'construction_obstruction': False, 'existence_or_stability_proved': False}}


if __name__ == '__main__':
    output = ROOT / 'evidence' / (STEM + '.json')
    assert not output.exists(), 'Preserve immutable receipts.'
    before = {p: digest(ROOT / p) for p in INPUTS}
    receipt = {'started_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256_before': before, 'python': sys.executable,
        'ordinary_argument_utf8': (ROOT / INPUTS[0]).read_text(encoding='utf-8'),
        'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8')}
    try:
        receipt['result'] = execute()
        receipt['status'] = 'supporting_arithmetic_passed_full_response_unresolved'
    except BaseException as exc:
        receipt['status'] = 'failed'
        receipt['failure'] = repr(exc)
        raise
    finally:
        receipt['input_sha256_after'] = {p: digest(ROOT / p) for p in INPUTS}
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    assert before == receipt['input_sha256_after']
    print(json.dumps({'status': receipt['status'],
        'constants': {k:v['ball'] for k,v in receipt['result']['constants'].items()}}, indent=2))
