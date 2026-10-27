"""Original-norm arithmetic for a signed source map; no regulator certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-signed-weak-source-transfer'
runtime = json.loads((ROOT / 'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/' + STEM + '.md',
    'checks/current_signed_weak_source_transfer.py',
    'notes/current-structured-weak-response-attempt.md',
    'notes/current-curved-paired-heat-bridge.md',
    'notes/current-continuum-energy-dual-estimator-test.md',
    'notes/current-relative-trace-lead-test.md',
    'evidence/current-relative-trace-lead-test.json',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-complete-polynomial-weak-transfer.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(value):
    d = value['decimal_enclosure']
    return arb(int(d['mid']), int(d['rad'])) * arb(10) ** int(d['exponent'])


def enc(value):
    m,r,e = value.mid_rad_10exp(60)
    return {'ball': value.str(35), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def prior(stem):
    run = json.loads((ROOT / 'evidence' / (stem+'.json')).read_bytes())['runs'][-1]
    assert run['status'] == 'passed'
    return run['result']


def execute():
    ctx.prec = 512
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for package in runtime['packages'].values():
        for name,expected in package['files'].items():
            assert digest(Path(package['root']) / name) == expected
    original = prior('current-coupled-truncation-energy-error')
    CL = dec(original['global_shape_L2_constant'])
    CD = dec(original['global_shape_Laplacian_constant'])
    gamma = dec(original['parameters']['gamma'])
    epsilon = dec(prior('current-complete-polynomial-weak-transfer')['constants']['epsilon'])
    assert 0 < epsilon < arb('1.227e-74')
    CE = (CL+CD)/(1-epsilon).sqrt()
    CS = 1/arb('.99').sqrt()+4/gamma
    lead = json.loads((ROOT / INPUTS[6]).read_bytes())
    assert lead['status'] == 'supporting_arithmetic_passed_full_response_unresolved'
    Dquarter = dec(lead['result']['candidate_targets']['one_quarter']
                   ['original_bolt_UV_coefficient_per_A_phi'])
    allocation = dec(lead['result']['parameters']
                     ['hypothetical_UV_allocation_other_unknown_debits_zero'])
    physical = Dquarter*CE
    exponents = {}
    for name,k,s in [('available_global',5,Fraction(3,2)),
                     ('optimistic_C3_grant',3,Fraction(2)),
                     ('required_C3_target',3,Fraction(19,8))]:
        power = 2*(s+1)/k-1
        exponents[name] = {'test_derivatives': k, 'strong_time_power': str(s),
                           'best_of_two_time_power': str(power),
                           'power_integrable_after_angular_weight': power>1}
    assert exponents['available_global']['best_of_two_time_power'] == '0'
    assert exponents['optimistic_C3_grant']['best_of_two_time_power'] == '1'
    assert exponents['required_C3_target']['best_of_two_time_power'] == '5/4'
    assert 2*(Fraction(2)+1)/Fraction(8,3)-1 == Fraction(5,4)
    values = {'C_L_original': CL, 'C_D_original': CD, 'epsilon': epsilon,
        'gamma_original': gamma, 'C_S_L2_to_supremum_upper_expression': CS,
        'C_E_L2_signed_source_to_original_energy_dual': CE,
        'D_quarter_original_UV_error_per_A_phi': Dquarter,
        'original_bolt_UV_error_per_B_phi': physical,
        'conditional_B_phi_ceiling_all_unknown_debits_zero': allocation/physical,
        'hypothetical_UV_allocation_not_free_margin': allocation}
    return {'constants': {k: enc(v) for k,v in values.items()},
        'exact_exponent_tests': exponents,
        'actual_B_phi': None, 'attained_time_exponent': None,
        'actual_signed_g_tau_evaluation': None, 'signed_Zr': None,
        'scope': {'norm_arrow_proved_by_ordinary_note': True,
            'interpolation_estimator_failure_proved_by_ordinary_note': True,
            'full_regulator_removal_certified': False,
            'physical_budget_improvement': 0,
            'conditional_ceiling_is_a_proved_constant': False,
            'construction_obstruction': False,
            'nonlinear_existence_or_stability_proved': False}}


if __name__ == '__main__':
    output = ROOT / 'evidence' / (STEM+'.json')
    assert not output.exists(), 'Preserve immutable receipts.'
    before = {p:digest(ROOT/p) for p in INPUTS}
    receipt = {'started_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256_before': before, 'python': sys.executable,
        'ordinary_argument_utf8': (ROOT/INPUTS[0]).read_text(encoding='utf-8'),
        'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8')}
    try:
        receipt['result'] = execute()
        receipt['status'] = 'supporting_arithmetic_passed_full_remainder_unresolved'
    except BaseException as exc:
        receipt['status'] = 'failed'
        receipt['failure'] = repr(exc)
        raise
    finally:
        receipt['input_sha256_after'] = {p:digest(ROOT/p) for p in INPUTS}
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    assert before == receipt['input_sha256_after']
    print(json.dumps({'status': receipt['status'],
        'constants': {k:v['ball'] for k,v in receipt['result']['constants'].items()},
        'exponents': receipt['result']['exact_exponent_tests']},indent=2))
