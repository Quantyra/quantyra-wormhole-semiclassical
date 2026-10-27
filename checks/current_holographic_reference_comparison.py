"""Physical bolt budget of an exact LLL reference; no signed massive solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-holographic-reference-comparison'
OWN = ['notes/' + STEM + '.md', 'checks/current_holographic_reference_comparison.py']
INPUTS = OWN + [
    'notes/stationary-paired-dirac-phase.md',
    'notes/coupled-fermion-domain.md',
    'notes/current-canonical-heat-restoration-assembly.md',
    'notes/current-physical-york-dual-assembly.md',
    'notes/current-angular-resummed-quantum-operator.md',
    'notes/current-admission-bottleneck-estimator-test.md',
    'notes/current-joint-radius-composition-test.md',
    'notes/current-psi-signed-action-factorization-test.md',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-joint-radius-composition-test.json',
    'evidence/current-admission-bottleneck-estimator-test.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(x):
    if not isinstance(x, dict):
        return arb(x.strip('[]'))
    d = x['decimal_enclosure']
    return arb(d['mid'] + ' +/- ' + d['rad']) * arb(10) ** int(d['exponent'])


def enc(x):
    m, r, e = x.mid_rad_10exp(50)
    return {'ball': x.str(35), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def receipt(stem):
    run = json.loads((ROOT / 'evidence' / (stem + '.json')).read_bytes())['runs'][-1]
    assert run['status'] == 'passed'
    return run['result']


def execute():
    ctx.prec = 512
    rt = json.loads((ROOT / 'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
    assert digest(Path(sys.executable)) == rt['python_sha256']
    for pkg in rt['packages'].values():
        for p, h in pkg['files'].items():
            assert digest(Path(pkg['root']) / p) == h
    energy = receipt('current-coupled-truncation-energy-error')
    joint = receipt('current-joint-radius-composition-test')['constants']
    admission = receipt('current-admission-bottleneck-estimator-test')['results']
    B2 = dec(energy['parameters']['B2'])
    CD = dec(energy['global_shape_Laplacian_constant'])
    Kret = dec(admission['Kret_original_mixed_point_conversion'])
    CJ = dec(joint['attained_joint_bolt_operator_upper'])
    paff = dec(joint['p_aff'])
    qop = arb(169) / 500
    f = qop / (1 + qop)
    cL = arb(2000) ** 4 / (6 * arb.pi() * B2)
    KL = Kret * cL * CD * f ** 3
    reference = paff * KL
    old = paff * CJ
    remainder = paff * (CJ + KL)
    budget = arb('1.667e-9')
    assert old.overlaps(dec(joint['attained_joint_physical_upper']))
    assert reference > 0 and reference < old
    assert reference > budget and remainder > budget
    assert KL / CJ < arb('1e-20')
    return {
        'precision_bits': 512,
        'outcome': 'finite_excessive_bolt_remainder_estimator_full_A_comparison_unresolved',
        'constants': {k: enc(v) for k, v in {
            'B_radius_squared': B2, 'c_L_in_original_normalization': cL,
            'accepted_CD': CD, 'accepted_Kret': Kret,
            'reference_bolt_operator_upper_KL': KL,
            'accepted_complete_bolt_operator_upper_CJ': CJ,
            'original_P32_affine_prefactor': paff,
            'reference_P32_bolt_upper': reference,
            'reference_upper_over_whole_bolt_allowance': reference / budget,
            'reference_upper_fraction_of_complete_upper': KL / CJ,
            'complete_remainder_P32_bolt_upper': remainder,
            'remainder_estimator_overrun': remainder / budget,
            'whole_bolt_allowance_minus_remainder_upper': budget - remainder,
            'optimistic_upper_expression_subtraction_not_a_lower_bound': old - reference,
            'optimistic_subtraction_overrun_not_actual_correction': (old - reference) / budget,
        }.items()},
        'analytic_results': {
            'constant_source_angle_derivative_on_identical_disks': 'exactly zero',
            'reference_response': 'exact existing kernel-free LLL Polyakov cocycle',
            'fourth_order_logarithmic_coefficient_removed': 'exactly zero',
            'normalized_unfiltered_remainder_quadratic_sequence_limit': 1,
            'restricted_F3_or_affine_map_obstruction': False,
        },
        'errors_and_unknowns': {
            'LLL_reference_response_model_error_in_fixed_chart': 0,
            'arithmetic_error': 'included in interval enclosures',
            'signed_massive_action_and_dual_norm_error': None,
            'full_A_remainder_comparison_and_phase_compatibility': None,
            'original_Y_comparison': None,
            'actual_mass_column': None,
            'other_coupled_columns_and_application_errors': None,
            'other_bolt_debits': None,
        },
        'scope': {
            'full_angular_or_thermal_spectrum_truncated': False,
            'reference_replaces_actual_boundary_domain': False,
            'original_radius_cost_retained': True,
            'signed_massive_response_evaluated': False,
            'new_physical_allowance_closed': False,
            'complete_original_A_norm_certified': False,
            'construction_obstruction': False,
            'goal_complete': False,
        },
    }


if __name__ == '__main__':
    suffix = '-replay' if '--replay' in sys.argv else ''
    out = ROOT / 'evidence' / (STEM + suffix + '.json')
    assert not out.exists(), 'Preserve prior execution receipts'
    hashes = {p: digest(ROOT / p) for p in INPUTS}
    run = {'started_utc': datetime.now(timezone.utc).isoformat(),
           'input_sha256_before': hashes,
           'frozen_sources_utf8': {p: (ROOT / p).read_text(encoding='utf-8') for p in OWN}}
    run['result'] = execute()
    run['input_sha256_after'] = {p: digest(ROOT / p) for p in INPUTS}
    assert run['input_sha256_after'] == hashes
    run['arithmetic_status'] = 'completed_with_outward_intervals'
    run['scientific_status'] = 'remainder_estimator_excessive_full_physical_comparison_unresolved'
    run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    out.write_bytes((json.dumps(run, indent=2) + '\n').encode())
    print(json.dumps({'scientific_status': run['scientific_status'],
        'results': {k: v['ball'] for k, v in run['result']['constants'].items()
                    if k not in ['B_radius_squared', 'accepted_CD', 'accepted_Kret']},
        'scope': run['result']['scope']}, indent=2))
