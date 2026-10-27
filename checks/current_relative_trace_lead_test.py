"""Budget arithmetic for the ordinary full-tower trace-hypothesis test.

Does not certify a heat remainder, evaluate a spectrum or prove stability.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-relative-trace-lead-test'
runtime = json.loads((ROOT / 'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/' + STEM + '.md',
    'checks/current_relative_trace_lead_test.py',
    'notes/current-structured-weak-response-attempt.md',
    'notes/current-signed-observation-regularization-test.md',
    'notes/current-weak-coordinate-transfer-test.md',
    'notes/current-canonical-heat-restoration-assembly.md',
    'notes/stationary-paired-dirac-phase.md',
    'evidence/current-coupled-truncation-energy-error.json',
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
    return {'ball': value.str(35), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def execute():
    ctx.prec = 512
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for package in runtime['packages'].values():
        for name, expected in package['files'].items():
            assert digest(Path(package['root']) / name) == expected
    run = json.loads((ROOT / INPUTS[7]).read_bytes())['runs'][-1]
    assert run['status'] == 'passed'
    par = run['result']['parameters']
    beta, C, B2 = [dec(par[k]) for k in ['beta', 'C', 'B2']]
    previous = json.loads((ROOT / INPUTS[8]).read_bytes())['result']
    K3 = dec(previous['physical_third_functional_error_coefficient'])
    K3_upper = arb('8.1e-45')
    assert 0 < K3 < K3_upper
    q = 2000**4
    assert (arb(q)+1)/C > 324
    tau_star = arb('.001')
    gross = arb('1.667e-9')
    known_partial_upper = arb('6.447981645e-11')
    conditional_UV_allocation = gross-known_partial_upper
    normalized_weight = 1/(B2*beta)
    rows = {}
    for name, alpha in [('one_quarter', arb(1)/4),
                        ('one_half', arb(1)/2), ('one', arb(1))]:
        # Complete angular heat bound, no finite angular sample.
        whole_UV = 4*C*tau_star**alpha/alpha
        physical = K3_upper*normalized_weight*whole_UV
        # The decreasing-sum bound includes every positive angular level.
        all_mass_sum_upper = 2*C**(1+alpha)*(
            (q+2)/arb(q+1)**(1+alpha)+1/(alpha*arb(q+1)**alpha))
        full_gamma_upper = (1+alpha).gamma()*all_mass_sum_upper
        N = 1000  # An example tail index, NOT a physical truncation.
        tail_upper = 2*C**(1+alpha)/(alpha*arb(N*(N+q))**alpha)
        rows[name] = {k: enc(v) for k, v in {
            'alpha_candidate_not_proved': alpha,
            'whole_angular_short_time_coefficient_per_A_phi': whole_UV,
            'original_bolt_UV_coefficient_per_A_phi': physical,
            'conditional_A_phi_ceiling_all_unknown_debits_zero':
                conditional_UV_allocation/physical,
            'all_level_mass_power_sum_upper_expression': all_mass_sum_upper,
            'all_level_gamma_weight_upper_expression': full_gamma_upper,
            'tail_above_1000_mass_power_sum_upper_expression': tail_upper,
        }.items()}
    # Exact multiplicity sum in the noncompactness lower bound.
    counts = {str(N): 2*q*N+2*N*(N+1) for N in [1, 1000, 1000000]}
    assert counts['1'] == 2*(q+2)
    assert counts['1000'] == sum(2*(q+2*n) for n in range(1, 1001))
    return {
        'parameters': {k: enc(v) for k,v in {
            'beta': beta, 'C': C, 'B2': B2,
            'actual_prior_K3_interval': K3,
            'K3_upper_used': K3_upper, 'tau_star': tau_star,
            'normalized_angular_weight': normalized_weight,
            'known_partial_debit_upper': known_partial_upper,
            'hypothetical_UV_allocation_other_unknown_debits_zero':
                conditional_UV_allocation,
            'angular_Sp_divergence_coefficient_2C_before_dividing_by_p': 2*C,
        }.items()},
        'candidate_targets': rows,
        'exact_angular_multiplicity_counts': counts,
        'actual_A_phi': None, 'attained_alpha': None,
        'signed_Zr': None, 'actual_phi_energy_norm': None,
        'scope': {
            'ordinary_domain_and_noncompactness_argument_in_note': True,
            'numerical_continuum_proof': False,
            'regulator_removal_certified': False,
            'signed_response_evaluated': False,
            'physical_budget_improvement': 0,
            'conditional_ceiling_is_an_attained_bound': False,
            'construction_obstruction': False,
            'existence_or_stability_proved': False,
        },
    }


if __name__ == '__main__':
    output = ROOT / 'evidence' / (STEM + '.json')
    assert not output.exists(), 'Preserve immutable receipts.'
    before = {p: digest(ROOT / p) for p in INPUTS}
    receipt = {
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256_before': before, 'python': sys.executable,
        'ordinary_argument_utf8': (ROOT / INPUTS[0]).read_text(encoding='utf-8'),
        'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8'),
        'primary_sources': [
            {'url': 'https://arxiv.org/html/2404.18422v3',
             'scope': 'Theorems 1 and 40 and full relative Schatten hypotheses'},
            {'url': 'https://arxiv.org/pdf/1802.06789v1',
             'scope': 'Section 2.3 smooth conformal factor and boundary limitation'},
            {'url': 'https://arxiv.org/pdf/2504.10099v3',
             'scope': 'Section II flat boundary-free scalar-potential setting'},
        ],
    }
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
        output.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    assert before == receipt['input_sha256_after']
    print(json.dumps({'status': receipt['status'],
        'parameters': {k:v['ball'] for k,v in receipt['result']['parameters'].items()},
        'targets': {name: {k:v['ball'] for k,v in row.items()}
                   for name,row in receipt['result']['candidate_targets'].items()}}, indent=2))
