"""S092 checks for the four inputs named by S091.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-endpoint-inputs.md. This checker does not march the
rule-D clock, does not call the RN series, and does not edit the pinned
checkers it calls.
"""
import argparse
import hashlib
import inspect
import json
import platform
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import sympy as s

import coupled_fermion_domain as coupled
import finite_local_electromagnetic_input as em
import source_clock_joint as joint
import source_clock_radiation as rad


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'notes' / 'source-clock-endpoint-inputs.md'
CURVED_NOTE = ROOT / 'notes' / 'curved-fermion-boundary.md'
COUPLED_NOTE = ROOT / 'notes' / 'coupled-fermion-domain.md'
MOUTH_NOTE = ROOT / 'notes' / 'finite-mouth-determinant.md'
EM_NOTE = ROOT / 'notes' / 'finite-local-electromagnetic-input.md'
PINNED = (
    ROOT / 'checks' / 'curved_fermion_boundary.py',
    ROOT / 'checks' / 'coupled_fermion_domain.py',
    ROOT / 'checks' / 'finite_local_electromagnetic_input.py',
    ROOT / 'checks' / 'source_clock_joint.py',
    ROOT / 'checks' / 'source_clock_radiation.py',
)
STENCIL = s.Rational(1, 1000)
DISPLAYED_AXIAL = '2.9452431128603523e-31'


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def flat(path):
    return ' '.join(path.read_text(encoding='utf-8').split())


def require_sentence(path, sentence):
    assert sentence in flat(path), sentence


def require_pre_fit_note():
    for sentence in (
        'before any checker run',
        'Assumption. A number for `c_theta(mu0)` would be the value of the running coupling at one numerical scale',
        'The definite increment between those endpoints is not that value',
        'Assumption. The additional invariants are the quantities the coupled note says vanish on the product',
        'The collar mismatch is not adopted as a two-collar theorem',
        'Assumption. An equality setting `D_j` would be a finite positive root',
        'The trace gap is not that root',
        'Assumption. A measured parallel coefficient would be a number for the second electric derivative',
        'The axial illustration is not that coefficient',
        'The S085 march is not repeated',
        'The S091 running coupling is not recomputed as a new result',
        'The S091 collar ratio is not recomputed as a new result',
        'The S091 trace derivative is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def stencil(at, center):
    samples = [
        at(center + sign * scale * STENCIL)
        for sign, scale in ((-2, 1), (-1, 1), (1, 1), (2, 1))
    ]
    return s.simplify(
        (-samples[3] + 8 * samples[2] - 8 * samples[1] + samples[0]) / (12 * STENCIL)
    )


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def contraction(matrix, gradient):
    trace = s.trace(matrix)
    hat = matrix - trace * s.eye(3) / 3
    return s.expand((gradient.T * hat * gradient)[0])


def scale_increment():
    """The increment between two scales does not assign c_theta(mu0)."""
    require_sentence(CURVED_NOTE, 'b_theta=4/[5(4pi)^2]')
    require_sentence(CURVED_NOTE, 'c_theta(mu) multiplying I_theta has beta function -b_theta')
    curved = flat(CURVED_NOTE)
    assert 'c_theta(1)' not in curved
    assert 'c_theta(mu0)=' not in curved
    beta_magnitude = s.Rational(4, 5) / (4 * s.pi)**2
    zero(s.factor(beta_magnitude) - 1 / (20 * s.pi**2))
    constant, scale, origin = s.symbols('c_theta_mu0 mu mu0', positive=True)
    running = constant - beta_magnitude * s.log(scale / origin)
    increment = s.simplify(running.subs(scale, s.E) - running.subs(scale, 1))
    zero(increment + beta_magnitude)
    assert constant not in increment.free_symbols
    assert origin not in increment.free_symbols
    assert abs(s.N(increment, 40) + s.N(beta_magnitude, 40)) < s.Integer(10)**(-20)
    endpoint = s.simplify(running.subs(scale, 1))
    assert constant in endpoint.free_symbols
    assert not endpoint.is_number
    declared = joint.charge_identities()
    assert declared['sigma'] == '2*c_theta*Lhat^{uu}'
    signature = list(inspect.signature(rad.march).parameters)
    for name in ('c_theta', 'mu'):
        assert name not in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        c_theta_set_to_zero=False,
        integration_constant_set_to_zero=False,
        sigma_replaced_by_increment=False,
        named_numerical_scale_in_note=False,
        increment_from_1_to_E=str(s.factor(increment)),
        endpoint_at_mu_1=str(endpoint),
        endpoint_is_a_number=False,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def collar_mismatch():
    """A nonzero Lhat mismatch is not a two-collar coefficient."""
    info = coupled.analytic()
    assert info['product_Delta_a4_unsmeared'] == 0
    require_sentence(COUPLED_NOTE, 'are not covered by (8).')
    require_sentence(COUPLED_NOTE, 'K_ab=0, E_;r=R_;r=Omega_tau,r=0')
    require_sentence(CURVED_NOTE, 'Unequal collar perturbations')
    require_sentence(CURVED_NOTE, 'Lhat_ab=L_ab-L gamma_ab/3')
    entry = s.symbols('x')
    gradient = s.Matrix([1, 0, 0])
    right = s.diag(0, 1, 0)

    def mismatch(value):
        left = s.diag(value, 0, 0)
        return contraction(left, gradient) - contraction(right, gradient)

    value = s.simplify(mismatch(1))
    zero(value - 1)
    exact = s.diff(mismatch(entry), entry)
    center = s.Integer(1)
    zero(stencil(mismatch, center) - exact.subs(entry, center))
    keys = sorted(info)
    for name in ('K_ab', 'E_r', 'Omega', 'two_collar'):
        assert not any(name in key for key in keys)
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        mismatch_adopted_as_theorem=False,
        positive_quadratic_mass_used_as_result=False,
        mismatch_at_predeclared_collars=str(value),
        product_unsmeared=info['product_Delta_a4_unsmeared'],
        product_return_contains_mismatch=False,
        shipped_return_keys=keys,
    )


def trace_gap():
    """The free-diagonal gap has no finite root that assigns D_j."""
    require_sentence(MOUTH_NOTE, 'Tr(R_0 W)=eta average_j coth(k(T+D_j)/2)/(2k).')
    require_sentence(MOUTH_NOTE, 'Its difference from eta/(2k)')
    eta, wave, clock, length = s.symbols('eta k T D_j', positive=True)
    argument = wave * (clock + length) / 2
    term = eta * s.coth(argument) / (2 * wave)
    free = eta / (2 * wave)
    gap = s.simplify(term - free)
    height = s.symbols('u', positive=True)
    assert s.simplify(s.coth(height) - 1) != 0
    assert s.limit(s.coth(height) - 1, height, s.oo) == 0
    point = {eta: 1, wave: 1, clock: 1, length: 1}
    exact = s.simplify(gap.subs(point))
    # coth stays as 1/tanh until it is rewritten in exponentials.
    reduced = s.simplify(s.expand(exact.rewrite(s.exp)))
    zero(reduced - 1 / (s.exp(2) - 1))
    assert abs(s.N(exact, 40) - s.N(1 / (s.exp(2) - 1), 40)) < s.Integer(10)**(-20)
    assert exact != 0
    declared = joint.charge_identities()['Lambda']
    # N is the channel count in the shipped string. Sympify would
    # otherwise bind that letter to the numerical function sympy.N.
    rhs = s.sympify(declared.split('=', 1)[1], locals={'N': s.Symbol('N')})
    zero(s.diff(rhs, s.Symbol('D_j')))
    signature = list(inspect.signature(rad.march).parameters)
    assert absent('D_j', ' '.join(signature))
    assert absent('D_j', declared)
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        delta_set_to_zero=False,
        channel_count_assigned=False,
        numeric_D_j=None,
        finite_positive_root=False,
        trace_gap_at_unit_point=str(s.factor(reduced)),
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def parallel_coefficient():
    """The axial illustration is not the parallel quadratic coefficient."""
    info = em.invariants()
    require_sentence(EM_NOTE, 'These axial shifts are not polar self-energy')
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    require_sentence(EM_NOTE, DISPLAYED_AXIAL)
    quadratic = s.sympify(info['quadratic_correction'].replace('^', '**'))
    field = s.Symbol('e_parallel')
    coupling = s.Symbol('cJ')
    unit = s.expand(quadratic.subs({s.Symbol('B0'): 1, s.Symbol('kappa'): 1}))
    second = s.diff(unit, field, 2)
    zero(second - 2 * coupling)
    step = s.Rational(1, 1000)
    # A quadratic has an exact second difference.
    sampled = (
        unit.subs(field, step) - 2 * unit.subs(field, 0) + unit.subs(field, -step)
    ) / step**2
    zero(s.expand(sampled - second))
    assert coupling in second.free_symbols
    assert not second.is_number
    level = s.Integer(2000)
    shape = s.sqrt(1 + 12 * s.pi * level**6)
    gamma = 6 * shape / (level**10 * (shape - 1))
    axial = 16 * s.pi * gamma * (1 + 1 / shape)
    displayed = s.Float(DISPLAYED_AXIAL, 40)
    assert abs(s.N(axial, 40) - displayed) / displayed < s.Integer(10)**(-12)
    example = second.subs(coupling, 1)
    assert example == 2
    assert s.simplify(example - axial) != 0
    signature = list(inspect.signature(rad.march).parameters)
    for name in ('cJ', 'chi', 'B0'):
        assert name not in signature
    assert 'kappa' in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        cJ_set_to_zero=False,
        conditional_axial_window_adopted=False,
        example_cJ_adopted_as_measurement=False,
        positive_quadratic_mass_used_as_result=False,
        parallel_second_derivative=str(s.factor(second)),
        parallel_coefficient_is_a_number=False,
        axial_ratio_equals_parallel_coefficient=False,
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    scale = scale_increment()
    collars = collar_mismatch()
    gap = trace_gap()
    magnetic = parallel_coefficient()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S091',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s091_running_coupling_recomputed=False,
        s091_collar_ratio_recomputed=False,
        s091_trace_derivative_recomputed=False,
        new_lambda_interval=False,
        area_budget_evaluated=False,
        area_budget='not evaluated; no replacement target was written before the test',
        zeta_s_set_to_zero=False,
        sigma_set_to_zero=False,
        c_theta_set_to_zero=False,
        cJ_set_to_zero=False,
        delta_set_to_zero=False,
        endpoint_energy_cancelled=False,
        positive_quadratic_mass_used_as_result=False,
        tasks=dict(
            scale_increment=scale,
            collar_mismatch=collars,
            trace_gap=gap,
            parallel_coefficient=magnetic,
        ),
        pinned_checker_sha256={
            path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in PINNED
        },
        versions=dict(python=platform.python_version(), sympy=s.__version__),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n',
    )
    print(json.dumps(dict(
        status=result['status'],
        scale_increment=scale['classification'],
        increment_from_1_to_E=scale['increment_from_1_to_E'],
        collar_mismatch=collars['classification'],
        mismatch_at_predeclared_collars=collars['mismatch_at_predeclared_collars'],
        trace_gap=gap['classification'],
        trace_gap_at_unit_point=gap['trace_gap_at_unit_point'],
        parallel_coefficient=magnetic['classification'],
        parallel_second_derivative=magnetic['parallel_second_derivative'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
