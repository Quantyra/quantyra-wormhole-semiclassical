"""S093 checks for the four inputs named by S092.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-scheme-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-scheme-inputs.md'
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
STEP = s.Rational(1, 1000)


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def flat(path):
    return ' '.join(path.read_text(encoding='utf-8').split())


def require_sentence(path, sentence):
    assert sentence in flat(path), sentence


def require_pre_fit_note():
    for sentence in (
        'before any checker run',
        'Assumption. A number for the endpoint value of `c_theta` would be the shipped scheme coefficient `d_theta`',
        '`d_theta` multiplies the normal term, not `I_theta`',
        'Assumption. The unequal-collar invariant would be the product trace `tr S^3`',
        'The product trace is not the collar mismatch',
        'Assumption. An equality setting one `D_j` would solve `Lbar = T + average_j D_j`',
        '`Lbar` does not assign that length',
        'Assumption. A measured `cJ` would be the `B0` derivative of the parallel susceptibility',
        'That slope is not a number',
        'The S085 march is not repeated',
        'The S092 increment is not recomputed as a new result',
        'The S092 collar mismatch is not recomputed as a new result',
        'The S092 trace gap is not recomputed as a new result',
        'The S092 second parallel derivative is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def scheme_coefficient():
    """d_theta is a number and does not assign the endpoint of c_theta."""
    require_sentence(CURVED_NOTE, 'b_theta=4/[5(4pi)^2]')
    require_sentence(CURVED_NOTE, 'd_theta=4/[15(4pi)^2]')
    require_sentence(
        CURVED_NOTE,
        'the finite J_theta choice and the finite c_theta must remain explicit',
    )
    d_theta = s.Rational(4, 15) / (4 * s.pi)**2
    zero(s.factor(d_theta) - 1 / (60 * s.pi**2))
    assert abs(s.N(d_theta, 40) - s.N(1 / (60 * s.pi**2), 40)) < s.Integer(10)**(-20)
    constant, origin = s.symbols('c_theta_mu0 mu0')
    endpoint = constant + s.log(origin) / (20 * s.pi**2)
    difference = s.simplify(endpoint - d_theta)
    assert constant in difference.free_symbols
    assert origin in difference.free_symbols
    assert not difference.is_number
    exact = s.diff(difference, constant)

    def at(value):
        return difference.subs(constant, value)

    center = s.Integer(1)
    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(exact - 1)
    zero(s.expand(sampled - exact))
    declared = joint.charge_identities()
    assert declared['sigma'] == '2*c_theta*Lhat^{uu}'
    signature = list(inspect.signature(rad.march).parameters)
    for name in ('c_theta', 'mu', 'd_theta'):
        assert name not in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        c_theta_set_to_zero=False,
        c_theta_set_equal_to_d_theta=False,
        sigma_replaced_by_d_theta=False,
        d_theta=str(s.factor(d_theta)),
        endpoint_minus_d_theta=str(s.factor(difference)),
        endpoint_is_a_number=False,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def product_trace():
    """The product trace tr(S^3) is not the unequal-collar mismatch."""
    require_sentence(COUPLED_NOTE, 'tr S=tr S^3=0')
    require_sentence(COUPLED_NOTE, 'are not covered by (8).')
    require_sentence(CURVED_NOTE, 'Unequal collar perturbations')
    angle, rate = s.symbols('theta theta_prime')
    chi = s.cos(angle) * coupled.chi0 - s.sin(angle) * coupled.kron(coupled.p1, coupled.I2)
    projector = (coupled.I4 + chi) / 2
    generator = coupled.C * projector * rate / 2
    trace = s.simplify(s.trigsimp(s.expand(s.trace(generator**3))))
    zero(trace)
    point = {angle: s.Rational(1, 2), rate: s.Rational(1, 3)}
    assert abs(s.N(trace.subs(point), 40)) < s.Integer(10)**(-20)
    assert s.simplify(trace - 1) != 0
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        trace_adopted_as_theorem=False,
        product_trace='0',
        trace_equals_mismatch=False,
        recorded_mismatch_recomputed=False,
    )


def mouth_length():
    """Lbar defines an average length and does not assign D_j."""
    require_sentence(MOUTH_NOTE, 'Let Lbar=T+average_j D_j.')
    require_sentence(MOUTH_NOTE, 'average D_j=O(delta log q)')
    clock, average = s.symbols('T Dbar')
    length = clock + average
    exact = s.diff(length, average)
    zero(exact - 1)
    center = s.Integer(1)

    def at(value):
        return length.subs(average, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact.subs(average, center)))
    point = length.subs({clock: 1, average: 1})
    zero(point - 2)
    assert abs(s.N(point, 40) - 2) < s.Integer(10)**(-20)
    assert average in s.simplify(length - clock).free_symbols
    declared = joint.charge_identities()['Lambda']
    rhs = s.sympify(declared.split('=', 1)[1], locals={'N': s.Symbol('N')})
    zero(s.diff(rhs, s.Symbol('D_j')))
    zero(s.diff(rhs, average))
    signature = list(inspect.signature(rad.march).parameters)
    joined = ' '.join(signature)
    assert absent('D_j', joined)
    assert absent('Lbar', joined)
    assert absent('D_j', declared)
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        delta_set_to_zero=False,
        channel_count_assigned=False,
        numeric_D_j=None,
        lbar_derivative_in_average='1',
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def susceptibility_slope():
    """The B0 slope of the parallel susceptibility does not measure cJ."""
    info = em.invariants()
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    susceptibility = s.sympify(
        info['parallel_electric_susceptibility'].replace('^', '**')
    )
    background = s.Symbol('B0')
    coupling = s.Symbol('cJ')
    grav = s.Symbol('kappa')
    unit = susceptibility.subs(grav, 1)
    slope = s.diff(unit, background).subs(background, 1)
    zero(s.expand(slope - 4 * coupling))
    center = s.Integer(1)

    def at(value):
        return unit.subs(background, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - slope))
    assert coupling in slope.free_symbols
    assert not slope.is_number
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
        positive_quadratic_mass_used_as_result=False,
        susceptibility_b0_derivative=str(s.factor(slope)),
        slope_is_a_number=False,
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    scheme = scheme_coefficient()
    trace = product_trace()
    length = mouth_length()
    slope = susceptibility_slope()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S092',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s092_increment_recomputed=False,
        s092_collar_mismatch_recomputed=False,
        s092_trace_gap_recomputed=False,
        s092_second_parallel_derivative_recomputed=False,
        new_lambda_interval=False,
        area_budget_evaluated=False,
        area_budget='not evaluated; no replacement target was written before the test',
        zeta_s_set_to_zero=False,
        sigma_set_to_zero=False,
        c_theta_set_to_zero=False,
        d_theta_set_to_zero=False,
        cJ_set_to_zero=False,
        delta_set_to_zero=False,
        endpoint_energy_cancelled=False,
        positive_quadratic_mass_used_as_result=False,
        tasks=dict(
            scheme_coefficient=scheme,
            product_trace=trace,
            mouth_length=length,
            susceptibility_slope=slope,
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
        scheme_coefficient=scheme['classification'],
        d_theta=scheme['d_theta'],
        product_trace=trace['classification'],
        product_trace_value=trace['product_trace'],
        mouth_length=length['classification'],
        lbar_derivative_in_average=length['lbar_derivative_in_average'],
        susceptibility_slope=slope['classification'],
        susceptibility_b0_derivative=slope['susceptibility_b0_derivative'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
