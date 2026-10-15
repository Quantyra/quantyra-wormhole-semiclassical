"""S094 checks for the four inputs named by S093.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-source-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-source-inputs.md'
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
        'Assumption. A number for `c_theta(mu0)` would be the finite coefficient of `integral sqrt(gamma) u^2`.',
        'That source coefficient is not `d_theta`.',
        'Assumption. The unequal-collar invariant would be the product trace `tr(SE)`.',
        'The product trace is not the collar mismatch.',
        'Assumption. A clock number for one `D_j` would be the bound constant `C` in the mouth estimate that multiplies `Lbar`.',
        'Independence is not an assignment.',
        'Assumption. A measured `cJ` would be the second `B0` derivative of the parallel susceptibility at gravitational `kappa = 1`.',
        'That second derivative can stay linear in `cJ`.',
        'The S085 march is not repeated',
        'The S093 scheme coefficient is not recomputed as a new result',
        'The S093 product trace is not recomputed as a new result',
        'The S093 mouth derivative is not recomputed as a new result',
        'The S093 susceptibility slope is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def source_coefficient():
    """The finite source coefficient is not a number for c_theta(mu0)."""
    require_sentence(
        CURVED_NOTE,
        'a dimension-one coefficient of integral sqrt(gamma)u^2',
    )
    require_sentence(
        CURVED_NOTE,
        'The finite coefficient of that source term, the finite J_theta choice and the finite c_theta must remain explicit',
    )
    assert re.search(
        r'source coefficient\s*=\s*-?(?:\d|\.\d)',
        flat(CURVED_NOTE),
    ) is None
    constant, origin = s.symbols('c_theta_mu0 mu0')
    endpoint = constant + s.log(origin) / (20 * s.pi**2)
    exact = s.diff(endpoint, constant)
    zero(exact - 1)
    center = s.Integer(1)

    def at(value):
        return endpoint.subs(constant, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    probe = endpoint.subs(origin, s.exp(1))
    numeric = (
        probe.subs(constant, center + STEP) - probe.subs(constant, center - STEP)
    ) / (2 * STEP)
    assert abs(s.N(numeric - 1, 40)) < s.Integer(10)**(-20)
    assert constant in endpoint.free_symbols
    assert origin in endpoint.free_symbols
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
        c_theta_set_equal_to_d_theta=False,
        source_coefficient_assigned_a_number=False,
        endpoint=str(endpoint),
        endpoint_derivative_in_constant=str(exact),
        endpoint_is_a_number=False,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def product_trace():
    """The product trace tr(SE) is not the unequal-collar mismatch."""
    require_sentence(
        COUPLED_NOTE,
        'tr(SE), tr(S^3), tr(S_:aa), tr((chi\')^2 S) vanish',
    )
    require_sentence(
        COUPLED_NOTE,
        'Unequal/background-dependent collars require their additional invariants and are not covered by (8).',
    )
    angle, rate = s.symbols('theta theta_prime')
    e0, e1, e2, e3 = s.symbols('e0 e1 e2 e3')
    chi = (
        s.cos(angle) * coupled.chi0
        - s.sin(angle) * coupled.kron(coupled.p1, coupled.I2)
    )
    projector = (coupled.I4 + chi) / 2
    generator = coupled.C * projector * rate / 2
    endomorphism = e0 * coupled.I2 + e1 * coupled.p1 + e2 * coupled.p2 + e3 * coupled.p3
    trace = s.trace(coupled.kron(generator, endomorphism))
    product = s.trace(generator) * s.trace(endomorphism)
    zero(s.expand(trace - product))
    reduced = s.simplify(s.trigsimp(s.expand(trace)))
    zero(reduced)
    point = {
        angle: s.Rational(1, 2),
        rate: s.Rational(1, 3),
        e0: s.Rational(2, 5),
        e1: s.Rational(1, 7),
        e2: s.Rational(3, 11),
        e3: s.Rational(4, 13),
    }
    assert abs(s.N(trace.subs(point), 40)) < s.Integer(10)**(-20)
    assert s.simplify(reduced - 1) != 0
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        trace_adopted_as_theorem=False,
        product_trace=str(reduced),
        trace_equals_mismatch=False,
        recorded_mismatch='1',
    )


def bound_constant():
    """The mouth bound constant C does not assign one D_j."""
    require_sentence(
        MOUTH_NOTE,
        '|D^j E_opt| <= C [1+sqrt(a)+a Lbar +a(1+log(1+a/epsilon)+|log mu|+|c|)], j=1,2. (2)',
    )
    require_sentence(
        MOUTH_NOTE,
        'The constant is independent of q, individual D_j, a and epsilon.',
    )
    coefficient, epsilon, mu, sea, clock, average = s.symbols(
        'a epsilon mu c T Dbar'
    )
    length = clock + average
    bracket = (
        1 + s.sqrt(coefficient)
        + coefficient * length
        + coefficient * (
            1 + s.log(1 + coefficient / epsilon)
            + s.Abs(s.log(mu)) + s.Abs(sea)
        )
    )
    exact = s.diff(bracket, average)
    zero(s.simplify(exact - coefficient))
    center = s.Integer(1)

    def at(value):
        return bracket.subs(average, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.simplify(sampled - exact))
    values = {
        coefficient: s.Rational(1, 4),
        epsilon: s.Rational(1, 5),
        mu: s.Integer(2),
        sea: s.Integer(3),
        clock: s.Integer(4),
    }
    numeric_bracket = bracket.subs(values)

    def numeric_at(value):
        return numeric_bracket.subs(average, value)

    numeric = (
        numeric_at(center + STEP) - numeric_at(center - STEP)
    ) / (2 * STEP)
    assert abs(s.N(numeric - values[coefficient], 40)) < s.Integer(10)**(-20)
    bound = s.Symbol('C')
    weighted = s.diff(bound * bracket, average)
    zero(s.simplify(weighted - bound * coefficient))
    assert bound in weighted.free_symbols
    assert not bound.is_number
    assert not weighted.is_number
    declared = joint.charge_identities()['Lambda']
    rhs = s.sympify(declared.split('=', 1)[1], locals={'N': s.Symbol('N')})
    clock_symbol = s.Symbol('D_j')
    zero(s.diff(rhs, clock_symbol))
    zero(s.simplify(rhs.subs(clock_symbol, 1) - rhs.subs(clock_symbol, 2)))
    signature = list(inspect.signature(rad.march).parameters)
    joined = ' '.join(signature)
    assert absent('D_j', joined)
    assert absent('Lbar', joined)
    assert absent('C', joined)
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        delta_set_to_zero=False,
        channel_count_assigned=False,
        numeric_D_j=None,
        bound_constant_is_a_number=False,
        bracket_derivative_in_average=str(exact),
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def second_susceptibility_derivative():
    """The second B0 derivative does not measure cJ."""
    info = em.invariants()
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    susceptibility = s.sympify(
        info['parallel_electric_susceptibility'].replace('^', '**')
    )
    background = s.Symbol('B0')
    coupling = s.Symbol('cJ')
    grav = s.Symbol('kappa')
    second = s.diff(susceptibility, background, 2)
    zero(s.expand(second - 4 * coupling * grav**4))
    assert background not in second.free_symbols
    at_unit = second.subs(grav, 1)
    zero(s.expand(at_unit - 4 * coupling))
    unit = susceptibility.subs(grav, 1)
    center = s.Integer(1)

    def at(value):
        return unit.subs(background, value)

    sampled = (
        at(center + STEP) - 2 * at(center) + at(center - STEP)
    ) / STEP**2
    zero(s.expand(sampled - at_unit))
    numeric_unit = unit.subs(coupling, s.Rational(3, 7))

    def numeric_at(value):
        return numeric_unit.subs(background, value)

    numeric = (
        numeric_at(center + STEP) - 2 * numeric_at(center) + numeric_at(center - STEP)
    ) / STEP**2
    assert abs(s.N(numeric - 4 * s.Rational(3, 7), 40)) < s.Integer(10)**(-18)
    assert coupling in at_unit.free_symbols
    assert not at_unit.is_number
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
        second_b0_derivative=str(s.factor(second)),
        second_b0_derivative_at_gravitational_kappa_1=str(s.factor(at_unit)),
        second_derivative_is_a_number=False,
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    source = source_coefficient()
    trace = product_trace()
    constant = bound_constant()
    slope = second_susceptibility_derivative()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S093',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s093_scheme_coefficient_recomputed=False,
        s093_product_trace_recomputed=False,
        s093_mouth_derivative_recomputed=False,
        s093_susceptibility_slope_recomputed=False,
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
            source_coefficient=source,
            product_trace=trace,
            bound_constant=constant,
            second_b0_derivative=slope,
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
        source_coefficient=source['classification'],
        endpoint_derivative_in_constant=source['endpoint_derivative_in_constant'],
        endpoint_is_a_number=source['endpoint_is_a_number'],
        product_trace=trace['classification'],
        product_trace_value=trace['product_trace'],
        bound_constant=constant['classification'],
        bracket_derivative_in_average=constant['bracket_derivative_in_average'],
        bound_constant_is_a_number=constant['bound_constant_is_a_number'],
        second_b0_derivative=slope['classification'],
        second_b0_derivative_at_gravitational_kappa_1=slope[
            'second_b0_derivative_at_gravitational_kappa_1'
        ],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
