"""S095 checks for the four inputs named by S094.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-choice-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-choice-inputs.md'
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
        'Assumption. A number for `c_theta(mu0)` would be the finite `J_theta` choice.',
        '`J_theta` is not that source coefficient.',
        'Assumption. The unequal-collar invariant would be the product trace `tr((chi\')^2 S)`.',
        'That trace is not `tr(SE)`.',
        'Assumption. A number for one `D_j` would be the local J derivative bound `O(a+a epsilon[1+|log(lambda_0/epsilon)|])`.',
        'A big-O statement is not an assignment.',
        'Assumption. A measured `cJ` would be the `cJ` derivative of the second `B0` derivative of the parallel susceptibility at gravitational `kappa = 1`.',
        'Removing `cJ` by differentiation does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S094 source coefficient is not recomputed as a new result',
        'The S094 product trace is not recomputed as a new result',
        'The S094 bracket derivative is not recomputed as a new result',
        'The S094 second `B0` derivative is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def j_theta_choice():
    """The finite J_theta choice does not assign c_theta(mu0)."""
    require_sentence(CURVED_NOTE, 'J_theta=integral sqrt(gamma)L u^2')
    require_sentence(CURVED_NOTE, 'Its coefficient is scheme dependent.')
    assert re.search(r'J_theta\s*=\s*-?(?:\d|\.\d)', flat(CURVED_NOTE)) is None
    scheme, functional = s.symbols('d_theta J_theta')
    term = -scheme * functional / 3
    exact = s.diff(term, functional)
    zero(s.simplify(exact + scheme / 3))
    center = s.Integer(1)

    def at(value):
        return term.subs(functional, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    assert scheme in exact.free_symbols
    assert functional not in exact.free_symbols
    assert not exact.is_number
    constant, origin = s.symbols('c_theta_mu0 mu0')
    endpoint = constant + s.log(origin) / (20 * s.pi**2)
    assert constant in endpoint.free_symbols
    assert not endpoint.is_number
    declared = joint.charge_identities()
    assert declared['sigma'] == '2*c_theta*Lhat^{uu}'
    signature = list(inspect.signature(rad.march).parameters)
    for name in ('J_theta', 'c_theta'):
        assert name not in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        c_theta_set_to_zero=False,
        c_theta_set_equal_to_d_theta=False,
        j_theta_assigned_a_number=False,
        j_theta_derivative=str(s.factor(exact)),
        derivative_assigns_c_theta=False,
        endpoint_keeps_c_theta_mu0=True,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def curvature_trace():
    """The product trace tr((chi')^2 S) is not the unequal-collar mismatch."""
    require_sentence(
        COUPLED_NOTE,
        "tr(SE), tr(S^3), tr(S_:aa), tr((chi')^2 S) vanish",
    )
    require_sentence(
        COUPLED_NOTE,
        'Unequal/background-dependent collars require their additional invariants and are not covered by (8).',
    )
    angle, rate = s.symbols('theta theta_prime')
    chi = (
        s.cos(angle) * coupled.chi0
        - s.sin(angle) * coupled.kron(coupled.p1, coupled.I2)
    )
    projector = (coupled.I4 + chi) / 2
    generator = coupled.C * projector * rate / 2
    chip = s.diff(chi, angle) * rate
    trace = s.simplify(s.trigsimp(s.expand(s.trace(chip * chip * generator))))
    zero(trace)
    point = {angle: s.Rational(1, 2), rate: s.Rational(1, 3)}
    assert abs(s.N(s.trace(chip * chip * generator).subs(point), 40)) < s.Integer(10)**(-20)
    assert s.simplify(trace - 1) != 0
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        trace_adopted_as_theorem=False,
        product_trace=str(trace),
        trace_equals_mismatch=False,
        recorded_mismatch='1',
    )


def local_bound():
    """The local J derivative bound does not assign one D_j."""
    require_sentence(
        MOUTH_NOTE,
        'O(a+a epsilon[1+|log(lambda_0/epsilon)|])',
    )
    coefficient, epsilon, scale = s.symbols('a epsilon lambda_0', positive=True)
    expression = coefficient + coefficient * epsilon * (
        1 + s.Abs(s.log(scale / epsilon))
    )
    assert not expression.is_number
    exact = s.diff(expression, scale)
    center = s.Integer(2)
    step = s.Symbol('h', positive=True)
    quotient = (
        expression.subs(scale, center + step) - expression.subs(scale, center - step)
    ) / (2 * step)
    limit = s.limit(quotient, step, 0)
    zero(s.simplify(limit - exact.subs(scale, center)))
    values = {coefficient: s.Rational(1, 4), epsilon: s.Integer(1)}
    numeric_expression = expression.subs(values)

    def numeric_at(value):
        return numeric_expression.subs(scale, value)

    numeric_step = s.Rational(1, 10**8)
    numeric = (
        numeric_at(center + numeric_step) - numeric_at(center - numeric_step)
    ) / (2 * numeric_step)
    expected = exact.subs(values).subs(scale, center)
    assert abs(s.N(numeric - expected, 40)) < s.Integer(10)**(-12)
    declared = joint.charge_identities()['Lambda']
    rhs = s.sympify(declared.split('=', 1)[1], locals={'N': s.Symbol('N')})
    clock_symbol = s.Symbol('D_j')
    zero(s.diff(rhs, clock_symbol))
    zero(s.simplify(rhs.subs(clock_symbol, 1) - rhs.subs(clock_symbol, 2)))
    signature = list(inspect.signature(rad.march).parameters)
    joined = ' '.join(signature)
    for name in ('D_j', 'a', 'epsilon', 'lambda_0'):
        assert absent(name, joined)
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        delta_set_to_zero=False,
        channel_count_assigned=False,
        big_o_evaluated_as_equality=False,
        numeric_D_j=None,
        bound_is_a_number=False,
        lambda_0_derivative=str(s.factor(exact)),
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def coupling_derivative():
    """Differentiating the second B0 slope in cJ does not measure cJ."""
    info = em.invariants()
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    susceptibility = s.sympify(
        info['parallel_electric_susceptibility'].replace('^', '**')
    )
    background = s.Symbol('B0')
    coupling = s.Symbol('cJ')
    grav = s.Symbol('kappa')
    second = s.diff(susceptibility, background, 2).subs(grav, 1)
    exact = s.diff(second, coupling)
    zero(s.expand(exact - 4))
    assert coupling not in exact.free_symbols
    assert exact.is_number
    center = s.Integer(1)

    def at(value):
        return second.subs(coupling, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    numeric = (
        at(s.Rational(3, 7) + STEP) - at(s.Rational(3, 7) - STEP)
    ) / (2 * STEP)
    assert abs(s.N(numeric - exact, 40)) < s.Integer(10)**(-18)
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
        cJ_derivative=str(s.factor(exact)),
        derivative_is_a_value_of_cJ=False,
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    choice = j_theta_choice()
    trace = curvature_trace()
    bound = local_bound()
    slope = coupling_derivative()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S094',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s094_source_coefficient_recomputed=False,
        s094_product_trace_recomputed=False,
        s094_bracket_derivative_recomputed=False,
        s094_second_b0_derivative_recomputed=False,
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
            j_theta_choice=choice,
            curvature_trace=trace,
            local_bound=bound,
            coupling_derivative=slope,
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
        j_theta_choice=choice['classification'],
        j_theta_derivative=choice['j_theta_derivative'],
        curvature_trace=trace['classification'],
        product_trace_value=trace['product_trace'],
        local_bound=bound['classification'],
        bound_is_a_number=bound['bound_is_a_number'],
        coupling_derivative=slope['classification'],
        cJ_derivative=slope['cJ_derivative'],
        derivative_is_a_value_of_cJ=slope['derivative_is_a_value_of_cJ'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
