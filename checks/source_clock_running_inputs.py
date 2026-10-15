"""S091 checks for the four inputs named by S090.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-running-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-running-inputs.md'
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


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def flat(path):
    return ' '.join(path.read_text(encoding='utf-8').split())


def require_sentence(path, sentence):
    assert sentence in flat(path), sentence


def require_pre_fit_note():
    for sentence in (
        'before any checker run',
        'Assumption. A renormalization condition is an equality that assigns one number to `c_theta',
        'The integration constant is not set to zero',
        'Assumption. The compared coefficient is the cross term between the sum of two equation-(8) integrals and the substituted square',
        'That ratio is not adopted as a two-collar theorem',
        'Assumption. An equality assigning `D_j` would be an equation `D_j`',
        'The trace derivative is not that equality',
        'Assumption. A measurement of `cJ` would be a number returned at `B0 = 1`',
        'The vanishing first variation is not that number',
        'The S085 march is not repeated',
        'The S090 substitution difference is not recomputed as a new result',
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


def running_coupling():
    """The beta function leaves the integration constant unfixed."""
    require_sentence(CURVED_NOTE, 'b_theta=4/[5(4pi)^2]')
    require_sentence(
        CURVED_NOTE,
        'c_theta(mu) multiplying I_theta has beta function -b_theta',
    )
    require_sentence(CURVED_NOTE, 'finite c_theta must remain explicit')
    beta_magnitude = s.Rational(4, 5) / (4 * s.pi)**2
    zero(s.factor(beta_magnitude) - 1 / (20 * s.pi**2))
    constant, scale, origin = s.symbols('c_theta_mu0 mu mu0', positive=True)
    running = constant - beta_magnitude * s.log(scale / origin)
    exact = s.diff(running, scale)
    zero(s.simplify(scale * exact + beta_magnitude))
    step = s.symbols('h', positive=True)
    quotient = scale * (
        running.subs(scale, scale + step) - running.subs(scale, scale - step)
    ) / (2 * step)
    zero(s.limit(quotient, step, 0) - scale * exact)
    base = {constant: s.Rational(1, 7), origin: 1}
    scale_value = s.Integer(2)
    step_value = s.Rational(1, 10**6)
    numeric = s.N(
        scale_value * (
            running.subs(base).subs(scale, scale_value + step_value)
            - running.subs(base).subs(scale, scale_value - step_value)
        ) / (2 * step_value),
        40,
    )
    assert abs(numeric + s.N(beta_magnitude, 40)) < s.Integer(10)**(-12)
    assert constant in running.free_symbols
    assert not running.is_number
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
        sigma_replaced_by_running_coupling=False,
        b_theta=str(s.factor(beta_magnitude)),
        beta_function=str(s.factor(-beta_magnitude)),
        running_coupling='c_theta(mu0) - b_theta*log(mu/mu0)',
        integration_constant='c_theta(mu0)',
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def collar_ratio():
    """The cross-to-diagonal ratio is not a two-collar theorem."""
    info = coupled.analytic()
    normal = s.sympify(info['product_Delta_a4_normal_smearing_without_4pi'])
    angle = s.Symbol('theta_prime')
    constant = s.simplify(normal.subs(angle, 1))
    zero(normal - constant * angle**2)
    left, right = s.symbols('theta_L theta_R')
    diagonal = constant * (left**2 + right**2)
    substituted = s.expand(constant * (left - right)**2)
    cross = s.expand(substituted - diagonal)
    zero(cross + 2 * constant * left * right)
    ratio = s.simplify(cross / diagonal)
    point = {left: 1, right: s.Rational(1, 2)}
    zero(ratio.subs(point) + s.Rational(4, 5))
    held = {right: s.Rational(1, 2)}
    center = s.Rational(1, 3)
    exact = s.diff(cross, left).subs(held).subs(left, center)

    def at(value):
        return cross.subs(held).subs(left, value)

    zero(stencil(at, center) - exact)
    keys = sorted(info)
    assert not any('bound' in key or 'two_collar' in key for key in keys)
    require_sentence(COUPLED_NOTE, 'are not covered by (8).')
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        ratio_adopted_as_theorem=False,
        positive_quadratic_mass_used_as_result=False,
        ratio_at_predeclared_point=str(s.factor(ratio.subs(point))),
        shipped_return_keys=keys,
        two_collar_error_bound_returned=False,
        error_bound_proved=False,
    )


def exterior_length():
    """A nonzero trace derivative does not assign D_j."""
    require_sentence(MOUTH_NOTE, 'arbitrary fixed D_j>0')
    require_sentence(MOUTH_NOTE, 'average D_j=O(delta log q)')
    eta, wave, clock, length = s.symbols('eta k T D_j', positive=True)
    term = eta * s.coth(wave * (clock + length) / 2) / (2 * wave)
    exact = s.diff(term, length)
    expected = -eta * s.csch(wave * (clock + length) / 2)**2 / 4
    zero(s.simplify(exact - expected))
    step = s.symbols('h', positive=True)
    quotient = (
        term.subs(length, length + step) - term.subs(length, length - step)
    ) / (2 * step)
    zero(s.limit(quotient, step, 0) - exact)
    point = {eta: 1, wave: 1, clock: 1, length: 1}
    zero(exact.subs(point) + 1 / (4 * s.sinh(1)**2))
    base = {eta: 1, wave: 1, clock: 1}
    step_value = s.Rational(1, 10**6)
    numeric = s.N(
        (
            term.subs(base).subs(length, 1 + step_value)
            - term.subs(base).subs(length, 1 - step_value)
        ) / (2 * step_value),
        40,
    )
    assert abs(numeric - s.N(exact.subs(point), 40)) < s.Integer(10)**(-12)
    assert exact.subs(point) != 0
    signature = list(inspect.signature(rad.march).parameters)
    declared = joint.charge_identities()['Lambda']
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
        average_D_j='O(delta log q)',
        trace_derivative='-eta*csch(k*(T+D_j)/2)**2/4',
        trace_derivative_at_unit_point=str(s.factor(exact.subs(point))),
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def magnetic_variation():
    """A vanishing first variation does not measure cJ."""
    info = em.invariants()
    assert info['magnetic_background_action_and_first_variation'] == 'zero'
    require_sentence(EM_NOTE, 'do not measure cJ')
    require_sentence(EM_NOTE, 'does not establish the physical value cJ=0')
    coupling, kappa, current = s.symbols('cJ kappa J')
    local = coupling * kappa**4 * current**2 / 16
    zero(local.subs(current, 0))
    exact = s.diff(local, current)
    zero(exact.subs(current, 0))
    center = s.Integer(0)
    held = {coupling: s.Rational(1, 3), kappa: 1}

    def at(value):
        return local.subs(held).subs(current, value)

    zero(stencil(at, center) - exact.subs(held).subs(current, center))
    quadratic = s.sympify(
        info['quadratic_correction'].replace('^', '**')
    )
    unit = s.expand(quadratic.subs({s.Symbol('B0'): 1, s.Symbol('kappa'): 1}))
    assert s.Symbol('cJ') in unit.free_symbols
    assert not unit.is_number
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
        magnetic_first_variation='zero',
        quadratic_at_unit_background=str(s.factor(unit)),
        cJ_remains=True,
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    running = running_coupling()
    ratio = collar_ratio()
    length = exterior_length()
    magnetic = magnetic_variation()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S090',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s089_l00_derivative_recomputed=False,
        s090_substitution_difference_recomputed=False,
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
            running_coupling=running,
            collar_ratio=ratio,
            exterior_length=length,
            magnetic_variation=magnetic,
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
        running_coupling=running['classification'],
        beta_function=running['beta_function'],
        collar_ratio=ratio['classification'],
        ratio_at_predeclared_point=ratio['ratio_at_predeclared_point'],
        exterior_length=length['classification'],
        trace_derivative_at_unit_point=length['trace_derivative_at_unit_point'],
        magnetic_variation=magnetic['classification'],
        quadratic_at_unit_background=magnetic['quadratic_at_unit_background'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
