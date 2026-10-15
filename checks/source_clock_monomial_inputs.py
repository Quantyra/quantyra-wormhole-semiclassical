"""S104 checks for the four inputs named by S103.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-monomial-inputs.md. This checker does not march the
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

import finite_local_electromagnetic_input as em
import source_clock_joint as joint
import source_clock_radiation as rad


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'notes' / 'source-clock-monomial-inputs.md'
CURVED_NOTE = ROOT / 'notes' / 'curved-fermion-boundary.md'
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
NUMERIC_STEP = s.Rational(1, 10**8)


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def flat(path):
    return ' '.join(path.read_text(encoding='utf-8').split())


def require_sentence(path, sentence):
    assert sentence in flat(path), sentence


def require_pre_fit_note():
    for sentence in (
        'before any checker run',
        'Assumption. A number for `c_theta(mu0)` would be the scale `mu` in `c_one*mu`.',
        'That scale is not a value of `c_theta(mu0)`.',
        'Assumption. The unequal-collar invariant would be the first-bracket term `b_theta (log f/2)(2/3)(R^2/f)w`.',
        'That bracket term is not a two-collar theorem.',
        'Assumption. A number for one `D_j` would be the width derivative `-2*a/h**2`.',
        'That derivative is not an assignment of `D_j`.',
        'Assumption. A measured `cJ` would be the susceptibility `1+4*B0**2`.',
        'The expression `1+4*B0**2` does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S103 scale derivative is not recomputed as a new result',
        'The S103 bracket derivative is not recomputed as a new result',
        'The S103 ceiling is not recomputed as a new result',
        'The S103 ratio is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def monomial_scale():
    """The coefficient derivative returns the scale, which does not assign c_theta."""
    require_sentence(CURVED_NOTE, 'it renormalizes a dimension-one coefficient')
    coefficient, scale = s.symbols('c_one mu')
    monomial = coefficient * scale
    exact = s.diff(monomial, coefficient)
    zero(s.expand(exact - scale))
    assert str(s.factor(exact)) == 'mu'
    assert str(s.factor(exact)) != 'c_one'
    assert not exact.is_number
    center = s.Integer(2)
    held = s.Rational(3, 7)

    def at(value):
        return monomial.subs(scale, held).subs(coefficient, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact.subs(scale, held)))
    assert abs(s.N(sampled - held, 40)) < s.Integer(10)**(-18)
    endpoint_coupling, origin = s.symbols('c_theta_mu0 mu0')
    endpoint = endpoint_coupling + s.log(origin) / (20 * s.pi**2)
    assert endpoint_coupling in endpoint.free_symbols
    assert not endpoint.is_number
    assert s.simplify(exact - endpoint_coupling) != 0
    declared = joint.charge_identities()
    assert declared['sigma'] == '2*c_theta*Lhat^{uu}'
    signature = list(inspect.signature(rad.march).parameters)
    joined = ' '.join(signature)
    assert absent('c_theta', joined)
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        c_theta_set_to_zero=False,
        c_theta_set_equal_to_b_theta=False,
        scale_is_a_number=False,
        scale_assigns_c_theta=False,
        coefficient_derivative='mu',
        derivative_assigns_c_theta=False,
        scale_derivative_recomputed=False,
        endpoint_keeps_c_theta_mu0=True,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def first_bracket():
    """The w derivative of the first bracket is not the collar mismatch."""
    require_sentence(CURVED_NOTE, 'b_theta (log f/2)(2/3)(R^2/f)w')
    require_sentence(CURVED_NOTE, 'b_theta=4/[5(4pi)^2]')
    require_sentence(
        CURVED_NOTE,
        'not a simultaneous large-flux error estimate, a proof for unequal driven collars',
    )
    scheme = s.Integer(4) / (5 * (4 * s.pi)**2)
    zero(scheme - 1 / (20 * s.pi**2))
    radius, lapse, weight = s.symbols('R f w')
    term = scheme * (s.log(lapse) / 2) * s.Rational(2, 3) * (radius**2 / lapse) * weight
    exact = s.diff(term, weight)
    expected = radius**2 * s.log(lapse) / (60 * s.pi**2 * lapse)
    zero(s.expand(exact - expected))
    stored = str(s.factor(exact))
    assert stored == 'R**2*log(f)/(60*pi**2*f)'
    assert not exact.is_number
    assert stored != '1'
    assert stored != '-R**2/(120*pi**2*f)'
    center = s.Rational(2, 5)
    held = {radius: s.Rational(3, 7), lapse: s.Rational(4, 9)}

    def at(value):
        return term.subs(held).subs(weight, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    target = exact.subs(held)
    zero(s.expand(sampled - target))
    assert abs(s.N(sampled - target, 40)) < s.Integer(10)**(-18)
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        d_theta_set_to_zero=False,
        two_collar_equation_introduced=False,
        derivative_adopted_as_theorem=False,
        rn_series_recomputed=False,
        bracket_derivative_recomputed=False,
        weight_derivative=stored,
        derivative_equals_mismatch=False,
        recorded_mismatch='1',
    )


def width_derivative():
    """The h derivative of the ceiling does not assign one D_j."""
    require_sentence(
        MOUTH_NOTE,
        'eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.',
    )
    amplitude, width = s.symbols('a h')
    ceiling = 2 * amplitude / width
    exact = s.diff(ceiling, width)
    expected = -2 * amplitude / width**2
    zero(s.expand(exact - expected))
    stored = str(s.factor(exact))
    assert stored == '-2*a/h**2'
    assert stored != '2*a/h'
    assert amplitude in exact.free_symbols
    assert width in exact.free_symbols
    assert not exact.is_number
    center = s.Rational(2, 5)
    held = s.Rational(3, 7)
    step_symbol = s.Symbol('step')

    def at_width(value):
        return ceiling.subs(amplitude, held).subs(width, value)

    central = (
        at_width(center + step_symbol) - at_width(center - step_symbol)
    ) / (2 * step_symbol)
    target = exact.subs(amplitude, held).subs(width, center)
    zero(s.expand(s.limit(central, step_symbol, 0) - target))
    sampled = (
        at_width(center + NUMERIC_STEP) - at_width(center - NUMERIC_STEP)
    ) / (2 * NUMERIC_STEP)
    assert abs(s.N(sampled - target, 40)) < s.Integer(10)**(-12)
    clock_symbol = s.Symbol('D_j')
    clock_exact = s.diff(exact, clock_symbol)
    zero(clock_exact)

    def at_clock(value):
        return exact.subs(clock_symbol, value)

    clock_sampled = (at_clock(s.Integer(2) + STEP) - at_clock(s.Integer(2) - STEP)) / (2 * STEP)
    zero(s.expand(clock_sampled - clock_exact))
    assert abs(s.N(clock_sampled, 40)) < s.Integer(10)**(-18)
    declared = joint.charge_identities()['Lambda']
    rhs = s.sympify(declared.split('=', 1)[1], locals={'N': s.Symbol('N')})
    zero(s.diff(rhs, clock_symbol))
    signature = list(inspect.signature(rad.march).parameters)
    joined = ' '.join(signature)
    for name in ('D_j', 'a', 'epsilon', 'h'):
        assert absent(name, joined)
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        delta_set_to_zero=False,
        a_set_to_zero=False,
        h_set_to_zero=False,
        derivative_read_as_D_j=False,
        derivative_is_a_number=False,
        ceiling_reported_as_result=False,
        numeric_D_j=None,
        width_derivative=stored,
        width_derivative_in_D_j='0',
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def substituted_susceptibility():
    """Writing 2 in place of cJ does not measure cJ."""
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    info = em.invariants()
    susceptibility = s.sympify(
        info['parallel_electric_susceptibility'].replace('^', '**')
    )
    names = {symbol.name: symbol for symbol in susceptibility.free_symbols}
    gravitational = names['kappa']
    background = names['B0']
    coefficient = names['cJ']
    reduced = s.expand(susceptibility.subs(gravitational, 1))
    zero(s.expand(reduced - (1 + 2 * coefficient * background**2)))
    substituted = s.expand(reduced.subs(coefficient, 2))
    expected = 1 + 4 * background**2
    zero(s.expand(substituted - expected))
    stored = str(s.factor(substituted))
    assert stored == '4*B0**2 + 1'
    assert not substituted.is_number
    assert background in substituted.free_symbols
    assert coefficient not in substituted.free_symbols
    slope = s.diff(substituted, coefficient)
    zero(s.expand(slope))
    center = s.Rational(2, 5)
    assert center != 1

    def at(value):
        return substituted.subs(coefficient, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - slope))
    assert abs(s.N(sampled, 40)) < s.Integer(10)**(-18)
    numeric = s.N(substituted.subs(background, center), 40)
    assert abs(numeric - s.Rational(41, 25)) < s.Integer(10)**(-18)
    signature = list(inspect.signature(rad.march).parameters)
    joined = ' '.join(signature)
    for name in ('cJ', 'chi', 'B0', 'e_parallel'):
        assert absent(name, joined)
    assert 'kappa' in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        cJ_set_to_zero=False,
        b0_set_to_one=False,
        conditional_axial_window_adopted=False,
        positive_quadratic_mass_used_as_result=False,
        ratio_recomputed=False,
        substituted_susceptibility=stored,
        susceptibility_is_cJ=False,
        susceptibility_measures_cJ=False,
        cJ_derivative_of_substituted='0',
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    scale = monomial_scale()
    bracket = first_bracket()
    width = width_derivative()
    susceptibility = substituted_susceptibility()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S103',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s103_scale_derivative_recomputed=False,
        s103_bracket_derivative_recomputed=False,
        s103_ceiling_recomputed=False,
        s103_ratio_recomputed=False,
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
            monomial_scale=scale,
            first_bracket=bracket,
            width_derivative=width,
            substituted_susceptibility=susceptibility,
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
        coefficient_derivative=scale['coefficient_derivative'],
        weight_derivative=bracket['weight_derivative'],
        width_derivative=width['width_derivative'],
        substituted_susceptibility=susceptibility['substituted_susceptibility'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
