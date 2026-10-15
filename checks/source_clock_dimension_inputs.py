"""S103 checks for the four inputs named by S102.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-dimension-inputs.md. This checker does not march the
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

import source_clock_joint as joint
import source_clock_radiation as rad


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'notes' / 'source-clock-dimension-inputs.md'
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
CLAUSE = 'it renormalizes a dimension-one coefficient'


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def flat(path):
    return ' '.join(path.read_text(encoding='utf-8').split())


def require_sentence(path, sentence):
    assert sentence in flat(path), sentence


def require_pre_fit_note():
    for sentence in (
        'before any checker run',
        'Assumption. A number for `c_theta(mu0)` would be the statement that the subtraction renormalizes a dimension-one coefficient.',
        'That statement is not a value of `c_theta(mu0)`.',
        'Assumption. The unequal-collar invariant would be the displayed second-bracket term `-d_theta (R^2/f)f\'/2`.',
        'That bracket term is not a two-collar theorem.',
        'Assumption. A number for one `D_j` would be the ceiling `2*a/h`.',
        'That ceiling is not an assignment of `D_j`.',
        'Assumption. A measured `cJ` would be the ratio `(2*B0**2)/B0**2`.',
        'The integer `2` does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S102 curvature derivative is not recomputed as a new result',
        'The S102 radial derivative is not recomputed as a new result',
        'The S102 integral strength is not recomputed as a new result',
        'The S102 factor derivative is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def dimension_coefficient():
    """The scale derivative of the dimension-one monomial does not assign c_theta."""
    require_sentence(CURVED_NOTE, CLAUSE)
    assert not any(character.isdigit() for character in CLAUSE)
    assert '=' not in CLAUSE
    coefficient, scale = s.symbols('c_one mu')
    monomial = coefficient * scale
    exact = s.diff(monomial, scale)
    zero(s.expand(exact - coefficient))
    assert str(s.factor(exact)) == 'c_one'
    assert not exact.is_number
    center = s.Integer(2)
    held = s.Rational(3, 7)

    def at(value):
        return monomial.subs(coefficient, held).subs(scale, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact.subs(coefficient, held)))
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
        clause_is_a_number=False,
        clause_assigns_c_theta=False,
        scale_derivative='c_one',
        derivative_assigns_c_theta=False,
        curvature_derivative_recomputed=False,
        endpoint_keeps_c_theta_mu0=True,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def second_bracket():
    """The f' derivative of the second bracket is not the collar mismatch."""
    require_sentence(CURVED_NOTE, '-d_theta (R^2/f)f\'/2')
    require_sentence(CURVED_NOTE, 'd_theta=4/[15(4pi)^2]')
    require_sentence(
        CURVED_NOTE,
        'not a simultaneous large-flux error estimate, a proof for unequal driven collars',
    )
    scheme = s.Integer(4) / (15 * (4 * s.pi)**2)
    zero(scheme - 1 / (60 * s.pi**2))
    radius, lapse, prime = s.symbols('R f f_prime')
    term = -scheme * (radius**2 / lapse) * prime / 2
    exact = s.diff(term, prime)
    expected = -radius**2 / (120 * s.pi**2 * lapse)
    zero(s.expand(exact - expected))
    stored = str(s.factor(exact))
    assert stored == '-R**2/(120*pi**2*f)'
    assert not exact.is_number
    assert radius in exact.free_symbols
    assert lapse in exact.free_symbols
    assert prime not in exact.free_symbols
    assert s.simplify(exact - 1) != 0
    assert stored != '-1/(60*pi**2)'
    center = s.Rational(2, 5)
    held = {radius: s.Rational(3, 7), lapse: s.Rational(4, 9)}

    def at(value):
        return term.subs(held).subs(prime, value)

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
        radial_derivative_recomputed=False,
        prime_derivative=stored,
        derivative_equals_mismatch=False,
        recorded_mismatch='1',
    )


def strength_ceiling():
    """The ceiling 2*a/h does not assign one D_j."""
    require_sentence(
        MOUTH_NOTE,
        'eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.',
    )
    amplitude, width = s.symbols('a h')
    ceiling = 2 * amplitude / width
    assert str(s.factor(ceiling)) == '2*a/h'
    assert amplitude in ceiling.free_symbols
    assert width in ceiling.free_symbols
    assert not ceiling.is_number
    clock_symbol = s.Symbol('D_j')
    exact = s.diff(ceiling, clock_symbol)
    zero(exact)
    center = s.Integer(2)

    def at(value):
        return ceiling.subs(clock_symbol, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    assert abs(s.N(sampled - exact, 40)) < s.Integer(10)**(-18)
    declared = joint.charge_identities()['Lambda']
    rhs = s.sympify(declared.split('=', 1)[1], locals={'N': s.Symbol('N')})
    zero(s.diff(rhs, clock_symbol))
    zero(s.simplify(rhs.subs(clock_symbol, 1) - rhs.subs(clock_symbol, 2)))
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
        ceiling_read_as_D_j=False,
        ceiling_is_a_number=False,
        strength_reported_as_result=False,
        numeric_D_j=None,
        ceiling='2*a/h',
        ceiling_derivative_in_D_j='0',
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def ratio_value():
    """Cancelling B0**2 in 2*B0**2 returns 2, which does not measure cJ."""
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    background = s.Symbol('B0')
    coefficient = s.Symbol('cJ')
    raw = (2 * background**2) / background**2
    exact_value = s.simplify(raw)
    zero(exact_value - 2)
    assert exact_value == 2
    assert exact_value.is_number
    assert coefficient not in raw.free_symbols
    slope = s.diff(raw, background)
    zero(s.expand(slope))
    center = s.Rational(2, 5)
    assert center != 1

    def at(value):
        return raw.subs(background, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - slope))
    assert abs(s.N(sampled, 40)) < s.Integer(10)**(-18)
    numeric = s.N(raw.subs(background, center), 40)
    assert abs(numeric - 2) < s.Integer(10)**(-18)
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
        factor_derivative_recomputed=False,
        ratio='2',
        ratio_is_cJ=False,
        ratio_measures_cJ=False,
        b0_derivative_of_ratio='0',
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    dimension = dimension_coefficient()
    bracket = second_bracket()
    ceiling = strength_ceiling()
    ratio = ratio_value()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S102',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s102_curvature_derivative_recomputed=False,
        s102_radial_derivative_recomputed=False,
        s102_integral_strength_recomputed=False,
        s102_factor_derivative_recomputed=False,
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
            dimension_coefficient=dimension,
            second_bracket=bracket,
            strength_ceiling=ceiling,
            ratio_value=ratio,
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
        scale_derivative=dimension['scale_derivative'],
        prime_derivative=bracket['prime_derivative'],
        ceiling=ceiling['ceiling'],
        ceiling_derivative=ceiling['ceiling_derivative_in_D_j'],
        ratio=ratio['ratio'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
