"""S102 checks for the four inputs named by S101.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-curvature-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-curvature-inputs.md'
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
CURVATURE = 'The a3 power subtraction is unchanged by curvature'


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def flat(path):
    return ' '.join(path.read_text(encoding='utf-8').split())


def require_sentence(path, sentence):
    assert sentence in flat(path), sentence


def require_pre_fit_note():
    for sentence in (
        'before any checker run',
        'Assumption. A number for `c_theta(mu0)` would be the statement that the a3 power subtraction is unchanged by curvature.',
        'That statement is not a value of `c_theta(mu0)`.',
        'Assumption. The unequal-collar invariant would be the leading term `-d_theta R` of the second bracket in (4).',
        'That leading term is not a two-collar theorem.',
        'Assumption. A number for one `D_j` would be the integral strength `eta=2a[1/(epsilon+h)-epsilon x^2/(1+b)]`.',
        'That strength is not an assignment of `D_j`.',
        'Assumption. A measured `cJ` would be the factor `2*B0**2` with `B0` kept symbolic.',
        'A factor `2*B0**2` does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S101 source derivative is not recomputed as a new result',
        'The S101 normal derivative is not recomputed as a new result',
        'The S101 decomposition is not recomputed as a new result',
        'The S101 reduced susceptibility is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def curvature_subtraction():
    """A vanishing curvature derivative does not assign c_theta(mu0)."""
    require_sentence(CURVED_NOTE, CURVATURE)
    assert any(character.isdigit() for character in CURVATURE)
    assert '=' not in CURVATURE
    assert CURVATURE.count('3') == 1
    subtraction, curvature = s.symbols('a3_power curvature')
    exact = s.diff(subtraction, curvature)
    zero(exact)
    center = s.Integer(2)

    def at(value):
        return subtraction.subs(curvature, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    assert abs(s.N(sampled - exact, 40)) < s.Integer(10)**(-18)
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
        sentence_assigns_c_theta=False,
        label_digit_is_c_theta=False,
        curvature_derivative='0',
        derivative_assigns_c_theta=False,
        source_derivative_recomputed=False,
        endpoint_keeps_c_theta_mu0=True,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def leading_term():
    """The radial derivative of -d_theta R is not the collar mismatch."""
    require_sentence(CURVED_NOTE, 'the second bracket has leading term -d_theta R')
    require_sentence(CURVED_NOTE, 'd_theta=4/[15(4pi)^2]')
    require_sentence(
        CURVED_NOTE,
        'not a simultaneous large-flux error estimate, a proof for unequal driven collars',
    )
    require_sentence(
        COUPLED_NOTE,
        'Unequal/background-dependent collars require their additional invariants and are not covered by (8).',
    )
    scheme = s.Integer(4) / (15 * (4 * s.pi)**2)
    zero(scheme - 1 / (60 * s.pi**2))
    radius = s.Symbol('R')
    term = -scheme * radius
    exact = s.diff(term, radius)
    zero(s.expand(exact + 1 / (60 * s.pi**2)))
    stored = str(s.factor(exact))
    assert stored == '-1/(60*pi**2)'
    assert exact.is_number
    assert s.simplify(exact - 1) != 0
    assert stored != 'u_h**2/(60*pi**2)'
    assert not term.is_number
    center = s.Integer(2)

    def at(value):
        return term.subs(radius, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    assert abs(s.N(sampled - exact, 40)) < s.Integer(10)**(-18)
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        d_theta_set_to_zero=False,
        two_collar_equation_introduced=False,
        derivative_adopted_as_theorem=False,
        rn_series_recomputed=False,
        normal_derivative_recomputed=False,
        radial_derivative=stored,
        derivative_equals_mismatch=False,
        recorded_mismatch='1',
    )


def integral_strength():
    """The integral strength eta does not assign one D_j."""
    require_sentence(
        MOUTH_NOTE,
        'eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.',
    )
    amplitude, offset, width, slope, shift = s.symbols('a epsilon h x b')
    strength = 2 * amplitude * (
        1 / (offset + width) - offset * slope**2 / (1 + shift)
    )
    for symbol in (amplitude, offset, width, slope, shift):
        assert symbol in strength.free_symbols
    assert not strength.is_number
    stored = str(s.simplify(strength))
    ceiling = 2 * amplitude / width
    assert s.simplify(strength - ceiling) != 0
    clock_symbol = s.Symbol('D_j')
    exact = s.diff(strength, clock_symbol)
    zero(exact)
    center = s.Integer(2)

    def at(value):
        return strength.subs(clock_symbol, value)

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
        epsilon_set_to_zero=False,
        h_set_to_zero=False,
        strength_read_as_D_j=False,
        strength_is_a_number=False,
        strength_equals_ceiling=False,
        decomposition_recomputed=False,
        numeric_D_j=None,
        strength=stored,
        strength_derivative_in_D_j='0',
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def factor_derivative():
    """The cJ derivative of 2*B0**2 does not measure cJ."""
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    background = s.Symbol('B0')
    coefficient = s.Symbol('cJ')
    factor = 2 * background**2
    assert str(s.factor(factor)) == '2*B0**2'
    assert not factor.is_number
    assert background in factor.free_symbols
    assert coefficient not in factor.free_symbols
    exact = s.diff(factor, coefficient)
    zero(exact)
    center = s.Rational(3, 7)
    held_background = s.Rational(2, 5)

    def at(value):
        return factor.subs(background, held_background).subs(coefficient, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    assert abs(s.N(sampled, 40)) < s.Integer(10)**(-18)
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
        reduced_susceptibility_recomputed=False,
        factor='2*B0**2',
        factor_is_a_number=False,
        cJ_derivative='0',
        derivative_measures_cJ=False,
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    curvature = curvature_subtraction()
    radial = leading_term()
    strength = integral_strength()
    factor = factor_derivative()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S101',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s101_source_derivative_recomputed=False,
        s101_normal_derivative_recomputed=False,
        s101_decomposition_recomputed=False,
        s101_reduced_susceptibility_recomputed=False,
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
            curvature_subtraction=curvature,
            leading_term=radial,
            integral_strength=strength,
            factor_derivative=factor,
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
        curvature_derivative=curvature['curvature_derivative'],
        radial_derivative=radial['radial_derivative'],
        strength=strength['strength'],
        strength_derivative=strength['strength_derivative_in_D_j'],
        factor=factor['factor'],
        cJ_derivative=factor['cJ_derivative'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
