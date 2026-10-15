"""S101 checks for the four inputs named by S100.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-heat-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-heat-inputs.md'
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
HEAT = (
    'No one-loop heat computation determines the apparatus implementing the source.'
)
UMBILIC = (
    'An umbilic boundary has no unsmeared derivative logarithm, but generally '
    'still has the second term in (3).'
)


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def flat(path):
    return ' '.join(path.read_text(encoding='utf-8').split())


def require_sentence(path, sentence):
    assert sentence in flat(path), sentence


def require_pre_fit_note():
    for sentence in (
        'before any checker run',
        'Assumption. A number for `c_theta(mu0)` would be the statement that no one-loop heat computation determines the apparatus implementing the source.',
        'That statement is not a value of `c_theta(mu0)`.',
        'Assumption. The unequal-collar invariant would be the statement that an umbilic boundary has no unsmeared derivative logarithm but generally still has the second term in (3).',
        'That second term is not a two-collar theorem.',
        'Assumption. A number for one `D_j` would be the decomposition `W=V_x-P`.',
        'That decomposition is not an assignment of `D_j`.',
        'Assumption. A measured `cJ` would be the parallel susceptibility `1+2*cJ*B0**2` at gravitational `kappa = 1` with `B0` kept symbolic.',
        'The expression `1+2*cJ*B0**2` does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S100 explicitness requirement is not recomputed as a new result',
        'The S100 omission gap is not recomputed as a new result',
        'The S100 interface profile is not recomputed as a new result',
        'The S100 coupling derivative is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def heat_sentence():
    """The heat-computation denial does not assign c_theta(mu0)."""
    require_sentence(CURVED_NOTE, HEAT)
    require_sentence(
        CURVED_NOTE,
        'dimension-one coefficient of integral sqrt(gamma)u^2',
    )
    assert not any(character.isdigit() for character in HEAT)
    assert '=' not in HEAT
    coefficient, density, square = s.symbols('c_one sqrt_gamma u2')
    source = coefficient * density * square
    exact = s.diff(source, coefficient)
    zero(s.expand(exact - density * square))
    assert str(s.factor(exact)) == 'sqrt_gamma*u2'
    assert not exact.is_number
    assert coefficient not in exact.free_symbols
    center = s.Integer(2)
    held = {density: s.Rational(3, 7), square: s.Rational(2, 5)}

    def at(value):
        return source.subs(held).subs(coefficient, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    expected = exact.subs(held)
    zero(s.expand(sampled - expected))
    assert abs(s.N(sampled - expected, 40)) < s.Integer(10)**(-18)
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
        sentence_is_a_number=False,
        sentence_assigns_c_theta=False,
        source_derivative='sqrt_gamma*u2',
        derivative_assigns_c_theta=False,
        assignment_derivative_recomputed=False,
        endpoint_keeps_c_theta_mu0=True,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def second_term():
    """The normal derivative of the second term in (3) is not the mismatch."""
    require_sentence(CURVED_NOTE, UMBILIC)
    require_sentence(CURVED_NOTE, 'd_theta(n_h sigma)u_h^2')
    require_sentence(CURVED_NOTE, 'd_theta=4/[15(4pi)^2]')
    require_sentence(
        COUPLED_NOTE,
        'Unequal/background-dependent collars require their additional invariants and are not covered by (8).',
    )
    scheme = s.Integer(4) / (15 * (4 * s.pi)**2)
    zero(scheme - 1 / (60 * s.pi**2))
    normal, amplitude = s.symbols('n_sigma u_h')
    term = scheme * normal * amplitude**2
    exact = s.diff(term, normal)
    zero(s.expand(exact - amplitude**2 / (60 * s.pi**2)))
    stored = str(s.factor(exact))
    assert stored == 'u_h**2/(60*pi**2)'
    assert not exact.is_number
    assert s.simplify(exact - 1) != 0
    assert stored != '4*L_ab*uu/5'
    assert stored != '1/(60*pi**2)'
    center = s.Rational(2, 5)
    held_amplitude = s.Rational(3, 7)

    def at(value):
        return term.subs(amplitude, held_amplitude).subs(normal, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    expected = exact.subs(amplitude, held_amplitude)
    zero(s.expand(sampled - expected))
    assert abs(s.N(sampled - expected, 40)) < s.Integer(10)**(-18)
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        d_theta_set_to_zero=False,
        two_collar_equation_introduced=False,
        derivative_adopted_as_theorem=False,
        shipped_trace_identified_with_derivative=False,
        omission_gap_recomputed=False,
        scheme_value_reported_as_result=False,
        normal_derivative=stored,
        derivative_equals_mismatch=False,
        recorded_mismatch='1',
    )


def potential_decomposition():
    """W = V_x - P does not assign one D_j."""
    require_sentence(MOUTH_NOTE, 'W=V_x-P.')
    require_sentence(MOUTH_NOTE, 'V_x(d)=a/p_x(d)^2.')
    require_sentence(MOUTH_NOTE, 'W=a(r-p_x)(r+p_x)/(p_x^2 r^2)')
    amplitude, radius, profile = s.symbols('a r p_x')
    difference = amplitude / profile**2 - amplitude / radius**2
    displayed = (
        amplitude * (radius - profile) * (radius + profile)
        / (profile**2 * radius**2)
    )
    zero(s.simplify(difference - displayed))
    assert amplitude in difference.free_symbols
    assert radius in difference.free_symbols
    assert profile in difference.free_symbols
    assert not difference.is_number
    stored = str(s.factor(s.together(difference)))
    assert stored != '64*a/(8*epsilon + 1)**2'
    assert stored != '0'
    clock_symbol = s.Symbol('D_j')
    exact = s.diff(difference, clock_symbol)
    zero(exact)
    center = s.Integer(2)

    def at(value):
        return difference.subs(clock_symbol, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    assert abs(s.N(sampled - exact, 40)) < s.Integer(10)**(-18)
    point = {
        amplitude: s.Rational(3, 7),
        radius: s.Rational(2, 5),
        profile: s.Rational(1, 4),
    }
    direct_gap = s.simplify(difference.subs(point) - displayed.subs(point))
    zero(direct_gap)
    assert abs(s.N(difference.subs(point) - displayed.subs(point), 40)) < s.Integer(10)**(-18)
    assert s.simplify(difference.subs(point)) != 0
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
        decomposition_read_as_D_j=False,
        decomposition_is_a_number=False,
        interface_profile_recomputed=False,
        numeric_D_j=None,
        decomposition=stored,
        decomposition_matches_display=True,
        decomposition_derivative_in_D_j='0',
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def reduced_susceptibility():
    """The reduced susceptibility at gravitational kappa = 1 does not measure cJ."""
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
    canonical = 1 + 2 * coefficient * background**2
    zero(s.expand(reduced - canonical))
    assert str(s.factor(canonical)) == '2*B0**2*cJ + 1'
    exact = s.diff(reduced, coefficient)
    zero(s.expand(exact - 2 * background**2))
    assert str(s.factor(exact)) == '2*B0**2'
    assert not reduced.is_number
    assert not exact.is_number
    assert background in exact.free_symbols
    assert coefficient not in exact.free_symbols
    assert str(s.factor(exact)) != '2'
    assert str(s.factor(exact)) != '4*B0'
    assert str(s.factor(exact)) != '4*B0*cJ'
    assert str(s.factor(canonical)) != '4*B0'
    center = s.Rational(3, 7)
    held_background = s.Rational(2, 5)

    def at(value):
        return reduced.subs(background, held_background).subs(coefficient, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    expected = exact.subs(background, held_background)
    zero(s.expand(sampled - expected))
    assert abs(s.N(sampled - expected, 40)) < s.Integer(10)**(-18)
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
        coupling_derivative_recomputed=False,
        susceptibility_solve_recomputed=False,
        reduced_susceptibility='1+2*cJ*B0**2',
        reduced_factored='2*B0**2*cJ + 1',
        reduced_is_a_number=False,
        cJ_derivative='2*B0**2',
        derivative_is_a_number=False,
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
    heat = heat_sentence()
    term = second_term()
    decomposition = potential_decomposition()
    susceptibility = reduced_susceptibility()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S100',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s100_explicitness_recomputed=False,
        s100_omission_gap_recomputed=False,
        s100_interface_profile_recomputed=False,
        s100_coupling_derivative_recomputed=False,
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
            heat_sentence=heat,
            second_term=term,
            potential_decomposition=decomposition,
            reduced_susceptibility=susceptibility,
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
        heat_sentence=heat['classification'],
        source_derivative=heat['source_derivative'],
        second_term=term['classification'],
        normal_derivative=term['normal_derivative'],
        potential_decomposition=decomposition['classification'],
        decomposition=decomposition['decomposition'],
        reduced_susceptibility=susceptibility['classification'],
        reduced=susceptibility['reduced_susceptibility'],
        cJ_derivative=susceptibility['cJ_derivative'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
