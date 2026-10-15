"""S105 checks for the four inputs named by S104.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-log-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-log-inputs.md'
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
        'Assumption. A number for `c_theta(mu0)` would be the logarithmic endpoint term `log(mu0)/(20*pi**2)`.',
        'That term is not a value of `c_theta(mu0)`.',
        'Assumption. The unequal-collar invariant would be the ratio of the first-bracket `w` derivative to the second-bracket `f\'` derivative.',
        'That ratio is not a two-collar theorem.',
        'Assumption. A number for one `D_j` would be the `a` derivative `-2/h**2`.',
        'That derivative is not an assignment of `D_j`.',
        'Assumption. A measured `cJ` would be the `B0` derivative `8*B0`.',
        'The expression `8*B0` does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S104 coefficient derivative is not recomputed as a new result',
        'The S104 weight derivative is not recomputed as a new result',
        'The S104 width derivative is not recomputed as a new result',
        'The S104 substituted susceptibility is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def logarithmic_endpoint():
    """The mu0 derivative of the log term does not assign c_theta."""
    require_sentence(CURVED_NOTE, 'b_theta=4/[5(4pi)^2]')
    scheme = s.Integer(4) / (5 * (4 * s.pi)**2)
    zero(scheme - 1 / (20 * s.pi**2))
    origin = s.Symbol('mu0')
    term = s.log(origin) / (20 * s.pi**2)
    exact = s.diff(term, origin)
    expected = 1 / (20 * s.pi**2 * origin)
    zero(s.expand(exact - expected))
    stored = str(s.factor(exact))
    assert stored == '1/(20*pi**2*mu0)'
    assert not exact.is_number
    assert origin in exact.free_symbols
    center = s.Rational(4, 9)
    assert center != 1
    step_symbol = s.Symbol('step')

    def at(value):
        return term.subs(origin, value)

    central = (at(center + step_symbol) - at(center - step_symbol)) / (2 * step_symbol)
    target = exact.subs(origin, center)
    zero(s.expand(s.limit(central, step_symbol, 0) - target))
    sampled = (at(center + NUMERIC_STEP) - at(center - NUMERIC_STEP)) / (2 * NUMERIC_STEP)
    assert abs(s.N(sampled - target, 40)) < s.Integer(10)**(-12)
    endpoint_coupling = s.Symbol('c_theta_mu0')
    endpoint = endpoint_coupling + term
    assert endpoint_coupling in endpoint.free_symbols
    assert not endpoint.is_number
    assert s.simplify(exact - endpoint_coupling) != 0
    assert s.simplify(term - endpoint_coupling) != 0
    declared = joint.charge_identities()
    assert declared['sigma'] == '2*c_theta*Lhat^{uu}'
    signature = list(inspect.signature(rad.march).parameters)
    assert absent('c_theta', ' '.join(signature))
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        c_theta_set_to_zero=False,
        c_theta_set_equal_to_b_theta=False,
        term_is_a_number=False,
        term_assigns_c_theta=False,
        mu0_derivative=stored,
        derivative_assigns_c_theta=False,
        coefficient_derivative_recomputed=False,
        endpoint_keeps_c_theta_mu0=True,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def bracket_ratio():
    """The ratio of the two bracket derivatives is not the collar mismatch."""
    require_sentence(CURVED_NOTE, 'b_theta (log f/2)(2/3)(R^2/f)w')
    require_sentence(CURVED_NOTE, '-d_theta (R^2/f)f\'/2')
    require_sentence(CURVED_NOTE, 'b_theta=4/[5(4pi)^2]')
    require_sentence(CURVED_NOTE, 'd_theta=4/[15(4pi)^2]')
    require_sentence(
        CURVED_NOTE,
        'not a simultaneous large-flux error estimate, a proof for unequal driven collars',
    )
    b_scheme = s.Integer(4) / (5 * (4 * s.pi)**2)
    d_scheme = s.Integer(4) / (15 * (4 * s.pi)**2)
    radius, lapse = s.symbols('R f')
    first = b_scheme * (s.log(lapse) / 2) * s.Rational(2, 3) * (radius**2 / lapse)
    second = -d_scheme * (radius**2 / lapse) / 2
    raw = first / second
    exact = s.simplify(raw)
    zero(s.expand(exact + 2 * s.log(lapse)))
    stored = str(s.factor(exact))
    assert stored == '-2*log(f)'
    assert stored != '1'
    assert not exact.is_number
    assert lapse in exact.free_symbols
    assert radius not in exact.free_symbols
    slope = s.diff(raw, radius)
    zero(s.expand(slope))
    center = s.Rational(3, 7)
    held_lapse = s.Rational(4, 9)
    assert held_lapse != 1

    def at(value):
        return raw.subs(lapse, held_lapse).subs(radius, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - slope))
    assert abs(s.N(sampled, 40)) < s.Integer(10)**(-18)
    numeric = s.N(raw.subs({radius: center, lapse: held_lapse}), 40)
    assert abs(numeric - s.N(-2 * s.log(held_lapse), 40)) < s.Integer(10)**(-18)
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        d_theta_set_to_zero=False,
        two_collar_equation_introduced=False,
        ratio_adopted_as_theorem=False,
        rn_series_recomputed=False,
        weight_derivative_recomputed=False,
        prime_derivative_recomputed=False,
        ratio=stored,
        ratio_equals_mismatch=False,
        recorded_mismatch='1',
        radius_derivative_of_ratio='0',
    )


def amplitude_derivative():
    """The a derivative of the width derivative does not assign one D_j."""
    require_sentence(
        MOUTH_NOTE,
        'eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.',
    )
    amplitude, width = s.symbols('a h')
    width_term = -2 * amplitude / width**2
    exact = s.diff(width_term, amplitude)
    expected = -2 / width**2
    zero(s.expand(exact - expected))
    stored = str(s.factor(exact))
    assert stored == '-2/h**2'
    assert stored != '-2*a/h**2'
    assert width in exact.free_symbols
    assert amplitude not in exact.free_symbols
    assert not exact.is_number
    center = s.Rational(3, 7)
    held = s.Rational(2, 5)
    assert held != 0

    def at(value):
        return width_term.subs(width, held).subs(amplitude, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    target = exact.subs(width, held)
    zero(s.expand(sampled - target))
    assert abs(s.N(sampled - target, 40)) < s.Integer(10)**(-18)
    clock_symbol = s.Symbol('D_j')
    clock_exact = s.diff(exact, clock_symbol)
    zero(clock_exact)

    def at_clock(value):
        return exact.subs(clock_symbol, value)

    clock_sampled = (
        at_clock(s.Integer(2) + STEP) - at_clock(s.Integer(2) - STEP)
    ) / (2 * STEP)
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
        width_derivative_recomputed=False,
        numeric_D_j=None,
        amplitude_derivative=stored,
        amplitude_derivative_in_D_j='0',
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def background_slope():
    """The B0 derivative of the substituted susceptibility does not measure cJ."""
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
    substituted = s.expand(reduced.subs(coefficient, 2))
    zero(s.expand(substituted - (1 + 4 * background**2)))
    exact = s.diff(substituted, background)
    zero(s.expand(exact - 8 * background))
    stored = str(s.factor(exact))
    assert stored == '8*B0'
    assert not exact.is_number
    assert background in exact.free_symbols
    assert coefficient not in exact.free_symbols
    center = s.Rational(2, 5)
    assert center != 1

    def at(value):
        return substituted.subs(background, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    target = exact.subs(background, center)
    zero(s.expand(sampled - target))
    assert abs(s.N(sampled - target, 40)) < s.Integer(10)**(-18)
    coupling_slope = s.diff(exact, coefficient)
    zero(s.expand(coupling_slope))
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
        substituted_susceptibility_recomputed=False,
        background_derivative=stored,
        derivative_is_cJ=False,
        derivative_measures_cJ=False,
        cJ_derivative_of_slope='0',
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    endpoint = logarithmic_endpoint()
    ratio = bracket_ratio()
    amplitude = amplitude_derivative()
    slope = background_slope()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S104',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s104_coefficient_derivative_recomputed=False,
        s104_weight_derivative_recomputed=False,
        s104_width_derivative_recomputed=False,
        s104_substituted_susceptibility_recomputed=False,
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
            logarithmic_endpoint=endpoint,
            bracket_ratio=ratio,
            amplitude_derivative=amplitude,
            background_slope=slope,
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
        mu0_derivative=endpoint['mu0_derivative'],
        ratio=ratio['ratio'],
        amplitude_derivative=amplitude['amplitude_derivative'],
        background_derivative=slope['background_derivative'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
