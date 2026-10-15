"""S100 checks for the four inputs named by S099.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-explicit-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-explicit-inputs.md'
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
EXPLICIT = (
    'the finite J_theta choice and the finite c_theta must remain explicit.'
)
OMISSION = (
    'Omitting that cross term would leave a spurious nonzero unsmeared '
    'contribution on an umbilic surface.'
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
        'Assumption. A number for `c_theta(mu0)` would be the requirement that the finite `c_theta` must remain explicit.',
        'That requirement is not a value of `c_theta(mu0)`.',
        'Assumption. The unequal-collar invariant would be the statement that omitting the cross term would leave a spurious nonzero unsmeared contribution on an umbilic surface.',
        'That omission statement is not a two-collar theorem.',
        'Assumption. A number for one `D_j` would be the cell profile `P(h)=a/(epsilon+h)^2` at `h=1/8`.',
        'A cell profile is not an assignment of `D_j`.',
        'Assumption. A measured `cJ` would be the `cJ` derivative of `4*B0*cJ` with `B0` kept symbolic.',
        'A factor `4*B0` does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S099 free-sector sentence is not recomputed as a new result',
        'The S099 curvature coefficient is not recomputed as a new result',
        'The S099 mouth width is not recomputed as a new result',
        'The S099 susceptibility slope is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def explicit_coupling():
    """The explicitness requirement does not assign c_theta(mu0)."""
    require_sentence(CURVED_NOTE, EXPLICIT)
    assert not any(character.isdigit() for character in EXPLICIT)
    assert '=' not in EXPLICIT
    coupling = s.Symbol('c_theta')
    assigned = s.Symbol('assigned')
    exact = s.diff(coupling, assigned)
    zero(exact)
    center = s.Integer(2)

    def at(value):
        return coupling.subs(assigned, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    assert abs(s.N(sampled - exact, 40)) < s.Integer(10)**(-18)
    endpoint_coupling, origin = s.symbols('c_theta_mu0 mu0')
    endpoint = endpoint_coupling + s.log(origin) / (20 * s.pi**2)
    assert endpoint_coupling in endpoint.free_symbols
    assert not endpoint.is_number
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
        assignment_derivative='0',
        endpoint_keeps_c_theta_mu0=True,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def omission_gap():
    """The umbilic omission statement is not the collar mismatch."""
    require_sentence(CURVED_NOTE, OMISSION)
    require_sentence(
        CURVED_NOTE,
        '(4/5)(L_ab-L gamma_ab/3)u^a u^b',
    )
    require_sentence(
        COUPLED_NOTE,
        'Unequal/background-dependent collars require their additional invariants and are not covered by (8).',
    )
    traceless, length, metric, bilinear = s.symbols('L_ab L gamma uu')
    full = s.Rational(4, 5) * (traceless - length * metric / 3) * bilinear
    omitted = s.Rational(4, 5) * (-length * metric / 3) * bilinear
    gap = s.simplify(full - omitted)
    exact = s.diff(gap, traceless)
    zero(s.expand(exact - s.Rational(4, 5) * bilinear))
    assert traceless in gap.free_symbols
    assert not gap.is_number
    assert s.simplify(gap - 1) != 0
    assert str(s.factor(gap)) != '-4'
    assert str(s.factor(gap)) != '-4*uu'
    center = s.Integer(2)

    def at(value):
        return gap.subs(traceless, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    numeric = sampled.subs(bilinear, s.Rational(3, 7))
    expected = exact.subs(bilinear, s.Rational(3, 7))
    assert abs(s.N(numeric - expected, 40)) < s.Integer(10)**(-18)
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        gap_adopted_as_theorem=False,
        shipped_trace_identified_with_gap=False,
        omission_gap=str(s.factor(gap)),
        gap_equals_mismatch=False,
        recorded_mismatch='1',
    )


def interface_profile():
    """The outer-interface profile does not assign one D_j."""
    require_sentence(MOUTH_NOTE, 'P(d)=a/(epsilon+d)^2')
    require_sentence(MOUTH_NOTE, 'Take h=1/8')
    amplitude, width_offset = s.symbols('a epsilon')
    clock_symbol = s.Symbol('D_j')
    profile = amplitude / (width_offset + s.Rational(1, 8))**2
    assert amplitude in profile.free_symbols
    assert width_offset in profile.free_symbols
    assert not profile.is_number
    exact = s.diff(profile, clock_symbol)
    zero(exact)
    center = s.Integer(2)

    def at(value):
        return profile.subs(clock_symbol, value)

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
        profile_read_as_D_j=False,
        profile_is_a_number=False,
        numeric_D_j=None,
        profile=str(s.simplify(profile)),
        profile_derivative_in_D_j='0',
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def coupling_derivative():
    """The cJ derivative of the symbolic slope does not measure cJ."""
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    info = em.invariants()
    susceptibility = s.sympify(
        info['parallel_electric_susceptibility'].replace('^', '**')
    )
    names = {symbol.name: symbol for symbol in susceptibility.free_symbols}
    gravitational = names['kappa']
    background = names['B0']
    coefficient = names['cJ']
    reduced = susceptibility.subs(gravitational, 1)
    slope = s.diff(reduced, background)
    zero(s.expand(slope - 4 * coefficient * background))
    exact = s.diff(slope, coefficient)
    zero(s.expand(exact - 4 * background))
    assert not exact.is_number
    assert background in exact.free_symbols
    assert coefficient not in exact.free_symbols
    assert str(s.factor(exact)) != '4'
    assert str(s.factor(exact)) != '4*B0*cJ'
    center = s.Rational(3, 7)
    held_background = s.Rational(2, 5)

    def at(value):
        return slope.subs(background, held_background).subs(coefficient, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    expected = exact.subs(background, held_background)
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
        cJ_derivative=str(s.factor(exact)),
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
    explicit = explicit_coupling()
    gap = omission_gap()
    profile = interface_profile()
    slope = coupling_derivative()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S099',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s099_free_sector_recomputed=False,
        s099_curvature_coefficient_recomputed=False,
        s099_mouth_width_recomputed=False,
        s099_susceptibility_slope_recomputed=False,
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
            explicit_coupling=explicit,
            omission_gap=gap,
            interface_profile=profile,
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
        explicit_coupling=explicit['classification'],
        sentence_is_a_number=explicit['sentence_is_a_number'],
        omission_gap=gap['classification'],
        gap=gap['omission_gap'],
        interface_profile=profile['classification'],
        profile=profile['profile'],
        coupling_derivative=slope['classification'],
        cJ_derivative=slope['cJ_derivative'],
        derivative_is_a_number=slope['derivative_is_a_number'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
