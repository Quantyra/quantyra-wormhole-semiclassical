"""S098 checks for the four inputs named by S097.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-relation-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-relation-inputs.md'
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
        'Assumption. A number for `c_theta(mu0)` would be the displayed scale relation `mu partial_mu Gamma_E=b_theta I_theta`.',
        'That scale relation is not a value of `c_theta(mu0)`.',
        'Assumption. The unequal-collar invariant would be the derivative of the curvature combination in `L_ab`.',
        'A derivative in `L_ab` is not a two-collar theorem.',
        'Assumption. A number for one `D_j` would be the statement that no physical four-dimensional matching coefficient follows from holding `c` fixed.',
        'Holding `c` fixed is not an assignment.',
        'Assumption. A measured `cJ` would be the `B0` derivative of `B0**2*cJ`.',
        'A factor `2*B0*cJ` does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S097 functional derivative is not recomputed as a new result',
        'The S097 curvature derivative is not recomputed as a new result',
        'The S097 self-energy derivative is not recomputed as a new result',
        'The S097 quadratic quotient is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def scale_relation():
    """The displayed scale relation does not assign c_theta(mu0)."""
    require_sentence(CURVED_NOTE, 'mu partial_mu Gamma_E=b_theta I_theta')
    scale, beta, functional = s.symbols('mu b_theta I_theta')
    rhs = beta * functional
    exact = s.diff(rhs, scale)
    zero(exact)
    assert not rhs.is_number
    center = s.Integer(2)

    def at(value):
        return rhs.subs(scale, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    probe = s.Rational(3, 7)
    numeric = (at(probe + STEP) - at(probe - STEP)) / (2 * STEP)
    assert abs(s.N(numeric - exact, 40)) < s.Integer(10)**(-18)
    coupling, origin = s.symbols('c_theta_mu0 mu0')
    endpoint = coupling + s.log(origin) / (20 * s.pi**2)
    assert coupling in endpoint.free_symbols
    assert not endpoint.is_number
    assert s.simplify(endpoint - exact) != 0
    declared = joint.charge_identities()
    assert declared['sigma'] == '2*c_theta*Lhat^{uu}'
    signature = list(inspect.signature(rad.march).parameters)
    joined = ' '.join(signature)
    for name in ('c_theta', 'b_theta', 'I_theta', 'mu'):
        assert absent(name, joined)
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        c_theta_set_to_zero=False,
        c_theta_set_equal_to_b_theta=False,
        mu_derivative=str(exact),
        rhs_is_a_number=False,
        derivative_assigns_c_theta=False,
        endpoint_keeps_c_theta_mu0=True,
        endpoint_equals_mu_derivative=False,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def curvature_slot():
    """The L_ab derivative is not the unequal-collar mismatch."""
    require_sentence(
        CURVED_NOTE,
        'Delta tr[(D chi)^2 S]=-2L u^2-4L_ab u^a u^b',
    )
    require_sentence(
        COUPLED_NOTE,
        'Unequal/background-dependent collars require their additional invariants and are not covered by (8).',
    )
    length, traceless, square, bilinear = s.symbols('L L_ab u2 uu')
    combination = -2 * length * square - 4 * traceless * bilinear
    exact = s.diff(combination, traceless)
    zero(s.simplify(exact + 4 * bilinear))
    assert bilinear in exact.free_symbols
    assert not exact.is_number
    assert s.simplify(exact - 1) != 0
    center = s.Integer(2)

    def at(value):
        return combination.subs(traceless, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    values = {bilinear: s.Rational(3, 7), square: s.Rational(2, 5), length: s.Rational(1, 4)}
    numeric = sampled.subs(values)
    expected = exact.subs(values)
    assert abs(s.N(numeric - expected, 40)) < s.Integer(10)**(-18)
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        derivative_adopted_as_theorem=False,
        shipped_trace_identified_with_derivative=False,
        derivative_in_L_ab=str(s.factor(exact)),
        derivative_equals_mismatch=False,
        recorded_mismatch='1',
    )


def fixed_coefficient():
    """Holding c fixed does not assign one D_j."""
    require_sentence(
        MOUTH_NOTE,
        'No claim about a physical four-dimensional matching coefficient follows from holding c fixed.',
    )
    held = s.Symbol('c')
    clock_symbol = s.Symbol('D_j')
    exact = s.diff(held, clock_symbol)
    zero(exact)
    assert not held.is_number
    center = s.Integer(2)

    def at(value):
        return held.subs(clock_symbol, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    assert abs(s.N(sampled - exact, 40)) < s.Integer(10)**(-18)
    declared = joint.charge_identities()['Lambda']
    rhs = s.sympify(declared.split('=', 1)[1], locals={'N': s.Symbol('N')})
    zero(s.diff(rhs, clock_symbol))
    zero(s.simplify(rhs.subs(clock_symbol, 1) - rhs.subs(clock_symbol, 2)))
    signature = list(inspect.signature(rad.march).parameters)
    joined = ' '.join(signature)
    for name in ('D_j', 'c'):
        assert absent(name, joined)
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        delta_set_to_zero=False,
        c_set_to_zero=False,
        channel_count_assigned=False,
        denial_read_as_equality=False,
        fixed_c_is_a_number=False,
        numeric_D_j=None,
        fixed_c_derivative_in_D_j='0',
        derivative_assigns_D_j=False,
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def quotient_slope():
    """The B0 derivative of the quadratic quotient does not measure cJ."""
    info = em.invariants()
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    correction = s.sympify(info['quadratic_correction'].replace('^', '**'))
    field = s.Symbol('e_parallel')
    grav = s.Symbol('kappa')
    coupling = s.Symbol('cJ')
    background = s.Symbol('B0')
    quotient = s.simplify((correction / field**2).subs(grav, 1))
    exact = s.diff(quotient, background)
    zero(s.simplify(exact - 2 * coupling * background))
    assert coupling in exact.free_symbols
    assert background in exact.free_symbols
    assert not exact.is_number
    probe_coupling = s.Rational(3, 7)
    probe_background = s.Rational(2, 5)
    probe_field = s.Rational(4, 3)
    prepared = (correction / field**2).subs({
        grav: 1,
        coupling: probe_coupling,
        field: probe_field,
    })

    def at(value):
        return prepared.subs(background, value)

    numeric = (at(probe_background + STEP) - at(probe_background - STEP)) / (2 * STEP)
    expected = exact.subs({coupling: probe_coupling, background: probe_background})
    assert abs(s.N(numeric - expected, 40)) < s.Integer(10)**(-18)
    signature = list(inspect.signature(rad.march).parameters)
    for name in ('cJ', 'chi', 'B0', 'e_parallel'):
        assert name not in signature
    assert 'kappa' in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        cJ_set_to_zero=False,
        conditional_axial_window_adopted=False,
        positive_quadratic_mass_used_as_result=False,
        b0_derivative=str(s.factor(exact)),
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
    relation = scale_relation()
    slot = curvature_slot()
    held = fixed_coefficient()
    slope = quotient_slope()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S097',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s097_functional_derivative_recomputed=False,
        s097_curvature_derivative_recomputed=False,
        s097_self_energy_derivative_recomputed=False,
        s097_quadratic_quotient_recomputed=False,
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
            scale_relation=relation,
            curvature_slot=slot,
            fixed_coefficient=held,
            quotient_slope=slope,
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
        scale_relation=relation['classification'],
        mu_derivative=relation['mu_derivative'],
        curvature_slot=slot['classification'],
        derivative_in_L_ab=slot['derivative_in_L_ab'],
        fixed_coefficient=held['classification'],
        fixed_c_derivative_in_D_j=held['fixed_c_derivative_in_D_j'],
        quotient_slope=slope['classification'],
        b0_derivative=slope['b0_derivative'],
        derivative_is_a_number=slope['derivative_is_a_number'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
