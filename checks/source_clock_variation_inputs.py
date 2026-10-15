"""S096 checks for the four inputs named by S095.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-variation-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-variation-inputs.md'
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
        'Assumption. A number for `c_theta(mu0)` would be the coefficient `-3` in the Weyl variation of `J_theta`.',
        'That coefficient is not a value of `c_theta(mu0)`.',
        'Assumption. The unequal-collar invariant would be the traced `S_:aa` term.',
        'Independence of `theta` is not a two-collar theorem.',
        'Assumption. A number for one `D_j` would be the statement that both derivatives of `E_cov` are `O(1)=o(q delta)`.',
        'A little-o comparison is not an assignment.',
        'Assumption. A measured `cJ` would be the ratio of the second `B0` derivative to its `cJ` derivative at gravitational `kappa = 1`.',
        'Returning the symbol `cJ` does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S095 `J_theta` derivative is not recomputed as a new result',
        'The S095 curvature trace is not recomputed as a new result',
        'The S095 local-bound derivative is not recomputed as a new result',
        'The S095 `cJ` derivative is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def weyl_coefficient():
    """The Weyl counting factor -3 does not assign c_theta(mu0)."""
    require_sentence(CURVED_NOTE, 'J_theta=integral sqrt(gamma)L u^2')
    require_sentence(
        CURVED_NOTE,
        'changes by -3 integral sqrt(gamma)(n omega)u^2',
    )
    assert re.search(r'c_theta\(mu0\)\s*=\s*-3', flat(CURVED_NOTE)) is None
    density = s.Symbol('density')
    variation = -3 * density
    exact = s.diff(variation, density)
    zero(exact + 3)
    center = s.Integer(2)

    def at(value):
        return variation.subs(density, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    probe = s.Rational(3, 7)
    numeric = (at(probe + STEP) - at(probe - STEP)) / (2 * STEP)
    assert abs(s.N(numeric - exact, 40)) < s.Integer(10)**(-18)
    assert exact.is_number
    assert density not in exact.free_symbols
    constant, origin = s.symbols('c_theta_mu0 mu0')
    endpoint = constant + s.log(origin) / (20 * s.pi**2)
    assert constant in endpoint.free_symbols
    assert not endpoint.is_number
    assert s.simplify(endpoint + 3) != 0
    declared = joint.charge_identities()
    assert declared['sigma'] == '2*c_theta*Lhat^{uu}'
    signature = list(inspect.signature(rad.march).parameters)
    assert 'c_theta' not in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        c_theta_set_to_zero=False,
        c_theta_set_equal_to_minus_three=False,
        weyl_coefficient='-3',
        density_derivative=str(s.factor(exact)),
        derivative_assigns_c_theta=False,
        endpoint_keeps_c_theta_mu0=True,
        endpoint_equals_minus_three=False,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def traced_term():
    """A vanishing theta derivative is not the unequal-collar mismatch."""
    require_sentence(
        CURVED_NOTE,
        'tr S=-2L is independent of theta, so the traced S_:aa term adds none.',
    )
    require_sentence(
        COUPLED_NOTE,
        "tr(SE), tr(S^3), tr(S_:aa), tr((chi')^2 S) vanish",
    )
    require_sentence(
        COUPLED_NOTE,
        'Unequal/background-dependent collars require their additional invariants and are not covered by (8).',
    )
    angle = s.Symbol('theta')
    length = s.Symbol('L')
    displayed = -2 * length
    displayed_exact = s.diff(displayed, angle)
    zero(displayed_exact)

    def displayed_at(value):
        return displayed.subs(angle, value)

    displayed_sampled = (
        displayed_at(s.Integer(1) + STEP) - displayed_at(s.Integer(1) - STEP)
    ) / (2 * STEP)
    zero(s.expand(displayed_sampled - displayed_exact))
    rate = s.Symbol('theta_prime')
    chi = (
        s.cos(angle) * coupled.chi0
        - s.sin(angle) * coupled.kron(coupled.p1, coupled.I2)
    )
    projector = (coupled.I4 + chi) / 2
    generator = coupled.C * projector * rate / 2
    raw = s.trace(generator)
    shipped = s.simplify(s.trigsimp(s.expand(raw)))
    shipped_exact = s.simplify(s.trigsimp(s.expand(s.diff(raw, angle))))
    zero(s.simplify(shipped_exact - s.diff(shipped, angle)))
    assert s.simplify(shipped - displayed) != 0
    assert s.simplify(displayed_exact - 1) != 0
    assert s.simplify(shipped_exact - 1) != 0
    center = s.Rational(1, 2)
    rate_value = s.Rational(1, 3)
    prepared = raw.subs(rate, rate_value)

    def shipped_at(value):
        return prepared.subs(angle, value)

    sampled = (
        shipped_at(center + STEP) - shipped_at(center - STEP)
    ) / (2 * STEP)
    expected = s.diff(prepared, angle).subs(angle, center)
    assert abs(s.N(s.simplify(sampled - expected), 40)) < s.Integer(10)**(-12)
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        trace_adopted_as_theorem=False,
        shipped_trace_identified_with_minus_2L=False,
        displayed_theta_derivative=str(displayed_exact),
        shipped_theta_derivative=str(s.factor(shipped_exact)),
        derivatives_equal_mismatch=False,
        recorded_mismatch='1',
    )


def little_o_bound():
    """The little-o comparison does not assign one D_j."""
    require_sentence(
        MOUTH_NOTE,
        'both derivatives of E_cov are O(1)=o(q delta).',
    )
    charge, width = s.symbols('q delta')
    product = charge * width
    exact = s.diff(product, charge)
    zero(s.simplify(exact - width))
    assert width in exact.free_symbols
    assert not exact.is_number
    center = s.Integer(2)

    def at(value):
        return product.subs(charge, value)

    sampled = (at(center + STEP) - at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    probe = s.Rational(3, 7)
    numeric = sampled.subs(width, probe)
    expected = exact.subs(width, probe)
    assert abs(s.N(numeric - expected, 40)) < s.Integer(10)**(-18)
    declared = joint.charge_identities()['Lambda']
    rhs = s.sympify(declared.split('=', 1)[1], locals={'N': s.Symbol('N')})
    clock_symbol = s.Symbol('D_j')
    zero(s.diff(rhs, clock_symbol))
    zero(s.simplify(rhs.subs(clock_symbol, 1) - rhs.subs(clock_symbol, 2)))
    signature = list(inspect.signature(rad.march).parameters)
    joined = ' '.join(signature)
    for name in ('D_j', 'q', 'delta'):
        assert absent(name, joined)
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        delta_set_to_zero=False,
        channel_count_assigned=False,
        little_o_evaluated_as_equality=False,
        little_o_is_a_number=False,
        numeric_D_j=None,
        q_delta_derivative=str(s.factor(exact)),
        derivative_assigns_D_j=False,
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def susceptibility_ratio():
    """The ratio of the two B0 derivatives returns the symbol cJ."""
    info = em.invariants()
    require_sentence(EM_NOTE, 'Magnetic background matching alone cannot')
    susceptibility = s.sympify(
        info['parallel_electric_susceptibility'].replace('^', '**')
    )
    background = s.Symbol('B0')
    coupling = s.Symbol('cJ')
    grav = s.Symbol('kappa')
    second = s.diff(susceptibility, background, 2).subs(grav, 1)
    slope = s.diff(second, coupling)
    ratio = s.simplify(second / slope)
    zero(s.simplify(ratio - coupling))
    assert coupling in ratio.free_symbols
    assert not ratio.is_number
    probe = s.Rational(3, 7)
    numeric = second.subs(coupling, probe) / slope
    expected = ratio.subs(coupling, probe)
    assert abs(s.N(numeric - expected, 40)) < s.Integer(10)**(-18)
    ratio_slope = s.diff(ratio, coupling)
    zero(s.simplify(ratio_slope - 1))
    center = s.Integer(2)
    sampled = (
        ratio.subs(coupling, center + STEP) - ratio.subs(coupling, center - STEP)
    ) / (2 * STEP)
    zero(s.expand(sampled - ratio_slope))
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
        ratio=str(s.factor(ratio)),
        ratio_is_a_number=False,
        ratio_measures_cJ=False,
        ratio_slope_stored_as_cJ_measurement=False,
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    coefficient = weyl_coefficient()
    traced = traced_term()
    bound = little_o_bound()
    ratio = susceptibility_ratio()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S095',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s095_j_theta_derivative_recomputed=False,
        s095_curvature_trace_recomputed=False,
        s095_local_bound_derivative_recomputed=False,
        s095_cJ_derivative_recomputed=False,
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
            weyl_coefficient=coefficient,
            traced_term=traced,
            little_o_bound=bound,
            susceptibility_ratio=ratio,
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
        weyl_coefficient=coefficient['classification'],
        density_derivative=coefficient['density_derivative'],
        endpoint_equals_minus_three=coefficient['endpoint_equals_minus_three'],
        traced_term=traced['classification'],
        displayed_theta_derivative=traced['displayed_theta_derivative'],
        shipped_theta_derivative=traced['shipped_theta_derivative'],
        little_o_bound=bound['classification'],
        q_delta_derivative=bound['q_delta_derivative'],
        little_o_is_a_number=bound['little_o_is_a_number'],
        susceptibility_ratio=ratio['classification'],
        ratio=ratio['ratio'],
        ratio_is_a_number=ratio['ratio_is_a_number'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
