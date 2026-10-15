"""S089 checks for the four inputs named by S088.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-declared-inputs.md. This checker does not march
the rule-D clock and does not edit the pinned checkers it calls.
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
import curved_fermion_boundary as curved
import finite_local_electromagnetic_input as em
import source_clock_joint as joint
import source_clock_radiation as rad


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'notes' / 'source-clock-declared-inputs.md'
CURVED_NOTE = ROOT / 'notes' / 'curved-fermion-boundary.md'
COUPLED_NOTE = ROOT / 'notes' / 'coupled-fermion-domain.md'
MOUTH_NOTE = ROOT / 'notes' / 'finite-mouth-determinant.md'
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
        'The candidate is the shipped mixed a4 curvature coefficient',
        '`sigma` is not replaced by that derivative',
        'Assumption. The folded coefficient is the shipped one-angle a4',
        'Assumption. The values named by the finite-mouth note',
        'Assumption. A measured susceptibility is the shipped parallel susceptibility',
        'The S085 march is not repeated',
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


def second_fundamental_form():
    names = s.symbols('L00 L11 L22 L01 L02 L12')
    matrix = s.Matrix([
        [names[0], names[3], names[4]],
        [names[3], names[1], names[5]],
        [names[4], names[5], names[2]],
    ])
    return names, matrix


def scalar_quadratic(matrix, vector):
    total = 0
    for i in range(3):
        for j in range(3):
            total += vector[i] * matrix[i, j] * vector[j]
    return s.expand(total)


def shipped_a4():
    if getattr(shipped_a4, 'value', None) is None:
        info = curved.clifford_checks()
        rows = info['coefficients']
        assert len(rows) == 2
        parsed = [s.expand(s.sympify(row['a4'])) for row in rows]
        zero(parsed[0] - parsed[1])
        shipped_a4.value = parsed[0]
    return shipped_a4.value


def curvature_contact():
    """Derivative of the shipped a4 curvature coefficient in one L component."""
    coefficient = shipped_a4()
    names, matrix = second_fundamental_form()
    gradient = s.Matrix(s.symbols('u0 u1 u2'))
    expected = s.Rational(4, 5) * (
        scalar_quadratic(matrix, gradient)
        - s.trace(matrix) * scalar_quadratic(s.eye(3), gradient) / 3
    )
    zero(coefficient - expected)
    derivative = s.diff(coefficient, names[0])
    point_u = {gradient[0]: 1, gradient[1]: 0, gradient[2]: 0}
    point_l = {
        names[1]: s.Rational(1, 3),
        names[2]: s.Rational(1, 5),
        names[3]: s.Rational(1, 7),
        names[4]: s.Rational(1, 11),
        names[5]: s.Rational(1, 13),
    }
    exact = s.factor(derivative.subs(point_u).subs(point_l))
    assert exact != 0
    center = s.Rational(2, 7)

    def at(value):
        return coefficient.subs(point_u).subs(point_l).subs(names[0], value)

    zero(stencil(at, center) - exact)
    scale = s.symbols('lam', real=True)
    umbilic = {
        names[0]: scale,
        names[1]: scale,
        names[2]: scale,
        names[3]: 0,
        names[4]: 0,
        names[5]: 0,
    }
    zero(coefficient.subs(umbilic))
    other = s.factor(coefficient.subs(point_u).subs({
        names[0]: 1,
        names[1]: 0,
        names[2]: 0,
        names[3]: 0,
        names[4]: 0,
        names[5]: 0,
    }))
    assert other != 0
    declared = joint.charge_identities()['Lambda']
    assert 'zeta_s' in declared and 'sigma' in declared
    assert s.Symbol('sigma') not in s.expand(derivative).free_symbols
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        sigma_set_to_zero=False,
        zeta_s_set_to_zero=False,
        a4=str(s.factor(coefficient)),
        l00_derivative=str(s.factor(derivative)),
        derivative_at_rational_gradient=str(exact),
        umbilic_value=0,
        nonumbilic_value=str(other),
        sigma_replaced_by_derivative=False,
        Lambda=declared,
        march_called=False,
    )


def folded_heat():
    """The shipped one-angle coefficient gains a rate product only by substitution."""
    info = coupled.analytic()
    assert info['product_Delta_a4_unsmeared'] == 0
    normal = s.sympify(info['product_Delta_a4_normal_smearing_without_4pi'])
    theta_p = s.Symbol('theta_prime')
    constant = s.simplify(normal.subs(theta_p, 1))
    zero(normal - constant * theta_p**2)
    assert constant != 0
    left, right = s.symbols('theta_L theta_R', real=True)
    substituted = s.expand(normal.subs(theta_p, left - right))
    product = substituted.coeff(left).coeff(right)
    zero(product - (-2 * constant))
    point = {right: s.Rational(1, 3)}
    center = s.Rational(2, 5)
    exact = s.diff(substituted, left).subs(point).subs(left, center)

    def at(value):
        return substituted.subs(point).subs(left, value)

    zero(stencil(at, center) - exact)
    coefficient = shipped_a4()
    names, matrix = second_fundamental_form()
    left_g = s.Matrix(s.symbols('p0:3'))
    right_g = s.Matrix(s.symbols('q0:3'))

    def evaluate(vector):
        return coefficient.subs({
            s.Symbol('u0'): vector[0],
            s.Symbol('u1'): vector[1],
            s.Symbol('u2'): vector[2],
        })

    polarization = s.expand(
        evaluate(left_g + right_g) - evaluate(left_g) - evaluate(right_g)
    ) / 2
    scale = s.symbols('lam', real=True)
    umbilic = {
        names[0]: scale,
        names[1]: scale,
        names[2]: scale,
        names[3]: 0,
        names[4]: 0,
        names[5]: 0,
    }
    zero(polarization.subs(umbilic))
    require_sentence(COUPLED_NOTE, 'are not covered by (8).')
    require_sentence(
        CURVED_NOTE,
        'Unequal collar perturbations are outside the following theorem.',
    )
    require_sentence(
        COUPLED_NOTE,
        'but not unequal jets, unequal clocks or arbitrary metric perturbations.',
    )
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        positive_quadratic_mass_used_as_result=False,
        shipped_unsmeared_a4=0,
        shipped_normal_smearing=info['product_Delta_a4_normal_smearing_without_4pi'],
        rate_product_coefficient=str(s.factor(product)),
        umbilic_polarization=0,
        folded_coefficient_proved=False,
    )


def mouth_scaling():
    """The published scaling depends on q, and the mean return stays a big-O."""
    channel = s.Symbol('q', positive=True)
    formulas = {
        'a': 1 / (4 * s.pi**2 * channel),
        'epsilon': channel**(-2),
        'lambda_0': 1 / (16 * s.pi * channel**3),
    }
    derivatives = {}
    gaps = {}
    center = s.Integer(4)
    for name, formula in formulas.items():
        exact = s.factor(s.diff(formula, channel))
        assert exact != 0
        step = s.symbols('h', positive=True)
        quotient = (
            formula.subs(channel, channel + step)
            - formula.subs(channel, channel - step)
        ) / (2 * step)
        zero(s.limit(quotient, step, 0) - exact)

        def at(value, formula=formula):
            return formula.subs(channel, value)

        exact_at = exact.subs(channel, center)
        gap = s.simplify((stencil(at, center) - exact_at) / exact_at)
        assert abs(float(gap)) < 1e-12
        derivatives[name] = str(exact)
        gaps[name] = str(gap)
    note = flat(MOUTH_NOTE)
    assert 'a=1/(4pi^2 q), epsilon=q^-2, delta=log^-2 q,' in note
    assert 'average D_j=O(delta log q), lambda_0=1/(16pi q^3),' in note
    signature = list(inspect.signature(rad.march).parameters)
    assert signature == ['eta', 'lam', 'kappa', 'z0', 'v0', 'method']
    declared = joint.charge_identities()['Lambda']
    for name in ('q', 'a', 'epsilon', 'D_j', 'lambda_0'):
        assert name not in signature
    assert re.search(r'(?<![A-Za-z0-9_])q(?![A-Za-z0-9_])', declared) is None
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        wormhole_determinant_computed=False,
        channel_count_assigned=False,
        q_derivatives=derivatives,
        finite_difference_relative_gaps=gaps,
        average_D_j='O(delta log q)',
        average_D_j_is_a_value=False,
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def numeric_background():
    """A numeric B0 leaves chi dependent on cJ and kappa."""
    info = em.invariants()
    susceptibility = s.sympify(
        info['parallel_electric_susceptibility'].replace('^', '**')
    )
    background = s.Symbol('B0')
    coupling = s.Symbol('cJ')
    kappa = s.Symbol('kappa')
    first = s.expand(susceptibility.subs(background, 1))
    second = s.expand(susceptibility.subs(background, s.Rational(1, 2)))
    assert coupling in first.free_symbols and kappa in first.free_symbols
    assert coupling in second.free_symbols and kappa in second.free_symbols
    assert s.simplify(first - second) != 0
    assert not first.is_number and not second.is_number
    ratio = s.simplify((first - 1) / (second - 1))
    zero(ratio - (1 / s.Rational(1, 2))**2)
    exact = s.diff(first, coupling)
    point = {kappa: s.Rational(3, 2)}
    exact_value = s.factor(exact.subs(point))
    center = s.Rational(1, 5)

    def at(value):
        return first.subs(point).subs(coupling, value)

    zero(stencil(at, center) - exact_value)
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
        chi_at_B0_1=str(s.factor(first)),
        chi_at_B0_one_half=str(s.factor(second)),
        chi_is_a_number=False,
        march_kappa_used_as_gravitational_kappa=False,
        cJ_derivative_at_B0_1=str(s.factor(exact)),
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    curvature = curvature_contact()
    heat = folded_heat()
    scaling = mouth_scaling()
    magnetic = numeric_background()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S088',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s087_product_throat_recomputed=False,
        s088_rn_series_recomputed=False,
        new_lambda_interval=False,
        area_budget_evaluated=False,
        area_budget='not evaluated; no replacement target was written before the test',
        zeta_s_set_to_zero=False,
        sigma_set_to_zero=False,
        cJ_set_to_zero=False,
        delta_set_to_zero=False,
        endpoint_energy_cancelled=False,
        positive_quadratic_mass_used_as_result=False,
        tasks=dict(
            curvature_contact=curvature,
            folded_heat=heat,
            mouth_scaling=scaling,
            numeric_background=magnetic,
        ),
        pinned_checker_sha256={
            path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in PINNED
        },
        versions=dict(
            python=platform.python_version(),
            sympy=s.__version__,
        ),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n',
    )
    print(json.dumps(dict(
        status=result['status'],
        curvature_contact=curvature['classification'],
        derivative_at_rational_gradient=curvature['derivative_at_rational_gradient'],
        folded_heat=heat['classification'],
        rate_product_coefficient=heat['rate_product_coefficient'],
        mouth_scaling=scaling['classification'],
        numeric_background=magnetic['classification'],
        chi_at_B0_1=magnetic['chi_at_B0_1'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
