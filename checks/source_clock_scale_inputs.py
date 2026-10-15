"""S090 checks for the four inputs named by S089.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-scale-inputs.md. This checker does not march the
rule-D clock and does not edit the pinned checkers it calls.
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
NOTE = ROOT / 'notes' / 'source-clock-scale-inputs.md'
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
        'Assumption. A renormalized value is a number for `c_theta`',
        '`sigma` is not replaced by `b_theta`',
        'Assumption. An error bound would compare the one-angle',
        'Assumption. The clock datum would be one declared positive `q`',
        'Assumption. A measurement is a number for `chi`',
        'The S085 march is not repeated',
        'The S089 `L00` derivative is not recomputed as a new result',
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


def scalar_quadratic(matrix, vector):
    total = 0
    for i in range(3):
        for j in range(3):
            total += vector[i] * matrix[i, j] * vector[j]
    return s.expand(total)


def curvature_form():
    names = s.symbols('L00 L11 L22 L01 L02 L12')
    matrix = s.Matrix([
        [names[0], names[3], names[4]],
        [names[3], names[1], names[5]],
        [names[4], names[5], names[2]],
    ])
    gradient = s.Matrix(s.symbols('u0 u1 u2'))
    bare = s.expand(
        scalar_quadratic(matrix, gradient)
        - s.trace(matrix) * scalar_quadratic(s.eye(3), gradient) / 3
    )
    return names, gradient, bare


def renormalized_coupling():
    """b_theta is a number, and c_theta remains inside sigma."""
    info = curved.clifford_checks()
    rows = info['coefficients']
    assert len(rows) == 2
    parsed = [s.expand(s.sympify(row['a4'])) for row in rows]
    zero(parsed[0] - parsed[1])
    coefficient = parsed[0]
    names, gradient, bare = curvature_form()
    zero(coefficient - s.simplify(coefficient.subs({
        names[0]: 1, names[1]: 0, names[2]: 0,
        names[3]: 0, names[4]: 0, names[5]: 0,
        gradient[0]: 1, gradient[1]: 0, gradient[2]: 0,
    })) * bare / bare.subs({
        names[0]: 1, names[1]: 0, names[2]: 0,
        names[3]: 0, names[4]: 0, names[5]: 0,
        gradient[0]: 1, gradient[1]: 0, gradient[2]: 0,
    }))
    sample = {
        names[0]: 1, names[1]: 0, names[2]: 0,
        names[3]: 0, names[4]: 0, names[5]: 0,
        gradient[0]: 1, gradient[1]: 0, gradient[2]: 0,
    }
    prefactor = s.simplify(coefficient.subs(sample) / bare.subs(sample))
    assert prefactor.is_number and prefactor != 0
    zero(coefficient - prefactor * bare)
    beta = s.factor(prefactor / (4 * s.pi)**2)
    assert beta.is_number
    assert s.Symbol('c_theta') not in beta.free_symbols
    declared = joint.charge_identities()
    assert declared['sigma'] == '2*c_theta*Lhat^{uu}'
    assert declared['sigma_set_to_zero'] is False
    coupling, collar = s.symbols('c_theta Lhat')
    sigma = 2 * coupling * collar
    exact = s.diff(sigma, coupling)
    center = s.Rational(1, 4)
    point = {collar: s.Rational(3, 5)}

    def at(value):
        return sigma.subs(point).subs(coupling, value)

    zero(stencil(at, center) - exact.subs(point))
    assert collar in exact.free_symbols
    signature = list(inspect.signature(rad.march).parameters)
    for name in ('c_theta', 'mu'):
        assert name not in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        c_theta_set_to_zero=False,
        sigma_set_to_zero=False,
        a4_prefactor=str(prefactor),
        b_theta=str(beta),
        sigma=declared['sigma'],
        c_theta_derivative=str(exact),
        sigma_replaced_by_b_theta=False,
        march_parameters=signature,
        march_called=False,
    )


def substitution_error():
    """The substitution difference is not a two-collar error bound."""
    info = coupled.analytic()
    assert info['product_Delta_a4_unsmeared'] == 0
    normal = s.sympify(info['product_Delta_a4_normal_smearing_without_4pi'])
    angle = s.Symbol('theta_prime')
    constant = s.simplify(normal.subs(angle, 1))
    zero(normal - constant * angle**2)
    left, right = s.symbols('theta_L theta_R')
    substituted = s.expand(constant * (left - right)**2)
    difference = s.expand(substituted - normal.subs(angle, left))
    assert difference != 0
    point = {right: s.Rational(1, 3)}
    center = s.Rational(2, 5)
    exact = s.diff(difference, left).subs(point).subs(left, center)

    def at(value):
        return difference.subs(point).subs(left, value)

    zero(stencil(at, center) - exact)
    keys = sorted(info)
    assert not any('bound' in key or 'two_collar' in key for key in keys)
    require_sentence(COUPLED_NOTE, 'are not covered by (8).')
    require_sentence(
        CURVED_NOTE,
        'Unequal collar perturbations are outside the following theorem.',
    )
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        positive_quadratic_mass_used_as_result=False,
        shipped_return_keys=keys,
        substitution_difference=str(s.factor(difference)),
        two_collar_error_bound_returned=False,
        error_bound_proved=False,
    )


def channel_data():
    """At one rational q the three formulas are numbers and D_j is not."""
    channel = s.Symbol('q', positive=True)
    formulas = {
        'a': 1 / (4 * s.pi**2 * channel),
        'epsilon': channel**(-2),
        'lambda_0': 1 / (16 * s.pi * channel**3),
    }
    illustrated = {}
    center = s.Integer(4)
    step = s.symbols('h', positive=True)
    for name, formula in formulas.items():
        value = s.factor(formula.subs(channel, center))
        assert value != 0
        quotient = (
            formula.subs(channel, channel + step)
            - formula.subs(channel, channel - step)
        ) / (2 * step)
        zero(s.limit(quotient, step, 0) - s.diff(formula, channel))
        illustrated[name] = str(value)
    note = flat(MOUTH_NOTE)
    assert 'average D_j=O(delta log q), lambda_0=1/(16pi q^3),' in note
    signature = list(inspect.signature(rad.march).parameters)
    declared = joint.charge_identities()['Lambda']
    for name in ('q', 'D_j'):
        assert name not in signature
        assert re.search(
            r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
            declared,
        ) is None
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        channel_count_assigned=False,
        illustrated_q='4',
        illustrated_values=illustrated,
        average_D_j='O(delta log q)',
        numeric_D_j=None,
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def measured_pair():
    """Numbers for B0 and kappa leave chi linear in cJ."""
    info = em.invariants()
    susceptibility = s.sympify(
        info['parallel_electric_susceptibility'].replace('^', '**')
    )
    coupling = s.Symbol('cJ')
    substituted = s.expand(
        susceptibility.subs({s.Symbol('B0'): 1, s.Symbol('kappa'): 1})
    )
    assert coupling in substituted.free_symbols
    assert s.Symbol('B0') not in substituted.free_symbols
    assert s.Symbol('kappa') not in substituted.free_symbols
    assert not substituted.is_number
    exact = s.diff(substituted, coupling)
    assert exact.is_number and exact != 0
    center = s.Rational(1, 5)

    def at(value):
        return substituted.subs(coupling, value)

    zero(stencil(at, center) - exact)
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
        chi_at_unit_B0_and_kappa=str(s.factor(substituted)),
        cJ_derivative=str(exact),
        chi_is_a_number=False,
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    coupling = renormalized_coupling()
    error = substitution_error()
    channel = channel_data()
    magnetic = measured_pair()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S089',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s087_product_throat_recomputed=False,
        s088_rn_series_recomputed=False,
        s089_l00_derivative_recomputed=False,
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
            renormalized_coupling=coupling,
            substitution_error=error,
            channel_data=channel,
            measured_pair=magnetic,
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
        renormalized_coupling=coupling['classification'],
        b_theta=coupling['b_theta'],
        a4_prefactor=coupling['a4_prefactor'],
        substitution_error=error['classification'],
        channel_data=channel['classification'],
        illustrated_values=channel['illustrated_values'],
        measured_pair=magnetic['classification'],
        chi_at_unit_B0_and_kappa=magnetic['chi_at_unit_B0_and_kappa'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
