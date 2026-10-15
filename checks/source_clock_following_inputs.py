"""S088 checks for the four inputs named by S087.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-following-inputs.md. This checker does not march
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

import curved_fermion_boundary as curved
import finite_local_electromagnetic_input as em
import source_clock_joint as joint
import source_clock_radiation as rad


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'notes' / 'source-clock-following-inputs.md'
CURVED_NOTE = ROOT / 'notes' / 'curved-fermion-boundary.md'
COUPLED_NOTE = ROOT / 'notes' / 'coupled-fermion-domain.md'
MOUTH_NOTE = ROOT / 'notes' / 'finite-mouth-determinant.md'
PINNED = (
    ROOT / 'checks' / 'curved_fermion_boundary.py',
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
        'The candidate is the shipped RN relative type-B series',
        'sigma` is not replaced by that zero',
        'Assumption. The missing invariant is the cross',
        'Assumption. The norm inputs are the symbols named by the finite-mouth',
        'Assumption. A nonzero background is the symbol `B0`',
        'The S085 march is not repeated',
        'No replacement area target is written before the test',
    ):
        require_sentence(NOTE, sentence)


def relative_gap(numeric, exact):
    scale = max(abs(exact), 1e-12)
    return abs(numeric - exact) / scale


def curved_boundary():
    """First derivative of the shipped RN type-B series at the boundary."""
    info = curved.weyl_and_asymptotic_checks()
    series = s.sympify(info['RN_typeB_relative_series'])
    zee = s.Symbol('z')
    assert zee in series.free_symbols
    value = s.simplify(series.subs(zee, 0))
    first = s.factor(s.diff(series, zee).subs(zee, 0))
    second = s.factor(s.diff(series, zee, 2).subs(zee, 0))
    zero(value)
    zero(first)
    assert second != 0
    sample = {symbol: s.Rational(index + 1, 3) for index, symbol in enumerate(
        series.free_symbols - {zee}
    )}
    center = s.Rational(1, 5)
    exact = s.diff(series, zee).subs(sample).subs(zee, center)
    samples = [
        series.subs(sample).subs(zee, center + sign * scale * STENCIL)
        for sign, scale in ((-2, 1), (-1, 1), (1, 1), (2, 1))
    ]
    numeric = s.simplify(
        (-samples[3] + 8 * samples[2] - 8 * samples[1] + samples[0]) / (12 * STENCIL)
    )
    zero(numeric - exact)
    declared = joint.charge_identities()
    lam = declared['Lambda']
    assert 'zeta_s' in lam and 'sigma' in lam
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        sigma_set_to_zero=False,
        zeta_s_set_to_zero=False,
        series=info['RN_typeB_relative_series'],
        boundary_value=0,
        first_derivative_at_boundary=0,
        second_derivative_at_boundary=str(second),
        sigma_replaced_by_boundary_zero=False,
        Lambda=lam,
        march_called=False,
    )


def angle_cross():
    """The two-angle cross is outside equation (8)."""
    left = s.Matrix(s.symbols('p0:3', real=True))
    right = s.Matrix(s.symbols('q0:3', real=True))

    def square(vec):
        return s.expand((vec.T * vec)[0])

    cross = s.expand((left.T * right)[0])
    difference = s.expand(((left - right).T * (left - right))[0])
    zero(difference - square(left) - square(right) + 2 * cross)
    component = right[1]
    fixed = {symbol: s.Rational(index + 1, 4) for index, symbol in enumerate(left)}
    fixed.update({
        symbol: s.Rational(index + 1, 5)
        for index, symbol in enumerate(right) if symbol != component
    })
    exact = s.diff(cross, component).subs(fixed)
    center = s.Rational(2, 5)
    step = s.Rational(1, 1000)

    def at(value):
        return cross.subs(fixed).subs(component, value)

    samples = [at(center + sign * scale * step) for sign, scale in (
        (-2, 1), (-1, 1), (1, 1), (2, 1),
    )]
    numeric = s.simplify(
        (-samples[3] + 8 * samples[2] - 8 * samples[1] + samples[0]) / (12 * step)
    )
    zero(numeric - exact)
    assert cross != 0
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
        cross_term='u_L·u_R',
        folded_coefficient_proved=False,
    )


def mouth_norm():
    """The mouth-note norm inputs are absent from the march."""
    signature = list(inspect.signature(rad.march).parameters)
    assert signature == ['eta', 'lam', 'kappa', 'z0', 'v0', 'method']
    declared = joint.charge_identities()['Lambda']
    absent = ('a', 'epsilon', 'D_j', 'lambda_0')
    for name in absent:
        assert name not in signature
        assert re.search(r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])', declared) is None
    require_sentence(
        MOUTH_NOTE,
        'Here a>0, epsilon, D_j and the exterior time normalization lambda_0',
    )
    require_sentence(
        MOUTH_NOTE,
        "control the corrected physical background's nonlocal determinant",
    )
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        wormhole_determinant_computed=False,
        missing_inputs=list(absent),
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def nonzero_background():
    """A symbolic nonzero B0 and the photon coefficient leave cJ free."""
    info = em.invariants()
    photon = info['four_null_photon_mixed_coefficient_of_J_squared']
    assert photon == '256'
    assert 'cJ' not in photon and 'B0' not in photon
    susceptibility = s.sympify(info['parallel_electric_susceptibility'].replace('^', '**'))
    assert s.Symbol('cJ') in susceptibility.free_symbols
    assert s.Symbol('B0') in susceptibility.free_symbols
    signature = list(inspect.signature(rad.march).parameters)
    for name in ('cJ', 'chi', 'B0'):
        assert name not in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        cJ_set_to_zero=False,
        conditional_axial_window_adopted=False,
        four_null_photon_coefficient=photon,
        cJ_remains_at_nonzero_B0=True,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    boundary = curved_boundary()
    cross = angle_cross()
    norm = mouth_norm()
    magnetic = nonzero_background()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S087',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s087_product_throat_recomputed=False,
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
            curved_boundary=boundary,
            angle_cross=cross,
            mouth_norm=norm,
            nonzero_background=magnetic,
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
        curved_boundary=boundary['classification'],
        second_derivative_at_boundary=boundary['second_derivative_at_boundary'],
        angle_cross=cross['classification'],
        mouth_norm=norm['classification'],
        nonzero_background=magnetic['classification'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
