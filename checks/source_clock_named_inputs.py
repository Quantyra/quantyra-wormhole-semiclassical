"""S087 checks for the four inputs named by S086.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-named-inputs.md. This checker does not march the
rule-D clock and does not edit the pinned checkers it calls.
"""
import argparse
import hashlib
import inspect
import json
import platform
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import sympy as s

import curved_fermion_boundary as curved
import finite_local_electromagnetic_input as em
import source_clock_joint as joint
import source_clock_radiation as rad


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'notes' / 'source-clock-named-inputs.md'
CURVED_NOTE = ROOT / 'notes' / 'curved-fermion-boundary.md'
COUPLED_NOTE = ROOT / 'notes' / 'coupled-fermion-domain.md'
MOUTH_NOTE = ROOT / 'notes' / 'finite-mouth-determinant.md'
NORM_NOTE = ROOT / 'notes' / 'nonlinear-current-determinant.md'
PINNED = (
    ROOT / 'checks' / 'source_clock_radiation.py',
    ROOT / 'checks' / 'source_clock_joint.py',
    ROOT / 'checks' / 'curved_fermion_boundary.py',
    ROOT / 'checks' / 'finite_local_electromagnetic_input.py',
)
STENCIL = 1e-6
DET_POINT = dict(z=0, v=s.Rational(1, 100), m=s.Rational(1, 2))


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def flat(path):
    return ' '.join(path.read_text(encoding='utf-8').split())


def require_sentence(path, sentence):
    text = flat(path)
    assert sentence in text, sentence


def require_pre_fit_note():
    text = flat(NOTE)
    for sentence in (
        'before any checker run',
        'The retained cutoff is the product throat',
        'gamma_tt = -r_e^2/epsilon^2',
        'N = r_e/epsilon',
        '`sigma` is not replaced by that limit',
        'Assumption. The jet-cross term',
        'Equation (8) of the coupled-domain note has no collar difference',
        'Assumption. The norm is the published Schatten bound',
        'The predeclared point remains `z = 0`, `v = 1/100`, `m = 1/2`',
        'a bound on a value does not become a bound on a derivative',
        'Assumption. The zero magnetic background',
        'The S085 march is not repeated',
        'No replacement area target is written before the test',
    ):
        assert sentence in text, sentence


def relative_gap(numeric, exact):
    scale = max(abs(exact), 1e-12)
    return abs(numeric - exact) / scale


def lapse_derivative():
    """Differentiate the shipped product-throat weight with respect to N."""
    info = rad.geometry_identities()
    density = s.sympify(info['integrated_density'])
    eps = s.Symbol('epsilon')
    assert eps in density.free_symbols
    re = s.Symbol('r_e', positive=True)
    lapse = s.Symbol('N', positive=True)
    weight = s.simplify(density.subs(eps, re / lapse))
    d_n = s.simplify(s.diff(weight, lapse))
    chain = s.simplify(
        s.diff(density, eps).subs(eps, re / lapse) * s.diff(re / lapse, lapse)
    )
    zero(d_n - chain)
    d_eps_0 = s.factor(s.diff(density, eps).subs(eps, 0))
    assert d_eps_0 != 0
    assert s.limit(weight, lapse, s.oo) == 0
    assert s.limit(d_n, lapse, s.oo) == 0
    sample = {
        s.Symbol('Sch'): s.Rational(1, 7),
        s.Symbol('rate'): s.Rational(2, 5),
        re: s.Rational(5, 2),
    }
    center = s.Rational(3, 2)
    exact = complex(d_n.subs(sample).subs(lapse, center))
    assert abs(exact.imag) < 1e-15

    def at(value):
        return complex(weight.subs(sample).subs(lapse, value))

    samples = [at(center + sign * scale * STENCIL) for sign, scale in (
        (-2, 1), (-1, 1), (1, 1), (2, 1),
    )]
    numeric = (-samples[3] + 8 * samples[2] - 8 * samples[1] + samples[0]) / (12 * STENCIL)
    gap = relative_gap(numeric.real, exact.real)
    assert gap < 1e-8, gap
    declared = joint.charge_identities()
    lam = declared['Lambda']
    assert 'zeta_s' in lam and 'sigma' in lam
    assert 'B0' not in lam
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        sigma_set_to_zero=False,
        zeta_s_set_to_zero=False,
        shipped_integrated_density=info['integrated_density'],
        epsilon_derivative_at_zero=str(d_eps_0),
        lapse_derivative=str(s.factor(d_n)),
        clock_boundary_limit_of_weight=0,
        clock_boundary_limit_of_lapse_derivative=0,
        finite_difference_relative_gap=gap,
        Lambda=lam,
        sigma_replaced_by_boundary_limit=False,
        march_called=False,
    )


def jet_cross():
    """Average the published quadratic jet. The fold stays unproved."""
    delta = s.symbols('delta', real=True)
    mean = s.Matrix(s.symbols('b0:3', real=True))
    jump = s.Matrix(s.symbols('j0:3', real=True))

    def square(vec):
        return s.expand((vec.T * vec)[0])

    left = mean + delta * jump / 2
    right = mean - delta * jump / 2
    average = s.expand((square(left) + square(right)) / 2)
    zero(s.simplify(average.subs(delta, 0) - square(mean)))
    second = s.factor(s.diff(average, delta, 2))
    zero(second - square(jump) / 2)
    assert second != 0
    point = {symbol: s.Rational(index + 1, 5) for index, symbol in enumerate(mean)}
    point.update({symbol: s.Rational(index + 2, 7) for index, symbol in enumerate(jump)})
    center = s.Rational(1, 3)
    exact = second.subs(point)
    step = s.Rational(1, 1000)
    samples = [
        average.subs(point).subs(delta, center + sign * scale * step)
        for sign, scale in ((-2, 1), (-1, 1), (0, 0), (1, 1), (2, 1))
    ]
    numeric = s.simplify(
        (-samples[4] + 16 * samples[3] - 30 * samples[2] + 16 * samples[1] - samples[0])
        / (12 * step**2)
    )
    zero(numeric - exact)
    float_step = 1e-3
    float_samples = [
        complex(average.subs(point).subs(delta, complex(center) + sign * scale * float_step))
        for sign, scale in ((-2, 1), (-1, 1), (0, 0), (1, 1), (2, 1))
    ]
    float_numeric = (
        -float_samples[4] + 16 * float_samples[3] - 30 * float_samples[2]
        + 16 * float_samples[1] - float_samples[0]
    ) / (12 * float_step**2)
    gap = relative_gap(float_numeric.real, float(exact))
    assert gap < 1e-8, gap
    clifford = curved.clifford_checks()
    assert clifford['arbitrary_symmetric_second_fundamental_form'] == 'PASS'
    require_sentence(
        CURVED_NOTE,
        'Unequal collar perturbations are outside the following theorem.',
    )
    require_sentence(
        COUPLED_NOTE,
        'Unequal/background-dependent collars require their additional invariants '
        'and are not covered by (8).',
    )
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        positive_quadratic_mass_used_as_result=False,
        second_derivative=str(second),
        finite_difference_relative_gap=gap,
        folded_coefficient_proved=False,
    )


def determinant_norm():
    """Apply the published Schatten bound to the predeclared stand-in."""
    z, v, m = s.symbols('z v m', positive=True)
    singular = (1 + z) * m / (1 + (1 + z) * v)
    radius = s.together(2 * singular**3)
    logdet3 = s.log(1 - singular**2) + singular**2
    bound = radius / (3 * (1 - radius ** s.Rational(1, 3)))
    point = {z: DET_POINT['z'], v: DET_POINT['v'], m: DET_POINT['m']}
    # z=0 is on the boundary of the positive symbol. Substitute after building.
    numeric_point = {z: s.Integer(0), v: DET_POINT['v'], m: DET_POINT['m']}
    radius_value = s.together(radius.subs(numeric_point))
    assert radius_value < 1
    value = s.Abs(logdet3.subs(numeric_point)).evalf(40)
    ceiling = bound.subs(numeric_point).evalf(40)
    assert value <= ceiling
    free = 2 * s.log(1 + (1 + z) * v)
    full = s.log((1 + (1 + z) * v) ** 2 - ((1 + z) * m) ** 2)
    d_free = s.simplify(s.diff(free, z).subs(numeric_point))
    d_full = s.simplify(s.diff(full, z).subs(numeric_point))
    assert d_full != d_free

    def at(z_value):
        vv = float(DET_POINT['v'])
        mm = float(DET_POINT['m'])
        factor = 1.0 + z_value
        return np.log((1.0 + factor * vv) ** 2 - (factor * mm) ** 2)

    samples = [at(sign * scale * STENCIL) for sign, scale in (
        (-2, 1), (-1, 1), (1, 1), (2, 1),
    )]
    numeric = (-samples[3] + 8 * samples[2] - 8 * samples[1] + samples[0]) / (12 * STENCIL)
    gap = relative_gap(numeric, float(d_full))
    assert gap < 1e-8, gap
    require_sentence(
        MOUTH_NOTE,
        "control the corrected physical background's nonlocal determinant",
    )
    require_sentence(
        NORM_NOTE,
        '|log det_3(I+X)| <= r/[3(1-r^(1/3))].',
    )
    return dict(
        classification='uncontrolled approximation',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        wormhole_determinant_computed=False,
        free_gaussian_recomputed=False,
        bound_controls_z_derivative=False,
        predeclared_point='z=0, v=1/100, m=1/2',
        schatten_radius=str(radius_value),
        schatten_hypothesis='r<1',
        logdet3_absolute_value=str(value),
        published_bound=str(ceiling),
        free_derivative=str(s.together(d_free)),
        interacting_derivative=str(s.together(d_full)),
        finite_difference_relative_gap=gap,
    )


def magnetic_background():
    """The shipped zero-background identities leave cJ free."""
    info = em.invariants()
    susceptibility = s.sympify(info['parallel_electric_susceptibility'].replace('^', '**'))
    correction = s.sympify(info['quadratic_correction'].replace('^', '**'))
    background = s.Symbol('B0')
    coefficient = s.Symbol('cJ')
    reduced = s.simplify(susceptibility.subs(background, 0) - 1)
    zero(reduced)
    assert coefficient not in reduced.free_symbols
    cleared = s.simplify(correction.subs(background, 0))
    zero(cleared)
    assert coefficient not in cleared.free_symbols
    declared = joint.charge_identities()
    lam = declared['Lambda']
    assert 'cJ' not in lam and 'B0' not in lam
    signature = list(inspect.signature(rad.march).parameters)
    assert signature == ['eta', 'lam', 'kappa', 'z0', 'v0', 'method']
    for absent in ('cJ', 'chi', 'B0', 'susceptibility'):
        assert absent not in signature
    assert info['magnetic_background_action_and_first_variation'] == 'zero'
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        cJ_set_to_zero=False,
        conditional_axial_window_adopted=False,
        susceptibility_at_zero_background='1',
        quadratic_correction_at_zero_background='0',
        cJ_constrained_by_zero_background=False,
        march_parameters=signature,
        march_called=False,
        Lambda=lam,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    lapse = lapse_derivative()
    jets = jet_cross()
    norm = determinant_norm()
    magnetic = magnetic_background()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S086',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
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
            lapse_derivative=lapse,
            jet_cross=jets,
            determinant_norm=norm,
            magnetic_background=magnetic,
        ),
        pinned_checker_sha256={
            path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in PINNED
        },
        versions=dict(
            python=platform.python_version(),
            sympy=s.__version__,
            numpy=np.__version__,
        ),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n',
    )
    print(json.dumps(dict(
        status=result['status'],
        lapse_derivative=lapse['classification'],
        jet_cross=jets['classification'],
        determinant_norm=norm['classification'],
        magnetic_background=magnetic['classification'],
        epsilon_derivative_at_zero=lapse['epsilon_derivative_at_zero'],
        schatten_radius=norm['schatten_radius'],
        opening_achieved=False,
        s074_gate_closed=False,
        s085_march_repeated=False,
    )), flush=True)


if __name__ == '__main__':
    main()
