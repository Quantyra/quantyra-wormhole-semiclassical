"""S086 checks for the four inputs named after the S085 clock obstructions.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-separated-inputs.md. This checker does not march the
rule-D clock, does not re-bracket Lambda, and does not edit the joint,
radiation, controller, allocation, curved-fermion, or finite-local checkers.
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
NOTE = ROOT / 'notes' / 'source-clock-separated-inputs.md'
CURVED_NOTE = ROOT / 'notes' / 'curved-fermion-boundary.md'
COUPLED_NOTE = ROOT / 'notes' / 'coupled-fermion-domain.md'
MOUTH_NOTE = ROOT / 'notes' / 'finite-mouth-determinant.md'
ANOMALY_NOTE = ROOT / 'notes' / 'spherical-anomaly-sources.md'

POINT = dict(
    N=s.Rational(3, 2),
    R=s.Rational(5, 2),
    theta_dot=s.Rational(2, 5),
    zeta_s=s.Rational(7, 5),
    c_theta=s.Rational(3, 5),
    Lhat=s.Rational(4, 5),
)
DET_POINT = dict(z=0, v=s.Rational(1, 100), m=s.Rational(1, 2))
STENCIL = 1e-6


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
        'Assumption. The coordinate component',
        'Error comparison. The collar tensor',
        'The clock charge depends on',
        'Hardest step. At',
        'Observable. At',
        'Call `contact_variations()`',
        'Assumption. Each collar',
        'Assumption. The free resolvent',
        'The predeclared point is `z = 0`, `v = 1/100`, `m = 1/2`',
        'Assumption. The matching condition is the shipped parallel',
        'The S085 march is not repeated',
    ):
        assert sentence in text, sentence


def stencil_derivative(samples, step):
    """Fourth-order derivative from f(x-2h), f(x-h), f(x+h), f(x+2h)."""
    return (-samples[3] + 8 * samples[2] - 8 * samples[1] + samples[0]) / (12 * step)


def relative_gap(numeric, exact):
    scale = max(abs(exact), 1e-12)
    return abs(numeric - exact) / scale


def separation():
    """Lapse response at fixed Lhat^{uu}, calibrated against published T_s."""
    N, R, th, zeta, coeff, lhat, sigma, ell, d_lhat = s.symbols(
        'N R theta_dot zeta_s c_theta Lhat sigma ell dLhat_dN',
        real=True, positive=True,
    )
    pi = s.pi
    lag_s = 2 * pi * zeta * R**2 * th**2 / N
    lag_th = 4 * pi * coeff * N * R**2 * lhat * th**2
    lag_sum = 2 * pi * R**2 * th**2 * (zeta + sigma) / N
    sigma_def = 2 * coeff * lhat
    round_metric = {N: 1, R: ell, sigma: sigma_def}
    zero((lag_s + lag_th - lag_sum).subs(round_metric))
    zero(lag_th.subs({N: 1, R: ell}) - 4 * pi * ell**2 * coeff * lhat * th**2)

    radius_gap = s.diff(lag_s + lag_th, R) - s.diff(lag_sum, R)
    zero(radius_gap.subs(round_metric))
    lapse_gap = s.simplify(
        (s.diff(lag_s + lag_th, N) - s.diff(lag_sum, N)).subs(round_metric)
    )
    zero(lapse_gap - 4 * pi * ell**2 * th**2 * sigma_def)
    assert s.simplify(lapse_gap) != 0
    assert zeta in (lag_s + lag_th).free_symbols
    assert lhat in lag_th.free_symbols

    def ratio(lag):
        return s.simplify(R * s.diff(lag, R) / (N * s.diff(lag, N)))

    ratio_s = ratio(lag_s)
    ratio_th = ratio(lag_th)
    assert s.simplify(ratio_s - ratio_th) != 0
    combined = s.simplify(
        ratio(lag_s + lag_th).subs({N: 1, lhat: sigma / (2 * coeff)})
    )
    clock_ratio = s.simplify(ratio(lag_sum).subs(N, 1))
    assert s.simplify(combined - clock_ratio) != 0
    shift = s.symbols('a', real=True)
    shifted = combined.subs({zeta: zeta + shift, sigma: sigma - shift})
    assert s.simplify(shifted - combined) != 0
    zero(shifted.subs(shift, 0) - combined)
    declared = joint.charge_identities()
    assert declared['Lambda'] == (
        'lambda_zeta+lambda_theta=(pi*ell**2/N)*(zeta_s+sigma)'
    )
    assert declared['sigma_set_to_zero'] is False
    coupled = (pi * ell**2 / s.symbols('number', positive=True)) * (zeta + sigma)
    zero(coupled.subs(zeta, zeta + shift).subs(sigma, sigma - shift) - coupled)

    # Published T_s. delta I = (1/2) integral sqrt(-gamma) T^{ij} delta gamma_ij.
    density = -th**2 / N**2
    t_tt = zeta * (th**2 - (-N**2) * density / 2)
    angular_trace = -zeta * density
    t_tt_up = t_tt / N**4
    measure = 4 * pi * N * R**2
    lapse_from_stress = s.Rational(1, 2) * measure * t_tt_up * (-2 * N)
    radius_from_stress = s.Rational(1, 2) * measure * (2 / R) * angular_trace
    zero(lapse_from_stress - s.diff(lag_s, N))
    zero(radius_from_stress - s.diff(lag_s, R))
    energy = s.simplify(
        (zeta * th**2 / (2 * N**2) * measure).subs({N: 1, R: ell})
    )
    zero(energy - 2 * pi * zeta * ell**2 * th**2)

    extra = s.diff(lag_th, lhat) * d_lhat
    assert d_lhat in extra.free_symbols
    assert s.simplify(extra) != 0

    values = {
        N: POINT['N'], R: POINT['R'], th: POINT['theta_dot'],
        zeta: POINT['zeta_s'], coeff: POINT['c_theta'], lhat: POINT['Lhat'],
    }

    def numeric_pair(analytic_n, analytic_r, builder):
        gaps = []
        for symbol, analytic in ((N, analytic_n), (R, analytic_r)):
            samples = []
            center = float(values[symbol])
            for sign, scale in ((-2, 1), (-1, 1), (1, 1), (2, 1)):
                moved = dict(values)
                moved[symbol] = center + sign * scale * STENCIL
                samples.append(float(builder.subs(moved)))
            # samples were appended in order -2,-1,+1,+2. Stencil wants -2,-1,+1,+2.
            numeric = stencil_derivative(samples, STENCIL)
            gaps.append(relative_gap(numeric, float(analytic.subs(values))))
        return max(gaps)

    source_gap = numeric_pair(s.diff(lag_s, N), s.diff(lag_s, R), lag_s)
    contact_gap = numeric_pair(s.diff(lag_th, N), s.diff(lag_th, R), lag_th)
    assert source_gap < 1e-8, source_gap
    assert contact_gap < 1e-8, contact_gap

    cases = curved.contact_variations()
    assert len(cases) == 8
    worst = max(case['discrepancy'] for case in cases)
    assert worst < 2e-7

    require_sentence(CURVED_NOTE, 'L_ab=-<nabla_a n,e_b>')
    return dict(
        classification='uncontrolled approximation',
        admissible_range=None,
        identified_physical_input=False,
        zeta_s_set_to_zero=False,
        sigma_set_to_zero=False,
        round_metric_actions_agree=True,
        round_metric_radius_derivatives_agree=True,
        source_response_ratio=str(ratio_s),
        fixed_Lhat_response_ratio=str(ratio_th),
        clock_sum_response_ratio=str(clock_ratio),
        fixed_Lhat_combined_ratio=str(combined),
        opposite_shift_changes_combined_ratio=True,
        lapse_gap='4*pi*ell**2*theta_dot**2*sigma',
        lhat_lapse_piece=str(s.factor(extra)),
        lhat_derivative_evaluated=False,
        source_finite_difference_relative_gap=source_gap,
        contact_finite_difference_relative_gap=contact_gap,
        shipped_contact_variation_max_discrepancy=worst,
        endpoint_energy_cancelled=False,
    )


def unequal_collars():
    """Average of two published energies. The jet cross term stays unevaluated."""
    zeta, th, ell_bar, delta = s.symbols(
        'zeta_s theta_dot ell_bar delta', real=True, positive=True,
    )
    pi = s.pi
    ell_l = ell_bar + delta / 2
    ell_r = ell_bar - delta / 2

    def energy(radius):
        return 2 * pi * zeta * radius**2 * th**2

    average = s.expand((energy(ell_l) + energy(ell_r)) / 2)
    matched = energy(ell_bar)
    zero(average.subs(delta, 0) - matched)
    zero((energy(ell_l) + energy(ell_r)).subs(delta, 0) - 2 * matched)
    second = s.diff(average, delta, 2)
    expanded_second = s.diff(s.expand(average), delta, 2)
    zero(second - expanded_second)
    assert s.simplify(second) != 0
    jet = s.Function('jet_cross')
    with_jet = average + jet(delta)
    gained = s.simplify(s.diff(with_jet, delta, 2) - second)
    assert gained != 0
    assert gained.has(s.Derivative)
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
        classification='uncontrolled approximation',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        two_collar_equation_introduced=False,
        matched_limit='average equals E_s(ell_bar) at delta=0',
        sum_of_energies_at_zero_difference='2*E_s',
        second_derivative=str(s.factor(second)),
        jet_cross_second_derivative='retained',
        positive_quadratic_mass_used_as_result=False,
    )


def interacting_determinant():
    """The free derivative does not bound the interacting derivative."""
    z, v, m = s.symbols('z v m', real=True)
    factor = 1 + z
    free = 2 * s.log(1 + factor * v)
    det = (1 + factor * v)**2 - (factor * m)**2
    interacting = s.log(det)
    d_free = s.simplify(s.diff(free, z))
    d_int = s.simplify(s.diff(interacting, z))
    assert m not in d_free.free_symbols
    assert m in d_int.free_symbols
    zero(d_int.subs(m, 0) - d_free)
    second = s.factor(s.diff(d_int, m, 2).subs(m, 0))
    assert second != 0
    point = {z: DET_POINT['z'], v: DET_POINT['v'], m: DET_POINT['m']}
    free_value = s.together(d_free.subs(point))
    int_value = s.together(d_int.subs(point))
    assert s.Abs(int_value) > s.Abs(free_value)

    def at(z_value):
        gg = 1.0 + z_value
        vv = float(DET_POINT['v'])
        mm = float(DET_POINT['m'])
        return np.log((1.0 + gg * vv)**2 - (gg * mm)**2)

    samples = [at(sign * scale * STENCIL) for sign, scale in (
        (-2, 1), (-1, 1), (1, 1), (2, 1),
    )]
    numeric = stencil_derivative(samples, STENCIL)
    gap = relative_gap(numeric, float(int_value))
    assert gap < 1e-8, gap
    require_sentence(
        MOUTH_NOTE,
        "control the corrected physical background's nonlocal determinant",
    )
    require_sentence(
        ANOMALY_NOTE,
        'It does not supply the independent nonlocal determinant',
    )
    return dict(
        classification='uncontrolled approximation',
        admissible_range=None,
        identified_physical_input=False,
        determinant_set_to_one=False,
        free_gaussian_recomputed=False,
        wormhole_determinant_computed=False,
        predeclared_point='z=0, v=1/100, m=1/2',
        free_derivative=str(free_value),
        interacting_derivative=str(int_value),
        second_m_derivative_at_zero=str(second),
        finite_difference_relative_gap=gap,
    )


def cJ_match():
    """Solve the shipped susceptibility. The benchmark does not supply it."""
    info = em.invariants()
    raw = info['parallel_electric_susceptibility']
    expr = s.sympify(raw.replace('^', '**'))
    cj = s.Symbol('cJ')
    chi = s.Symbol('chi')
    solved_list = s.solve(s.Eq(chi, expr), cj)
    assert len(solved_list) == 1
    solved = s.factor(solved_list[0])
    zero(expr.subs(cj, solved) - chi)
    assert s.Symbol('Lambda') not in solved.free_symbols
    assert chi in solved.free_symbols
    assert s.Symbol('B0') in solved.free_symbols
    declared = joint.charge_identities()
    assert 'cJ' not in declared['Lambda']
    zeta, sigma, ell, number = s.symbols('zeta_s sigma ell N', positive=True)
    coupled = (s.pi * ell**2 / number) * (zeta + sigma)
    assert cj not in coupled.free_symbols
    signature = list(inspect.signature(rad.march).parameters)
    assert signature == ['eta', 'lam', 'kappa', 'z0', 'v0', 'method']
    for absent in ('cJ', 'chi', 'B0', 'susceptibility'):
        assert absent not in signature
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        cJ_set_to_zero=False,
        conditional_axial_window_adopted=False,
        shipped_parallel_susceptibility=raw,
        matching_condition=str(solved),
        matching_data='parallel electric susceptibility and B0',
        data_supplied_by_rule_D_march=False,
        march_parameters=signature,
        quadratic_correction=info['quadratic_correction'],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    separated = separation()
    collars = unequal_collars()
    determinant = interacting_determinant()
    electromagnetic = cJ_match()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_after_S085',
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
        endpoint_energy_cancelled=False,
        positive_quadratic_mass_used_as_result=False,
        tasks=dict(
            separation=separated,
            unequal_collars=collars,
            nonlocal_interacting_determinant=determinant,
            cJ=electromagnetic,
        ),
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
        separation=separated['classification'],
        unequal_collars=collars['classification'],
        nonlocal_interacting_determinant=determinant['classification'],
        cJ=electromagnetic['classification'],
        opening_achieved=False,
        s074_gate_closed=False,
        s085_march_repeated=False,
    )), flush=True)


if __name__ == '__main__':
    main()
