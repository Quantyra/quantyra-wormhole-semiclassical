"""S123 checks for the four inputs named by S122.

The comparison, hardest step, observable, and test are fixed in
notes/source-clock-nineteenth-inputs.md. This checker does not march the
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
NOTE = ROOT / 'notes' / 'source-clock-nineteenth-inputs.md'
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
EXPECTED_PINNED = {
    'curved_fermion_boundary.py': '0aed0bba4cdf72d120edec075f13854b68914eed8779b782806cdd20169c1240',
    'coupled_fermion_domain.py': '692cfa606cc38c89f975766d2db76f665e218410f51418ccf7b1e8a8d6772a91',
    'finite_local_electromagnetic_input.py': '8d8c040a4ebdb4ed15d267a6bf5d5f1e793da94ba6363b40c3ff7c252100f1cd',
    'source_clock_joint.py': 'a17adbc5b1acb8327a89a0d1fb1524b7dbf88493318f8164bd9a5549d37eecb3',
    'source_clock_radiation.py': 'bbfe392c627630d570c04a7b1132314aab132cf1997ef606ae57e9410542b101',
}
MARCH = ('eta', 'lam', 'kappa', 'z0', 'v0', 'method')
STEP = s.Rational(1, 1000)
# Probe at h=2/5: step 1/10**20 on 12804747411456000/h**19 has error about 9.68e-13, under 1e-12.
NUMERIC_STEP = s.Rational(1, 10**20)


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def flat(path):
    return ' '.join(path.read_text(encoding='utf-8').split())


def require_sentence(path, sentence):
    assert sentence in flat(path), sentence


def require_pre_fit_note():
    for sentence in (
        'before any checker run',
        'Assumption. A number for `c_theta(mu0)` would be the nineteenth `mu0` derivative `320118685286400/(pi**2*mu0**19)`.',
        'That derivative is not a value of `c_theta(mu0)`.',
        'Assumption. The unequal-collar invariant would be the `f` derivative `711374856192000/f**18` of `-41845579776000/f**17`.',
        'That derivative is not a two-collar theorem.',
        'Assumption. A number for one `D_j` would be the `h` derivative `-243290200817664000/h**20`.',
        'That derivative is not an assignment of `D_j`.',
        'Assumption. A measured `cJ` would be the `B0` derivative `0` of the eighteenth derivative `0`.',
        'The derivative `0` does not measure `cJ`.',
        'The S085 march is not repeated',
        'The S122 eighteenth mu0 derivative is not recomputed as a new result',
        'The S122 f derivative is not recomputed as a new result',
        'The S122 h derivative is not recomputed as a new result',
        'The S122 eighteenth B0 derivative is not recomputed as a new result',
        'The S088 relative series is not differentiated again',
        'No replacement area target is written before the test',
        'A positive quadratic term is not used as a stability mass',
        '`cJ` is not set to zero',
    ):
        require_sentence(NOTE, sentence)


def absent(name, text):
    return re.search(
        r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
        text,
    ) is None


def agree_limit(term, symbol, exact, center):
    step_symbol = s.Symbol('step')

    def at(value):
        return term.subs(symbol, value)

    central = (at(center + step_symbol) - at(center - step_symbol)) / (2 * step_symbol)
    target = exact.subs(symbol, center)
    zero(s.expand(s.limit(central, step_symbol, 0) - target))
    sampled = (at(center + NUMERIC_STEP) - at(center - NUMERIC_STEP)) / (2 * NUMERIC_STEP)
    assert abs(s.N(sampled - target, 40)) < s.Integer(10)**(-12)


def nineteenth_mu0_derivative():
    """The nineteenth mu0 derivative does not assign c_theta."""
    require_sentence(CURVED_NOTE, 'b_theta=4/[5(4pi)^2]')
    scheme = s.Integer(4) / (5 * (4 * s.pi)**2)
    zero(scheme - 1 / (20 * s.pi**2))
    origin = s.Symbol('mu0')
    term = -17784371404800 / (s.pi**2 * origin**18)
    exact = s.diff(term, origin)
    expected = 320118685286400 / (s.pi**2 * origin**19)
    zero(s.expand(exact - expected))
    stored = str(s.factor(exact))
    assert stored == '320118685286400/(pi**2*mu0**19)'
    assert stored != '-17784371404800/(pi**2*mu0**18)'
    assert not exact.is_number
    assert origin in exact.free_symbols
    center = s.Rational(4, 9)
    assert center != 1
    agree_limit(term, origin, exact, center)
    endpoint_coupling = s.Symbol('c_theta_mu0')
    endpoint = endpoint_coupling + s.log(origin) / (20 * s.pi**2)
    assert endpoint_coupling in endpoint.free_symbols
    assert not endpoint.is_number
    assert s.simplify(exact - endpoint_coupling) != 0
    declared = joint.charge_identities()
    assert declared['sigma'] == '2*c_theta*Lhat^{uu}'
    signature = list(inspect.signature(rad.march).parameters)
    assert tuple(signature) == MARCH
    assert absent('c_theta', ' '.join(signature))
    return dict(
        classification='unavailable physical input',
        admissible_range=None,
        identified_physical_input=False,
        c_theta_set_to_zero=False,
        c_theta_set_equal_to_b_theta=False,
        derivative_is_a_number=False,
        derivative_assigns_c_theta=False,
        nineteenth_mu0_derivative=stored,
        s122_eighteenth_mu0_derivative_recomputed=False,
        endpoint_keeps_c_theta_mu0=True,
        sigma=declared['sigma'],
        march_parameters=signature,
        march_called=False,
    )


def eighteenth_lapse_derivative():
    """The f derivative of -41845579776000/f**17 is not the collar mismatch."""
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
    zero(b_scheme - 1 / (20 * s.pi**2))
    zero(d_scheme - 1 / (60 * s.pi**2))
    lapse = s.Symbol('f')
    term = -41845579776000 / lapse**17
    exact = s.diff(term, lapse)
    expected = 711374856192000 / lapse**18
    zero(s.expand(exact - expected))
    stored = str(s.factor(exact))
    assert stored == '711374856192000/f**18'
    assert stored != '1'
    assert stored != '-41845579776000/f**17'
    assert not exact.is_number
    assert lapse in exact.free_symbols
    center = s.Rational(4, 9)
    assert center != 1
    agree_limit(term, lapse, exact, center)
    return dict(
        classification='unresolved proof',
        admissible_range=None,
        identified_physical_input=False,
        delta_set_to_zero=False,
        d_theta_set_to_zero=False,
        two_collar_equation_introduced=False,
        derivative_adopted_as_theorem=False,
        rn_series_recomputed=False,
        s122_f_derivative_recomputed=False,
        f_derivative=stored,
        derivative_equals_mismatch=False,
        recorded_mismatch='1',
    )


def width_twentieth():
    """The h derivative of 12804747411456000/h**19 does not assign one D_j."""
    require_sentence(
        MOUTH_NOTE,
        'eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.',
    )
    width = s.Symbol('h')
    term = 12804747411456000 / width**19
    exact = s.diff(term, width)
    expected = -243290200817664000 / width**20
    zero(s.expand(exact - expected))
    stored = str(s.factor(exact))
    assert stored == '-243290200817664000/h**20'
    assert stored != '12804747411456000/h**19'
    assert width in exact.free_symbols
    assert not exact.is_number
    center = s.Rational(2, 5)
    assert center != 0
    agree_limit(term, width, exact, center)
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
    assert tuple(signature) == MARCH
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
        s122_h_derivative_recomputed=False,
        numeric_D_j=None,
        h_derivative=stored,
        h_derivative_in_D_j='0',
        clock_derivative_in_D_j='0',
        march_parameters=signature,
        march_called=False,
        Lambda=declared,
    )


def nineteenth_background_derivative():
    """The B0 derivative of the eighteenth derivative 0 does not measure cJ."""
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
    eighteenth = s.diff(substituted, background, 18)
    zero(s.expand(eighteenth))
    exact = s.diff(eighteenth, background)
    zero(s.expand(exact))
    stored = str(s.factor(exact))
    assert stored == '0'
    assert exact.is_number
    assert coefficient not in exact.free_symbols
    assert stored != 'cJ'
    center = s.Rational(2, 5)
    assert center != 1

    def at(value):
        return substituted.subs(background, value)

    def second_at(value):
        return (at(value + STEP) - 2 * at(value) + at(value - STEP)) / STEP**2

    def third_at(value):
        return (second_at(value + STEP) - second_at(value - STEP)) / (2 * STEP)

    def fourth_at(value):
        return (third_at(value + STEP) - third_at(value - STEP)) / (2 * STEP)

    def fifth_at(value):
        return (fourth_at(value + STEP) - fourth_at(value - STEP)) / (2 * STEP)

    def sixth_at(value):
        return (fifth_at(value + STEP) - fifth_at(value - STEP)) / (2 * STEP)

    def seventh_at(value):
        return (sixth_at(value + STEP) - sixth_at(value - STEP)) / (2 * STEP)

    def eighth_at(value):
        return (seventh_at(value + STEP) - seventh_at(value - STEP)) / (2 * STEP)

    def ninth_at(value):
        return (eighth_at(value + STEP) - eighth_at(value - STEP)) / (2 * STEP)

    def tenth_at(value):
        return (ninth_at(value + STEP) - ninth_at(value - STEP)) / (2 * STEP)

    def eleventh_at(value):
        return (tenth_at(value + STEP) - tenth_at(value - STEP)) / (2 * STEP)

    def twelfth_at(value):
        return (eleventh_at(value + STEP) - eleventh_at(value - STEP)) / (2 * STEP)

    def thirteenth_at(value):
        return (twelfth_at(value + STEP) - twelfth_at(value - STEP)) / (2 * STEP)

    def fourteenth_at(value):
        return (thirteenth_at(value + STEP) - thirteenth_at(value - STEP)) / (2 * STEP)

    def fifteenth_at(value):
        return (fourteenth_at(value + STEP) - fourteenth_at(value - STEP)) / (2 * STEP)

    def sixteenth_at(value):
        return (fifteenth_at(value + STEP) - fifteenth_at(value - STEP)) / (2 * STEP)

    def seventeenth_at(value):
        return (sixteenth_at(value + STEP) - sixteenth_at(value - STEP)) / (2 * STEP)

    def eighteenth_at(value):
        return (seventeenth_at(value + STEP) - seventeenth_at(value - STEP)) / (2 * STEP)

    sampled = (eighteenth_at(center + STEP) - eighteenth_at(center - STEP)) / (2 * STEP)
    zero(s.expand(sampled - exact))
    assert abs(s.N(sampled, 40)) < s.Integer(10)**(-18)
    coupling_slope = s.diff(exact, coefficient)
    zero(s.expand(coupling_slope))
    signature = list(inspect.signature(rad.march).parameters)
    assert tuple(signature) == MARCH
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
        s122_eighteenth_B0_derivative_recomputed=False,
        nineteenth_B0_derivative=stored,
        derivative_is_cJ=False,
        derivative_measures_cJ=False,
        cJ_derivative_of_nineteenth_derivative='0',
        march_kappa_used_as_gravitational_kappa=False,
        march_parameters=signature,
        march_called=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require_pre_fit_note()
    pinned = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in PINNED
    }
    assert pinned == EXPECTED_PINNED
    endpoint = nineteenth_mu0_derivative()
    lapse = eighteenth_lapse_derivative()
    width = width_twentieth()
    slope = nineteenth_background_derivative()
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_for_the_four_inputs_named_by_S122',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        successor_model_opened=False,
        s085_march_repeated=False,
        s088_rn_series_recomputed=False,
        s122_eighteenth_mu0_derivative_recomputed=False,
        s122_f_derivative_recomputed=False,
        s122_h_derivative_recomputed=False,
        s122_eighteenth_B0_derivative_recomputed=False,
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
            nineteenth_mu0_derivative=endpoint,
            eighteenth_lapse_derivative=lapse,
            width_twentieth=width,
            nineteenth_background_derivative=slope,
        ),
        pinned_checker_sha256=pinned,
        versions=dict(python=platform.python_version(), sympy=s.__version__),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n',
    )
    print(json.dumps(dict(
        status=result['status'],
        nineteenth_mu0_derivative=endpoint['nineteenth_mu0_derivative'],
        f_derivative=lapse['f_derivative'],
        h_derivative=width['h_derivative'],
        nineteenth_B0_derivative=slope['nineteenth_B0_derivative'],
        opening_achieved=False,
    )), flush=True)


if __name__ == '__main__':
    main()
