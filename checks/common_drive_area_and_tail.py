"""S136.4 conditional free-JT area and post-source tail identities.

Frozen before first run: exact symbolic checks, 60 decimal digits, independent
tail/Green quadratures within 1e-35. Examples are mathematical verification
profiles; no physical remainder assignment or preparation-box certificate.
"""
import hashlib
import json
import platform
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'evidence/common-drive-area-and-tail.json'
DPS = 60
QUADRATURE_TOL = '1e-35'
INPUTS = [
    'notes/common-drive-area-and-tail.md',
    'checks/common_drive_area_and_tail.py',
    'notes/jt-prescribed-source-constraint.md',
    'notes/laboratory-drive-clock-feedback.md',
    'notes/physical-clock-remainder.md',
]


def require_zero(expression):
    result = sp.trigsimp(sp.simplify(expression))
    if result != 0:
        result = sp.trigsimp(sp.simplify(sp.expand_trig(result)))
    if result != 0:
        raise AssertionError(('nonzero symbolic remainder', result))


def symbolic_checks():
    s, r, v = sp.symbols('s r v', real=True)
    A, R, B = sp.symbols('A R B', positive=True)
    phi = sp.Function('Phi')(s)
    require_zero(sp.diff(sp.sin(s)*phi, s)/sp.sin(s)
                 - sp.diff(phi, s)-sp.cot(s)*phi)
    require_zero(-sp.diff(sp.diff(phi, s)+sp.cot(s)*phi, s)
                 - (-sp.diff(phi, s, 2)-sp.cot(s)*sp.diff(phi, s)
                    + phi/sp.sin(s)**2))
    # K is continuous; its flux Phi_s+cot(s)*Phi jumps from 1 to 0.
    # Thus -d_s(flux K)=delta(s-r), including the correct sign.
    k_left = (1-sp.cos(s))/sp.sin(s)  # s<r
    k_right = (1-sp.cos(r))/sp.sin(s)  # s>r
    require_zero(k_left.subs(s, r)-k_right.subs(s, r))
    require_zero(sp.diff(k_left, s)+sp.cot(s)*k_left-1)
    require_zero(sp.diff(k_right, s)+sp.cot(s)*k_right)
    require_zero(sp.limit(sp.sin(s)*k_left, s, 0))
    require_zero(sp.diff(sp.sin(s)*k_right, s).subs(s, sp.pi/2))
    # Positive-denominator cross multiplication proves the half-angle
    # representation sin(s)/(1+cos(s))=tan(s/2), without asking SymPy to
    # select a half-angle simplification branch.
    require_zero((1-sp.cos(s))*(1+sp.cos(s))-sp.sin(s)**2)

    witness = -A*sp.sin(s)**4
    Ttt = sp.expand_trig(-sp.diff(witness, s, 2)
                        -sp.cot(s)*sp.diff(witness, s)
                        +witness/sp.sin(s)**2)
    Tss = -sp.cot(s)*sp.diff(witness, s)-witness/sp.sin(s)**2
    require_zero(Ttt-5*A*sp.sin(s)**2*(4*sp.cos(s)**2-1))
    require_zero(Tss-A*sp.sin(s)**2*(4*sp.cos(s)**2+1))
    require_zero(sp.integrate(Ttt, (s, 0, sp.pi)))
    require_zero(sp.diff(Tss, s)+sp.cot(s)*(Tss-Ttt))
    require_zero(sp.sin(s)**2*(-Ttt+Tss)
                 - (sp.sin(s)**2*sp.diff(witness, s, 2)-2*witness))
    require_zero(sp.limit(sp.sin(s)*witness, s, 0))
    require_zero(sp.limit(sp.diff(witness, s)+sp.cot(s)*witness, s, 0))
    require_zero(witness.subs(s, sp.pi/2)+A)
    resonant = R*v*sp.sin(v)/2
    require_zero(sp.diff(resonant, v, 2)+resonant-R*sp.cos(v))
    require_zero(resonant.subs(v, 0))
    require_zero(sp.diff(resonant, v).subs(v, 0))
    n = sp.symbols('n', integer=True, nonnegative=True)
    require_zero(resonant.subs(v, 3*sp.pi/2+2*sp.pi*n)
                 + R*(3*sp.pi/2+2*sp.pi*n)/2)

    eta, m, reference_bound = sp.Rational(1, 2000), sp.Rational(3, 200), sp.Rational(13, 10)
    qmax = 4*eta*reference_bound/m
    dmax = 3*qmax
    assert qmax == sp.Rational(13, 75)
    assert dmax == sp.Rational(13, 25) < 1
    assert 1-dmax == sp.Rational(12, 25)
    assert m-sp.Rational(1, 2) == -sp.Rational(97, 200)
    # Squaring positive sides proves 1/2-sqrt(47/200)>3/200.
    floor_square_gap = (sp.Rational(1, 2)-m)**2-sp.Rational(47, 200)
    assert floor_square_gap == sp.Rational(9, 40000) > 0
    assert m-sp.Rational(1, 200) == sp.Rational(1, 100)
    # Leading classical ratio is monotone for k>=2000, X>=k^3.
    classical_margin = sp.Integer(2000)**10/4-sp.Rational(97, 200)
    assert classical_margin > 0
    return {
        'spatial_operator_and_distributional_green_kernel': 'exact symbolic pass',
        'kernel_range_proof': 'K<=tan(s/2)<=1 on 0<s<=pi/2; K>=0',
        'zero_charge_witness_Ttt': '5*A*sin(s)^2*(4*cos(s)^2-1)',
        'zero_charge_witness_Tss': 'A*sin(s)^2*(4*cos(s)^2+1)',
        'zero_charge_witness_integral_boundary_and_conservation': 'exact symbolic pass',
        'resonant_forcing_and_unbounded_negative_subsequence': 'exact symbolic pass',
        'exact_constants': dict(eta=str(eta), qmax=str(qmax), Dmax=str(dmax),
                                gf_min=str(1-dmax), area_floor='-97/200',
                                free_clock_floor_square_gap=str(floor_square_gap),
                                tail_L1='1/200', future_clock_floor='>1/100',
                                classical_margin_at_k_2000=str(classical_margin)),
    }


def text_number(value):
    return mp.nstr(value, 55)


def tail_checks():
    def exponential_response(a, t):
        return (mp.exp(a*t)-mp.cos(t)-a*mp.sin(t))/(a*a+1)

    examples = [
        ('exp_decay', lambda x: mp.exp(-x)/1000, mp.mpf(1)/1000,
         lambda t: exponential_response(-1, t)/1000),
        ('modulated_decay', lambda x: mp.exp(-2*x)*(1+mp.cos(x))/1600,
         mp.mpf(1)/800/2,
         lambda t: (exponential_response(-2, t)
                    +mp.re(exponential_response(-2+1j, t)))/1600),
        ('signed_gaussian', lambda x: mp.exp(-x*x)*mp.sin(3*x)/2500,
         mp.sqrt(mp.pi)/5000, None),
    ]
    rows, max_error = [], mp.mpf(0)
    for name, residual, envelope_L1, closed in examples:
        assert envelope_L1 < mp.mpf(1)/200
        for t in [mp.mpf('0.2'), mp.mpf('1.7'), mp.mpf(8)]:
            direct = mp.quad(lambda x: mp.sin(t-x)*residual(x), [0, t], method='tanh-sinh')
            velocity = mp.quad(lambda x: mp.cos(t-x)*residual(x), [0, t], method='tanh-sinh')
            A = mp.quad(lambda x: mp.cos(x)*residual(x), [0, t], method='gauss-legendre')
            B = mp.quad(lambda x: mp.sin(x)*residual(x), [0, t], method='gauss-legendre')
            alternate = mp.sin(t)*A-mp.cos(t)*B
            alternate_velocity = mp.cos(t)*A+mp.sin(t)*B
            errors = [abs(direct-alternate), abs(velocity-alternate_velocity)]
            if closed is not None:
                errors += [abs(direct-closed(t)), abs(velocity-mp.diff(closed, t))]
            error = max(errors)
            assert error <= mp.mpf(QUADRATURE_TOL), (name, t, error)
            assert abs(direct) <= envelope_L1 and abs(velocity) <= envelope_L1
            max_error = max(max_error, error)
            rows.append(dict(profile=name, time=text_number(t),
                             delta_z=text_number(direct), delta_v=text_number(velocity),
                             rigorous_envelope_L1=text_number(envelope_L1),
                             independent_error=text_number(error)))
    return dict(rows=rows, maximum_error=text_number(max_error),
                profile_scope='smooth verification residuals, no physical assignment')


def green_checks():
    s = sp.symbols('s', real=True)
    cases = [
        ('negative_sine_squared', sp.Rational(3, 100), -sp.Rational(2, 7)*sp.sin(s)**2),
        ('zero_charge_witness', sp.Rational(1, 9), -sp.Rational(3, 5)*sp.sin(s)**4),
        ('mixed_smooth_profile', sp.Rational(1, 7), sp.Rational(4, 11)*sp.sin(s)**2*sp.cos(s)**2),
    ]
    rows, max_error = [], mp.mpf(0)
    for name, boundary, regular in cases:
        # Homogeneous B/sin(s) has zero source. Only regular terms are used
        # in the source/flux functions, avoiding cancellation at s=0.
        T = sp.trigsimp(-sp.diff(regular, s, 2)-sp.cot(s)*sp.diff(regular, s)
                       +regular/sp.sin(s)**2)
        flux = sp.trigsimp(sp.diff(regular, s)+sp.cot(s)*regular)
        require_zero(sp.diff(flux, s)+T)
        require_zero(flux.subs(s, sp.pi/2))
        stress = sp.lambdify(s, T, 'mpmath')
        flux_value = sp.lambdify(s, flux, 'mpmath')
        phi = sp.lambdify(s, boundary/sp.sin(s)+regular, 'mpmath')
        bv = mp.mpf(str(boundary.p))/mp.mpf(str(boundary.q))
        for sigma in [mp.pi/13, mp.pi/5, mp.pi/2]:
            # Direct Green integral, split at the derivative discontinuity.
            first = mp.quad(lambda x: (1-mp.cos(x))*stress(x), [0, sigma], method='tanh-sinh')
            second = (1-mp.cos(sigma))*mp.quad(stress, [sigma, mp.pi/2], method='tanh-sinh')
            green = (bv+first+second)/mp.sin(sigma)
            # Independent first-order spatial reconstruction and quadrature.
            nested_reduced = (bv+mp.quad(lambda x: mp.sin(x)*flux_value(x),
                                        [0, sigma], method='gauss-legendre'))/mp.sin(sigma)
            expected = phi(sigma)
            error = max(abs(green-expected), abs(nested_reduced-expected), abs(green-nested_reduced))
            assert error <= mp.mpf(QUADRATURE_TOL), (name, sigma, error)
            max_error = max(max_error, error)
            rows.append(dict(profile=name, sigma=text_number(sigma),
                             reconstructed_phi=text_number(green), independent_error=text_number(error)))
    return dict(rows=rows, maximum_error=text_number(max_error),
                profile_scope='smooth metric-equation verification sources, no matter realization')


def normalized_hash(path):
    raw = (ROOT/path).read_bytes().replace(b'\r\n', b'\n').replace(b'\r', b'\n')
    return hashlib.sha256(raw).hexdigest()


def main():
    mp.mp.dps = DPS
    result = dict(
        status='pass', milestone='E003/S136.4',
        scientific_verdict='conditional free-JT area/tail control; bounded persistent residual is insufficient',
        frozen_protocol=dict(decimal_digits=DPS, independent_quadrature_tolerance='1e-35',
                             eta='1/2000', driving_clock_target='3/200',
                             completion_target='7', endpoint_energy_target='-3/200',
                             postcompletion_tail_L1='1/200'),
        symbolic=symbolic_checks(), tail=tail_checks(), spatial_green=green_checks(),
        claim_limits=dict(full_preparation_rectangle_certified=False,
                          driving_and_endpoint_targets_assumed_not_proved=True,
                          area_bound_for_arbitrary_clock_residual=False,
                          free_area_requires_consistent_free_JT_r_zero=True,
                          physical_tail_bound_derived=False,
                          renormalized_extremal_area_or_full_physical_area_certified=False,
                          physical_opening_or_stability_claimed=False, Lean_run=False),
        input_hash_policy='SHA-256 of UTF-8 source bytes after CRLF/CR normalization to LF; no other normalization',
        inputs_sha256={path: normalized_hash(path) for path in INPUTS},
        runtime=dict(python=sys.version, platform=platform.platform(),
                     mpmath=mp.__version__, sympy=sp.__version__),
    )
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=result['status'], milestone=result['milestone'],
                          tail_maximum_error=result['tail']['maximum_error'],
                          green_maximum_error=result['spatial_green']['maximum_error'],
                          receipt=str(OUTPUT)), indent=2))


if __name__ == '__main__':
    main()
