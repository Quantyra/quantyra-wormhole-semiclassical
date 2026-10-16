"""S137.1 exact rational pressure bracket and reference-centering identities."""
from fractions import Fraction as F
from pathlib import Path
import datetime
import hashlib
import json
import math
import platform
import subprocess
import sys
import traceback

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
NAME = 'quantum-reference-clock-centering'
INPUTS = ['notes/'+NAME+'.md', 'checks/quantum_reference_clock_centering.py',
          'notes/uniform-dirac-reference.md', 'checks/uniform_dirac_reference.py',
          'notes/physical-clock-remainder.md', 'notes/spherical-anomaly-sources.md']
DPS = 60
TOL = '1e-35'
CELLS = 64


class R:
    """Closed rational interval; no binary floating operations in acceptance."""
    def __init__(self, lo, hi=None):
        if isinstance(lo, R):
            self.lo, self.hi = lo.lo, lo.hi
        else:
            self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(a, b):
        b = R(b)
        return R(a.lo+b.lo, a.hi+b.hi)
    __radd__ = __add__

    def __neg__(a):
        return R(-a.hi, -a.lo)

    def __sub__(a, b):
        return a+-R(b)

    def __rsub__(a, b):
        return R(b)+-a

    def __mul__(a, b):
        b = R(b)
        v = [a.lo*b.lo, a.lo*b.hi, a.hi*b.lo, a.hi*b.hi]
        return R(min(v), max(v))
    __rmul__ = __mul__

    def __truediv__(a, b):
        b = R(b)
        assert not b.lo <= 0 <= b.hi
        return a*R(1/b.hi, 1/b.lo)

    def __rtruediv__(a, b):
        return R(b)/a

    def __pow__(a, n):
        assert isinstance(n, int) and n >= 0
        out = R(1)
        for _ in range(n):
            out = out*a
        return out

    def compact(a, denominator=10**20):
        lo = (a.lo.numerator*denominator)//a.lo.denominator
        hi = -((-a.hi.numerator*denominator)//a.hi.denominator)
        return R(F(lo, denominator), F(hi, denominator))

    def record(a):
        return {'rational': [str(a.lo), str(a.hi)],
                'decimal_display_only': [float(a.lo), float(a.hi)]}


def exponential(x):
    x = F(x)
    assert x >= 0
    terms = max(180, math.ceil(2*x)+40)
    term = F(1)
    total = term
    for k in range(1, terms+1):
        term *= x/k
        total += term
    first = term*x/(terms+1)
    ratio = x/(terms+2)
    assert ratio < 1
    return R(total, total+first/(1-ratio)).compact()


def pi_bounds():
    def atan(x, n):
        s = sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n)), F(0))
        other = s+(-1)**n*x**(2*n+1)/F(2*n+1)
        return R(min(s, other), max(s, other))
    return (16*atan(F(1, 5), 48)-4*atan(F(1, 239), 16)).compact()


def constants():
    pi = pi_bounds()
    assert F(314, 100) < pi.lo < pi.hi < F(22, 7)
    q = 2000**4
    logq = R(F(30403, 1000), F(30404, 1000))
    assert exponential(logq.lo).hi < q < exponential(logq.hi).lo
    log100 = R(F('4.60517018'), F('4.60517019'))
    assert exponential(log100.lo).hi < 100 < exponential(log100.hi).lo
    H = sum((F(1, j) for j in range(1, 101)), F(0))
    gamma_precise = R(H)-log100-F(1, 200)+R(0, F(1, 120000))
    gamma = R(F(577, 1000), F(578, 1000))
    assert gamma.lo < gamma_precise.lo and gamma_precise.hi < gamma.hi
    return pi, q, logq, gamma


def flat_pressure(pi, logq, gamma):
    values = [R(1)]
    for j in range(1, CELLS+1):
        values.append((1+2/(exponential(F(CELLS, j))-1)).compact())
    # Monotonicity follows analytically from the coth derivative. Nearly
    # constant cells may have overlapping rational value enclosures.
    assert all(values[j].lo <= values[j+1].hi for j in range(CELLS))
    ir = R(sum((v.lo for v in values[:-1]), F(0))/CELLS,
           sum((v.hi for v in values[1:]), F(0))/CELLS)
    uv = R(-F(1, 720), -F(1, 720)+F(1, 60480))
    integral = (uv+ir)/(8*pi**2)
    value = -(integral-1/(8*pi**2)+(gamma-logq)/(48*pi**2))/2
    return value, {'J_ir': ir.record(), 'J_uv': uv.record(),
                   'pressure': value.record(), 'monotone_cells': CELLS}


def sphere_poly(t, h, pi):
    out = 1/(2*pi*t)
    for j in range(1, 4):
        for v in range(j+1):
            out += sp.bernoulli(2*j)*(-1)**(2*j-v)*t**(2*j-v-1)*h**v/(
                2*pi*sp.factorial(2*j-2*v)*sp.factorial(v))
    return sp.expand(out)


def sf(value):
    value = sp.Rational(value)
    return F(int(value.p), int(value.q))


def full_kernel_error(q):
    t, h, r = sp.symbols('t h r', positive=True)
    pi = sp.Symbol('p', positive=True)
    S3 = sphere_poly(t, h, pi)
    HP3 = 1/(2*pi*t)+r**2*t/(120*pi)-r**3*t**2/(252*pi)
    product = sp.expand(HP3*S3)
    expected = [1/(4*pi**2), -h/(24*pi**2), (15-3*h*h+3*r*r)/(720*pi**2)]
    assert all(sp.simplify(product.coeff(t, power)-c) == 0
               for power, c in zip([-2, -1, 0], expected))
    difference = sp.expand(product-S3.subs(h, 0)/(2*pi*t))
    difference -= sum(difference.coeff(t, j)*t**j for j in [-2, -1, 0])
    tau, rmax, hi = F(1, 4), F(31, 10000), F(1, q)
    uv_poly = F(0)
    for powers, coeff in sp.Poly(sp.expand(difference.subs(pi, 3)), t, h, r).terms():
        pt, ph, pr = powers
        assert pt >= 1
        uv_poly += abs(sf(coeff))*tau**pt*hi**ph*rmax**pr/pt
    # Absolute multivariate coefficients first: setting h=1 first could
    # cancel opposite signs and would not majorize the full interval.
    poly_bound = sum(abs(sf(c))*tau**p[0] for p, c in
                     sp.Poly(sp.expand(t*S3.subs(pi, 3)), t, h).terms())
    assert poly_bound < 1
    assert F(1, 6)+tau*tau/F(360) < 1

    def sphere_error(inverse_q):
        result = F(0)
        for j in range(4):
            for v in range(4-j):
                coefficient = F(math.factorial(7)*math.comb(3-j, v)*4**v*math.factorial(v),
                                math.factorial(j)*math.factorial(7-2*j)*69984)
                power = 5-j-v
                result += coefficient*inverse_q**(j+v)*tau**(power-1)/(power-1)
        return result

    uv_remainders = sphere_error(hi)+sphere_error(F(0))+rmax**4*tau**2/2880
    # These are strict upper bounds using pi>3, log4<2, |gamma-logq|<31.
    angular_tail = (F(1, 6)/tau**2+rmax**2/360)*(F(17, 6)/q+F(65, 3)/q**2)
    curvature_tail = rmax/48
    dp1, dp2 = F(1, 216*q), (rmax*rmax+hi*hi)/2160
    subtraction = 3*dp1+2*dp2
    finite_terms = dp1+31*dp2
    total = (uv_poly+uv_remainders+angular_tail+curvature_tail+subtraction+finite_terms)/2
    assert total < F(1, 10000)
    return total, {k: str(v) for k, v in dict(sphere_polynomial_majorant=poly_bound,
               UV_polynomial=uv_poly,
               UV_remainders=uv_remainders, angular_IR=angular_tail,
               curvature_IR=curvature_tail, subtraction=subtraction,
               finite_terms=finite_terms, pressure_error=total).items()}


def candidate_bounds(pi, q, pressure):
    k = 2000
    assert q == k**4
    c_lower, c_upper = F(6*k**3), F(123*k**3, 20)
    assert c_lower*c_lower < 1+12*pi.lo*k**6
    assert (6*k**3+1)**2 < 1+12*pi.lo*k**6  # xi=(C-1)/6 > k^3
    assert 1+12*pi.hi*k**6 < c_upper*c_upper
    assert c_upper/q < F(31, 10000)
    ratio = 2*pi*k**5/R(c_lower, c_upper)
    remainder = F(5, 10**43)
    mu_lower = F(160, 41*k**10)
    mu_upper = F(4, k**10)
    displacement = -(ratio*pressure)+R(-remainder/mu_lower, remainder/mu_lower)
    assert displacement.hi < -100000
    eps_lower = F(4, k**8)
    sphere_down = eps_lower*pressure.lo-remainder
    radius_down = sphere_down/2
    assert radius_down > F(2, 10**28)
    assert radius_down > F('5.46875e-29')
    return {'mu_bounds': R(mu_lower, mu_upper).record(),
            'epsilon_over_mu': ratio.record(), 'normalized_area_translation': displacement.record(),
            'sphere_squared_fractional_decrease_lower': str(sphere_down),
            'sphere_radius_fractional_decrease_lower': str(radius_down),
            'reference_remainder': str(remainder),
            'raw_K_100000_hypothesis_for_stationary_product': False}


def symbolic_centering():
    X, xs, G, a, L, rm, phi = sp.symbols('X xs G a L rm phi', positive=True)
    V = 2+6*X/L**2-2*rm**2/X
    V0 = V.subs(X, xs)
    Js = 2/a**2-sp.diff(V, X).subs(X, xs)
    ell = (-V0+Js*(X-xs))/(4*G)
    effective = V+4*G*ell
    assert sp.simplify(effective.subs(X, xs)) == 0
    assert sp.simplify(sp.diff(effective, X).subs(X, xs)-2/a**2) == 0
    assert sp.simplify(4*G*sp.diff(ell, X)-Js) == 0
    assert sp.simplify(sp.diff(ell, X)-Js/(4*G)) == 0
    U = effective.subs(X, xs+2*G*phi)/(4*G)-phi/a**2
    assert sp.factor(U+2*G*rm**2*phi**2/(xs**2*(xs+2*G*phi))) == 0
    C, k, pi = sp.symbols('C k pi', positive=True)
    eps = 4/k**8*(1+1/C)
    mu = 24/((C-1)*k**7)
    assert sp.factor((eps/mu-2*pi*k**5/C)*6*C*k-(C*C-1-12*pi*k**6)) == 0
    b, alpha, lam = sp.symbols('b alpha lam', positive=True)
    rate = sp.exp(alpha)*lam/a**2
    assert sp.simplify(b*lam/G-b*a**2*sp.exp(-alpha)/G*rate) == 0
    return {'potential_tangent': 'exact', 'joint_Ward_reference_subtraction': 'exact',
            'nonlinear_rational_remainder': 'exact', 'normalization_with_exterior_lapse': 'exact'}


def numerical_centering():
    rows = []
    for b, a, G, rm, L, phi, alpha in [
        ('1.4', '.8', '.02', '1.9', '2.1', '.2', '.03'),
        ('2.1', '.6', '.015', '3.4', '1.7', '-.4', '-.02'),
        ('1.1', '.4', '.01', '2.8', '3.2', '.05', '.1')]:
        b, a, G, rm, L, phi, alpha = map(mp.mpf, [b, a, G, rm, L, phi, alpha])
        x0, x = b*b, b*b+2*G*phi
        V = lambda xx: 2+6*xx/L**2-2*rm**2/xx
        Js = 2/a**2-mp.diff(V, x0)
        effective = V(x)-V(x0)+Js*(x-x0)
        direct = effective/(4*G)-phi/a**2
        rational = -2*G*rm**2*phi**2/(b**4*x)
        error = abs(direct-rational)
        assert error < mp.mpf(TOL)
        # Independent coefficient from finite near-horizon metric/area jets.
        width = mp.mpf('.013')
        dt_dU = mp.exp(alpha)*width/a**2
        area_leading = b*width/G
        phi_r = b*a*a*mp.exp(-alpha)/G
        norm_error = abs(area_leading-phi_r*dt_dU)
        assert norm_error < mp.mpf(TOL)
        rows.append({'potential_error': str(error), 'normalization_error': str(norm_error)})
    return rows


def hashes():
    return {p: hashlib.sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in INPUTS}


def main():
    mp.mp.dps = DPS
    began = datetime.datetime.now(datetime.timezone.utc).isoformat()
    source = hashes()
    log = ROOT/'evidence'/f'{NAME}.log'
    result = {'started_utc': began, 'story': 'S137.1', 'inputs_sha256': source,
              'hash_policy': 'SHA256 with CRLF-to-LF normalization',
              'execution': {'executable': sys.executable, 'argv': sys.orig_argv,
                            'cwd': str(Path.cwd()), 'command': subprocess.list2cmdline(sys.orig_argv)},
              'protocol': {'rational_cells': CELLS, 'numeric_digits': DPS, 'numeric_tolerance': TOL,
                           'pressure_target': ['.028', '.030'], 'raw_area_target': '<-100000'},
              'runtime': {'python': sys.version, 'platform': platform.platform(),
                          'sympy': sp.__version__, 'mpmath': mp.__version__}}
    def emit(message):
        with log.open('a', encoding='utf-8', newline='\n') as f:
            f.write(message+'\n')
        print(message, flush=True)
    emit('RUN '+began+'; exact rational acceptance; frozen protocol '+json.dumps(result['protocol']))
    try:
        pi, q, logq, gamma = constants()
        flat, flat_data = flat_pressure(pi, logq, gamma)
        err, err_data = full_kernel_error(q)
        pressure = flat+R(-err, err)
        assert F(28, 1000) < pressure.lo and pressure.hi < F(30, 1000), pressure.record()
        result.update(constants={'pi': pi.record(), 'logq': logq.record(), 'gamma': gamma.record()},
                      flat_comparison=flat_data, infinite_kernel_errors=err_data,
                      pressure=pressure.record(), candidate=candidate_bounds(pi, q, pressure),
                      symbolic=symbolic_centering(), numerical=numerical_centering())
        assert hashes() == source, 'source changed during validation'
        result['status'] = 'PASS_RATIONAL_REFERENCE_COMPARISON_AND_EXACT_CENTERING'
        result['claims'] = {'specified_reference_raw_jet_window_excluded': True,
                            'centered_dynamical_jet_window_excluded': False,
                            'matched_quantum_exterior_constructed': False,
                            'physical_residual_or_opening_certified': False, 'Lean_run': False}
        emit('Pressure interval '+json.dumps(pressure.record()))
        emit('Normalized area translation '+json.dumps(result['candidate']['normalized_area_translation']))
        emit('TERMINAL PASS exit 0; all frozen inputs match')
    except Exception:
        result['status'] = 'FAIL'
        result['traceback'] = traceback.format_exc()
        emit(result['traceback'])
        emit('TERMINAL FAIL exit 1; retained assertion, no target change')
        raise
    finally:
        with (ROOT/'evidence'/f'{NAME}.json').open('w', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, indent=2)
            f.write('\n')


if __name__ == '__main__':
    main()
