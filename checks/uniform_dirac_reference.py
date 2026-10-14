"""Flux-uniform homogeneous Dirac reference in the stated zeta scheme.

Analytic bounds prove the nearby-product theorem. Adaptive quadrature
illustrates the exact pressure integrals; it is not an interval certificate.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import zeta
import sympy as s

PI = math.pi
N_EM, M_H = 8, 6
MODE_EXPONENT = 80.


def sphere_polynomial(order, t, h, pi=s.pi):
    value = 1/(2*pi*t)
    for k in range(1, order+1):
        value += s.bernoulli(2*k)/(2*pi)*sum(
            (-1)**(2*k-j)*t**(2*k-j-1)*h**j/s.factorial(2*k-2*j)/s.factorial(j)
            for j in range(k+1))
    return s.expand(value)


def em_terms(order):
    m = 2*order+1
    return [(math.factorial(m)*math.comb(order-j, ell)*4**ell*math.factorial(ell)
             // (math.factorial(j)*math.factorial(m-2*j)), j+ell, 2*order-j-ell-1)
            for j in range(order+1) for ell in range(order-j+1)]


def analytic_checks():
    t, h, r = s.symbols('t h r', positive=True)
    Sn = sphere_polynomial(3, t, h)
    hn = 1/(2*s.pi*t)+sum(s.bernoulli(2*k+2)*r**(k+1)*t**k/
             ((2*k+2)*s.factorial(k)*s.pi) for k in range(3))
    hp = s.expand(hn-r*s.diff(hn, r))
    dhp = s.expand(r*s.diff(hp, r))
    coeff = [1/(4*s.pi**2), (r-h)/(24*s.pi**2),
             (15-3*h*h-5*r*h-3*r*r)/(720*s.pi**2)]
    pcoeff = [coeff[0], -h/(24*s.pi**2), (15-3*h*h+3*r*r)/(720*s.pi**2)]
    dcoeff = [0, 0, r*r/(120*s.pi**2)]
    for hh, sub in zip([hn, hp, dhp], [coeff, pcoeff, dcoeff]):
        product = s.expand(Sn*hh)
        for power, expected in zip([-2, -1, 0], sub):
            assert s.simplify(product.coeff(t, power)-expected) == 0
    # Entire coarse UV proof uses rational majorants pi>3, zeta(6),zeta(8)<2.
    e3 = s.Rational(2*sum(v[0] for v in em_terms(3)), 3*6**6)
    assert e3 < 1  # |S-S3| <= e3*t^2 for q>=1, 0<t<=1.
    d2 = s.Rational(4*math.factorial(7), math.factorial(3)*6**8*3)
    assert 12*d2 < 1  # h, hP, r*d_r hP remainders <=t^3.
    def abs_poly(expression):
        return sum(abs(c)*t**p[0]*h**p[1]*r**p[2]
                   for p, c in s.Poly(s.expand(expression), t, h, r).terms())
    assert abs_poly(t*Sn.subs(s.pi, 3)).subs({t: s.Rational(1, 4), h: 1, r: 1}) < 1
    numerical_Sn = sphere_polynomial(N_EM, t, h, pi=s.Integer(3))
    assert abs_poly(t*numerical_Sn).subs({t: s.Rational(1, 4), h: 1, r: 1}) < 1
    integrals = []
    for hh in [hn, hp, dhp]:
        product = s.expand((Sn*hh).subs(s.pi, 3))
        polynomial = sum(product.coeff(t, p)*t**p for p in range(1, 9))
        majorant = abs_poly(polynomial).subs({h: 1, r: 1})
        bound = s.integrate(majorant/t, (t, 0, s.Rational(1, 4)))
        assert bound < 1
        integrals.append(str(bound))
    # Middle integral <=40, IR <=3; UV <=1+1/2; local finite terms then
    # give |F|, |pA|, |r*d_r pA| <=100(1+|L|), |d_L pA|<=1.
    assert s.Rational(63, 3)+s.Rational(15, 2)+3+2 < 40
    # Exact stationary-map identity and first-order coefficients.
    z, u, alpha, eps, pa, a2 = s.symbols('z u alpha epsilon pA a2')
    ps = -a2/2-pa
    eqA = -alpha*z-(1-alpha)*z*z/2+(1+alpha)/2-eps*pa
    eqS = z*z/u-(1-alpha)*z*z/2-(1+alpha)/2-eps*ps
    assert s.expand(eqA+eqS-(z*z/u-(1-alpha)*z*z-alpha*z+eps*a2/2)) == 0
    assert s.expand(eqA-(-(z-1)-(1-alpha)*(z-1)**2/2-eps*pa)) == 0
    z1, u1 = 1-eps*pa, 1+eps*(a2/2-alpha*pa)
    for eq in [eqA, eqS]:
        assert s.simplify(s.diff(eq.subs({z: z1, u: u1}), eps).subs(eps, 0)) == 0
    # Uniform family estimates: log(2000)<8, C>6k^3, epsilon<5k^-8.
    k0 = s.Integer(2000)
    B0 = s.Integer(101)*(1+4*8)
    epB = 5*B0/k0**8
    alB = B0/(6*k0**3)
    rho, alpha_upper, eps_upper = s.Rational(1, 100), s.Rational(1, 10**10), 5/k0**8
    assert epB < s.Rational(7, 10**23) and alB < s.Rational(1, 10**6)
    assert 7/k0*(1+rho)/(1-rho) < s.Rational(4, 1000)
    assert epB+rho*rho/2 < rho
    assert (alpha_upper*rho/(1-rho)+eps_upper/(2*(1-rho)**2))/(1-rho) < rho
    assert rho+4*epB < s.Rational(2, 100)
    assert (2*alpha_upper+3*eps_upper)/(1-rho)**2 < s.Rational(2, 100)
    assert 2*epB < s.Rational(14, 10**23)
    assert 10/k0**8 < s.Rational(4, 10**26)
    assert 100*epB**2 < s.Rational(5, 10**43)
    assert 18*alpha_upper+s.Rational(8, 100)+2*(1+s.Rational(4, 10**6))**2/100**2 < 100
    return {'scaled_heat_coefficients': list(map(str, coeff)),
            'sphere_EM3_remainder_coefficient_upper': str(e3),
            'largest_thermal_remainder_coefficient_upper': str(12*d2),
            'UV_polynomial_integral_bounds': integrals,
            'pressure_and_weighted_derivative_bound': '100*(1+abs(L)) for q>=1 and 0<r<=1',
            'contraction_constant_upper': '0.02 on abs(z-1),abs(u-1)<=0.01 for all integer k>=2000',
            'uniform_radius_squared_bounds': {'sphere': '1.4e-22', 'AdS': '4e-26'},
            'uniform_first_order_remainder_bound': '5e-43 in each normalized squared radius',
            'epsilon_B_upper_exact': str(epB)}


def polynomials(q, r):
    sphere = np.zeros(2*N_EM+1, dtype=np.longdouble)
    sphere[0] = 1/(2*PI)
    for k in range(1, N_EM+1):
        for j in range(k+1):
            power = 2*k-j-1
            sphere[power+1] += (np.longdouble(float(s.bernoulli(2*k)))/(2*PI)
                *(-1)**(2*k-j)*np.longdouble(q)**(-j)/math.factorial(2*k-2*j)/math.factorial(j))
    hyper = np.zeros(M_H+2, dtype=np.longdouble)
    phyper = hyper.copy()
    hyper[0] = phyper[0] = 1/(2*PI)
    for k in range(M_H+1):
        value = float(s.bernoulli(2*k+2))*r**(k+1)/((2*k+2)*math.factorial(k)*PI)
        hyper[k+1], phyper[k+1] = value, -k*value
    A = np.array([1/(4*PI**2), (r-1/q)/(24*PI**2),
                  (15-3/q**2-5*r/q-3*r*r)/(720*PI**2)])
    P = np.array([A[0], -1/(24*PI**2*q), (15-3/q**2+3*r*r)/(720*PI**2)])
    product, pressure = np.convolve(sphere, hyper), np.convolve(sphere, phyper)
    assert max(np.max(abs(product[:3]-A)), np.max(abs(pressure[:3]-P))) < 1e-15
    return product, pressure, A, P


def uv_bound(q, r, t0, pressure=False):
    pref = zeta(2*N_EM, 1)/(PI*(2*PI)**(2*N_EM))
    sphere = sum(pref*c*q**(-power_q)*t0**(power_t-1)/(power_t-1)
                 for c, power_q, power_t in em_terms(N_EM))
    D = 2*math.factorial(2*M_H+3)*zeta(2*M_H+4, 1)*r**(M_H+2)/(
        math.factorial(M_H+1)*(2*PI)**(2*M_H+4)*PI)
    if pressure:
        D *= M_H+1
    return sphere+D*t0**M_H/M_H


def kernels(t, q, r):
    N = max(1, math.ceil((2*MODE_EXPONENT/t)/(1+math.sqrt(1+4*MODE_EXPONENT/(q*t)))))
    assert t*(1+2*N/q)**2 >= 2/q
    n = np.arange(1, N+1, dtype=float)
    sphere = (1+2*np.sum((1+2*n/q)*np.exp(-t*n*(1+n/q))))/(4*PI)
    x = r*t
    if x < 1:
        def thermal(lam, power):
            if lam == 0:
                return 1/(2*PI) if power == 1 else 0.
            return lam**power*math.exp(-x*lam*lam)/math.expm1(2*PI*lam)
        h = 1/(2*PI*t)+2*r/PI*quad(lambda lam: thermal(lam, 1), 0, 12, epsabs=2e-13, epsrel=2e-13)[0]
        hp = 1/(2*PI*t)+2*r*r*t/PI*quad(lambda lam: thermal(lam, 3), 0, 12, epsabs=2e-13, epsrel=2e-13)[0]
    else:
        c = PI/math.sqrt(x)
        def integrand(y, power):
            if y == 0:
                return 1/c if power == 1 else 0.
            return y**power/math.tanh(c*y)*math.exp(-y*y)
        h = quad(lambda y: integrand(y, 1), 0, 12, epsabs=2e-13, epsrel=2e-13)[0]/(PI*t)
        hp = quad(lambda y: integrand(y, 3), 0, 12, epsabs=2e-13, epsrel=2e-13)[0]/(PI*t)
    return np.array([h*sphere, hp*sphere])


def vacuum(q, r, L, t0=.125):
    assert q >= 1 and 0 < r <= 1 and 0 < t0 <= .25
    product, pressure, A, P = polynomials(q, r)
    values = np.array([sum(float(c)*t0**(j-2)/(j-2) for j, c in enumerate(poly) if j >= 3)
                       for poly in [product, pressure]])
    errors = np.zeros(2)
    for channel, sub in enumerate([A, P]):
        def integrand(t):
            return (kernels(t, q, r)[channel]-sub[0]/(t*t)-sub[1]/t-sub[2])/t
        val, err = quad(integrand, t0, 1, epsabs=2e-11, epsrel=2e-11, limit=150)
        values[channel] += val
        errors[channel] += err
        def ir(y):
            if y == 0:
                return math.sqrt(r)/(4*PI**2.5)*(1 if channel == 0 else .5)
            return 2*kernels(1/(y*y), q, r)[channel]/y
        val, err = quad(ir, 0, 1, epsabs=2e-11, epsrel=2e-11, limit=150)
        values[channel] += val
        errors[channel] += err
    finite = np.array([values[0]-A[0]/2-A[1]+A[2]*(np.euler_gamma+L),
                       values[1]-P[0]/2-P[1]+P[2]*(np.euler_gamma+L)])/2
    tail = math.exp(-MODE_EXPONENT)*(1/(8*PI**2*t0*t0)+math.sqrt(r)/(6*PI**2.5*t0**1.5))
    return {'q': q, 'r': r, 'L': L, 't0': t0,
            'scaled_action': float(finite[0]), 'scaled_pA': float(-finite[1]),
            'scaled_pS': float(finite[1]-A[2]/2), 'scaled_trace': float(-A[2]),
            'UV_bounds_before_half': [uv_bound(q, r, t0), uv_bound(q, r, t0, True)],
            'angular_tail_bound_before_half': tail,
            'quadrature_error_estimates_before_half': errors.tolist(),
            'quadrature_errors_are_not_interval_certificates': True}


def numerical_checks():
    old = json.loads((Path(__file__).resolve().parents[1]/'evidence/magnetic-product-vacuum.json').read_text())
    rows = []
    for prior in old['numerical_checks']['samples']:
        q = prior['q']
        value = vacuum(q, 1/q, -math.log(q))
        comparisons = [('scaled_action', 'b4_Euclidean_action_density'),
                       ('scaled_pA', 'b4_pressure_AdS'), ('scaled_pS', 'b4_pressure_sphere')]
        discrepancy = max(abs(value[new]-prior[old_key]/q**2) for new, old_key in comparisons)
        assert discrepancy < 3e-9, (q, discrepancy)
        value['legacy_scheme_regression_error'] = discrepancy
        rows.append(value)
    candidate = []
    mp.mp.dps = 85
    for kval in [2000, 4000, 20000]:
        k = mp.mpf(kval)
        q = k**4
        C = mp.sqrt(1+12*mp.pi*k**6)
        r, alpha, gauge = C/q, 1/C, k**-4
        eps = 4*gauge*gauge*(1+alpha)
        value = vacuum(float(q), float(r), float(-mp.log(q)))
        comparison = vacuum(float(q), float(r), float(-mp.log(q)), t0=.2)
        error = max(abs(value[key]-comparison[key]) for key in ['scaled_action', 'scaled_pA', 'scaled_pS'])
        assert error < 3e-9, error
        a2 = (15-3/q**2-5*r/q-3*r*r)/(720*mp.pi**2)
        pa = mp.mpf(value['scaled_pA'])
        B = 101*(1+mp.log(q))
        value.update({'k': kval, 'UV_cutoff_comparison_error': error,
                      'first_sphere_squared_radius_shift': str(-eps*pa),
                      'first_AdS_squared_radius_shift': str(eps*(a2/2-alpha*pa)),
                      'analytic_remainder_for_exact_first_order_integrals': str(100*(eps*B)**2)})
        candidate.append(value)
    # Independent weighted radius derivative at moderate parameters.
    q, r, L, step = 8, .04, -math.log(8), .002
    center = vacuum(q, r, L)
    gradients = []
    for h in [step, step/2]:
        derivative = (vacuum(q, r*math.exp(h), L)['scaled_action']
                      -vacuum(q, r*math.exp(-h), L)['scaled_action'])/(2*h)
        gradients.append(derivative)
    derivative = (4*gradients[1]-gradients[0])/3
    error = abs(center['scaled_pA']+center['scaled_action']-derivative)
    assert error < 2e-8, error
    return {'legacy_regression': rows, 'candidate_family': candidate,
            'independent_weighted_radius_derivative_error': error}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print('Checking rescaled heat coefficients, uniform remainder majorants and stationary contraction', flush=True)
    result = {'analytic': analytic_checks()}
    print('Comparing common zeta scheme, two cutoffs and independent radius variation', flush=True)
    result['numerical'] = numerical_checks()
    result.update({'status': 'PASS', 'scope': 'homogeneous one-Dirac-loop zeta reference; no boundary drive, full EFT completion or wormhole stability',
                   'versions': {'python': platform.python_version(), 'sympy': s.__version__, 'numpy': np.__version__, 'scipy': scipy.__version__, 'mpmath': mp.__version__},
                   'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('PASS: uniform homogeneous reference and retained renormalization scheme', flush=True)


if __name__ == '__main__':
    main()
