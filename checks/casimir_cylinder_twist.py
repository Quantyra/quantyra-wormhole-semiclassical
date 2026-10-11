"""Static spatial Casimir response on a magnetic product cylinder.

Exact angular identities, an independent Jacobi-basis check, Matsubara-loop
checks, and numerical evaluation of analytic all-level remainder bounds.
This does not certify a dynamical kernel or wormhole stability.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import expit, kv
import sympy as sp


def angular_data(q):
    weights = [6*q*(q-1)/(5*(q+1)*(q+3)), 12*q/(5*(q+2)*(q+3))]
    shifts = [(q-2)/2, (2*q+1)/2]
    lambdas = [q+1, 2*(q+2)]
    return weights, shifts, lambdas


def mixed_loop(mass, length, wave):
    """Exact spatial-vertex mixed LLL/massive finite-size loop coefficient."""
    def integrand(x):
        p = x/length
        energy = math.hypot(p, mass)
        den = (mass*mass+wave*wave)**2+4*wave*wave*p*p
        low = 4*mass*mass*p*expit(-x)/den
        high = 2*mass*mass*(mass*mass+wave*wave+2*p*p)*expit(-length*energy)/(den*energy)
        return (low-high)/(math.pi*length)
    return quad(integrand, 0, np.inf, epsabs=2e-14, epsrel=2e-12)


def mixed_error_bound(mass, length, wave):
    leading = math.pi/(3*mass**2*length**2)
    power = leading*(2*wave**2/mass**2+wave**4/mass**4
                     +14*math.pi**2*wave**2/(5*mass**4*length**2))
    exponential = 2/math.pi*(kv(2, mass*length)+wave**2/mass**2*kv(0, mass*length))
    return power+exponential


def geometric_moments(r):
    return [1/(1-r), 1/(1-r)**2, (1+r)/(1-r)**3,
            (1+4*r+r*r)/(1-r)**4, (1+11*r+11*r*r+r**3)/(1-r)**5]


def heavy_bound(q, radius, length, wave, aa, hh):
    """Bound the entire massive/massive bubble plus contact contribution.

    Uses M_n>=M_1+(n-1)/R, M_n<=(n+q/2)/R and a Gaussian
    upper bound on K_nu; the infinite polynomial-geometric sum is exact.
    Floating values evaluate the analytic bound, not certified intervals.
    """
    x = length*math.sqrt(q+1)/radius
    r = math.exp(-length/radius)
    moments = geometric_moments(r)
    poly = np.polynomial.Polynomial
    deg = poly([q, 2])
    mass_upper = poly([q/2, 1])/radius
    vertex_bound = poly([abs(aa)+3*q*abs(hh)/radius**2, 3*abs(hh)/radius**2])
    def summed(p):
        return sum(c*moments[i] for i, c in enumerate(p.coef))
    def envelope(nu):
        return math.sqrt(math.pi/(2*x))*math.exp(-x+nu*nu/(2*x))
    bubble = 4/(math.pi*abs(wave))*envelope(1)*summed(deg*vertex_bound**2*mass_upper)
    contact = 6*hh*hh/(5*math.pi*radius**2)*(envelope(0)+envelope(2))*summed(deg*mass_upper**2)
    return bubble+contact


def symbolic_checks():
    u, m, p, r, k, mass = sp.symbols('u m p r k M', real=True)
    q = sp.Symbol('q', positive=True, integer=True)
    checks = {}
    def zero(name, value):
        assert sp.factor(value) == 0, name
        checks[name] = 'exact_zero'
    H = lambda f: sp.expand(-(1-u*u)*sp.diff(f, u, 2)+((q+1)*u-2*m)*sp.diff(f, u))
    y = sp.legendre(2, u)
    vertex = sp.diff(y, u)*(m-q*u/2)+(1-u*u)*sp.diff(y, u, 2)/4
    zero('metric_vertex_as_angular_operator', vertex-3*y/2+H(y)/2)
    moments = {0: sp.S.One, 1: 2*m/(q+1)}
    for n in range(2, 5):
        moments[n] = sp.factor((2*m*moments[n-1]+(n-1)*moments[n-2])/(q+n))
    def mean(f):
        return sp.factor(sum(c*moments[n] for (n,), c in sp.Poly(sp.expand(f), u).terms()))
    centered = y-mean(y)
    lambdas = [q+1, 2*(q+2)]
    modes = [(lambdas[1]*centered-H(centered))/(lambdas[1]-lambdas[0]),
             (H(centered)-lambdas[0]*centered)/(lambdas[1]-lambdas[0])]
    index = sp.Symbol('a', integer=True, nonnegative=True)
    def total(f):
        return sp.factor(sum(c*sp.summation(((q-1)/2-index)**n, (index, 0, q-1))
                             for (n,), c in sp.Poly(f, m).terms()))
    weights, shifts, _ = angular_data(q)
    for n, mode in enumerate(modes):
        zero('angular_eigenmode_'+str(n+1), H(mode)-lambdas[n]*mode)
        zero('angular_weight_'+str(n+1), total(mean(mode*mode))-weights[n])
    zero('angular_orthogonality', mean(modes[0]*modes[1]))
    zero('contact_angular_weight', total(mean((1-u*u)*sp.diff(y, u)**2))-6*q/5)
    gtt = sp.factor(sum(weights[n]*shifts[n]**2/lambdas[n] for n in range(2)))
    deficit = 18*q*(2*q**3+3*q**2+q+3)/(5*(q+1)**2*(q+2)**2)
    zero('finite_flux_deficit', 6*q/5-4*gtt-deficit)
    shifted = sp.Poly(sp.expand(sp.together(deficit-2).as_numer_denom()[0].subs(q, p+3)), p)
    assert all(c > 0 for c in shifted.all_coeffs())
    checks['deficit_greater_than_two_for_q_ge_three'] = 'positive_shifted_polynomial'
    assert sp.limit(deficit, q, sp.oo) == sp.Rational(36, 5)

    den = (mass*mass+k*k)**2+4*k*k*p*p
    A = -2*k*(mass*mass+k*k+4*p*p)/den
    B = 4*mass*mass*p*p/den
    C = -A
    D = -2*mass*mass*(mass*mass+k*k+2*p*p)/den
    numerator = 2*(p*p-r*(r+k))
    zero('Matsubara_partial_fractions', numerator-(A*r+B)*((r+k)**2+p*p+mass*mass)
         -(C*(r+k)+D)*(r*r+p*p))

    # Independent matrix Schur expansion fixes the sign of the kinetic term.
    sx, sy, sz = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)
    v0, v1, p0, p1 = sp.symbols('v0 v1 p0 p1', real=True)
    vp, pp = sx*v0+sy*v1, sx*p0+sy*p1
    expected = 2*vp*(v0*p0+v1*p1)-(v0*v0+v1*v1)*pp
    for i, value in enumerate(vp*pp*vp-expected):
        zero('Clifford_Schur_entry_'+str(i), value)

    # An affine metric path needs a quadratic coframe correction.
    R, sn, v = sp.symbols('R sn v', nonzero=True, real=True)
    co0 = sp.diag(1, 1, R, R*sn)
    co1, co2 = sp.zeros(4), sp.zeros(4)
    co1[3, 1] = R*sn*v
    co2[1, 1] = -R*R*sn*sn*v*v/2
    assert co0.T*co2+co2.T*co0+co1.T*co1 == sp.zeros(4)
    fr2 = co0.inv()*co1*co0.inv()*co1*co0.inv()-co0.inv()*co2*co0.inv()
    zero('affine_metric_contact_sign', fr2[1, 1]-R*R*sn*sn*v*v/2)
    # l=1 plus the pulled-back monopole is a periodic coordinate change.
    # Its full angular vertex is a scalar in each m sector, so the
    # nonzero-frequency bubble vanishes; only the affine-path contact remains.
    h = sp.Symbol('h', real=True)
    pulled_back_charge = -q*h*u/(2*R**2)
    l1_metric_vertex = m-q*u/2
    zero('l1_coordinate_vertex_cancellation',
         pulled_back_charge-h*l1_metric_vertex/R**2+h*m/R**2)
    rr = sp.Rational(1, 64)
    assert (1+11*rr+11*rr**2+rr**3)/(1-rr)**5 < 2
    checks['uniform_geometric_tail_bound'] = 'exact_rational_check'
    return checks


def independent_mode_checks():
    u = sp.Symbol('u')
    cases = 0
    for q in [3, 4, 7, 12]:
        got = [sp.zeros(2), sp.zeros(2)]
        for a in range(q):
            alpha, beta = a, q-1-a
            m = sp.Rational(q-1, 2)-a
            y = sp.legendre(2, u)
            t = sp.diff(y, u)*(m-q*u/2)+(1-u*u)*sp.diff(y, u, 2)/4
            # Raw beta-weight integrals, independent of the H polynomial inverse.
            density = (1-u)**alpha*(1+u)**beta
            def integral(f):
                poly = sp.Poly(sp.expand(f*density), u)
                return sum(c*sp.Rational(2, n+1) for (n,), c in poly.terms() if n % 2 == 0)
            normalization = integral(1)
            for n in [1, 2]:
                jacobi = sp.jacobi(n, alpha, beta, u)
                norm = integral(jacobi*jacobi)*normalization
                xy, xt = integral(y*jacobi), integral(t*jacobi)
                got[n-1] += sp.Matrix([[xy*xy, xy*xt], [xy*xt, xt*xt]])/norm
                cases += 1
        weights, shifts, _ = angular_data(sp.Integer(q))
        for i in range(2):
            expected = weights[i]*sp.Matrix([[1, -shifts[i]], [-shifts[i], shifts[i]**2]])
            assert (got[i]-expected).applyfunc(sp.factor) == sp.zeros(2)
    return cases


def numerical_checks():
    residuals = []
    for length, mass, p, mode in [(2., .7, .4, 1), (3., 1.2, .9, 2), (5., 2., .3, 1)]:
        k = 2*math.pi*mode/length
        energy = math.hypot(mass, p)
        den = (mass*mass+k*k)**2+4*k*k*p*p
        B = 4*mass*mass*p*p/den
        D = -2*mass*mass*(mass*mass+k*k+2*p*p)/den
        n = np.arange(-20000, 20000)
        rr = 2*math.pi*(n+.5)/length
        f = 2*(p*p-rr*(rr+k))/((rr*rr+p*p)*((rr+k)**2+p*p+mass*mass))
        # Subtract the leading tail and restore its exact antiperiodic sum.
        discrete = np.sum(f+2/(rr*rr))/length-length/2
        continuous = B/(2*p)+D/(2*energy)
        expected = -B*expit(-length*p)/p-D*expit(-length*energy)/energy
        error = abs(discrete-continuous-expected)
        assert error < 2e-11
        residuals.append(error)
    error_bound_ratios = []
    # Include short circles as well as the theorem's long-circle regime.
    for mass, length, mode in [(.7, 2., 1), (1.2, 3., 2), (2., 5., 1),
                              (2., 16., 1), (3., 32., 2)]:
        wave = 2*math.pi*mode/length
        value, error = mixed_loop(mass, length, wave)
        deficit = math.pi/(3*mass**2*length**2)-value
        bound = mixed_error_bound(mass, length, wave)
        assert -error <= deficit <= bound+error
        error_bound_ratios.append(deficit/bound)
    samples = []
    for q in [3, 8, 16, 64]:
        for length in [16., 32.]:
            radius = 1.
            k = 2*math.pi/length
            weights, shifts, lambdas = angular_data(q)
            matrix = np.zeros((2, 2))
            quadrature_error = 0.
            for weight, shift, lam in zip(weights, shifts, lambdas):
                mass = math.sqrt(lam)
                value, error = mixed_loop(mass, length, k)
                assert value <= math.pi/(3*mass*mass*length*length)+error
                vec = np.array([1., shift])
                matrix += weight*value*np.outer(vec, vec)
                quadrature_error += weight*error*shift*shift
            matrix[1, 1] -= math.pi*(6*q/5)/(12*length*length)
            gtt = sum(weights[i]*shifts[i]**2/lambdas[i] for i in range(2))
            leading = math.pi*(4*gtt-6*q/5)/(12*length*length)
            bound = heavy_bound(q, radius, length, k, 0., 1.)
            assert matrix[1, 1]+bound+quadrature_error < 0
            assert bound < abs(leading)/128
            samples.append({'q': q, 'length_over_radius': length,
                            'mixed_plus_contact_hh': matrix[1, 1],
                            'leading_upper_bound_hh': leading,
                            'all_higher_sector_absolute_bound_hh': bound,
                            'quadrature_error_estimate_hh': quadrature_error,
                            'full_hh_coefficient_strictly_negative': True})
    baseline, error = quad(lambda x: math.log1p(math.exp(-x)), 0, np.inf, epsabs=1e-13)
    assert abs(-2*baseline/math.pi+math.pi/6) < 1e-12
    return {'independent_Matsubara_cases': len(residuals), 'max_Matsubara_residual': max(residuals),
            'mixed_error_bound_cases': len(error_bound_ratios),
            'max_mixed_error_to_bound_ratio': max(error_bound_ratios),
            'Casimir_energy_unit_circle': -2*baseline/math.pi, 'samples': samples}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = {'status': 'passed_static_cylinder_twist_response_checks',
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'sympy_version': sp.__version__, 'scipy_version': scipy.__version__,
              'symbolic_checks': symbolic_checks(), 'independent_Jacobi_mode_cases': independent_mode_checks(),
              'numerical_checks': numerical_checks(),
              'scope': 'Free Dirac finite-size static spatial response on R_t x S1 x S2, l=2, q>=3',
              'analytic_all_level_bound': 'See note; infinite massive sector bounded, not discarded at a numerical cutoff',
              'dynamical_kernel_complete': False, 'wormhole_stability_proved': False,
              'lean_certified': False}
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
