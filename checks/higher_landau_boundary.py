"""S072: exact RN exterior matching and conditional massive-channel bound.

The finite-cut theorem assumes a bounded inner retarded load. Candidate
thermal spectra are leading matched estimates, not a uniform error theorem.
"""
from pathlib import Path
import hashlib
import json
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]


def H(X, u):
    return X + 2*X*u + (1+3*X)*u*u


def action(X, u):
    C = 1+6*X
    delta = 1-u
    return mp.log((mp.sqrt(C*H(X, u))+C-(1+4*X)*delta)
                  / (delta*(mp.sqrt(C*X)+2*X)))/mp.sqrt(C)


def matching_checks():
    rows = []
    for X in map(mp.mpf, ['.1', '1', '100', '8.186613663405e9']):
        C = 1+6*X
        Lam = 2*C/(mp.sqrt(C*X)+2*X)
        derivative_errors = []
        quadrature_errors = []
        for u in map(mp.mpf, ['.1', '.5', '.9', '.99']):
            integrand = 1/((1-u)*mp.sqrt(H(X, u)))
            derivative_errors.append(abs(mp.diff(lambda v: action(X, v), u)/integrand-1))
            independent = mp.quad(lambda v: 1/((1-v)*mp.sqrt(H(X, v))), [0, u])
            quadrature_errors.append(abs(action(X, u)/independent-1))
        delta = mp.mpf('1e-20')
        asymptotic_error = abs(mp.sqrt(C)*action(X, 1-delta)+mp.log(delta)-mp.log(Lam))
        assert max(derivative_errors+quadrature_errors) < mp.mpf('1e-55')
        assert asymptotic_error < 2*delta
        rows.append(dict(X=float(X),max_relative_error=float(max(derivative_errors+quadrature_errors)),
                         near_horizon_constant_error=float(asymptotic_error)))
    return rows


def ir_spectral(M, T, omega):
    y = omega/(2*mp.pi*T)
    A = M+mp.mpf('.5')
    return (mp.pi*T)**(2*M)*abs(mp.gamma(A+1j*y))**2/mp.gamma(A)**2*mp.cosh(mp.pi*y)


def inner_checks():
    rows = []
    for M in map(mp.mpf, ['.7', '1.3', '3.2', '18.04813031688602']):
        # Independent exact Poincare-AdS2 Hankel solution and its flux.
        z = mp.mpf('1e-6')
        ratio = mp.hankel1(M-mp.mpf('.5'), z)/mp.hankel1(M+mp.mpf('.5'), z)
        flux = 2/(mp.pi*z*abs(mp.hankel1(M+mp.mpf('.5'), z))**2)
        leading = mp.pi/mp.gamma(M+mp.mpf('.5'))**2*(z/2)**(2*M)
        assert abs(mp.im(ratio)/flux-1) < mp.mpf('1e-55')
        assert abs(flux/leading-1) < mp.mpf('1e-8')
        errors = []
        for y in map(mp.mpf, ['.2', '2', '10']):
            T = mp.mpf('.03')
            exact = (1j*(mp.pi*T)**(2*M)*mp.gamma(mp.mpf('.5')-M)
                     /mp.gamma(mp.mpf('.5')+M)
                     *mp.gamma(mp.mpf('.5')+M-1j*y)
                     /mp.gamma(mp.mpf('.5')-M-1j*y))
            positive = ir_spectral(M, T, 2*mp.pi*T*y)
            errors.append(abs(mp.im(exact)/positive-1))
        assert max(errors) < mp.mpf('1e-55')
        rows.append(dict(M=float(M),bessel_asymptotic_relative_error=float(abs(flux/leading-1)),
                         thermal_identity_relative_error=float(max(errors))))
    for y in [0, 1, 20]:
        assert abs(ir_spectral(mp.mpf(0), mp.mpf('.1'), 2*mp.pi*mp.mpf('.1')*y)-1) < mp.mpf('1e-55')
    for M in [mp.mpf('.5'), mp.mpf('1.5')]:
        center = ir_spectral(M, mp.mpf('.1'), mp.mpf('.2'))
        assert center > 0
        for sign in [-1, 1]:
            assert abs(ir_spectral(M+sign*mp.mpf('1e-8'), mp.mpf('.1'), mp.mpf('.2'))/center-1) < mp.mpf('1e-6')
    return rows


def outer_checks():
    rows = []
    for X, lam, delta, load in [(1., 3., .1, .04+.08j), (5., 8., .05, -.07+.04j),
                                (50., 12., .2, .02+.02j)]:
        R = .1
        uc = 1-delta
        Wmin = lam*min(np.sqrt(X), delta*np.sqrt(H(X, uc)))
        w = .1*Wmin
        assert w*(1+R*R) < 2*Wmin*R
        def rhs(u, state):
            a, ell, integ = state
            b = np.exp(ell)
            den = (1-u)**2*H(X, u)
            W = lam*(1-u)*np.sqrt(H(X, u))
            return [(2*W*a-w*(1+a*a-b*b))/den,
                    (2*W-2*w*a)/den, -a/den]
        sol = solve_ivp(rhs, (uc, 0), [load.real, np.log(load.imag), 0.],
                        method='DOP853', rtol=2e-12, atol=2e-13, max_step=.002)
        assert sol.success
        a, ell, integ = sol.y[:, -1]
        radius = np.hypot(sol.y[0], np.exp(sol.y[1]))
        assert max(radius) <= R
        S = float(lam*action(mp.mpf(X), mp.mpf(uc)))
        identity_error = abs(ell-np.log(load.imag)+2*S-2*w*integ)
        assert identity_error < 2e-10
        tau_upper = (1/delta-1)/X
        b = np.exp(ell)
        rho = 2*b/((1+a)**2+b*b)
        bound = 2*R/(1-R)**2*np.exp(-2*S+2*w*R*tau_upper)
        assert 0 < rho <= bound
        rows.append(dict(X=X,lambda_n=lam,delta=delta,omega_b=w,
                         largest_cut_ratio=float(max(radius)),
                         log_flux_identity_error=float(identity_error),rho=float(rho),bound=float(bound)))
    return rows


def tower_bracket(a, q):
    return 1+2/q+2*(1/a+1/a**2)+4/q*(1/a+3/a**2+6/a**3+6/a**4)


def tower_checks():
    rows = []
    for q in [1, 10, 1000]:
        for a in [4., 10., 30.]:
            n = np.arange(1, 100001, dtype=float)
            # Normalize by exp(-a) to keep the comparison well scaled.
            partial = np.sum((1+2*n/q)*np.exp(-a*(np.sqrt(n)-1)))
            bound = tower_bracket(a, q)
            assert partial <= bound
            rows.append(dict(q=q,a=a,normalized_partial_sum=float(partial),analytic_bound=bound))
    return rows


def candidate(ki):
    k = mp.mpf(ki); q = k**4
    X = (mp.sqrt(1+12*mp.pi*k**6)-1)/6
    C = 1+6*X
    Tb = (2*q/k)*mp.mpf('.5')*C/(2*mp.pi*k**10*X)
    Omega = mp.sqrt(q*q/(C*(1+3*X)))
    Gamma = 4*mp.pi/k*Omega**2
    Lam = 2*C/(mp.sqrt(C*X)+2*X)
    M = mp.sqrt((q+1)/C)
    delta = k**-4; R = mp.mpf('.1'); uc = 1-delta
    S = mp.sqrt(q+1)*action(X, uc)
    wmax = 2*mp.pi*Tb*(2*Gamma)
    Wmin = mp.sqrt(q+1)*min(mp.sqrt(X), delta*mp.sqrt(H(X, uc)))
    correction = 2*wmax*R*(1/delta-1)/X
    assert wmax*(1+R*R) < 2*Wmin*R
    a = 2*S
    log_bound = mp.log(2*R/(1-R)**2)-a+correction
    tower_log_bound = log_bound+mp.log(2*tower_bracket(a, q))
    spectra = []
    for E in [Omega, 2*Gamma]:
        A = M+mp.mpf('.5')
        logrho = (mp.log(2)+2*M*mp.log(mp.pi*Tb/(C*Lam))
                  +mp.log(mp.cosh(mp.pi*E))+2*mp.re(mp.loggamma(A+1j*E))
                  -2*mp.loggamma(A))
        spectra.append(dict(global_frequency=float(E),log10_leading_rho=float(logrho/mp.log(10))))
    return dict(k=ki,q=int(q),X=float(X),Tb=float(Tb),M1=float(M),Lambda=float(Lam),
                cut_delta=float(delta),assumed_inner_radius=float(R),
                band_max_global_frequency=float(2*Gamma),
                omega_b_over_Wmin_b=float(wmax/Wmin),
                logarithmic_outer_frequency_correction=float(correction),
                log10_conditional_first_level_bound=float(log_bound/mp.log(10)),
                log10_conditional_tower_bound_per_LLL_channel=float(tower_log_bound/mp.log(10)),
                leading_thermal_spectral_estimates=spectra)


def main():
    mp.mp.dps = 75
    result = dict(date='2026-10-13',status='PASS',python=platform.python_version(),
                  numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  matching=matching_checks(),inner_normalization=inner_checks(),
                  finite_frequency_exterior=outer_checks(),tower_sum=tower_checks(),
                  candidate=[candidate(k) for k in [2000,4000,8000]],
                  limitations=[
                      'The exact finite-cut theorem assumes |f/g| <= 0.1 on the inner cut for every included mode and frequency',
                      'No uniform matching remainder establishes that condition for the evolving physical throat',
                      'Thermal spectra are leading initial-RN matching estimates, not certified bounds',
                      'Boundary spectral suppression does not bound bulk gauge polarization, all-frequency switching, or complete stress'])
    (ROOT/'evidence/higher-landau-boundary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: exact radial integral, Bessel/thermal normalization, finite-frequency cut propagation and tower bound.')
    print(json.dumps(result['candidate'][0],indent=2))


if __name__ == '__main__':
    main()
