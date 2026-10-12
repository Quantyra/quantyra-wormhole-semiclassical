"""Global determinant, current contact and gauge reweighting checks."""
import os
for _key in ['OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS']:
    os.environ[_key] = '1'
import hashlib
import json
import platform
import runpy
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.special import erfc

ROOT = Path(__file__).resolve().parents[1]
OLD = runpy.run_path(str(ROOT/'checks/driven_fermion_sea.py'))
T = np.pi
RNG = np.random.default_rng(20261012)


def axes(q):
    m = np.arange(q)-(q-1)/2
    jp = np.diag(np.sqrt(np.arange(1, q)*(q-np.arange(1, q))), -1)
    return np.diag(m), (jp+jp.T)/2


def transfer_log(U, lengths, beta, cutoff=60):
    value = 0.
    for n in range(cutoff):
        r = np.exp(-(2*n+1)*np.pi*lengths/beta)
        value += 4*(np.linalg.slogdet(np.eye(len(r))+r[:, None]*U)[1]-np.sum(np.log1p(r)))
    # Each difference is at most 2q r/(1-r) in magnitude; geometric Matsubara tail.
    rnext = np.exp(-(2*cutoff+1)*np.pi*min(lengths)/beta)
    tail = 8*len(lengths)*rnext/((1-rnext)*(1-np.exp(-2*np.pi*min(lengths)/beta)))
    return float(value), float(tail)


def unitary_checks():
    rows = []
    for q in [1, 2, 3, 7, 16]:
        for _ in range(8):
            Z = RNG.normal(size=(q, q))+1j*RNG.normal(size=(q, q))
            U = expm(1j*(Z+Z.conj().T))
            r = np.exp(RNG.uniform(-12, -.001, q))
            actual = np.linalg.slogdet(np.eye(q)+r[:, None]*U)[1]
            bound = np.sum(np.log1p(r))
            assert actual <= bound+2e-13
            rows.append(dict(q=q, log_ratio=float(actual-bound)))
    return rows


def static_checks():
    rows = []
    for q in [2, 3]:
        lengths = T+np.linspace(.2, 2.5, q)
        beta = 1.7
        Z, X = axes(q)
        for a, b in [(.3, .2), (2., -1.), (6., 4.)]:
            def source(s):
                f = np.cos(np.pi*s/T)**2
                return f*(a*Z+b*np.sin(2*np.pi*s/T)*X)
            def rhs(s, u):
                return (-1j*source(s)@u.reshape(q, q)).ravel()
            sol = solve_ivp(rhs, [-T/2, T/2], np.eye(q, dtype=complex).ravel(),
                            method='DOP853', rtol=2e-12, atol=2e-13)
            assert sol.success
            U = sol.y[:, -1].reshape(q, q)
            assert np.linalg.norm(U.conj().T@U-np.eye(q)) < 2e-11
            exact, tail = transfer_log(U, lengths, beta)
            counter = beta/(2*np.pi)*quad(lambda s: np.trace(source(s)@source(s)).real,
                                         -T/2, T/2, epsabs=1e-11)[0]
            values = []
            for N in [4, 12, 32]:
                p, AA = OLD['matrices'](q, N, lengths)
                eig = np.linalg.eigvalsh(np.diag(p)+a*AA[0]+b*AA[1])
                raw = 2*np.sum(np.logaddexp(beta*eig/2, -beta*eig/2)-
                               np.logaddexp(beta*p/2, -beta*p/2))
                value = raw-counter
                assert value <= 2e-10
                values.append(dict(cutoff=N, covariant_log_ratio=float(value),
                                   uncorrected_log_ratio=float(raw), error=float(abs(value-exact))))
            assert values[-1]['error'] < 6e-8
            assert values[-1]['error'] < values[0]['error']/100
            assert exact < 0 and tail < 1e-100
            rows.append(dict(q=q, amplitudes=[a, b], transfer_log_ratio=exact, transfer_tail_bound=tail,
                             contact=float(counter), spectral_refinements=values))
    return rows


def time_functions(t, beta, scale):
    x = 2*np.pi*t/beta
    return scale*np.array([np.cos(x), .6*np.sin(2*x), .5*np.sin(x), .4*np.cos(2*x)])


def dynamic_case(q, N, scale, steps, beta=.7):
    lengths = T+np.linspace(.2, 2.5, q)
    p, AA = OLD['matrices'](q, N, lengths)
    tq = q*(q*q-1)/12
    dt = beta/steps
    U = np.eye(len(p), dtype=complex)
    trace_h, contact, holder = 0., 0., 0.
    reference = 2*np.sum(np.logaddexp(-beta*p, 0))
    for n in range(steps):
        a, b, c, d = time_functions((n+.5)*dt, beta, scale)
        h = np.diag(p)+a*AA[0]+b*AA[1]
        charge = c*AA[1]+d*AA[0]
        e, W = np.linalg.eigh(h)
        half = (W*np.exp(-dt*e/2))@W.conj().T
        U = half@expm(-1j*dt*charge)@half@U
        scalar = tq*(a*a*3*T/8+b*b*5*T/32)/(2*np.pi)
        contact += dt*scalar
        trace_h += dt*np.trace(h).real
        static = 2*np.sum(np.logaddexp(-beta*e, 0))+beta*np.sum(e)-reference-beta*scalar
        holder += dt/beta*static
    logdelta = 2*np.linalg.slogdet(np.eye(len(p))+U)[1]+trace_h-reference-contact
    assert logdelta <= holder+2e-9 and holder < 0
    return float(logdelta), float(holder), p, AA, float(contact)


def dynamic_checks():
    rows = []
    for q in [2, 3]:
        for scale in [.7, 3.]:
            for N in [3, 6]:
                ref = []
                for steps in [32, 64, 128]:
                    value, holder, p, AA, contact = dynamic_case(q, N, scale, steps)
                    ref.append(dict(steps=steps, log_ratio=value, holder_bound=holder))
                beta = .7
                def rhs(t, u):
                    a, b, c, d = time_functions(t, beta, scale)
                    h = np.diag(p)+a*AA[0]+b*AA[1]
                    charge = c*AA[1]+d*AA[0]
                    return (-(h+1j*charge)@u.reshape(len(p), len(p))).ravel()
                sol = solve_ivp(rhs, [0, beta], np.eye(len(p), dtype=complex).ravel(),
                                method='DOP853', rtol=2e-11, atol=2e-12)
                assert sol.success
                U = sol.y[:, -1].reshape(len(p), len(p))
                ode = 2*np.linalg.slogdet(np.eye(len(p))+U)[1]-2*np.sum(np.logaddexp(-beta*p, 0))-contact
                errors = [abs(r['log_ratio']-ode) for r in ref]
                assert errors[-1] < errors[0]/8 and errors[-1] < .001
                rows.append(dict(q=q, cutoff=N, scale=scale, time_refinements=ref,
                                 adaptive_ODE_log_ratio=float(ode), ODE_errors=list(map(float, errors))))
    return rows


def gaussian_checks():
    rows, amplitude_rows = [], []
    beta = 2.3
    for q in [2, 5, 12]:
        lengths = T+np.linspace(.2, 2.5, q)
        E = axes(q)[1]/np.sqrt(q*(q*q-1)/12)
        F = T/2
        rcoef = np.sqrt(2)*beta/(np.pi*np.sqrt(q))*beta*5*T/16
        quadratic = 0.
        for n in range(60):
            r = np.exp(-(2*n+1)*np.pi*lengths/beta)
            k = r/(1+r)
            quadratic += 2*F*F*np.sum(k[:, None]*(1-k[None, :])*abs(E)**2)
        def logweight(a):
            return transfer_log(expm(-1j*a*F*E), lengths, beta)[0]
        for a in [0., .5, 5., 20., 100.]:
            value = logweight(a)
            assert value <= 2e-13
            amplitude_rows.append(dict(q=q, amplitude=a, log_ratio=value,
                                       nonlinear_log_weight=float(value+quadratic*a*a)))
        for variance in [.002, .02, .05]:
            kappa = quadratic*variance
            rho = rcoef*2*np.sqrt(2/np.pi)*variance**1.5
            threshold = (np.sqrt(rho)/rcoef)**(1/3)
            probbad = erfc(threshold/np.sqrt(2*variance))
            eta = kappa+np.sqrt(rho)*(1+1/(3*(1-rho**(1/6))))
            z, w = np.polynomial.hermite.hermgauss(64)
            w = w/np.sqrt(np.pi)
            a = np.sqrt(2*variance)*z
            determinant = np.array([np.exp(logweight(x)) for x in a])
            normalizer = np.dot(w, determinant)
            tv = np.dot(w, abs(determinant/normalizer-1))/2
            normalizer_quad = quad(lambda x: np.exp(logweight(np.sqrt(variance)*x)-x*x/2)/np.sqrt(2*np.pi),
                                   -10, 10, epsabs=2e-11)[0]
            assert abs(normalizer-normalizer_quad) < 2e-10
            assert 0 <= 1-normalizer <= eta and tv <= 1-normalizer+1e-13
            assert probbad <= np.sqrt(rho)
            zg = (1+2*quadratic*variance)**(-.5)
            bad_integral = 2*quad(lambda x: np.exp(logweight(np.sqrt(variance)*x)-x*x/2)/np.sqrt(2*np.pi),
                                  threshold/np.sqrt(variance), np.inf, epsabs=2e-11)[0]
            weighted_tail = bad_integral/zg
            assert weighted_tail <= np.exp(kappa)*np.sqrt(rho)+1e-12
            rows.append(dict(q=q, variance=variance, kappa=float(kappa), rho=float(rho), eta=float(eta),
                             normalizer=float(normalizer), total_variation=float(tv),
                             normalizer_quadrature_error=float(abs(normalizer-normalizer_quad)),
                             bare_bad_probability=float(probbad), weighted_bad_integral=float(weighted_tail),
                             weighted_bad_bound=float(np.exp(kappa)*np.sqrt(rho))))
    return dict(gaussian_cases=rows, large_amplitude_cases=amplitude_rows,
                scope='One real Gaussian spatial rotation mode with exact finite-temperature determinant; no physical covariance claim.')


def global_det3_checks():
    rows = []
    beta = 2.3
    for q in [2, 3]:
        lengths = T+np.linspace(.2, 2.5, q)
        E = axes(q)[1]/np.sqrt(q*(q*q-1)/12)
        F = T/2
        coefficient = 0.
        for n in range(60):
            ev = np.exp(-(2*n+1)*np.pi*lengths/beta)
            k = ev/(1+ev)
            coefficient += 2*F*F*np.sum(k[:, None]*(1-k[None, :])*abs(E)**2)
        for amplitude in [.5, 3.]:
            exact = transfer_log(expm(-1j*amplitude*F*E), lengths, beta)[0]+coefficient*amplitude**2
            refinements = []
            for N in [8, 16, 32]:
                p = np.concatenate([2*np.pi*(np.arange(-N, N)+.5)/L for L in lengths])
                i = np.repeat(np.arange(q), 2*N)
                A = amplitude*E[i[:, None], i[None, :]]*OLD['fourier'](p[None, :]-p[:, None], 0)
                A /= np.sqrt(lengths[i[:, None]]*lengths[i[None, :]])
                value = 0.
                for n in range(2*N):
                    w = (2*n+1)*np.pi/beta
                    X = A/(p[:, None]+1j*w)
                    value += 4*(np.linalg.slogdet(np.eye(len(p))+X)[1]-np.trace(X).real+np.trace(X@X).real/2)
                refinements.append(dict(spatial_cutoff=N, temporal_cutoff=2*N,
                                        log_det3=float(value), absolute_error=float(abs(value-exact))))
            errors = [x['absolute_error'] for x in refinements]
            assert errors[-1] < 1e-5 and errors[0]/errors[1] > 3 and errors[1]/errors[2] > 3
            rows.append(dict(q=q, amplitude=amplitude, transfer_plus_quadratic=float(exact), refinements=refinements))
    return rows


def covariance_scaling():
    rows = []
    for q in [2, 3, 100, 1000000]:
        K, tq = np.pi*q**4, q*(q*q-1)/12
        common = 2*tq/K
        e2 = 3*tq*(q*q+2)/(K*(q+1)**2)
        ir, cutoff = .1, 10.
        def scalar(p):
            aw, av, d = K/6, K/3, 2*K/3
            kw, kv = 1/(q+1), q/(q+1)
            return tq*(d+(av*kw*kw+aw*kv*kv)*p*p)/(p*p*(aw*av*p*p+d*(aw+av)))
        for p in [.1, 1., 10.]:
            exact = common/p**2+(e2-common)/(p*p+6)
            assert abs(scalar(p)-exact) < 1e-12
        integrated = quad(lambda logp: np.exp(2*logp)*scalar(np.exp(logp))/(2*np.pi),
                          np.log(ir), np.log(cutoff), epsabs=1e-13)[0]
        formula = common/(2*np.pi)*np.log(cutoff/ir)+(e2-common)/(4*np.pi)*np.log((cutoff**2+6)/(ir**2+6))
        assert abs(integrated-formula) < 2e-13
        rows.append(dict(q=q, one_axis_bare_variance=float(formula), integral_error=float(abs(integrated-formula))))
    scales = []
    from scipy.special import digamma
    kappa_return = 1+np.arcsinh(np.sqrt(3))/(2*np.sqrt(3))
    AT = 4*np.pi*max(2/np.pi**2, 2/np.pi*(1/(2*np.pi)+1/T))*(1+T*T/4)
    B = 8*AT+1/np.pi
    for exponent in [4, 8, 16, 32]:
        q = 10.**exponent
        logq = np.log(q)
        tq, K = q*(q*q-1)/12, np.pi*q**4
        e0, er = 1/(4*np.pi*q), 3*tq*(q*q+2)/(K*(q+1)**2)
        tau, ir, cutoff = logq**2, 1/logq**2, q**.25
        sigma2 = (e0+3*er)/(2*np.pi)*np.log(cutoff/ir)
        Lbar = T+kappa_return*(2*digamma(2*q+1)-digamma(q+1)+np.euler_gamma)/logq**2
        kappa = B/2*tau*T*sigma2
        rho = 4*np.sqrt(2)/(np.pi*np.sqrt(q))*Lbar*tau*T*3**.75*sigma2**1.5
        eta = kappa+np.sqrt(rho)*(1+1/(3*(1-rho**(1/6)))) if rho < 1 else None
        scales.append(dict(q_power_of_ten=exponent, covariance_trace_bound=float(sigma2),
                           kappa=float(kappa), rho=float(rho), eta=None if eta is None else float(eta),
                           nontrivial_eta=bool(eta is not None and eta < 1)))
    return dict(bare_reference=rows, asymptotic_substitutions=scales,
                limitations='Flat annulus reference, fixed zero modes and Gaussian gauge action; not physical matching.')


def pure_gauge_checks():
    # A_mu=partial_mu theta: zero field strength, identity temporal/spatial holonomy.
    beta, q = .3, 2
    lengths = T+np.array([.2, 2.5])
    E = axes(q)[1]/np.sqrt(q*(q*q-1)/12)
    rows = []
    for amplitude in [.4, 1.]:
        refinements = []
        for N in [4, 8, 16, 32]:
            p = np.concatenate([2*np.pi*(np.arange(-N, N)+.5)/L for L in lengths])
            i = np.repeat(np.arange(q), 2*N)
            gap = p[None, :]-p[:, None]
            theta = E[i[:, None], i[None, :]]*OLD['fourier'](gap, 0)
            theta /= np.sqrt(lengths[i[:, None]]*lengths[i[None, :]])
            derivative = -1j*gap*theta
            def rhs(t, u):
                x = 2*np.pi*t/beta
                h = np.diag(p)+amplitude*np.sin(x)*derivative
                h = h+1j*amplitude*2*np.pi/beta*np.cos(x)*theta
                return (-h@u.reshape(len(p), len(p))).ravel()
            sol = solve_ivp(rhs, [0, beta], np.eye(len(p), dtype=complex).ravel(),
                            method='DOP853', rtol=1e-11, atol=1e-12)
            assert sol.success
            logdelta = 2*np.linalg.slogdet(np.eye(len(p))+sol.y[:, -1].reshape(len(p), len(p)))[1]
            logdelta -= 2*np.sum(np.logaddexp(-beta*p, 0))+beta*amplitude**2*np.pi/(8*T)
            refinements.append(dict(cutoff=N, log_ratio=float(logdelta), exact_log_ratio=0.))
        assert abs(refinements[-1]['log_ratio']) < 3e-6
        assert abs(refinements[0]['log_ratio']) > 100*abs(refinements[-1]['log_ratio'])
        rows.append(dict(amplitude=amplitude, refinements=refinements))
    return rows


def main():
    result = dict(status='passed_global_inner_determinant_and_regulated_gauge_tail_not_physical_state',
                  unitary=unitary_checks(), static=static_checks(), dynamic=dynamic_checks(),
                  gauge_marginal=gaussian_checks(), global_det3=global_det3_checks(), covariance=covariance_scaling(),
                  pure_gauge=pure_gauge_checks(),
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
                  input_sha256={},
                  physical_exterior=False, photon_cutoff_removed=False, physical_state_constructed=False)
    for name in ['notes/determinant-tail-control.md', 'notes/nonlinear-current-determinant.md',
                 'notes/propagating-photon-polarization.md', 'notes/driven-fermion-sea.md',
                 'checks/driven_fermion_sea.py', 'notes/rotational-current-response.md']:
        result['input_sha256'][name] = hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
    (ROOT/'evidence/determinant-tail-control.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], static_cases=len(result['static']),
                         dynamic_cases=len(result['dynamic']), gaussian_cases=len(result['gauge_marginal']['gaussian_cases']))))


if __name__ == '__main__':
    main()
