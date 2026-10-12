"""Independent finite-matrix and kernel checks; no interacting-state claim."""
import os
for _name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_name] = '1'
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import block_diag
from scipy.special import digamma, k1

ROOT = Path(__file__).resolve().parents[1]
T = np.pi
GAMMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]]])
KAPPA = 1+np.arcsinh(np.sqrt(3))/(2*np.sqrt(3))
RNG = np.random.default_rng(20261011)


def color(q):
    j = (q-1)/2
    m = np.arange(q)-j
    plus = np.diag(np.sqrt((j-m[:-1])*(j+m[:-1]+1)), -1)
    tq = q*(q*q-1)/12
    return np.array([np.eye(q)/np.sqrt(q), (plus+plus.T)/(2*np.sqrt(tq)),
                     (plus-plus.T)/(2j*np.sqrt(tq)), np.diag(m)/np.sqrt(tq)])


def potential(phi, E):
    C = np.einsum('ma,aij->mij', phi, E)
    return 1j*sum(np.kron(C[mu], GAMMA[mu]) for mu in range(2))


def color_checks():
    rows = []
    for q in [2, 3, 4, 8, 15]:
        E = color(q)
        assert np.max(abs(np.einsum('aij,bji->ab', E, E)-np.eye(4))) < 2e-14
        assert np.max(abs(sum(A@A for A in E)-4*np.eye(q)/q)) < 2e-14
        for _ in range(12):
            phi = RNG.normal(size=(2, 4))
            weights = np.exp(RNG.uniform(-5, 5, q))
            V = potential(phi, E)
            lam, U = np.linalg.eigh(-1j*V)
            cube = (U*abs(lam)**3)@U.conj().T
            actual = np.dot(np.repeat(weights, 2), np.diag(cube)).real
            bound = 16*np.sqrt(2)*np.mean(weights)/np.sqrt(q)*np.linalg.norm(phi)**3
            op_bound = 2*np.sqrt(2/q)*np.linalg.norm(phi)
            assert max(abs(lam)) <= op_bound*(1+2e-14)
            assert actual <= bound*(1+2e-14)
            assert abs(np.trace(V@V.conj().T).real-2*np.sum(phi**2)) < 2e-12
            rows.append(dict(q=q, weighted_ratio=float(actual/bound),
                             operator_ratio=float(max(abs(lam))/op_bound)))
    return rows


def coefficient(beta, L, N):
    modes = (2*np.arange(-N, N)+1)*np.pi
    return np.sum(((modes[:, None]/beta)**2+(modes[None, :]/L)**2)**(-1.5))/(beta*L)


def poisson_coefficient(beta, L, N=300):
    # First Poisson-sum the longer circle: exponentially convergent correction.
    short, long = sorted([beta, L])
    w = (2*np.arange(N)+1)*np.pi/short
    ell = np.arange(1, 60)
    correction = np.sum((-1.)**ell[None, :]*ell[None, :]*long /
                        (np.pi*w[:, None])*k1(w[:, None]*ell[None, :]*long))
    return short/(4*np.pi)+4*correction/short


def kernel_checks():
    rows = []
    for beta, L in [(.2, 3.), (1., 1.), (2., 7.), (13., 2.), (60., 3.)]:
        exact = poisson_coefficient(beta, L)
        bound = min(beta, L)/(4*np.pi)
        approximations = [coefficient(beta, L, N) for N in [40, 80, 160]]
        errors = [exact-x for x in approximations]
        assert 0 < exact <= bound*(1+1e-13)
        assert all(x > 0 for x in errors)
        assert errors[2] < .6*errors[1] < .36*errors[0]
        assert abs(poisson_coefficient(L, beta)-exact) < 1e-14
        rows.append(dict(beta=beta, L=L, poisson_value=float(exact), bound=float(bound),
                         square_cutoffs=[40, 80, 160], sums=list(map(float, approximations)),
                         cutoff_errors=list(map(float, errors))))
    for p in [.04, .7, 3.]:
        actual = quad(lambda w: (w*w+p*p)**(-1.5)/np.pi, 0, np.inf,
                      epsabs=1e-10, epsrel=1e-12)[0]
        assert abs(actual-1/(np.pi*p*p)) < 2e-9
    return rows


def spatial(k, which):
    def window(x):
        return T*np.sinc(x*T/(2*np.pi))
    def f(x):
        return .5*window(x)+.25*(window(x+2*np.pi/T)+window(x-2*np.pi/T))
    return f(k) if which == 0 else f(k)+.3*(f(k+2*np.pi/T)-f(k-2*np.pi/T))/(2j)


def profile(t, s, beta, amplitude):
    x = 2*np.pi*t/beta
    f = np.cos(np.pi*s/T)**2
    g = f*(1+.3*np.sin(2*np.pi*s/T))
    return amplitude*np.array([[0, f*np.cos(x), 0, .3*f*np.sin(2*x)],
                                [.2*g*np.cos(2*x), 0, .4*g*np.sin(x), g]])


def dirac_matrix(q, N, lengths, beta, amplitude):
    E = color(q)
    basis = [(i, n, (2*r+1)*np.pi/lengths[i])
             for i in range(q) for n in range(-N, N) for r in range(-N, N)]
    labels = np.array([b[0] for b in basis])
    ns = np.array([b[1] for b in basis])
    ps = np.array([b[2] for b in basis])
    k = ps[None, :]-ps[:, None]
    shift = ns[:, None]-ns[None, :]
    norm = np.sqrt(lengths[labels, None]*lengths[None, labels])
    # (mu, color, spatial profile, Fourier coefficients of periodic source)
    terms = [(0, 1, 0, {1: .5, -1: .5}),
             (0, 3, 0, {2: .3/(2j), -2: -.3/(2j)}),
             (1, 3, 1, {0: 1}), (1, 2, 1, {1: .4/(2j), -1: -.4/(2j)}),
             (1, 0, 1, {2: .1, -2: .1})]
    C = [np.zeros_like(k, dtype=complex) for _ in range(2)]
    for mu, a, which, coefficients in terms:
        temporal = sum(coef*(shift == m) for m, coef in coefficients.items())
        C[mu] += amplitude*E[a][labels[:, None], labels[None, :]]*spatial(k, which)*temporal/norm
    V = 1j*sum(np.kron(C[mu], GAMMA[mu]) for mu in range(2))
    G = block_diag(*[-1j*(GAMMA[0]*((2*n+1)*np.pi/beta)+GAMMA[1]*p) /
                    (((2*n+1)*np.pi/beta)**2+p*p) for _, n, p in basis])
    assert np.linalg.norm(V+V.conj().T) < 1e-12
    return G@V, V, basis


def integrate_bounds(q, lengths, beta, amplitude, N, order=32):
    E = color(q)
    zs, ws = np.polynomial.legendre.leggauss(order)
    cs = np.array([coefficient(beta, L, N) for L in lengths])
    weighted, integral = 0., 0.
    for zt, wt in zip(zs, ws):
        for zs_, ws_ in zip(zs, ws):
            phi = profile(beta*(zt+1)/2, T*zs_/2, beta, amplitude)
            V = potential(phi, E)
            lam, U = np.linalg.eigh(-1j*V)
            diag_cube = (abs(U)**2)@(abs(lam)**3)
            weight = wt*ws_*beta*T/4
            weighted += weight*np.dot(np.repeat(cs, 2), diag_cube)
            integral += weight*np.linalg.norm(phi)**3
    r = 4*np.sqrt(2)/(np.pi*np.sqrt(q))*np.mean(np.minimum(beta, lengths))*integral
    alternative = np.sqrt(2)*beta/(np.pi*np.sqrt(q))*integral
    return float(weighted), float(min(r, alternative)), float(integral)


def independent_matrix(q, lengths, beta, amplitude, basis):
    # Quadrature of the original local multiplication operator; no analytic Fourier entries.
    z, w = np.polynomial.legendre.leggauss(28)
    E = color(q)
    result = np.zeros((2*len(basis), 2*len(basis)), dtype=complex)
    for zt, wt in zip(z, w):
        t = beta*(zt+1)/2
        for zs, ws in zip(z, w):
            s = T*zs/2
            B = np.zeros((2*q, 2*len(basis)), dtype=complex)
            for b, (i, n, p) in enumerate(basis):
                v = np.exp(1j*((2*n+1)*np.pi*t/beta+p*s))/np.sqrt(beta*lengths[i])
                B[2*i:2*i+2, 2*b:2*b+2] = v*np.eye(2)
            result += wt*ws*beta*T/4*(B.conj().T@potential(profile(t, s, beta, amplitude), E)@B)
    return result


def determinant_checks():
    rows = []
    for q in [2, 3, 4]:
        lengths = T+np.linspace(.2, 3.1, q)**1.5
        beta, amplitude = 2.3, .25
        for N in [1, 2, 3]:
            X, V, basis = dirac_matrix(q, N, lengths, beta, amplitude)
            singular = np.linalg.svd(X, compute_uv=False)
            s3, op = float(np.sum(singular**3)), float(singular[0])
            weighted, r, integral = integrate_bounds(q, lengths, beta, amplitude, N)
            assert s3 <= weighted*(1+2e-6) and weighted <= r*(1+2e-6)
            assert r < 1 and op < 1
            sign, logabs = np.linalg.slogdet(np.eye(len(X))+X)
            direct = logabs+1j*np.angle(sign)-np.trace(X)+np.trace(X@X)/2
            power = X@X
            series = 0j
            for n in range(3, 15):
                power = power@X
                series += (-1)**(n+1)*np.trace(power)/n
            tail = s3*op**12/(15*(1-op))
            assert abs(direct-series) <= tail+3e-13
            assert abs(direct) <= r/(3*(1-r**(1/3)))
            row = dict(q=q, cutoff=N, matrix_dimension=len(X), schatten_cube=s3,
                       operator_norm=op, coincident_bound=weighted, source_bound=r,
                       source_L3_integral=integral, logdet3_real=float(direct.real),
                       logdet3_imag=float(direct.imag), series_error=float(abs(direct-series)),
                       trace_series_tail_bound=float(tail),
                       nonnormality=float(np.linalg.norm(X@X.conj().T-X.conj().T@X)))
            if q == 2 and N == 1:
                error = np.linalg.norm(V-independent_matrix(q, lengths, beta, amplitude, basis))
                assert error < 1e-12
                row['independent_spacetime_matrix_error'] = float(error)
            rows.append(row)
    return rows


def gaussian_and_scaling():
    # Finite Gaussian examples only, not the covariance of the physical action.
    rows = []
    for rank in [1, 3, 8]:
        covariance = np.diag(np.linspace(.2, 1., rank))
        samples = RNG.normal(size=(120000, rank))@np.sqrt(covariance)
        estimate = np.mean(np.linalg.norm(samples, axis=1)**3)
        bound = 3**.75*np.trace(covariance)**1.5
        assert estimate < bound
        rows.append(dict(rank=rank, empirical_third_moment=float(estimate), bound=float(bound)))
    scales = []
    for exponent in [4, 8, 16, 32, 64]:
        q = 10.**exponent
        logq = np.log(q)
        # H_n=digamma(n+1)+EulerGamma; q here is integral in mathematical formula.
        h = 2*digamma(2*q+1)-digamma(q+1)+np.euler_gamma
        Lbar_bound = T+KAPPA*h/logq**2
        duration = logq**2
        common = 4*np.sqrt(2)/(np.pi*np.sqrt(q))*Lbar_bound*duration*T
        rotating_amplitude = logq**3/np.sqrt(48*np.pi)*np.sqrt(1-1/q**2)
        scales.append(dict(q_power_of_ten=exponent, Lbar_bound=float(Lbar_bound),
                           unit_profile_r=float(common),
                           moving_profile_r=float(common*rotating_amplitude**3),
                           assumed_sigma_squared_logq_rho=float(common*3**.75*logq**1.5)))
    return dict(finite_gaussian_examples=rows, asymptotic_substitutions=scales,
                actual_gauge_covariance_derived=False, weighted_tail_control=False)


def rotational_reference():
    rows = []
    for q in [2, 3, 10, 100]:
        K, tq = np.pi*q**4, q*(q*q-1)/12
        aw, av, d = K/6, K/3, 2*K/3
        k = np.array([1., q])/(q+1)
        e2 = tq*(k[0]**2/aw+k[1]**2/av)
        assert abs(e2-3*tq*(q*q+2)/(K*(q+1)**2)) < 1e-14
        errors, ratios = [], []
        for p in [.003, .03, .3, 3., 30.]:
            B = p*p*np.diag([aw, av])+d*np.array([[1., -1.], [-1., 1.]])
            c0 = tq*k@np.linalg.solve(B, k)
            schur = c0/(1+c0/np.pi)
            direct = tq*k@np.linalg.solve(B+(tq/np.pi)*np.outer(k, k), k)
            assert abs(direct-schur) <= 2e-9
            upper = e2/(p*p+e2/np.pi)
            assert 0 < direct <= upper*(1+2e-9)
            errors.append(float(abs(direct-schur)))
            ratios.append(float(direct/upper))
        cutoff = q**2
        integrated = quad(lambda t: 3*e2/(4*np.pi)/(t+e2/np.pi), 0, cutoff**2,
                          epsabs=1e-11, epsrel=1e-11)[0]
        variance = 3*e2/(4*np.pi)*np.log1p(np.pi*cutoff**2/e2)
        assert abs(integrated-variance) < 1e-10
        rows.append(dict(q=q, effective_coupling_squared=float(e2), cutoff=cutoff,
                         rotational_variance_bound=float(variance),
                         maximum_schur_error=max(errors), covariance_bound_ratios=ratios,
                         integrated_bound_error=float(abs(integrated-variance))))
    return dict(model='Flat infinite-plane, R=1, transverse quadratic rotation reference; no actual finite returns or global modes', cases=rows)


def main():
    result = dict(status='passed_nonlinear_inner_current_bound_not_coupled_state',
                  color=color_checks(), kernels=kernel_checks(), determinants=determinant_checks(),
                  gaussian_and_scaling=gaussian_and_scaling(), rotational_reference=rotational_reference(),
                  limits='Finite tests support the ordinary proof; no gauge covariance, weighted-tail, stress, real-time or physical matching conclusion.',
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
                  input_sha256={})
    for name in ['notes/nonlinear-current-determinant.md', 'notes/driven-fermion-sea.md',
                 'notes/angular-current-coset.md', 'notes/rotational-current-response.md']:
        result['input_sha256'][name] = hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
    (ROOT/'evidence/nonlinear-current-determinant.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], color_cases=len(result['color']),
                         kernel_cases=len(result['kernels']), determinant_cases=len(result['determinants']),
                         maximum_series_error=max(x['series_error'] for x in result['determinants']))))


if __name__ == '__main__':
    main()
