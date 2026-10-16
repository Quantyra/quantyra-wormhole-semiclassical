"""S133 fixed-metric paired-angle map; no physical gap, state or clock claim."""
import os
for _thread_name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_thread_name] = '1'
import hashlib
import json
import platform
import sys
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.special import gamma, gammaincc

ROOT = Path(__file__).resolve().parents[1]
TOL = dict(identity=2e-11, local_difference=2e-7, determinant_first=2e-6,
           determinant_second=2e-5, series_rounding=5e-11,
           inequality_relative=5e-10, omitted_normal_minimum=1e-4)
STEPS = [2e-4, 1e-4]
RNG = np.random.default_rng(2026101501)
I2 = np.eye(2)
P1 = np.array([[0, 1], [1, 0]], complex)
P2 = np.array([[0, -1j], [1j, 0]])
P3 = np.diag([1., -1.])


def kron(a, b, c):
    return np.kron(np.kron(a, b), c)


GAMMA = [kron(P3, P1, I2), -kron(I2, P3, I2),
         -kron(I2, P2, P1), -kron(I2, P2, P2)]
GEN = 1j*kron(P2, P1, I2)
EYE = np.eye(8)


def check_identity(a):
    error = float(np.linalg.norm(a))
    assert error < TOL['identity'], error
    return error


def cutoff(r):
    # Smooth plateau, smooth transition, then zero before the collar join.
    u = (r-.2)/.6
    if u <= 0:
        return 1., 0.
    if u >= 1:
        return 0., 0.
    a, b = np.exp(-1/u), np.exp(-1/(1-u))
    fraction = a/(a+b)
    derivative = a*b*(1/u**2+1/(1-u)**2)/(a+b)**2/.6
    return 1-fraction, -derivative


KAPPA, HORIZON = .08, 2.


def theta(t):
    return .3+.1*np.sin(t)+KAPPA*(.2*t+.02*np.sin(2*t))


def rate(t):
    return .1*np.cos(t)+KAPPA*(.2+.04*np.cos(2*t))


def direction(t):
    return np.sin(.7*t)+.1


def primitive(t):
    return (1-np.cos(.7*t))/.7+.1*t


def rotation(r, t):
    return expm(cutoff(r)[0]*theta(t)*GEN/2)


def insertion(r, t):
    f, fp = cutoff(r)
    return (GAMMA[0]@GEN*fp*theta(t)+GAMMA[1]@GEN*f*rate(t))/2


def domain_and_transfer():
    for i in range(4):
        check_identity(GEN@GAMMA[i]+GAMMA[i]@GEN)
        for j in range(4):
            check_identity(GAMMA[i]@GAMMA[j]+GAMMA[j]@GAMMA[i]-2*(i==j)*EYE)
    check_identity(GEN.conj().T+GEN)
    check_identity(GEN@GEN+EYE)
    p0 = (EYE+GAMMA[0])/2
    rows = []
    for angle in [-.7, .1, 1.1]:
        u = expm(angle*GEN/2)
        p = u@p0@u.conj().T
        p_adj = EYE-GAMMA[0]@p@GAMMA[0]
        error = check_identity(p_adj-u.conj().T@(EYE-p0)@u)
        v = RNG.normal(size=8)+1j*RNG.normal(size=8)
        psi, adjpsi = u@p0@v, u.conj().T@(EYE-p0)@v
        check_identity((EYE-p)@psi)
        check_identity((EYE-p_adj)@adjpsi)
        check_identity((EYE-p0)@u.conj().T@psi)
        check_identity(p0@u@adjpsi)
        rows.append(dict(angle=angle, adjoint_projector_error=error))
    jets = []
    for r, t in [(.1, .3), (.35, .7), (.55, 1.3), (.9, 1.7)]:
        u = rotation(r, t)
        f, fp = cutoff(r)
        psi = RNG.normal(size=8)+1j*RNG.normal(size=8)
        derivatives = RNG.normal(size=(4, 8))+1j*RNG.normal(size=(4, 8))
        d_alpha = [fp*theta(t), f*rate(t), 0., 0.]
        d_u = [GEN@u*x/2 for x in d_alpha]
        direct = u@sum(GAMMA[i]@(d_u[i]@psi+u@derivatives[i]) for i in range(4))
        expected = sum(GAMMA[i]@derivatives[i] for i in range(4))+insertion(r, t)@psi
        algebra_error = check_identity(direct-expected)
        eps = 1e-5
        numerical_dr = (rotation(r+eps, t)-rotation(r-eps, t))/(2*eps)
        numerical_dt = (rotation(r, t+eps)-rotation(r, t-eps))/(2*eps)
        finite = u@(GAMMA[0]@numerical_dr+GAMMA[1]@numerical_dt)
        fd_error = float(np.linalg.norm(finite-insertion(r, t)))
        assert fd_error < TOL['local_difference'], fd_error
        jets.append(dict(r=r, t=t, algebra_error=algebra_error, finite_difference_error=fd_error))
    r, t = .45, .7
    f, fp = cutoff(r)
    omitted = float(np.linalg.norm(GAMMA[0]@GEN*fp*theta(t)/2))
    assert omitted > TOL['omitted_normal_minimum']
    # Pure algebraic singular-spectrum test; not a numerical boundary solver.
    a = RNG.normal(size=(8, 8))+1j*RNG.normal(size=(8, 8))
    u = rotation(.4, .9)
    transferred = u@a@u
    singular_error = check_identity(np.linalg.svd(a, compute_uv=False)-np.linalg.svd(transferred, compute_uv=False))
    check_identity(transferred.conj().T@transferred-u.conj().T@a.conj().T@a@u)
    return dict(boundary_domains=rows, local_jets=jets,
                omitted_normal_term_norm=omitted, unitary_singular_value_error=singular_error)


def realized_angle():
    rows = []
    for r, t in [(.1, .3), (.35, .7), (.55, 1.3)]:
        f, fp = cutoff(r)
        w_prime = KAPPA*(GAMMA[1]@GEN*f*direction(t)+GAMMA[0]@GEN*fp*primitive(t))/2
        ceiling = abs(KAPPA)/2*(abs(f)+HORIZON*abs(fp))*max(abs(direction(x)) for x in np.linspace(0, HORIZON, 1001))
        assert np.linalg.norm(w_prime, 2) <= ceiling*(1+TOL['inequality_relative'])
        # Independent derivative through the realized angle and its rate.
        def perturbed(epsilon):
            return (GAMMA[0]@GEN*fp*(theta(t)+epsilon*KAPPA*primitive(t))
                    +GAMMA[1]@GEN*f*(rate(t)+epsilon*KAPPA*direction(t)))/2
        fd = (perturbed(1e-5)-perturbed(-1e-5))/(2e-5)
        error = check_identity(fd-w_prime)
        rows.append(dict(r=r, t=t, insertion_norm=float(np.linalg.norm(w_prime, 2)),
                         window_ceiling=float(ceiling), source_difference_error=error))
    # Completion identity d(eta*Ph)=eta*h+eta'*Ph at an exterior point,
    # where h=0 and Ph is its late constant. Example eta is a smooth local
    # profile; analytic completion bound, not a physical time contour.
    constant, tt = .4, 2.5
    eta = lambda x: np.exp(-(x-2.)**2)
    derivative = (eta(tt+1e-5)*constant-eta(tt-1e-5)*constant)/(2e-5)
    exact = -2*(tt-2)*eta(tt)*constant
    assert abs(derivative-exact) < TOL['local_difference']
    return dict(kappa=KAPPA, H=HORIZON, local_cases=rows,
                completion_derivative_error=float(abs(derivative-exact)),
                completion_and_Euclidean_state_are_additional_data=True)


def s5norm(a):
    singular = np.linalg.svd(a, compute_uv=False)
    return float(np.sum(singular**5)**.2), float(singular[0])


def f5(k):
    # Eigenvalue log for the analytic branch in ||K||op<1, independently
    # compared with the finite regularized determinant subtraction.
    result = np.sum(np.log1p(np.linalg.eigvals(k)))
    power = np.eye(len(k), dtype=complex)
    for n in range(1, 5):
        power = power@k
        result -= (-1)**(n+1)*np.trace(power)/n
    return result


def df5(k, h):
    return np.trace(np.linalg.matrix_power(k, 4)@np.linalg.solve(np.eye(len(k))+k, h))


def d2f5(k, j, h):
    r = np.linalg.inv(np.eye(len(k))+k)
    ordered = sum(np.linalg.matrix_power(k, l)@j@np.linalg.matrix_power(k, 3-l) for l in range(4))
    return np.trace((ordered@r-np.linalg.matrix_power(k, 4)@r@j@r)@h)


def determinant_cases():
    rows = []
    for dim, target in [(4, .25), (8, .55), (12, .8)]:
        k = RNG.normal(size=(dim, dim))+1j*RNG.normal(size=(dim, dim))
        k *= target/s5norm(k)[1]
        h = RNG.normal(size=(dim, dim))+1j*RNG.normal(size=(dim, dim))
        j = RNG.normal(size=(dim, dim))+1j*RNG.normal(size=(dim, dim))
        h *= .2/s5norm(h)[1]
        j *= .15/s5norm(j)[1]
        s, a = s5norm(k)
        hh, _ = s5norm(h)
        jj, _ = s5norm(j)
        d, d2 = df5(k, h), d2f5(k, j, h)
        first_bound, second_bound = s**4*hh/(1-a), s**3*(4-3*a)*hh*jj/(1-a)**2
        assert abs(d) <= first_bound
        assert abs(d2) <= second_bound
        check_identity(d2-d2f5(k, h, j))
        differences = []
        for step in STEPS:
            first = (f5(k+step*h)-f5(k-step*h))/(2*step)
            second = (f5(k+step*h+step*j)-f5(k+step*h-step*j)
                      -f5(k-step*h+step*j)+f5(k-step*h-step*j))/(4*step**2)
            errors = [float(abs(first-d)), float(abs(second-d2))]
            assert errors[0] < TOL['determinant_first'], errors
            assert errors[1] < TOL['determinant_second'], errors
            differences.append(dict(step=step, first_error=errors[0], mixed_error=errors[1]))
        tails = []
        for degree in [5, 8, 12]:
            powers = [np.linalg.matrix_power(k, n) for n in range(degree+1)]
            series = sum((-1)**(n+1)*np.trace(powers[n])/n for n in range(5, degree+1))
            first_series = sum((-1)**(n+1)*np.trace(powers[n-1]@h) for n in range(5, degree+1))
            mixed_series = sum((-1)**(n+1)*sum(np.trace(powers[l]@j@powers[n-2-l]@h) for l in range(n-1)) for n in range(5, degree+1))
            errors = list(map(float, [abs(f5(k)-series), abs(d-first_series), abs(d2-mixed_series)]))
            bounds = [s**5*a**(degree-4)/((degree+1)*(1-a)),
                      s**4*hh*a**(degree-4)/(1-a),
                      s**3*hh*jj*a**(degree-4)*(degree/(1-a)+a/(1-a)**2)]
            assert all(error <= bound*(1+TOL['inequality_relative'])+TOL['series_rounding'] for error, bound in zip(errors, bounds))
            tails.append(dict(degree=degree, errors=errors, bounds=bounds))
        # Full finite determinant reconstructs all orders 1..4, not zero.
        sign, logabs = np.linalg.slogdet(np.eye(dim)+k)
        independent = logabs+1j*np.angle(sign)
        lower = sum((-1)**(n+1)*np.trace(np.linalg.matrix_power(k, n))/n for n in range(1, 5))
        reconstruction = float(abs(independent-f5(k)-lower))
        assert reconstruction < TOL['identity']
        rows.append(dict(dimension=dim, operator_norm=a, first_absolute=float(abs(d)),
                         first_bound=first_bound, mixed_absolute=float(abs(d2)), mixed_bound=second_bound,
                         differences=differences, tails=tails, full_determinant_reconstruction_error=reconstruction,
                         lower_order_sum_absolute=float(abs(lower))))
    return rows


def heat_bound():
    rows = []
    for gap in [.4, 1., 2.]:
        c = 3.
        analytic = c/gamma(2.5)*(2+np.exp(gap**2)*gap**(-5)*gamma(2.5)*gammaincc(2.5, gap**2))
        direct = c/gamma(2.5)*(quad(lambda t: t**(-.5), 0, 1)[0]
                                  +quad(lambda t: t**1.5*np.exp(-gap**2*(t-1)), 1, np.inf)[0])
        assert abs(analytic-direct) < 1e-9
        rows.append(dict(assumed_C=c, assumed_gap=gap, B5_power5=float(analytic), integration_error=float(abs(analytic-direct))))
    return dict(reference_constants_are_illustrative_not_physical=True, cases=rows)


def conditional_geometry_interface():
    dim = 8
    def matrix(scale):
        value = RNG.normal(size=(dim, dim))+1j*RNG.normal(size=(dim, dim))
        return value*scale/np.linalg.norm(value, 2)
    d0, d1, d2 = 2*np.eye(dim)+matrix(.3), matrix(.2), matrix(.15)
    w0, w1, w2 = matrix(.3), matrix(.1), matrix(.08)
    r = np.linalg.inv(d0)
    rp = -r@d1@r
    rpp = 2*r@d1@r@d1@r-r@d2@r
    kp = rp@w0+r@w1
    kpp = rpp@w0+2*rp@w1+r@w2
    q, _ = s5norm(r)
    b1, b2 = np.linalg.norm(d1@r, 2), np.linalg.norm(d2@r, 2)
    wn = [np.linalg.norm(w, 2) for w in [w0, w1, w2]]
    first_bound = q*(b1*wn[0]+wn[1])
    second_bound = q*((2*b1*b1+b2)*wn[0]+2*b1*wn[1]+wn[2])
    assert s5norm(kp)[0] <= first_bound
    assert s5norm(kpp)[0] <= second_bound
    def at(z):
        return np.linalg.solve(d0+z*d1+z*z*d2/2, w0+z*w1+z*z*w2/2)
    rows = []
    for step in STEPS:
        first_error = float(np.linalg.norm((at(step)-at(-step))/(2*step)-kp))
        second_error = float(np.linalg.norm((at(step)-2*at(0)+at(-step))/step**2-kpp))
        assert first_error < TOL['local_difference']
        assert second_error < TOL['determinant_second']
        rows.append(dict(step=step, first_error=first_error, second_error=second_error))
    return dict(first_S5_norm=s5norm(kp)[0], first_bound=float(first_bound),
                second_S5_norm=s5norm(kpp)[0], second_bound=float(second_bound), differences=rows,
                actual_wormhole_graph_domain_or_relative_derivative_bounds_derived=False)


def spectral_and_lapse():
    c = .2
    fourier = []
    for cutoff_n in [2, 4, 8, 16]:
        grid = np.indices((2*cutoff_n+1,)*4, dtype=float).reshape(4, -1)-cutoff_n
        norm = np.sqrt(np.sum(grid*grid, axis=0))
        norm = norm[norm > 0]
        sums = {str(p):float(8*np.sum((c/norm)**p)) for p in [3, 4, 5]}
        fourier.append(dict(maxnorm_cutoff=cutoff_n, mode_count=len(norm), schatten_power_sums=sums,
                            analytic_S5_tail_ceiling=4096*c**5/cutoff_n))
    blocks = []
    for radius in [1, 2, 4, 8]:
        grid = np.indices((radius,)*4, dtype=float).reshape(4, -1)+radius
        norm = np.sqrt(np.sum(grid*grid, axis=0))
        assert len(norm) == radius**4
        for p in [3, 4, 5]:
            actual = float(8*np.sum((c/norm)**p))
            lower = 8*c**p/4**p*radius**(4-p)
            upper = 8*c**p/2**p*radius**(4-p)
            assert lower <= actual <= upper*(1+1e-13)
            blocks.append(dict(radius=radius, p=p, sum=actual, lower=lower, upper=upper))
    mode_errors = []
    w = c*GAMMA[1]@GEN
    for k in [np.array([1., 0, 0, 0]), np.array([2., -1, 3, 1]), np.array([9., 1, 2, -3])]:
        d = 1j*sum(GAMMA[i]*k[i] for i in range(4))
        singular = np.linalg.svd(np.linalg.solve(d, w), compute_uv=False)
        mode_errors.append(check_identity(singular-c/np.linalg.norm(k)))
    lapse = []
    n_lapse = 1.5
    for temporal in [1, 2, 4, 8, 16]:
        d1 = 1j*temporal*GAMMA[1]
        dn = d1/n_lapse
        relative = np.linalg.solve(d1, dn-d1)
        error = check_identity(relative-(1/n_lapse-1)*EYE)
        lapse.append(dict(temporal_mode=temporal, minimum_singular_value=float(np.linalg.svd(relative, compute_uv=False)[-1]), identity_error=error))
    return dict(source_strength=c, Fourier_partial_sums=fourier, positive_dyadic_blocks=blocks,
                Fourier_block_matrix_errors=mode_errors, lapse_N=n_lapse, lapse_modes=lapse,
                noncompact_lapse_proved_by_infinite_temporal_sequence=True,
                S3_S4_divergence_proved_by_infinite_disjoint_blocks=True,
                S5_convergence_proved_by_analytic_shell_tail=True,
                physical_wormhole_spectral_theorem=False)


def main():
    names = ['checks/four_dimensional_boundary_response.py', 'notes/four-dimensional-boundary-response.md',
             'notes/coupled-fermion-domain.md', 'notes/curved-fermion-boundary.md',
             'notes/source-clock-radiation.md', 'notes/source-clock-joint.md',
             'notes/determinant-response-control.md']
    result = dict(story='S133', date='2026-10-15', status='passed_conditional_4D_source_transfer_metric_obstruction',
                  tolerances=TOL, finite_difference_steps=STEPS, seed=2026101501,
                  domains_and_transfer=domain_and_transfer(), realized_angle=realized_angle(),
                  determinant_response=determinant_cases(), heat_bound=heat_bound(), spectral_and_lapse=spectral_and_lapse(),
                  conditional_geometry_interface=conditional_geometry_interface(),
                  execution=dict(executable=sys.executable, script=str(Path(__file__).resolve()), cwd=str(Path.cwd())),
                  hash_policy='SHA-256 of file bytes after CRLF-to-LF normalization; no other changes',
                  input_sha256={name:hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for name in names},
                  source_review=[dict(url='https://arxiv.org/html/2308.00871v2', scope='Section 2.1.2 equations 2.8 through 2.12: physical paired boundary data', observed_equation='Psi_+^R+h Psi_+^L=0; Psi_-^L+h Psi_-^R=0'),
                                 dict(url='https://arxiv.org/html/hep-th/0306138v3', scope='Equation 2.21 and section 5.3 mixed heat expansion; physical gap not inferred', observed_equation='Tr(f exp(-tD)) asymptotic sum_k t^((k-n)/2) a_k(f,D)'),
                                 dict(url='https://arxiv.org/html/2007.12834v3', scope='Equation 1.3 regularized determinant, A=-K and k=5', observed_equation='det_k(I-A)=det((I-A) exp(sum_(ell=1)^(k-1) A^ell/ell))'),
                                 dict(url='https://ems.press/content/serial-article-files/31026', scope='Published theorem 6.1, p.748: closed-manifold negative-order singular-value law; no boundary-exterior theorem imported', observed_equation='P in Psi^(-m), m>0, p=n/m; lim j^(1/p) mu_j(P) equals principal-symbol integral power 1/p')],
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  claims=dict(physical_gap_or_heat_constant_derived=False, physical_Euclidean_state=False,
                              full_clock_metric_response_controlled=False, det5_lower_orders_set_to_zero=False,
                              transfer_phase_inferred=False, interacting_gauge_measure_controlled=False,
                              s074_gate_closed=False, opening_achieved=False, novelty_claim=False, lean_invoked=False))
    result['summary'] = dict(domain_cases=len(result['domains_and_transfer']['boundary_domains']),
                             local_jet_cases=len(result['domains_and_transfer']['local_jets']),
                             determinant_cases=len(result['determinant_response']), trace_tail_checks=27,
                             maximum_first_difference_error=max(d['first_error'] for r in result['determinant_response'] for d in r['differences']),
                             maximum_mixed_difference_error=max(d['mixed_error'] for r in result['determinant_response'] for d in r['differences']))
    (ROOT/'evidence/four-dimensional-boundary-response.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=result['status'], **result['summary'], execution=result['execution'],
                          hash_policy=result['hash_policy'], input_sha256=result['input_sha256']), indent=2))


if __name__ == '__main__':
    main()
