"""S135.1 exact correlated free covariance, retiming and finite-width errors."""
import os
for _thread_name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_thread_name] = '1'
import hashlib
import json
import platform
import sys
from pathlib import Path
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import roots_legendre

ROOT = Path(__file__).resolve().parents[1]
TAU, P0 = .25, np.pi
T0 = P0+TAU
MAPS = [(1.02, 0.), (1., .002)]
RATIOS = [.97, 1., 1.04]
TOL = dict(point_split=2e-8, finite_identity='1e-45', work_and_time=2e-9,
           Fourier_rounding=2e-10, Fourier_refinement=2e-10, inequality_relative=2e-8)
SEPARATIONS = ['1e-4', '1e-5', '1e-6']


def step(x):
    x = np.clip(x, 0., 1.)
    return x**3*(10-15*x+6*x*x)


def step_primitive(t):
    x = np.clip(np.asarray(t)/TAU, 0., 1.)
    middle = TAU*(2.5*x**4-3*x**5+x**6)
    return np.where(np.asarray(t)<=0, 0., np.where(np.asarray(t)<TAU, middle, np.asarray(t)-TAU/2))


def base_rate(u, eta):
    return 4*eta*(step(np.asarray(u)/TAU)-step((np.asarray(u)-P0)/TAU))


def base_angle(u, eta):
    return 4*eta*(step_primitive(u)-step_primitive(np.asarray(u)-P0))


def clock(u, s, b):
    return s*np.asarray(u)+b*np.asarray(u)**3


def inverse_clock(t, s, b):
    if b == 0:
        return np.asarray(t)/s
    return 2*np.sqrt(s/(3*b))*np.sinh(np.arcsinh(3*np.asarray(t)/(2*s)*np.sqrt(3*b/s))/3)


def rate(t, eta, s, b):
    u = inverse_clock(t, s, b)
    return base_rate(u, eta)/(s+3*b*u*u)


def angle(t, eta, s, b):
    return base_angle(inverse_clock(t, s, b), eta)


def returned(t, period, eta, s, b):
    if t <= 0:
        return 0.
    return float(sum(rate(t-j*period, eta, s, b) for j in range(int(np.floor(t/period))+1)))


def pieces(fun, lower, upper, cuts):
    knots = sorted(set([lower, upper]+[float(x) for x in cuts if lower < x < upper]))
    return sum(quad(fun, left, right, epsabs=2e-12, epsrel=2e-12, limit=200)[0]
               for left, right in zip(knots[:-1], knots[1:]))


def source_cuts(s, b):
    return [float(clock(u, s, b)) for u in [0., TAU, P0, T0]]


def period_quantities(period, eta, s, b):
    end = float(clock(T0, s, b))
    count = int(np.ceil(end/period))+1
    def periodization(v):
        return float(sum(rate(v+j*period, eta, s, b) for j in range(count)))
    cuts = [cut-j*period for cut in source_cuts(s, b) for j in range(count)]
    mean = 4*eta*P0/period
    squared = pieces(lambda v: periodization(v)**2, 0., period, cuts)
    excess = pieces(lambda v: (periodization(v)-mean)**2/(4*np.pi), 0., period, cuts)
    assert abs(excess-(squared-(4*eta*P0)**2/period)/(4*np.pi)) < TOL['work_and_time']
    return excess, squared


def point_split():
    mp.mp.dps = 65
    pi, ii = mp.pi, mp.j
    ident, gen = mp.eye(2), mp.matrix([[0, -1], [1, 0]])
    def smooth(x):
        x = min(mp.mpf(1), max(mp.mpf(0), x))
        return x**3*(10-15*x+6*x*x)
    def prim(u):
        if u <= 0:
            return mp.mpf(0)
        if u >= mp.mpf('.25'):
            return u-mp.mpf('.125')
        x = 4*u
        return mp.mpf('.25')*(mp.mpf('2.5')*x**4-3*x**5+x**6)
    def th0(u, eta):
        return 4*eta*(prim(u)-prim(u-pi))
    def q0(u, eta):
        return 4*eta*(smooth(4*u)-smooth(4*(u-pi)))
    def rotation(a):
        return mp.cos(a)*ident+mp.sin(a)*gen
    def background(d, period):
        return ident/(4*period*ii*mp.sin(pi*d/(2*period)))+ii*gen/(4*period*mp.cos(pi*d/(2*period)))
    rows = []
    for s0, b0 in MAPS:
        s, b = mp.mpf(str(s0)), mp.mpf(str(b0))
        eta, period = mp.mpf('.02'), pi*mp.mpf('1.04')
        tmap = lambda u: s*u+b*u**3
        y = lambda u: s+3*b*u*u
        def inv(t):
            if b == 0:
                return t/s
            return 2*mp.sqrt(s/(3*b))*mp.sinh(mp.asinh(3*t/(2*s)*mp.sqrt(3*b/s))/3)
        def accum(t, per=period, reference=False):
            if t <= 0:
                return mp.mpf(0)
            return sum(th0(t-j*per if reference else inv(t-j*per), eta)
                       for j in range(int(mp.floor(t/per))+1))
        def dsum(t):
            if t <= 0:
                return mp.mpf(0)
            return sum(q0(inv(t-j*period), eta)/y(inv(t-j*period))
                       for j in range(int(mp.floor(t/period))+1))
        for u in map(mp.mpf, ['.12', '.8', '3.25']):
            t = tmap(u)
            dd = dsum(t)
            p = (dd*dd-pi*dd/period)/(8*pi)
            vacuum = -pi/(96*period*period)
            sch = 6*b/y(u)-mp.mpf('1.5')*(6*b*u/y(u))**2
            expected = y(u)**2*(vacuum+p)-sch/(48*pi)
            reference_d = sum(q0(u-j*pi, eta) for j in range(max(0, int(mp.floor(u/pi))+1)))
            reference_p = (reference_d**2-reference_d)/(8*pi)
            reference_expected = -1/(96*pi)+reference_p
            point_errors, global_errors, same_errors, different_errors = [], [], [], []
            for separation in SEPARATIONS:
                eps = mp.mpf(separation)
                def covariance(h, mode):
                    left, right = u+h/2, u-h/2
                    if mode == 'reference':
                        delta = left-right
                        mat = rotation(accum(left, pi, True))*background(delta, pi)*rotation(accum(right, pi, True)).T
                        return mat[0, 0]
                    tl, tr = tmap(left), tmap(right)
                    jacobian = mp.sqrt(y(left)*y(right))
                    base = background(tl-tr, period)
                    if mode == 'same_reference':
                        return jacobian*base[0, 0]
                    mat = rotation(accum(tl))*base*rotation(accum(tr)).T
                    return jacobian*mat[0, 0]
                def regular(h, mode):
                    return covariance(h, mode)-1/(2*pi*ii*h)
                actual = -ii/2*(regular(eps, 'actual')-regular(-eps, 'actual'))/(2*eps)
                ref = -ii/2*(regular(eps, 'reference')-regular(-eps, 'reference'))/(2*eps)
                same = -ii/2*((covariance(eps, 'actual')-covariance(eps, 'same_reference'))
                             -(covariance(-eps, 'actual')-covariance(-eps, 'same_reference')))/(2*eps)
                def global_regular(h):
                    mat = rotation(accum(t+h/2))*background(h, period)*rotation(accum(t-h/2)).T
                    return mat[0, 0]-1/(2*pi*ii*h)
                glob = -ii/2*(global_regular(eps)-global_regular(-eps))/(2*eps)
                point_errors.append(float(abs(actual-expected)))
                global_errors.append(float(abs(glob-vacuum-p)))
                same_errors.append(float(abs(same-y(u)**2*p)))
                different_errors.append(float(abs((actual-ref)-(expected-reference_expected))))
            assert max(point_errors[-1], global_errors[-1], same_errors[-1], different_errors[-1]) < TOL['point_split']
            finite_identity = y(u)**2*(vacuum+p)-sch/(48*pi)-expected
            assert abs(finite_identity) < mp.mpf(TOL['finite_identity'])
            assert abs(inv(t)-u) < mp.mpf(TOL['finite_identity'])
            rows.append(dict(s=s0, b=b0, lab_time=float(u), Schwarzian=float(sch),
                             canonical_lab_stress=float(expected), same_map_state_difference=float(y(u)**2*p),
                             different_map_vacuum_Schwarzian_contribution=float(y(u)**2*vacuum+1/(96*pi)-sch/(48*pi)),
                             lab_point_split_errors=point_errors, global_point_split_errors=global_errors,
                             same_map_point_split_errors=same_errors, different_map_errors=different_errors))
    return rows


def fourier_coefficients(period, eta, s, b, nodes):
    x, w = roots_legendre(nodes)
    frequencies = 2*np.pi*np.arange(1, 129)/period
    result = np.zeros(128, complex)
    for lower, upper in [(0., TAU), (TAU, P0), (P0, T0)]:
        u = (upper-lower)*x/2+(upper+lower)/2
        weights = (upper-lower)*w/2*base_rate(u, eta)
        result += np.exp(-1j*frequencies[:, None]*clock(u, s, b)[None, :])@weights
    return result


def state_and_bounds():
    rows = []
    for eta in [.0005, .02]:
        for s, b in MAPS:
            for ratio in RATIOS:
                period = ratio*P0
                excess, squared = period_quantities(period, eta, s, b)
                coefficients = fourier_coefficients(period, eta, s, b, 768)
                coarse = fourier_coefficients(period, eta, s, b, 512)
                refinement = float(np.max(abs(coefficients-coarse)))
                assert refinement < TOL['Fourier_refinement']
                theta_final = 4*eta*P0
                second_clock_max = 6*b*T0
                l1_majorant = (30*eta/TAU/s**2+3*second_clock_max*8*eta/s**3
                               +6*b*theta_final/s**3+3*second_clock_max**2*theta_final/s**4)
                parseval = []
                for harmonics in [64, 128]:
                    partial = float(np.sum(abs(coefficients[:harmonics])**2)/(2*np.pi*period))
                    tail = l1_majorant**2*period**3/(96*np.pi**5*harmonics**3)
                    assert partial <= excess+TOL['Fourier_rounding']
                    assert abs(excess-partial) <= tail+TOL['Fourier_rounding']
                    parseval.append(dict(harmonics=harmonics, partial_excitation=partial, tail_bound=tail,
                                         discrepancy=float(abs(excess-partial))))
                end = float(clock(T0, s, b))
                common_end = max(end, ratio*T0)
                u_max = float(inverse_clock(common_end, s, b))
                epsilon0 = u_max*abs(1-s/ratio)+b*u_max**3/ratio
                epsilon1 = abs(1/s-1/ratio)+3*b*u_max*u_max/s**2
                bound_rate = (15*eta/(2*TAU))/s*epsilon0+4*eta*epsilon1
                multiplicity = int(np.ceil(common_end/period))
                excitation_bound = multiplicity*common_end*bound_rate**2/(4*np.pi)
                assert excess <= excitation_bound*(1+TOL['inequality_relative'])
                # Compare the actual global rate with the period-adapted
                # reference on its full support, not just a finite Fourier set.
                adapted_cuts = source_cuts(s, b)+[ratio*x for x in [0., TAU, P0, T0]]
                rate_difference_squared = pieces(lambda t: float(rate(t, eta, s, b)-base_rate(t/ratio, eta)/ratio)**2,
                                                  0., common_end, adapted_cuts)
                l2_bound = multiplicity*rate_difference_squared/(4*np.pi)
                assert excess <= l2_bound*(1+TOL['inequality_relative'])+1e-15
                rows.append(dict(eta=eta, s=s, b=b, P_over_pi=ratio, final_angle=theta_final,
                                 actual_source_end=end, final_excitation=excess,
                                 final_energy_shift=-theta_final/(4*period)+squared/(4*np.pi),
                                 continuum_Fourier_refinement_error=refinement, Fourier= parseval,
                                 q_second_L1_majorant=l1_majorant, epsilon0=epsilon0, epsilon1=epsilon1,
                                 overlap_multiplicity=multiplicity, exact_rate_L2_excitation_bound=l2_bound,
                                 retiming_width_excitation_bound=excitation_bound))
    return rows


def work_and_compression():
    rows = []
    eta, s, b, period = .02, 1., .002, 1.04*P0
    cuts0 = source_cuts(s, b)
    def energy(t):
        cuts = [cut+j*period for cut in cuts0 for j in range(int(np.ceil(max(t, 0)/period))+1)]
        squared = pieces(lambda v: returned(v, period, eta, s, b)**2, t-period, t, cuts)
        return float(-angle(t, eta, s, b)/(4*period)+squared/(4*np.pi))
    for t in [.12, .8, period+.12, 4.]:
        q = float(rate(t, eta, s, b))
        predicted = (q*q+2*q*returned(t-period, period, eta, s, b)-np.pi*q/period)/(4*np.pi)
        errors = []
        for step_size in [1e-5, 5e-6]:
            numeric = (energy(t+step_size)-energy(t-step_size))/(2*step_size)
            errors.append(abs(numeric-predicted))
        assert max(errors) < TOL['work_and_time']
        rows.append(dict(global_time=t, work_derivative=predicted, finite_difference_errors=errors))
    compression = []
    norm_base_squared = 16*eta**2*(P0-50*TAU/231)
    for speed in [2., 4., 8.]:
        excess, _ = period_quantities(P0, eta, 1/speed, 0.)
        exact = (speed*norm_base_squared-(4*eta*P0)**2/P0)/(4*np.pi)
        assert abs(excess-exact) < TOL['work_and_time']
        compression.append(dict(speed=speed, positive_clock_rate=1/speed, final_excitation=excess,
                                exact_compression_excitation=exact))
    return dict(global_work=rows, positive_clock_insufficient=compression,
                unbounded_excitation_proved_by_exact_linear_speed_formula=True)


def main():
    names = ['notes/retimed-returned-quantum-stress.md', 'checks/retimed_returned_quantum_stress.py',
             'notes/returned-switch-capture.md', 'notes/prescribed-source-clock-variation.md',
             'notes/exterior-lll-delayed-boundary.md', 'notes/transparent-clock-flux.md']
    result = dict(story='S135.1', date='2026-10-15', status='passed_exact_retimed_free_quantum_stress_and_excitation',
                  tolerances=TOL, point_split_decimal_precision=65, point_split_separations=SEPARATIONS,
                  point_split=point_split(), quantum_states=state_and_bounds(), work_and_compression=work_and_compression(),
                  execution=dict(executable=sys.executable, script=str(Path(__file__).resolve()), cwd=str(Path.cwd())),
                  hash_policy='SHA-256 of file bytes after CRLF-to-LF normalization; no other changes',
                  input_sha256={name:hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for name in names},
                  environment=dict(python=platform.python_version(), mpmath=mp.__version__, numpy=np.__version__, scipy=scipy.__version__),
                  source_review=dict(url='https://arxiv.org/html/1804.00491v3',
                                     scope='Appendix C.1-C.3, especially C.207-C.211: real-fermion spectrum, vacuum energy, c=N/2 and AdS geometric completion',
                                     observed_equations=['cos(pi omega)=-sin(pi epsilon)', 'E=-(1+12 epsilon(1-epsilon))/48', 'c=N/2'],
                                     malformed_HTML_C206_not_imported=True),
                  claims=dict(actual_physical_clock_error_supplied=False, evolving_optical_geometry_matched=False,
                              laboratory_Schwarzian_added_to_global_JT_charge=False,
                              different_map_subtraction_replaced_by_same_map=False,
                              interacting_state_constructed=False, finite_contacts_matched=False,
                              s074_gate_closed=False, opening_achieved=False, novelty_claim=False, lean_invoked=False))
    result['summary'] = dict(point_split_cases=len(result['point_split']), state_cases=len(result['quantum_states']),
                             Fourier_tail_comparisons=2*len(result['quantum_states']),
                             maximum_finest_lab_point_split_error=max(r['lab_point_split_errors'][-1] for r in result['point_split']),
                             maximum_fourier_refinement_error=max(r['continuum_Fourier_refinement_error'] for r in result['quantum_states']))
    (ROOT/'evidence/retimed-returned-quantum-stress.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=result['status'], **result['summary'], execution=result['execution'],
                          hash_policy=result['hash_policy'], input_sha256=result['input_sha256']), indent=2))


if __name__ == '__main__':
    main()
