"""S134.2 exact free LLL cut response; no interacting/state claim."""
import os
for _threads in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_threads] = '1'
import hashlib
import json
import platform
import sys
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
TOL = dict(identity=3e-11, characteristics=2e-10, source_and_clock=2e-9,
           inequality_relative=2e-8, incorrect_model_minimum=1e-4)
DELAYS = [np.array([.07, .19]), np.array([.11, .31])]
SCALES = [1., .5, .25]
P1 = np.array([[0, 1], [1, 0]], complex)
P2 = np.array([[0, -1j], [1j, 0]])
P3 = np.diag([1., -1.])
EYE = np.eye(2)
JROT = np.array([[0., -1.], [1., 0.]])
WIDTH = 1.7


def assert_small(value, tolerance=None):
    error = float(np.linalg.norm(value))
    assert error < (TOL['identity'] if tolerance is None else tolerance), error
    return error


def theta(t):
    return .23+.12*np.sin(.7*t)


def physical_source(t):
    th = theta(t)
    s = -1j*np.cos(th)*EYE-np.sin(th)*P1
    th1, th2 = .084*np.cos(.7*t), -.0588*np.sin(.7*t)
    sp = -1j*th1*P1@s
    spp = (-1j*th2*P1-th1**2*EYE)@s
    return s, sp, spp


def general_source(t):
    # A general unitary comparison tests the theorem beyond commuting
    # one-angle sources; it is not substituted for the paired coupling.
    s, sp, spp = physical_source(t)
    phi, phi1, phi2 = .1*np.sin(.4*t), .04*np.cos(.4*t), -.016*np.sin(.4*t)
    u = np.cos(phi)*EYE-1j*np.sin(phi)*P3
    up, upp = -1j*phi1*P3@u, (-1j*phi2*P3-phi1**2*EYE)@u
    return s@u, sp@u+s@up, spp@u+2*sp@up+s@upp


def packet(t):
    if abs(t) >= WIDTH:
        return np.zeros(2, complex), np.zeros(2, complex), np.zeros(2, complex)
    denominator = 1-(t/WIDTH)**2
    f = np.exp(-1/denominator)
    l1 = -2*t/WIDTH**2/denominator**2
    l2 = -2/WIDTH**2/denominator**2-8*t*t/WIDTH**4/denominator**3
    fp, fpp = f*l1, f*(l1*l1+l2)
    frequencies = np.array([.7, -1.1])
    vector = np.array([1., .3+.4j])*np.exp(-1j*frequencies*t)
    return f*vector, (fp-1j*frequencies*f)*vector, (fpp-2j*frequencies*fp-frequencies**2*f)*vector


def response(t, delays, source=physical_source):
    return np.array([sum(source(t-delays[i])[0][i, j]*packet(t-delays[i]-delays[j])[0][j]
                         for j in range(2)) for i in range(2)])


def direct_characteristics(t, delays, source):
    output = []
    for i in range(2):
        scattering_time = t-delays[i]
        incident_at_boundary = np.array([packet(scattering_time-travel)[0][j] for j, travel in enumerate(delays)])
        boundary_return = source(scattering_time)[0]@incident_at_boundary
        output.append(boundary_return[i])
    return np.array(output)


def norm2(function, limit, tolerance=1e-10):
    points = sorted(set([-WIDTH, WIDTH]+[float(x) for x in limit[2:]]))
    left, right = limit[:2]
    points = [x for x in points if left < x < right]
    value, error = quad(lambda t: float(np.vdot(function(t), function(t)).real),
                        left, right, points=points, epsabs=tolerance, epsrel=tolerance, limit=300)
    return np.sqrt(max(0., value)), error


def projector_and_spectrum():
    e_plus, e_minus = np.array([1, 1j])/np.sqrt(2), np.array([1, -1j])/np.sqrt(2)
    incident = np.column_stack([np.r_[e_minus, [0., 0.]], np.r_[[0., 0.], e_plus]])
    returning = np.column_stack([np.r_[e_plus, [0., 0.]], np.r_[[0., 0.], e_minus]])
    rows = []
    for angle in [-.4, .2, .9]:
        chi = np.cos(angle)*np.kron(P3, P1)-np.sin(angle)*np.kron(P1, EYE)
        s = -1j*np.cos(angle)*EYE-np.sin(angle)*P1
        projected = returning.conj().T@chi@incident
        error = assert_small(projected-s)
        forbidden = (np.eye(4)-chi)/2
        incoming_amplitude = np.array([.3+.4j, -.7+.1j])
        solved = np.linalg.lstsq(forbidden@returning, -forbidden@incident@incoming_amplitude, rcond=None)[0]
        assert_small(solved-s@incoming_amplitude)
        assert_small(s.conj().T@s-EYE)
        rephase = np.diag([1., 1j])
        assert_small(rephase@(1j*s)@rephase.conj().T-expm(angle*JROT))
        rows.append(dict(angle=angle, projector_scattering_error=error))
    spectral = []
    angle = .21
    s = -1j*np.cos(angle)*EYE-np.sin(angle)*P1
    for delays in DELAYS:
        central = np.pi-np.sum(delays)
        period = central+np.sum(delays)
        for omega in [.3, .8, 1.7]:
            delay = np.diag(np.exp(1j*omega*delays))
            actual = np.linalg.det(EYE-np.exp(1j*omega*central)*delay@s@delay@P1)
            exact = 2*np.exp(1j*omega*period)*(np.cos(omega*period)+np.sin(angle))
            spectral.append(dict(delays=delays.tolist(), omega=omega, determinant_error=assert_small(actual-exact)))
        for n in range(6):
            omega = ((n+.5)*np.pi+(-1)**n*angle)/period
            hb = np.tan(angle/2)
            assert_small((1+hb*hb)*np.cos(omega*period)+2*hb)
            delay = np.diag(np.exp(1j*omega*delays))
            assert_small(np.linalg.det(EYE-np.exp(1j*omega*central)*delay@s@delay@P1))
    return dict(projectors=rows, spectral_determinants=spectral,
                static_reflection_phase_retained=True, total_optical_width_not_central_width_is_period=True)


def canonical_and_clocks():
    canonical = []
    field = np.array([.3+.2j, -.4+.7j])
    for radius in [3., 5., 11.]:
        f = 1+radius**2/4-.8/radius+.2/radius**2
        original = f**(-.25)*field/radius
        radial_norm = radius**2/np.sqrt(f)*np.vdot(original, original).real
        expected = np.vdot(field, field).real/f
        proper_field = f**(-.25)*field
        assert_small(radial_norm-expected)
        assert_small(proper_field-radius*original)
        scaling = .27
        assert_small(np.vdot(field/np.sqrt(scaling), field/np.sqrt(scaling)).real*scaling-np.vdot(field, field).real)
        canonical.append(dict(radius=radius, f=f, radial_measure_error=float(abs(radial_norm-expected))))
    alpha, beta = np.array([.9, 1.3]), np.array([.01, .02])
    clock = lambda i, u: alpha[i]*u+beta[i]*u**3
    clockp = lambda i, u: alpha[i]+3*beta[i]*u**2
    inverse = lambda i, target: brentq(lambda u: clock(i, u)-target, -20., 20., xtol=1e-14)
    rows = []
    delays = DELAYS[0]
    for u in [-.7, .2, 1.1]:
        for i in range(2):
            t = clock(i, u)
            expected = np.sqrt(clockp(i, u))*response(t, delays)[i]
            transported = 0j
            derivative_errors = []
            for j in range(2):
                target = t-delays[i]-delays[j]
                v = inverse(j, target)
                vp = clockp(i, u)/clockp(j, v)
                input_proper = np.sqrt(clockp(j, v))*packet(clock(j, v))[0][j]
                transported += np.sqrt(vp)*physical_source(t-delays[i])[0][i, j]*input_proper
                h = 1e-5
                numeric = (inverse(j, clock(i, u+h)-delays[i]-delays[j])
                           -inverse(j, clock(i, u-h)-delays[i]-delays[j]))/(2*h)
                derivative_errors.append(assert_small(numeric-vp, TOL['source_and_clock']))
            error = assert_small(transported-expected, TOL['source_and_clock'])
            rows.append(dict(u=u, port=i, half_density_error=error, inverse_clock_derivative_errors=derivative_errors))
    return dict(canonical_radial=canonical, unequal_clock_transport=rows)


def delayed_packets():
    norms = [norm2(lambda t, n=n: packet(t)[n], [-WIDTH, WIDTH], 1e-11)[0] for n in range(3)]
    rows = []
    for source_name, source, derivative_bound, second_bound in [
        ('actual_paired_angle', physical_source, .084, .0588+.084**2),
        ('general_unitary_comparison_only', general_source, .124, .0588+.084**2+2*.084*.04+.016+.04**2)]:
        for base in DELAYS:
            for scale in SCALES:
                delays = scale*base
                ell = max(delays)
                limits = [-WIDTH-1, WIDTH+2*ell+1]
                breaks = [sign*WIDTH+delays[i]+delays[j] for sign in [-1, 1] for i in range(2) for j in range(2)]
                limits += breaks
                def zero_error(t):
                    return response(t, delays, source)-source(t)[0]@packet(t)[0]
                def first_error(t):
                    s, sp, _ = source(t)
                    value, derivative, _ = packet(t)
                    length = np.diag(delays)
                    approximation = s@value-length@sp@value-(length@s+s@length)@derivative
                    return response(t, delays, source)-approximation
                actual_norm, _ = norm2(lambda t: response(t, delays, source), limits, 1e-11)
                coarse_norm, _ = norm2(lambda t: response(t, delays, source), limits, 1e-9)
                norm_error = abs(actual_norm-norms[0])
                assert norm_error < TOL['source_and_clock'], norm_error
                assert abs(actual_norm-coarse_norm) < TOL['source_and_clock']
                error0, _ = norm2(zero_error, limits, 1e-11)
                error1, _ = norm2(first_error, limits, 1e-11)
                bound0 = ell*(derivative_bound*norms[0]+2*norms[1])
                bound1 = ell*ell/2*(second_bound*norms[0]+4*derivative_bound*norms[1]+4*norms[2])
                assert error0 <= bound0*(1+TOL['inequality_relative'])
                assert error1 <= bound1*(1+TOL['inequality_relative'])
                characteristic_error = max(assert_small(response(t, delays, source)-direct_characteristics(t, delays, source), TOL['characteristics']) for t in [-.9, -.2, .4, 1., 1.8])
                commutator = float(np.linalg.norm(np.diag(delays)@source(.6)[0]-source(.6)[0]@np.diag(delays)))
                assert commutator > TOL['incorrect_model_minimum']
                rows.append(dict(source=source_name, delays=delays.tolist(), scale=scale,
                                 L2_norm_error=float(norm_error), characteristic_error=characteristic_error,
                                 delay_source_commutator=commutator, zeroth_error=float(error0), zeroth_bound=float(bound0),
                                 first_error=float(error1), first_bound=float(bound1), quadrature_refinement_error=float(abs(actual_norm-coarse_norm))))
    source_noncommutator = float(np.linalg.norm(general_source(.1)[0]@general_source(1.2)[0]-general_source(1.2)[0]@general_source(.1)[0]))
    assert source_noncommutator > 1e-5
    delays = DELAYS[1]
    wrong_time = max(float(np.linalg.norm(response(t, delays)-np.array([sum(physical_source(t)[0][i, j]*packet(t-delays[i]-delays[j])[0][j] for j in range(2)) for i in range(2)]))) for t in np.linspace(-1, 2, 41))
    one_leg = max(float(np.linalg.norm(response(t, delays)-physical_source(t)[0]@np.array([packet(t-delays[j])[0][j] for j in range(2)]))) for t in np.linspace(-1, 2, 41))
    assert wrong_time > TOL['incorrect_model_minimum']
    assert one_leg > TOL['incorrect_model_minimum']
    return dict(packet_H2_norms=list(map(float, norms)), cases=rows,
                general_source_time_commutator=source_noncommutator,
                incorrect_source_time_discrepancy=wrong_time, omitted_leg_discrepancy=one_leg)


def work_and_frequency():
    work = []
    for t in [-.5, .2, .9]:
        a, ap, _ = packet(t)
        s, sp, _ = physical_source(t)
        b, bp = s@a, sp@a+s@ap
        energy = lambda field, derivative: (1j/2)*(np.vdot(field, derivative)-np.vdot(derivative, field))
        actual = energy(b, bp)-energy(a, ap)
        exact = 1j*np.vdot(a, s.conj().T@sp@a)
        work.append(dict(time=t, source_work=float(exact.real), identity_error=assert_small(actual-exact)))
    s = physical_source(.3)[0]
    frequency = []
    for ell in [.1, .01, .001]:
        omega = np.pi/(2*ell)
        returned = np.exp(2j*omega*ell)*s
        distance = float(np.linalg.norm(returned-s, 2))
        assert_small(distance-2)
        frequency.append(dict(ell=ell, frequency=omega, operator_distance=distance))
    return dict(classical_work=work, no_all_frequency_smallness=frequency,
                classical_work_not_unsubtracted_quantum_stress=True)


def main():
    names = ['notes/exterior-lll-delayed-boundary.md', 'checks/exterior_lll_delayed_boundary.py',
             'notes/coupled-fermion-domain.md', 'notes/boundary-vertex-matching.md',
             'notes/returned-switch-capture.md', 'notes/switching-radiation-capture.md']
    result = dict(story='S134.2', date='2026-10-15', status='passed_exact_free_LLL_delayed_cut_response',
                  tolerances=TOL, delays=[x.tolist() for x in DELAYS], scales=SCALES,
                  projector_and_spectrum=projector_and_spectrum(), canonical_and_clocks=canonical_and_clocks(),
                  delayed_packets=delayed_packets(), work_and_frequency=work_and_frequency(),
                  execution=dict(executable=sys.executable, script=str(Path(__file__).resolve()), cwd=str(Path.cwd())),
                  hash_policy='SHA-256 of file bytes after CRLF-to-LF normalization; no other changes',
                  input_sha256={name:hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for name in names},
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  source_review=dict(url='https://arxiv.org/html/2308.00871v2',
                                     scope='Equations 2.2 through 2.14 and 2.23 through 2.25',
                                     observed_equations=['Psi=e^(-sigma/2) psi tensor eta/R',
                                                         'Psi_+^R+h Psi_+^L=0; Psi_-^L+h Psi_-^R=0',
                                                         't=C tau/L_th; f=1+r^2/ell_AdS^2-2GM/r+r_e^2/r^2']),
                  claims=dict(finite_cuts_new_physical_boundary=False, instantaneous_source_used=False,
                              higher_angular_modes_controlled=False, interacting_exterior_eliminated=False,
                              full_quantum_influence_action_derived=False, physical_RN_optical_width_matched=False,
                              moving_geometry_controlled=False, quantum_stress_error_bounded=False,
                              s074_gate_closed=False, opening_achieved=False, novelty_claim=False, lean_invoked=False))
    result['summary'] = dict(projector_cases=len(result['projector_and_spectrum']['projectors']),
                             spectral_cases=len(result['projector_and_spectrum']['spectral_determinants']),
                             delayed_packet_cases=len(result['delayed_packets']['cases']),
                             maximum_packet_norm_error=max(x['L2_norm_error'] for x in result['delayed_packets']['cases']),
                             maximum_characteristic_error=max(x['characteristic_error'] for x in result['delayed_packets']['cases']),
                             maximum_first_Taylor_error_ratio=max(x['first_error']/x['first_bound'] for x in result['delayed_packets']['cases']))
    (ROOT/'evidence/exterior-lll-delayed-boundary.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=result['status'], **result['summary'], execution=result['execution'],
                          hash_policy=result['hash_policy'], input_sha256=result['input_sha256']), indent=2))


if __name__ == '__main__':
    main()
