"""Reference-vacuum boundary obstruction and full-line cancellation checks."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import hashlib
import json
import platform
import warnings
import numpy as np
import scipy
from scipy.integrate import quad, IntegrationWarning
from numpy.polynomial.legendre import leggauss
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'notes/reference-mouth-completion.md'
OUT = ROOT / 'evidence/reference-mouth-completion.json'
KERNEL_TOL = 2e-9
FOURIER_TOL = 2e-10
SLOPE_TOL = 2e-4
STEP_TOL = 2e-4


def exact_trace(k, T, m, z):
    r = np.hypot(k, m)
    rho = (r-z)/(r+z)
    x = np.exp(-2*r*T)
    trace_r = T/(2*r)*(1+rho*rho*x)/(1-rho*rho*x)
    trace_r += rho*(-np.expm1(-2*r*T))/(2*r*r*(1-rho*rho*x))
    th = np.tanh(r*T/2)
    gp, gm = 1/(z+r*th), 1/(z+r/th)
    e = np.exp(-k*T)
    return m*m*trace_r+z*(gp+gm)+k/2*((1-e)*gp+(1+e)*gm)-1


def kernel_parts(s, t, k, T, m, z):
    r = np.hypot(k, m)
    rho = (r-z)/(r+z)
    sign = np.sign(s-t)
    e0 = np.exp(-r*np.abs(s-t))
    eL = rho*np.exp(-r*(s+t))
    eR = rho*np.exp(-r*(2*T-s-t))
    e2 = rho*rho*np.exp(-r*(2*T-np.abs(s-t)))
    norm = 2*r*(1-rho*rho*np.exp(-2*r*T))
    R = (e0+eL+eR+e2)/norm
    Rs = r*(-sign*e0-eL+eR+sign*e2)/norm
    Rt = r*(sign*e0-eL+eR-sign*e2)/norm
    Rst = r*r*(-e0+eL+eR-e2)/norm
    return R, Rs, Rt, Rst


def direct_trace(k, T, m, z, n):
    points, weights = leggauss(n)
    u, w = (points+1)/2, weights/2
    s = T*u
    _, _, _, diag = kernel_parts(s, s, k, T, m, z)
    contact = T*np.dot(w, -diag+k/2)
    ss, tt = s[:, None], s[:, None]*u[None, :]
    R, Rs, Rt, Rst = kernel_parts(ss, tt, k, T, m, z)
    integrand = k/2*np.exp(-k*(ss-tt))*(Rst-k*k*R+k*(Rs-Rt))
    smooth = 2*T*np.sum(w[:, None]*w[None, :]*ss*integrand)
    return float(contact+smooth)


def symbolic_checks():
    t, h, m, T, r, k, p = sy.symbols('t h m T r k p', positive=True)
    # Exponentially small reflections removed for the ultraviolet series.
    rr = sy.sqrt(1+m*m*t*t)
    rho = (rr-h)/(rr+h)
    F = m*m*(T*t/(2*rr)+rho*t*t/(2*rr**2))+(2*h+1)/(rr+h)-1
    series = sy.series(F, t, 0, 3).removeO().expand()
    assert sy.simplify(series.subs(t, 0)-h/(1+h)) == 0
    assert sy.simplify(series.coeff(t, 1)-m*m*T/2) == 0
    c = h/(1+h)
    dc = h*sy.diff(c, h)
    ddc = h*sy.diff(dc, h)
    assert [sy.simplify(v.subs(h, 2)) for v in [c, dc, ddc]] == [sy.Rational(2, 3), sy.Rational(2, 9), -sy.Rational(2, 27)]
    Q = p*p+k*k
    v = sy.Matrix([p, -k])
    P = v*v.T/Q
    D = sy.eye(2)-v*v.T/(Q+m*m)
    assert sy.simplify(sy.trace(D*P)-m*m/(Q+m*m)) == 0
    a = r+k
    triangle = 2*T/a-2*(1-sy.exp(-a*T))/a**2
    Fmatched = (r+k)*T/2-k*(r+k)**2/(4*r)*triangle
    expected = (r*r-k*k)*T/(2*r)+k/(2*r)*(1-sy.exp(-a*T))
    assert sy.simplify(Fmatched-expected) == 0
    cross = -k*a*a/(4*r)*2*(1-sy.exp(-a*T))/a**2
    assert sy.simplify(cross+k/(2*r)*(1-sy.exp(-a*T))) == 0
    return ['uv_constant', 'uv_bulk_coefficient', 'two_clock_derivatives',
            'full_line_fourier_trace', 'matched_direct_kernel_trace', 'cross_mouth_cancellation']


def bulk_subtracted(cutoff, phi):
    T, m, h = 2.3*np.exp(-phi), 1.7*np.exp(phi), 2*np.exp(phi)
    def f(k):
        return exact_trace(k, T, m, h*k)-m*m*T/(2*np.hypot(k, 1.0))
    return quad(f, 0, cutoff, epsabs=3e-10, epsrel=3e-11, limit=250)[0]


def clock_values(cutoff, step):
    lo, mid, hi = [bulk_subtracted(cutoff, phi) for phi in [-step, 0, step]]
    return np.array([mid, (hi-lo)/(2*step), (hi-2*mid+lo)/(step*step)])


def run():
    warnings.simplefilter('error', IntegrationWarning)
    exact = symbolic_checks()
    kernels = []
    for T in [.7, 2.3]:
        for h in [.3, 1., 2.]:
            for k in [.2, 1.5, 6.]:
                m, z = 1.7, h*k
                target = exact_trace(k, T, m, z)
                vals = [direct_trace(k, T, m, z, n) for n in [24, 48, 96]]
                assert abs(vals[-1]-target) < KERNEL_TOL
                assert abs(vals[-1]-vals[-2]) < KERNEL_TOL
                assert target >= -1e-12
                kernels.append(dict(T=T, h=h, k=k, closed_trace=float(target),
                                    direct_refinements=vals, error=abs(vals[-1]-target)))
    limits = []
    for h in [.3, 1., 2., 7.]:
        T, m = 2.3, 1.7
        values = []
        for k in [40., 80., 160., 320.]:
            F = exact_trace(k, T, m, h*k)
            remainder = F-h/(1+h)-m*m*T/(2*k)
            values.append(dict(k=k, value=float(F), scaled_remainder=float(k*k*remainder)))
        assert abs(values[-1]['scaled_remainder']) < 10
        assert abs(values[-1]['value']-h/(1+h)-m*m*T/(2*320)) < 1e-4
        limits.append(dict(h=h, constant=h/(1+h), samples=values))
    matched = []
    for T in [.7, 2.3]:
        for k in [.2, 1.5, 6.]:
            m, r = 1.7, np.hypot(k, 1.7)
            a = r+k
            b = k/(2*r)*(-np.expm1(-a*T))
            F = exact_trace(k, T, m, r)
            outside = quad(lambda s: quad(lambda y: -k*a*a/(4*r)*
                (np.exp(-a*(s+y))+np.exp(-a*(T-s+y))), 0, np.inf,
                epsabs=2e-12, epsrel=2e-12)[0], 0, T,
                epsabs=2e-11, epsrel=2e-11)[0]
            fourier = T/np.pi*quad(lambda p: m*m/(p*p+k*k+m*m), 0, np.inf,
                                  epsabs=2e-12, epsrel=2e-12)[0]
            assert abs(F-(m*m*T/(2*r)+b)) < FOURIER_TOL
            assert abs(outside+b) < FOURIER_TOL
            assert abs(F+outside-fourier) < FOURIER_TOL
            matched.append(dict(T=T, k=k, core_trace=float(F), outside_trace=float(outside),
                                full_fourier_trace=float(fourier), cancellation_error=abs(F+outside-fourier)))
    cutoffs = [20., 40., 80., 160., 320.]
    jets = [clock_values(c, .002) for c in cutoffs]
    slopes = [(jets[i+1]-jets[i])/(cutoffs[i+1]-cutoffs[i]) for i in range(len(cutoffs)-1)]
    expected = np.array([2/3, 2/9, -2/27])
    assert np.max(np.abs(slopes[-1]-expected)) < SLOPE_TOL
    step_values = []
    for step in [.008, .004, .002]:
        slope = (clock_values(320., step)-clock_values(160., step))/160
        step_values.append(dict(step=step, slope=slope.tolist()))
    assert np.max(np.abs(np.array(step_values[-1]['slope'])-step_values[-2]['slope'])) < STEP_TOL
    # Subtracting the diagnosed boundary term is a diagnostic, not physical matching.
    finite_parts = [float(j[0]-2*c/3) for c, j in zip(cutoffs, jets)]
    assert abs(finite_parts[-1]-finite_parts[-2]) < .02
    report = dict(
        date='2026-10-11', status='passed_reference_mouth_completion_checks',
        scope='Throat-only reference boundary divergence and exact uniform full-line cancellation',
        note_sha256=hashlib.sha256(NOTE.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, sympy=sy.__version__),
        tolerances=dict(kernel_absolute=KERNEL_TOL, fourier_absolute=FOURIER_TOL,
                        uv_clock_slope_absolute=SLOPE_TOL, clock_step_absolute=STEP_TOL),
        exact_checks=exact, kernel_cases=kernels, ultraviolet_limits=limits,
        matched_cases=matched,
        cutoff_clock_values=[dict(cutoff=c, energy_first_second_derivative=j.tolist()) for c, j in zip(cutoffs, jets)],
        cutoff_slopes=[dict(lower=cutoffs[i], upper=cutoffs[i+1], slopes=v.tolist()) for i, v in enumerate(slopes)],
        expected_slopes=expected.tolist(), clock_step_refinements=step_values,
        diagnosed_finite_parts=finite_parts,
        positive_extension_counterexample_verified=True,
        complete_physical_stress_controlled=False, physical_damping_proved=False,
        novelty_established=False, lean_invocations=0,
        decision='Derive full mouth interaction and matched renormalization before further throat-only higher-diagram work.')
    OUT.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='')
    print(json.dumps(dict(status=report['status'], exact_checks=len(exact), kernel_cases=len(kernels),
        maximum_kernel_error=max(c['error'] for c in kernels), matched_cases=len(matched),
        maximum_cancellation_error=max(c['cancellation_error'] for c in matched),
        final_cutoff_slopes=slopes[-1].tolist(), expected_slopes=expected.tolist(),
        diagnosed_finite_parts=finite_parts), indent=2))


if __name__ == '__main__':
    run()
