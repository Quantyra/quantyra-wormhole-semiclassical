"""Full-line reference Ward trace, renormalization and static derivatives."""
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
from scipy.linalg import eigh
from scipy.special import ive, kve

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT/'notes/full-line-reference-stress.md'
OUT = ROOT/'evidence/full-line-reference-stress.json'
MATRIX_TOL = 2e-10
DERIVATIVE_TOL = 5e-5
TAIL_TOL = 2e-6
ENERGY_TOL = 2e-8


def barrier_trace(k, T, m2, z=None):
    r = np.sqrt(k*k+m2)
    z = k if z is None else z
    rho = (r-z)/(r+z)
    x = np.exp(-2*r*T)
    return m2*(T/(2*r)*(1+rho*rho*x)/(1-rho*rho*x)
               +rho*(-np.expm1(-2*r*T))/(2*r*r*(1-rho*rho*x)))


def barrier_logdet(k, T, m2):
    r = np.sqrt(k*k+m2)
    rho = (r-k)/(r+k)
    return (r-k)*T+2*np.log(r+k)-np.log(4*r*k)+np.log1p(-rho*rho*np.exp(-2*r*T))


def barrier_energy(T, m2, mu):
    # Stable local subtraction avoids cancellation of two 1/k terms.
    def integrand(k):
        r, u = np.sqrt(k*k+m2), np.hypot(k, mu)
        rho = m2/(r+k)**2
        x = np.exp(-2*r*T)
        local = T/2*(mu*mu-m2)/(r*u*(r+u))
        image = T/r*rho*rho*x/(1-rho*rho*x)
        image += rho*(-np.expm1(-2*r*T))/(2*r*r*(1-rho*rho*x))
        return m2*(local+image)
    return quad(integrand, 0, np.inf, epsabs=3e-11, epsrel=3e-11, limit=200)[0]


def exterior_data(k, u2, a):
    nu = np.sqrt(.25+u2)
    x = k*a
    kv = kve(nu, x)
    kd = -.5*(kve(nu-1, x)+kve(nu+1, x))
    z = -1/(2*a)-k*kd/kv
    logamp = .5*np.log(2*k*a/np.pi)+np.log(kv)
    return nu, z, logamp


def tailed_logdet(k, T, m2, u2, a):
    _, z, logamp = exterior_data(k, u2, a)
    r = np.sqrt(k*k+m2)
    rho = (r-z)/(r+z)
    # Jost solutions on both infinite tails normalized to exp(-k|s|).
    return 2*logamp+(r-k)*T+2*np.log(r+z)-np.log(4*r*k)+np.log1p(-rho*rho*np.exp(-2*r*T))


def tailed_green_trace(k, T, m2, u2, a):
    nu, z, _ = exterior_data(k, u2, a)
    r = np.sqrt(k*k+m2)
    th = np.tanh(r*T)
    left_log_derivative = (r*th+z)/(1+z*th/r)
    x = k*a
    iv, kv = ive(nu, x), kve(nu, x)
    ider = .5*(ive(nu-1, x)+ive(nu+1, x))
    kder = -.5*(kve(nu-1, x)+kve(nu+1, x))
    reflected = -((.5/a-left_log_derivative)*iv+k*ider)/((.5/a-left_log_derivative)*kv+k*kder)
    def integrand(t):
        y = k*t
        kval = kve(nu, y)
        return u2/t*(ive(nu, y)*kval+reflected*np.exp(-2*k*(t-a))*kval*kval)
    tail = quad(integrand, a, np.inf, epsabs=3e-11, epsrel=3e-11, limit=200)[0]
    return barrier_trace(k, T, m2, z)+2*tail


def matrix_checks():
    rng = np.random.default_rng(20261011)
    rows = []
    for n in [4, 7]:
        for k in [.2, 1., 4.]:
            for replicate in range(2):
                derivative = rng.normal(size=(n, n))
                L = np.vstack([derivative, -1j*k*np.eye(n)])
                A = (L.conj().T@L).real
                eigen, basis = eigh(A)
                invsqrt = (basis/np.sqrt(eigen))@basis.T
                v = rng.uniform(.1, 3., n)
                f = rng.uniform(-.8, .8, n)
                h = rng.uniform(-.4, .4, n)
                V = np.diag(v)
                V1, V2 = np.diag(v*f), np.diag(v*(f*f+h))
                B, B1, B2 = [invsqrt@X@invsqrt for X in [V, V1, V2]]
                C = np.linalg.inv(np.eye(n)+B)
                R = np.linalg.inv(A+V)
                U = L@invsqrt
                P = U@U.conj().T
                D = np.eye(2*n)-L@R@L.conj().T
                F = float(np.trace(V@R))
                ward_error = abs(np.trace(D@P)-F)
                assert ward_error < MATRIX_TOL
                assert np.linalg.norm(P@P-P) < MATRIX_TOL
                assert np.linalg.norm(B@B1-B1@B) > 1e-5
                first = np.trace(C@C@B1)
                second = np.trace(C@C@B2)-2*np.trace(C@B1@C@C@B1)
                independent_first = np.trace(V1@R-V@R@V1@R)
                independent_second = np.trace(V2@R-2*V1@R@V1@R-V@R@V2@R+2*V@R@V1@R@V1@R)
                assert abs(first-independent_first) < MATRIX_TOL
                assert abs(second-independent_second) < MATRIX_TOL
                c1, c2 = np.max(np.abs(f)), np.max(np.abs(f*f+h))
                assert abs(first) <= c1*F+MATRIX_TOL
                assert abs(second) <= (c2+2*c1*c1)*F+MATRIX_TOL
                steps = []
                for step in [.002, .001, .0005]:
                    vals = []
                    for phi in [-step, 0, step]:
                        vv = np.diag(v*np.exp(phi*f+phi*phi*h/2))
                        vals.append(float(np.trace(vv@np.linalg.inv(A+vv))))
                    fd1 = (vals[2]-vals[0])/(2*step)
                    fd2 = (vals[2]-2*vals[1]+vals[0])/(step*step)
                    steps.append(dict(step=step, first=fd1, second=fd2))
                assert abs(steps[-1]['first']-first) < DERIVATIVE_TOL
                assert abs(steps[-1]['second']-second) < DERIVATIVE_TOL
                rows.append(dict(n=n, k=k, replicate=replicate, ward_error=float(ward_error),
                                 F=F, first=float(first), second=float(second), steps=steps))
    return rows


def run():
    warnings.simplefilter('error', IntegrationWarning)
    matrices = matrix_checks()
    barriers = []
    for m2 in [.3, 2., 7.]:
        for k in [.02, .2, 1., 4., 12.]:
            target = barrier_trace(k, 1.7, m2)
            # Holomorphic coupling derivative, independent Jost construction.
            jost = barrier_logdet(k, 1.7, m2*(1+1j*1e-20)).imag/1e-20
            assert abs(jost-target) < 2e-10
            barriers.append(dict(m2=m2, k=k, green=float(target), jost_derivative=float(jost), error=float(abs(jost-target))))
    tails = []
    for u2 in [.35, 2.]:
        for k in [.03, .3, 2., 8.]:
            T, m2, a = 1.5, 1.3**2, .8
            target = tailed_green_trace(k, T, m2, u2, a)
            vals = []
            for step in [.002, .001, .0005]:
                lo = tailed_logdet(k, T, m2*(1-step), u2*(1-step), a)
                hi = tailed_logdet(k, T, m2*(1+step), u2*(1+step), a)
                vals.append(float((hi-lo)/(2*step)))
            assert abs(vals[-1]-target) < TAIL_TOL
            assert abs(vals[-1]-vals[-2]) < TAIL_TOL
            tails.append(dict(u2=u2, k=k, full_green_trace=float(target), jost_refinements=vals, error=abs(vals[-1]-target)))
    energy = []
    for T, m2, mu in [(1.7, .3, .7), (2.3, 2., 1.), (.9, 7., .2)]:
        base = barrier_energy(T, m2, mu)
        I = m2*T
        vals = []
        for phi in [-.2, -.05, .1, .3]:
            actual = barrier_energy(T*np.exp(-phi), m2*np.exp(2*phi), mu)
            predicted = np.exp(phi)*(base-I*phi/2)
            assert abs(actual-predicted) < ENERGY_TOL
            vals.append(dict(phi=phi, actual=actual, predicted=float(predicted), error=abs(actual-predicted)))
        jets = []
        for step in [.008, .004, .002]:
            lo = barrier_energy(T*np.exp(step), m2*np.exp(-2*step), mu)
            hi = barrier_energy(T*np.exp(-step), m2*np.exp(2*step), mu)
            jets.append(dict(step=step, first=(hi-lo)/(2*step), second=(hi-2*base+lo)/(step*step)))
        expected = [base-I/2, base-I]
        assert max(abs(jets[-1]['first']-expected[0]), abs(jets[-1]['second']-expected[1])) < DERIVATIVE_TOL
        energy.append(dict(T=T, m2=m2, mu=mu, energy=base, dilation_checks=vals, derivative_refinements=jets, expected_derivatives=expected))
    angular = []
    for q in [2, 3, 8, 64, 256]:
        ell = np.arange(1, q, dtype=float)
        beta = np.cumprod((q-ell)/(q+ell))
        kap = ell*(ell+1)
        tau = (2*ell+1)*q*beta/(4*np.pi)
        assert np.all(beta <= np.exp(-kap/q)*(1+1e-13))
        s0, s1, shalf = np.sum(tau), np.sum(tau/kap), np.sum(tau/np.sqrt(kap))
        hq = np.sum(1/np.arange(1, q+1, dtype=float))
        assert abs(s0-q*(q-1)/(4*np.pi)) < 1e-10*max(1,s0)
        assert abs(s1-q*(hq-1)/(4*np.pi)) < 1e-10*max(1,s1)
        bound = 3*q**1.5/(8*np.sqrt(np.pi))
        assert shalf <= bound
        angular.append(dict(q=q, sum_tau=float(s0), inverse_mass_sum=float(s1), inverse_sqrt_mass_sum=float(shalf), bound=bound))
    report = dict(date='2026-10-11', status='passed_full_line_reference_stress_checks',
        scope='Complete order-g^2 unwrapped full-line current reference in the prescribed spherical photon model',
        note_sha256=hashlib.sha256(NOTE.read_bytes()).hexdigest(), script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        tolerances=dict(matrix_absolute=MATRIX_TOL, static_derivative_absolute=DERIVATIVE_TOL, full_tail_absolute=TAIL_TOL, energy_absolute=ENERGY_TOL),
        matrix_cases=matrices, barrier_cases=barriers, inverse_square_tail_cases=tails, energy_cases=energy, angular_cases=angular,
        artificial_mouth_subtraction_used=False, exterior_domain_truncated=False,
        physical_two_mouth_matching_proved=False, complete_finite_return_diagram_controlled=False,
        higher_order_remainder_controlled=False, physical_damping_proved=False, novelty_established=False, lean_invocations=0)
    OUT.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='')
    print(json.dumps(dict(status=report['status'], matrix_cases=len(matrices), barrier_cases=len(barriers), tail_cases=len(tails),
        max_matrix_error=max(v['ward_error'] for v in matrices), max_barrier_error=max(v['error'] for v in barriers),
        max_tail_error=max(v['error'] for v in tails), energies=[v['energy'] for v in energy], angular_cases=len(angular)),indent=2))


if __name__ == '__main__':
    run()
