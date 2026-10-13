"""S072 common-clock charged midpoint matching and boundary conventions.

Checks compare Matsubara sums with a continuous spectral integral and
the exact barrier transfer with a separately integrated Riccati equation.
These are projected-model checks, not an interacting wormhole simulation.
"""
from pathlib import Path
import hashlib
import json
import warnings

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import IntegrationWarning, quad, solve_ivp

warnings.simplefilter("error", IntegrationWarning)
ROOT = Path(__file__).resolve().parents[1]


def boundary_conventions():
    J = np.array([[0., -1.], [1., 0.]])
    unit = np.eye(2)
    rows = []
    for hb in (0., 1.e-4, .01, .2, .7, .99):
        theta = 2*np.arctan(hb)
        scattering = np.linalg.solve(unit-hb*J, unit+hb*J)
        rotation = np.cos(theta)*unit+np.sin(theta)*J
        error = float(np.max(abs(scattering-rotation)))
        lam_b = .5*np.arctan(2*hb/abs(1-hb*hb))
        branch_error = abs(theta-2*lam_b)
        assert error < 5.e-15 and branch_error < 5.e-15
        assert np.max(abs(scattering.T@scattering-unit)) < 5.e-15
        rows.append(dict(h_boundary=hb, theta=theta,
                         matrix_error=error, branch_error=branch_error))

    eta = 1/2000
    tau = .25
    end = np.pi+tau

    def S(t):
        v = np.clip(t/tau, 0., 1.)
        return 10*v**3-15*v**4+6*v**5

    def primitive(t):
        if t <= 0:
            return 0.
        if t >= tau:
            return t-tau/2
        v = t/tau
        return tau*(2.5*v**4-3*v**5+v**6)

    def theta(t):
        return 4*eta*(primitive(t)-primitive(t-np.pi))

    def rate(t, contact):
        answer = 4*eta*(S(t)-S(t-np.pi))
        return answer/np.cos(theta(t)/2)**2 if contact else answer

    fourier = []
    for contact in (False, True):
        value = []
        for weight in ("cos", "sin"):
            value.append(sum(quad(lambda t: rate(t, contact), a, b,
                                  weight=weight, wvar=2, epsabs=2.e-15,
                                  epsrel=2.e-12)[0]
                             for a, b in zip((0, tau, np.pi),
                                             (tau, np.pi, end))))
        fourier.append(complex(*value))
    assert abs(fourier[0]) < 2.e-16
    assert abs(fourier[1]) > 1.e-10
    final_angle = 4*np.pi*eta
    return dict(branches=rows, shaped_switch=dict(
        eta=eta, tau=tau, final_angle=final_angle,
        final_boundary_parameter=np.tan(final_angle/2),
        contact_relative_correction=2*np.tan(final_angle/2)/final_angle-1,
        angle_rate_fourier_magnitude=abs(fourier[0]),
        contact_rate_fourier_magnitude=abs(fourier[1])))


def constant_delta(s, mass):
    return mass*mass/(np.sqrt(s*s+mass*mass)+s)


def barrier_delta(s, mass, length):
    kval = np.sqrt(s*s+mass*mass)
    tangent = np.tanh(kval*length)
    return mass*mass*tangent/(kval+s*tangent)


def midpoint_sum(delta, beta, q, K0, vmax, cutoff=200001, epsilon=0):
    assert cutoff % 2 == 1
    n = np.arange(1, cutoff+1, 2, dtype=float)
    freq = 2*np.pi*n/beta
    lam = q/np.pi
    answer = 4/(beta*lam)*np.sum(delta(freq)*np.exp(-epsilon*n)/freq**2)
    tail = min(beta*K0/(np.pi**2*lam*cutoff),
               beta*beta*vmax/(8*np.pi**3*lam*cutoff**2))
    return float(answer), float(tail)


def spectral_constant(mass, beta, q):
    """Independent real-frequency physical-minus-reference covariance.

    Set w=m cosh(v) in the high-frequency integral. Truncation at v=30
    has absolute error below 2e-26/(pi lambda); quadrature is independent
    of the Matsubara sum.
    """
    lam = q/np.pi
    scale = beta*mass/4
    low, low_err = quad(lambda y: np.tanh(scale*y)/y if y else scale,
                        0, 1, epsabs=2.e-13, epsrel=2.e-13)

    def high(v):
        tangent = np.tanh(v)
        # Stable tanh(v)-1, including v beyond floating-point saturation.
        deficit = -2*np.exp(-2*v)/(1+np.exp(-2*v))
        return deficit*tangent*np.tanh(scale*np.cosh(v))

    upper, upper_err = quad(high, 0, 30, epsabs=2.e-13, epsrel=2.e-13)
    value = beta*mass/(8*lam)+(upper-low)/(np.pi*lam)
    error = (low_err+upper_err+2.e-26)/(np.pi*lam)
    return float(value), float(error)


def spectral_checks():
    rows = []
    for beta in (.25, 1.7, 6.):
        for mass in (.05, .8, 3.):
            q = 7
            discrete, tail = midpoint_sum(lambda s: constant_delta(s, mass),
                                          beta, q, mass, mass*mass)
            continuous, quadrature = spectral_constant(mass, beta, q)
            bound = beta*mass/(8*(q/np.pi))
            difference = abs(discrete-continuous)
            assert difference <= tail+quadrature+3.e-14
            assert -3.e-14 <= discrete <= continuous+quadrature+3.e-14
            assert 0 <= continuous <= bound
            rows.append(dict(beta=beta, mass=mass, q=q,
                             matsubara_partial=discrete,
                             spectral_integral=continuous,
                             absolute_difference=difference,
                             matsubara_tail_bound=tail,
                             quadrature_error_estimate=quadrature,
                             theorem_upper_bound=bound))
    return rows


def barrier_checks():
    rows = []
    for mass, length in ((.2, .3), (.8, 1.4), (3., .7)):
        K0 = mass*np.tanh(mass*length)
        frequencies = np.r_[0., np.geomspace(1.e-6, 1.e4, 100)]
        differences = barrier_delta(frequencies, mass, length)
        assert np.all(differences >= 0)
        assert np.max(differences) <= K0*(1+2.e-14)
        dynamic = []
        for s in (0., .01, .5, 2., 7.):
            reference = s+barrier_delta(s, mass, length)
            solution = solve_ivp(lambda x, y: [mass*mass+s*s-y[0]**2],
                                 (-length, 0), [s], method="DOP853",
                                 rtol=3.e-12, atol=3.e-13,
                                 max_step=length/60)
            assert solution.success
            error = abs(solution.y[0, -1]-reference)
            assert error < 3.e-11
            dynamic.append(dict(s=s, exact=reference, riccati_error=error))
        beta, q = 1.7, 7
        value, tail = midpoint_sum(lambda s: barrier_delta(s, mass, length),
                                   beta, q, K0, mass*mass)
        assert 0 <= value and value+tail <= beta*K0/(8*(q/np.pi))
        rows.append(dict(mass=mass, length=length, K0=K0,
                         riccati=dynamic, matching_partial=value,
                         tail_bound=tail))
    return rows


def regulator_and_clock_checks():
    beta, mass, q = 1.3, .9, 7
    lam = q/np.pi
    kappa = 2*np.pi/beta
    rows = []
    last_difference = -1.
    for epsilon in (.4, .1, .025, .00625):
        n = np.arange(1, 200002, 2, dtype=float)
        s = kappa*n
        cutoff = np.exp(-epsilon*n)
        K = s+constant_delta(s, mass)
        full = 4/(beta*lam)*np.sum(K/s**2*cutoff)
        ref_sum = 2/q*np.sum(cutoff/n)
        ref_closed = np.log(1/np.tanh(epsilon/2))/q
        difference, _ = midpoint_sum(lambda v: constant_delta(v, mass),
                                      beta, q, mass, mass*mass,
                                      epsilon=epsilon)
        assert abs(ref_sum-ref_closed) < 2.e-15
        assert abs(full-ref_closed-difference) < 2.e-15
        assert difference > last_difference
        last_difference = difference
        # beta_t=2pi, G_single,t(n)=kappa G_single,U(kappa n).
        global_midpoint = 4/(2*np.pi)*np.sum(kappa*K/(lam*s*s)*cutoff)
        assert abs(global_midpoint-full) < 2.e-15
        rows.append(dict(epsilon=epsilon, physical_cutoff=epsilon/kappa,
                         W_full=full, W_reference=ref_closed,
                         finite_difference=difference,
                         wrong_single_boundary_difference=full/2-ref_closed,
                         clock_error=abs(global_midpoint-full)))
    wrong = [r["wrong_single_boundary_difference"] for r in rows]
    assert all(a > b for a, b in zip(wrong, wrong[1:]))
    return rows


def candidate_checks():
    mp.mp.dps = 70
    rows = []
    for kval in (2000, 4000, 10000):
        k = mp.mpf(kval)
        X = (mp.sqrt(1+12*mp.pi*k**6)-1)/6
        G = mp.mpf(1)
        L = k**5
        b = L*mp.sqrt(X)
        C = 1+6*X
        N = 2*k**4
        eta = 1/k
        bare_inertia = 4*mp.pi*b*k**8
        for z0 in (mp.mpf(1)/4, mp.mpf(7)/25):
            phi_r = b**3/(G*C)
            beta = 2*mp.pi*phi_r/(N*eta*z0)
            direct = beta/(2*bare_inertia)
            reduced = X/(8*k*z0*C)
            assert abs(direct/reduced-1) < mp.mpf("1e-65")
            upper = 1/(192*k*z0)
            assert direct/4 < upper <= 1/(48*k)
            rows.append(dict(k=kval, z0=float(z0),
                             x_bare=float(direct),
                             midpoint_exponent_upper=float(direct/4),
                             uniform_exponent_upper=float(1/(48*k)),
                             matched_amplitude_lower=float(mp.exp(-direct/4)),
                             substitution_relative_error=float(abs(direct/reduced-1))))
    return rows


def main():
    result = dict(story="S072", date="2026-10-13", status="PASS",
                  scope="Free boundary map and initial projected electric Gaussian midpoint normalization",
                  versions=dict(numpy=np.__version__, scipy=scipy.__version__,
                                mpmath=mp.__version__),
                  boundary=boundary_conventions(),
                  constant_potential_spectral=spectral_checks(),
                  finite_barrier=barrier_checks(),
                  regulator_and_clock=regulator_and_clock_checks(),
                  candidate=candidate_checks(),
                  limitations=["No interacting switching-stress or evolved-state bound",
                               "Compact charge-spin weights and angular/coset matching unresolved",
                               "Cutoff removal is within the projected massless Gaussian model",
                               "No physical preparation, full stability or finite-energy transport proof"])
    result["artifacts"] = [dict(path=path, sha256=hashlib.sha256(
        (ROOT/path).read_bytes()).hexdigest()) for path in
        ("notes/boundary-vertex-matching.md", "checks/boundary_vertex_matching.py")]
    output = ROOT/"evidence/boundary-vertex-matching.json"
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps(dict(status="PASS", spectral_cases=len(result["constant_potential_spectral"]),
                          riccati_cases=sum(len(r["riccati"]) for r in result["finite_barrier"]),
                          candidate_cases=len(result["candidate"]),
                          maximum_spectral_difference=max(r["absolute_difference"] for r in
                                                          result["constant_potential_spectral"]),
                          maximum_riccati_error=max(t["riccati_error"] for r in
                                                   result["finite_barrier"] for t in r["riccati"]),
                          k2000_exponent_bound=result["candidate"][0]["midpoint_exponent_upper"],
                          output=str(output)), indent=2))


if __name__ == "__main__":
    main()
