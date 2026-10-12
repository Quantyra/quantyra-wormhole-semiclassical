"""Check the specified untwisted gluing models; no full wormhole claim."""
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
import sympy as s
from scipy.integrate import quad, solve_bvp
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
KAPPA = 1 + np.arcsinh(np.sqrt(3)) / (2 * np.sqrt(3))


def exact_checks():
    checks = []

    def zero(name, expression):
        assert s.simplify(expression) == 0, name
        checks.append(name)

    a, t, Q = s.symbols('a t Q', positive=True)
    R = (a + Q) / 2
    u2 = (Q**2 - a**2) / (4*t)
    F = R + 4*t*R**2/Q**2 - 4*t*R*u2/Q**2
    simple = R*(a**2 + 4*t*R)/Q**2
    zero('arc integrand identity', F-simple)
    z = a/Q
    log_derivative = s.diff(simple, Q)*Q*(1-z*z)/simple
    negative = -(1-z)*(2*z+a*a/(a*a+2*t*(a+Q)))
    zero('global arc monotonicity identity', log_derivative-negative)
    zero('axis endpoint integrand', simple.subs({a: 1-t, Q: 1+t})-1)

    A, B, phi = s.symbols('A B phi', positive=True)
    quadratic = A*phi**2+B*(phi+s.pi)**2
    zero('weighted phase minimizer', s.diff(quadratic, phi).subs(phi, -s.pi*B/(A+B)))
    zero('weighted minimum', quadratic.subs(phi, -s.pi*B/(A+B))-s.pi**2*A*B/(A+B))
    energy = -s.pi*(A+B)/6+quadratic/(2*s.pi)
    zero('minimum energy decomposition', energy.subs(phi, -s.pi*B/(A+B))
         +s.pi*(A+B)/24+s.pi*(A-B)**2/(8*(A+B)))
    q = s.symbols('q', positive=True)
    rho = energy.subs(phi, -s.pi*B/(A+B))+q*s.pi/24
    zero('optimized support coefficient', -8*rho/(s.pi*q)
         -((A+B)/(3*q)-s.Rational(1, 3)+(A-B)**2/(q*(A+B))))

    r, x, J, e, k = s.symbols('r x J e k', real=True)
    f = (r-1)**2/r**2
    u0 = (r-1)**2*(r+2)/r
    u1 = (2*r-1)*(5*r*r-5*r+2)/(30*r*(r-1)**3)
    u = e*u0+u1-2*J/r
    zero('RN inhomogeneous equation', s.diff(f*s.diff(u, r), r)-(2/r**2+4/r**4)*u-12*J/r**4)
    zero('RN Wronskian', f*(u0*s.diff(u1, r)-s.diff(u0, r)*u1)+1)
    H = x*(1+x*x)
    V = k*H-2*J*x
    w = -2*k*H-2*J*x
    zero('odd throat Maxwell', (1+x*x)*s.diff(V, x, 2)-2*(V-w))
    zero('odd throat rotation', (1+x*x)*s.diff(w, x, 2)+4*(V-w))
    zero('odd throat conserved charge', s.diff(w+2*V, x)+6*J)
    Oa = s.diag(-1, -1, -1)
    Fmap = s.diag(1, 1, -1)
    assert Oa.T*Fmap == s.diag(-1, -1, 1)
    assert Fmap*s.Matrix([0, 0, 1]) == s.Matrix([0, 0, -1])
    assert Fmap*s.Matrix([1, 0, 0]) == s.Matrix([1, 0, 0])
    checks.extend(['relative half turn', 'preserved magnetic equality', 'reversed electric parity'])
    return checks


def g_of_u(u):
    def integrand(eta):
        a = 1-eta*eta
        Q = np.hypot(a, 2*u*eta)
        R = (a+Q)/2
        return np.sqrt(R*(a*a+4*eta*eta*R))/Q
    points = [1-v*u for v in [1, 10, 100, 1000] if 0 < 1-v*u < 1]
    return quad(integrand, 0, 1, points=points, epsabs=2e-11, epsrel=2e-11)[0]


def geometry_checks():
    rows = []
    for u in [0.0001, 0.003, 0.03, 0.2, 0.6, 0.95]:
        # Independent Euclidean arclength of the prolate-parametrized curve.
        def arc(eta):
            a = 1-eta*eta
            Q = np.hypot(a, 2*u*eta)
            mu = (a+Q)/(2*u)
            mup = 2*eta*(u-mu)/Q
            rh = np.sqrt((mu*mu-1)*(1-eta*eta))/2
            rhp = (mu*mup*(1-eta*eta)-eta*(mu*mu-1))/(4*rh)
            zp = (mu+eta*mup)/2
            return 2*np.hypot(rhp, zp)
        direct = quad(arc, 0, 1, epsabs=1e-8, epsrel=2e-10)[0]
        g = g_of_u(u)
        rel = abs(direct-g/u)/direct
        assert rel < 2e-9
        assert 1 <= g <= KAPPA
        rows.append(dict(u=u, g=g, arclength=direct, relative_error=rel))
    grid = np.geomspace(1e-6, 1, 181)
    gs = np.array([g_of_u(u) for u in grid])
    assert np.all(np.diff(gs) < 0)
    assert np.all(np.diff(gs/grid) < 0)
    return dict(independent_arc=rows, monotonicity_grid=len(grid), small_u_g=float(gs[0]), kappa=KAPPA)


def wrapped(phi, theta):
    return (phi+theta+np.pi) % (2*np.pi)-np.pi


def minimize_phase(theta, weights):
    # Enumerate every smooth branch; this is independent of the two-class formula.
    points = np.unique(np.r_[0., 2*np.pi, (np.pi-theta) % (2*np.pi)])
    best = float('inf')
    for lo, hi in zip(points[:-1], points[1:]):
        if hi-lo < 1e-13:
            continue
        sol = minimize_scalar(lambda phi: np.dot(weights, wrapped(phi, theta)**2),
                              bounds=(lo, hi), method='bounded', options={'xatol': 1e-12})
        best = min(best, sol.fun, np.dot(weights, wrapped(lo, theta)**2),
                   np.dot(weights, wrapped(hi, theta)**2))
    return best


def phase_checks():
    rows = []
    for q in range(1, 25):
        m = np.arange(q)-(q-1)/2
        theta = (np.pi*m) % (2*np.pi)
        minimum = minimize_phase(theta, np.ones(q))
        expected = np.pi**2*(q//2)*((q+1)//2)/q
        assert abs(minimum-expected) < 2e-10
        chi = 1-4*minimum/(q*np.pi**2)
        assert abs(chi-(q % 2)/q**2) < 1e-11
        rows.append(dict(q=q, minimum_phase_square=minimum, best_equal_length_support=chi))
    # Check regularized determinant, independently of the spectral polynomial.
    winding = []
    for L, angle in [(0.7, -2.3), (1.1, 0.3), (2.7, 1.7)]:
        integral = -2/np.pi*quad(lambda x: np.log(abs(1+np.exp(-x*L+1j*angle))),
                                0, 80/L, epsabs=2e-12)[0]
        target = -np.pi/(6*L)+angle**2/(2*np.pi*L)
        assert abs(integral-target) < 3e-12
        winding.append(dict(L=L, phase=angle, absolute_error=abs(integral-target)))
    delayed = []
    for q in [7, 8, 31, 32, 127, 128, 511, 512]:
        ratio = np.log(q)**-2/np.pi  # d/T in the retained family.
        us = (np.arange(q)+0.5)/q
        delays = ratio*np.array([g_of_u(u)/u for u in us])
        lengths = 1+delays
        theta = np.pi*(np.arange(q) % 2)
        b = lengths**-2
        B, Db = b.sum(), b[::2].sum()-b[1::2].sum()
        chi = B/(3*q)-1/3+Db*Db/(q*B)
        upper = -(1-(1+ratio)**-2)/3+3*(1+2*KAPPA*ratio)**2/q**2
        numeric = 4*B/(3*q)-1/3-4*minimize_phase(theta, b)/(np.pi**2*q)
        w = 1/lengths
        W, Dw = w.sum(), w[::2].sum()-w[1::2].sum()
        Emin = -np.pi*W/24-np.pi*Dw*Dw/(8*W)
        numeric_E = -np.pi*W/6+minimize_phase(theta, w)/(2*np.pi)
        assert abs(numeric-chi) < 2e-11
        assert abs(Emin-numeric_E) < 2e-10
        assert abs(Db) <= 1 and chi <= upper
        phi_energy = -np.pi*w[1::2].sum()/W
        actual_chi = 4*B/(3*q)-1/3-4*np.dot(b, wrapped(phi_energy, theta)**2)/(np.pi**2*q)
        assert actual_chi <= chi+1e-13
        # Direct loop derivative at fixed phase matches the local stress.
        eps = 1e-5
        a2 = wrapped(phi_energy, theta)**2
        coeff = -np.pi/6+a2/(2*np.pi)
        derivative = (np.sum(coeff/(lengths+eps))-np.sum(coeff/(lengths-eps)))/(2*eps)
        stress = np.sum(coeff/(lengths*lengths))
        assert abs(derivative+stress)/max(1, abs(stress)) < 3e-9
        delayed.append(dict(q=q, d_over_T=ratio, B=float(B), alternating_B=float(Db),
                            best_support=float(chi), rigorous_upper_bound=float(upper),
                            support_at_energy_minimum=float(actual_chi),
                            energy_minimum=float(Emin), phase_stress_difference=float(chi-actual_chi)))
    return dict(equal_lengths=rows, winding_integral=winding, point_dipole_delays=delayed)


def basis(r):
    z = r-1
    return (z*z*(r+2)/r, 2*r-2/r**2,
            (2*r-1)*(5*r*r-5*r+2)/(30*r*z**3),
            -(5*z**4+15*z**3+21*z*z+13*z+3)/(15*r*r*z**4))


def matched_checks():
    rows = []
    rout = 21.
    for lam in [0.01, 0.003, 0.001, 0.0003, 0.0001]:
        R = lam**(-1/3)
        z0, rs = lam*R, 1+lam*R
        f = z0*z0/(rs*rs)
        H, g, gp = R*(1+R*R), 1+3*R*R, 6*R
        u0, u0p, u1, u1p = basis(rs)
        matrix = np.array([[u1, -g], [f*u1p, -lam*(1+R*R)*gp]])
        bE, kE = np.linalg.solve(matrix, [-u0, -f*u0p])
        bJ, kJ = np.linalg.solve(matrix, [2*(1/rs-1), -2*f/rs**2])

        def frame(e, J, b, k):
            def integrand(t):
                z = np.exp(t)
                r = 1+z
                return z*(-6*J-2*(e*basis(r)[0]+b*basis(r)[2]-2*J/r))/r**4
            integral = quad(integrand, np.log(z0), np.log(rout-1),
                            epsabs=3e-13, epsrel=3e-13)[0]
            return lam*(-2*k*H-2*J*R)+integral
        wE, wJ = frame(1, 0, bE, kE), frame(0, 1, bJ, kJ)
        J = -wE/wJ
        b, k = bE+J*bJ, kE+J*kJ
        assert abs(frame(1, J, b, k)) < 1e-11

        def rhs(t, y, param):
            z = np.exp(t)
            r = 1+z
            return np.vstack((y[1], (2+4/r**2)*y[0]-(2/r-1)*y[1]+12*param[0]/r**2,
                              z*(-6*param[0]-2*y[0])/r**4))

        def bc(left, right, param):
            j = param[0]
            ou0, ou0p, ou1, ou1p = basis(rout)
            kin = (left[0]+2*j)/g
            return np.array([left[1]-z0/f*lam*(1+R*R)*gp*kin,
                             left[2]-lam*(-2*kin*H-2*j*R),
                             right[1]/(rout-1)-ou0p-(right[0]-ou0+2*j/rout)*ou1p/ou1-2*j/rout**2,
                             right[2]])
        refinements = []
        for size, tol in [(120, 2e-7), (360, 2e-9)]:
            mesh = np.linspace(np.log(z0), np.log(rout-1), size)
            rr = 1+np.exp(mesh)
            # Initial data omit spin and inner load; this linear BVP solves for J.
            guess = np.vstack((basis(rr)[0], (rr-1)*basis(rr)[1], np.zeros(size)))
            sol = solve_bvp(rhs, bc, mesh, guess, p=[0.], tol=tol, max_nodes=16000)
            assert sol.success, sol.message
            test = np.linspace(mesh[0], mesh[-1], 401)
            rad = 1+np.exp(test)
            exact = basis(rad)[0]+b*basis(rad)[2]-2*J/rad
            got = sol.sol(test)
            error = np.max(abs(got[0]-exact)/np.maximum(1, abs(exact)))
            jerr = abs(sol.p[0]-J)
            kerr = abs((got[0, 0]+2*sol.p[0])/g-k)/max(abs(k), 1e-20)
            assert max(error, jerr, kerr) < 3e-6
            refinements.append(dict(nodes=int(sol.x.size), tolerance=tol,
                                    flux_relative_error=float(error), charge_error=float(jerr),
                                    throat_coefficient_relative_error=float(kerr)))
        assert refinements[-1]['charge_error'] < refinements[0]['charge_error']
        w0out = 2/rout-(2*rout-1)/rout**4
        wKNout = (2*rout-1)/rout**4
        asym = -(1-w0out)/(1-wKNout)
        rows.append(dict(lam=lam, R=R, r_out=rout, J_over_e=float(J),
                         asymptotic_J_over_e=asym, k=float(k), b=float(b),
                         electric_k_over_lambda_squared=float(kE/lam**2),
                         zero_charge_outer_rotation=float(wE), refinements=refinements))
    return rows


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    result = dict(status='passed_specified_gluing_support_and_odd_classical_matching',
                  python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, sympy=s.__version__,
                  exact_checks=exact_checks(), geometry=geometry_checks(), phases=phase_checks(), matching=matched_checks(),
                  script_sha256=digest(Path(__file__)),
                  input_sha256={n: digest(ROOT/n) for n in ['notes/gluing-support-compatibility.md',
                      'notes/global-mouth-trial.md', 'notes/binary-mouth-drive.md', 'notes/stationary-electric-dipole.md',
                      'notes/lll-return-holonomy.md', 'notes/collective-clock-response.md', 'notes/rotational-current-response.md']},
                  compensating_rotational_connection_excluded_from_tested_model=True,
                  finite_flux_physical_return_proved=False, complete_quantum_state_selected=False,
                  smooth_semiclassical_junction=False, novelty_established=False, lean_invocations=0, goal_complete=False)
    (ROOT/'evidence/gluing-support-compatibility.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], exact_checks=len(result['exact_checks']),
                         phase_cases=len(result['phases']['equal_lengths'])+len(result['phases']['point_dipole_delays']),
                         matching_cases=len(result['matching']),
                         max_refined_charge_error=max(x['refinements'][-1]['charge_error'] for x in result['matching']))))


if __name__ == '__main__':
    main()
