"""Check the proper-mouth static completion and isolated-patch Weyl action.

No numerical sample certifies the complete transparent influence action.
Run with Python, NumPy, SciPy, SymPy and mpmath; no Lean invocation.
"""
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT/'notes/proper-clock-anomaly.md'
OUT = ROOT/'evidence/proper-clock-anomaly.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integrate(fun, lo, hi):
    return quad(fun, lo, hi, epsabs=2e-11, epsrel=2e-11, limit=200)[0]


def symbolic_checks():
    x, e, k = sy.symbols('x e k', positive=True)
    b = sy.sqrt(1-e*e*x*x)
    a = sy.acos(e*x)
    E = -(b/e-x*a)/(12*sy.pi)
    D = lambda f: x*sy.diff(f, x)
    checks = {
        'core_force': sy.simplify(D(E)-x*a/(12*sy.pi)),
        'core_stiffness': sy.simplify(D(D(E))-(x*a-e*x*x/b)/(12*sy.pi)),
        'moving_domain_term': sy.simplify(x*sy.tan(a)**2*e**3*x*k*k/(2*b**3)-e*k*k/(2*b)),
        'boundary_term': sy.simplify(-2*x*sy.tan(a)*(e*k/b)**2+2*e*k*k/b),
    }
    v, vp, rp = sy.symbols('v vp rp', positive=True)
    checks['covariant_square'] = sy.expand(vp*vp/v+2*vp*rp+v*rp*rp-v*(vp/v+rp)**2)
    xp, sp, spp = sy.symbols('xp sp spp', real=True)
    speed = sy.sqrt(x*x-sp*sp)
    for side in [-1, 1]:
        tangent = sy.Matrix([x, side*sp])/speed
        outward = sy.Matrix([sp, side*x])/speed
        accel = (tangent.diff(x)*xp+tangent.diff(sp)*spp)/speed
        curvature = -outward[0]*accel[0]+outward[1]*accel[1]
        checks[f'outward_boundary_curvature_{side}'] = sy.simplify(
            curvature*speed-(x*spp-sp*xp)/(x*x-sp*sp))
    for name, value in checks.items():
        assert value == 0, (name, value)
    return [{'identity': name, 'residual': str(value)} for name, value in checks.items()]


def static_metrics():
    rows = []
    for amp in [.07, .3, .65]:
        def fields(s):
            lam = amp*(math.cos(s)+.3*math.sin(2*s))
            lp = amp*(-math.sin(s)+.6*math.cos(2*s))
            rho = .4*math.sin(s)-.12*math.cos(3*s)
            rp = .4*math.cos(s)+.36*math.sin(3*s)
            rpp = -.4*math.sin(s)+1.08*math.cos(3*s)
            v = math.exp(lam)
            # N=exp(lam+rho), dy=ds/v. The optical second derivative
            # integrates to zero, including both geometric contributions.
            npp = -amp*(math.cos(s)+1.2*math.sin(2*s))+rpp
            return v, lp, rp, npp
        J = integrate(lambda s: 1/fields(s)[0], 0, 2*math.pi)
        canonical = -math.pi/(6*J)-integrate(lambda s: fields(s)[0]*fields(s)[1]**2, 0, 2*math.pi)/(24*math.pi)
        weyl = integrate(lambda s: fields(s)[0]*(2*fields(s)[1]*fields(s)[2]+fields(s)[2]**2), 0, 2*math.pi)/(24*math.pi)
        covariant = -math.pi/(6*J)-integrate(lambda s: fields(s)[0]*(fields(s)[1]+fields(s)[2])**2, 0, 2*math.pi)/(24*math.pi)
        stress = -math.pi/(6*J)+integrate(lambda s: fields(s)[0]*(2*(fields(s)[1]*(fields(s)[1]+fields(s)[2])+fields(s)[3])-(fields(s)[1]+fields(s)[2])**2), 0, 2*math.pi)/(24*math.pi)
        err = max(abs(canonical-weyl-covariant), abs(stress-covariant))
        assert err < 2e-12
        rows.append(dict(amplitude=amp, canonical_energy=canonical, weyl_action=weyl,
                         covariant_energy=covariant, integrated_stress_energy=stress, error=err))
    return rows


def core_checks():
    rows = []
    mp.mp.dps = 45
    for e in [.01, .04, .1, .2]:
        x = 1.1
        a = math.acos(e*x)
        raw = 2*integrate(lambda sig: x*math.tan(sig)**2, 0, a)
        expected = 2*(math.sqrt(1-e*e*x*x)/e-x*a)
        assert abs(raw-expected) < 2e-10
        ee, xx = mp.mpf(str(e)), mp.mpf(str(x))
        def energy(phi):
            y = xx*mp.exp(phi)
            return -(mp.sqrt(1-ee*ee*y*y)/ee-y*mp.acos(ee*y))/(12*mp.pi)
        numerical_stiffness = float(mp.diff(energy, 0, 2))
        stiffness = (x*a-e*x*x/math.sqrt(1-e*e*x*x))/(12*math.pi)
        bound = e*x*x/(6*math.pi*math.sqrt(15/16))
        assert abs(numerical_stiffness-stiffness) < 2e-14
        assert abs(stiffness-x/24) <= bound
        rows.append(dict(epsilon=e, x=x, spatial_integral=raw, formula=expected,
                         integration_error=abs(raw-expected), stiffness=stiffness,
                         independent_stiffness=numerical_stiffness, leading_error=abs(stiffness-x/24),
                         analytic_bound=bound))
    return rows


def joins():
    rows = []
    mp.mp.dps = 35
    def correction(z):
        b = mp.sqrt(1-z*z)
        def integrand(u):
            den = mp.cos(z*u)-b*mp.sin(z*u)/z
            denp = -z*mp.sin(z*u)-b*mp.cos(z*u)
            rin = mp.log(1-u)-mp.log(den)
            rp = -1/(1-u)-denp/den
            v = u+1
            chi = 1-10*v**3+15*v**4-6*v**5
            chip = -30*v*v+60*v**3-30*v**4
            inner = -denp/den
            blended = 1/(1-u)+chip*rin+chi*rp
            return blended**2-inner**2
        # Two mouths and c=1, with epsilon=z and x=1 for tabulation.
        return -mp.quad(integrand, [-1, 0])/(12*mp.pi)
    for e in [.02, .06, .12, .24]:
        ee = mp.mpf(str(e))
        fun = lambda phi: correction(ee*mp.exp(phi))/ee
        values = [float(mp.diff(fun, 0, k)) for k in range(3)]
        scaled = [v/e for v in values]
        assert max(abs(v) for v in scaled) < 1
        rows.append(dict(epsilon=e, energy_force_stiffness_correction=values,
                         normalized_by_epsilon_x_squared=scaled,
                         note='Samples of the analytic O(epsilon x^2) join bound; not a uniform numerical proof.'))
    return rows


def return_checks():
    rows = []
    mp.mp.dps = 40
    for e in [.003, .03, .12, .24]:
        for D in [0, .1, 1, 100, 1e8]:
            ee, dd = mp.mpf(str(e)), mp.mpf(str(D))
            def difference(phi):
                x = mp.exp(phi)
                return -mp.pi*x/6*(1/(2*mp.acos(ee*x)+dd*x)-1/(mp.pi+dd*x))
            values = [float(mp.diff(difference, 0, k)) for k in range(3)]
            assert max(abs(v/e) for v in values) < 1
            rows.append(dict(epsilon=e, return_length=D, energy_force_stiffness=values,
                             normalized_by_epsilon=[v/e for v in values]))
    return rows


def moving_patch():
    rows = []
    omega, amp = 1.3, .2
    period = 2*math.pi/omega
    for e in [.1, .05, .025, .0125]:
        def data(t):
            X = math.exp(amp*math.cos(omega*t))
            Xp = -amp*omega*math.sin(omega*t)*X
            Xpp = (-amp*omega**2*math.cos(omega*t)+(amp*omega*math.sin(omega*t))**2)*X
            b = math.sqrt(1-e*e*X*X)
            s = math.acos(e*X)
            sp = -e*Xp/b
            spp = -e*(Xpp/b+e*e*X*Xp*Xp/b**3)
            x = math.sqrt(X*X+sp*sp)
            xp = (X*Xp+sp*spp)/x
            beta = sp/x
            betap = (spp*x-sp*xp)/(x*x)
            rho = -math.log(math.cos(s))
            bulk = x*(math.tan(s)-s)-1/e
            direct = bulk+2*rho*betap/(1-beta*beta)
            ibp = bulk-2*x*beta*math.tan(s)*math.atanh(beta)
            bx = math.sqrt(1-e*e*x*x)
            approx = (bx-1)/e-x*math.acos(e*x)-3*e*(xp/x)**2/(2*bx)
            lapse_error = abs((x*x-sp*sp)/math.cos(s)**2-1/e**2)*e**2
            return direct, ibp, approx, lapse_error, x, xp, s
        direct = integrate(lambda t: data(t)[0], 0, period)/(12*math.pi)
        ibp = integrate(lambda t: data(t)[1], 0, period)/(12*math.pi)
        approx = integrate(lambda t: data(t)[2], 0, period)/(12*math.pi)
        lapse_error = max(data(t)[3] for t in np.linspace(0, period, 101))
        assert abs(direct-ibp) < 3e-11
        assert lapse_error < 5e-13
        error = abs(ibp-approx)
        assert error < e**3
        rows.append(dict(epsilon=e, period=period, direct_action=direct, integrated_by_parts_action=ibp,
                         boundary_integration_error=abs(direct-ibp), two_derivative_action=approx,
                         derivative_remainder=error, remainder_over_epsilon_cubed=error/e**3,
                         induced_lapse_relative_error=lapse_error))
    for old, new in zip(rows, rows[1:]):
        order = math.log(old['derivative_remainder']/new['derivative_remainder'], 2)
        assert 2.8 < order < 3.2, order
        new['observed_remainder_order'] = order
    return rows


def main():
    start = time.monotonic()
    result = dict(status='passed_proper_clock_anomaly_checks', date='2026-10-11',
        exact_checks=symbolic_checks(), static_metric_cases=static_metrics(), core_cases=core_checks(),
        join_cases=joins(), return_cases=return_checks(), moving_patch_cases=moving_patch(),
        static_completion_scope='specified free-CFT two-dimensional AdS core and controlled static joins',
        full_transparent_influence_action_derived=False, physical_four_dimensional_matching_closed=False,
        physical_clock_inertia_derived=False, physical_damping_proved=False, novelty_established=False,
        lean_invocations=0, note_sha256=sha(NOTE), script_sha256=sha(Path(__file__)),
        dependency_hashes={name: sha(ROOT/name) for name in ['notes/collar-conservative-response.md',
            'notes/collective-clock-response.md', 'notes/weak-coupling-window.md']},
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                      sympy=sy.__version__, mpmath=mp.__version__))
    result['elapsed_seconds'] = time.monotonic()-start
    OUT.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: result[k] for k in ['status', 'elapsed_seconds', 'physical_clock_inertia_derived']}))
    print(json.dumps({'maximum_static_identity_error': max(v['error'] for v in result['static_metric_cases']),
        'maximum_induced_lapse_error': max(v['induced_lapse_relative_error'] for v in result['moving_patch_cases']),
        'remainder_orders': [v['observed_remainder_order'] for v in result['moving_patch_cases'][1:]]}))


if __name__ == '__main__':
    main()
