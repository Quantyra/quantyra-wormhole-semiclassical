"""Local fermion coset and two-frequency current return checks."""
import hashlib
import json
import platform
import runpy
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import sici

ROOT = Path(__file__).resolve().parents[1]
geometry = runpy.run_path(str(ROOT/'checks/gluing_support_compatibility.py'))
KAPPA = geometry['KAPPA']
CF = 25*KAPPA/16


def generators(q):
    m = np.arange(q)-(q-1)/2
    jp = np.diag(np.sqrt((np.arange(1, q))*(q-np.arange(1, q))), -1)
    return [(jp+jp.T)/2, (jp-jp.T)/(2j), np.diag(m)], jp


def inner(a, b):
    return sum(np.conjugate(value)*b.get(key, 0) for key, value in a.items())


def add(*terms):
    out = {}
    for coeff, state in terms:
        for key, value in state.items():
            out[key] = out.get(key, 0)+coeff*value
    return {key: value for key, value in out.items() if abs(value) > 1e-14}


def car_bilinear(bits, creation, annihilation):
    if not ((bits >> annihilation) & 1):
        return None
    sign = (-1)**((bits & ((1 << annihilation)-1)).bit_count())
    bits ^= 1 << annihilation
    if (bits >> creation) & 1:
        return None
    sign *= (-1)**((bits & ((1 << creation)-1)).bit_count())
    return bits | (1 << creation), sign


def current(state, X, mode, q):
    out = {}
    entries = [(i, j, X[i, j]) for i, j in zip(*np.nonzero(X))]
    for bits, amp in state.items():
        for level in range(4):
            source = level+mode
            if not 0 <= source < 4:
                continue
            for i, j, coef in entries:
                value = car_bilinear(bits, level*q+int(i), source*q+int(j))
                if value is not None:
                    key, sign = value
                    out[key] = out.get(key, 0)+amp*coef*sign
        if mode == 0:
            out[bits] = out.get(bits, 0)-2*np.trace(X)*amp
    return {key: value for key, value in out.items() if abs(value) > 1e-13}


def fock_checks():
    rows = []
    for q in [2, 3, 4, 5, 6, 8]:
        Ja, _ = generators(q)
        k = q*(q*q-1)//6
        vac = {(1 << (2*q))-1: 1.}
        # Four NS levels -3/2,-1/2,1/2,3/2 suffice at total level two.
        free = {}
        for level, coefficient in [(2, -.5), (3, .5)]:
            for i in range(q):
                value = car_bilinear(next(iter(vac)), level*q+i, (level-2)*q+i)
                key, sign = value
                free[key] = coefficient*sign
        gauge = add(*[(1/(k+2), current(current(vac, X, -1, q), X, -1, q)) for X in Ja])
        electric = add((1/(2*q), current(current(vac, np.eye(q), -1, q), np.eye(q), -1, q)))
        coset = add((1, free), (-1, gauge), (-1, electric))
        c = Fraction(q-4)+Fraction(6, k+2)
        tests = [abs(inner(free, free)-q/2), abs(inner(gauge, gauge)-3*k/(2*(k+2))),
                 abs(inner(electric, electric)-.5), abs(inner(coset, coset)-float(c)/2),
                 abs(inner(electric, gauge)), abs(inner(coset, gauge)), abs(inner(coset, electric))]
        annihilation = []
        for X in [np.eye(q), *Ja]:
            for mode in [0, 1, 2]:
                state = current(coset, X, mode, q)
                annihilation.append(abs(inner(state, state)))
        assert max(tests+annihilation) < 2e-11, (q, tests, annihilation)
        # The quadrupole is an affine primary of angular spin two.
        if q >= 3:
            X = Ja[2]@Ja[2]-np.eye(q)*(q*q-1)/12
            casimir = sum(J@(J@X-X@J)-(J@X-X@J)@J for J in Ja)
            assert np.linalg.norm(casimir-6*X) < 1e-10
            assert max(abs(np.trace(J@X)) for J in Ja) < 1e-11
        rows.append(dict(q=q, k=k, c_coset=str(c), norm_central_error=float(max(tests)),
                         maximum_current_annihilation_norm_squared=float(max(annihilation)),
                         coset_level_two_states=len(coset),
                         quadrupole_coset_weight=str(1-Fraction(6, k+2)) if q >= 3 else None))
    return rows


def leakage(q, delays, w1, w2):
    n = np.arange(q)
    m = n-(q-1)/2
    Tq = q*(q*q-1)/12
    phase_z = np.exp(1j*(w2-w1)*delays)
    alpha = np.dot(m*m, phase_z)/Tq
    singlet = abs(np.dot(m, phase_z))**2/(q*Tq)
    longitudinal = max(0., float(1-abs(alpha)**2-singlet))
    weights = (n[1:])*(q-n[1:])
    phases = np.exp(1j*(w2*delays[:-1]-w1*delays[1:]))
    average = np.dot(weights, phases)/weights.sum()
    transverse = max(0., float(1-abs(average)**2))
    return longitudinal, transverse


def bounds(q, delta, w1, w2):
    x = np.arange(1, q, dtype=float)
    eta = abs(w2-w1)
    C = 2*KAPPA*eta*delta*q
    sep = 6/(q*q-1)*np.sum(x*np.minimum(4, (C/x)**2))
    A = CF*max(abs(w1), abs(w2))*delta*q
    B = 4*A/3
    adjacent = 6/(q*q-1)*(2*B+2*np.sqrt(2*B))
    xx = np.arange(1, q+1, dtype=float)
    longitudinal = 3*q/(q*q-1)*np.sum(np.minimum(4, (C/xx)**2))
    return min(1., longitudinal), min(1., 2*(sep+adjacent)), min(1., adjacent)


def projection_checks():
    rows = []
    for q in [2, 3, 4, 7, 12, 19]:
        Ja, jp = generators(q)
        Xs = [np.eye(q), *Ja]
        delay = .17*KAPPA/((np.arange(q)+.5)/q)
        for w1, w2 in [(0.8, 0.8), (0.3, 1.2), (-.4, .7)]:
            D1, D2 = np.diag(np.exp(1j*w1*delay)), np.diag(np.exp(1j*w2*delay))
            results = []
            for X in [Ja[2], jp, jp.T, Ja[0], Ja[1]]:
                M = D1.conj().T@X@D2
                proj = sum(np.trace(Y.conj().T@M)/np.trace(Y.conj().T@Y)*Y for Y in Xs)
                results.append(float(np.linalg.norm(M-proj, 'fro')**2/np.linalg.norm(X, 'fro')**2))
            lz, lp = leakage(q, delay, w1, w2)
            lm = leakage(q, delay, w2, w1)[1]
            scalar = [lz, lp, lm, (lp+lm)/2, (lp+lm)/2]
            assert np.max(abs(np.array(results)-scalar)) < 1e-12
            assert max(leakage(q, np.ones(q)*.7, w1, w2)) < 1e-13
            rows.append(dict(q=q, omega1=w1, omega2=w2, longitudinal=results[0],
                             plus=results[1], minus=results[2], x=results[3], y=results[4]))
    return rows


def delay_cases():
    rows = []
    for kind, qs in [('inverse_tail', [64, 256, 1024, 4096, 16384, 65536]),
                     ('point_dipole', [32, 128, 512, 2048])]:
        for q in qs:
            delta = q**(-.4)
            u = (np.arange(q)+.5)/q
            f = KAPPA/u if kind == 'inverse_tail' else np.array([geometry['g_of_u'](x)/x for x in u])
            delay = delta*f
            jump_bound = CF*delta*q/((np.arange(q-1)+.5)*(np.arange(q-1)+1.5))
            jump_ratio = float(np.max(abs(np.diff(delay))/jump_bound))
            assert jump_ratio <= 1+1e-10 and np.all(np.diff(delay) < 0)
            for w1, w2 in [(.8, .8), (.3, 1.2), (-.4, .7)]:
                lz, lp = leakage(q, delay, w1, w2)
                lm = leakage(q, delay, w2, w1)[1]
                bz, bp, ba = bounds(q, delta, w1, w2)
                assert lz <= bz+1e-13 and max(lp, lm) <= bp+1e-13
                if w1 == w2:
                    assert lz < 1e-13 and lp <= ba+1e-13
                rows.append(dict(kind=kind, q=q, delta=delta, omega1=w1, omega2=w2,
                                 longitudinal=lz, transverse=lp, transverse_minus=lm,
                                 transverse_cartesian=(lp+lm)/2, maximum_delay_jump_ratio=jump_ratio,
                                 upper_longitudinal=bz,
                                 upper_transverse=bp, upper_adjacent=ba,
                                 longitudinal_over_delta=lz/delta,
                                 predicted_longitudinal_limit=3*np.pi*KAPPA*abs(w2-w1),
                                 transverse_over_delta=lp/delta))
    return rows


def moments(alpha):
    si, ci = sici(alpha)
    E = np.exp(1j*alpha)
    C0 = E-alpha*(np.pi/2-si)-1j*alpha*ci
    C1 = (E+1j*alpha*C0)/2
    C2 = (E+1j*alpha*C1)/3
    return np.array([C0, C1, C2])


def continuum_checks():
    independent = []
    for alpha in [.02, .1, .4, 1.3]:
        direct = []
        for n in range(3):
            real = quad(lambda x: x**(-n-2), 1, np.inf, weight='cos', wvar=alpha, epsabs=3e-12)[0]
            imag = quad(lambda x: x**(-n-2), 1, np.inf, weight='sin', wvar=alpha, epsabs=3e-12)[0]
            direct.append(real+1j*imag)
        error = float(np.max(abs(moments(alpha)-direct)))
        assert error < 2e-11
        independent.append(dict(alpha=alpha, error=error))
    asymptotics = []
    for delta in [.02, .005, .001, .0002, .00005, .00001]:
        alpha = KAPPA*.9*delta
        C0, C1, C2 = moments(alpha)
        az = 3*(C0-4*C1+4*C2)
        beta = C0-2*C1
        lz = 1-abs(az)**2-3*abs(beta)**2
        lp = 1-abs(6*(C1-C2))**2
        asymptotics.append(dict(delta=delta, longitudinal_over_delta=float(lz/delta),
                                longitudinal_limit=3*np.pi*KAPPA*.9,
                                transverse_over_delta=float(lp/delta)))
    assert abs(asymptotics[-1]['longitudinal_over_delta']/asymptotics[-1]['longitudinal_limit']-1) < .004
    return dict(independent_moments=independent, continuum_asymptotics=asymptotics)


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    result = dict(status='passed_local_current_coset_and_two_frequency_return_not_state_closure',
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        fock=fock_checks(), projections=projection_checks(), delays=delay_cases(), continuum=continuum_checks(),
        script_sha256=digest(Path(__file__)),
        input_sha256={name: digest(ROOT/name) for name in ['notes/angular-current-coset.md',
          'checks/gluing_support_compatibility.py', 'notes/gluing-support-compatibility.md',
          'notes/rotational-current-response.md', 'notes/angular-current-renormalization.md',
          'notes/singlet-stress-continuum.md', 'notes/gauged-mouth-interaction.md']},
        global_coset_factorization_proved=False, physical_quantum_state_selected=False,
        full_renormalized_stress_bound=False, novel_physical_endpoint_accepted=False,
        lean_invocations=0, goal_complete=False)
    (ROOT/'evidence/angular-current-coset.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], fock_cases=len(result['fock']),
        matrix_projection_cases=len(result['projections']), delay_cases=len(result['delays']),
        maximum_fock_error=max(x['norm_central_error'] for x in result['fock']))))


if __name__ == '__main__':
    main()
