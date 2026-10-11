"""Check the free-channel return-holonomy criterion in the accompanying note.

No actual exterior, interacting vacuum, or four-dimensional solution is
computed. The regulator and determinant checks are independent numerical
checks of the ordinary finite-channel argument, not certified error bounds.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
import sympy as sp


def angles(value):
    return (np.asarray(value) + np.pi) % (2*np.pi) - np.pi


def energy(theta, phase=0.0, length=1.0):
    a = angles(np.asarray(theta)+phase)
    return -len(a)*np.pi/(6*length) + np.dot(a, a)/(2*np.pi*length)


def minimize_phase(theta):
    """Enumerate the exact quadratic pieces, not a local optimizer/grid."""
    theta = np.asarray(theta)
    cuts = np.unique(np.r_[0.0, (np.pi-theta) % (2*np.pi), 2*np.pi])
    candidates = list(cuts)
    for left, right in zip(cuts[:-1], cuts[1:]):
        middle = (left+right)/2
        vertex = middle - np.mean(angles(theta+middle))
        candidates.append(float(np.clip(vertex, left, right)))
    phase = min(candidates, key=lambda p: energy(theta, p))
    d2 = float(np.dot(angles(theta+phase), angles(theta+phase)))
    return {'phase': float(phase), 'energy_L1': float(energy(theta, phase)),
            'D2': d2, 'chi': 1-4*d2/(len(theta)*np.pi**2)}


def spin_matrices(q):
    j = (q-1)/2
    m = np.arange(q)-j
    raising = np.zeros((q, q))
    for k in range(q-1):
        raising[k+1, k] = np.sqrt((j-m[k])*(j+m[k]+1))
    return np.diag(m), (raising+raising.T)/2


def determinant_energy(unitary):
    q = len(unitary)
    def integrand(xi):
        sign, logabs = np.linalg.slogdet(np.eye(q)+np.exp(-xi)*unitary)
        return -2*logabs/np.pi
    result, error = quad(integrand, 0.0, np.inf, epsabs=2e-11, epsrel=2e-11)
    return result, error


def regulated_spectrum(generator, cutoff):
    """Dirac H=diag(K,-K); vacuum energy is -tr|K|, not -tr|K|/2.

    Regulate physical energies by exp(-cutoff*|H|), then subtract the
    infinite-line integral with the identical regulator. The extra finite
    index truncation is taken deep into the exponential tail.
    """
    q = len(generator)
    limit = int(np.ceil(45/(2*np.pi*cutoff)))+q
    total = 0.0
    for n in range(-limit, limit+1):
        k = 2*np.pi*(n+0.5)*np.eye(q)+generator
        vals = np.abs(np.linalg.eigvalsh(k))
        total -= float(np.sum(vals*np.exp(-cutoff*vals)))
    return total+q/(np.pi*cutoff**2)


def verify():
    q, L, a, d, k = sp.symbols('q L a d k', positive=True)
    n = sp.symbols('n', integer=True, nonnegative=True)
    j = (q-1)/2
    e = -q*sp.pi/(6*L)+d/(2*sp.pi*L)
    chi = 1-4*d/(q*sp.pi**2)
    m2 = sp.summation((n-j)**2, (n, 0, q-1))
    roots_d2 = sp.pi**2*(q*q-1)/(3*q)
    identities = {
        'Bernoulli_single_channel': 2*sp.pi/L*((sp.Rational(1,2)+a/(2*sp.pi))**2-(sp.Rational(1,2)+a/(2*sp.pi))+sp.Rational(1,6)) - (-sp.pi/(6*L)+a*a/(2*sp.pi*L)),
        'reflecting_reference_subtraction': e/L+q*sp.pi/(24*L**2)+q*sp.pi*chi/(8*L**2),
        'spin_second_moment': m2-q*(q*q-1)/12,
        'root_phase_spread': (2*sp.pi/q)**2*m2-roots_d2,
        'cyclic_energy': e.subs(d, roots_d2)+sp.pi/(6*q*L),
        'cyclic_support': chi.subs(d, roots_d2)-(-sp.Rational(1,3)+4/(3*q*q)),
        'rigid_rotation_support': chi.subs(d, a*a*m2)-(1-a*a*(q*q-1)/(3*sp.pi**2)),
    }
    for name, value in identities.items():
        assert sp.simplify(value) == 0, name

    cyclic = []
    for count in [1, 2, 3, 4, 8, 16, 64]:
        theta = 2*np.pi*np.arange(count)/count
        result = minimize_phase(theta)
        assert abs(result['energy_L1']+np.pi/(6*count)) < 2e-13
        assert abs(result['chi']-(-1/3+4/(3*count**2))) < 2e-14
        cyclic.append({'q': count, **result})

    # Direct matrices test conjugacy invariance and the winding determinant.
    determinant = []
    spectral = []
    rng = np.random.default_rng(20261010)
    for count, scale in [(3, 2*np.pi/3), (5, 0.23), (7, 0.7)]:
        jz, jx = spin_matrices(count)
        generator = scale*jx
        unitary = expm(1j*generator)
        theta = np.angle(np.linalg.eigvals(unitary))
        direct, quad_error = determinant_energy(unitary)
        expected = energy(theta)
        assert abs(direct-expected) < 3e-10
        basis, _ = np.linalg.qr(rng.normal(size=(count,count))+1j*rng.normal(size=(count,count)))
        rotated = basis@unitary@basis.conj().T
        rotated_energy, _ = determinant_energy(rotated)
        assert abs(rotated_energy-direct) < 2e-10
        determinant.append({'q': count, 'rotation': scale, 'energy_L1': expected,
                            'determinant_residual': abs(direct-expected),
                            'quadrature_error_estimate': quad_error,
                            'basis_change_residual': abs(rotated_energy-direct)})
        cutoffs = [0.12, 0.06, 0.03, 0.015]
        regulated = [regulated_spectrum(generator, eps) for eps in cutoffs]
        extrap = [(regulated[i]-20*regulated[i+1]+64*regulated[i+2])/45 for i in range(2)]
        residuals = [abs(x-expected) for x in extrap]
        assert residuals[1] < 3e-7 and residuals[1] < residuals[0]/10
        spectral.append({'q': count, 'cutoffs': cutoffs, 'regulated_energies': regulated,
                         'extrapolated_energies': extrap, 'residuals': residuals})

    # A rigid sphere rotation already realizes the centered root spectrum.
    rotations = []
    for count in [3, 8, 32, 128, 1024]:
        delta = 2*np.pi/count
        theta = delta*(np.arange(count)-(count-1)/2)
        result = minimize_phase(theta)
        fidelity = np.cos(delta/2)**(2*(count-1))
        assert abs(result['chi']-(-1/3+4/(3*count**2))) < 8e-14
        case = {'q': count, 'rotation': delta,
                          'rotation_over_width_scale': delta*np.sqrt(count),
                          'equatorial_coherent_state_fidelity': float(fidelity),
                          'chi': result['chi']}
        if count <= 32:
            jz, jx = spin_matrices(count)
            _, vectors = np.linalg.eigh(jx)
            equator = vectors[:, -1]
            matrix_fidelity = abs(np.vdot(equator, expm(1j*delta*jz)@equator))**2
            case['coherent_state_matrix_residual'] = abs(matrix_fidelity-fidelity)
            assert case['coherent_state_matrix_residual'] < 2e-14
        rotations.append(case)

    # The positive-support bound is uniform in q at delta=c/q.
    slow = []
    for count in [3, 16, 64, 256]:
        c = 0.2
        theta = c/count*(np.arange(count)-(count-1)/2)
        expected_chi = 1-c*c*(1-1/count**2)/(3*np.pi**2)
        actual_chi = 1-4*float(np.dot(theta,theta))/(count*np.pi**2)
        assert abs(actual_chi-expected_chi)<2e-15
        assert minimize_phase(theta)['chi'] >= expected_chi-2e-14
        slow.append({'q': count, 'q_times_rotation': c, 'chi_at_common_phase_zero': actual_chi})

    # Numerical checks of sufficient chord/phase bounds; the note proves them.
    bounds = []
    for count in [2, 3, 8, 20]:
        theta = rng.uniform(-0.9, 0.9, count)
        chord2 = float(np.sum(4*np.sin(theta/2)**2))
        a2 = float(np.dot(theta, theta))
        measured = 1-4*a2/(count*np.pi**2)
        assert measured >= 1-chord2/count-2e-14
        assert measured >= 1-4*0.9**2/np.pi**2-2e-14
        bounds.append({'q': count, 'chi': measured, 'Frobenius_lower_bound': 1-chord2/count})
    return {'symbolic': {name: 'exact_zero' for name in identities},
            'cyclic_phase_minima': cyclic, 'matrix_determinant_checks': determinant,
            'independent_regulated_spectra': spectral, 'rigid_rotation_cases': rotations,
            'slow_rotation_cases': slow, 'sufficient_bound_cases': bounds}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = {'status': 'passed_LLL_return_holonomy_cross_checks',
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'checks': verify(),
              'model': 'q free complex massless 2D Dirac channels, equal optical length, fixed energy-independent unitary return and common U(1) phase',
              'reference': 'LLL state difference relative to q identical conformally transformed reflecting strips; negligible exterior optical length',
              'support_observable': 'Angular trace / spherical average only',
              'arbitrary_unitary_support_preservation_refuted_in_model': True,
              'actual_MMP_return_map_derived': False, 'actual_MMP_support_refuted': False,
              'interacting_vacuum_computed': False, 'full_support_stress_controlled': False,
              'physical_stability_proved': False, 'novelty_established': False,
              'full_numerical_error_certified': False, 'lean_certified': False, 'lean_invocations': 0}
    args.output.write_bytes((json.dumps(result, indent=2)+'\n').encode())
    print(json.dumps({'status': result['status'], 'script_sha256': result['script_sha256'],
                      'spectral_residuals': [x['residuals'] for x in result['checks']['independent_regulated_spectra']]}))


if __name__ == '__main__':
    main()
