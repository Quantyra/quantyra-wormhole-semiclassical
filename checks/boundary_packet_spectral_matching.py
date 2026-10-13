"""S072 boundary packet CAR, radial variation and continuum matching null.

No four-dimensional state or stress is inferred from the finite checks.
The localization and signed-moment counterexample are ordinary proofs.
"""
from pathlib import Path
import hashlib
import json
import warnings

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad, IntegrationWarning
from scipy.linalg import eigh, expm
from scipy.special import digamma, polygamma

warnings.simplefilter('error', IntegrationWarning)
ROOT = Path(__file__).resolve().parents[1]


def annihilators(D):
    result = []
    for mode in range(D):
        c = np.zeros((1 << D, 1 << D), dtype=complex)
        for state in range(1 << D):
            if (state >> mode) & 1:
                sign = (-1)**((state & ((1 << mode)-1)).bit_count())
                c[state ^ (1 << mode), state] = sign
        result.append(c)
    return result


def second_quantized(h, cs):
    result = np.zeros_like(cs[0])
    for a, b in zip(*np.nonzero(h)):
        result += h[a, b]*(cs[a].conj().T@cs[b])
    return result


def spin(q):
    j = (q-1)/2
    m = np.arange(q)-j
    plus = np.diag(np.sqrt(j*(j+1)-m[:-1]*(m[:-1]+1)), -1)
    return [(plus+plus.T)/2, (plus-plus.T)/(2j), np.diag(m)]


def channel(O, packets):
    return sum((c.conj().T@O@c+c@O@c.conj().T for c in packets),
               np.zeros_like(O))/len(packets)


def thermal_spectrum(H, packets):
    energies, vectors = eigh(H)
    probs = np.exp(-(energies-min(energies))); probs /= sum(probs)
    rho = (vectors*probs)@vectors.conj().T
    sqrt_rho = (vectors*np.sqrt(probs))@vectors.conj().T
    delta = energies[None, :]-energies[:, None]
    measure = np.zeros(H.shape, dtype=float)
    for c in packets:
        matrix = vectors.conj().T@c@vectors
        measure += (probs[:, None]+probs[None, :])*abs(matrix)**2/len(packets)
    C = float(np.sum(measure/np.cosh(delta/2)))
    S = float(2*np.sum(measure*delta*np.tanh(delta/2)))
    direct_C = float(np.trace(sqrt_rho@channel(sqrt_rho, packets)).real)
    direct_S = float(2*np.trace(rho@(channel(H, packets)-H)).real)
    assert abs(np.sum(measure)-1) < 2.e-12
    assert abs(C-direct_C) < 2.e-12
    assert abs(S-direct_S) < 2.e-11
    assert 1-C <= S/8+2.e-12
    return dict(C=C, S=S, spectral_mass=float(np.sum(measure)),
                direct_C_error=abs(C-direct_C), direct_S_error=abs(S-direct_S)), rho


def packet_case(q, orbitals, beta, x, a):
    D = q*orbitals
    cs = annihilators(D)
    unit = np.eye(1 << D)
    raw = (1+.13*np.arange(orbitals))*np.exp(.37j*np.arange(orbitals))
    f = raw/np.linalg.norm(raw)
    packets = [sum((f[r].conj()*cs[r*q+m] for r in range(orbitals)),
                   np.zeros_like(cs[0])) for m in range(q)]
    car_error = max(float(np.max(abs(c@d.conj().T+d.conj().T@c-(unit if m == n else 0))))
                    for m, c in enumerate(packets) for n, d in enumerate(packets))
    assert car_error < 3.e-15
    Qt = second_quantized(np.eye(D), cs)-D/2*unit
    Qf = sum((c.conj().T@c for c in packets), np.zeros_like(unit, dtype=complex))-q/2*unit
    ts = spin(q)
    Jt = [second_quantized(np.kron(np.eye(orbitals), t), cs) for t in ts]
    Jf = [second_quantized(t, packets) for t in ts]
    K = x*(Qt@Qt)+a*sum((J@J for J in Jt), np.zeros_like(Qt))
    mixed = x*(Qt@Qf)+a*sum((J@F for J, F in zip(Jt, Jf)), np.zeros_like(Qt))
    lam = x+a*(q*q-1)/4
    compact_error = float(np.max(abs(channel(K, packets)-K-lam*unit+4*mixed/q)))
    assert compact_error < 4.e-12

    radial = np.diag(np.linspace(-1.2, 2.1, orbitals)).astype(complex)
    for r in range(orbitals-1):
        radial[r, r+1] = .17+.27j
        radial[r+1, r] = .17-.27j
    h = np.kron(radial, np.eye(q))
    P = np.kron(np.outer(f, f.conj()), np.eye(q))
    R = beta*second_quantized(h, cs)
    formula = beta/q*(np.trace(P@h)*unit-second_quantized(P@h+h@P, cs))
    radial_error = float(np.max(abs(channel(R, packets)-R-formula)))
    assert radial_error < 5.e-12
    omitted_cross = beta/q*second_quantized(P@h@(np.eye(D)-P)+(np.eye(D)-P)@h@P, cs)
    assert np.linalg.norm(omitted_cross) > .01
    spectrum, rho = thermal_spectrum(K+R, packets)
    mixed_mean = float(np.trace(rho@mixed).real)
    radial_mean = float(np.trace(rho@formula).real)
    target_S = 2*lam-8*mixed_mean/q+2*radial_mean
    assert abs(target_S-spectrum['S']) < 2.e-11
    return dict(q=q, orbitals=orbitals, beta=beta, x=x, a=a,
                CAR_error=car_error, compact_identity_error=compact_error,
                radial_identity_error=radial_error,
                omitted_cross_operator_norm=float(np.linalg.norm(omitted_cross)),
                mixed_charge_expectation=mixed_mean,
                radial_channel_expectation=radial_mean,
                total_first_moment_error=abs(target_S-spectrum['S']), **spectrum)


def stable_sech(x):
    z = np.exp(-abs(x))
    return 2*z/(1+z*z)


def localization():
    rows = []
    beta = 2*np.pi
    for ratio in (1., .2, .05, .01, .001, .0001, .00001):
        epsilon = beta*ratio
        d = 2*ratio
        integral, err = quad(lambda v: np.exp(-d*v)*stable_sech(v),
                             0, np.inf, epsabs=2.e-13, epsrel=2.e-13)
        closed = (digamma((d+3)/4)-digamma((d+1)/4))/2
        assert abs(integral-closed) < 2.e-12
        # Subtract the analytic large-frequency part before quadrature.
        correction, _ = quad(lambda v: v*np.exp(-d*v)*2*np.exp(-2*v)/(1+np.exp(-2*v)),
                             0, np.inf, epsabs=2.e-13, epsrel=2.e-13)
        S = 4/d*(1-d*d*correction)
        S_closed = 4/d*(1-d*d/8*(polygamma(1, (d+2)/4)-polygamma(1, (d+4)/4)))
        assert abs(S-S_closed) < 2.e-12*max(1, S)
        C = d*integral
        N = 1/(np.pi*epsilon)
        boundary = N*C/2
        target = 1/(2*beta)
        regulator_bound = 4*float(mp.catalan)*epsilon/(np.pi*beta**2)
        assert 0 <= target-boundary <= regulator_bound+2.e-14
        assert 0 < C <= 1 and 1-C <= S/8+2.e-12
        rows.append(dict(epsilon_over_beta=ratio, epsilon=epsilon, beta=beta,
                         canonical_norm=N, C=C, S=S, boundary_midpoint=boundary,
                         point_limit=target, regulator_error_bound=regulator_bound,
                         C_asymptotic_ratio=C/(np.pi*ratio),
                         S_asymptotic_ratio=S*ratio/2,
                         digamma_error=abs(integral-closed), quadrature_error=err))
    return rows


def moment_counterexample():
    mp.mp.dps = 65
    F = lambda e: e*mp.tanh(e/2)
    r = F(1)/F(3)
    s = 1/mp.cosh(mp.mpf(1)/2)
    reference = 1-r+r/mp.cosh(mp.mpf(3)/2)
    gap = reference-s
    closed_gap = 4*(1-s)**2*(s+2)/(3*(4-s*s))
    assert gap > 0 and abs(gap-closed_gap) < mp.mpf('1e-60')
    cs = annihilators(2)
    packet = [float(mp.sqrt(1-r))*cs[0]+float(mp.sqrt(r))*cs[1]]
    ref, _ = thermal_spectrum(second_quantized(np.diag([0., 3.]), cs), packet)
    candidate, _ = thermal_spectrum(second_quantized(np.eye(2), cs), packet)
    assert abs(ref['S']-candidate['S']) < 2.e-14
    assert abs(ref['C']-candidate['C']-float(gap)) < 2.e-14
    assert candidate['C'] < ref['C']
    embedded = []
    nu = .3
    for epsilon in (.2, .05, .01, .001):
        d = 2*epsilon
        bulkC = d*(digamma((d+3)/4)-digamma((d+1)/4))/2
        bulkS = 4/d*(1-d*d/8*(polygamma(1, (d+2)/4)-polygamma(1, (d+4)/4)))
        p = nu*epsilon
        Cref = (1-p)*bulkC+p*ref['C']
        Cnew = (1-p)*bulkC+p*candidate['C']
        Sref = (1-p)*bulkS+p*ref['S']
        Snew = (1-p)*bulkS+p*candidate['S']
        physical_gap = (Cref-Cnew)/(2*np.pi*epsilon)
        expected = nu*float(gap)/(2*np.pi)
        assert abs(physical_gap-expected) < 3.e-15
        assert abs(Sref-Snew) < 2.e-12
        embedded.append(dict(epsilon=epsilon, spectral_mixing_weight=p,
                             subtracted_first_moment=Snew-Sref,
                             physical_midpoint_loss=physical_gap))
    return dict(weight_r=float(r), reference_C=float(reference), C=float(s),
                exact_positive_gap=float(gap), common_first_moment=float(2*F(1)),
                matrix_C_gap_error=abs(ref['C']-candidate['C']-float(gap)),
                matrix_first_moment_difference=candidate['S']-ref['S'],
                common_UV_tail=embedded,
                conclusion='Equal reference-subtracted first moments do not bound the midpoint loss.')


def dressed_recoil():
    rows = []
    beta, omega, level, g = 1.7, .8, .6, .45
    fermion = np.array([[0., 1.], [0., 0.]])
    occupation = fermion.T@fermion
    for modes in (32, 48):
        lowering = np.diag(np.sqrt(np.arange(1, modes)), 1)
        X = (lowering+lowering.T)/np.sqrt(2*omega)
        P = 1j*np.sqrt(omega/2)*(lowering.T-lowering)
        HB = omega*np.diag(np.arange(modes)+.5)
        V = expm(1j*g*X)
        b = np.kron(fermion, V)
        Hbath = np.kron(np.eye(2), HB)
        H = beta*(Hbath+level*np.kron(occupation, np.eye(modes)))
        spectrum, rho = thermal_spectrum(H, [b])
        variation = channel(Hbath, [b])-Hbath
        expected = g*np.kron(2*occupation-np.eye(2), P)+g*g/2*np.eye(2*modes)
        defect = variation-expected
        thermal_squared_error = float(np.trace(rho@(defect.conj().T@defect)).real)
        assert thermal_squared_error < 2.e-13
        target_C = stable_sech(beta*level/2)*np.exp(-g*g/(2*omega)*np.tanh(beta*omega/4))
        target_S = 2*beta*level*np.tanh(beta*level/2)+beta*g*g
        assert abs(spectrum['C']-target_C) < 2.e-12
        assert abs(spectrum['S']-target_S) < 2.e-11
        car_error = float(np.max(abs(b@b.conj().T+b.conj().T@b-np.eye(2*modes))))
        assert car_error < 3.e-14
        rows.append(dict(oscillator_modes=modes, beta=beta, omega=omega,
                         fermion_energy=level, dressing_g=g,
                         CAR_error=car_error, recoil_thermal_squared_error=thermal_squared_error,
                         midpoint_formula_error=abs(spectrum['C']-target_C),
                         first_moment_formula_error=abs(spectrum['S']-target_S), **spectrum))
    assert abs(rows[0]['C']-rows[1]['C']) < 2.e-12
    return rows


def main():
    packets = [packet_case(q, r, beta, x, a) for q, r, beta, x, a in
               ((2, 2, 1.2, .07, .02), (3, 2, .7, .12, .035),
                (2, 3, 1.7, .04, .017), (1, 3, 1.2, .08, .02))]
    result = dict(date='2026-10-13', story='S072', status='PASS',
                  decision='Raw and reference-subtracted CAR first moments cannot alone certify the point boundary transfer.',
                  versions=dict(numpy=np.__version__, scipy=scipy.__version__),
                  packet_cases=packets, free_localization=localization(),
                  moment_counterexample=moment_counterexample(), dressed_recoil=dressed_recoil(),
                  limitations=[
                      'Packet identities hold at a canonical spatial regulator; gravitational dressing remains unmatched.',
                      'Free KMS localization is a required normalization comparator, not the interacting RN spectral density.',
                      'The counterexample rejects an inference, not the candidate wormhole.',
                      'The actual prepared weighted spectrum and driven two-time covariance/stress remain open.'])
    result['artifacts'] = [dict(path=p, sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest())
                           for p in ('notes/boundary-packet-spectral-matching.md',
                                     'checks/boundary_packet_spectral_matching.py')]
    output = ROOT/'evidence/boundary-packet-spectral-matching.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status='PASS', packet_cases=len(packets),
                         localization_cases=len(result['free_localization']),
                         maximum_radial_identity_error=max(r['radial_identity_error'] for r in packets),
                         maximum_spectral_first_moment_error=max(r['direct_S_error'] for r in packets),
                         same_moment_midpoint_gap=result['moment_counterexample']['exact_positive_gap'],
                         physical_common_tail_gap=result['moment_counterexample']['common_UV_tail'][0]['physical_midpoint_loss'],
                         output=str(output)), indent=2))


if __name__ == '__main__':
    main()
