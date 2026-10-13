"""Fock traces, group integration and reservoir limit for S072.

Only the stated global Casimir Hamiltonian is tested. No radial angular
interaction, gravitational gluing, driven stress or stability is certified.
"""
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
import scipy
from scipy.linalg import eigh
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'checks'))
from fermionic_compact_transfer import operators


def su_heat(cos_half_angle, a, time):
    """Full SU(2) cover kernel; spin l=n/2, character U_n(cos(theta/2))."""
    c = np.asarray(cos_half_angle)
    nmax = int(math.ceil(2 * math.sqrt(55 / (a * time)))) + 6
    before = np.ones_like(c)
    answer = before.copy()
    current = 2 * c
    for n in range(1, nmax + 1):
        answer += (n + 1) * math.exp(-time * a * n * (n + 2) / 4) * current
        before, current = current, 2 * c * current - before
    return answer


def electric_kernel(alpha, x, time=1):
    nmax = int(math.ceil(math.sqrt(55 / (x * time)))) + 3
    n = np.arange(-nmax, nmax + 1)
    return (np.exp(-time * x * n * n)[:, None]
            * np.exp(1j * n[:, None] * alpha)).sum(axis=0).real


def electric_insert(alpha, x, t):
    nmax = int(math.ceil(math.sqrt(55 / (x * min(t, 1-t))))) + 3
    n = np.arange(-nmax, nmax + 1)
    return (np.exp(-x * ((1-t) * (n-1)**2 + t * n*n))[:, None]
            * np.exp(1j * n[:, None] * alpha)).sum(axis=0)


def angular_insert(q, a, t, theta_g, order):
    """Integrate h angle and axis; no Clebsch-Gordan formula is used."""
    nodes, weights = leggauss(order)
    theta_h = math.pi * (nodes + 1)
    wh = weights * np.sin(theta_h / 2)**2
    axis = nodes[None, :]
    ch, sh = np.cos(theta_h[:, None]/2), np.sin(theta_h[:, None]/2)
    aa, abar = ch + 1j*sh*axis, ch - 1j*sh*axis
    minus_bb = -sh*sh*(1-axis*axis)
    diagonals = []
    for p in range(q):
        s = q-1-p
        diagonals.append(sum(math.comb(p, n)*math.comb(s, n)
                             * aa**(p-n)*abar**(s-n)*minus_bb**n
                             for n in range(min(p, s)+1)))
    common = wh[:, None] * (weights[None, :]/2)
    common *= su_heat(ch, a, 1-t)
    result = np.empty((len(theta_g), q), dtype=complex)
    for i, angle in enumerate(theta_g):
        cos_product = ch*math.cos(angle/2) + sh*math.sin(angle/2)*axis
        weight = common * su_heat(cos_product, a, t)
        result[i] = [np.sum(weight * diagonal) for diagonal in diagonals]
    return result


def group_ratio(q, energies, x, a, t=.5, tag=0, order=48, nalpha=96,
                identical_reservoir=None):
    """Exact group-integral formula, evaluated on the U(1)xSU(2) cover.

    The fermion determinant imposes the center quotient automatically.
    identical_reservoir=(energy, count) avoids constructing a large Fock space.
    """
    nodes, weights = leggauss(order)
    theta = math.pi*(nodes+1)
    wtheta = weights*np.sin(theta/2)**2
    alpha = 2*math.pi*np.arange(nalpha)/nalpha
    m = np.arange(q)-(q-1)/2
    eigen = np.exp(1j*(alpha[:, None, None]+theta[None, :, None]*m))
    z = math.exp(-energies[tag])
    A = (1+z*eigen)/(1+z)
    f = np.prod(abs(A)**2, axis=2)
    # det(A)^* adj(A): remains regular at z=1 and det(A)=0.
    B = np.empty_like(A)
    det_conj = np.prod(A.conj(), axis=2)
    for s in range(q):
        B[:, :, s] = det_conj*np.prod(np.delete(A, s, axis=2), axis=2)
    reservoir = np.ones_like(f)
    for s, u in enumerate(energies):
        if s != tag:
            zz = math.exp(-u)
            reservoir *= np.prod(abs((1+zz*eigen)/(1+zz))**2, axis=2)
    if identical_reservoir is not None:
        u, count = identical_reservoir
        zz = math.exp(-u)
        reservoir *= np.prod(abs((1+zz*eigen)/(1+zz))**2, axis=2)**count
    ks = su_heat(np.cos(theta/2), a, 1)
    ku = electric_kernel(alpha, x)
    denominator = np.sum(reservoir*f*ku[:, None]*ks[None, :]*wtheta)/nalpha
    angular = angular_insert(q, a, t, theta, max(48, min(order, 88)))
    eu = electric_insert(alpha, x, t)
    numerator = np.sum(reservoir*eu[:, None]
                       * np.sum(B*angular[None, :, :], axis=2)/q*wtheta)/nalpha
    assert abs(numerator.imag) < 2e-10*max(1, abs(numerator.real))
    return float(numerator.real/denominator), float(denominator)


def fock_ratio(q, energies, x, a, t=.5, tag=0):
    copies = len(energies)
    D, basis, Q, Js, J2 = operators(q, copies)
    u_modes = np.tile(np.repeat(energies, q), 2)
    kinetic = np.array([sum(u_modes[s] for s in range(D) if state >> s & 1)
                        for state in basis])
    H = (x*(Q@Q)+a*J2).toarray()+np.diag(kinetic)
    energy, vec = eigh(H)
    Z = np.exp(-energy).sum()
    weights = np.exp(-t*energy[:, None]-(1-t)*energy[None, :])/Z
    b = 0.
    for s in range(tag*q, (tag+1)*q):
        c = np.zeros_like(H)
        for state in basis:
            if state >> s & 1:
                c[state ^ (1 << s), state] = (-1)**((state & ((1 << s)-1)).bit_count())
        transformed = vec.conj().T@c@vec
        b += np.sum(weights*abs(transformed)**2).real/q
    z = math.exp(-energies[tag])
    Zfree = math.prod((1+math.exp(-u))**(2*q) for u in energies)
    return float(b/(z**(1-t)/(1+z))), float(Z/Zfree)


def allowed(q, Q, twice_l):
    return twice_l % 2 == ((q-1)*Q) % 2


def compact_ratio(q, x, a, t=.5):
    """Independent allowed-representation/fusion sum on the faithful group."""
    qmax = int(math.ceil(math.sqrt(65/(x*min(t, 1-t)))))+q+3
    lmax = int(math.ceil(2*math.sqrt(65/(a*min(t, 1-t)))))+q+3
    Z = numerator = 0.
    for Q in range(-qmax, qmax+1):
        for l2 in range(lmax+1):
            if not allowed(q, Q, l2):
                continue
            E = x*Q*Q+a*l2*(l2+2)/4
            Z += (l2+1)**2*math.exp(-E)
            for lp2 in range(abs(l2-(q-1)), l2+q, 2):
                assert allowed(q, Q+1, lp2)
                Ep = x*(Q+1)**2+a*lp2*(lp2+2)/4
                numerator += (l2+1)*(lp2+1)*math.exp(-(1-t)*E-t*Ep)
    return numerator/(q*Z)


def heat_constants(q, x, a, t):
    qmax = int(math.ceil(math.sqrt(80/(x*t))))+3
    lmax = int(math.ceil(2*math.sqrt(80/(a*t))))+6
    P = L = 0.
    j = (q-1)/2
    for Q in range(-qmax, qmax+1):
        for l2 in range(lmax+1):
            if allowed(q, Q, l2):
                term = (l2+1)**2*math.exp(-t*(x*Q*Q+a*l2*(l2+2)/4))
                P += term
                L += 2*(abs(Q)+l2/(2*j))*term
    return P, L


def certificate(q, x, a, V, t=.5):
    S = q*(q*q-1)/12
    inv_c = 48*math.e*math.pi**4*math.sqrt(q)*S**1.5
    eps = math.sqrt((math.log(inv_c)+4*math.log(V))/V)
    assert eps < .25
    tail = min(1., inv_c*V*V*math.exp(-V*eps*eps))
    P, L1 = heat_constants(q, x, a, 1)
    _, Lt = heat_constants(q, x, a, t)
    U = eps*(Lt/P+(2*q-1)/2)+2*tail
    W = eps*L1/P+eps*eps/4+tail
    error = (U+W)/(1-W) if W < 1 else None
    return dict(q=q, x=x, a=a, t=t, V=V, epsilon=eps,
                tail_bound=tail, numerator_error=U, denominator_error=W,
                absolute_ratio_error=error)


def algebra_checks():
    rng = np.random.default_rng(2072)
    maximum = 0.
    for q in [2, 3, 4, 9]:
        m = np.arange(q)-(q-1)/2
        for z in [0., .15, .999999, 1.]:
            for close in [False, True]:
                alpha, theta = rng.uniform(-.01, .01, 2) if close else rng.uniform(-math.pi, math.pi, 2)
                r = np.exp(1j*(alpha+m*theta))
                aa = (1+z*r)/(1+z)
                ff = np.prod(abs(aa)**2)
                B = np.array([np.prod(aa.conj())*np.prod(np.delete(aa, n)) for n in range(q)])
                assert np.max(abs(B)) <= 1+2e-14
                assert 1-ff <= .25*np.sum(abs(r-1)**2)+2e-14
                assert np.max(abs(B-1)) <= (2*q-1)/2*np.max(abs(r-1))+2e-14
                if min(abs(aa)) > 1e-10:
                    maximum = max(maximum, float(np.max(abs(B-ff/aa))))
    assert maximum < 3e-14
    return dict(cases=32, cofactor_identity_error=maximum,
                zero_energy_singular_inverse_avoided=True)


def main():
    rows = []
    cases = [(2, [0.], .3, .4, .5),
             (2, [.7, 1.4], .2, .25, .2),
             (2, [.7, 1.4], .2, .25, .8),
             (3, [.4], .3, .2, .5),
             (3, [0.], .15, .3, .2),
             (4, [.9], .2, .15, .5)]
    for q, energies, x, a, t in cases:
        exact, Z = fock_ratio(q, energies, x, a, t)
        group, Zg = group_ratio(q, energies, x, a, t)
        err = max(abs(exact-group), abs(Z-Zg))
        assert err < 3e-9, (q, energies, t, exact, group, Z, Zg)
        row = dict(q=q, energies=energies, x=x, a=a, t=t,
                   fock_ratio=exact, group_ratio=group, error=err)
        rows.append(row)
        print('fock/group', row, flush=True)
    limits = []
    for q, t in [(2, .2), (2, .5), (3, .5)]:
        x, a, u = .3, .25, .8
        target = compact_ratio(q, x, a, t)
        previous = float('inf')
        for count in [16, 64, 256]:
            value, _ = group_ratio(q, [u], x, a, t, order=144, nalpha=384,
                                   identical_reservoir=(.5, count))
            error = abs(value-target)
            assert error < previous
            previous = error
            row = dict(q=q, t=t, remaining_multiplets=count,
                       ratio=value, compact_limit=target, error=error)
            limits.append(row)
            print('reservoir', row, flush=True)
        refined, _ = group_ratio(q, [u], x, a, t, order=192, nalpha=512,
                                  identical_reservoir=(.5, 256))
        assert abs(refined-value) < 3e-10, (refined, value)
        limits[-1]['quadrature_refinement_error'] = abs(refined-value)
        if t == .5:
            lam = x+a*(q*q-1)/4
            assert 0 <= 1-target <= lam/4
    certs = [certificate(q, .3, .25, V, t)
             for q, t in [(2, .2), (2, .5), (3, .5)] for V in [1e8, 1e12, 1e16]]
    for row in certs:
        assert 0 < row['absolute_ratio_error'] < 1
        row['meets_two_percent_target'] = row['absolute_ratio_error'] < .02
        if row['V'] >= 1e12:
            assert row['meets_two_percent_target']
    # Retain the informative failure: the first trial volume did not meet
    # the requested error budget. Its bound was not invalidated or loosened.
    assert all(not row['meets_two_percent_target'] for row in certs if row['V'] == 1e8)
    receipt = dict(status='PASS', scope='Global electric/angular Casimir fermion reservoir only',
                   algebra=algebra_checks(), finite_fock_group=rows,
                   reservoir_convergence=limits, explicit_certificates=certs,
                   versions=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    path = ROOT/'evidence/angular-reservoir-correlator.json'
    path.write_bytes((json.dumps(receipt, indent=2)+'\n').encode())
    print('PASS', len(rows), 'Fock/group;', len(limits), 'reservoir;', len(certs), 'certificates', flush=True)


if __name__ == '__main__':
    main()
