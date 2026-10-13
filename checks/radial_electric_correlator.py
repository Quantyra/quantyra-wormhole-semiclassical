"""S072 radial Gauss profile, charge lattice and boundary thermal correlator.

Independent canonical traces and radial/spectral representations test the
ordinary proof. This is not a four-dimensional driven simulation.
"""
from pathlib import Path
import hashlib
import itertools
import json
import math

import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import eigh, solve
from scipy.integrate import solve_ivp
from scipy.special import polygamma

ROOT = Path(__file__).resolve().parents[1]


def lattice(n, c, x, t, cutoff=10, dual=False):
    p = np.array(list(itertools.product(range(-cutoff, cutoff+1), repeat=n)), float)
    total = p.sum(axis=1)
    norm = (p*p).sum(axis=1)
    if dual:
        energy = np.pi**2*((norm-total**2/n)/c+total**2/(n*n*x))
        weights = np.exp(-energy)
        return float(weights@np.cos(2*np.pi*t*p[:, 0])/sum(weights))
    energy = c*(norm-total**2/n)+x*total**2
    shifted = p.copy()
    shifted[:, 0] += t
    q = shifted.sum(axis=1)
    shifted_energy = c*((shifted*shifted).sum(axis=1)-q*q/n)+x*q*q
    return float(np.exp(-shifted_energy).sum()/np.exp(-energy).sum())


def rotor_theta(x, t):
    k = np.arange(-100, 101, dtype=float)
    return float(np.exp(-x*(k+t)**2).sum()/np.exp(-x*k*k).sum())


def lattice_checks():
    rows = []
    for n in (1, 2, 3):
        for c, x in ((.4, .8), (1.3, 2.1), (3., 4.)):
            for t in (.23, .5):
                direct = lattice(n, c, x, t, 12)
                refined = lattice(n, c, x, t, 16)
                dual = lattice(n, c, x, t, 8, True)
                assert abs(direct-refined) < 2.e-12
                assert abs(refined-dual) < 3.e-12
                rows.append(dict(n=n, c=c, x=x, t=t, direct=refined,
                                 dual=dual, error=abs(refined-dual)))
    limits = []
    # Includes a strong compact effect so a missing theta factor is visible.
    for n in (1, 2, 3):
        for x in (.8, 4.):
            target = rotor_theta(x, .5)
            values = [lattice(n, c, x, .5, 9, True) for c in (2., 1., .3, .15)]
            assert abs(values[-1]-target) < 4.e-9
            assert abs(values[-1]-target) <= abs(values[0]-target)+3.e-15
            limits.append(dict(n=n, x=x, c_values=[2., 1., .3, .15],
                               values=values, limit=target))
    return dict(finite_lattices=rows, optical_limits=limits)


def oscillator_factor(v, t):
    return (-np.expm1(-v*t))*(-np.expm1(-v*(1-t)))/(-np.expm1(-v))


def canonical_trace(n, c, x, beta, omega, displacement, t, dim):
    """Thermal matrices of H_Q in a common, unshifted oscillator basis.

    The boundary ladder is the identity on that coordinate, not the
    identity on the differently centered oscillator eigenstates.
    """
    a = np.diag(np.sqrt(np.arange(1, dim)), 1)
    coordinate = (a+a.T)/np.sqrt(2*omega)
    h0 = np.diag(omega*(np.arange(dim)+.5))
    vectors = list(itertools.product(range(-6, 7), repeat=n))
    totals = range(-6*n, 6*n+2)
    factors, traces = {}, {}
    for Q in totals:
        h = h0-omega**2*displacement*Q*coordinate
        h = h+.5*omega**2*(displacement*Q)**2*np.eye(dim)
        e, u = eigh(h)
        factors[Q] = ((u*np.exp(-beta*(1-t)*e))@u.T,
                      (u*np.exp(-beta*t*e))@u.T)
        traces[Q] = float(np.exp(-beta*e).sum())
    denom, numer = 0., 0.
    charge_denom, charge_numer = 0., 0.
    for vec in vectors:
        v = np.asarray(vec, float)
        Q = int(v.sum())
        energy = c*(v@v-Q*Q/n)+x*Q*Q
        plus = v.copy()
        plus[0] += 1
        energy_plus = c*(plus@plus-(Q+1)**2/n)+x*(Q+1)**2
        weight = math.exp(-(1-t)*energy-t*energy_plus)
        denom += math.exp(-energy)*traces[Q]
        # trace(M N), using their symmetry, avoids extra matrix products.
        numer += weight*np.sum(factors[Q][0]*factors[Q+1][1])
        charge_denom += math.exp(-energy)
        charge_numer += weight
    direct = numer/denom
    osc = math.exp(-displacement**2*omega*oscillator_factor(beta*omega, t)/2)
    formula = charge_numer/charge_denom*osc
    return direct, formula


def trace_checks():
    rows = []
    for n, c, x in ((1, .8, .9), (2, .9, .7)):
        for t in (.23, .5):
            a, pred = canonical_trace(n, c, x, 1.3, 1.7, .37, t, 48)
            b, _ = canonical_trace(n, c, x, 1.3, 1.7, .37, t, 64)
            assert abs(a-b) < 2.e-11 and abs(b-pred) < 2.e-11
            charge = math.exp(-t*(1-t)*(c*(1-1/n)+x))*lattice(n, c, x, t, 12)
            expected = charge*math.exp(-.37**2*1.7*oscillator_factor(1.3*1.7, t)/2)
            assert abs(b-expected) < 3.e-11
            rows.append(dict(n=n, t=t, matrix_trace=b, prediction=expected,
                             error=abs(b-expected), refinement=abs(a-b)))
    return rows


def constant_spectrum():
    rows = []
    for length, mass, beta in ((.8, .6, 1.2), (2.1, 1.7, .7), (4., .25, 2.)):
        lam = 3/np.pi
        K0 = mass/np.tanh(mass*length)
        k = np.pi*np.arange(1, 160001)/length
        w = np.sqrt(k*k+mass*mass)
        j = np.arange(1, 160001, dtype=float)
        s = 2*np.pi*j/beta
        sv = np.sqrt(s*s+mass*mass)
        diff = sv/np.tanh(sv*length)-s/np.tanh(s*length)
        for t in (.19, .37, .5):
            osc = np.sum(k*k/(length*lam*w**3)*oscillator_factor(beta*w, t)
                         -1/(length*lam*k)*oscillator_factor(beta*k, t))
            canonical = beta*(K0-1/length)/(2*lam)*t*(1-t)+osc
            mats = 2/(beta*lam)*np.sum((1-np.cos(2*np.pi*j*t))*diff/s**2)
            # Separate convergent 1/k^3 tails; tolerance exceeds both tails.
            assert abs(canonical-mats) < 2.e-10
            assert 0 <= mats <= beta*(K0-1/length)/(2*lam)*t*(1-t)
            rows.append(dict(length=length, mass=mass, beta=beta, t=t,
                             canonical=canonical, matsubara=mats,
                             error=abs(canonical-mats)))
    return rows


def fem(length, count, potential):
    nodes = np.linspace(-length, 0, count+1)
    h = length/count
    mass, action = np.zeros((count+1, count+1)), np.zeros((count+1, count+1))
    points, weights = np.polynomial.legendre.leggauss(6)
    shape = np.stack(((1-points)/2, (1+points)/2))
    for j in range(count):
        local_m = h*np.array([[2., 1.], [1., 2.]])/6
        local_a = np.array([[1., -1.], [-1., 1.]])/h
        xx = (nodes[j]+nodes[j+1])/2+h*points/2
        local_a += (shape*(h/2*weights*potential(xx)))@shape.T
        mass[j:j+2, j:j+2] += local_m
        action[j:j+2, j:j+2] += local_a
    # Zero inner value, leave outer value in the matrix.
    mass, action = mass[1:, 1:], action[1:, 1:]
    z = np.r_[-solve(action[:-1, :-1], action[:-1, -1], assume_a='pos'), 1.]
    K0 = float(z@action@z)
    vals, vecs = eigh(action[:-1, :-1], mass[:-1, :-1])
    coeffs = vecs.T@(mass@z)[:-1]
    contact = float(z@mass@z-coeffs@coeffs)
    return mass, action, K0, vals, coeffs, contact


def radial_checks():
    rows = []
    for label, potential in (
        ('localized', lambda x: .8*np.exp(1.2*x)*(1+.2*np.sin(2*x)**2)),
        ('nonconstant', lambda x: .3+.15*(x+1.5)**2)):
        length, beta, lam = 3., 1.7, 3/np.pi
        previous = None
        for count in (40, 80, 160):
            current = fem(length, count, potential)
            free = fem(length, count, lambda x: np.zeros_like(x))
            M, H, K0, vals, coeffs, contact = current
            M0, H0, k00, vals0, coeffs0, contact0 = free
            assert abs(contact-contact0) < 2.e-13 and contact > 0
            errors = []
            diffs = []
            continuum_errors = []
            for s in (0., .1, .8, 3., 11.):
                A = H+s*s*M
                primal = A[-1, -1]-A[-1, :-1]@solve(A[:-1, :-1], A[:-1, -1], assume_a='pos')
                spectral = K0+s*s*contact+np.sum(s*s*vals*coeffs**2/(s*s+vals))
                assert abs(primal-spectral) < 3.e-11
                errors.append(float(abs(primal-spectral)))
                # Independently integrate continuum Riccati from a small
                # positive distance after the inner Dirichlet singularity.
                eps = 1.e-6
                initial = 1/eps+(potential(-length)+s*s)*eps/3
                ode = solve_ivp(lambda xx, y: [potential(xx)+s*s-y[0]**2],
                                (-length+eps, 0), [initial], rtol=2.e-10, atol=2.e-11)
                assert ode.success
                exact = float(ode.y[0, -1])
                continuum_errors.append(float(abs(primal-exact)))
                free_exact = 1/length if s == 0 else s/np.tanh(s*length)
                diffs.append(exact-free_exact)
            assert min(diffs) > 0 and np.max(np.diff(diffs)) < 2.e-10
            # Spectral formula versus a separately evaluated Matsubara sum.
            j = np.arange(1, 20002, 2, dtype=float)
            s = 2*np.pi*j/beta
            dk = np.full_like(s, K0-k00)
            for start in range(0, len(s), 500):
                ss = s[start:start+500, None]**2
                dk[start:start+500] += np.sum(ss*vals*coeffs**2/(ss+vals), axis=1)
                dk[start:start+500] -= np.sum(ss*vals0*coeffs0**2/(ss+vals0), axis=1)
            mats = 4/(beta*lam)*np.sum(dk/s**2)
            # A finite mesh has a nonzero constant delta-K at high s.
            # Sum that tail exactly, then bound the remaining 1/s^4 tail
            # by the two separate positive oscillator measures.
            dk_inf = K0-k00+np.sum(vals*coeffs**2)-np.sum(vals0*coeffs0**2)
            next_odd = j[-1]+2
            tail2 = float(polygamma(1, next_odd/2))/4
            tail4 = float(polygamma(3, next_odd/2))/96
            mats += beta*dk_inf/(np.pi**2*lam)*tail2
            remainder = beta**3/(4*np.pi**4*lam)*tail4*(
                np.sum(vals**2*coeffs**2)+np.sum(vals0**2*coeffs0**2))
            w, w0 = np.sqrt(vals), np.sqrt(vals0)
            osc = np.sum(w*coeffs**2*np.tanh(beta*w/4)-w0*coeffs0**2*np.tanh(beta*w0/4))/(2*lam)
            canonical = beta*(K0-k00)/(8*lam)+osc
            assert abs(canonical-mats) <= remainder+3.e-11
            rows.append(dict(profile=label, elements=count, K0=K0,
                             mesh_contact=contact, contact_difference=abs(contact-contact0),
                             DtN_spectral_error=max(errors), continuum_delta_K=diffs,
                             DtN_continuum_error=max(continuum_errors),
                             midpoint=canonical, matsubara=mats,
                             thermal_error=abs(canonical-mats),
                             matsubara_remainder_bound=float(remainder),
                             refinement=None if previous is None else abs(canonical-previous)))
            previous = canonical
        last = rows[-3:]
        assert last[2]['refinement'] < .4*last[1]['refinement']
        assert last[2]['DtN_continuum_error'] < .4*last[1]['DtN_continuum_error']
    return rows


def candidate_checks():
    mp.mp.dps = 60
    rows = []
    for k in (2000, 10000, 1000000):
        x_bound = mp.mpf(1)/(12*k)
        gaussian = -mp.expm1(-x_bound/4)
        log_compact = mp.log(4)-mp.pi**2/x_bound-mp.log1p(-mp.exp(-8*mp.pi**2/x_bound))
        assert gaussian < mp.mpf(1)/(48*k)
        rows.append(dict(k=k, x_upper=str(x_bound), loss_main=str(gaussian),
                         log_compact_correction_upper=str(log_compact)))
    return rows


def main():
    result = dict(date='2026-10-13', status='PASS',
                  scope='Fixed-background spherical electric LLL neutral-inner-cut Gibbs ensemble',
                  lattices=lattice_checks(), shifted_thermal_traces=trace_checks(),
                  constant_potential_spectra=constant_spectrum(),
                  radial_profile_checks=radial_checks(), candidate=candidate_checks(),
                  limitations=['No horizon black-hole degeneracy matching',
                               'No angular or massive interactions',
                               'No driven stress or geometry',
                               'No physical creation or full stability claim'],
                  versions=dict(numpy=np.__version__, scipy=scipy.__version__, mpmath=mp.__version__))
    result['source_hashes'] = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
                             ['notes/radial-electric-correlator.md', 'checks/radial_electric_correlator.py']}
    path = ROOT/'evidence/radial-electric-correlator.json'
    path.write_bytes((json.dumps(result, indent=2)+'\n').encode('utf-8'))
    print('PASS: 18 charge lattices, 6 optical limits, 4 shifted thermal traces,')
    print('9 independent constant-potential comparisons, 6 nonconstant radial meshes, 3 candidate bounds')
    print('maximum canonical trace error', max(r['error'] for r in result['shifted_thermal_traces']))
    print('maximum continuum spectral error', max(r['error'] for r in result['constant_potential_spectra']))
    print('maximum mesh thermal error', max(r['thermal_error'] for r in result['radial_profile_checks']))


if __name__ == '__main__':
    main()
