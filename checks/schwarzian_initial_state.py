"""S072: exact disk-Schwarzian initial moments and preparation marginals.

Independent Gamma spectral, Bessel matrix-element and length-wavefunction
checks. These floating-point checks supplement the ordinary bounds; they
are not interval certificates or a quantum capture theorem.
"""
from pathlib import Path
import hashlib
import json
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.special import loggamma, gammaln

ROOT = Path(__file__).resolve().parents[1]


def moments(a):
    a = mp.mpf(a)
    s0 = 2*mp.quad(lambda v: mp.exp(-a*v*v)/mp.cosh(mp.pi*v), [0, 1, 4, mp.inf])
    s2 = 2*mp.quad(lambda v: v*v*mp.exp(-a*v*v)/mp.cosh(mp.pi*v), [0, 1, 4, mp.inf])
    w = lambda v: 2/mp.pi if not v else 2*v/mp.sinh(mp.pi*v)
    w0 = 2*mp.quad(lambda v: w(v)*mp.exp(-a*v*v), [0, 1, 4, mp.inf])
    w2 = 2*mp.quad(lambda v: v*v*w(v)*mp.exp(-a*v*v), [0, 1, 4, mp.inf])
    half = (1+2*a/mp.pi**2)*s0-4*a*a*s2/mp.pi**2
    one = (1+6*a/mp.pi**2)*w0-4*a*a*w2/mp.pi**2
    p2 = 3/(8*a)+mp.pi**2*(1-one)/(16*a*a)
    return half, one, p2


def source_spectral(a, delta, count):
    # Direct source (4.10), with h=4a, 2C=1, divided by (pi/h)^(2Delta).
    # No u,v factorization or closed midpoint formula enters this check.
    extent = np.pi/(4*a)+12/np.sqrt(4*a)
    nodes, weights = np.polynomial.legendre.leggauss(count)
    k, w = extent*(nodes+1)/2, extent*weights/2
    log_mu = np.log(k)+2*np.pi*k+np.log1p(-np.exp(-4*np.pi*k))
    k1, k2 = k[:, None], k[None, :]
    log_vertex = 2*np.real(loggamma(delta+1j*(k1+k2))
                         +loggamma(delta+1j*(k1-k2)))-gammaln(2*delta)
    log_z = 1.5*np.log(np.pi/(4*a))+np.pi**2/(4*a)
    log_integrand = (log_mu[:, None]+log_mu[None, :]-2*a*(k1*k1+k2*k2)
                     +log_vertex-np.log(2*np.pi**2)-log_z
                     -2*delta*np.log(np.pi/(4*a)))
    return float(w@np.exp(log_integrand)@w)


def length_wavefunction(a, count):
    # Insert K_(2ik)(z)=int exp[-z cosh(t)] cos(2kt) dt in the HH wavefunction.
    # Perform the Gaussian k integral first, then quadrature over t and length.
    # This independently checks normalization, both moments, and int |psi'|^2.
    nodes, weights = np.polynomial.legendre.leggauss(count)
    ell, well = 26*nodes+14, 26*weights  # [-12,40]
    tmax = 16*np.sqrt(a)
    t, wt = tmax*(nodes+1)/2, tmax*weights/2
    log_z = 1.5*np.log(np.pi/(4*a))+np.pi**2/(4*a)
    pref = np.sqrt(2)/np.pi*np.sqrt(np.pi/(2*a))*np.exp(np.pi**2/(2*a)-log_z/2)/(2*a)
    kernel = pref*np.exp(-t*t/(2*a))*(np.pi*np.cos(np.pi*t/a)-t*np.sin(np.pi*t/a))
    z = 4*np.exp(-ell/2)
    zcosh = z[:, None]*np.cosh(t)[None, :]
    ex = np.exp(-zcosh)
    psi = ex@(wt*kernel)
    derivative = (zcosh*ex/2)@(wt*kernel)
    observable = z*2*a/np.pi
    return dict(a=a, nodes=count, length_interval=[-12, 40], t_max=float(tmax),
                norm=float(well@(psi*psi)),
                half=float(well@(observable*psi*psi)),
                one=float(well@(observable**2*psi*psi)),
                p2=float(well@(derivative*derivative)),
                endpoint_amplitude_max=float(max(abs(psi[0]), abs(psi[-1]))))


@mp.workdps(40)
def bessel_checks():
    rows = []
    for dd, kk, jj in [('.5', '.25', '.7'), ('1', '.8', '1.3'), ('1.2', '.4', '1.1')]:
        d, k, j = map(mp.mpf, (dd, kk, jj))
        integrand = lambda z: z**(2*d-1)*mp.besselk(2j*k, z)*mp.besselk(2j*j, z)
        direct = mp.quad(integrand, [0, mp.mpf('.1'), 1, 8, 30, mp.inf])
        exact = 2**(2*d-3)*abs(mp.gamma(d+1j*(k+j))*mp.gamma(d+1j*(k-j)))**2/mp.gamma(2*d)
        relative = abs(direct/exact-1)
        assert relative < mp.mpf('1e-30')
        ell = mp.mpf('.3')
        f = lambda x: mp.besselk(2j*k, 4*mp.exp(-x/2))
        residual = abs(-mp.diff(f, ell, 2)+4*mp.exp(-ell)*f(ell)-k*k*f(ell))
        assert residual < mp.mpf('1e-35')
        rows.append(dict(delta=dd, k=kk, j=jj, matrix_relative_error=float(relative),
                         schrodinger_residual=float(residual)))
        print('Bessel case', dd, kk, jj, 'PASS', flush=True)
    return rows


def small_a_checks():
    rows = []
    gamma = mp.mpf(1)/4-2/mp.pi**2
    c1 = 6/mp.pi**2-mp.mpf(1)/2
    for astr in ['1e-12', '1e-8', '.0001', '.01', '.05', '.1']:
        a = mp.mpf(astr)
        half, one, p2 = moments(a)
        rh, r1 = (half-1+gamma*a)/(a*a), (one-1-c1*a)/(a*a)
        centered = one-2*half+1
        assert -3/(2*mp.pi**2) <= rh <= mp.mpf(5)/32*(1+2*a/mp.pi**2)
        assert abs(rh) <= mp.mpf(1)/5
        assert -5/mp.pi**2 <= r1 <= mp.mpf(1)/2*(1+6*a/mp.pi**2)
        assert abs(r1) <= mp.mpf(3)/5
        assert 0 < 1-half <= gamma*a+a*a/5
        assert 0 < centered <= 2*a/mp.pi**2+a*a
        assert 0 < p2 <= mp.pi**2/(32*a)+mp.mpf(5)/16
        rows.append(dict(a=astr, half=mp.nstr(half, 45), one=mp.nstr(one, 45),
                         half_remainder_over_a2=mp.nstr(rh, 30),
                         one_remainder_over_a2=mp.nstr(r1, 30),
                         nominal_centered_second=mp.nstr(centered, 30),
                         p2=mp.nstr(p2, 30)))
    return rows


def candidate_checks():
    rows = []
    gamma = mp.mpf(1)/4-2/mp.pi**2
    zc, width, speed = mp.mpf(53)/200, mp.mpf(3)/200, mp.mpf(1)/2000
    for k in [2000, 4000, 10000, 1000000]:
        M = mp.mpf(2)*k**3  # N eta, with N=2 k^4, eta=1/k
        a = mp.pi/(4*M*zc)
        half, one, p2 = moments(a)
        displacement_bound = zc**2*(2*a/mp.pi**2+a*a)
        velocity_bound = (mp.pi**2/(32*a)+mp.mpf(5)/16)/(M*M)
        assert 0 < zc**2*(one-2*half+1) <= displacement_bound
        assert 0 < p2/(M*M) <= velocity_bound
        # Separate Born marginals. Their sum is NOT a classical joint capture probability.
        ztail, vtail = displacement_bound/(width*width), velocity_bound/(speed*speed)
        assert ztail < mp.mpf('1.172e-8') and vtail < mp.mpf('2.602e-5')
        all_interval_a = mp.pi/(4*M*(mp.mpf(1)/4))
        all_midpoint_loss = gamma*all_interval_a+all_interval_a**2/5
        assert all_midpoint_loss < mp.mpf('9.3e-12')
        # Robertson bound with [z,v]=i z/(2M).
        variance_z = zc**2*(one-half*half)
        assert variance_z*p2/(M*M) >= (zc*half)**2/(16*M*M)
        rows.append(dict(k=k, M=str(2*k**3), nominal_z=str(zc), a=mp.nstr(a, 30),
                         mean_relative_loss=mp.nstr(1-half, 30),
                         loss_upper_all_nominal_z_in_capture_interval=mp.nstr(all_midpoint_loss, 30),
                         displacement_second_bound=mp.nstr(displacement_bound, 30),
                         velocity_second_bound=mp.nstr(velocity_bound, 30),
                         separate_redshift_tail_bound=mp.nstr(ztail, 30),
                         separate_velocity_tail_bound=mp.nstr(vtail, 30)))
    return rows


def main():
    mp.mp.dps = 70
    bessel = bessel_checks()
    print('Bessel matrix elements and canonical eigen-equation: PASS', flush=True)
    small = small_a_checks()
    spectral = []
    for a in [.05, .1, .5, 2.]:
        exact = moments(str(a))
        for delta, answer in zip([.5, 1.], exact[:2]):
            coarse, fine = (source_spectral(a, delta, n) for n in [640, 960])
            error, refinement = abs(fine-float(answer)), abs(fine-coarse)
            assert error < 3e-10 and refinement < 3e-10, (a, delta, error, refinement)
            spectral.append(dict(a=a, delta=delta, coarse=coarse, fine=fine,
                                 exact=float(answer), absolute_error=error, refinement=refinement))
    print('Direct two-energy Gamma integrals and small-a bounds: PASS', flush=True)
    lengths = []
    for a in [1., 2., 4.]:
        expected = [1., *map(float, moments(str(a))) ]
        coarse, fine = (length_wavefunction(a, n) for n in [256, 512])
        keys = ['norm', 'half', 'one', 'p2']
        error = max(abs(fine[key]-val) for key, val in zip(keys, expected))
        refinement = max(abs(fine[key]-coarse[key]) for key in keys)
        assert error < 2e-9 and refinement < 2e-9, (a, error, refinement)
        fine.update(max_absolute_error=error, refinement=refinement)
        lengths.append(fine)
    print('Length wavefunction norm, first/second moment and momentum: PASS', flush=True)
    candidates = candidate_checks()
    result = dict(date='2026-10-13', status='PASS', python=platform.python_version(),
                  mpmath=mp.__version__, numpy=np.__version__, scipy=scipy.__version__,
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  precision_decimal_digits=mp.mp.dps, bessel_decimal_digits=40,
                  bessel=bessel, small_a=small, direct_spectral=spectral,
                  length_wavefunction=lengths, candidates=candidates,
                  numerical_scope='Finite floating-point quadratures with refinement, not interval certification.',
                  physical_scope=['Disk Schwarzian/pure two-boundary JT initial Hartle-Hawking state.',
                                  'Separate initial redshift and velocity marginals; no joint classical probability.',
                                  'No driven quantum capture, matter determinant, charged horizon gluing or topology preparation.'])
    target = ROOT/'evidence/schwarzian-initial-state.json'
    target.write_bytes((json.dumps(result, indent=2)+'\n').encode('utf8'))
    print('PASS: exact pure-JT initial-state moments and separate preparation marginals.')
    print('Worst candidate:', json.dumps(candidates[0]))
    print('Maximum Gamma spectral error:', max(r['absolute_error'] for r in spectral))
    print('Maximum length-wavefunction error:', max(r['max_absolute_error'] for r in lengths))


if __name__ == '__main__':
    main()
