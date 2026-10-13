"""S072 full charged-transition midpoint in a radial fermion reservoir.

This checks the global electric charging model, not a substitution of a
screened capacitance for the complete radial gauge Hamiltonian.
"""
from pathlib import Path
import hashlib
import json
import math

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
CSTAR = math.e*math.sqrt(2)*math.pi**3.5/32


def reservoir(us):
    probabilities = np.array([1.])
    for u in us:
        p = expit(-u)
        v = p*(1-p)
        probabilities = np.convolve(probabilities, [v, 1-2*v, v])
    charges = np.arange(len(probabilities))-(len(probabilities)-1)//2
    assert abs(sum(probabilities)-1) < 2.e-12
    assert np.max(abs(probabilities-probabilities[::-1])) < 2.e-14
    return charges, probabilities


def formula(us, tag, x):
    remaining = np.delete(np.array(us, dtype=float), tag)
    Q, P = reservoir(remaining)
    f = np.array([P@np.exp(-x*(Q+s)**2) for s in (0, .5, 1)])
    z = math.exp(-us[tag])
    D = (1+z*z)*f[0]+2*z*f[2]
    b = math.sqrt(z)*(1+z)*math.exp(-x/4)*f[1]/D
    b0 = math.sqrt(z)/(1+z)
    ratio = b/b0
    lower, upper = math.exp(-x/2), math.exp(x/4)/math.cosh(x/2)
    assert lower-2.e-12 <= ratio <= upper+2.e-12
    for s, value in zip((0, .5, 1), f):
        assert math.exp(-x*s*s)*f[0]-2.e-14 <= value <= f[0]+2.e-14
    frequencies, weights = [], []
    for a in (0, 1):
        frequencies.extend(us[tag]+x*(2*(Q-a)+1))
        weights.extend(P*z**a*(np.exp(-x*(Q-a)**2)
                       +z*np.exp(-x*(Q+1-a)**2))/D)
    frequencies, weights = np.array(frequencies), np.array(weights)
    assert min(weights) >= 0 and abs(sum(weights)-1) < 3.e-13
    # Stable evaluation of 1/(2 cosh(v/2)).
    kernel = np.exp(-abs(frequencies)/2)/(1+np.exp(-abs(frequencies)))
    spectral_b = float(weights@kernel)
    assert abs(spectral_b-b) < 3.e-13
    return dict(b=float(b), b0=b0, ratio=float(ratio), lower=lower, upper=upper,
                spectral_error=abs(spectral_b-b), f=f.tolist())


def direct_trace(us, x, seed_charge=0):
    modes = 2*len(us)
    states = np.arange(1 << modes)
    occupations = ((states[:, None] >> np.arange(modes)) & 1)
    charges = occupations[:, ::2].sum(axis=1)-occupations[:, 1::2].sum(axis=1)
    energies = occupations@np.repeat(us, 2)+x*(charges+seed_charge)**2
    shift = min(energies)
    boltzmann = np.exp(-(energies-shift))
    Z = sum(boltzmann)
    values = []
    for tag in range(len(us)):
        empty = states[occupations[:, 2*tag] == 0]
        filled = empty | (1 << (2*tag))
        value = sum(np.exp(-(energies[empty]+energies[filled]-2*shift)/2))/Z
        values.append(float(value))
    return values


def small_fock_checks():
    cases = []
    spectra = [[0.], [.2, 1.3], [0., .7, 4.], [.1, .5, 1., 3., 9.]]
    for us in spectra:
        for x in (.00004, .03, .4, 2., 6.):
            direct = direct_trace(us, x)
            errors, ratios = [], []
            for tag, value in enumerate(direct):
                answer = formula(us, tag, x)
                errors.append(abs(value-answer['b']))
                ratios.append(answer['ratio'])
            assert max(errors) < 3.e-13
            cases.append(dict(energies=us, x=x, fock_states=4**len(us),
                              maximum_trace_error=max(errors), ratios=ratios))
    # The paired neutral-seed hypothesis has real content.
    biased = direct_trace([1.], 1., seed_charge=4)[0]
    biased_ratio = biased/(1/(2*math.cosh(.5)))
    assert biased_ratio < math.exp(-.5)/5
    return cases, dict(seed_charge=4, x=1., u=1., ratio=biased_ratio,
                       inapplicable_neutral_lower=math.exp(-.5))


def fourier_checks():
    cases = []
    for us, x in [([.2, .7, 3.], .07), ([0., 1., 4., 8.], .8), ([2.], 3.)]:
        Q, P = reservoir(us)
        ps = expit(-np.array(us))
        vs = ps*(1-ps)
        def phi(alpha):
            return float(np.prod(1-4*vs*np.sin(alpha/2)**2))
        errors = []
        for s in (0., .5, 1., 1.7):
            cutoff = math.sqrt(4*x*65)
            value = quad(lambda alpha: math.exp(-alpha*alpha/(4*x))*phi(alpha)
                         *math.cos(s*alpha)/math.sqrt(math.pi*x),
                         0, cutoff, epsabs=2.e-13, epsrel=2.e-13, limit=300)[0]
            # Even transform: 2/sqrt(4 pi x) = 1/sqrt(pi x).
            expected = float(P@np.exp(-x*(Q+s)**2))
            errors.append(abs(value-expected))
        assert max(errors) < 2.e-12
        cases.append(dict(energies=us, x=x, maximum_fourier_error=max(errors)))
    return cases


def theta_data(x):
    mp.mp.dps = 65
    xx = mp.mpf(str(x))
    cutoff = int(mp.ceil(mp.sqrt(170/xx)))+2
    t0 = mp.fsum(mp.exp(-xx*n*n) for n in range(-cutoff, cutoff+1))
    th = mp.fsum(mp.exp(-xx*(mp.mpf(n)+mp.mpf('.5'))**2)
                 for n in range(-cutoff, cutoff+1))
    s0 = mp.fsum(n*n*mp.exp(-xx*n*n) for n in range(-cutoff, cutoff+1))
    sh = mp.fsum((mp.mpf(n)+mp.mpf('.5'))**2
                 *mp.exp(-xx*(mp.mpf(n)+mp.mpf('.5'))**2)
                 for n in range(-cutoff, cutoff+1))
    dual0 = mp.jtheta(3, 0, mp.exp(-mp.pi**2/xx))
    dualh = mp.jtheta(4, 0, mp.exp(-mp.pi**2/xx))
    error = abs(th/t0-dualh/dual0)
    assert error < mp.mpf('1e-55')
    tail_bound = 4*mp.exp(-mp.pi**2/xx)/(-mp.expm1(-8*mp.pi**2/xx))
    assert 1-th/t0 <= tail_bound+mp.mpf('1e-60')
    return dict(C=float(mp.exp(-xx/4)*th/t0), mean0=float(s0/t0),
                meanhalf=float(sh/th), poisson_error=float(error),
                compact_deficit_upper=float(tail_bound))


def concentration(us):
    ps = expit(-np.array(us))
    vs = ps*(1-ps)
    V = float(sum(vs))
    def phi(alpha):
        # log1p handles a broad spectrum and thermodynamic concentration.
        factors = 4*vs*np.sin(alpha/2)**2
        if np.any(factors >= 1):
            return 0.
        return float(np.exp(np.log1p(-factors).sum()))
    area = quad(phi, 0, math.pi, epsabs=1.e-13, epsrel=1.e-12)[0]
    second = quad(lambda alpha: alpha*alpha*phi(alpha), 0, math.pi,
                  epsabs=1.e-13, epsrel=1.e-12)[0]/area
    universal = min(math.pi**2, CSTAR/V) if V >= .5 else math.pi**2
    assert second <= universal+2.e-12
    Q, P = reservoir(us)
    ratio = P/P[len(us)]
    assert max(ratio) <= 1+2.e-12
    assert np.max(1-ratio-Q*Q*second/2) < 2.e-11
    return V, second, universal


def reservoir_limit_checks():
    cases = []
    for M in (8, 32, 128, 512, 2048):
        us = np.full(M, 1.)
        V, m2, universal = concentration(us)
        for x in (.05, .2, 1.5):
            td = theta_data(x)
            dn = m2/2*(td['meanhalf']+.25)
            dd = m2/2*(td['mean0']+.5)
            udn = universal/2*(td['meanhalf']+.25)
            udd = universal/2*(td['mean0']+.5)
            ratios = []
            for tag_u in (.2, 2., 8.):
                result = formula(np.r_[tag_u, us], 0, x)
                relative = result['ratio']/td['C']
                assert relative >= 1-dn-2.e-12
                if dd < 1:
                    assert relative <= 1/(1-dd)+2.e-12
                ratios.append(relative)
            cases.append(dict(reservoir_pairs=M, V=V, x=x,
                              measured_m2=m2, universal_m2=universal,
                              lower=1-dn, upper=1/(1-dd) if dd < 1 else None,
                              universal_lower=1-udn,
                              universal_upper=1/(1-udd) if udd < 1 else None,
                              mode_ratios_to_rotor=ratios,
                              maximum_relative_error=max(abs(r-1) for r in ratios)))
    for x in (.05, .2, 1.5):
        selected = [c for c in cases if c['x'] == x]
        errors = [c['maximum_relative_error'] for c in selected]
        assert all(a > b for a, b in zip(errors[:-1], errors[1:]))
        assert errors[-1] < .0002
    return cases


def continuum_field(ell, q, x, cutoff, epsilon, grid):
    # beta=1. Spatial momentum cutoff/smearing, not shifted-frequency filtering.
    count = int(math.ceil(cutoff*ell/(2*math.pi)))
    u = 2*math.pi*(np.arange(count)+.5)/ell
    z = np.exp(-u)
    alpha = 2*math.pi*np.arange(grid)/grid-math.pi
    v = z/(1+z)**2
    factors = 1-4*v[:, None]*np.sin(alpha[None, :]/2)**2
    logphi = q*np.log(factors).sum(axis=0)
    phi = np.exp(logphi)
    charges = np.arange(-int(math.ceil(math.sqrt(70/x)))-2,
                        int(math.ceil(math.sqrt(70/x)))+3)
    fourier = np.exp(-1j*charges[:, None]*alpha[None, :])
    kernels = [np.exp(-x*(charges+s)**2)@fourier for s in (0, .5, 1)]
    f = np.array([(phi[None, :]/factors@kernel/grid).real for kernel in kernels])
    D = (1+z*z)*f[0]+2*z*f[2]
    ratio = (1+z)**2*math.exp(-x/4)*f[1]/D
    assert np.min(ratio) >= math.exp(-x/2)-2.e-11
    assert np.max(ratio) <= math.exp(x/4)/math.cosh(x/2)+2.e-11
    b0_modes = np.exp(-u/2)/(1+np.exp(-u))
    b = float(2/ell*np.sum(np.exp(-epsilon*u)*b0_modes*ratio))
    free = float(2/ell*np.sum(np.exp(-epsilon*u)*b0_modes))
    return dict(ell=ell, species=q, x=x, cutoff=cutoff, epsilon=epsilon,
                angular_grid=grid, field_midpoint=b, free_midpoint=free,
                field_ratio=b/free, lowest_mode_ratio=float(ratio[0]),
                reservoir_V=float(q*sum(v)), modes_per_charge=count)


def continuum_checks():
    cases = []
    for x in (.05, .4):
        target = theta_data(x)['C']/2
        for ell in (2., 4., 8., 24., 72., 216.):
            row = continuum_field(ell, 2, x, 32., 0., 2048)
            refined = continuum_field(ell, 2, x, 40., 0., 4096)
            row['cutoff_and_grid_refinement_error'] = abs(row['field_midpoint']-refined['field_midpoint'])
            row['infinite_reservoir_point_value'] = target
            row['absolute_limit_error'] = abs(row['field_midpoint']-target)
            row['refined_absolute_limit_error'] = abs(refined['field_midpoint']-target)
            assert row['cutoff_and_grid_refinement_error'] < 8.e-8
            cases.append(row)
        errors = [c['refined_absolute_limit_error'] for c in cases if c['x'] == x]
        # At large ell the UV tail, not volume error, sets the observed floor.
        assert errors[-1] < 5.e-8 and errors[-1] < errors[0]/1000
    smeared = []
    for epsilon in (.1, .03, .01, 0.):
        row = continuum_field(216., 2, .4, 40., epsilon, 4096)
        assert row['field_ratio'] >= math.exp(-.2)-2.e-12
        smeared.append(row)
    assert all(a['field_midpoint'] < b['field_midpoint']
               for a, b in zip(smeared[:-1], smeared[1:]))
    return cases, smeared


def lattice_checks():
    cases = []
    for ell in (2., 4., 8., 24.):
        c = math.pi/ell
        cutoff = int(math.ceil(math.sqrt(70/c)))+1
        integers = np.arange(-cutoff, cutoff+1)
        # Separate product-versus-charge check of the Fock-space branching.
        levels = 2*c*(np.arange(int(math.ceil(80/(2*c))))+.5)
        z = np.exp(-levels)
        v = z/(1+z)**2
        theta0 = sum(np.exp(-c*integers**2))
        jacobi_errors = []
        for alpha in (.17, 1.1, 2.8):
            product = np.prod(1-4*v*np.sin(alpha/2)**2)
            charge_sum = np.sum(np.exp(-c*integers**2)*np.cos(integers*alpha))/theta0
            jacobi_errors.append(abs(product-charge_sum))
        assert max(jacobi_errors) < 2.e-13
        free_point = float(2/ell*np.sum(np.exp(-levels/2)/(1+np.exp(-levels))))
        free_charge = theta_data(c)['C']
        for q in (1, 2, 3):
            grid = np.array(np.meshgrid(*([integers]*q), indexing='ij')).reshape(q, -1)
            shifted = grid.astype(float)
            shifted[0] += .5
            for x in (.05, .4):
                Z = np.exp(-c*(grid*grid).sum(axis=0)-x*grid.sum(axis=0)**2).sum()
                th = np.exp(-c*(shifted*shifted).sum(axis=0)-x*shifted.sum(axis=0)**2).sum()
                charged = math.exp(-(c+x)/4)*th/Z
                predicted = free_point*charged/free_charge
                direct = continuum_field(ell, q, x, 40., 0., 4096)
                error = abs(predicted-direct['field_midpoint'])
                assert error < 2.e-9
                cases.append(dict(ell=ell, species=q, x=x, lattice_cutoff=cutoff,
                                  free_point=free_point, lattice_point=predicted,
                                  transition_point=direct['field_midpoint'],
                                  midpoint_error=error,
                                  jacobi_product_error=max(jacobi_errors)))
    return cases


def main():
    small, biased = small_fock_checks()
    fourier = fourier_checks()
    limits = reservoir_limit_checks()
    continuum, smeared = continuum_checks()
    lattice = lattice_checks()
    theta = [dict(x=x, **theta_data(x)) for x in (.00004, .03, .4, 2., 6.)]
    candidate = dict(k_min=2000, electric_loss_upper=1/(24*2000),
                     condition='The coefficient of the total-charge Hamiltonian must be matched; no angular or stress bound.')
    artifacts = ['notes/electric-reservoir-midpoint.md', 'checks/electric_reservoir_midpoint.py']
    receipt = dict(status='PASS', story='S072', date='2026-10-13',
        scope='Global electric charging plus the full free radial fermion reservoir',
        numpy=np.__version__, scipy=scipy.__version__, concentration_constant=CSTAR,
        fock_cases=small, biased_seed_counterexample=biased, fourier_cases=fourier,
        reservoir_limit_cases=limits, continuum_cases=continuum,
        finite_circle_lattice_cases=lattice,
        spatial_regulator_cases=smeared, theta_cases=theta, candidate=candidate,
        artifacts=[dict(path=p, sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()) for p in artifacts],
        limitations=['The actual radial Coulomb kernel is not merely a total-charge square.',
                     'A screened inertia and the same explicit fermions cannot be combined without matching.',
                     'The angular sector, evolving state, stress and geometry are not computed.',
                     'The infinite-reservoir compact factor replaces the Gaussian zero mode; it is not an extra independent rotor.'])
    output = ROOT/'evidence/electric-reservoir-midpoint.json'
    output.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status='PASS', fock_cases=len(small),
        maximum_trace_error=max(c['maximum_trace_error'] for c in small),
        maximum_fourier_error=max(c['maximum_fourier_error'] for c in fourier),
        reservoir_cases=len(limits), continuum_cases=len(continuum),
        lattice_cases=len(lattice), maximum_lattice_error=max(c['midpoint_error'] for c in lattice),
        maximum_continuum_refinement=max(c['cutoff_and_grid_refinement_error'] for c in continuum),
        candidate_loss_upper=candidate['electric_loss_upper']), indent=2))


if __name__ == '__main__':
    main()
