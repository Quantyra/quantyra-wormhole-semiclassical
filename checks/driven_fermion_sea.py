"""Finite CAR, unequal-return evolution and continuum form-bound checks."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigh, expm

ROOT = Path(__file__).resolve().parents[1]
T = np.pi
KAPPA = 1+np.arcsinh(np.sqrt(3))/(2*np.sqrt(3))
MT = max(2/np.pi**2, 2/np.pi*(1/(2*np.pi)+1/T))
AT = 4*np.pi*MT*(1+T*T/4)
PEAK = 2*(5/6)**2.5/np.sqrt(6)


def restriction_constant(lengths):
    return 1+T/(2*np.pi)*(2+np.log1p(max(lengths)/np.pi))


def fourier(k, which):
    def window(x):
        return T*np.sinc(x*T/(2*np.pi))
    def f1(x):
        return .5*window(x)+.25*(window(x+2*np.pi/T)+window(x-2*np.pi/T))
    return f1(k) if which == 0 else (f1(k+2*np.pi/T)-f1(k-2*np.pi/T))/(2j)


def matrices(q, N, lengths):
    m = np.arange(q)-(q-1)/2
    jp = np.diag(np.sqrt(np.arange(1, q)*(q-np.arange(1, q))), -1)
    axes = [np.diag(m), (jp+jp.T)/2]
    channel = np.repeat(np.arange(q), 2*N)
    freq = np.concatenate([2*np.pi*(np.arange(-N, N)+.5)/L for L in lengths])
    gap = freq[None, :]-freq[:, None]
    norm = np.sqrt(np.array(lengths)[channel[:, None]]*np.array(lengths)[channel[None, :]])
    potentials = [X[channel[:, None], channel[None, :]]*fourier(gap, a)/norm
                  for a, X in enumerate(axes)]
    assert max(np.linalg.norm(A-A.conj().T) for A in potentials) < 1e-11
    # Reversing chirality conjugates these matrices. A fixed channel phase
    # implements the conjugation, so both chirality energies are equal here.
    phase = (-1.)**channel
    assert max(np.linalg.norm(A.conj()-phase[:, None]*A*phase[None, :])
               for A in potentials) < 1e-11
    return freq, potentials


def pair_form(freq, A):
    p, h = np.flatnonzero(freq > 0), np.flatnonzero(freq < 0)
    return float(np.sum(abs(A[np.ix_(p, h)])**2/(freq[p, None]-freq[None, h])))


def diagonal_form(freq, A):
    values = []
    for mask in [freq > 0, freq < 0]:
        z = np.flatnonzero(mask)
        scaled = A[np.ix_(z, z)]/np.sqrt(abs(freq[z, None]*freq[None, z]))
        values.append(float(max(abs(np.linalg.eigvalsh(scaled)))))
    return max(values)


def projection_checks():
    rows = []
    for q in [2, 3, 5]:
        lengths = T+np.linspace(.1, 7, q)
        for N in [4, 12, 32]:
            freq, potentials = matrices(q, N, lengths)
            tq = q*(q*q-1)/12
            for a, A in enumerate(potentials):
                pair = pair_form(freq, A)
                upper_pair = AT*tq*T*([3/8, 5/32][a])
                diag = diagonal_form(freq, A)
                upper_diag = restriction_constant(lengths)*(q-1)/2
                assert pair <= upper_pair and diag <= upper_diag
                rows.append(dict(q=q, modes_per_sign=N, profile=a, pair_form=pair,
                                 upper_pair=upper_pair, diagonal_form=diag, upper_diagonal=upper_diag))
    # Spatial quadrature independent of the analytic Fourier matrix formula.
    errors = []
    for k in [-9.1, -2., 0., .7, 4.8]:
        for a in [0, 1]:
            def f(s):
                base = np.cos(np.pi*s/T)**2
                return base if a == 0 else base*np.sin(2*np.pi*s/T)
            direct = quad(lambda s: f(s)*np.cos(k*s), -T/2, T/2, epsabs=1e-12)[0]
            direct += 1j*quad(lambda s: f(s)*np.sin(k*s), -T/2, T/2, epsabs=1e-12)[0]
            errors.append(abs(direct-fourier(k, a)))
    assert max(errors) < 2e-12
    return dict(cases=rows, maximum_fourier_quadrature_error=float(max(errors)))


def covariance_from_fock(state, size):
    # gamma_ij=<c_j^dagger c_i>, with no quasifree assumption on state.
    vectors = []
    for a in range(size):
        out = np.zeros_like(state)
        for bits in range(len(state)):
            if (bits >> a) & 1:
                sign = (-1)**((bits & ((1 << a)-1)).bit_count())
                out[bits ^ (1 << a)] += sign*state[bits]
        vectors.append(out)
    V = np.array(vectors)
    return V@V.conj().T


def algebra_checks():
    rng = np.random.default_rng(531110)
    rows = []
    size = 8
    freq = np.array([-.13, -.7, -2.1, -5., .07, .6, 1.7, 6.])
    negative = freq < 0
    P = np.diag(negative.astype(float))
    for case in range(36):
        if case < 12:
            state = rng.normal(size=2**size)+1j*rng.normal(size=2**size)
            state /= np.linalg.norm(state)
            gamma = covariance_from_fock(state, size)
            kind = 'general_Fock_vector'
        else:
            W = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
            W = expm(.12j*(W+W.conj().T))
            occupancy = negative.astype(float) if case < 24 else rng.uniform(0, 1, size)
            gamma = W@np.diag(occupancy)@W.conj().T
            kind = 'Slater' if case < 24 else 'mixed_covariance'
        A = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
        A = (A+A.conj().T)/2
        energy = float(np.trace(np.diag(freq)@(gamma-P)).real)
        p, h = np.flatnonzero(freq > 0), np.flatnonzero(freq < 0)
        weighted_cross = np.sum((freq[p, None]-freq[None, h])*abs(gamma[np.ix_(p, h)])**2)
        a, b = pair_form(freq, A), diagonal_form(freq, A)
        actual = abs(np.trace(A@(gamma-P)))
        bound = b*energy+2*np.sqrt(a*energy)
        assert weighted_cross <= energy+1e-11 and actual <= bound+1e-11
        rows.append(dict(kind=kind, energy=energy, weighted_cross=float(weighted_cross),
                         expectation=float(actual), bound=float(bound)))
    return rows


def pulse(u):
    h1 = np.sin(np.pi*u)**4
    return h1, h1*np.sin(2*np.pi*u)


def bound_data(q, lengths, epsilon, finite_a=None):
    j, tq = (q-1)/2, q*(q*q-1)/12
    bstar = restriction_constant(lengths)*epsilon*j*(1+PEAK)
    B = restriction_constant(lengths)*epsilon*j*(2+4*PEAK)
    # Both chiralities; finite_a optionally gives independent discrete pair forms.
    a1, a2 = ([2*AT*tq*T*3/8, 2*AT*tq*T*5/32] if finite_a is None
              else [2*x for x in finite_a])
    root_a = epsilon*(np.sqrt(a1)+PEAK*np.sqrt(a2))
    W = epsilon*(2*np.sqrt(a1)+4*PEAK*np.sqrt(a2))
    upper = 2*np.exp(2*B)*(2*root_a+np.sqrt(2)*W)**2
    return dict(b_star=float(bstar), integrated_b=float(B), a_star=float(root_a**2),
                integrated_root_a_dot=float(W), energy_bound=float(upper))


def evolution(q, N, lengths, epsilon, duration, steps):
    freq, As = matrices(q, N, lengths)
    eig = [eigh(A, driver='evd') for A in As]
    neg = np.flatnonzero(freq < 0)
    pos = np.flatnonzero(freq > 0)
    Y = np.eye(len(freq), dtype=complex)
    dt = duration/steps
    free = np.exp(-.5j*dt*freq)
    def kick(state, which, strength):
        val, vec = eig[which]
        return vec@(np.exp(-1j*strength*val)[:, None]*(vec.conj().T@state))
    maximum_energy = 0.
    work, old_power = 0., 0.
    for step in range(steps):
        h1, h2 = pulse((step+.5)/steps)
        Y = free[:, None]*Y
        Y = kick(Y, 0, .5*dt*epsilon*h1)
        Y = kick(Y, 1, dt*epsilon*h2)
        Y = kick(Y, 0, .5*dt*epsilon*h1)
        Y = free[:, None]*Y
        # Use transition blocks, avoiding subtraction of two filled-sea energies.
        ph, hp = Y[np.ix_(pos, neg)], Y[np.ix_(neg, pos)]
        pp, hh = ph@ph.conj().T, -hp@hp.conj().T
        cross = ph@Y[np.ix_(neg, neg)].conj().T
        energy = 2*float(np.dot(freq[pos], pp.diagonal()).real+np.dot(freq[neg], hh.diagonal()).real)
        maximum_energy = max(maximum_energy, energy)
        u = (step+1)/steps
        h1, h2 = pulse(u)
        dh1 = 4*np.pi*np.sin(np.pi*u)**3*np.cos(np.pi*u)
        dh2 = dh1*np.sin(2*np.pi*u)+2*np.pi*h1*np.cos(2*np.pi*u)
        derivative = epsilon/duration*(dh1*As[0]+dh2*As[1])
        power = 2*float((np.trace(derivative[np.ix_(pos, pos)]@pp)
                        +np.trace(derivative[np.ix_(neg, neg)]@hh)
                        +2*np.trace(derivative[np.ix_(neg, pos)]@cross).real).real)
        work += dt*(old_power+power)/2
        old_power = power
    gram_error = float(np.linalg.norm(Y.conj().T@Y-np.eye(len(freq))))
    bounds = bound_data(q, lengths, epsilon)
    tighter = bound_data(q, lengths, epsilon, [pair_form(freq, A) for A in As])
    assert bounds['b_star'] < .25
    assert gram_error < 2e-9 and maximum_energy < min(bounds['energy_bound'], tighter['energy_bound'])
    return dict(final_energy=energy, maximum_energy=maximum_energy, gram_error=gram_error,
                integrated_work=work, work_error=abs(work-energy),
                continuum_bound=bounds, finite_pair_bound=tighter)


def dynamics_checks():
    rows = []
    for q, duration in [(3, .7), (3, 3.), (4, 1.5)]:
        lengths = T+np.linspace(.2, 5., q)
        for N in [6, 12, 18]:
            values = [evolution(q, N, lengths, .018, duration, steps) for steps in [32, 64, 128]]
            energies = [v['final_energy'] for v in values]
            difference = [abs(energies[i+1]-energies[i]) for i in range(2)]
            assert difference[1] < .4*difference[0]+2e-13, (q, N, duration, energies, difference)
            assert values[-1]['work_error'] < 2e-8*max(1, values[-1]['final_energy'])
            rows.append(dict(q=q, modes_per_sign=N, duration=duration,
                             refinements=values, time_differences=difference))
    # The largest spectral truncations should have resolved these smooth pulses.
    for start in [0, 3, 6]:
        values = [row['refinements'][-1]['final_energy'] for row in rows[start:start+3]]
        assert abs(values[-1]-values[-2]) < .015*values[-1]+1e-9
    return rows


def independent_evolution():
    rows = []
    for q, duration in [(3, 3.), (4, 1.5)]:
        N, epsilon = 6, .018
        lengths = T+np.linspace(.2, 5., q)
        freq, As = matrices(q, N, lengths)
        size = len(freq)
        def rhs(t, flat):
            U = flat.reshape(size, size)
            h1, h2 = pulse(t/duration)
            return (-1j*(freq[:, None]*U+epsilon*(h1*As[0]+h2*As[1])@U)).ravel()
        sol = solve_ivp(rhs, (0, duration), np.eye(size, dtype=complex).ravel(),
                        method='DOP853', rtol=2e-11, atol=2e-13, t_eval=[duration])
        assert sol.success
        U = sol.y[:, -1].reshape(size, size)
        p, h = np.flatnonzero(freq > 0), np.flatnonzero(freq < 0)
        energy = 2*(np.dot(freq[p], np.sum(abs(U[np.ix_(p, h)])**2, axis=1))
                    -np.dot(freq[h], np.sum(abs(U[np.ix_(h, p)])**2, axis=1)))
        compare = evolution(q, N, lengths, epsilon, duration, 128)
        error = abs(energy-compare['final_energy'])
        assert error < 2e-9
        rows.append(dict(q=q, duration=duration, ode_energy=float(energy),
                         split_energy=compare['final_energy'], absolute_error=float(error), nfev=sol.nfev))
    return rows


def scales():
    rows = []
    for exponent in [4, 8, 12, 16, 24, 40, 60]:
        q = 10.**exponent
        delta = np.log(q)**-2
        epsilon = np.log(q)**3/(2*np.sqrt(np.pi)*q**1.5)
        lengths = [T, T+2*KAPPA*delta*q]
        b = bound_data(q, lengths, epsilon)
        valid = b['b_star'] <= .25
        rows.append(dict(q=q, delta=float(delta), omega_ell=float(epsilon), **b,
                         smallness_condition=valid,
                         bound_over_q_delta=float(b['energy_bound']/(q*delta)) if valid else None))
    assert rows[-1]['bound_over_q_delta'] < 1e-30
    return rows


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    result = dict(status='passed_driven_free_sea_energy_not_interacting_stress',
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        algebra=algebra_checks(), projections=projection_checks(), dynamics=dynamics_checks(),
        independent_dynamics=independent_evolution(), scales=scales(),
        script_sha256=sha(Path(__file__)),
        input_sha256={p: sha(ROOT/p) for p in ['notes/driven-fermion-sea.md',
            'notes/propagating-photon-polarization.md', 'notes/moving-reflection-load.md',
            'notes/angular-current-coset.md']},
        continuum_state_energy_proved_in_defined_model=True,
        full_quantum_gauge_state=False,local_renormalized_stress_bound=False,
        physical_return_matching=False,novel_physical_endpoint=False,lean_invocations=0,goal_complete=False)
    (ROOT/'evidence/driven-fermion-sea.json').write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], algebra_cases=len(result['algebra']),
        projection_cases=len(result['projections']['cases']), dynamics_cases=len(result['dynamics']),
        maximum_gram_error=max(v['gram_error'] for r in result['dynamics'] for v in r['refinements']))))


if __name__ == '__main__':
    main()
