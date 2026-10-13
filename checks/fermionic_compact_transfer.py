"""S072 exact fermionic compact transfer, branching and driven joint state.

Finite CAR checks do not stand in for continuum stress or a physical
black-hole state. The dimension-uniform bound is proved in the note.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter, defaultdict
import hashlib
import json

import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import block_diag, eigh
from scipy.sparse import coo_matrix, diags, csr_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]


def spin(q):
    j = (q-1)/2
    m = np.arange(q)-j
    plus = np.diag(np.sqrt(j*(j+1)-m[:-1]*(m[:-1]+1)), -1)
    return [(plus+plus.T)/2, (plus-plus.T)/(2j), np.diag(m)]


def one_particle(q, copies=1):
    tj = spin(q)
    generators = [block_diag(*([t]*copies+[-t.T]*copies)) for t in tj]
    charges = np.r_[np.ones(q*copies), -np.ones(q*copies)]
    return charges, generators


def one_body(matrix, basis, lookup):
    rows, cols, values = [], [], []
    nonzero = list(zip(*np.nonzero(matrix)))
    for col, state in enumerate(basis):
        for a, b in nonzero:
            a, b = int(a), int(b)
            if not (state >> b) & 1:
                continue
            sign = (-1)**((state & ((1 << b)-1)).bit_count())
            middle = state ^ (1 << b)
            if (middle >> a) & 1:
                continue
            sign *= (-1)**((middle & ((1 << a)-1)).bit_count())
            target = middle | (1 << a)
            rows.append(lookup[target]); cols.append(col)
            values.append(sign*matrix[a, b])
    return coo_matrix((values, (rows, cols)),
                      shape=(len(basis), len(basis)), dtype=complex).tocsr()


def operators(q, copies=1):
    charges, ts = one_particle(q, copies)
    D = len(charges)
    basis = list(range(1 << D)); lookup = {b: i for i, b in enumerate(basis)}
    Q = one_body(np.diag(charges), basis, lookup)
    Js = [one_body(t, basis, lookup) for t in ts]
    J2 = sum((J@J for J in Js), csr_matrix(Q.shape, dtype=complex))
    return D, basis, Q, Js, J2


def majoranas(D):
    states = np.arange(1 << D)
    for a in range(D):
        permutation = states ^ (1 << a)
        signs = np.array([(-1)**((int(s) & ((1 << a)-1)).bit_count())
                          for s in states], dtype=complex)
        yield permutation, signs
        yield permutation, 1j*signs*(1-2*((states >> a) & 1))


def channel(A, D):
    answer = np.zeros_like(A, dtype=complex)
    for permutation, phase in majoranas(D):
        answer += phase.conj()[:, None]*A[np.ix_(permutation, permutation)]*phase[None, :]
    return answer/(2*D)


def branch_check(q, copies, Q, Js, J2):
    charges = np.rint(Q.diagonal().real).astype(int)
    weights = np.rint(2*Js[2].diagonal().real).astype(int)
    counts = Counter(zip(charges, weights))
    branches = []
    total = 0
    J2 = J2.toarray()
    for charge in sorted(set(charges)):
        indices = np.flatnonzero(charges == charge)
        eigenvalues = np.linalg.eigvalsh(J2[np.ix_(indices, indices)])
        spins = np.rint(np.sqrt(1+4*np.maximum(eigenvalues, 0))-1).astype(int)
        spectral = Counter(spins)
        for twice_spin in range(max(abs(weights[indices]))+1):
            multiplicity = counts[(charge, twice_spin)]-counts[(charge, twice_spin+2)]
            assert multiplicity >= 0
            if multiplicity:
                dimension = multiplicity*(twice_spin+1)
                assert spectral[twice_spin] == dimension
                assert (twice_spin-(q-1)*charge) % 2 == 0
                total += dimension
                branches.append(dict(charge=int(charge), twice_spin=int(twice_spin),
                                     multiplicity=int(multiplicity)))
    assert total == len(charges)
    # This center-allowed branch is absent even though no rotor forbids it.
    if q == 4 and copies == 1:
        assert counts[(4, 0)] == 1 and counts[(4, 2)] == 0
    return branches


def thermal_case(q, copies, x, a):
    D, basis, Q, Js, J2 = operators(q, copies)
    K = (x*(Q@Q)+a*J2).toarray()
    Cj = (q*q-1)/4
    lam = x+a*Cj
    E_K = channel(K, D)
    identity_error = float(np.max(abs(E_K-K-lam*np.eye(len(K))+4*K/D)))
    assert identity_error < 2.e-12
    energies, vectors = eigh(K)
    weights = np.exp(-energies); probabilities = weights/sum(weights)
    W = np.zeros_like(K.real)
    for permutation, phase in majoranas(D):
        transformed = vectors.conj().T@(phase.conj()[:, None]*vectors[permutation])
        W += abs(transformed)**2/(2*D)
    assert max(np.max(abs(W.sum(axis=0)-1)), np.max(abs(W-W.T))) < 2.e-12
    midpoint = np.sqrt(probabilities[:, None]*probabilities[None, :])*W
    delta = energies[:, None]-energies[None, :]
    S = float(np.sum((probabilities[None, :]-probabilities[:, None])*delta*W))
    mean = float(probabilities@energies)
    rhs = 2*lam-8*mean/D
    assert abs(S-rhs) < 2.e-12
    C = float(np.sum(midpoint))
    root_rho = (vectors*np.sqrt(probabilities))@vectors.conj().T
    direct = float(np.trace(root_rho@channel(root_rho, D)).real)
    assert abs(C-direct) < 2.e-12
    M2 = float(np.sum(midpoint*delta**2))
    assert 0 <= 1-C+2.e-12 <= S/8+4.e-12
    assert M2 <= S+2.e-12 and S <= 2*lam+2.e-12
    times = []
    for tau in (0., .2, 1., 3.):
        value = float(np.sum(midpoint*np.cos(delta*tau)))
        bound = 1-S/8-S*tau*tau/2
        assert value >= bound-2.e-12
        times.append(dict(tau=tau, correlation=value, lower_bound=bound))
    return dict(q=q, copies=copies, D=D, x=x, a=a, lambda_value=lam,
                identity_error=identity_error, mean_K=mean, S=S,
                affinity=C, affinity_upper_loss=S/8, second_moment=M2,
                direct_affinity_error=abs(C-direct), times=times,
                branches=branch_check(q, copies, Q, Js, J2)), C


def spin_environment():
    J = [np.kron(t, np.eye(3)) for t in spin(2)]
    B = [np.kron(np.eye(2), t) for t in spin(3)]
    total = [a+b for a, b in zip(J, B)]
    J2 = sum(t@t for t in total)
    vals, vecs = eigh(J2)
    low = vecs[:, abs(vals-.75) < 1.e-12]
    mixed = sum(t@j for t, j in zip(total, J))
    projected = low.conj().T@mixed@low
    error = float(np.max(abs(projected+.25*np.eye(2))))
    assert error < 1.e-12
    return dict(fermion_spin=.5, environment_spin=1, total_spin=.5,
                mixed_product=-.25, matrix_error=error)


def drive_case(q, x, a, final_amplitude=1., max_step=.025, tolerance=2.e-12):
    charges, ts = one_particle(q)
    D = len(charges)
    basis = [sum(1 << i for i in c) for c in combinations(range(2*D), D)]
    lookup = {b: i for i, b in enumerate(basis)}
    zeros = np.zeros((D, D))
    QL = one_body(block_diag(np.diag(charges), zeros), basis, lookup)
    QR = one_body(block_diag(zeros, np.diag(charges)), basis, lookup)
    JL = [one_body(block_diag(t, zeros), basis, lookup) for t in ts]
    JR = [one_body(block_diag(zeros, t), basis, lookup) for t in ts]
    KL = x*(QL@QL)+a*sum((J@J for J in JL), csr_matrix(QL.shape))
    KR = x*(QR@QR)+a*sum((J@J for J in JR), csr_matrix(QR.shape))
    K = KL+KR
    hopping = np.block([[zeros, 1j*np.eye(D)], [-1j*np.eye(D), zeros]])
    X = one_body(hopping, basis, lookup)
    amplitudes = {0: 1.+0j}
    for mode in range(D):
        updated = defaultdict(complex)
        for state, value in amplitudes.items():
            for target_mode, coefficient in ((mode, 1), (mode+D, -1j)):
                sign = (-1)**((state & ((1 << target_mode)-1)).bit_count())
                updated[state | (1 << target_mode)] += coefficient*sign*value/np.sqrt(2)
        amplitudes = updated
    initial = np.array([amplitudes.get(b, 0j) for b in basis])
    assert np.linalg.norm(X@initial-D*initial) < 2.e-12
    initial = expm_multiply(-KL/2, initial)
    initial /= np.linalg.norm(initial)
    _, expected_C = thermal_case(q, 1, x, a)
    actual_C = float(np.vdot(initial, X@initial).real/D)
    assert abs(actual_C-expected_C) < 2.e-12
    initial_energy = float(np.vdot(initial, K@initial).real)
    tau = .25
    end = np.pi+tau

    def S(t):
        z = np.clip(t/tau, 0., 1.)
        return 10*z**3-15*z**4+6*z**5

    def primitive(t):
        if t <= 0:
            return 0.
        if t >= tau:
            return t-tau/2
        z = t/tau
        return tau*(2.5*z**4-3*z**5+z**6)

    def amplitude(t):
        return final_amplitude*(primitive(t)-primitive(t-np.pi))/np.pi

    def rate(t):
        return final_amplitude*(S(t)-S(t-np.pi))/np.pi

    def rhs(t, y):
        psi = y[:-1]
        Xpsi = X@psi
        return np.r_[-1j*(K@psi-amplitude(t)*Xpsi),
                     -rate(t)*np.vdot(psi, Xpsi).real]

    samples = np.linspace(0, end, 81)
    times, states = [], []
    current = np.r_[initial, 0j]
    for left, right in zip((0, tau, np.pi), (tau, np.pi, end)):
        selected = np.r_[samples[(samples >= left) & (samples < right)], right]
        solution = solve_ivp(rhs, (left, right), current, method='DOP853',
                             rtol=tolerance, atol=tolerance/100,
                             max_step=max_step, t_eval=selected)
        assert solution.success
        current = solution.y[:, -1]
        times.extend(solution.t); states.extend(solution.y.T)
    energy_errors, constraints, affinities, norms = [], [], [], []
    for t, column in zip(times, states):
        psi, work = column[:-1], column[-1].real
        energy = float(np.vdot(psi, K@psi-amplitude(t)*(X@psi)).real)
        energy_errors.append(abs(energy-initial_energy-work))
        constraint = np.linalg.norm((QL+QR)@psi)**2
        constraint += sum(np.linalg.norm((l+r)@psi)**2 for l, r in zip(JL, JR))
        constraints.append(float(constraint))
        norms.append(abs(np.vdot(psi, psi).real-1))
        affinities.append(float(np.vdot(psi, X@psi).real/D))
    assert max(energy_errors) < 3.e-9, (q, x, a, max(energy_errors), max(norms), energy)
    assert max(constraints) < 1.e-20
    assert max(norms) < 2.e-10
    return dict(q=q, D=D, dimension=len(basis), x=x, a=a,
                final_amplitude=final_amplitude, end=end,
                initial_affinity=actual_C, minimum_affinity=min(affinities),
                final_affinity=affinities[-1], affinity_change=affinities[-1]-actual_C,
                initial_energy=initial_energy, final_energy=energy,
                maximum_work_error=max(energy_errors),
                maximum_gauss_constraint=max(constraints),
                maximum_norm_error=max(norms), max_step=max_step,
                tolerance=tolerance)


def candidate():
    mp.mp.dps = 70
    rows = []
    for integer in (2000, 4000, 10000):
        k = mp.mpf(integer); q = k**4
        X = (mp.sqrt(1+12*mp.pi*k**6)-1)/6
        for z in (mp.mpf(1)/4, mp.mpf(7)/25):
            electric = X/(8*k*z*(1+6*X))
            angular = mp.pi*(q*q-1)/(8*q/k*z*(1+6*X)*(X+mp.mpf(5)/6))
            upper = 1/(12*k*z)
            assert electric+angular < upper <= 1/(3*k)
            rows.append(dict(k=integer, z0=float(z), lambda_upper=float(electric+angular),
                             initial_loss_upper=float((electric+angular)/4),
                             uniform_loss_upper=float(1/(12*k))))
    return rows


def main():
    thermal = [thermal_case(q, copies, x, a)[0] for q, copies, x, a in
               ((2, 1, .08, .012), (3, 1, .2, .05),
                (4, 1, .08, .012), (2, 2, .2, .05))]
    driven = [drive_case(q, x, a) for q, x, a in
              ((2, .08, .012), (2, .2, .05), (3, .2, .05))]
    refined = drive_case(2, .2, .05, max_step=.0125, tolerance=2.e-13)
    refinement = max(abs(refined[key]-driven[1][key]) for key in
                     ('final_energy', 'final_affinity', 'minimum_affinity'))
    assert refinement < 2.e-9
    result = dict(date='2026-10-13', story='S072', status='PASS',
                  scope='Exact CAR compact-charge thermal theorem and finite driven joint-state diagnostics',
                  versions=dict(numpy=np.__version__, scipy=scipy.__version__),
                  thermal=thermal, mixed_spin_counterexample=spin_environment(),
                  driven=driven, refined=refined, refinement_difference=refinement,
                  candidate=candidate(), limitations=[
                      'No identification of the physical radial or black-hole Hamiltonian with K',
                      'Uniform canonical-mode channel is not an unsmeared boundary transfer operator',
                      'Driven finite charge model has no continuum stress or JT geometry',
                      'Physical mixed charges and residual Hamiltonian terms remain unbounded'])
    result['artifacts'] = [dict(path=p, sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest())
                           for p in ('notes/fermionic-compact-transfer.md',
                                     'checks/fermionic_compact_transfer.py')]
    output = ROOT/'evidence/fermionic-compact-transfer.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status='PASS', thermal_cases=len(thermal), driven_cases=len(driven),
                         maximum_operator_error=max(r['identity_error'] for r in thermal),
                         maximum_work_error=max(r['maximum_work_error'] for r in driven),
                         maximum_gauss_error=max(r['maximum_gauss_constraint'] for r in driven),
                         refinement_difference=refinement,
                         k2000_conditional_loss=result['candidate'][0]['initial_loss_upper'],
                         output=str(output)), indent=2))


if __name__ == '__main__':
    main()
