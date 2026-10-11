"""Finite independent checks of the angular current route and its limits."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import sph_harm_y
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT/'notes/angular-current-renormalization.md'
OUT = ROOT/'evidence/angular-current-renormalization.json'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def exact_checks():
    rt, I = sy.sqrt, sy.I
    Jx = sy.Matrix([[0, 1, 0], [1, 0, 1], [0, 1, 0]])/rt(2)
    Jy = sy.Matrix([[0, -I, 0], [I, 0, -I], [0, I, 0]])/rt(2)
    Jz = sy.diag(1, 0, -1)
    E = [Jx/rt(2), Jy/rt(2), Jz/rt(2),
         (Jx*Jx-Jy*Jy)/rt(2), (Jx*Jy+Jy*Jx)/rt(2),
         (Jx*Jz+Jz*Jx)/rt(2), (Jy*Jz+Jz*Jy)/rt(2), (3*Jz*Jz-2*sy.eye(3))/rt(6)]
    assert sy.Matrix([[sy.simplify(sy.trace(a*b)) for b in E] for a in E]) == sy.eye(8)
    a, b = sy.symbols('a b', real=True)
    G = [a]*3+[b]*5
    beta = sy.zeros(8)
    for i in range(8):
        for j in range(8):
            comm = -I*(E[i]*E[j]-E[j]*E[i])
            f = sy.Matrix([sy.simplify(sy.trace(e*comm)) for e in E])
            beta += G[i]*G[j]*(f*f.T)/2
    beta = sy.simplify(beta)
    expected = sy.diag(*([(a*a+5*b*b)/2]*3+[3*a*b]*5))
    assert sy.simplify(beta-expected) == sy.zeros(8)
    assert sy.simplify(expected[0, 0].subs(b, a/15)/a) == sy.Rational(23, 45)*a
    n = sy.symbols('n', integer=True, positive=True)
    assert sy.simplify((2*n+1)/(n*n*(n+1)**2)-1/n**2+1/(n+1)**2) == 0
    return [dict(name='q3_orthonormal_basis', exact=True),
            dict(name='q3_complete_one_loop_tensor', spin_one=str(expected[0, 0]), spin_two=str(expected[3, 3])),
            dict(name='q3_tree_ray_failure', spin_one_log_rate='23*a/45', spin_two_log_rate='3*a'),
            dict(name='frobenius_sum_telescopes', exact=True)]


def angular_basis(q, refinement=1):
    ntheta, nphi = refinement*(2*q+8), refinement*(4*q+8)
    u, wu = np.polynomial.legendre.leggauss(ntheta)
    phi = 2*np.pi*np.arange(nphi)/nphi
    th, ph = np.meshgrid(np.arccos(u), phi, indexing='ij')
    th, ph = th.ravel(), ph.ravel()
    w = np.repeat(wu, nphi)*2*np.pi/nphi
    j = (q-1)/2
    eta = np.array([math.sqrt(q*math.comb(q-1, k)/(4*np.pi))*np.sin(th/2)**k
                    *np.cos(th/2)**(q-1-k)*np.exp(1j*(j-k)*ph) for k in range(q)]).T
    E, degree, sigma, beta_l = [], [], [], []
    b = 1.
    for l in range(1, q):
        b *= (q-l)/(q+l)
        beta_l.append(b)
        harmonics = [sph_harm_y(l, 0, th, ph).real]
        for m in range(1, l+1):
            h = sph_harm_y(l, m, th, ph)
            harmonics.extend([np.sqrt(2)*h.real, np.sqrt(2)*h.imag])
        for h in harmonics:
            T = eta.conj().T@((w*h)[:, None]*eta)
            E.append(T/np.sqrt(q*b/(4*np.pi)))
            degree.append(l-1)
            sigma.append(q*b/(4*np.pi*l*(l+1)))
    E = np.asarray(E)
    gram = np.einsum('aij,bji->ab', E, E)
    error = float(np.max(abs(gram-np.eye(q*q-1))))
    assert error < 2e-12, (q, error)
    return E, np.array(degree), np.array(sigma), np.array(beta_l), error


def structure(E):
    comm = -1j*(np.einsum('aij,bjk->abik', E, E)-np.einsum('bij,ajk->abik', E, E))
    f = np.einsum('cji,abij->abc', E, comm)
    assert np.max(abs(f.imag)) < 1e-12
    assert np.max(abs(np.einsum('abc,cij->abij', f, E)-comm)) < 2e-12
    return f.real, comm


def angular_flow_checks():
    rows, bases = [], {}
    for q in [2, 3, 4, 5, 7]:
        E, degree, sigma, beta_l, gram_error = angular_basis(q)
        f, comm = structure(E)
        d = q*q-1
        flat = f.reshape(d*d, d)
        beta = .5*flat.T@((sigma[:, None]*sigma[None, :]).reshape(-1, 1)*flat)
        predicted_diagonal = np.zeros(d)
        C = np.zeros((q-1, q-1, q-1))
        for c in range(q-1):
            out = np.where(degree == c)[0]
            predicted_diagonal[out] = np.diag(beta)[out].mean()
            for a in range(q-1):
                ia = np.where(degree == a)[0]
                for b in range(q-1):
                    ib = np.where(degree == b)[0]
                    C[c, a, b] = .5*np.sum(f[np.ix_(ia, ib, out)]**2)/len(out)
        closure_error = float(np.max(abs(beta-np.diag(predicted_diagonal))))
        assert closure_error < 2e-13
        kap = (sum(1/k for k in range(1, q+1))-1)/(4*np.pi)
        Q = np.einsum('a,aij,ajk->ik', sigma, E, E)
        Qdot = np.einsum('ab,aij,bjk->ik', beta, E, E)
        comm_square = .5*np.einsum('a,b,abij,abjk->ik', sigma, sigma, comm, comm)
        assert np.max(abs(Q-kap*np.eye(q))) < 2e-13
        assert np.max(abs(Qdot-comm_square)) < 2e-13
        eigen_dot = np.linalg.eigvalsh(Qdot)
        assert eigen_dot.min() >= -1e-13 and eigen_dot.max() <= 2*kap*kap+1e-13
        ray_rate = float(np.dot(sigma, np.diag(beta))/np.dot(sigma, sigma))
        ray_error = float(np.linalg.norm(np.diag(beta)-ray_rate*sigma)/np.linalg.norm(beta))
        if q == 2:
            assert ray_error < 1e-12
            assert np.max(abs(np.diag(beta)-2*sigma*sigma)) < 1e-13
        else:
            assert ray_error > .01
        sig0 = np.array([sigma[degree == k][0] for k in range(q-1)])
        rhs = lambda t, s: np.einsum('cab,a,b->c', C, s, s)
        times = np.linspace(0, .3/kap, 61)
        coarse = solve_ivp(rhs, (0, times[-1]), sig0, t_eval=times, rtol=2e-8, atol=1e-11)
        fine = solve_ivp(rhs, (0, times[-1]), sig0, t_eval=times, rtol=2e-11, atol=1e-14)
        assert coarse.success and fine.success
        ref_error = float(np.max(abs(fine.y-coarse.y)))
        assert ref_error < 2e-8
        mult = 2*np.arange(1, q)+1
        kappa = mult@fine.y/q
        upper = kap/(1-2*kap*times)
        assert np.max(kappa-upper) < 1e-12 and np.min(np.diff(fine.y, axis=1)) > -1e-12
        K = sum(s*np.kron(e, e.conj()) for s, e in zip(sigma, E))
        pair_norm = float(np.max(abs(np.linalg.eigvalsh(K))))
        assert abs(pair_norm-kap) < 2e-12
        rows.append(dict(q=q, current_dimension=d, quadrature_gram_error=gram_error,
                         one_loop_closure_error=closure_error, tree_ray_relative_failure=ray_error,
                         initial_kappa=kap, initial_kappa_derivative=float(eigen_dot.mean()),
                         riccati_derivative_ratio=float(eigen_dot.max()/(2*kap*kap)),
                         pair_norm_error=abs(pair_norm-kap), flow_refinement_error=ref_error,
                         flow_end=dict(t=float(times[-1]), coefficients=fine.y[:, -1].tolist(),
                                       kappa=float(kappa[-1]), analytic_upper=float(upper[-1]))))
        bases[q] = (E, sigma)
    return rows, bases


def bilinear(state, create, annihilate):
    if not (state >> annihilate) & 1:
        return None
    sign = -1 if (state & ((1 << annihilate)-1)).bit_count() % 2 else 1
    state ^= 1 << annihilate
    if (state >> create) & 1:
        return None
    if (state & ((1 << create)-1)).bit_count() % 2:
        sign = -sign
    return state | (1 << create), sign


def fock_vacuum(q, N, E, sigma):
    # Exact CAR signs in a finite fermion sea, with NS momenta
    # -N+1/2,...,N-1/2 in each chirality. No current-algebra norm is used
    # to construct these amplitudes.
    def index(chiral, momentum_slot, flavor):
        return (chiral*2*N+momentum_slot)*q+flavor
    sea = sum(1 << index(c, s, i) for c in range(2) for s in range(N) for i in range(q))
    result, energies = {}, {}
    for n in range(1, N+1):
        for e, coupling in zip(E, sigma):
            right = {}
            for p in range(n):
                for i in range(q):
                    for j in range(q):
                        op = bilinear(sea, index(1, N+p, i), index(1, N+p-n, j))
                        if op:
                            state, sign = op
                            right[state] = right.get(state, 0)+sign*e[i, j]
            for state, amplitude in right.items():
                for p in range(n):
                    for i in range(q):
                        for j in range(q):
                            op = bilinear(state, index(0, N+p, i), index(0, N+p-n, j))
                            if op:
                                final, sign = op
                                result[final] = result.get(final, 0)+coupling*amplitude*sign*e[i, j]
                                if final in energies:
                                    assert energies[final] == 2*n
                                energies[final] = 2*n
    norm = float(sum(abs(v)**2 for v in result.values()))
    weighted = float(sum(abs(v)**2/(energies[k]+1) for k, v in result.items()))
    G2 = float(sigma@sigma)
    expected = G2*N*(N+1)*(2*N+1)/6
    expected_weighted = G2*sum(n*n/(2*n+1) for n in range(1, N+1))
    assert abs(norm-expected) < 3e-12 and abs(weighted-expected_weighted) < 3e-12
    return dict(q=q, cutoff=N, fock_modes=4*q*N, nonzero_vacuum_amplitudes=int(sum(abs(v)>1e-14 for v in result.values())),
                norm_squared=norm, form_dual_norm_squared=weighted, current_algebra_norm_error=abs(norm-expected),
                current_algebra_form_error=abs(weighted-expected_weighted))


def asymptotic_checks():
    rows = []
    for q in [2, 8, 32, 128, 512, 2048, 8192, 65536]:
        b, total = 1., 0.
        for l in range(1, q):
            b *= (q-l)/(q+l)
            total += (2*l+1)*b*b/(l*l*(l+1)**2)
        assert 0 < total < 1
        H = sum(1/k for k in range(1, q+1))
        kap = (H-1)/(4*np.pi*np.pi*q*q)
        t = 5*np.log(q)
        assert 2*kap*t < 1
        delta = q**(-1/16)
        N = math.isqrt(q)
        G2 = q**(-2)/(16*np.pi**4)*total
        bound = np.sqrt(G2*N*(N+1)*(2*N+1)/6)/delta
        rows.append(dict(q=q, normalized_frobenius_sum=total, running_fraction_bound=1/(1-2*kap*t)-1,
                         coefficient_norm_over_delta=q*kap/(1-2*kap*t)/delta,
                         running_increment_over_delta=2*q*kap*kap*t/(1-2*kap*t)/delta,
                         finite_cutoff_vacuum_bound=float(bound), scaled_vacuum_bound=float(bound*q**(3/16))))
    assert rows[-1]['normalized_frobenius_sum'] > .999
    assert rows[-1]['coefficient_norm_over_delta'] < rows[2]['coefficient_norm_over_delta']
    return rows


def photon_checks():
    s = sy.symbols('s', real=True)
    rows = []
    for l in [1, 2, 4]:
        for n in [0, 1, 2]:
            a = l+1
            f = sy.cos(s)**a*sy.gegenbauer(n, a, sy.sin(s))
            residual = sy.trigsimp(-sy.diff(f, s, 2)+l*(l+1)/sy.cos(s)**2*f-(n+a)**2*f)
            assert sy.simplify(sy.expand_trig(residual)) == 0, (l, n, residual)
            fn = sy.lambdify(s, f, 'numpy')
            dfn = sy.lambdify(s, sy.diff(f, s), 'numpy')
            norm = quad(lambda x: fn(x)**2, -np.pi/2, np.pi/2, epsabs=1e-10)[0]
            rq = quad(lambda x: dfn(x)**2+l*(l+1)*fn(x)**2/np.cos(x)**2,
                      -np.pi/2, np.pi/2, epsabs=1e-10)[0]/norm
            assert abs(rq-(n+a)**2) < 2e-9
            rows.append(dict(l=l, radial_level=n, symbolic_residual='0', frequency=n+a, rayleigh_error=abs(rq-(n+a)**2)))
    return rows


def main():
    exact = exact_checks()
    angular, bases = angular_flow_checks()
    receipt = dict(status='passed_angular_current_checks', date='2026-10-11', exact_checks=exact,
                   angular_flow_cases=angular, fock_vacuum_cases=[fock_vacuum(q, N, *bases[q]) for q in [2, 3] for N in [1, 2, 3, 4]],
                   asymptotic_cases=asymptotic_checks(), photon_cases=photon_checks(),
                   one_loop_positive_flow_bound_proved=True, single_contact_coupling_closed=False,
                   bare_contact_cutoff_uniform_form_bound=False, full_nonspherical_response_controlled=False,
                   higher_loop_remainder_controlled=False, physical_damping_proved=False, physical_state_constructed=False,
                   novelty_established=False, lean_invocations=0,
                   note_sha256=sha(NOTE), script_sha256=sha(Path(__file__)),
                   dependency_hashes={f: sha(ROOT/f) for f in ['notes/weak-coupling-window.md', 'notes/singlet-static-curvature.md']})
    OUT.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'], counts={k: len(receipt[k]) for k in
                     ['exact_checks', 'angular_flow_cases', 'fock_vacuum_cases', 'asymptotic_cases', 'photon_cases']},
                     ray_failures=[v['tree_ray_relative_failure'] for v in angular],
                     output=str(OUT))), flush=True)


if __name__ == '__main__':
    main()
