"""Exact and finite checks of the projected quadratic photon comparator."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
from scipy.linalg import block_diag
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT/'notes/propagating-photon-polarization.md'
OUT = ROOT/'evidence/propagating-photon-polarization.json'
DEP = ROOT/'checks/angular_current_renormalization.py'
spec = importlib.util.spec_from_file_location('angular', DEP)
angular = importlib.util.module_from_spec(spec)
spec.loader.exec_module(angular)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def exact_checks():
    w, p, a, k, t, g, Q = sy.symbols('w p alpha kappa tau g Q', nonzero=True)
    z, s = w*w-p*p, g*g*t/sy.pi
    v = sy.Matrix([p, w])
    chi = t/(sy.pi*z)*sy.Matrix([[p*p, p*w], [p*w, p*p]])
    S = sy.diag(0, t/sy.pi)
    Sig = sy.simplify(chi+S)
    assert sy.simplify(Sig-t/(sy.pi*z)*v*v.T) == sy.zeros(2)
    assert sy.simplify(Sig*sy.Matrix([w, -p])) == sy.zeros(2, 1)
    C = sy.Matrix([[k+a*p*p, a*w*p], [a*w*p, a*w*w-k]])
    assert sy.factor((C-g*g*Sig).det()-k*(a*z-k-s)) == 0
    dressed = sy.simplify(Sig+g*g*Sig*(C-g*g*Sig).inv()*Sig)
    assert sy.simplify(dressed-Sig*(a*z-k)/(a*z-k-s)) == sy.zeros(2)
    h = sy.symbols('h', real=True)
    mass = k/a*sy.Matrix([[h*h, -h], [-h, 1]])
    prop = sy.factor((Q*Q*sy.eye(2)+mass).inv()[0, 0])
    assert sy.factor(prop-(Q*Q+k/a)/(Q*Q*(Q*Q+k*(1+h*h)/a))) == 0
    assert sy.factor(prop-1/(1+h*h)/Q**2-h*h/(1+h*h)/(Q*Q+k*(1+h*h)/a)) == 0
    return [dict(name=n, exact=True) for n in ['local_contact_and_Ward_identity',
        'Sommerfield_mass_shift', 'dressed_current_residue', 'independent_positive_two_scalar_propagator']]


def frame_constant(Tmin, Tmax):
    mt = max(2/np.pi**2, 2/np.pi*(1/(2*np.pi)+1/Tmin))
    return mt, 4*np.pi*mt*(1+Tmax*Tmax/4)


def pair_energies(Li, Lj, cutoff):
    ni, nj = int(cutoff*Li/(2*np.pi)+1), int(cutoff*Lj/(2*np.pi)+1)
    E = 2*np.pi*((np.arange(ni)[:, None]+.5)/Li+(np.arange(nj)[None, :]+.5)/Lj)
    return E[E <= cutoff]


def frame_checks():
    rows = []
    T = 2.6
    mt, A = frame_constant(T, T)
    for Li, Lj in [(2.61, 2.61), (2.9, 7.3), (3., 91.), (43., 91.)]:
        E = pair_energies(Li, Lj, 20.)
        count_ratio = max(np.count_nonzero(E <= R)/(Li*Lj)*np.pi**2/R**2
                          for R in np.geomspace(.04, 20., 70))
        bin_ratio = max(float(np.sum(1/E[(E >= k)&(E < k+1)])/(Li*Lj*mt)) for k in range(20))
        assert count_ratio <= 1+1e-12 and bin_ratio <= 1+1e-12
        F0 = T*np.sinc(E*T/(2*np.pi))
        p = np.pi/T
        F1 = .5*T*(np.sinc((E-p)*T/(2*np.pi))+np.sinc((E+p)*T/(2*np.pi)))
        ratios = [float(np.sum(abs(F)**2/E)/(Li*Lj*A*norm)) for F, norm in [(F0, T), (F1, T/2)]]
        assert max(ratios) <= 1+1e-12
        rows.append(dict(Li=Li, Lj=Lj, pair_count=len(E), count_bound_ratio=float(count_ratio),
                         unit_bin_bound_ratio=bin_ratio, frame_bound_ratios=ratios))
    return rows


def fock_checks():
    rows = []
    for q in [2, 3]:
        basis, deg, sigma, _, _ = angular.angular_basis(q)
        tau = sigma*(deg+1)*(deg+2)
        matrices = basis*np.sqrt(tau)[:, None, None]
        co0 = np.arange(1, q*q)/(q*q)
        co1 = np.cos(np.arange(q*q-1))
        A0 = np.einsum('a,aij->ij', co0, matrices)
        A1 = np.einsum('a,aij->ij', co1, matrices)
        lengths = 2.6+np.array([.3, 1.2, 4.7][:q])
        for n in [1, 2, 3]:
            modes = [(c, j, r, 2*np.pi*(r+.5)/lengths[j])
                     for c in [1, -1] for j in range(q) for r in range(-n, n)]
            sea = sum(1 << k for k, (_, _, r, _) in enumerate(modes) if r < 0)
            energies = [m[3] for m in modes]
            base = sum(e for e in energies if e < 0)
            amps = {}
            def smear(k, i, j, c):
                # Real test functions: indicator [-1.3,.2] and cos(pi s/2.6).
                f0 = 1.5*np.sinc(k*1.5/(2*np.pi))*np.exp(-.55j*k)
                p = np.pi/2.6
                f1 = 1.3*(np.sinc((k-p)*1.3/np.pi)+np.sinc((k+p)*1.3/np.pi))
                return (A0[i, j]*f0+c*A1[i, j]*f1)/np.sqrt(lengths[i]*lengths[j])
            for ci, (c, i, ri, ei) in enumerate(modes):
                for aj, (d, j, rj, ej) in enumerate(modes):
                    if c != d or ri < 0 or rj >= 0:continue
                    op = angular.bilinear(sea, ci, aj)
                    if op:
                        state, sign = op
                        amps[state] = amps.get(state, 0)+sign*smear(c*(ei-ej), i, j, c)
            w = 1+.23j
            actual = 0j
            for state, amp in amps.items():
                energy = sum(e for k, e in enumerate(energies) if (state >> k)&1)-base
                actual += abs(amp)**2*(1/(w-energy)-1/(w+energy))
            expected = 0j
            for i in range(q):
                for j in range(q):
                    for r in range(n):
                        for s in range(n):
                            energy = 2*np.pi*((r+.5)/lengths[i]+(s+.5)/lengths[j])
                            for c in [1, -1]:
                                expected += abs(smear(c*energy, i, j, c))**2*(1/(w-energy)-1/(w+energy))
            error = float(abs(actual-expected))
            assert error < 2e-12
            rows.append(dict(q=q, ns_modes_per_sign=n, excited_states=len(amps),
                             explicit_CAR_susceptibility_error=error))
    return rows


def photon(N, a, l, w):
    T = 2*a
    h = T/(N+1)
    nodes = -a+h*np.arange(1, N+1)
    edges = -a+h*(np.arange(N+1)+.5)
    G = np.zeros((N+1, N))
    for j in range(N):
        G[j, j] = 1/h
        G[j+1, j] = -1/h
    k = l*(l+1)
    alpha = np.diag(np.cos(nodes)**2)
    H = G.T@G+k*np.diag(1/np.cos(nodes)**2)
    B, Bs = np.hstack([G.T, w*np.eye(N)]), np.vstack([G, -w*np.eye(N)])
    CE = k*np.eye(2*N+1)+Bs@alpha@B
    DE = (np.eye(2*N+1)-Bs@np.linalg.solve(H-w*w*np.eye(N), B))/k
    U = np.r_[np.ones(N+1), -1j*np.ones(N)]
    CM = U[:, None]*CE*U[None, :]
    D = DE/U[:, None]/U[None, :]
    err = float(np.linalg.norm(CM@D-np.eye(2*N+1), 2))
    bound = 2/(1-abs(w)**2/4)
    norm = float(np.linalg.norm(k*D, 2))
    assert err < 2e-11 and norm <= bound
    assert np.linalg.eigvalsh(H)[0] >= 4
    return D, CM, nodes, edges, h, dict(N=N, l=l, w=[w.real, w.imag],
        inverse_error=err, weighted_norm=norm, norm_bound=float(bound), smallest_H=float(np.linalg.eigvalsh(H)[0]))


def pair_current(Li, Lj, cut, w, nodes, edges, h):
    E = pair_energies(Li, Lj, cut)
    # Exact integrals against normalized constant cells; no point-sampling UV alias.
    u0 = np.sqrt(h)*np.sinc(E*h/(2*np.pi))[None, :]*np.exp(1j*edges[:, None]*E)
    u1 = np.sqrt(h)*np.sinc(E*h/(2*np.pi))[None, :]*np.exp(1j*nodes[:, None]*E)
    vr, vl = np.vstack([u0, u1]), np.vstack([u0.conj(), -u1.conj()])
    ans = np.zeros((len(nodes)+len(edges),)*2, dtype=complex)
    for v in [vr, vl]:
        ans += (v/(w-E)[None, :])@v.conj().T-(v.conj()/(w+E)[None, :])@v.T
    return ans/(Li*Lj)


def coupled(q, N, cut, equal=False):
    a, g2, w = 1.3, .002, 1+.15j
    E, degrees, sigmas, betas, gram_err = angular.angular_basis(q)
    ls = degrees+1
    kappas = ls*(ls+1)
    taus = sigmas*kappas
    Ts = E*np.sqrt(taus)[:, None, None]
    photon_cache = {l: photon(N, a, l, w) for l in range(1, q)}
    _, _, nodes, edges, h, _ = photon_cache[1]
    size = 2*N+1
    L = 2*a+np.array(([1.2]*q) if equal else [.3, 1.2, 4.7][:q])
    Sig = np.zeros(((q*q-1)*size,)*2, dtype=complex)
    for i in range(q):
        for j in range(q):
            tij = Ts[:, i, j]
            Pij = pair_current(L[i], L[j], cut, w, nodes, edges, h)
            Sig += np.kron(tij[:, None]*tij.conj()[None, :], Pij)
    contact = np.r_[np.zeros(N+1), np.ones(N)]
    Sig += np.diag(np.repeat(taus/np.pi, size)*np.tile(contact, q*q-1))
    D0 = block_diag(*[photon_cache[int(l)][0] for l in ls])
    C0 = block_diag(*[photon_cache[int(l)][1] for l in ls])
    Dg = np.linalg.inv(C0-g2*Sig)
    mhalf = np.sqrt(np.repeat(kappas, size))
    phalf = np.sqrt(np.repeat(taus, size))
    sj = Sig/phalf[:, None]/phalf[None, :]
    _, AT = frame_constant(2*a, 2*a)
    BJ = 4*AT*(1+abs(w)/w.imag)+1/np.pi
    Cph = 2/(1-abs(w)**2/4)
    eps = g2*Cph*BJ*max(taus/kappas)
    assert eps < 1
    normSig = float(np.linalg.norm(sj, 2))
    errD = float(np.linalg.norm(mhalf[:, None]*(Dg-D0)*mhalf[None, :], 2))
    boundD = Cph*eps/(1-eps)
    assert normSig <= BJ and errD <= boundD
    diag_blocks = block_diag(*[sj[i*size:(i+1)*size, i*size:(i+1)*size] for i in range(q*q-1)])
    mixing = float(np.linalg.norm(sj-diag_blocks, 2))
    # q=2 has only one diagonal traceless vertex; pair symmetry also
    # cancels the x/y cross term. Unequal returns need q>=3 for this test.
    assert mixing < 2e-12 if equal or q == 2 else mixing > 1e-5
    # A fixed continuum smooth A_1 source projected to the node cells.
    probe = np.zeros((q*q-1)*size)
    p = np.pi/(2*a)
    probe[N+1:size] = np.sqrt(h)*np.cos(p*nodes)*np.sinc(p*h/(2*np.pi))*np.sqrt(1/a)
    val = probe.conj()@(Dg-D0)@probe
    return dict(q=q, N=N, cutoff=cut, equal_returns=equal, angular_gram_error=gram_err,
        current_norm=normSig, current_bound=float(BJ), epsilon=float(eps),
        photon_correction_norm=errD, photon_correction_bound=float(boundD),
        angular_mixing_norm=mixing, probe_correction=[float(val.real), float(val.imag)],
        photon_checks=[row[5] for row in photon_cache.values()])


def run():
    exact, frames, fock = exact_checks(), frame_checks(), fock_checks()
    photons = [photon(N, 1.3, l, w)[5] for N in [12, 24, 48]
               for l in [1, 2, 4] for w in [1+.15j, .8+.4j]]
    cases = [coupled(q, 12, 64., eq) for q in [2, 3] for eq in [False, True]]
    uv = [coupled(2, 48, cut) for cut in [16., 32., 64., 128.]]
    mesh = [coupled(2, N, 128.) for N in [12, 24, 48, 96]]
    def difference(rows):
        vs = [complex(*r['probe_correction']) for r in rows]
        return [abs(vs[i]-vs[i-1])/abs(vs[i]) for i in range(1, len(vs))]
    uv_error, mesh_error = difference(uv), difference(mesh)
    assert uv_error[-1] < .01, uv_error
    assert mesh_error[-1] < .01, mesh_error
    scales = []
    for q in [16, 256, 4096, 65536, 2**24]:
        delta = q**(-1/16)
        # The asymptotic estimate uses a fixed band |w|<=1.5 and x_min=.9.
        _, AT = frame_constant(2., 3.5)
        BJ = 4*AT*(1+1.5/delta)+1/np.pi
        Cph = 2/(1-1.5**2/(4*.9**2))
        eps = q**(-2)*Cph*BJ*q*(q-1)/(q+1)/(8*np.pi)
        scales.append(dict(q=q, delta=delta, epsilon_bound=float(eps),
                           scaled_epsilon=float(eps*q**(15/16))))
    r = dict(status='passed_propagating_photon_checks', date='2026-10-11',
        note_sha256=sha(NOTE), script_sha256=sha(Path(__file__)),
        dependency_hashes={str(DEP.relative_to(ROOT)).replace('\\', '/'): sha(DEP)},
        scope='Free projected current polarization and quadratic photon resummation with dual Dirichlet comparator boundaries',
        current_frequency_cutoff_removed_in_ordinary_proof=True,
        full_interacting_stress_controlled=False, physical_matching_proved=False,
        higher_current_remainder_controlled=False, physical_damping_proved=False,
        novelty_established=False, lean_invocations=0,
        exact_checks=exact, frame_cases=frames, fock_cases=fock, photon_cases=photons,
        angular_cases=cases, ultraviolet_refinement=uv, spatial_refinement=mesh,
        ultraviolet_probe_relative_changes=uv_error, spatial_probe_relative_changes=mesh_error,
        probe_acceptance_tolerance=.01, asymptotic_cases=scales)
    OUT.write_text(json.dumps(r, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=r['status'], uv_error=uv_error, mesh_error=mesh_error,
        max_inverse_error=max(x['inverse_error'] for x in photons),
        counts={k:len(r[k]) for k in ['exact_checks','frame_cases','fock_cases','photon_cases','angular_cases','ultraviolet_refinement','spatial_refinement','asymptotic_cases']})))


if __name__ == '__main__':
    run()
