"""Supplementary scalar ODE, return-matrix and variational graph checks.

No Lean, cloud invocation, physical damping claim or numerical proof.
Run from any directory; the receipt is tied to this file and its note.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import scipy
import sympy as sp
from numpy.linalg import det, norm, solve, svd
from scipy.integrate import solve_bvp
from scipy.linalg import eigh, null_space
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
J = np.array([[0., 1.], [1., 0.]])
I2 = np.eye(2)


def require(value, message):
    if not value:
        raise AssertionError(message)


def relerr(a, b):
    return float(norm(a-b)/max(1., norm(b)))


def exact_checks():
    R, lapse, g, q, E, phi = sp.symbols('R lapse g q E phi', positive=True)
    e2 = g**2*lapse**2/(4*sp.pi*R**2)
    electric = E**2/(2*e2)+E*sp.sqrt(q/sp.pi)*phi
    stationary = -e2*sp.sqrt(q/sp.pi)*phi
    residual = sp.simplify(electric.subs(E, stationary)+q*e2*phi**2/(2*sp.pi))
    require(residual == 0, 'electric square normalization')
    nu = sp.Rational(15, 16)
    exponent = sp.simplify(-1-5*(nu-1))
    require(exponent == -sp.Rational(11, 16), 'finite-resolution exponent')
    n, r = sp.symbols('n r', integer=True, positive=True)
    pair_sum = sp.simplify(sp.summation(r*(n-r), (r, 1, n-1))-n*(n*n-1)/6)
    require(pair_sum == 0, 'boson pair sum')
    return [dict(identity='electric square', residual=str(residual)),
            dict(identity='relative boundary-readout exponent at nu=15/16', value=str(exponent)),
            dict(identity='ordered positive boson pair sum', residual=str(pair_sum))]


def scalar_load(z, T, mass, tol=2e-9):
    """Two independent prescribed-endpoint columns, not a transfer shooting solve."""
    mesh = np.linspace(0., T, 65)
    real_initial = np.vstack((1-mesh/T, -np.ones_like(mesh)/T,
                              mesh/T, np.ones_like(mesh)/T))
    initial = np.vstack((real_initial, np.zeros_like(real_initial)))

    def ode(x, y):
        coefficient = z*z+mass(x)**2
        c = y[:4]+1j*y[4:]
        rhs = np.vstack((c[1], coefficient*c[0], c[3], coefficient*c[2]))
        return np.vstack((rhs.real, rhs.imag))

    def bc(left, right):
        lc, rc = left[:4]+1j*left[4:], right[:4]+1j*right[4:]
        residual = np.array([lc[0]-1, rc[0], lc[2], rc[2]-1])
        return np.concatenate((residual.real, residual.imag))

    result = solve_bvp(ode, bc, mesh, initial, tol=tol, max_nodes=18000)
    require(result.success, result.message)
    def solution(x):
        values = result.sol(x)
        return values[:4]+1j*values[4:]

    ends = solution(np.array([0., T]))
    load = np.array([[-ends[1, 0], -ends[3, 0]], [ends[1, 1], ends[3, 1]]])
    scattering = solve(load+z*I2, z*I2-load)
    return load, scattering, SimpleNamespace(x=result.x, sol=solution)


def analytic_load(z, T, m):
    a = np.sqrt(complex(z*z+m*m))
    return a*np.array([[1/np.tanh(a*T), -1/np.sinh(a*T)],
                       [-1/np.sinh(a*T), 1/np.tanh(a*T)]])


def ode_checks():
    T = 1.3
    constant_rows = []
    for m in [0., 1., 4., 10.]:
        for z in [.2-.7j, .03-3.2j, 1.1-.3j]:
            load, scattering, result = scalar_load(z, T, lambda x: m+0*x)
            exact = analytic_load(z, T, m)
            error = relerr(load, exact)
            require(error < 2e-8, ('constant ODE', m, z, error))
            require(norm(scattering, 2) <= 1+2e-9, 'passive constant load')
            if m == 0:
                require(relerr(scattering, np.exp(-z*T)*J) < 2e-8, 'free propagation')
            constant_rows.append(dict(m=m, z=[z.real, z.imag], relative_error=error,
                                      scattering_norm=float(norm(scattering, 2)), nodes=len(result.x)))
    variable_rows = []
    nodes, weights = leggauss(320)
    xs, weights = T*(nodes+1)/2, T*weights/2
    b = np.array([.7+.2j, -.3+.4j])
    for minimum in [1.5, 5., 12.]:
        mass = lambda x: minimum*(1+.3*np.sin(np.pi*x/T)**2+.2*x/T)
        for z in [.1-.8j, .01-2.1j]:
            load, scattering, result = scalar_load(z, T, mass)
            finer, _, refined = scalar_load(z, T, mass, tol=2e-10)
            y = result.sol(xs)
            f, fp = b[0]*y[0]+b[1]*y[2], b[0]*y[1]+b[1]*y[3]
            lhs = float(np.real(np.vdot(b, load@b)/z))
            rhs = float(z.real*np.dot(weights, abs(f)**2)
                        +(1/z).real*np.dot(weights, abs(fp)**2+mass(xs)**2*abs(f)**2))
            identity_error = abs(lhs-rhs)/max(1., abs(rhs))
            refinement_error = relerr(load, finer)
            require(identity_error < 2e-8, 'passivity energy identity')
            require(refinement_error < 2e-8, 'variable ODE refinement')
            require(norm(scattering, 2) <= 1+2e-9, 'passive variable load')
            variable_rows.append(dict(minimum_mass=minimum, z=[z.real, z.imag],
                identity_error=identity_error, refinement_error=refinement_error,
                scattering_norm=float(norm(scattering, 2)), nodes=len(result.x), refined_nodes=len(refined.x)))
    reflection_rows = []
    for minimum in [2., 5., 12.]:
        mass = lambda x: minimum*(1+.3*np.sin(np.pi*x/T)**2+.2*x/T)
        for fraction in [.05, .2, .7]:
            k = minimum*fraction
            load, scattering, _ = scalar_load(-1j*k, T, mass)
            a = np.sqrt(minimum**2-k**2)
            lower = a*np.tanh(a*T/2)
            observed = float(norm(scattering+I2, 2))
            bound = 2*k/lower
            require(np.min(eigh(load.real, eigvals_only=True)) >= lower-2e-8, 'DtN lower bound')
            unitary_error = relerr(scattering.conj().T@scattering, I2)
            require(unitary_error < 2e-8, 'real-frequency unitary')
            require(observed <= bound+2e-8, 'reflection bound')
            reflection_rows.append(dict(minimum_mass=minimum, k=k, norm_error_from_reflection=observed,
                                        bound=bound, unitary_error=unitary_error))
    return constant_rows, variable_rows, reflection_rows


def full_return(q, T, lengths, z, singlet):
    v = np.ones((q, 1))/np.sqrt(q)
    P, N = v@v.T, np.eye(q)-v@v.T
    s0 = np.exp(-z*T)*J
    throat = np.kron(s0, N)+np.kron(singlet, P)
    free = np.kron(s0, np.eye(q))
    exterior = np.kron(J, np.diag(np.exp(-z*lengths)))
    identity = np.eye(2*q)
    G, G0 = solve(identity-exterior@throat, identity), solve(identity-exterior@free, identity)
    U = np.kron(I2, v)
    C = singlet-s0
    # Woodbury in an order which does not invert possibly singular C.
    woodbury = G0+G0@exterior@U@solve(I2-C@U.T@G0@exterior@U, C@U.T@G0)
    return G, G0, woodbury, throat, exterior


def return_checks():
    rows, determinant_rows = [], []
    T = 1.3
    for q in [2, 5, 17]:
        lengths = .4+3*np.linspace(0., 1., q)**1.3
        for z in [.3-.8j, .015-2.4j, .004-4.7j]:
            load = analytic_load(z, T, 4.)
            sm = solve(load+z*I2, z*I2-load)
            for kind, singlet in [('finite_constant_mass', sm), ('reflecting', -I2)]:
                G, G0, woodbury, throat, exterior = full_return(q, T, lengths, z, singlet)
                diff = G-G0
                singular = svd(diff, compute_uv=False)
                tail = float(norm(singular[2:]))
                r = np.exp(-z.real*min(lengths))
                trace_norm, bound = float(sum(singular)), float(4*r/(1-r)**2)
                error = relerr(woodbury, G)
                require(error < 2e-11, 'Woodbury comparison')
                require(tail < 2e-10*max(1., singular[0]), 'rank two')
                require(trace_norm <= bound*(1+2e-11), 'trace norm bound')
                rows.append(dict(q=q, z=[z.real, z.imag], kind=kind, rank_tail=tail,
                    woodbury_error=error, trace_norm=trace_norm, trace_norm_bound=bound))
            _, _, _, throat, exterior = full_return(q, T, lengths, z, -I2)
            t, ds = np.exp(-z*T), np.exp(-z*lengths)
            product = 1.+0j
            parity_errors = []
            P = np.ones((q, q))/q
            for sigma in [-1, 1]:
                direct = det(np.eye(q)-t*np.diag(ds)+(t+sigma)*np.diag(ds)@P)
                reduced = np.prod(1-t*ds)*(1+(t+sigma)/q*np.sum(ds/(1-t*ds)))
                parity_errors.append(float(abs(direct-reduced)/max(1., abs(direct))))
                product *= reduced
            error = float(abs(product-det(np.eye(2*q)-exterior@throat))/max(1., abs(product)))
            require(max(parity_errors+[error]) < 2e-11, 'secular determinant')
            determinant_rows.append(dict(q=q, z=[z.real, z.imag], parity_errors=parity_errors,
                                         full_determinant_error=error))
    # Unequal lengths: a free pole is no longer a pole of the coupled system.
    counterexample = []
    for eta in [.01, .001, .0001]:
        k = 2*np.pi/3
        G, G0, _, _, _ = full_return(2, 1., np.array([2., 3.7]), eta-1j*k, -I2)
        counterexample.append(dict(eta=eta, difference_norm=float(norm(G-G0, 2)),
                                   coupled_norm=float(norm(G, 2))))
    require(counterexample[-1]['difference_norm'] > 50*counterexample[0]['difference_norm'],
            'real-axis nonuniformity example')
    require(max(x['coupled_norm'] for x in counterexample) < 20., 'shifted pole example')
    return rows, determinant_rows, counterexample


def assemble_circles(T, lengths, n):
    q, per = len(lengths), 2*n
    K, M = np.zeros((q*per, q*per)), np.zeros((q*per, q*per))
    KT, MT = np.zeros_like(K), np.zeros_like(M)
    for j, D in enumerate(lengths):
        for edge in range(per):
            h = (T if edge < n else D)/n
            ix = np.array([j*per+edge, j*per+(edge+1) % per])
            ke = np.array([[1., -1.], [-1., 1.]])/h
            me = h*np.array([[2., 1.], [1., 2.]])/6
            K[np.ix_(ix, ix)] += ke
            M[np.ix_(ix, ix)] += me
            if edge < n:
                KT[np.ix_(ix, ix)] += ke
                MT[np.ix_(ix, ix)] += me
    C = np.zeros((2, q*per))
    CI = np.zeros((n+1, q*per))
    for j in range(q):
        C[0, j*per], C[1, j*per+n] = 1., 1.
        for s in range(n+1):
            CI[s, j*per+s] = 1.
    return K, M, C, CI, KT, MT


def constrained_spectrum(K, M, C):
    Q = null_space(C)
    return eigh(Q.T@K@Q, Q.T@M@Q, eigvals_only=True)


def finite_element_checks():
    rows, equal_rows = [], []
    T = 1.3
    for lengths in [np.array([.7, 1.4, 2.8]), np.array([1.4]*3)]:
        for n in [8, 16, 32]:
            K, M, C, CI, _, _ = assemble_circles(T, lengths, n)
            free = eigh(K, M, eigvals_only=True)
            cut = constrained_spectrum(K, M, C)
            limiting = constrained_spectrum(K, M, CI)
            h = T/n
            KD = (2*np.eye(n-1)-np.eye(n-1, k=1)-np.eye(n-1, k=-1))/h
            MD = h*(4*np.eye(n-1)+np.eye(n-1, k=1)+np.eye(n-1, k=-1))/6
            interval = eigh(KD, MD, eigvals_only=True)
            split = np.sort(np.concatenate([limiting, interval]))
            split_error = float(np.max(abs(split-cut)/(1+abs(cut))))
            lower_violation = float(max(0., np.max((free[:-2]-cut)/(1+abs(cut)))))
            upper_violation = float(max(0., np.max((cut-free[2:])/(1+abs(cut)))))
            require(split_error < 2e-9, 'independent form splitting')
            require(max(lower_violation, upper_violation) < 2e-9, 'codimension two interlacing')
            require(sum(abs(free) < 1e-7) == len(lengths), 'free kernel')
            require(sum(abs(limiting) < 1e-7) == len(lengths)-1, 'reflecting graph kernel')
            # Positive-frequency counts with the discrete Dirichlet interval.
            for omega in [.8, 2., 4.5, 7.]:
                nf = int(sum((free > 1e-7) & (free <= omega**2)))
                ni = int(sum((limiting > 1e-7) & (limiting <= omega**2)))
                nd = int(sum(interval <= omega**2))
                require(nf-nd-1 <= ni <= nf-nd+1, 'positive count correction')
            rows.append(dict(lengths=lengths.tolist(), cells_per_edge=n,
                split_error=split_error, lower_violation=lower_violation, upper_violation=upper_violation,
                free_kernel=int(sum(abs(free) < 1e-7)), limiting_kernel=int(sum(abs(limiting) < 1e-7))))
            if np.all(lengths == lengths[0]):
                expected = []
                for mode in range(1, 15):
                    expected += [(2*np.pi*mode/(T+lengths[0]))**2]*(2*(len(lengths)-1))
                    expected += [(np.pi*mode/lengths[0])**2]
                expected = np.sort(expected)[:12]
                observed = limiting[limiting > 1e-7][:12]
                error = float(max(abs(observed-expected)/expected))
                equal_rows.append(dict(cells_per_edge=n, first_twelve_relative_error=error))
    require(equal_rows[2]['first_twelve_relative_error'] < .013, 'equal-length limiting spectrum')
    require(equal_rows[0]['first_twelve_relative_error']/equal_rows[1]['first_twelve_relative_error'] > 3.7,
            'first mesh refinement')
    require(equal_rows[1]['first_twelve_relative_error']/equal_rows[2]['first_twelve_relative_error'] > 3.7,
            'second mesh refinement')
    return rows, equal_rows


def pair_weights(frequencies, A, B):
    products = np.outer(frequencies, frequencies)
    difference = B-products*A
    return difference**2/(8*products)


def band_weight(output_frequency):
    x = (output_frequency-4.8)/1.8
    out = np.zeros_like(x)
    mask = abs(x) < 1
    out[mask] = np.exp(-x[mask]**2/(1-x[mask]**2))
    return out


def stress_checks():
    """Independent Fourier modes, local quadratic forms and vacuum-energy Hessian."""
    T, L = 1.3, 3.2
    nx, wx = leggauss(220)
    xs, weights = T*(nx+1)/2, T*wx/2
    modes, derivatives, frequencies = [], [], []
    for m in range(1, 8):
        k = 2*np.pi*m/L
        modes.extend([np.sqrt(2/L)*np.cos(k*xs), np.sqrt(2/L)*np.sin(k*xs)])
        derivatives.extend([-np.sqrt(2/L)*k*np.sin(k*xs), np.sqrt(2/L)*k*np.cos(k*xs)])
        frequencies.extend([k, k])
    modes, derivatives, frequencies = np.array(modes), np.array(derivatives), np.array(frequencies)
    A, B = (modes*weights)@modes.T, (derivatives*weights)@derivatives.T
    W = pair_weights(frequencies, A, B)
    output = frequencies[:, None]+frequencies[None, :]
    fourier = []
    for n in range(1, 9):
        observed = float(sum(W[np.isclose(output, 2*np.pi*n/L)]))
        expected = 2*(n*n-1)*np.sin(np.pi*n*T/L)**2/(3*L*L*n)
        error = abs(observed-expected)
        require(error < 2e-11, 'free-fermion weights from boson pairs')
        fourier.append(dict(harmonic=n, observed=observed, expected=float(expected), error=error))

    rows, contacts = [], []
    for lengths in [np.array([.7, 1.4, 2.8]), np.array([1.4]*3)]:
        shape = 'equal' if np.all(lengths == lengths[0]) else 'unequal'
        for n in [12, 24, 48, 96]:
            K, M, _, CI, KT, MT = assemble_circles(T, lengths, n)
            one_row = dict(shape=shape, cells_per_edge=n)
            kinds = ['free', 'reflecting']+(['finite_mass'] if n == 24 else [])
            for kind in kinds:
                Q = null_space(CI) if kind == 'reflecting' else np.eye(len(K))
                # Constant singlet mass 4, with its vertex retained in V.
                potential = 16*np.kron(np.ones((3, 3))/3, MT[:2*n, :2*n]) if kind == 'finite_mass' else np.zeros_like(K)
                full_stiffness, throat_stiffness = K+potential, KT+potential
                eig, basis = eigh(Q.T@full_stiffness@Q, Q.T@M@Q)
                # The output band is (3,6.6); each positive contributing mode is below 6.6.
                selected = (eig > 1e-7) & (eig < 6.6**2)
                ws, U = np.sqrt(eig[selected]), Q@basis[:, selected]
                A, B = U.T@MT@U, U.T@throat_stiffness@U
                W = pair_weights(ws, A, B)
                value = float(np.sum(W*band_weight(ws[:, None]+ws[None, :])))
                one_row[kind+'_band_weight'] = value
                whole = pair_weights(ws, U.T@M@U, U.T@full_stiffness@U)
                require(np.max(abs(whole)) < 2e-19, 'whole Hamiltonian creates no pairs')
                if kind != 'free' and n == 24:
                    # Finite oscillator family: H(phi)=H(0)+(exp(phi)-1)V.
                    omega2 = np.diag(ws**2)
                    def vacuum_energy(phi):
                        shift = np.expm1(phi)
                        kinetic = np.eye(len(ws))+shift*A
                        potential = omega2+shift*B
                        es, vs = eigh(kinetic)
                        square = (vs*np.sqrt(es))@vs.T
                        return .5*sum(np.sqrt(eigh(square@potential@square, eigvals_only=True)))
                    contact = float(.25*sum(ws*np.diag(A)+np.diag(B)/ws))
                    static = float(-2*np.sum(W/(ws[:, None]+ws[None, :])))
                    predicted = contact+static
                    fd = []
                    for h in [.002, .001]:
                        observed = (vacuum_energy(h)-2*vacuum_energy(0)+vacuum_energy(-h))/h**2
                        fd.append(dict(step=h, energy_hessian=float(observed), error=float(abs(observed-predicted))))
                    require(fd[-1]['error'] < 3e-5*max(1., abs(predicted)), 'contact plus retarded static response')
                    contacts.append(dict(shape=shape, kind=kind, oscillators=len(ws), bare_contact=contact,
                                         static_susceptibility=static, predicted_hessian=predicted, differences=fd))
            one_row['naive_q_minus_one_band_weight'] = 2*one_row['free_band_weight']/3
            one_row['relative_naive_error'] = abs(one_row['reflecting_band_weight']-one_row['naive_q_minus_one_band_weight'])/one_row['reflecting_band_weight']
            if shape == 'equal':
                expected = 0.
                circle = T+lengths[0]
                for harmonic in range(2, 15):
                    wn = 2*(harmonic*harmonic-1)*np.sin(np.pi*harmonic*T/circle)**2/(3*circle*circle*harmonic)
                    expected += 2*wn*float(band_weight(np.array([2*np.pi*harmonic/circle]))[0])
                one_row['continuum_equal_band_weight'] = float(expected)
                one_row['relative_continuum_error'] = abs(one_row['reflecting_band_weight']-expected)/expected
                require(one_row['relative_naive_error'] < 2e-9, 'equal delays: q-1 response')
            rows.append(one_row)
    unequal = [x for x in rows if x['shape'] == 'unequal']
    equal = [x for x in rows if x['shape'] == 'equal']
    require(unequal[-1]['relative_naive_error'] > .03, 'unequal-delay counterexample')
    changes = [abs(unequal[i+1]['reflecting_band_weight']-unequal[i]['reflecting_band_weight']) for i in range(len(unequal)-1)]
    require(changes[-2]/changes[-1] > 3.5, 'unequal stress mesh convergence')
    require(changes[-1]/unequal[-1]['reflecting_band_weight'] < .025, 'unequal stress mesh resolution')
    require(equal[-1]['relative_continuum_error'] < .02, 'equal stress continuum convergence')
    return fourier, rows, contacts


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    exact = exact_checks()
    constant, variable, reflection = ode_checks()
    returns, determinants, counterexample = return_checks()
    forms, spectra = finite_element_checks()
    fourier, stress, contacts = stress_checks()
    result = dict(status='passed_schwinger_return_feedback_checks',
        scope='defined Gaussian graph; supplementary floating point verification',
        script_sha256=sha(Path(__file__)), note_sha256=sha(ROOT/'notes/schwinger-return-feedback.md'),
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, sympy=sp.__version__,
        exact_checks=exact, constant_load_cases=constant, variable_load_cases=variable,
        reflection_cases=reflection, return_cases=returns, determinant_cases=determinants,
        real_axis_counterexample=counterexample, finite_element_cases=forms, equal_length_refinements=spectra,
        free_stress_harmonics=fourier, coupled_stress_refinements=stress, contact_term_cases=contacts,
        oscillator_stress_formula_derived=True,
        physical_damping_proved=False, interacting_stress_response_controlled=False,
        full_physical_model_closed=False, novelty_established=False, lean_invocations=0)
    path = ROOT/'evidence/schwinger-return-feedback.json'
    path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], constant_load_cases=len(constant), variable_load_cases=len(variable),
        reflection_cases=len(reflection), return_cases=len(returns), determinant_cases=len(determinants),
        finite_element_cases=len(forms), max_ode_error=max(x['relative_error'] for x in constant),
        max_splitting_error=max(x['split_error'] for x in forms), equal_length_refinements=spectra,
        stress_refinements=stress, contact_term_cases=contacts), indent=2))


if __name__ == '__main__':
    main()
