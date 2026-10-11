"""Supplementary potential, propagation and source-insertion checks.

All-mode and unbounded-operator conclusions use the accompanying proof.
Finite matrices here test identities and constants, not a quantum vacuum.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_bvp
from scipy.linalg import expm
from scipy.special import eval_gegenbauer, roots_legendre
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('coupled_maxwell_matching', ROOT/'checks/coupled_maxwell_matching.py')
cm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cm)


def text_hash(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').replace('\r\n', '\n').encode()).hexdigest()


def norm(a):
    return float(np.linalg.norm(a, 2))


def divided_basis(n, c):
    return np.array([eval_gegenbauer(l-1, 1.5, c)/np.sqrt(2*l*(l+1)/(2*l+1))
                     for l in range(1, n+1)])


def throat(ell, length, core, tolerance):
    def ode(t, y):
        mu = ell*(ell+1)*cm.profile(t, length, core)**2
        return np.array([y[2], y[3], mu*y[0], mu*y[1]])
    def bc(a, b):
        return np.array([a[0]-1, a[1], b[0], b[1]-1])
    edge = np.linspace(0, .5, 80)
    grid = np.unique(np.r_[edge, np.linspace(.5, length-.5, 320), length-edge])
    guess = np.array([1-grid/length, grid/length, -np.ones_like(grid)/length, np.ones_like(grid)/length])
    sol = solve_bvp(ode, bc, grid, guess, tol=tolerance, max_nodes=20000)
    assert sol.success, (ell, length, core, sol.message)
    block = np.array([[-sol.y[2, 0], -sol.y[3, 0]], [sol.y[2, -1], sol.y[3, -1]]])
    return block, sol


def field_case(p, n, profile, maps, order):
    e = cm.exterior(p, n, order)
    LT = np.zeros_like(e['LE'])
    for l in range(1, n+1):
        orient = np.diag([1., (-1.)**(l+1)])
        ids = [l-1, n+l-1]
        LT[np.ix_(ids, ids)] = orient@maps[l][0]@orient/p
    forcing = e['LE']@e['v']-e['f']
    h = np.linalg.solve(e['LE']+LT, forcing)
    coef = np.linalg.solve(e['T'], h-e['v'])
    c = np.linspace(-1, 1, 601)
    basis = divided_basis(n, c)
    hv = h.reshape(2, n)@basis
    fv = forcing.reshape(2, n)@basis
    kv = (h-e['v']).reshape(2, n)@basis
    hmax = float(np.max(abs(hv)))
    comparison = hmax/(p*np.max(abs(fv)))
    assert comparison <= 448/445+1e-7, comparison
    length, core = profile
    t = np.unique(np.r_[np.linspace(0, .5, 91), np.linspace(.5, length-.5, 251),
                        length-np.linspace(0, .5, 91)])
    radial = []
    for l in range(1, n+1):
        endpoints = np.array([h[l-1], (-1.)**(l+1)*h[n+l-1]])
        radial.append(endpoints@maps[l][1].sol(t)[:2])
    bv = np.array(radial).T@basis
    maximum_ratio = float(np.max(abs(bv))/hmax)
    assert maximum_ratio <= 1+1e-8
    optical = cm.profile(t, length, core)[:, None]*abs(bv)*np.sqrt(1-c*c)[None, :]/p
    # Independent spatial evaluation of the outgoing two-center multipoles.
    rr = []
    zz = []
    for center in [-.5, .5]:
        rho = p*np.geomspace(1, 4/p, 81)
        rr.append((rho[:, None]*np.sqrt(1-c*c)).ravel())
        zz.append((center+rho[:, None]*c).ravel())
    r = np.concatenate(rr)
    z = np.concatenate(zz)
    rho = np.array([np.hypot(r, z-center) for center in [-.5, .5]])
    keep = np.min(rho, axis=0) >= p*(1-1e-12)
    r, z, rho = r[keep], z[keep], rho[:, keep]
    U = np.zeros_like(r)
    for side, center in enumerate([-.5, .5]):
        ci = (z-center)/rho[side]
        div = divided_basis(n, ci)
        for l in range(1, n+1):
            U += coef[side*n+l-1]*(p/rho[side])**l*div[l-1]/rho[side]**2
    envelope = np.max(abs(kv))/p**2*np.sum((p/rho)**3, axis=0)
    harmonic_ratio = float(np.max(abs(U)/envelope))
    assert harmonic_ratio <= 1+1e-6, harmonic_ratio
    return dict(p=p, cutoff=n, length_over_p=length, core_scaled_w=core,
                boundary_potential_scaled=hmax/p**2, boundary_inverse_ratio=comparison,
                throat_maximum_ratio=maximum_ratio, exterior_harmonic_ratio=harmonic_ratio,
                exterior_optical_potential_scaled=float(np.max(r*abs(U))/p),
                throat_optical_potential_scaled=float(np.max(optical)/p))


def hermitian(rng, n):
    a = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
    a = a+a.conj().T
    return a/norm(a)


def resolvent_derivatives(R, W, Z):
    return [R, -R@W@R, R@W@R@Z@R+R@Z@R@W@R]


def evolution_derivatives(H, t, W, Z):
    n = len(H)
    A = -1j*t*H
    def ordered(directions):
        k = len(directions)+1
        block = np.kron(np.eye(k), A)
        for i, source in enumerate(directions):
            block[i*n:(i+1)*n, (i+1)*n:(i+2)*n] = -1j*t*source
        return expm(block)[:n, (k-1)*n:k*n]
    return [expm(A), ordered([W]), ordered([W, Z])+ordered([Z, W])]


def finite_difference(fn, step):
    first = (fn(step, 0)-fn(-step, 0))/(2*step)
    second = (fn(step, step)-fn(step, -step)-fn(-step, step)+fn(-step, -step))/(4*step**2)
    return [first, second]


def operator_checks():
    rng = np.random.default_rng(10531011)
    evolution, resolvents, ports = [], [], []
    for n, epsilon in itertools.product([6, 12], [.01, .07]):
        H = 2*hermitian(rng, n)
        V = epsilon*hermitian(rng, n)
        W, Z = hermitian(rng, n), hermitian(rng, n)
        assert norm(V[:n//2, n//2:]) > .05*epsilon
        for t in [.2, 2., 7.]:
            values = [evolution_derivatives(H+lam*V, t, W, Z) for lam in [0, 1]]
            ratios = [norm(values[1][j]-values[0][j])/(epsilon*t**(j+1)) for j in range(3)]
            assert max(ratios) <= 1+1e-10
            evals0, vecs0 = np.linalg.eigh(H)
            evals1, vecs1 = np.linalg.eigh(H+V)
            unitary = lambda e, v, time: (v*np.exp(-1j*time*e))@v.conj().T
            integrals = []
            for order in [32, 64]:
                nodes, weights = roots_legendre(order)
                integral = sum(weight*unitary(evals1, vecs1, t*(1-node)/2)@V@
                               unitary(evals0, vecs0, t*(1+node)/2)
                               for node, weight in zip(nodes, weights))*(-1j*t/2)
                integrals.append(integral)
            identity_error = norm(values[1][0]-values[0][0]-integrals[-1])
            assert identity_error < 1e-11
            fd_errors = []
            for step in [.002/t, .001/t]:
                fd = finite_difference(lambda x, y: expm(-1j*t*(H+x*W+y*Z)), step)
                fd_errors.append(max(norm(fd[j-1]-values[0][j])/(1+norm(values[0][j])) for j in [1, 2]))
            assert fd_errors[-1] < 3e-6
            assert fd_errors[-1] < .4*fd_errors[0]+3e-8
            evolution.append(dict(dimension=n, epsilon=epsilon, time=t, bound_ratios=ratios,
                                  duhamel_identity_error=identity_error,
                                  integral_refinement=norm(integrals[0]-integrals[1]),
                                  source_difference_errors=fd_errors))
        B = rng.normal(size=(n, 3))+1j*rng.normal(size=(n, 3))
        B /= norm(B)
        for nu in [.03, .2, .9]:
            z = .17+1j*nu
            values = [resolvent_derivatives(np.linalg.inv(H+lam*V-z*np.eye(n)), W, Z) for lam in [0, 1]]
            ratios = [norm(values[1][j]-values[0][j])/(math.factorial(j+1)*epsilon/nu**(j+2)) for j in range(3)]
            assert max(ratios) <= 1+1e-10
            R0, R1 = values[0][0], values[1][0]
            error = norm(R1-R0+R1@V@R0)
            assert error < 1e-10
            resolvents.append(dict(dimension=n, epsilon=epsilon, nu=nu, bound_ratios=ratios,
                                    resolvent_identity_error=error))
            port_values = []
            defects = []
            for lam in [0, 1]:
                A = (nu-.17j)*np.eye(n)+1j*(H+lam*V)+B@B.conj().T/2
                Q = np.linalg.inv(A)
                P = Q@B
                T = np.eye(3)-B.conj().T@P
                defects.append(norm(T.conj().T@T+2*nu*P.conj().T@P-np.eye(3)))
                assert norm(P) <= 1/np.sqrt(2*nu)+1e-10
                derivatives = resolvent_derivatives(Q, 1j*W, 1j*Z)
                port_values.append([T]+[-B.conj().T@derivatives[j]@B for j in [1, 2]])
            ratios = [norm(port_values[1][j]-port_values[0][j])/(math.factorial(j+1)*epsilon/(2*nu**(j+1)))
                      for j in range(3)]
            assert max(ratios) <= 1+1e-10
            assert max(defects) < 1e-11
            ports.append(dict(dimension=n, epsilon=epsilon, nu=nu, bound_ratios=ratios,
                              current_identity_error=max(defects)))
    return evolution, resolvents, ports


def attenuation_checks():
    # Independent scalar barrier tests; these round necks are not asserted
    # to have the Euclidean seam collars or to solve Einstein equations.
    cases = []
    for kappa in [.5, .1, .01, .001]:
        edge = np.sqrt(1-kappa*kappa)
        def ode(t, v):
            return np.array([v[1], 2*v[0]/(t*t+kappa*kappa)])
        def bc(a, b):
            return np.array([a[0]-1, b[0]-1])
        grids = [np.linspace(-edge, edge, 240), kappa*np.linspace(-10, 10, 120)]
        grid = np.unique(np.concatenate(grids))
        grid = grid[(grid >= -edge) & (grid <= edge)]
        solutions = []
        for tolerance in [1e-8, 1e-10]:
            guess = np.array([np.ones_like(grid), np.zeros_like(grid)])
            sol = solve_bvp(ode, bc, grid, guess, tol=tolerance, max_nodes=24000)
            assert sol.success, (kappa, sol.message)
            solutions.append(sol.sol(grid)[0])
        exact = grid*grid+kappa*kappa
        error = float(np.max(abs(solutions[-1]-exact)))
        assert error < 1e-8
        y = np.sqrt(exact)
        weighted = float(np.max(abs(solutions[-1])/y))
        assert weighted <= 1+1e-8
        cases.append(dict(kappa=kappa, p_max_w=1/kappa, optical_potential_bound=weighted,
                          exact_solution_error=error,
                          bvp_refinement=float(np.max(abs(solutions[0]-solutions[1])))))
    return cases


def verify():
    r, z, x, c, p, gamma = sp.symbols('r z x c p gamma', positive=True)
    U = sp.Function('U')(r, z)
    a = r*r*U
    assert sp.simplify(sp.diff(a, r, 2)-sp.diff(a, r)/r+sp.diff(a, z, 2)
                       -r*r*(sp.diff(U, r, 2)+3*sp.diff(U, r)/r+sp.diff(U, z, 2))) == 0
    b = sp.Function('b')(x, c)
    a = (1-c*c)*b
    assert sp.expand(sp.diff(a, c, 2)-((1-c*c)*sp.diff(b, c, 2)-4*c*sp.diff(b, c)-2*b)) == 0
    barrier = 8*gamma/7+(1-8*gamma/7)*(p/r)**3
    assert sp.simplify(sp.diff(barrier, r, 2)+4*sp.diff(barrier, r)/r) == 0
    assert sp.simplify(barrier.subs(r, p)-1) == 0
    assert sp.simplify(barrier.subs(r, 2*p)-(sp.Rational(1, 8)+gamma)) == 0
    assert sp.simplify(-p*sp.diff(barrier, r).subs(r, p)-2-(1-24*gamma/7)) == 0
    assert 1-sp.Rational(24, 7)/512 == sp.Rational(445, 448)
    gradient_barrier = 1-sp.exp(-8*(r-p)/p)
    assert sp.simplify(sp.diff(gradient_barrier, r, 2)+4*sp.diff(gradient_barrier, r)/r
                       -8*sp.exp(-8*(r-p)/p)/p**2*(-8+4*p/r)) == 0
    re, delta = sp.symbols('re delta', positive=True)
    f = (1-re/r)**2
    y = r*r/(r-re)
    dx = lambda value: sp.simplify(f*sp.diff(value, r))
    assert sp.simplify(dx(y)-(1-2*re/r)) == 0
    assert sp.simplify(y*dx(dx(y))-2*re*(r-re)/r**2) == 0
    # Difference from the bound 1/2 is an explicit nonnegative square.
    assert sp.simplify(sp.Rational(1, 2)-y*dx(dx(y))-(r-2*re)**2/(2*r*r)) == 0
    ads_y = sp.cos(delta*x)/delta
    assert sp.simplify(ads_y*sp.diff(ads_y, x, 2)+sp.cos(delta*x)**2) == 0
    # Exact static l=1 test-Maxwell solution on the prescribed RN metric.
    rn_mode = r*r-re*re
    assert sp.simplify(sp.diff(f*sp.diff(rn_mode, r), r)-2*rn_mode/r**2) == 0
    scale_powers = []
    for j in range(3):
        epsilon_power = 1-sp.Rational(3, 2)
        assert epsilon_power == -sp.Rational(1, 2)
        scale_powers.append(dict(source_order=j, q_power=str(epsilon_power),
                                 resolvent_log_power=4*(j+2), return_log_power=4*(j+1),
                                 finite_time_log_power=2*(j+1)))
    profiles = [(2., .2), (6., 1.), (20., .05), (1000., .005)]
    fields, refinements = [], []
    for profile in profiles:
        maps = {l: throat(l, *profile, 1e-10) for l in range(1, 9)}
        coarse = throat(1, *profile, 1e-8)
        change = norm(coarse[0]-maps[1][0])
        assert change < 1e-7
        refinements.append(dict(length_over_p=profile[0], dtn_refinement=change))
        for radius in [.1, .03, .003]:
            low = field_case(radius, 6, profile, maps, 160)
            high = field_case(radius, 8, profile, maps, 240)
            change = max(abs(high[k]-low[k]) for k in ['boundary_potential_scaled',
                'exterior_optical_potential_scaled', 'throat_optical_potential_scaled'])
            assert change < 1e-5
            high['angular_refinement'] = change
            fields.append(high)
    evolution, resolvents, ports = operator_checks()
    attenuation = attenuation_checks()
    inputs = ['notes/maxwell-dirac-comparison.md', 'notes/coupled-maxwell-matching.md',
              'notes/loaded-operator-realization.md', 'checks/coupled_maxwell_matching.py']
    return dict(status='passed_maxwell_dirac_comparison_checks', exact_lifts_checked=2,
                exact_barrier_coefficient='445/448', field_cases=fields, throat_refinements=refinements,
                scale_powers=scale_powers,
                rn_ads_optical_comparison_checked=True, attenuation_cases=attenuation,
                evolution_cases=evolution, resolvent_cases=resolvents, loaded_port_cases=ports,
                optical_condition='sup(w)<=C/p OR y=1/w satisfies y_double_prime<=2/y',
                physical_rn_neck_violates_fixed_cap=True,
                physical_global_metric_matching_open=True,
                potential_error_order='q p', selected_error_order='q^(-1/2)',
                physical_metric_solved=False, quantum_vacuum_response_proved=False,
                novelty_established=False, supplementary_numerics=True, lean_invocations=0,
                input_sha256={v: text_hash(ROOT/v) for v in inputs},
                versions=dict(numpy=np.__version__, scipy=scipy.__version__, sympy=sp.__version__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    receipt = verify()
    receipt['script_sha256'] = text_hash(Path(__file__))
    Path(args.output).write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: receipt[k] if k == 'status' else len(receipt[k])
                      for k in ['status', 'field_cases', 'evolution_cases', 'resolvent_cases', 'loaded_port_cases']}))
