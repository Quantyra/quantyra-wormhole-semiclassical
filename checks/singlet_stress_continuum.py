"""Independent checks for the uniform Gaussian singlet absorption limit.

The analytic note proves the uniform bound. Numerical cases are supplementary.
The first-order energy representation is bookkeeping, not a new physical field.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
import sympy as sp
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.linalg import eigh, null_space

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('schwinger_return_feedback', ROOT/'checks/schwinger_return_feedback.py')
DEPENDENCY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DEPENDENCY)
assemble_circles, pair_weights, band_weight = DEPENDENCY.assemble_circles, DEPENDENCY.pair_weights, DEPENDENCY.band_weight
PI = np.pi
KAPPA = 1+np.arcsinh(np.sqrt(3))/(2*np.sqrt(3))
S = 8*np.log(2)-13*PI**2/24
CSTAR = KAPPA*S/(3*PI**2)
WP = Poly([693/512])*Poly([1, 0, -1])**5
BP = Poly([0, 0, 0, 0, 0, 126, -420, 540, -315, 70])


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def bump(t, order=0):
    t = np.asarray(t)
    return np.where(abs(t) < 1, WP.deriv(order)(t), 0.)


def cutoff(t, order=0):
    t = np.asarray(t)
    middle = (t > .5) & (t < 1.)
    return np.where(middle, 2**order*BP.deriv(order)(np.clip(2*t-1, 0., 1.)),
                    np.where(t >= 1., 1., 0.) if order == 0 else 0.)


def exact_checks():
    x, y, a, b, c, cp = sp.symbols('x y a b c cp', real=True)
    # Independent formal differentiation of chi(y)[f(x)+iy f'(x)].
    f, chi = sp.Function('f')(x), sp.Function('chi')(y)
    ext = chi*(f+sp.I*y*sp.diff(f, x))
    expected = sp.I*sp.diff(chi, y)*f/2-y*sp.diff(chi, y)*sp.diff(f, x)/2+sp.I*y*chi*sp.diff(f, x, 2)/2
    require(sp.simplify((sp.diff(ext, x)+sp.I*sp.diff(ext, y))/2-expected) == 0, 'Cauchy derivative')
    w = sp.Rational(693, 512)*(1-x*x)**5
    require(sp.integrate(w, (x, -1, 1)) == 1, 'kernel normalization')
    require(all(sp.diff(w, x, n).subs(x, 1) == 0 for n in range(5)), 'kernel C4 boundary')
    omega, nu = sp.symbols('omega nu', positive=True)
    require(sp.simplify((b-omega*nu*a)**2/(8*omega*nu)-omega*nu*(a-b/(omega*nu))**2/8) == 0, 'pair normalization')
    exponent = sp.Rational(1, 16)
    return [dict(identity='order-one antiholomorphic derivative', passed=True),
            dict(identity='even normalized C4 kernel', integral=1),
            dict(identity='bounded-projection pair normalization', passed=True),
            dict(identity='relative powers at delta=q^(-1/16)', boundary=str(-1+5*exponent),
                 infrared=str(-exponent), discrete=str(-1+2*exponent), continuum=str(-2*exponent))]


def energy_representation_checks():
    rows = []
    T, lengths = 1.3, np.array([.7, 1.4, 12.8])
    for n in [12, 24]:
        K, M, C, CI, KT, MT = assemble_circles(T, lengths, n)
        ev, U = eigh(M)
        sqrtM = (U*np.sqrt(ev))@U.T
        invsqrtM = (U/np.sqrt(ev))@U.T
        gradients = np.zeros_like(K)
        throat_edges = []
        for j, length in enumerate(lengths):
            for edge in range(2*n):
                index = j*2*n+edge
                h = (T if edge < n else length)/n
                gradients[index, index] = -1/np.sqrt(h)
                gradients[index, j*2*n+(edge+1) % (2*n)] = 1/np.sqrt(h)
                throat_edges.append(edge < n)
        require(np.linalg.norm(gradients.T@gradients-K) < 1e-10, 'independent incidence stiffness')
        QI, QC = null_space(CI), null_space(C)
        reflect_eig, reflect_basis = eigh(QI.T@K@QI, QI.T@M@QI)
        keep = (reflect_eig > 1e-7) & (reflect_eig < 6.6**2)
        wr, ur = np.sqrt(reflect_eig[keep]), QI@reflect_basis[:, keep]
        reflecting = float(np.sum(pair_weights(wr, ur.T@MT@ur, ur.T@KT@ur)*band_weight(wr[:, None]+wr[None, :])))
        for mass in [0., 4., 16.]:
            potential = mass**2*np.kron(np.ones((3, 3))/3, MT[:2*n, :2*n])
            pe, pu = eigh(potential)
            retained = pe > 1e-9
            mass_factor = (pu[:, retained]*np.sqrt(pe[retained])).T
            gradient = np.vstack((gradients, mass_factor))
            eig, basis = eigh(K+potential, M)
            keep = (eig > 1e-7) & (eig < 6.6**2)
            ws, u = np.sqrt(eig[keep]), basis[:, keep]
            plus = np.vstack((sqrtM@u, (gradient@u)/ws))/np.sqrt(2)
            minus = np.vstack((sqrtM@u, -(gradient@u)/ws))/np.sqrt(2)
            xscalar = invsqrtM@MT@invsqrtM
            xgradient = np.concatenate((np.array(throat_edges, float), np.ones(len(mass_factor))))
            overlap = minus[:len(K)].T@xscalar@plus[:len(K)]
            overlap += minus[len(K):].T@(xgradient[:, None]*plus[len(K):])
            A, B = u.T@MT@u, u.T@(KT+potential)@u
            predicted = (A-B/np.outer(ws, ws))/2
            pair_error = float(np.max(abs(overlap-predicted)))
            require(pair_error < 3e-10, 'energy representation matrix elements')
            d = gradient@invsqrtM
            residual = np.vstack((d.T@plus[len(K):], d@plus[:len(K)]))-plus*ws
            residual_norm = float(np.linalg.norm(residual)/max(1., np.linalg.norm(plus*ws)))
            require(residual_norm < 1e-9, 'first-order eigenvectors')
            W = pair_weights(ws, A, B)
            output = ws[:, None]+ws[None, :]
            infrared_cases = []
            for epsilon in [.75, 1.5]:
                omitted = float(np.sum(W*band_weight(output)*((ws[:, None] <= epsilon) | (ws[None, :] <= epsilon))))
                bound = 6.6*epsilon*sum(ws <= epsilon)  # max of the unnormalized previous band weight is 1.
                require(omitted <= bound+1e-9, 'Bessel soft-pair bound')
                infrared_cases.append(dict(epsilon=epsilon, omitted=omitted, bound=float(bound)))
            ce, cb = eigh(QC.T@(K+potential)@QC, QC.T@M@QC)
            ck = (ce > 1e-7) & (ce < 6.6**2)
            wc, uc = np.sqrt(ce[ck]), QC@cb[:, ck]
            cut_value = float(np.sum(pair_weights(wc, uc.T@MT@uc, uc.T@(KT+potential)@uc)*band_weight(wc[:, None]+wc[None, :])))
            cut_error = abs(cut_value-reflecting)
            require(cut_error < 3e-9, 'isolated massive singlet has no throat-energy pair response')
            rows.append(dict(cells_per_edge=n, mass=mass, pair_matrix_error=pair_error,
                first_order_residual=residual_norm, cut_response=cut_value, reflecting_response=reflecting,
                cut_response_error=cut_error, infrared_cases=infrared_cases))
    return rows


def chi_and_derivative(y):
    r = np.clip(y-1, 0., 1.)
    polynomial = Poly([0, 0, 0, 10, -15, 6])
    return np.where(y <= 1, 1., 1-polynomial(r)), np.where(y <= 1, 0., -polynomial.deriv()(r))


def cauchy_scalar(lam, center, width, n):
    nodes, weights = leggauss(n)
    split_x = sorted(set([center-width, center+width]+([lam] if center-width < lam < center+width else [])))
    total = 0.
    for left, right in zip(split_x[:-1], split_x[1:]):
        xs, wx = (left+right)/2+(right-left)*nodes/2, (right-left)*weights/2
        fx, fp, fpp = (bump((xs-center)/width, j)/width**j for j in range(3))
        for low, high in zip([0., .01, .1, 1.], [.01, .1, 1., 2.]):
            ys, wy = (low+high)/2+(high-low)*nodes/2, (high-low)*weights/2
            chi, chip = chi_and_derivative(ys)
            antiholo = (.5j*(ys*chi)[:, None]*fpp[None, :]
                        +.5j*chip[:, None]*fx[None, :]-.5*(ys*chip)[:, None]*fp[None, :])
            values = 2*np.real(antiholo/(lam-xs[None, :]-1j*ys[:, None]))/PI
            total += float(wy@values@wx)
    return total


def cauchy_checks():
    rng = np.random.default_rng(71812)
    A = np.diag([-2.4, -1.6, -.4, .7, 1.8, 2.3])
    vectors = rng.normal(size=(6, 2))*.15
    B = A+vectors@vectors.T
    raw = rng.normal(size=(6, 6))
    Q, _ = np.linalg.qr(raw)
    X = Q[:, :3]@Q[:, :3].T
    comparisons = []
    for z in [.1+.2j, 1.8+.03j, -1.2+.7j]:
        diff = np.linalg.inv(A-z*np.eye(6))-np.linalg.inv(B-z*np.eye(6))
        values = np.linalg.svd(diff, compute_uv=False)
        require(np.linalg.norm(values[2:]) < 2e-12, 'rank-two resolvent algebra')
        require(sum(values) <= 4/abs(z.imag)+1e-10, 'rank trace bound')
        comparisons.append(dict(z=[z.real, z.imag], rank_tail=float(np.linalg.norm(values[2:])), trace_norm=float(sum(values)), bound=4/abs(z.imag)))
    rows = []
    for H in [A, B]:
        ev, basis = eigh(H)
        exact_f, exact_g = bump((ev+2)/1.3), bump((ev-2)/1.3)
        exact_trace = float(np.sum(exact_f[:, None]*exact_g[None, :]*abs(basis.T@X@basis)**2))
        refinements = []
        for n in [64, 128]:
            f = np.array([cauchy_scalar(lam, -2., 1.3, n) for lam in ev])
            g = np.array([cauchy_scalar(lam, 2., 1.3, n) for lam in ev])
            F, G = (basis*f)@basis.T, (basis*g)@basis.T
            integral_trace = float(np.trace(F@X@G@X))
            error = float(max(np.max(abs(f-exact_f)), np.max(abs(g-exact_g))))
            refinements.append(dict(nodes_per_piece=n, scalar_error=error, pair_trace=integral_trace,
                                    pair_trace_error=abs(integral_trace-exact_trace)))
        require(refinements[-1]['scalar_error'] < 2e-6, 'Cauchy integral against exact functional calculus')
        require(refinements[-1]['pair_trace_error'] < 2e-6, 'two-frequency trace formula')
        rows.append(dict(exact_trace=exact_trace, refinements=refinements))
    return comparisons, rows


def derivative_norms(h, n):
    nodes, weights = leggauss(n)
    norms = np.zeros((3, 3))
    for t, wt in zip(nodes, weights):
        omega = 1+h*t
        breaks = sorted(set([0., omega]+[p for p in [h/2, h, omega-h, omega-h/2] if 0 < p < omega]))
        for left, right in zip(breaks[:-1], breaks[1:]):
            r = (left+right)/2+(right-left)*nodes/2
            s = omega-r
            integration = wt*h*(right-left)*weights/2
            def gd(x):
                ratio = x/h
                return [x*cutoff(ratio), cutoff(ratio)+ratio*cutoff(ratio, 1),
                        (2*cutoff(ratio, 1)+ratio*cutoff(ratio, 2))/h]
            gr, gs = gd(r), gd(s)
            for i in range(3):
                for j in range(3):
                    value = np.zeros_like(r)
                    for p in range(i+1):
                        for v in range(j+1):
                            k = i+j-p-v
                            value += .5*math.comb(i, p)*math.comb(j, v)*gr[p]*gs[v]*bump(t, k)/h**(k+1)
                    norms[i, j] += np.dot(integration, abs(value))
    return norms


def derivative_checks():
    rows = []
    # Polynomial coefficient sums provide genuine (very conservative) supremum bounds.
    p1, p2 = sum(abs(BP.deriv(1).coef)), sum(abs(BP.deriv(2).coef))
    GB = [1.25, 1+2*p1, 4*p1+4*p2]
    WB = [sum(abs(WP.deriv(k).coef)) for k in range(5)]
    constant = 0.
    for i in range(3):
        for j in range(3):
            constant += 1.25*sum(math.comb(i, p)*math.comb(j, v)*GB[p]*GB[v]*WB[i+j-p-v]
                                 for p in range(i+1) for v in range(j+1))
    for h in [.2, .1, .05, .025]:
        first, second = derivative_norms(h, 40), derivative_norms(h, 80)
        change = float(np.linalg.norm(second-first)/np.linalg.norm(second))
        scaled = float(h**4*np.sum(second))
        require(change < .004, 'cutoff derivative quadrature refinement')
        require(scaled <= constant, 'explicit h^-4 bound')
        rows.append(dict(h=h, derivative_l1_norms=second.tolist(), relative_refinement_change=change,
                         h4_scaled_sum=scaled, conservative_uniform_constant=constant))
    return rows


def soft_count_checks():
    rows = []
    for q in [100, 10000, 1000000]:
        for delta, epsilon in [(.2, .2), (.1, .1), (.04, .04)]:
            u = (np.arange(q)+.5)/q
            counts = 2*np.floor(epsilon*(PI+KAPPA*delta/u)/(2*PI)).astype(int)
            observed = int(sum(counts))
            J = min(q, max(0, int(np.floor(q*KAPPA*delta*epsilon/(2*PI-epsilon*PI)+.5))))
            bound = 0. if J == 0 else epsilon/PI*(PI*J+KAPPA*delta*q*(2+np.log(max(1, 2*J-1))))
            require(observed <= bound+1e-8, 'soft-mode harmonic bound')
            rows.append(dict(q=q, delta=delta, epsilon=epsilon, count=observed, eligible_channels=J, bound=float(bound)))
    return rows


def continuum_density(omega, nmax=30000):
    n = np.arange(2, nmax+1, dtype=float)
    value = KAPPA*omega**2*np.sin(omega*PI/2)**2/3*np.sum((1-1/n**2)/(2*PI*n-PI*omega)**2)
    # Uniform, conservative positive tail on Omega in [.75,1.25].
    bound = KAPPA*1.25**2/(3*PI**2*nmax)
    return float(value), float(bound)


def free_discrete(q, delta):
    h = delta
    maximum = int(np.floor((1+h)*(PI+2*KAPPA*delta*q)/(2*PI)))
    value, terms = 0., 0
    for n in range(2, maximum+1):
        ua = KAPPA*delta*(1-h)/(2*PI*n-PI*(1-h))
        ub = KAPPA*delta*(1+h)/(2*PI*n-PI*(1+h))
        first = max(0, int(np.ceil(q*ua-.5)))
        last = min(q-1, int(np.floor(q*ub-.5)))
        if first > last:
            continue
        u = (np.arange(first, last+1)+.5)/q
        omega = 2*PI*n/(PI+KAPPA*delta/u)
        weight = (n*n-1)/(6*PI**2*n**3)*omega**2*np.sin(omega*PI/2)**2
        value += float(np.sum(weight*bump((omega-1)/h)/h))
        terms += len(u)
    return PI*value/(q*delta), terms, maximum


def continuum_checks():
    rows = []
    for q, delta in [(20000, .08), (200000, .04), (2000000, .02)]:
        nodes, weights = leggauss(72)
        vals = [continuum_density(1+delta*t) for t in nodes]
        continuum = float(np.dot(weights*bump(nodes), [v[0] for v in vals]))
        tail = max(v[1] for v in vals)
        discrete, terms, maximum = free_discrete(q, delta)
        require(abs(discrete-continuum) < .015*CSTAR+tail, 'finite-channel continuum comparison')
        require(abs(continuum-CSTAR) < .04*CSTAR+tail, 'small-band coefficient comparison')
        rows.append(dict(q=q, delta=delta, q_delta_squared=q*delta**2, q_delta_fifth=q*delta**5,
            discrete_normalized_response=discrete, continuum_normalized_response=continuum,
            exact_limit=CSTAR, continuum_tail_bound=tail, contributing_channel_harmonics=terms, maximum_harmonic=maximum,
            scope='free discrete-to-continuum check, not a small interacting-error parameter point',
            relative_discrete_difference=abs(discrete-continuum)/CSTAR,
            relative_continuum_difference=abs(continuum-CSTAR)/CSTAR))
    return rows


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    exact = exact_checks()
    energy = energy_representation_checks()
    resolvents, cauchy = cauchy_checks()
    derivatives = derivative_checks()
    soft = soft_count_checks()
    continuum = continuum_checks()
    result = dict(status='passed_singlet_stress_continuum_checks',
        scope='supplementary checks for the analytic Gaussian-model absorption limit',
        script_sha256=sha(Path(__file__)), note_sha256=sha(ROOT/'notes/singlet-stress-continuum.md'),
        dependency_sha256=sha(ROOT/'checks/schwinger_return_feedback.py'),
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, sympy=sp.__version__,
        exact_checks=exact, energy_representation_cases=energy, resolvent_cases=resolvents,
        cauchy_cases=cauchy, derivative_cases=derivatives, soft_count_cases=soft, continuum_cases=continuum,
        model_uniform_absorption_argument=True, physical_damping_proved=False,
        renormalized_real_response_controlled=False, physical_model_closed=False,
        novelty_established=False, lean_invocations=0)
    (ROOT/'evidence/singlet-stress-continuum.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], counts={k:len(v) for k,v in result.items() if isinstance(v,list)},
                         coefficient=CSTAR, continuum_cases=continuum,
                         maximum_pair_error=max(x['pair_matrix_error'] for x in energy),
                         maximum_cut_error=max(x['cut_response_error'] for x in energy)), indent=2))


if __name__ == '__main__':
    main()
