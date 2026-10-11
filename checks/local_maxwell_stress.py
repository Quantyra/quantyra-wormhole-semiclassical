"""Direct tensor and supplementary radial checks of local Maxwell stress.

This does not solve the joined Einstein-Maxwell/quantum system.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_bvp, solve_ivp
from scipy.special import roots_legendre
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def text_hash(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').replace('\r\n', '\n').encode()).hexdigest()


def tensor_checks():
    t, x, c, az = sp.symbols('t x c az', real=True)
    coords = [t, x, c, az]
    N, R = sp.Function('N')(x), sp.Function('R')(x)
    diagonal = [-N**2, N**2, R**2/(1-c*c), R**2*(1-c*c)]
    metric = sp.diag(*diagonal)
    inverse = sp.diag(*[1/v for v in diagonal])
    gam = {}
    for a in range(4):
        for b in range(4):
            for d in range(4):
                value = (sp.diff(metric[a, d], coords[b])+sp.diff(metric[a, b], coords[d])
                         -sp.diff(metric[b, d], coords[a]))/(2*diagonal[a])
                if value != 0:
                    gam[a, b, d] = sp.simplify(value)
    Gm = lambda a, b, d: gam.get((a, b, d), 0)
    ricci = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            value = sum(sp.diff(Gm(k, a, b), coords[k])-sp.diff(Gm(k, a, k), coords[b])
                        +sum(Gm(k, a, b)*Gm(j, k, j)-Gm(j, a, k)*Gm(k, b, j)
                             for j in range(4)) for k in range(4))
            ricci[a, b] = sp.simplify(value)
    scalar = sp.simplify(sum(ricci[k, k]/diagonal[k] for k in range(4)))
    einstein = ricci-metric*scalar/2
    rho, pressure = einstein[0, 0]/N**2, einstein[1, 1]/N**2
    assert sp.simplify(sp.diff(R**2, x, 2)-2*N**2-N**2*R**2*(pressure-rho)) == 0
    assert sp.simplify(sp.diff(R, x, 2)-2*sp.diff(N, x)*sp.diff(R, x)/N
                       +N**2*R*(rho+pressure)/2) == 0
    q, coupling = sp.symbols('q coupling', positive=True)
    b = sp.Function('b')(x, c)
    a = (1-c*c)*b
    F = sp.zeros(4)
    F[1, 3], F[3, 1] = sp.diff(a, x), -sp.diff(a, x)
    F[2, 3], F[3, 2] = -q/2+sp.diff(a, c), q/2-sp.diff(a, c)
    raised = inverse*F*inverse
    Fsq = sum(F[i, j]*raised[i, j] for i in range(4) for j in range(4))
    stress = (F*inverse*F.T-metric*Fsq/4)/coupling**2
    br = (q/2-sp.diff(a, c))/R**2
    bt2 = (1-c*c)*sp.diff(b, x)**2/(N**2*R**2)
    assert sp.simplify(stress[0, 0]/N**2-(br**2+bt2)/(2*coupling**2)) == 0
    assert sp.simplify(stress[1, 1]/N**2-(-br**2+bt2)/(2*coupling**2)) == 0
    assert sp.simplify((stress[0, 0]+stress[1, 1])/N**2-bt2/coupling**2) == 0
    assert sp.simplify(stress[1, 2]+br*sp.diff(b, x)/coupling**2) == 0
    averages = []
    for ell in range(1, 9):
        Y = (1-c*c)*sp.diff(sp.legendre(ell, c), c)
        assert sp.integrate(sp.diff(Y, c), (c, -1, 1)) == 0
        averages.append(ell)
    ell, mass = sp.symbols('ell mass', positive=True)
    L = ell*(ell+1)
    scalar_mass = sp.Matrix([[L+2, -4*L], [-1, L]])
    vector_mass = sp.Matrix([[L+2, 2], [2*L, L]])
    S = sp.diag(1, 2*sp.sqrt(L))
    symmetric = sp.simplify(S*scalar_mass*S.inv())
    assert sp.simplify(symmetric-symmetric.T) == sp.zeros(2)
    characteristic = (mass-ell*(ell-1))*(mass-(ell+1)*(ell+2))
    for M in [scalar_mass, vector_mass]:
        assert sp.simplify((mass*sp.eye(2)-M).det()-characteristic) == 0
    box, trace, pi, em = sp.symbols('box trace pi em')
    eq1 = trace/2-(box-4)*pi/4-2*em
    eq2 = trace-pi-(box-2)*em
    assert sp.expand(4*eq1-2*eq2+(box-6)*(pi-2*em)) == 0
    radius, re, rr, alpha = sp.symbols('radius re rr alpha', positive=True)
    fn = (1-re/radius)**2
    dx = lambda value: sp.simplify(fn*sp.diff(value, radius))
    assert sp.simplify(dx(dx(radius**2))-2*fn+2*fn*re**2/radius**2) == 0
    rnmode = radius**2-re**2
    assert sp.simplify(sp.diff(fn*sp.diff(rnmode, radius), radius)-2*rnmode/radius**2) == 0
    f = 1+rr**2
    b_even = 1+rr*sp.atan(rr)
    assert sp.simplify(sp.diff(f*sp.diff(b_even, rr), rr)-2*b_even) == 0
    assert sp.simplify(sp.diff(f*sp.diff(rr, rr), rr)-2*rr) == 0
    I2 = (sp.atan(rr)+rr/f)/2
    assert sp.simplify(sp.diff(I2, rr)-1/f**2) == 0
    cc = 1+rr*sp.atan(rr)
    ss = -rr**2-rr*(rr**2+3)*sp.atan(rr)+sp.log(f)
    C, SS = 1+alpha*cc, f+alpha*ss
    condition = sp.diff(SS, rr)*sp.diff(C, rr)+SS*sp.diff(C, rr, 2)
    first = sp.expand(condition).coeff(alpha, 1)
    assert sp.simplify(first-2*(1+rr*sp.atan(rr))) == 0
    regularity = []
    numeric = sp.lambdify((rr, alpha), [C, SS/f, condition], 'numpy')
    for av in [1e-2, 1e-4, 1e-8]:
        grid = np.r_[0., np.geomspace(1e-6, av**-.5, 2000)]
        cv, sv, iv = numeric(grid, av)
        assert np.min(cv) > 0 and np.min(sv) > .8 and np.max(iv) < 2
        regularity.append(dict(alpha=av, overlap_endpoint=av**-.5,
                               max_area_comparison_lhs=float(np.max(iv)),
                               minimum_S_over_f=float(np.min(sv))))
    return dict(direct_einstein_identities=2, direct_maxwell_identities=4,
                zero_flux_angular_modes=averages, mass_matrices_checked=2,
                physical_l1_mass_squared=6, reference_overlap_samples=regularity)


def radial_checks():
    cases = []
    for hmin in [1e-4, 1e-8]:
        for hmax in [.01, 1., 100.]:
            lo, hi = np.log(hmin), np.log(hmax)
            normalization = hmax*(2+hmax)
            exact = lambda u: np.exp(u)*(2+np.exp(u))/normalization
            exact_derivative = lambda u: 2*np.exp(u)*(1+np.exp(u))/normalization
            def ode(u, y):
                h = np.exp(u)
                return np.array([y[1], 2*y[0]-(1-h)/(1+h)*y[1]])
            def bc(a, b):
                return np.array([a[0]-exact(lo), b[0]-1])
            solutions = []
            grid = np.linspace(lo, hi, 1201)
            for count, tolerance in [(120, 1e-8), (240, 1e-10)]:
                mesh = np.linspace(lo, hi, count)
                guess = np.array([np.linspace(exact(lo), 1, count), np.ones(count)/(hi-lo)])
                sol = solve_bvp(ode, bc, mesh, guess, tol=tolerance, max_nodes=16000)
                assert sol.success, (hmin, hmax, sol.message)
                solutions.append(sol.sol(grid))
            h = np.exp(grid)
            values = solutions[-1]
            er = np.max(abs(values[0]-exact(grid)))
            ed = np.max(abs(values[1]-exact_derivative(grid)))
            assert er < 2e-9 and ed < 2e-9
            # Physical field coefficients normalized by B_ext, independently
            # reconstructed from the solved potential and derivative.
            radial = normalization*values[0]/(1+h)**2
            transverse = normalization*values[1]/(2*(1+h)**2)
            expected_r = h*(2+h)/(1+h)**2
            expected_t = h/(1+h)
            ferr = float(max(np.max(abs(radial-expected_r)), np.max(abs(transverse-expected_t))))
            assert ferr < 1e-6
            assert max(np.max(abs(radial)), np.max(abs(transverse))) <= 1+1e-7
            cases.append(dict(hmin=hmin, hmax=hmax, potential_error=float(er),
                              derivative_error=float(ed), physical_field_error=ferr,
                              solution_refinement=float(np.max(abs(solutions[0]-values)))))
    return cases


def focusing_checks():
    cases = []
    for end in [.1, 1., 10., 1000.]:
        f_end = 1+end*end
        base = 1+end*np.arctan(end)
        deriv = np.arctan(end)+end/f_end
        amplitude = 1/max(2*base, np.sqrt(f_end)*deriv)
        def mean_null(rho):
            f = 1+rho*rho
            ap = amplitude*(np.arctan(rho)+rho/f)
            return (2/3)*f*ap*ap
        integrals = []
        for order in [64, 128]:
            nodes, weights = roots_legendre(order)
            u = (nodes+1)*np.arcsinh(end)/2
            rho = np.sinh(u)
            # Integrate the actual angular field, not the analytic average.
            ca, wa = roots_legendre(order//2)
            f = 1+rho*rho
            ap = amplitude*(np.arctan(rho)+rho/f)
            transverse = -np.sqrt(f)[:, None]*ap[:, None]*np.sqrt(1-ca*ca)[None, :]
            averaged = transverse**2@wa/2
            integral = np.dot(weights, averaged*np.cosh(u)/f)*np.arcsinh(end)/2
            integrals.append(float(integral))
        # Independent evolution of the focusing equation, kappa=re=g=1.
        def ode(rho, y):
            return np.array([y[1], -mean_null(rho)/(2*(1+rho*rho))])
        sol = solve_ivp(ode, (0, end), [0., 0.], method='DOP853', rtol=1e-10, atol=1e-12)
        assert sol.success
        slope = float(sol.y[1, -1])
        discrepancy = abs(slope+integrals[-1]/2)
        assert discrepancy < 1e-9
        I1 = np.arctan(end)
        I2 = (I1+end/(1+end*end))/2
        assert 1 <= I1/I2 <= 2
        assert 0 <= integrals[-1]/I2 <= 2
        assert slope <= 0
        assert abs(integrals[-1]-integrals[0]) < 1e-10
        cases.append(dict(overlap_endpoint=end, exact_field_amplitude=amplitude,
                          normalized_maxwell_integral=integrals[-1], focusing_slope=slope,
                          bound_ratio=integrals[-1]/(2*I2), weight_ratio=float(I1/I2),
                          ode_quadrature_error=discrepancy,
                          quadrature_refinement=abs(integrals[-1]-integrals[0])))
    return cases


def verify():
    exact = tensor_checks()
    radial = radial_checks()
    focusing = focusing_checks()
    # q/g^2=q^3 and re/d=(16 pi delta q^3)^-1, delta=log(q)^-2.
    powers = dict(component_q_power=3-2*3, component_log_power=2*2,
                  null_q_power=3-4*3, null_log_power=2*4)
    assert powers['component_q_power'] == -3 and powers['null_q_power'] == -9
    inputs = ['notes/local-maxwell-stress.md', 'notes/maxwell-dirac-comparison.md',
              'notes/coupled-maxwell-matching.md', 'notes/casimir-reference-matching.md']
    return dict(status='passed_local_maxwell_stress_checks', exact=exact,
                rn_radial_cases=radial, focusing_cases=focusing, scale_powers=powers,
                local_field_bound='C q/d^2 under inner area and local regularity conditions',
                direct_null_stress_nonnegative=True, monopole_cross_term_cancels=True,
                full_einstein_response_solved=False, full_second_order_monopole_response=False,
                actual_quantum_state_derived=False, bulk_mass_gap_closes_boundary_modes=False,
                novelty_established=False, lean_invocations=0, supplementary_numerics=True,
                versions=dict(numpy=np.__version__, scipy=scipy.__version__, sympy=sp.__version__),
                input_sha256={p: text_hash(ROOT/p) for p in inputs})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    receipt = verify()
    receipt['script_sha256'] = text_hash(Path(__file__))
    Path(args.output).write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'], rn_cases=len(receipt['rn_radial_cases']),
                         focusing_cases=len(receipt['focusing_cases']))))
