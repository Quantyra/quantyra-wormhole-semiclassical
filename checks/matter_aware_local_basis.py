"""First-order local coefficient map with charged stress and current retained.

The map is a bulk EFT field redefinition, not a complete renormalization
prescription. Angular contact tests do not bound a continuum many-body vacuum.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import math
import platform

import mpmath as mp
import numpy as np
import sympy as s


def zero(value):
    value = s.expand(value)
    assert value == 0, str(value)


def symbolic_map():
    A, B, C, D, kap2, cosm = s.symbols('A B C D kappa2 Lambda')
    tau, t2, st, I, I4 = s.symbols('tau T2 ST I I4')
    R0 = 4*cosm-kap2*tau
    Ric0sq = 4*cosm**2-2*cosm*kap2*tau+kap2**2*(I4-I**2/4+2*st+t2)
    Ric0X = cosm*I+kap2*(I4-I**2/4+st-tau*I/4)
    direct = A*R0**2+B*Ric0sq+kap2*C*R0*I+kap2*D*Ric0X
    reduced = (16*A+4*B)*cosm**2-(8*A+2*B)*cosm*kap2*tau
    reduced += (4*C+D)*cosm*kap2*I
    reduced += kap2**2*((B+D)*(I4-I**2/4)+A*tau**2+B*t2+(2*B+D)*st-(C+D/4)*tau*I)
    zero(direct-reduced)
    # Off-shell residual with E=tr(E_mu_nu), E2=E_mu_nu E^mu_nu.
    E, E2, R0E, EX = s.symbols('E E2 Ric0E EX')
    residual = A*(-2*R0*E+E**2)+B*(2*R0E-R0*E+E2)
    residual += kap2*(-C*I*E+D*(EX-I*E/2))
    contraction = A*(-2*R0+E)*E+B*(2*R0E-R0*E+E2)
    contraction += kap2*(D*EX-(C+D/2)*I*E)
    zero(residual-contraction)
    c = s.symbols('c1:10')
    aa, bb, cc, dd = c[0]-c[2], c[1]+4*c[2], c[3]-c[5]/3, c[4]+2*c[5]
    i2, i4 = c[6]-(bb+dd)/4, c[7]+bb+dd
    # Independent source regression: Cheung--Liu--Remmen eq. 95, D=4.
    r1, r2, r3, r4 = s.symbols('r1:5')
    shifts = [-r1/4-r2/2, r1/2, 0, r1/8-r3/4-r4/2,
              -r1/2+r3/2, 0, r3/8, -r3/2, 0]
    repl = {x: x+dx for x, dx in zip(c, shifts)}
    for invariant in [c[2], c[5], i2, i4]:
        zero(invariant.subs(repl, simultaneous=True)-invariant)
    j, div = s.symbols('j divF')
    zero(div**2-j**2-(div-j)*(div+j))
    # Regression to the preceding isolated J^2 direction.
    cj = s.symbols('cJ')
    jj = {x: 0 for x in c}
    jj.update({c[6]: -cj/8, c[7]: cj/4})
    zero(i2.subs(jj)+cj/8)
    zero(i4.subs(jj)-cj/4)
    return {'A': str(aa), 'B': str(bb), 'U': str(cc), 'V': str(dd),
            'I_squared_coefficient': str(i2), 'I4_coefficient': str(i4),
            'matter_contacts': 'kappa^4[A*tau^2+B*T2+(2B+V)*S.T-(U+V/4)*tau*I]+kappa^2*c9*j^2',
            'off_shell_residual': 'E_mu_nu H^mu_nu+kappa^2*c9*(divF-j).(divF+j)',
            'source_free_field_redefinition_invariants': 'c3, c6, c7-(B+V)/4, c8+B+V',
            'order': 'linear in the local coefficients; boundary/source transformations required'}


def tensor_map_checks():
    # Direct algebraic curvature contractions, independent of the scalar map.
    # Sum of simple bivector squares obeys all algebraic Riemann identities.
    rng = np.random.default_rng(20261013)
    eta = [-1, 1, 1, 1]
    metric = s.diag(*eta)
    dot = lambda U, V: sum(eta[i]*eta[j]*U[i, j]*V[i, j] for i in range(4) for j in range(4))
    trace = lambda U: sum(eta[i]*U[i, i] for i in range(4))
    rows = []
    for case in range(8):
        curv = s.MutableDenseNDimArray.zeros(4, 4, 4, 4)
        for n in range(7):
            u, v = [list(map(int, rng.integers(-3, 4, size=4))) for _ in range(2)]
            form = s.Matrix(4, 4, lambda i, j: u[i]*v[j]-u[j]*v[i])
            weight = s.Rational((-1)**n, n+1)
            for i, j, k, l in itertools.product(range(4), repeat=4):
                curv[i, j, k, l] += weight*form[i, j]*form[k, l]
        Ric = s.Matrix(4, 4, lambda j, l: sum(eta[i]*curv[i, j, i, l] for i in range(4)))
        R = trace(Ric)
        Riem2 = sum(eta[i]*eta[j]*eta[k]*eta[l]*curv[i, j, k, l]**2 for i, j, k, l in itertools.product(range(4), repeat=4))
        Euler = Riem2-4*dot(Ric, Ric)+R**2
        fraw = s.Matrix(rng.integers(-3, 4, size=(4, 4)).tolist())
        F = fraw-fraw.T
        X = F*metric*F.T
        I, I4 = trace(X), s.trace((F*metric)**4)
        zero(dot(X, X)-I4)
        S = X-metric*I/4
        rawT = s.Matrix(rng.integers(-4, 5, size=(4, 4)).tolist())
        T = (rawT+rawT.T)/2
        tau = trace(T)
        kap2, cosm = s.Rational(case+1, 7), s.Rational(2-case, 11)
        Ric0 = cosm*metric+kap2*(S+T-metric*tau/2)
        R0 = trace(Ric0)
        Eten = Ric-metric*R/2+cosm*metric-kap2*(S+T)
        E = trace(Eten)
        zero(R-R0+E)
        for value in Ric-Ric0-Eten+metric*E/2:
            zero(value)
        RF, WF = 0, 0
        for i, j, k, l in itertools.product(range(4), repeat=4):
            # Riem = Weyl + (g wedge Ric)/2 - R(g wedge g)/6.
            ricpart = (metric[i, k]*Ric[j, l]-metric[i, l]*Ric[j, k]
                       -metric[j, k]*Ric[i, l]+metric[j, l]*Ric[i, k])/2
            scalarpart = -R*(metric[i, k]*metric[j, l]-metric[i, l]*metric[j, k])/6
            ff = eta[i]*eta[j]*eta[k]*eta[l]*F[i, j]*F[k, l]
            RF += curv[i, j, k, l]*ff
            WF += (curv[i, j, k, l]-ricpart-scalarpart)*ff
        zero(RF-WF-2*dot(Ric, X)+R*I/3)
        div = s.Matrix(list(map(int, rng.integers(-3, 4, size=4))))
        current = s.Matrix(list(map(int, rng.integers(-3, 4, size=4))))
        square = lambda v: (v.T*metric*v)[0]
        # Every coefficient direction separately, plus a mixed direction.
        for choice in range(10):
            coeff = [s.Integer(int(i == choice)) for i in range(9)] if choice < 9 else [s.Rational(i-4, i+1) for i in range(9)]
            c1, c2, c3, c4, c5, c6, c7, c8, c9 = coeff
            A, B, C, D = c1-c3, c2+4*c3, c4-c6/3, c5+2*c6
            original = c1*R**2+c2*dot(Ric, Ric)+c3*Riem2
            original += kap2*(c4*R*I+c5*dot(Ric, X)+c6*RF+c9*square(div))
            original += kap2**2*(c7*I**2+c8*I4)
            reduced = c3*Euler+kap2*c6*WF+(16*A+4*B)*cosm**2
            reduced += -(8*A+2*B)*cosm*kap2*tau+(4*C+D)*cosm*kap2*I
            reduced += kap2**2*((c7-(B+D)/4)*I**2+(c8+B+D)*I4)
            reduced += kap2**2*(A*tau**2+B*dot(T, T)+(2*B+D)*dot(S, T)-(C+D/4)*tau*I)
            reduced += kap2*c9*square(current)
            H = A*(-2*R0+E)*metric+B*(2*Ric0-R0*metric+Eten)
            H += kap2*(D*X-(C+D/2)*I*metric)
            residual = dot(Eten, H)+kap2*c9*((div-current).T*metric*(div+current))[0]
            zero(original-reduced-residual)
        rows.append({'case': case, 'coefficient_directions': 10, 'exact_off_shell_residual': True})
    return rows


def beta(q, ell):
    return s.Rational(math.factorial(q-1)*math.factorial(q), math.factorial(q-1-ell)*math.factorial(q+ell))


def differential_basis_checks():
    r = s.symbols('r', positive=True)
    rows = []
    for powers in [(1, 0, 1, 1), (2, -1, 1, 0), (-1, 1, 2, 1)]:
        diag = [(-1 if i == 0 else 1)*r**(2*powers[i]) for i in range(4)]
        inv = [1/value for value in diag]
        metric = s.diag(*diag)
        diff = lambda value, index: s.diff(value, r) if index == 1 else s.S.Zero
        Gamma = s.MutableDenseNDimArray.zeros(4, 4, 4)
        for a, b, c in itertools.product(range(4), repeat=3):
            Gamma[a, b, c] = inv[a]*(diff(metric[a, c], b)+diff(metric[a, b], c)-diff(metric[b, c], a))/2
        Riem = s.MutableDenseNDimArray.zeros(4, 4, 4, 4)
        for a, b, c, d in itertools.product(range(4), repeat=4):
            Riem[a, b, c, d] = s.expand(diag[a]*(diff(Gamma[a, b, d], c)-diff(Gamma[a, b, c], d)
                +sum(Gamma[a, c, e]*Gamma[e, b, d]-Gamma[a, d, e]*Gamma[e, b, c] for e in range(4))))
        Ric = s.Matrix(4, 4, lambda b, d: sum(inv[a]*Riem[a, b, a, d] for a in range(4)))
        F = s.Matrix([[0, r*r, 0, 0], [-r*r, 0, r**3, r+1], [0, -r**3, 0, 2], [0, -r-1, -2, 0]])
        # This two-form is closed, including its nonzero constant 23 flux.
        for a, b, c in itertools.product(range(4), repeat=3):
            zero(diff(F[b, c], a)+diff(F[c, a], b)+diff(F[a, b], c))
        raised = s.diag(*inv)*F*s.diag(*inv)
        volume = r**sum(powers)
        div = s.Matrix([s.diff(volume*raised[1, b], r)/volume for b in range(4)])
        grad = s.MutableDenseNDimArray.zeros(4, 4, 4)
        for a, b, c in itertools.product(range(4), repeat=3):
            grad[a, b, c] = s.expand(diff(F[b, c], a)-sum(Gamma[e, a, b]*F[e, c]+Gamma[e, a, c]*F[b, e] for e in range(4)))
        grad2 = sum(inv[a]*inv[b]*inv[c]*grad[a, b, c]**2 for a, b, c in itertools.product(range(4), repeat=3))
        X = F*s.diag(*inv)*F.T
        RicX = sum(inv[a]*inv[b]*Ric[a, b]*X[a, b] for a, b in itertools.product(range(4), repeat=2))
        RFF = sum(Riem[a, b, c, d]*raised[a, b]*raised[c, d] for a, b, c, d in itertools.product(range(4), repeat=4))
        boundary_radial = inv[1]*sum(raised[b, c]*grad[1, b, c] for b, c in itertools.product(range(4), repeat=2))
        boundary_radial -= 2*sum(raised[1, b]*diag[b]*div[b] for b in range(4))
        divergence = s.diff(volume*boundary_radial, r)/volume
        zero(grad2-2*(div.T*metric*div)[0]+2*RicX-RFF-divergence)
        rows.append({'metric_powers': powers, 'closed_F_and_covariant_derivative_identity': 'exact for all r>0',
                     'boundary_divergence': str(s.expand(divergence))})
    return rows


def angular_exact_checks():
    z = s.symbols('z')
    rows = []
    for q in range(1, 17):
        bs = [beta(q, ell) for ell in range(q)]
        zero(sum((2*ell+1)*b for ell, b in enumerate(bs))-q)
        zero(sum((2*ell+1)*b*b for ell, b in enumerate(bs))-s.Rational(q*q, 2*q-1))
        polynomial = sum((2*ell+1)*bs[ell]*s.legendre(ell, z)/q for ell in range(q))
        zero(polynomial-((1+z)/2)**(q-1))
        rows.append({'q': q, 'linear_sum': q, 'square_sum': str(s.Rational(q*q, 2*q-1))})
    return rows


def angular_quadrature():
    rows = []
    for q in [1, 2, 3, 4, 6, 8]:
        for refinement in [1, 2]:
            x, w = np.polynomial.legendre.leggauss(refinement*(q+4))
            phi = 2*np.pi*np.arange(refinement*(4*q+4))/(refinement*(4*q+4))
            theta, phi = np.meshgrid(np.arccos(x), phi, indexing='ij')
            theta, phi = theta.ravel(), phi.ravel()
            weights = np.repeat(w, refinement*(4*q+4))*2*np.pi/(refinement*(4*q+4))
            eta = np.array([math.sqrt(q*math.comb(q-1, k)/(4*np.pi))*np.sin(theta/2)**k
                            *np.cos(theta/2)**(q-1-k)*np.exp(1j*((q-1)/2-k)*phi) for k in range(q)]).T
            gram = eta.conj().T@(weights[:, None]*eta)
            products = np.einsum('ni,nj->nij', eta, eta).reshape(len(weights), q*q)
            pair = products.conj().T@(weights[:, None]*products)
            conjugate = pair.reshape(q, q, q, q).transpose(0, 3, 2, 1).reshape(q*q, q*q)
            eval_pair = np.linalg.eigvalsh(pair)
            eval_conj = np.linalg.eigvalsh(conjugate)
            predicted_pair = np.array([0.]*(q*q-(2*q-1))+[q*q/(4*np.pi*(2*q-1))]*(2*q-1))
            predicted_conj = np.sort(np.concatenate([np.full(2*l+1, q*float(beta(q, l))/(4*np.pi)) for l in range(q)]))
            errors = {'gram': float(np.max(abs(gram-np.eye(q)))),
                      'pair_spectrum': float(np.max(abs(eval_pair-predicted_pair))),
                      'conjugate_spectrum': float(np.max(abs(eval_conj-predicted_conj)))}
            assert max(errors.values()) < 8e-12, (q, refinement, errors)
            rows.append({'q': q, 'refinement': refinement, 'errors': errors})
    return rows


def candidate_and_vacuum():
    mp.mp.dps = 85
    k = mp.mpf(2000)
    q, gauge = k**4, k**-4
    C = mp.sqrt(1+12*mp.pi*k**6)
    gamma = 6*C/(k**10*(C-1))
    alpha = 8*mp.pi*gamma*gauge**2/C  # alpha9 per unit c9, before physical matching
    frobenius2 = alpha**2*q**4/(4*mp.pi**4*(2*q-1))
    rows = []
    for N in [4, 16, 64, 256, 1024]:
        vacnorm = sum(s.Integer(n)**2 for n in range(1, N+1))
        zero(vacnorm-s.Rational(N*(N+1)*(2*N+1), 6))
        weighted = sum(s.Rational(n*n, 2*n+1) for n in range(1, N+1))
        assert weighted > s.Rational(N*N, 6)
        rows.append({'N': N, 'vacuum_column_norm_over_G2': str(vacnorm),
                     'weighted_norm_over_G2_N2': float(weighted/(N*N))})
    return {'candidate_k2000_per_unit_c9': {
                'alpha9': str(alpha), 'conjugate_pair_coefficient_norm': str(abs(alpha)*q/(4*mp.pi)),
                'maximum_current_coupling': str(abs(alpha)*q/(2*mp.pi**2)),
                'current_coupling_Frobenius_norm': str(mp.sqrt(frobenius2)),
                'physical_coefficient_assumed': False},
            'vacuum_cases': rows,
            'continuum_decision': 'For nonzero G, weighted vacuum column grows as ||G||_F^2 N^2/4; no cutoff-uniform free-vacuum form estimate follows from the finite angular norm.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print('Checking full displayed coefficient map and off-shell tensor residuals', flush=True)
    result = {'map': symbolic_map(), 'tensor_cases': tensor_map_checks()}
    result['differential_basis_cases'] = differential_basis_checks()
    print('Checking exact angular identities and independent LLL quadrature spectra', flush=True)
    result['angular_exact'] = angular_exact_checks()
    result['angular_quadrature'] = angular_quadrature()
    result['candidate_and_vacuum'] = candidate_and_vacuum()
    result.update({'status': 'PASS', 'scope': 'displayed parity-even bosonic four-derivative span with induced charged-matter contacts; not a complete Dirac/boundary basis or renormalized opening bound',
                   'versions': {'python': platform.python_version(), 'sympy': s.__version__, 'numpy': np.__version__, 'mpmath': mp.__version__},
                   'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('PASS: matter-aware local map, angular spectra and contact-estimate limitation', flush=True)


if __name__ == '__main__':
    main()
