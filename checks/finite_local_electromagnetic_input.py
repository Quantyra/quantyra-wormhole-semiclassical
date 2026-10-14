"""S074 finite local input invisible to the magnetic reference background.

The selected J^2 operator, J=F*F, is an EFT coefficient, not a prediction
of the previous nonlocal bubble. All-orders rational formulas refer only
to this specified local truncation; an EFT use retains its stated order.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform

import mpmath as mp
import numpy as np
import sympy as s


def zero(expr):
    value = s.factor(s.expand(expr))
    assert value == 0, str(value)


def zeros(matrix):
    for value in matrix:
        zero(value)


def invariants():
    e = s.Matrix(s.symbols('e1:4'))
    b = s.Matrix(s.symbols('b1:4'))
    F = s.zeros(4)
    for i in range(3):
        F[0, i+1], F[i+1, 0] = e[i], -e[i]
        for j in range(3):
            F[i+1, j+1] = sum(s.LeviCivita(i, j, k)*b[k] for k in range(3))
    metric = s.diag(-1, 1, 1, 1)
    # Choose the oriented contravariant Levi-Civita density epsilon^0123=+1.
    dual = s.Matrix(4, 4, lambda i, j: sum(s.LeviCivita(i, j, k, l)*F[k, l]/2 for k in range(4) for l in range(4)))
    J = s.expand(sum(F[i, j]*dual[i, j] for i in range(4) for j in range(4)))
    I = s.expand(sum(F[i, j]**2*metric[i, i]*metric[j, j] for i in range(4) for j in range(4)))
    I4 = s.expand(s.trace((F*metric)**4))
    zero(J-4*e.dot(b))
    zero(I-2*(b.dot(b)-e.dot(e)))
    zero(J*J-(4*I4-2*I*I))
    zeros(F*dual.T-s.eye(4)*J/4)
    coeff, kap2, B0, t = s.symbols('cJ kappa2 B0 t', real=True)
    perturb = {e[i]: t*e[i] for i in range(3)}
    perturb.update({b[i]: t*b[i]+(B0 if i == 2 else 0) for i in range(3)})
    local = coeff*kap2**2*J**2/16
    series = s.expand(local.subs(perturb, simultaneous=True))
    zero(series.coeff(t, 0))
    zero(series.coeff(t, 1))
    zero(series.coeff(t, 2)-coeff*kap2**2*B0*B0*e[2]**2)
    # Direct constitutive response from the six independent F components.
    lag = -I/4+local
    response = s.hessian(lag, list(e)).subs({**{v: 0 for v in e}, b[0]: 0, b[1]: 0, b[2]: B0})
    zeros(response-s.diag(1, 1, 1+2*coeff*kap2**2*B0**2))
    # Four transverse, null photon configurations with incoming covectors
    # k1=(1,0,0,1), k2=(1,0,0,-1), k3=-k1, k4=-k2.
    # Polarizations x,y,x,y. F3=-F1, F4=-F2, and momenta sum to zero.
    t1, t2, t3, t4 = s.symbols('t1 t2 t3 t4')
    k1, k2 = s.Matrix([1, 0, 0, 1]), s.Matrix([1, 0, 0, -1])
    pol1, pol2 = s.Matrix([0, 1, 0, 0]), s.Matrix([0, 0, 1, 0])
    photons = []
    for momentum, pol in [(k1, pol1), (k2, pol2), (-k1, pol1), (-k2, pol2)]:
        ff = momentum*pol.T-pol*momentum.T
        zero((momentum.T*metric*momentum)[0])
        zero((momentum.T*metric*pol)[0])
        zeros(momentum.T*metric*ff)
        photons.append(ff)
    totalF = sum((weight*ff for weight, ff in zip([t1, t2, t3, t4], photons)), s.zeros(4))
    wave_e = s.Matrix([t1-t3, t2-t4, 0])
    wave_b = s.Matrix([t2-t4, t1-t3, 0])
    zeros(s.Matrix([totalF[0, i] for i in range(1, 4)])-wave_e)
    zeros(s.Matrix([totalF[2, 3], totalF[3, 1], totalF[1, 2]])-wave_b)
    waveJ = 4*wave_e.dot(wave_b)
    four_coefficient = s.expand(waveJ**2).coeff(t1).coeff(t2).coeff(t3).coeff(t4)
    assert four_coefficient == 256
    # The counterterm stress alone is -L_J g, since J depends on inverse
    # metric only through 1/sqrt(-g). Its Ward identity follows from
    # F_(nu rho) dualF^(mu rho)=delta_nu^mu J/4 and dF=0.
    dj = s.Matrix(s.symbols('dJ0:4'))
    divH = -coeff*kap2**2*dual.T*dj/2
    ward = F*divH+coeff*kap2**2*J*dj/8
    zeros(ward)
    return {'F_dual_F': str(J), 'basis_relation': 'J^2=4 Tr((F g^-1)^4)-2(F^2)^2',
            'magnetic_background_action_and_first_variation': 'zero',
            'quadratic_correction': 'cJ*kappa^4*B0^2*e_parallel^2',
            'parallel_electric_susceptibility': '1+2*cJ*kappa^4*B0^2',
            'four_null_photon_mixed_coefficient_of_J_squared': str(four_coefficient),
            'constitutive_and_stress_Ward_identity': 'exact, using dF=0'}


def dualization():
    R, kap2 = s.symbols('R kappa2', positive=True)
    sig, coeff = s.symbols('sigma cJ', real=True)
    ds = s.Matrix(s.symbols('sigma_t sigma_x sigma_c'))
    metric = s.diag(-1, 1, 1)
    K = 2*s.pi*R*s.Matrix([[R**2/(2*kap2)+sig**2, sig], [sig, 1]])
    Ki = s.Matrix([[2*kap2/R**2, -2*kap2*sig/R**2], [-2*kap2*sig/R**2, 1+2*kap2*sig**2/R**2]])/(2*s.pi*R)
    zeros(K*Ki-s.eye(2))
    w = s.Matrix([sig, 1])
    zeros(Ki*w-s.Matrix([0, 1/(2*s.pi*R)]))
    A = s.kronecker_product(K, metric)
    Ai = s.kronecker_product(Ki, metric)
    u = s.kronecker_product(w, ds)
    rank = 4*s.pi*coeff*kap2**2/R
    v = Ai*u
    eta = 2*coeff*kap2**2*(ds.T*metric*ds)[0]/R**2
    zero(rank*(u.T*Ai*u)[0]-eta)
    inverse = Ai-rank*v*v.T/(1+eta)
    zeros((A+rank*u*u.T)*inverse-s.eye(6))
    # Dual scalar action obtained by eliminating six independent curvature
    # variables from 1/2 v^T A v - dPhi^T v, not by naive substitution.
    dphi = s.Matrix(s.symbols('pB_t pB_x pB_c pA_t pA_x pA_c'))
    correction = s.factor(-(dphi.T*(inverse-Ai)*dphi)[0]/2)
    expected = coeff*kap2**2*((ds.T*metric*dphi[3:6, 0])[0])**2/(2*s.pi*R**3*(1+eta))
    zero(correction-expected)
    # A three-dimensional horizontal two-form plus d sigma wedge dphi.
    f01, f02, f12 = s.symbols('f01 f02 f12')
    v_horizontal = s.Matrix([f12, -f02, f01])
    J = 4*(f01*ds[2]-f02*ds[1]+f12*ds[0])/R
    zero(J-4*(ds.T*v_horizontal)[0]/R)
    return {'rank_one_update': 'A=K tensor gamma + (4*pi*cJ*kappa^4/R) u u^T, u=(sigma,1) tensor d sigma',
            'eta': str(eta), 'inverse': 'Sherman-Morrison identity verified exactly for all 6 curvature variables',
            'dual_scalar_correction': str(correction),
            'determinant_ratio': 'det(A_updated)/det(A)=1+eta; finite Gaussian multiplier=(1+eta)^(-1/2)',
            'scope': 'local bosonic axisymmetric reduction on a prescribed even background; fixed flux sectors and boundary Legendre term retained'}


def spectrum_algebra():
    C, L, theta = s.symbols('C L theta', positive=True)
    off = s.sqrt(2*(C+1)*(L-2))
    M0 = s.Matrix([[L+2*(C+1), off], [off, L-2]])/C
    M = s.Matrix([[L/(1+theta)+2*(C+1), off], [off, L-2]])/C
    r = theta/(1+theta)
    zeros(M-(M0-r*L*s.diag(1, 0)/C))
    zero(M.det()-L*(L-2)/((1+theta)*C**2))
    zero(s.trace(M)-(2+L*(2+theta)/((1+theta)*C)))
    zero(M0.inv()[0, 0]-C/L)
    projector = M0.inv()*L*s.diag(1, 0)/C
    zeros(projector**2-projector)
    zero(s.trace(projector)-1)
    dtheta = 2*(C+1)/C+2/((1+theta)*C)
    zero(dtheta-(2+4/C-2*r/C))
    # All-m harmonics: Maxwell radial electric energy is a scalar sphere
    # harmonic square. Check the m=0 norm algebra for l=1..10 explicitly.
    c = s.symbols('c')
    for ell in range(1, 11):
        val = ell*(ell+1)
        A = (1-c*c)*s.diff(s.legendre(ell, c), c)
        wa = s.Rational(2*val, 2*ell+1)
        zero(s.integrate(s.diff(A, c)**2, (c, -1, 1))/wa-val)
    return {'all_l_ge_2_mass_matrix': [[str(x) for x in row] for row in M.tolist()],
            'determinant': 'L*(L-2)/[(1+theta)*C^2]',
            'dipole_mass': str(dtheta),
            'relative_update_rank': 1,
            'relative_update_nonzero_eigenvalue': '-theta/(1+theta)',
            'comparison': 'min(1,1/(1+theta))*M0 <= M(theta) <= max(1,1/(1+theta))*M0 for theta>-1',
            'EFT_linear_relative_shift_bound': '|theta|; exact selected-model remainder <=theta^2/(1-|theta|) for |theta|<1'}


def published_master_comparison():
    # Nomura--Yoshida--Soda 2004.07560v2, equations 70-72 and 83.
    # Their action contains -calL, with calF=F^2/4, calG=F*F/4.
    calF, calG, coeff = s.symbols('calF calG cJ', real=True)
    kap2, charge, radius, C, L = s.symbols('kappa2 P b C L', positive=True)
    theta = s.symbols('theta', real=True)
    calL = calF-coeff*kap2**2*calG**2
    LF = s.diff(calL, calF)
    LGG = s.diff(calL, calG, 2)
    denominator = LF-charge**2*LGG/radius**4
    zero(denominator-(1+2*coeff*kap2**2*charge**2/radius**4))
    # The derivative terms in the ordinary master potential also vanish.
    zero(s.diff(calL, calF, 2))
    zero(s.diff(calL, calF, 3))
    # Swap their gravitational/electromagnetic order to ours.
    V = s.Matrix([[L/(radius**2*denominator)+2*kap2*charge**2/radius**4,
                   charge*s.sqrt(2*kap2*(L-2))/radius**3],
                  [charge*s.sqrt(2*kap2*(L-2))/radius**3,
                   (L-2)/radius**2]])
    # At the extremal product kappa^2 P^2/b^2=C+1, a^2=b^2/C.
    replacement = {coeff: theta*radius**4/(2*kap2**2*charge**2)}
    product = ((radius**2/C)*V.subs(replacement)).subs(charge, radius*s.sqrt((C+1)/kap2))
    expected = s.Matrix([[L/(1+theta)+2*(C+1), s.sqrt(2*(C+1)*(L-2))],
                         [s.sqrt(2*(C+1)*(L-2)), L-2]])/C
    zeros(s.simplify(product-expected))
    dipole = (radius**2/C)*(2/(radius**2*denominator)+2*kap2*charge**2/radius**4)
    dipole = dipole.subs(replacement).subs(charge**2, radius**2*(C+1)/kap2)
    zero(dipole-(2/(1+theta)+2*(C+1))/C)
    return {'source': 'https://arxiv.org/html/2004.07560v2',
            'equations': ['70', '71', '72', '83'],
            'source_lagrangian': 'calL=calF-cJ*kappa^4*calG^2; action contains -calL',
            'source_constitutive_denominator': '1+2*cJ*kappa^4*P^2/b^4=1+theta',
            'all_l_ge_2_product_matrix': 'exact match after swapping field order',
            'dipole_product_mass': 'exact match',
            'novelty': 'application and normalization check; the general NLED perturbation dependence is established prior work'}


def numerical_dual_checks():
    rows = []
    # Both Lorentzian and positive Euclidean Gaussian matrices. No continuum
    # quantum determinant is inferred from these finite integrals.
    for signature in [-1., 1.]:
        for coeff in [-.02, .015, .07]:
            R, kap2, sig = 1.3, .8, -.6
            metric = np.array([[signature*1.7, .08, -.04], [.08, 1.2, .06], [-.04, .06, .9]])
            mi = np.linalg.inv(metric)
            K = 2*np.pi*R*np.array([[R*R/(2*kap2)+sig*sig, sig], [sig, 1.]])
            Ai = np.kron(np.linalg.inv(K), mi)
            A = np.kron(K, metric)
            ds = np.array([.21, -.35, .72])
            u = np.kron([sig, 1.], ds)
            rank = 4*np.pi*coeff*kap2**2/R
            eta = 2*coeff*kap2**2*(ds@mi@ds)/(R*R)
            assert 1+eta > 0
            updated = A+rank*np.outer(u, u)
            predicted = Ai-rank*np.outer(Ai@u, Ai@u)/(1+eta)
            error = float(np.max(np.abs(np.linalg.inv(updated)-predicted)))
            ratio = np.linalg.det(updated)/np.linalg.det(A)
            assert error < 2e-13 and abs(ratio-(1+eta)) < 2e-13
            gaussian = None
            if signature == 1:
                assert np.linalg.eigvalsh(updated)[0] > 0
                gaussian = float(np.prod(np.diag(np.linalg.cholesky(A)))/np.prod(np.diag(np.linalg.cholesky(updated))))
                assert abs(gaussian-1/np.sqrt(1+eta)) < 2e-13
            rows.append({'signature': signature, 'cJ': coeff, 'inverse_error': error, 'determinant_ratio_error': float(abs(ratio-1-eta)), 'Euclidean_Gaussian_ratio': gaussian})
    return rows


def mass_data(C, ell, theta):
    L = mp.mpf(ell*(ell+1))
    off = mp.sqrt(2*(C+1)*(L-2))/C
    matrix = mp.matrix([[(L/(1+theta)+2*(C+1))/C, off], [off, (L-2)/C]])
    trace = 2+L*(2+theta)/((1+theta)*C)
    det = L*(L-2)/((1+theta)*C*C)
    heavy = (trace+mp.sqrt(trace*trace-4*det))/2
    return matrix, [det/heavy, heavy]


def mass_checks():
    mp.mp.dps = 110
    rows = []
    for C in map(mp.mpf, ['1', '4', '100', '1e12']):
        for ell in [2, 3, 8, 100]:
            M0, base = mass_data(C, ell, mp.mpf(0))
            for theta in map(mp.mpf, ['-.2', '-.01', '.01', '.2']):
                M, mass = mass_data(C, ell, theta)
                lower, upper = min(1, 1/(1+theta)), max(1, 1/(1+theta))
                for old, new in zip(base, mass):
                    assert lower*old <= new <= upper*old
                # Direct high-precision eigenvalues versus rationalized formula.
                eigen, U = mp.eigsy(M)
                relative = max(abs(eigen[p]/mass[p]-1) for p in range(2))
                assert relative < mp.mpf('1e-65')
                rows.append({'C': str(C), 'l': ell, 'theta': str(theta), 'light_mass_ratio': mp.nstr(mass[0]/base[0], 20), 'direct_eigenvalue_relative_error': mp.nstr(relative, 8)})
    k = mp.mpf(2000)
    C = mp.sqrt(1+12*mp.pi*k**6)
    gamma = 6*C/(k**10*(C-1))
    factor = 16*mp.pi*gamma*(1+1/C)
    candidate = {'C': str(C), 'G_over_a2': str(gamma), 'theta_per_unit_cJ': str(factor),
                 'cJ_absolute_limit_for_one_percent_two_sided_mass_comparison': str(1/(101*factor)),
                 'positive_cJ_for_theta_0p01': str(mp.mpf('.01')/factor), 'cases': []}
    for coeff in [mp.mpf(0), mp.mpf(1), mp.mpf('.01')/factor]:
        theta = factor*coeff
        base = mass_data(C, 3, mp.mpf(0))[1]
        masses = mass_data(C, 3, theta)[1]
        candidate['cases'].append({'cJ': str(coeff), 'theta': str(theta), 'l3_light_relative_shift': str(masses[0]/base[0]-1), 'magnetic_background_shift_from_this_operator': '0 exactly'})
    return {'parameter_cases': rows, 'candidate_k2000': candidate}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print('Checking covariant invariants, background blindness and constitutive Ward identity', flush=True)
    result = {'invariants': invariants()}
    print('Checking constrained dual action, measure factor and all-angular mass comparison', flush=True)
    result['dualization'] = dualization()
    result['spectrum'] = spectrum_algebra()
    result['published_master_comparison'] = published_master_comparison()
    result['finite_dual_checks'] = numerical_dual_checks()
    result['mass_checks'] = mass_checks()
    result.update({'status': 'PASS', 'scope': 'one specified finite local EFT coefficient; full S074 action, counterterms, state and opening remain incomplete', 'versions': {'python': platform.python_version(), 'sympy': s.__version__, 'numpy': np.__version__, 'mpmath': mp.__version__}, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('PASS: finite coefficient is independent of magnetic background; axial response and admissible range verified', flush=True)


if __name__ == '__main__':
    main()
