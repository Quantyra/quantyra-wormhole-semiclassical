"""S074: spacelike axial reduction and a mixed l=1/l=3 cubic invariant.

This is a bulk classical Hessian on prescribed even backgrounds. It is not
the full constrained quantum action, a loop bound, or a stability certificate.
Run with --output PATH. No Lean or external compute is invoked.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform

import mpmath as mp
import numpy as np
import sympy as s
from sympy.physics.wigner import wigner_3j

ROOT = Path(__file__).resolve().parents[1]
C, lam = s.symbols('C lambda', positive=True)
c = s.symbols('c', real=True)
ss = 1-c*c
Y = s.legendre(2, c)
V = ss*s.diff(Y, c)
P = s.sqrt((C+1)/lam)


def zero(expr):
    value = canonical(s.expand(expr))
    assert value == 0, str(value)


def canonical(expr):
    # Factor square-root arguments before simplification. With C>0,
    # sqrt(C^2+2C+1)=C+1; do not introduce a numerical zero tolerance.
    expr = expr.replace(lambda x: x.is_Pow and x.exp.is_Rational and x.exp.q == 2,
                        lambda x: s.Pow(s.factor(x.base), x.exp))
    return s.simplify(expr)


def zeros(mat):
    for value in mat:
        zero(value)


def moment(expr):
    poly = s.Poly(s.cancel(expr), c)
    return s.simplify(sum(coef*s.Rational(2, n[0]+1)
                          for n, coef in poly.terms() if n[0] % 2 == 0))


def angular(mat):
    return mat.applyfunc(lambda x: canonical(moment(x)/(2*s.pi)))


def curvature_scalar(g, dg, ddg):
    n = g.rows
    gi = g.inv()
    dgi = [-gi*d*gi for d in dg]
    gam, dgam = {}, {}
    for i in range(n):
        for j in range(n):
            for k in range(n):
                gam[i, j, k] = sum(gi[i, z]*(dg[j][z, k]+dg[k][z, j]-dg[z][j, k])/2 for z in range(n))
                for v in range(n):
                    dgam[v, i, j, k] = sum((dgi[v][i, z]*(dg[j][z, k]+dg[k][z, j]-dg[z][j, k])+gi[i, z]*(ddg[v][j][z, k]+ddg[v][k][z, j]-ddg[v][z][j, k]))/2 for z in range(n))
    ric = s.zeros(n)
    for i in range(n):
        for j in range(n):
            ric[i, j] = sum(dgam[k, k, i, j]-dgam[j, k, i, k]+sum(gam[k, i, j]*gam[z, k, z]-gam[z, i, k]*gam[k, j, z] for z in range(n)) for k in range(n))
    return s.expand(s.trace(gi*ric))


def four_dimensional_action():
    # Independent metric jets for all three connection curvatures, with
    # arbitrary first and second c derivatives of the even diagonal metric.
    a, b, h, r = s.symbols('a b h r', positive=True)
    ap, bp, hp, rp, app, bpp, hpp, rpp = s.symbols('ap bp hp rp app bpp hpp rpp')
    f01, f02, f12 = s.symbols('f01 f02 f12')
    f = s.Matrix([[0, f01, f02], [-f01, 0, f12], [-f02, -f12, 0]])
    dconn = f/2  # partial_M B_N; B=0 at the test point
    g = s.diag(-a, b, h, r*r)
    dg = [s.zeros(4) for _ in range(4)]
    ddg = [[s.zeros(4) for _ in range(4)] for _ in range(4)]
    dg[2] = s.diag(-ap, bp, hp, 2*r*rp)
    ddg[2][2] = s.diag(-app, bpp, hpp, 2*rp*rp+2*r*rpp)
    for u in range(3):
        for m in range(3):
            dg[u][m, 3] = dg[u][3, m] = r*r*dconn[u, m]
        for v in range(3):
            for m in range(3):
                ddg[u][v][m, 3] = ddg[u][v][3, m] = 2*r*rp*(int(u == 2)*dconn[v, m]+int(v == 2)*dconn[u, m])
                for n in range(3):
                    ddg[u][v][m, n] += r*r*(dconn[u, m]*dconn[v, n]+dconn[v, m]*dconn[u, n])
    ric4 = curvature_scalar(g, dg, ddg)
    base = g[:3, :3]
    # Build the even base jets independently of the connection.
    bdg = [s.zeros(3), s.zeros(3), s.diag(-ap, bp, hp)]
    bddg = [[s.zeros(3) for _ in range(3)] for _ in range(3)]
    bddg[2][2] = s.diag(-app, bpp, hpp)
    ric3 = curvature_scalar(base, bdg, bddg)
    fsq = sum(f[i, j]**2/(base[i, i]*base[j, j]) for i in range(3) for j in range(3))
    lapr = (rpp+(ap/(2*a)+bp/(2*b)-hp/(2*h))*rp)/h
    zero(ric4-ric3+r*r*fsq/4+2*lapr/r)

    sig, kap, gauge = s.symbols('sigma kappa2 gauge', positive=True)
    gs = s.symbols('g01 g02 g12')
    gg = s.Matrix([[0, gs[0], gs[1]], [-gs[0], 0, gs[2]], [-gs[1], -gs[2], 0]])
    dsig = s.symbols('sigma0:3')
    em = s.zeros(4)
    em[:3, :3] = gg+sig*f
    for i in range(3):
        em[i, 3], em[3, i] = dsig[i], -dsig[i]
    emsq = sum(em[i, j]**2/(g[i, i]*g[j, j]) for i in range(4) for j in range(4))
    zero(emsq-sum((gg[i, j]+sig*f[i, j])**2/(base[i, i]*base[j, j]) for i in range(3) for j in range(3))-2*sum(dsig[i]**2/base[i, i] for i in range(3))/(r*r))
    K = 2*s.pi*r*s.Matrix([[r*r/(2*kap)+sig*sig/gauge**2, sig/gauge**2], [sig/gauge**2, 1/gauge**2]])
    Ki = s.Matrix([[2*kap/(r*r), -2*kap*sig/(r*r)], [-2*kap*sig/(r*r), gauge**2+2*kap*sig*sig/(r*r)]])/(2*s.pi*r)
    zeros(K*Ki-s.eye(2))
    # First-order Lorentzian dualization: +Phi dF, then integrate by parts.
    # Flat oriented frame density: 1/2(F01 K F01+F02 K F02-F12 K F12)
    # -d0Phi F12+d1Phi F02-d2Phi F01.
    v0, v1, v2 = [s.Matrix(s.symbols(f'v{i}_0:2')) for i in range(3)]
    f01v, f02v, f12v = Ki*v2, -Ki*v1, -Ki*v0
    dual = ((f01v.T*K*f01v+f02v.T*K*f02v-f12v.T*K*f12v)/2-v0.T*f12v+v1.T*f02v-v2.T*f01v)[0]
    zero(dual-(v0.T*Ki*v0-v1.T*Ki*v1-v2.T*Ki*v2)[0]/2)
    return {'Ricci_identity': 'R4=R3-R^2 F_B^2/4-2 Laplacian(R)/R', 'Maxwell_identity': 'F4^2=(dA+sigma dB)^2+2(d sigma)^2/R^2', 'K_inverse_and_Lorentzian_dualization': 'exact', 'jet_scope': 'all three F_B components; arbitrary diagonal even metric c-jets through second order'}


def harmonic_action():
    norms = []
    for ell in range(1, 9):
        L = ell*(ell+1)
        pl = s.legendre(ell, c)
        aa, xx = ss*s.diff(pl, c), ss**2*s.diff(pl, c, 2)
        wa = s.Rational(2*L, 2*ell+1)
        zero(s.diff(aa, c)+L*pl)
        zero(s.diff(xx, c)+(L-2)*aa)
        zero(moment(aa**2/ss)-wa)
        zero(moment(xx**2/ss**2)-wa*(L-2))
        norms.append({'l': ell, 'wA': str(wa), 'wX': str(wa*(L-2))})
    L = s.symbols('L', positive=True)
    matrix = s.Matrix([[L+2*(C+1), s.sqrt(2*(C+1)*(L-2))], [s.sqrt(2*(C+1)*(L-2)), L-2]])/C
    zero(s.trace(matrix)-(2+2*L/C))
    zero(matrix.det()-L*(L-2)/C**2)
    profiles = []
    for ell in [1, 3]:
        L = ell*(ell+1)
        pl = s.legendre(ell, c)
        wa = s.Rational(2*L, 2*ell+1)
        profiles.append((ss*s.diff(pl, c)/s.sqrt(wa/(2*s.pi)), s.S.Zero))
        if ell >= 2:
            profiles.append((s.S.Zero, ss**2*s.diff(pl, c, 2)/s.sqrt(lam*wa*(L-2)/s.pi)))
    A = s.Matrix([[v[0] for v in profiles]])
    X = s.Matrix([[v[1] for v in profiles]])
    Ap, Xp = A.diff(c), X.diff(c)
    D = Xp-P*A
    E = P*V*Ap
    Z = P*V*A
    T0 = A.T*A/ss+2*lam*X.T*X/ss**2
    U0 = Ap.T*Ap/C+2*lam*D.T*D/(C*ss)
    Tc = -4*lam*Y*X.T*X/ss**2
    Ta = 2*lam*(X.T*Z+Z.T*X)/ss**2
    mats = {
        'K0': angular(T0), 'M0': angular(U0),
        'Kchi': angular(Tc), 'Kalpha': angular(Ta), 'Kh': angular(-Y*T0),
        'Vchi': angular(-2*Y*Ap.T*Ap/C-8*lam*Y*D.T*D/(C*ss)),
        'Valpha': angular(2*lam*(D.T*E+E.T*D)/(C*ss)),
        # Second-order coefficients, without Taylor factorials.
        'Kchichi': angular(4*lam*Y**2*X.T*X/ss**2),
        'Kchialpha': angular(-4*lam*Y*(X.T*Z+Z.T*X)/ss**2),
        'Kalphaalpha': angular(2*lam*Z.T*Z/ss**2),
        'KJ': angular(Y**2*T0/2),
        'Khchi': angular(-Y*Tc), 'Khalpha': angular(-Y*Ta),
        'Vchichi': angular(2*Y**2*Ap.T*Ap/C+16*lam*Y**2*D.T*D/(C*ss)),
        'Vchialpha': angular(-8*lam*Y*(D.T*E+E.T*D)/(C*ss)),
        'Valphaalpha': angular(2*lam*E.T*E/(C*ss)),
        'VJ': angular(-Y**2*U0/2),
    }
    zeros(mats['K0']-s.eye(3))
    expected = s.Matrix([[2+4/C, 0, 0], [0, 2+14/C, 2*s.sqrt(5*(C+1))/C], [0, 2*s.sqrt(5*(C+1))/C, 10/C]])
    zeros(mats['M0']-expected)
    zeros(mats['Kchi'])
    zeros(mats['Kh']-s.Matrix([[s.Rational(1, 5), -3*s.sqrt(14)/35, 0], [-3*s.sqrt(14)/35, -s.Rational(1, 5), 0], [0, 0, 0]]))
    return mats, (A, X, Ap, Xp), norms


def mixed_stress_identity():
    # Unequal masses and frequencies; curved AdS2 B=1+x^2, arbitrary jets.
    x, w, v, d, n = s.symbols('x w v d n', real=True)
    f, g = s.Function('f')(x), s.Function('g')(x)
    B = 1+x*x
    substitutions = {s.diff(f, x, 2): (d-w*w/B)*f/B-s.diff(B, x)*s.diff(f, x)/B,
                     s.diff(g, x, 2): (n-v*v/B)*g/B-s.diff(B, x)*s.diff(g, x)/B}
    dot = w*v*f*g/B+B*s.diff(f, x)*s.diff(g, x)
    tt = -w*v*f*g/B**2+dot/(2*B)
    tx = s.I*(w*f*s.diff(g, x)+v*g*s.diff(f, x))/2
    xx = B**2*s.diff(f, x)*s.diff(g, x)-B*dot/2
    dt = -s.I*(w+v)*tt+s.diff(tx, x)+s.diff(B, x)*tx/B
    dx = -s.I*(w+v)*tx+s.diff(xx, x)+B*s.diff(B, x)*tt/2-s.diff(B, x)*xx/(2*B)
    zero(s.expand(dt-s.I*(d*v+n*w)*f*g/(2*B)).subs(substitutions))
    zero(s.expand(dx-B*(d*f*s.diff(g, x)+n*g*s.diff(f, x))/2).subs(substitutions))
    # Divergence once more yields ((d+n) Box(fg)-(d-n)^2 fg)/4.
    zero((d+n)*dot/2+d*n*f*g-((d+n)*((d+n)*f*g+2*dot)-(d-n)**2*f*g)/4)
    return 'HessTF chi * grad(f) grad(g) integrates to [(d+n) Box(chi)-(d-n)^2 chi] f g /4, plus boundary terms'


def invariant(mats):
    d = 2+4/C
    Mp = s.Matrix([[2+6/C, s.sqrt(12*(C+1))/C], [s.sqrt(12*(C+1))/C, 6/C]])
    W = s.Matrix([[2*(1+C/3)/lam, 2*P], [2*P, 6]])
    inverse_scale = s.diag(s.sqrt(3*lam), 1)
    zeros(inverse_scale*W*inverse_scale-C*Mp)
    Ma = mats['M0'][1:3, 1:3]
    vc, va, hr, kr = [mats[name][0:1, 1:3] for name in ['Vchi', 'Valpha', 'Kh', 'Kalpha']]
    tc, ta = s.Matrix([s.sqrt(3), 0]), s.Matrix([0, 1/s.sqrt(C+1)])
    B = tc*vc+ta*va-C*tc*hr*(d*s.eye(2)-Ma)**2/6-ta*kr*(d*s.eye(2)+Ma)/2+Mp*(C*tc*hr*(d*s.eye(2)+Ma)/6+ta*kr/2)
    B = B.applyfunc(canonical)
    printed_B = s.Matrix([[-2*s.sqrt(42)*(2*C*C+22*C+27)/(35*C), -4*(C+2)*s.sqrt(210*(C+1))/(35*C)], [-12*(C+8)*s.sqrt(14*(C+1))/(35*C), -6*s.sqrt(70)*(2*C+1)/(35*C)]])
    zeros(B-printed_B)
    q = s.factor(s.trace(Mp.inv()**2*B*B.T)*5/(4*s.pi*C))
    zero(q-(53*C*C+830*C+944)/(280*s.pi*C))
    # Regression against the previous independently derived diagonal source.
    z = s.symbols('z')
    direct_chi = mats['Vchi'][0, 0]+C*mats['Kh'][0, 0]*2*d*z/6
    direct_a = mats['Valpha'][0, 0]
    # S3=-1/2 q1^2*(...); Box chi=(2+6/C)chi+6 lambda P alpha/C.
    coeff_chi = s.simplify(-(direct_chi.subs(z, 2+6/C))/2)
    coeff_alpha = s.simplify(-(direct_a+C*mats['Kh'][0, 0]*2*d*(6*(C+1)/C)/6)/2)
    zero(coeff_chi+2*(C*C+11*C+6)/(15*C))
    zero(coeff_alpha+2*(C+1)*(C+8)/(5*C))
    angular_rows = []
    reference = wigner_3j(1, 2, 3, 0, 0, 0)**2
    for m2 in range(-2, 3):
        total = sum(wigner_3j(1, 2, 3, m1, m2, m3)**2 for m1 in range(-1, 2) for m3 in range(-3, 4))
        zero(total-s.Rational(1, 5))
        zero(total/reference-s.Rational(7, 3))
        angular_rows.append({'m2': m2, 'sum_of_squares': str(total), 'relative_to_m0': str(total/reference)})
    return B, {'B': [[str(x) for x in row] for row in B.tolist()], 'sum_g_squared_over_mu_m0': 'lambda*(53*C^2+830*C+944)/(280*pi*C)', 'all_m_fixed_external_m': '(G/a^2)*(53+830/C+944/C^2)/15', 'uniform_bound_C_ge_1': '(1827/15)*(G/a^2)', 'diagonal_source_coefficients_chi_alpha': [str(coeff_chi), str(coeff_alpha)], 'angular_orthogonality': angular_rows}


def numerical_contacts(mats, profiles):
    # Independent Gauss integration of the unexpanded nonlinear 3D action.
    # Check its second-order Taylor action, including q!=0 base shear.
    nodes, weights = np.polynomial.legendre.leggauss(100)
    rows = []
    for cv in [1., 4., 100.]:
        lv = .37
        pv = np.sqrt((cv+1)/lv)
        pf = [s.lambdify(c, m.subs({C: cv, lam: lv}), 'numpy') for m in profiles]
        numeric = {k: np.array(m.subs({C: cv, lam: lv}), dtype=float) for k, m in mats.items()}
        direction = np.array([.21, -.13, .17, -.09])
        errors = []
        for eps in [.01, .005, .0025]:
            chi, alpha, p, shear = eps*direction
            J = p*p-shear*shear
            baseh = np.array([[p, shear], [shear, p]])
            eta = np.diag([-1., 1.])
            hh = eta@baseh@eta
            kinetic = np.zeros((2, 2, 3, 3))
            mass = np.zeros((3, 3))
            for cc, wt in zip(nodes, weights):
                aa, xx, ap, xp = [np.asarray(f(cc), dtype=float).reshape(3) for f in pf]
                sv, yy, vv = 1-cc*cc, (3*cc*cc-1)/2, 3*cc*(1-cc*cc)
                gg = eta+yy*baseh
                detroot = np.sqrt(-np.linalg.det(gg))
                xeff = xx+pv*alpha*vv*aa
                deff = xp-pv*aa+pv*alpha*vv*ap
                T = np.outer(aa, aa)/sv+2*lv*np.exp(-2*chi*yy)*np.outer(xeff, xeff)/sv**2
                U = np.exp(-2*chi*yy)*np.outer(ap, ap)/cv+2*lv*np.exp(-4*chi*yy)*np.outer(deff, deff)/(cv*sv)
                kinetic += wt*np.einsum('ab,ij->abij', detroot*np.linalg.inv(gg), T)/(2*np.pi)
                mass += wt*detroot*U/(2*np.pi)
            kn = numeric['K0']+chi*numeric['Kchi']+alpha*numeric['Kalpha']+chi*chi*numeric['Kchichi']+chi*alpha*numeric['Kchialpha']+alpha*alpha*numeric['Kalphaalpha']+J*numeric['KJ']
            hn = numeric['Kh']+chi*numeric['Khchi']+alpha*numeric['Khalpha']
            kt = np.einsum('ab,ij->abij', eta, kn)+np.einsum('ab,ij->abij', hh, hn)
            mt = numeric['M0']+chi*numeric['Vchi']+alpha*numeric['Valpha']+chi*chi*numeric['Vchichi']+chi*alpha*numeric['Vchialpha']+alpha*alpha*numeric['Valphaalpha']+J*numeric['VJ']
            err = max(np.max(np.abs(kinetic-kt)), np.max(np.abs(mass-mt)))
            assert err < 200*eps**3, (cv, eps, err)
            errors.append(float(err))
        ratios = [errors[i]/errors[i+1] for i in range(2)]
        assert min(ratios) > 7.8 and max(ratios) < 8.2, ratios
        rows.append({'C': cv, 'amplitudes': [.01, .005, .0025], 'max_absolute_remainders': errors, 'halving_ratios': ratios})
    return rows


def eigenbasis_comparison(B):
    mp.mp.dps = 90
    bfun = s.lambdify(C, B, 'mpmath')
    rows = []
    kval = mp.mpf(2000)
    candidate = mp.sqrt(1+12*mp.pi*kval**6)
    for cv in [mp.mpf(1), mp.mpf(4), mp.mpf(100), mp.mpf('1e12'), candidate]:
        rp = mp.sqrt(cv*cv+12*cv+12)
        ra = mp.sqrt(cv*cv+24*cv+24)
        pp, pm = (cv+6+rp)/cv, 24/(cv*(cv+6+rp))
        ap, am = (cv+12+ra)/cv, 120/(cv*(cv+12+ra))
        cp = mp.sqrt((rp+cv)/(2*rp))
        sp = mp.sqrt(12*(cv+1))/mp.sqrt(2*rp*(rp+cv))
        ca = mp.sqrt((ra+cv+2)/(2*ra))
        sa = mp.sqrt(20*(cv+1))/mp.sqrt(2*ra*(ra+cv+2))
        pol = [(pm, mp.matrix([-sp, cp])), (pp, mp.matrix([cp, sp]))]
        axial = [(am, mp.matrix([-sa, ca])), (ap, mp.matrix([ca, sa]))]
        d = 2+4/cv
        total = mp.mpf(0)
        branch_rows = []
        for nu, ea in axial:
            gs = []
            for z, ep in pol:
                cc = mp.sqrt(3)*ep[0]/mp.sqrt(4*mp.pi*cv*z/5)
                ac = ep[1]/mp.sqrt((cv+1)*4*mp.pi*cv*z/5)
                fc = [-12*mp.sqrt(14)*(2*cv+5)/(35*cv)-3*mp.sqrt(14)*cv*((d+nu)*z-(d-nu)**2)/210, -24*mp.sqrt(70*(cv+1))/(35*cv)]
                fa = [-12*mp.sqrt(14)*(cv+1)/(35*cv), 24*mp.sqrt(70*(cv+1))/(35*cv)+3*mp.sqrt(70*(cv+1))*(z-d-nu)/35]
                direct = cc*mp.fdot(ea, fc)+ac*mp.fdot(ea, fa)
                matrix = (ep.T*bfun(cv)*ea)[0]/mp.sqrt(4*mp.pi*cv*z/5)
                assert abs(direct-matrix) < mp.mpf('1e-60')*max(1, abs(direct)), (cv, direct, matrix)
                gs.append(direct)
                total += direct**2/z
            branch_rows.append({'mu_axial': mp.nstr(nu, 22), 'g_over_sqrtlambda': [mp.nstr(x, 22) for x in gs]})
        expected = (53*cv*cv+830*cv+944)/(280*mp.pi*cv)
        rel = abs(total/expected-1)
        assert rel < mp.mpf('1e-60'), (cv, rel)
        rows.append({'C': mp.nstr(cv, 22), 'branches': branch_rows, 'relative_norm_discrepancy': mp.nstr(rel, 8)})
    a2 = kval**10*(candidate-1)/(6*candidate)
    q_all = (53+830/candidate+944/candidate**2)/(15*a2)
    return {'precision_digits': 90, 'cases': rows, 'k2000_G_over_a2': mp.nstr(1/a2, 22), 'k2000_all_m_norm': mp.nstr(q_all, 22), 'scope': 'classical on-shell cubic coefficient norm, not the mixed loop integral or opening error'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print('Checking original curvature, Maxwell and first-order dual action', flush=True)
    result = {'action': four_dimensional_action()}
    print('Checking harmonics and cubic/quartic prescribed-background action', flush=True)
    mats, profiles, norms = harmonic_action()
    result['harmonics'] = norms
    result['matrices'] = {name: [[str(x) for x in row] for row in mat.tolist()] for name, mat in mats.items()}
    print('Checking curved mixed stress and exact inverse-mass coupling norm', flush=True)
    result['mixed_stress'] = mixed_stress_identity()
    B, result['invariant'] = invariant(mats)
    result['nonlinear_contact_quadrature'] = numerical_contacts(mats, profiles)
    result['independent_eigenbasis'] = eigenbasis_comparison(B)
    result.update({'status': 'PASS', 'scope': 'classical bulk axisymmetric axial Hessian, its prescribed polar contacts, and mixed on-shell cubic norm; full S074 gate remains open', 'versions': {'python': platform.python_version(), 'sympy': s.__version__, 'numpy': np.__version__, 'mpmath': mp.__version__}, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('PASS: action, harmonic normalization, mixed vertex and exact coupling invariant', flush=True)


if __name__ == '__main__':
    main()
