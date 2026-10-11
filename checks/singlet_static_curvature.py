"""Reproduce S053 singlet static checks; no Lean, fitting or network calls."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy.linalg as la
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss
import sympy as sy
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT/'notes/singlet-static-curvature.md'
OUT = ROOT/'evidence/singlet-static-curvature.json'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def exact_checks():
    x, e, Q = sy.symbols('x e Q', positive=True)
    b = sy.sqrt(1-e**2*x**2)
    I = Q*b/e
    d = lambda f: sy.simplify(x*sy.diff(f, x))
    checks = {
        'integrated_profile_first_derivative': d(I)+Q*e*x**2/b,
        'integrated_profile_second_derivative': d(d(I))+Q*e*(2*x**2/b+e**2*x**4/b**3),
        'cut_angle_first_derivative': d(sy.acos(e*x))+e*x/b,
        'cut_angle_second_derivative': d(d(sy.acos(e*x)))+e*x/b**3,
    }
    for k, v in checks.items():
        assert sy.simplify(v) == 0, k
    # Independent convergent integrals establish the zero-mode and scheme signs.
    rows = [dict(name=k, residual='0') for k in checks]
    for a in [.001, .2, 3.]:
        observed = quad(lambda u: np.log1p(a/u**2), 0, np.inf, epsabs=1e-10)[0]
        assert abs(observed-np.pi*np.sqrt(a)) < 2e-9
        rows.append(dict(name='zero_mode_integral', a=a, error=abs(observed-np.pi*np.sqrt(a))))
    for m1, m2 in [(.2, 3.), (2., .05)]:
        observed = quad(lambda k: 1/np.sqrt(k*k+m1*m1)-1/np.sqrt(k*k+m2*m2),
                        0, np.inf, epsabs=1e-11)[0]
        assert abs(observed-np.log(m2/m1)) < 1e-10
        rows.append(dict(name='subtraction_scale_sign', mu1=m1, mu2=m2,
                         error=abs(observed-np.log(m2/m1))))
    exponents = dict(integrated_mass=sy.Rational(-3, 4), energy=sy.Rational(-3, 8),
                     equilibrium_shift=sy.Rational(-11, 8), static_relative=sy.Rational(-21, 16),
                     recentering_relative=sy.Rational(-19, 16), dynamic_relative=sy.Rational(-11, 16))
    assert exponents['energy']-1+sy.Rational(3, 16) == exponents['recentering_relative']
    rows.append(dict(name='working_exponents', values={k: str(v) for k, v in exponents.items()}))
    return rows


def derivative_data(A, Ap, App, V, Vp, Vpp):
    R, S = la.inv(A), la.inv(A+V)
    L = np.linalg.slogdet(A+V)[1]-np.linalg.slogdet(A)[1]
    Z = np.trace(S@V)
    D1 = np.trace(S@Vp-(R-S)@Ap)
    D2 = np.trace(S@Vpp-(R-S)@App+R@Ap@R@Ap-S@Ap@S@Ap
                  -2*S@Ap@S@Vp-S@Vp@S@Vp)
    Rp = -R@Ap@R
    Rpp = 2*R@Ap@R@Ap@R-R@App@R
    born = [np.trace(R@V), np.trace(Rp@V+R@Vp), np.trace(Rpp@V+2*Rp@Vp+R@Vpp)]
    return np.array([L, D1, D2]), float(Z), np.array(born)


def matrix_checks():
    rng = np.random.default_rng(731041)
    rows = []
    for n in [4, 9, 17]:
        mats = []
        for _ in range(4):
            U = rng.normal(size=(n, n))
            mats.append(U@U.T/n+.1*np.eye(n))
        K1, K2, K3, W = mats
        vectors = rng.normal(size=(n, 3))/np.sqrt(n)
        alpha = np.array([-2., .3, 1.7])
        for scale in [1e-3, .2, 4.]:
            for k in [.1, 1., 20.]:
                def data(phi):
                    A = np.exp(phi)*K1+np.exp(-phi)*K2+K3+k*k*W
                    Ap, App = np.exp(phi)*K1-np.exp(-phi)*K2, np.exp(phi)*K1+np.exp(-phi)*K2
                    vv = [scale*(vectors*(alpha**j*np.exp(alpha*phi)))@vectors.T for j in range(3)]
                    return A, Ap, App, *vv
                vals, Z, born = derivative_data(*data(.13))
                assert abs(vals[1]) <= 4*Z+1e-10
                assert abs(vals[2]) <= 28*Z+1e-10
                h = .002
                samples = {i: derivative_data(*data(.13+i*h))[0][0] for i in [-2, -1, 0, 1, 2]}
                fd1 = (samples[-2]-8*samples[-1]+8*samples[1]-samples[2])/(12*h)
                fd2 = (-samples[-2]+16*samples[-1]-30*samples[0]+16*samples[1]-samples[2])/(12*h*h)
                error = max(abs(fd1-vals[1]), abs(fd2-vals[2]))
                assert error < 2e-7, (n, scale, k, error)
                row = dict(dimension=n, strength=scale, k=k, logdet=float(vals[0]),
                           trace_fraction=Z, derivatives=vals[1:].tolist(), finite_difference_error=error,
                           first_bound_ratio=abs(vals[1])/(4*Z), second_bound_ratio=abs(vals[2])/(28*Z))
                if born[0] <= .5:
                    row['uv_remainder_over_trace_squared'] = ((vals-born)/born[0]**2).tolist()
                rows.append(row)
    return rows


# The exact Fourier/Bernoulli expansion removes catastrophic subtraction of
# 1/(k^2 L) when the zero oscillator is separated from the determinant.
t = sy.symbols('t')
BERNOULLI = [np.array([float(v) for v in sy.Poly(sy.bernoulli(2*r, t)/sy.factorial(2*r), t).all_coeffs()])
             for r in range(1, 7)]


def green_regular(k, distance, lengths):
    answer = np.zeros_like(distance)
    for L in lengths:
        if k*L < .5:
            u = distance/L
            term = sum(k**(2*r)*L**(2*r+1)*np.polyval(coef, u) for r, coef in enumerate(BERNOULLI))
        else:
            term = (np.exp(-k*distance)+np.exp(-k*(L-distance)))/(2*k*(-np.expm1(-k*L)))-1/(k*k*L)
        answer += term/len(lengths)
    return answer


def fredholm(k, T, returns, y, weights, m2):
    lengths = T+np.asarray(returns)
    distance = T*np.abs(y[:, None]-y[None, :])
    p = np.sqrt(T*weights*m2)
    reg = p[:, None]*green_regular(k, distance, lengths)*p[None, :]
    ch = la.cho_factor(np.eye(len(y))+reg, lower=True, check_finite=False)
    c0 = np.mean(1/lengths)
    L = 2*np.log(np.diag(ch[0])).sum()+np.log1p(c0/k**2*(p@la.cho_solve(ch, p, check_finite=False)))
    return float(L)


def transfer_logdet(k, T, D, m2):
    a = np.sqrt(k*k+m2)
    X, Y = a*T, k*D
    r = (a/k+k/a)/2
    B = (1+r)/2*(1+np.exp(-2*(X+Y)))+(1-r)/2*(np.exp(-2*X)+np.exp(-2*Y))-2*np.exp(-(X+Y))
    return m2*T/(a+k)+np.log(B)-2*np.log(-np.expm1(-k*(T+D)))


def kernel_checks():
    mp.mp.dps = 60
    rows = []
    for L in [3.4, 6., 15.]:
        for k in [1e-5, .03, 1., 10.]:
            d = np.linspace(0, min(3., L), 11)
            calc = green_regular(k, d, [L])
            refs = []
            for dd in d:
                kk, LL, DD = map(mp.mpf, map(str, [k, L, dd]))
                refs.append(float(mp.cosh(kk*(LL/2-DD))/(2*kk*mp.sinh(kk*LL/2))-1/(kk*kk*LL)))
            error = float(np.max(abs(calc-refs)))
            assert error < 2e-12, (L, k, error)
            diag = green_regular(k, np.array([0.]), [L])[0]
            assert 0 <= diag <= L/12+1e-12
            rows.append(dict(length=L, k=k, regular_green_error=error, regular_diagonal=float(diag)))
    for k in [1e-5, .1, 1., 20.]:
        T, D, m2 = np.pi, 1.7, .03
        a = mp.sqrt(mp.mpf(str(k))**2+mp.mpf(str(m2)))
        kk, TT, DD = map(mp.mpf, map(str, [k, T, D]))
        ref = float(mp.log((2*mp.cosh(a*TT)*mp.cosh(kk*DD)+(a/kk+kk/a)*mp.sinh(a*TT)*mp.sinh(kk*DD)-2)
                          /(2*mp.cosh(kk*(TT+DD))-2)))
        error = abs(transfer_logdet(k, T, D, m2)-ref)
        assert error < 3e-12, error
        row = dict(k=k, transfer_high_precision_error=error, reference=ref, nystrom=[])
        for n in [32, 64, 128]:
            g, w = leggauss(n)
            val = fredholm(k, T, [D], (g+1)/2, w/2, np.full(n, m2))
            row['nystrom'].append(dict(nodes=n, value=val, error=abs(val-ref)))
        assert row['nystrom'][-1]['error'] < 3e-6
        assert row['nystrom'][-1]['error'] < row['nystrom'][0]['error']*.15
        # Replicating identical returns changes q but leaves the singlet
        # determinant unchanged. This checks the 1/q normalization directly.
        replicated = fredholm(k, T, [D]*11, (g+1)/2, w/2, np.full(n, m2))
        row['eleven_channel_replication_error'] = abs(replicated-val)
        assert row['eleven_channel_replication_error'] < 2e-12
        rows.append(row)
    return rows


def profile_checks():
    rows = []
    rho = np.linspace(-1, 1, 257)
    for z in [1e-5, .001, .03, .1, .25]:
        b, a = np.sqrt(1-z*z), np.arccos(z)
        ap, app = -z/b, -z/b**3
        p = 1+ap/a+2*rho*ap*np.tan(a*rho)
        pp = app/a-(ap/a)**2+2*rho*app*np.tan(a*rho)+2*rho**2*ap**2/np.cos(a*rho)**2
        np1 = 1+rho*ap*np.tan(a*rho)
        np2 = rho*app*np.tan(a*rho)+rho**2*ap**2/np.cos(a*rho)**2
        assert np.max(abs(p)) <= 3 and np.max(abs(p*p+pp)) <= 16
        assert np.max(abs(np1)) <= 1+1e-10 and np.max(abs(np2)) <= 3
        rows.append(dict(z=z, max_first_relative=float(np.max(abs(p))),
                         max_second_relative=float(np.max(abs(p*p+pp))),
                         max_lapse_log_first=float(np.max(abs(np1))), max_lapse_log_second=float(np.max(abs(np2)))))
    return rows


def energy(phi, n, case, knodes=72, K=40.):
    g, w = leggauss(n)
    y, weights, x = (g+1)/2, w/2, np.exp(phi)
    strength = case['I0']
    if case['kind'] == 'proper':
        e = case['epsilon']
        ac = np.arccos(e*x)
        T = 2*ac/x
        mu2 = strength*e/(2*np.sqrt(1-e*e))
        m2 = mu2*x*x/np.cos(ac*(2*y-1))**2
        Iexact = 2*mu2*np.sqrt(1-e*e*x*x)/e
    else:
        T = np.pi/x
        shape = np.ones(n) if case['kind'] == 'constant' else 1+.6*np.cos(2*np.pi*y)
        m2 = strength/np.pi*x*x*shape
        Iexact = strength*x
    Inum = T*(weights@m2)
    lengths = T+np.asarray(case['returns'])
    a = Inum*np.mean(1/lengths)
    # Integrate the zero mode analytically, so the numerical quadrature sees
    # a finite infrared remainder rather than a logarithmic singularity.
    total = K*np.log1p(a/K**2)+2*np.sqrt(a)*np.arctan(K/np.sqrt(a))
    kg, kw = leggauss(knodes)
    for lo, hi in [(0., 1.), (1., K)]:
        ks, ws = lo+(kg+1)*(hi-lo)/2, kw*(hi-lo)/2
        total += sum(wk*(fredholm(k, T, case['returns'], y, weights, m2)-np.log1p(a/k**2)
                         -Inum/(2*np.sqrt(k*k+1))) for k, wk in zip(ks, ws))
    total /= 2*np.pi
    Lmin = float(np.min(lengths))
    tail = (Inum**2/8/K/np.tanh(K*Lmin/2)**2
            +Inum*np.exp(-K*Lmin)/(K*Lmin*(-np.expm1(-K*Lmin)))
            +Inum/2*np.log((K+np.sqrt(K*K+1))/(2*K)))/(2*np.pi)
    return dict(energy=float(total), integrated_mass=float(Inum), integrated_mass_exact=Iexact,
                integrated_mass_error=abs(Inum-Iexact), continuum_energy_tail_bound=float(tail))


def energy_checks():
    cases = [dict(name='equal_constant', kind='constant', I0=.02, returns=[1.7]),
             dict(name='unequal_constant', kind='constant', I0=.02, returns=[.3, 1.2, 4.7]),
             dict(name='unequal_variable', kind='variable', I0=.02, returns=[.3, 1.2, 4.7]),
             dict(name='proper_epsilon_0.1', kind='proper', epsilon=.1, I0=.02, returns=[.3, 1.2, 4.7]),
             dict(name='proper_epsilon_0.03', kind='proper', epsilon=.03, I0=.004, returns=[.3, 1.2, 4.7])]
    result = []
    for case in cases:
        row = dict(case=case, mesh=[])
        for n in [32, 64, 128]:
            samples = {i: energy(i*.01, n, case) for i in [-4, -2, -1, 0, 1, 2, 4]}
            e = {i: r['energy'] for i, r in samples.items()}
            derivatives = []
            for stride in [2, 1]:
                h = stride*.01
                d1 = (e[-2*stride]-8*e[-stride]+8*e[stride]-e[2*stride])/(12*h)
                d2 = (-e[-2*stride]+16*e[-stride]-30*e[0]+16*e[stride]-e[2*stride])/(12*h*h)
                derivatives.append(dict(step=h, first=d1, second=d2))
            change = max(abs(derivatives[0][key]-derivatives[1][key]) for key in ['first', 'second'])
            assert change < 2e-7, (case, n, change)
            row['mesh'].append(dict(nodes=n, reference=samples[0], derivatives=derivatives,
                                    step_refinement_change=change, max_integrated_mass_error=max(v['integrated_mass_error'] for v in samples.values())))
        def summary(mesh):
            return np.array([mesh['reference']['energy'], mesh['derivatives'][-1]['first'], mesh['derivatives'][-1]['second']])
        changes = [float(np.max(abs(summary(b)-summary(a)))) for a, b in zip(row['mesh'][:-1], row['mesh'][1:])]
        assert changes[-1] < 2e-6 and changes[-1] < changes[0]*.4, (case, changes)
        row['successive_mesh_changes'] = changes
        refined = energy(0, 128, case, knodes=108)
        row['frequency_quadrature_change'] = abs(refined['energy']-row['mesh'][-1]['reference']['energy'])
        assert row['frequency_quadrature_change'] < 2e-8
        if case['name'] == 'equal_constant':
            m2, T, D = case['I0']/np.pi, np.pi, case['returns'][0]
            # Independent transfer determinant, independently adaptive quadrature
            # to infinity (including its actual tail, rather than the tail bound).
            f = lambda k: transfer_logdet(k, T, D, m2)-case['I0']/(2*np.sqrt(k*k+1))
            exact = (quad(f, 0, 1, epsabs=2e-11)[0]+quad(f, 1, np.inf, epsabs=2e-11)[0])/(2*np.pi)
            err = abs(exact-refined['energy'])
            row['independent_transfer_energy'] = exact
            row['transfer_energy_difference'] = err
            assert err < refined['continuum_energy_tail_bound']+2e-7, (err, refined)
            # Compare derivatives over the same finite frequency interval;
            # no derivative tail certificate is inferred from the energy tail.
            def transfer_energy(phi):
                x = np.exp(phi)
                f = lambda k: transfer_logdet(k, np.pi/x, D, m2*x*x)-case['I0']*x/(2*np.sqrt(k*k+1))
                return (quad(f, 0, 1, epsabs=2e-12)[0]+quad(f, 1, 40., epsabs=2e-12)[0])/(2*np.pi)
            h = .01
            te = {i: transfer_energy(i*h) for i in [-2, -1, 0, 1, 2]}
            td = [(te[-2]-8*te[-1]+8*te[1]-te[2])/(12*h),
                  (-te[-2]+16*te[-1]-30*te[0]+16*te[1]-te[2])/(12*h*h)]
            nd = row['mesh'][-1]['derivatives'][-1]
            errors = [abs(td[0]-nd['first']), abs(td[1]-nd['second'])]
            assert max(errors) < 2e-8, errors
            row['independent_transfer_derivatives'] = td
            row['transfer_derivative_errors'] = errors
        result.append(row)
        print(json.dumps(dict(case=case['name'], mesh_changes=changes,
                              energy=row['mesh'][-1]['reference']['energy'])), flush=True)
    return result


def main():
    receipt = dict(status='passed_singlet_static_curvature_checks', date='2026-10-11',
                   exact_checks=exact_checks(), matrix_cases=matrix_checks(), kernel_cases=kernel_checks(),
                   profile_cases=profile_checks(), energy_cases=energy_checks(),
                   gaussian_static_curvature_bound_proved=True, canonical_gaussian_decay_transfer_proved=True,
                   physical_damping_proved=False, nonspherical_interaction_bound_proved=False,
                   physical_state_constructed=False, novelty_established=False, lean_invocations=0,
                   numerical_scope='Finite samples and independent discretizations; uniform constants follow from the ordinary proof. Tail bound applies to continuum energy, not derivative or spatial discretization error.',
                   note_sha256=sha(NOTE), script_sha256=sha(Path(__file__)),
                   dependency_hashes={f: sha(ROOT/f) for f in ['notes/singlet-complex-response.md', 'notes/finite-channel-decay.md', 'notes/schwinger-return-feedback.md']})
    OUT.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'], output=str(OUT))), flush=True)


if __name__ == '__main__':
    main()
