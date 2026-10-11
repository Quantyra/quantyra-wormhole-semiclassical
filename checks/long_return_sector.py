"""Independent geometry/operator checks for a local shrinking-sector theorem.

No physical full-return, state or damping certificate is produced.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.special import roots_legendre
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.])


def main(output):
    R, Z, v, eta = sp.symbols('R Z v eta', real=True)
    # Integrate derivatives of the nonsingular representation, without
    # subtracting nearly equal monopole fluxes at small u.
    integrand = R**2 / (2*(R**2+(Z+v)**2)**sp.Rational(3, 2))
    derivatives = [integrand, sp.diff(integrand, R), sp.diff(integrand, Z),
                   sp.diff(integrand, R, 2), sp.diff(integrand, R, Z),
                   sp.diff(integrand, Z, 2)]
    calc = sp.lambdify((R, Z, v), derivatives, 'numpy', cse=True)
    nodes, weights = roots_legendre(24)

    def field(r, z, u):
        values = np.array(calc(r, z, u*nodes/2)) @ (weights/2)
        flux, fr, fz, frr, frz, fzz = values
        bvec = np.array([fz, -fr])/r
        jac = np.array([[frz/r-fz/r**2, fzz/r],
                        [-frr/r+fr/r**2, -frz/r]])
        mass = (flux-1)/r
        hessian = np.array([[frr/r-2*fr/r**2+2*(flux-1)/r**3,
                             frz/r-fz/r**2],
                            [frz/r-fz/r**2, fzz/r]])
        return mass, bvec, jac, hessian

    def mass_only(r, z, u):
        # This vectorized evaluation retains the exact rescaled dipole.
        r = np.asarray(r)
        z = np.asarray(z)
        flux = np.sum(r[..., None]**2 /
                      (2*(r[..., None]**2+(z[..., None]+u*nodes/2)**2)**1.5)
                      * weights/2, axis=-1)
        return (flux-1)/r

    a = 1-eta**2
    b = (a+sp.sqrt(a*a+4*v*v*eta*eta))/2
    cr = sp.sqrt((b*b-v*v)*(1-eta*eta))/2
    cz = b*eta/2
    curve = sp.lambdify((v, eta), [cr, cz, sp.diff(cr, eta), sp.diff(cz, eta)],
                       'numpy', cse=True)

    def basic(u, e):
        r, z, dr, dz = np.array(curve(u, e), float)
        mass, bvec, jac, hess = field(r, z, u)
        strength = np.linalg.norm(bvec)
        tangent = bvec/strength
        normal = np.array([-tangent[1], tangent[0]])
        speed = np.hypot(dr, dz)
        assert np.linalg.norm(tangent-np.array([dr, dz])/speed) < 2e-12
        assert abs(mass) < 2e-12
        bs = float(tangent @ jac @ tangent)
        return dict(r=r, z=z, b=strength, a=bs/(4*strength),
                    c=float(normal @ hess @ normal),
                    k=float(normal @ jac @ tangent/strength),
                    speed=speed, normal=normal)

    def geometry(u, e, step):
        g = basic(u, e)
        nearby = [basic(u, e+j*step) for j in [-2, -1, 1, 2]]
        for key in ['a', 'c', 'k']:
            g[key+'s'] = (nearby[0][key]-8*nearby[1][key]
                         +8*nearby[2][key]-nearby[3][key])/(12*step*g['speed'])
        # Same fixed rescaled tube width in all numerical cases.
        g['width'] = 0.012
        return g

    exact = {}
    def check(name, expr):
        assert sp.simplify(expr) == 0, name
        exact[name] = '0'
    f = Z/sp.sqrt(R*R+Z*Z)
    check('flux_integral_derivative', sp.diff(f, Z)/2-integrand.subs(v, 0))
    limit_flux = R**2/(2*(R*R+Z*Z)**sp.Rational(3, 2))
    # On the curve, r^2+z^2=(1-eta^2)^2/4. Use positive a explicitly.
    aa = sp.symbols('a', positive=True)
    check('limit_curve_level', (aa**3/4)/(2*(aa**2/4)**sp.Rational(3, 2))-1)
    q, u, r, zz = sp.symbols('q u r zz', positive=True)
    vu = sp.symbols('V')
    check('mass_dilation_parameter', (q*(u*vu-u)/(r/u))/u-q*u*(vu-1)/r)
    n = sp.symbols('n', positive=True)
    sn = (n*n-1)/(n*n*(2*n-1)**2)
    check('spectral_partial_fraction', sn-(8/(2*n-1)-3/(2*n-1)**2-4/n-1/n**2))
    hh, pp, jj, fp = sp.symbols('h P J Fprime')
    vec = sp.Matrix([pp, jj])
    advective = -sp.I*fp*(sp.diag(1,-1)*vec/hh-vec)
    check('moving_envelope_upper', advective[0]+sp.I*fp*(1/hh-1)*pp)
    check('moving_envelope_lower', advective[1]-sp.I*fp*(1/hh+1)*jj)

    flux_cases = []
    for u0 in [.1, .01, .001]:
        for e in [-.7, -.35, 0., .35, .7]:
            r0, z0, _, _ = curve(u0, e)
            # Direct monopole expression is independent of integral evaluation.
            direct = ((z0+u0/2)/np.hypot(r0, z0+u0/2)
                      -(z0-u0/2)/np.hypot(r0, z0-u0/2))/(2*u0)
            value = mass_only(r0, z0, u0)*r0+1
            error = abs(value-direct)
            assert error < 2e-13
            rlim = (1-e*e)**1.5/2
            zlim = e*(1-e*e)/2
            lim = r0*r0/(2*(r0*r0+z0*z0)**1.5)
            flux_cases.append(dict(u=u0, eta=e, flux_error=float(error),
                scaled_curve_error=float(np.hypot(r0-rlim, z0-zlim)/u0**2),
                scaled_flux_error=float(abs(value-lim)/u0**2)))

    def residual(g, u0, Q, wave, order):
        # Resolve the C1 cutoff at its actual breakpoints. Global Hermite
        # quadrature converges poorly when nodes straddle that transition.
        scale = np.sqrt(Q*g['b'])
        edge = min(g['width']*scale, 12.)
        joins = sorted({-edge, edge, 0.,
                        max(-edge, -g['width']*scale/2),
                        min(edge, g['width']*scale/2)})
        nodes0, weights0 = roots_legendre(order)
        yn, wn = [], []
        for left, right in zip(joins[:-1], joins[1:]):
            yy = (left+right)/2+(right-left)*nodes0/2
            yn.extend(yy)
            wn.extend(weights0*(right-left)/2*np.exp(-yy**2)/np.sqrt(np.pi))
        y, wt = np.array(yn)/scale, np.array(wn)
        b, c, a, k = [g[key] for key in ['b', 'c', 'a', 'k']]
        h = 1-k*y
        hs = -g['ks']*y
        assert np.min(h) > .85
        m = mass_only(g['r']+g['normal'][0]*y,
                      g['z']+g['normal'][1]*y, u0)
        P = 1-Q*c*y**3/6
        J = a*y
        v0 = np.array([P, J])
        G = a*(1-2*Q*b*y*y)
        dn = np.array([-Q*c*y*y/2-Q*b*y*P, a-Q*b*y*J])
        ds = np.array([-Q*g['cs']*y**3/6+P*(G+1j*wave),
                       g['as']*y+J*(G+1j*wave)])
        # Direct first-order Pauli operator, without the canceled formula.
        direct = (1j*SX@dn-1j*SZ@ds/h+1j*SZ@v0*hs/(2*h*h)
                  +SY@v0*(Q*m)-wave*v0)
        Xi = -1j*wave
        dm = m-b*y
        r3 = dm-c*y*y/2
        upper = (G*(-k*y+Q*c*y**3/6)/h+Q*g['cs']*y**3/(6*h)
                 +Xi*(1/h-1)*P+hs*P/(2*h*h)-Q*a*y*dm)
        lower = (Q*r3-Q*Q*c*y**3*dm/6+g['as']*y/h+J*G/h
                 -Xi*J*(1/h+1)-hs*J/(2*h*h))
        mismatch = float(np.max(abs(direct-1j*np.array([upper, lower]))))
        assert mismatch < 2e-8
        chi = np.ones_like(y)
        chp = np.zeros_like(y)
        transition = abs(y) > g['width']/2
        ang = np.pi*(abs(y[transition])/g['width']-.5)
        chi[transition] = np.cos(ang)**2
        chp[transition] = -np.pi/g['width']*np.sin(2*ang)*np.sign(y[transition])
        residual = chi*direct+1j*SX@v0*chp
        norm = float(np.sqrt(np.sum(wt*np.sum(abs(residual)**2, axis=0))))
        packet_norm = float(np.sum(wt*chi**2*np.sum(abs(v0)**2, axis=0)))
        # An uncanceled leading Gaussian checks the need for both correctors.
        raw = np.array([np.ones_like(y), np.zeros_like(y)])
        raw_dn = np.array([-Q*b*y, np.zeros_like(y)])
        raw_ds = np.array([G+1j*wave, np.zeros_like(y)])
        raw_res = (1j*SX@raw_dn-1j*SZ@raw_ds/h+1j*SZ@raw*hs/(2*h*h)
                   +SY@raw*(Q*m)-wave*raw)
        raw_res = chi*raw_res+1j*SX@raw*chp
        return norm, packet_norm, mismatch, float(np.sqrt(np.sum(wt*np.sum(abs(raw_res)**2, axis=0))))

    sections = []
    for u0 in [.01, .001, .0001]:
        for e in [-.6, 0., .6]:
            coarse = geometry(u0, e, 2e-4)
            g = geometry(u0, e, 1e-4)
            refinement = max(abs(g[k]-coarse[k])/max(1, abs(g[k]))
                             for k in ['as', 'cs', 'ks'])
            assert refinement < 1e-7
            for wave in [1., 4.]:
                rows = []
                for Q in [10000.5, 40000.5, 160000.5, 640000.5]:
                    # u=10^-k and Q half-integral give integer q=Q/u and
                    # physical half-integral azimuthal label q-Q.
                    q0 = round(Q/u0)
                    assert abs(q0*u0-Q) < 2e-9
                    assert (q0-Q) % 1 == .5
                    val = residual(g, u0, Q, wave, 96)
                    fine = residual(g, u0, Q, wave, 160)
                    err = abs(val[0]-fine[0])
                    assert err < 2e-7*max(1, fine[0]), (u0, e, Q, err)
                    assert fine[0] < fine[3]
                    rows.append(dict(q=q0, Q=Q, residual=fine[0],
                        sqrt_Q_residual=float(np.sqrt(Q)*fine[0]),
                        packet_norm_squared=fine[1], operator_identity_error=fine[2],
                        leading_residual=fine[3], quadrature_error=err))
                ratio = rows[-1]['sqrt_Q_residual']/rows[-2]['sqrt_Q_residual']
                assert ratio < 1.02
                sections.append(dict(u=u0, eta=e, wave=wave,
                    derivative_refinement=refinement, b=g['b'],
                    curvature=g['k'], rows=rows, last_scaled_ratio=ratio))

    tails = []
    total = 8*np.log(2)-13*np.pi**2/24
    for N in [2, 8, 32, 128, 512]:
        ns = np.arange(2, N+1, dtype=float)
        partial = np.sum((ns*ns-1)/(ns*ns*(2*ns-1)**2))
        omitted = total-partial
        bound = 1/(2*(2*N-1))
        assert 0 < omitted < bound
        tails.append(dict(N=N, omitted=float(omitted), bound=bound,
                          relative_bound=bound/total))
    nu = F(15, 16)
    exponents = dict(discrete_resolution=1-2*nu,
                    conditional_residual_accumulation=1-F(3, 2)*nu,
                    frozen_rotation=2-F(5, 2)*nu,
                    nominal_overlap=F(1, 2)-nu, radiation_lifetime=2-4*nu)
    assert all(x < 0 for x in exponents.values())
    expected = [F(-7, 8), F(-13, 32), F(-11, 32), F(-7, 16), F(-7, 4)]
    assert list(exponents.values()) == expected
    out = dict(status='passed_local_shrinking_sector_checks',
        exact_checks=exact, flux_cases=flux_cases, section_cases=sections,
        spectral_tail_cases=tails, rational_exponents={k:str(v) for k,v in exponents.items()},
        local_uniform_packet_argument=True, global_return_proved=False,
        four_dimensional_state_constructed=False, physical_damping_proved=False,
        novelty_established=False, lean_invocations=0,
        versions=dict(numpy=np.__version__, scipy=scipy.__version__, sympy=sp.__version__),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        note_sha256=hashlib.sha256((ROOT/'notes/long-return-sector.md').read_bytes()).hexdigest())
    Path(output).write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=out['status'], flux_cases=len(flux_cases),
                         section_cases=len(sections), exponents=out['rational_exponents']), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    main(parser.parse_args().output)
