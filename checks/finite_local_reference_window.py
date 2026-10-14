"""Exact finite-local product identities and conditional radius bounds; S074."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import mpmath as mp
import sympy as s


def zero(value):
    assert s.simplify(value) == 0, value


def algebra_checks():
    x, y, magnetic, kap = s.symbols('x y magnetic kappa2', positive=True)
    cs = s.symbols('c1:10')
    c1, c2, c3, c4, c5, c6, c7, c8, c9 = cs
    eta = [-1, 1, 1, 1]
    metric = s.diag(*eta)
    F = s.zeros(4)
    F[2, 3], F[3, 2] = magnetic, -magnetic
    Fup = metric*F*metric

    def riem(a, b, c, d):
        if all(i < 2 for i in [a, b, c, d]):
            curvature = -x
        elif all(i >= 2 for i in [a, b, c, d]):
            curvature = y
        else:
            return s.S.Zero
        return curvature*(metric[a, c]*metric[b, d]-metric[a, d]*metric[b, c])

    ric = s.Matrix(4, 4, lambda a, b: sum(eta[i]*riem(i, a, i, b) for i in range(4)))
    scalar = s.trace(metric*ric)
    ric2 = sum(eta[a]*eta[b]*ric[a, b]**2 for a in range(4) for b in range(4))
    riem2 = sum(eta[a]*eta[b]*eta[c]*eta[d]*riem(a, b, c, d)**2
                for a in range(4) for b in range(4) for c in range(4) for d in range(4))
    inv = sum(F[a, b]*Fup[a, b] for a in range(4) for b in range(4))
    X = F*metric*F.T
    ricX = sum(eta[a]*eta[b]*ric[a, b]*X[a, b] for a in range(4) for b in range(4))
    rFF = sum(riem(a, b, c, d)*Fup[a, b]*Fup[c, d]
              for a in range(4) for b in range(4) for c in range(4) for d in range(4))
    i4 = s.trace((F*Fup)**2)
    actual = [scalar, ric2, riem2, ricX, rFF, i4]
    expected = [2*(y-x), 2*(x*x+y*y), 4*(x*x+y*y), y*inv, 2*y*inv, inv**2/2]
    for left, right in zip(actual, expected):
        zero(left-right)
    tensor_density = (c1*scalar**2+c2*ric2+c3*riem2
                      +kap*(c4*scalar*inv+c5*ricX+c6*rFF)
                      +kap**2*(c7*inv**2+c8*i4))
    I = s.symbols('I', positive=True)
    Ag, K, H = 4*c1+2*c2+4*c3, 2*c4+c5+2*c6, c7+c8/2
    density = Ag*(x*x+y*y)-8*c1*x*y+kap*I*(-2*c4*x+K*y)+kap**2*H*I*I
    zero(tensor_density-density.subs(I, inv))
    pa = Ag*(y*y-x*x)+kap*K*y*I+kap**2*H*I*I
    ps = Ag*(x*x-y*y)+2*kap*c4*x*I-2*kap*K*y*I-3*kap**2*H*I*I
    zero(density-x*s.diff(density, x)-pa)
    zero(density-y*s.diff(density, y)-2*I*s.diff(density, I)-ps)

    # Independent fixed-flux action differentiation, starting with the tensors.
    va, vb, flux = s.symbols('a_squared b_squared flux', positive=True)
    sub = {x: 1/va, y: 1/vb, magnetic: flux/vb}
    volume_action = va*vb*tensor_density.subs(sub)
    target = {x: 1/va, y: 1/vb, I: 2*flux**2/vb**2}
    zero(s.diff(volume_action, va)/vb-pa.subs(target))
    zero(s.diff(volume_action, vb)/va-ps.subs(target))
    zero(2*(pa+ps)+2*kap*I*(-2*c4*x+K*y)+4*kap**2*H*I**2)

    # Directions invisible to these bulk background equations, not to all physics.
    t = s.symbols('t')
    for direction in [{c1:t, c2:-4*t, c3:t}, {c7:-t/8, c8:t/4}, {c9:t}]:
        values = {c: direction.get(c, 0) for c in cs}
        zero(pa.subs(values)); zero(ps.subs(values))

    z, u, alpha, lam, eps, p, heat = s.symbols('z u alpha lambda epsilon A a2', positive=True)
    aa, cc, kk, hh = s.symbols('Ag c4 K H')
    g = aa*(alpha**2-z*z/u**2)+2*alpha*(1+alpha)*kk/z+4*(1+alpha)**2*hh/z**2
    sigma = 4*(1+alpha)*cc/(u*z**2)-2*alpha*(1+alpha)*kk/z**3-8*(1+alpha)**2*hh/z**4
    # Scale x0=1; kappa^2 I=2(1+alpha)/z^2 follows from the classical product.
    normalized = {x:1/u, y:alpha/z, I:2*(1+alpha)/(kap*z**2)}
    compact_pa = pa.subs({c1:(aa-2*c2-4*c3)/4, c7:hh-c8/2, c5:kk-2*c4-2*c6}).subs(c4, cc)
    compact_ps = ps.subs({c1:(aa-2*c2-4*c3)/4, c7:hh-c8/2, c5:kk-2*c4-2*c6}).subs(c4, cc)
    zero(z*z*compact_pa.subs(normalized)-g)
    zero((compact_pa+compact_ps).subs(normalized)-sigma)
    eq1 = -alpha*z-(1-alpha)*z*z/2+(1+alpha)/2-eps*p-lam*g
    eq2 = z*z/u-(1-alpha)*z*z/2-(1+alpha)/2-eps*(-p-heat/2)-lam*(z*z*sigma-g)
    tz = 1-eps*p-lam*g-(1-alpha)*(z-1)**2/2
    denominator = 1-alpha+alpha/z-eps*heat/(2*z*z)+lam*sigma
    zero(eq1-(tz-z))
    zero(eq1+eq2-z*z*(1/u-denominator))
    g0, sigma0 = g.subs({z:1, u:1}), sigma.subs({z:1, u:1})
    dz = -eps*p-lam*g0
    du = eps*(heat/2-alpha*p)-lam*(sigma0+alpha*g0)
    zero(-dz-eps*p-lam*g0)
    zero((1+alpha)*dz-du-eps*(-p-heat/2)-lam*(sigma0-g0))
    return {'product_tensor_contractions': [str(s.simplify(v)) for v in actual],
            'fixed_flux_variations': 'both exact for all eight nonderivative coefficients',
            'trace_scaling_identity': 'PASS',
            'invisible_direction_regressions': ['Euler', 'J^2', 'D.D'],
            'normalized_pressures_and_Einstein_map': 'PASS',
            'first_order_radius_equations': 'PASS'}


def bound_checks():
    lo, hi, cap = s.Rational(99, 100), s.Rational(101, 100), s.Rational(1, 4)
    weights = [8, 32, 4, 128]  # |Ag|, |c4|, |K|, |H|
    tables = {
        'g': [cap**2+(hi/lo)**2, 0, 2*cap*(1+cap)/lo, 4*(1+cap)**2/lo**2],
        'sigma': [0, 4*(1+cap)/lo**3, 2*cap*(1+cap)/lo**3, 8*(1+cap)**2/lo**4],
        'g_derivative_row': [2*hi/lo**2+2*hi**2/lo**3, 0,
                             2*cap*(1+cap)/lo**2, 8*(1+cap)**2/lo**3],
        'sigma_derivative_row': [0, 12*(1+cap)/lo**4,
                                 6*cap*(1+cap)/lo**4, 32*(1+cap)**2/lo**5]}
    for row in tables.values():
        for coefficient, weight in zip(row, weights):
            assert 0 <= coefficient <= weight, (coefficient, weight)
    rho, alpha = s.Rational(1, 100), s.Rational(1, 10**10)
    eB, finite = s.Rational(7, 10**23), s.Rational(1, 10000)
    eps = eB/100  # B>=100; deliberately coarse, uniform bound.
    h = alpha*rho/lo+eps/(2*lo**2)+finite
    assert h < 1-lo
    assert eB+finite+rho**2/2 < rho
    assert h/lo < rho
    assert 1/lo**2 < 2 and 2/lo**3 < 3
    rows = [rho+4*eB+finite, (2*alpha+3*eps+finite)/lo**2]
    assert all(row < s.Rational(2, 100) for row in rows)
    assert 1/(1-rho/2) < 2 and (1+4*alpha)/lo < 2
    assert 2*(eB+finite) < s.Rational(2001, 10**7)
    assert 4/lo**3 < 8  # derivative row of a2(r)/z^2
    # z remainder <=10 D0^2; u remainder below100 D0^2.
    u_remainder = 10*alpha+4*alpha/lo+s.Rational(8,100)+2+(1+4*alpha)**2/lo
    assert u_remainder < 100
    assert 8*10+32+4*5+128*s.Rational(3,2) == 324
    # C>6k^3, alpha<1e-10 imply C/(C-1)<1+1e-9.
    lambda_upper = 48*s.Rational(22,7)*(1+s.Rational(1,10**9))/s.Integer(2000)**10
    assert lambda_upper*324*(2*10**24) < finite
    return {'coefficient_order': ['abs(Ag)', 'abs(c4)', 'abs(K)', 'abs(H)'],
            'weights': weights,
            'rational_coefficient_majorants': {key:list(map(str,row)) for key,row in tables.items()},
            'denominator_change_upper': str(h),
            'contraction_rows_upper': list(map(str, rows)),
            'u_quadratic_remainder_coefficient_upper': str(u_remainder),
            'sufficient_condition': 'lambda*N4<=1e-4 for every k>=2000',
            'squared_radius_bound': 'each <=2*(epsilon*B+lambda*N4)<0.0002001',
            'first_order_remainder': 'each <=100*(epsilon*B+lambda*N4)^2',
            'simple_uniform_coefficient_bound': 'abs(c_i)<=2e24 for i=1,...,8 is sufficient; c9 untested by this background'}


def numerical_checks():
    mp.mp.dps = 90
    family = []
    for kval in [2000, 4000, 20000]:
        k = mp.mpf(kval)
        C = mp.sqrt(1+12*mp.pi*k**6)
        lam = 48*mp.pi*C/(k**10*(C-1))
        family.append({'k':kval, 'lambda':str(lam),
                       'N4_limit':str(1/(10000*lam)),
                       'common_coefficient_limit':str(1/(10000*324*lam))})
    alpha = 1/mp.sqrt(1+12*mp.pi*mp.mpf(2000)**6)
    directions = [[1,0,0,0], [-1,0,0,0], [0,1,0,0], [0,-1,0,0],
                  [0,0,1,0], [0,0,-1,0], [0,0,0,1], [0,0,0,-1],
                  [1,-2,3,-4], [-4,3,-2,1], [1,1,1,1], [-1,-1,-1,-1]]
    records = []
    for direction in directions:
        norm = sum(w*abs(v) for w,v in zip([8,32,4,128], direction))
        for scale in [mp.mpf('0.0001'), mp.mpf('0.00001')]:
            ag, cc, kk, hh = [scale*v/norm for v in direction]  # lambda times coefficients

            def pressures(z,u):
                # Differentiate the homogeneous action itself, not the map's g/sigma.
                def action(av,bv):
                    x, y, ki = 1/av, alpha/bv, 2*(1+alpha)/bv**2
                    # The omitted xy term has identically zero radius derivatives.
                    density = ag*(x*x+y*y)+ki*(-2*cc*x+kk*y)+hh*ki*ki
                    return av*bv*density
                return (mp.diff(lambda av:action(av,z),u)/z,
                        mp.diff(lambda bv:action(u,bv),z)/u)

            def equations(z,u):
                pa, ps = pressures(z,u)
                return (-alpha*z-(1-alpha)*z*z/2+(1+alpha)/2-z*z*pa,
                        z*z/u-(1-alpha)*z*z/2-(1+alpha)/2-z*z*ps)

            z,u = mp.findroot(equations, (1,1), tol=mp.mpf('1e-80'), maxsteps=30)
            pa0,ps0 = pressures(mp.mpf(1),mp.mpf(1))
            dz1,du1 = -pa0, -(pa0+ps0+alpha*pa0)
            residual = max(abs(v) for v in equations(z,u))
            remainder = max(abs(z-1-dz1), abs(u-1-du1))
            assert residual < mp.mpf('1e-70')
            assert max(abs(z-1),abs(u-1)) <= 2*scale
            assert remainder <= 100*scale*scale
            records.append({'direction':direction, 'lambda_N4':str(scale),
                            'sphere_squared_radius_shift':str(z-1),
                            'AdS_squared_radius_shift':str(u-1),
                            'maximum_equation_residual':str(residual),
                            'remainder_over_lambda_N4_squared':str(remainder/(scale*scale))})
    return {'candidate_scales':family, 'local_only_root_cases':records,
            'root_scope':'epsilon=0; independent high-precision local-action variation. Full Dirac extension proved by analytic bounds, not these numerical roots.',
            'maximum_residual':str(max(mp.mpf(r['maximum_equation_residual']) for r in records)),
            'maximum_remainder_ratio':str(max(mp.mpf(r['remainder_over_lambda_N4_squared']) for r in records))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    paths = ['evidence/uniform-dirac-reference.json', 'notes/uniform-dirac-reference.md',
             'notes/matter-aware-local-basis.md', 'notes/finite-local-electromagnetic-input.md']
    deps = {name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in paths}
    prior = json.loads((root/paths[0]).read_text())
    assert prior['status']=='PASS'
    assert prior['analytic']['pressure_and_weighted_derivative_bound']=='100*(1+abs(L)) for q>=1 and 0<r<=1'
    print('Checking product contractions, fixed-flux variation and stationary identities',flush=True)
    result = {'algebra':algebra_checks()}
    print('Checking exact coefficient majorants, uniform contraction and quadratic remainder',flush=True)
    result['analytic_bounds'] = bound_checks()
    print('Checking independent local-action roots in 24 signed coefficient cases',flush=True)
    result['numerical'] = numerical_checks()
    result.update({'status':'PASS',
                   'scope':'conditional homogeneous S0+L4+one-Dirac-loop reference; no physical coefficient matching or full EFT remainder',
                   'dependencies_sha256':deps,
                   'versions':{'python':platform.python_version(),'sympy':s.__version__,'mpmath':mp.__version__},
                   'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('PASS: finite local coefficient window and nonlinear reference-radius bounds',flush=True)


if __name__=='__main__':
    main()
