"""S070: explicit material and collisionless radiation-wall comparisons.

All accepted identities concern a specified probe matter model. This does
not compute a gravitational instanton determinant or a nucleation rate.
"""
import hashlib
import json
import platform
from pathlib import Path
import mpmath as mp
import sympy as s


def same(a, b):
    assert s.simplify(a-b) == 0, s.simplify(a-b)


def symbolic_checks():
    e, vt, vx, vy = s.symbols('e vt vx vy',real=True)
    rho, b, C = s.symbols('rho b C',positive=True)
    # Embedding-induced metric and comoving scalar density, derived directly.
    eta = s.diag(-1,1,1)
    dz = s.Matrix([vt,vx,vy])
    h = eta+e**2*dz*dz.T
    det = s.factor(h.det())
    same(det,-(1+e**2*(vx*vx+vy*vy-vt*vt)))
    B = h.inv()[1:3,1:3]
    same(B.det(),(1-e**2*vt*vt)/(1+e**2*(vx*vx+vy*vy-vt*vt)))
    lag = -rho*(1-e**2*vt*vt)**s.Rational(3,4)*(1+e**2*(vx*vx+vy*vy-vt*vt))**s.Rational(-1,4)
    quadratic = s.expand(s.diff(lag,e,2).subs(e,0)/2)
    same(quadratic,rho*vt*vt/2+rho*(vx*vx+vy*vy)/4)
    density = C*b**s.Rational(3,2)
    pressure = b*s.diff(density,b)-density
    same(pressure,density/2)
    same(s.diff(pressure,b)/s.diff(density,b),s.Rational(1,2))

    # Static sphere: positive kinetic term, negative angular gradient term.
    R,f0,f2,x,v,q,J = s.symbols('R f0 f2 x v q J',positive=True)
    ff = f0+2*f0*x/R+f2*x*x/2
    u = s.sqrt(ff)/(R+x)
    same(s.diff(u,x).subs(x,0),0)
    u2 = s.simplify(s.diff(u,x,2).subs(x,0))
    radial_growth = f0*f0/R**2-f0*f2/2
    same(-R*f0**s.Rational(3,2)*u2,radial_growth)
    X = f0+e*e*v*v/f0
    Y = e*e*q/(R*R)
    exact_euclidean = C/R*X**s.Rational(3,4)*(X+Y)**s.Rational(-1,4)
    same(s.diff(exact_euclidean,e,2).subs(e,0)/2,
         C*v*v/(2*R*f0**s.Rational(3,2))-C*q/(4*R**3*s.sqrt(f0)))
    # A tangential compression J only mixes through u'(R), which vanishes.
    same(s.diff(u*J**s.Rational(3,2),x,J).subs(x,0),0)

    # Fixed angular-momentum Routh reduction of a massless particle.
    z = s.symbols('z',real=True)
    pr = J*z/(R*s.sqrt(f0)*s.sqrt(1-z*z))
    H = J*s.sqrt(f0)/(R*s.sqrt(1-z*z))
    same(H*H,f0*f0*pr*pr+f0*J*J/R**2)
    same(pr*(f0*z)-H,-J*s.sqrt(f0)*s.sqrt(1-z*z)/R)
    # Independent Hamiltonian Hessian yields the same radial instability.
    p = s.symbols('p',real=True)
    hh = s.sqrt(ff*ff*p*p+ff*J*J/(R+x)**2)
    aa = s.diff(hh,p,2).subs({x:0,p:0})
    bb = s.diff(hh,x,2).subs({x:0,p:0})
    same(-aa*bb,radial_growth)

    # Explicit extremal RN-AdS geometry. a is the extremal horizon radius;
    # zeta=a^2/L_AdS^2. Derive from the general metric, not printed F.8.
    r,a,zeta = s.symbols('r a zeta',positive=True)
    fr = 1+zeta*r*r/a**2-2*a*(1+2*zeta)/r+a*a*(1+3*zeta)/r**2
    rr = 2*a*(1+3*zeta)
    same((r*s.diff(fr,r)-2*fr).subs(r,rr),0)
    fR = s.factor(fr.subs(r,rr))
    growth = s.factor((fr*fr/r**2-fr*s.diff(fr,r,2)/2).subs(r,rr))
    expected = (1+2*zeta)*(1+6*zeta)**4/(32*a*a*(1+3*zeta)**4)
    same(growth,expected)
    same(growth.subs(zeta,0),1/(32*a*a))
    # The retrieved paper's F.8 differs by a factor R and by dimensions.
    printed_f8 = (1+2*zeta)*(1+6*zeta)**4/(64*a**3*(1+3*zeta)**5)
    same(growth/printed_f8,rr)

    # Local fixed-energy thermodynamics does not reverse the area term.
    entropy,area,energy,k = s.symbols('entropy area energy k',positive=True)
    E = k*entropy**s.Rational(3,2)/s.sqrt(area)
    temp = s.diff(E,entropy)
    press = -s.diff(E,area)
    same(press,E/(2*area))
    S_at_E = (energy*s.sqrt(area)/k)**s.Rational(2,3)
    same(s.diff(S_at_E,area).subs(energy,E),press/temp)
    return {'flat_normal_L2':str(quadratic),'surface_pressure_ratio':'1/2',
            'sound_speed_squared':'1/2','RNAdS_static_radius':str(rr),
            'RNAdS_f_at_wall':str(fR),'RNAdS_radial_growth_squared':str(growth),
            'printed_F8_discrepancy_factor':str(rr)}


def finite_difference_checks():
    """Integrate the nonlinear material action, independently of its Hessian."""
    mp.mp.dps = 70
    rows = []
    for zeta in [mp.mpf(0),mp.mpf('0.1')]:
        R = 2*(1+3*zeta)
        def ff(r):
            return 1+zeta*r*r-2*(1+2*zeta)/r+(1+3*zeta)/(r*r)
        u = lambda r: mp.sqrt(ff(r))/r
        u2 = mp.diff(u,R,2)
        for ell in (0,2,5):
            def action(eps):
                def integrand(z):
                    P = mp.legendre(ell,z)
                    dP = mp.diff(lambda w:mp.legendre(ell,w),z)
                    rr = R+eps*P
                    grad2 = eps*eps*(1-z*z)*dP*dP
                    return u(rr)*(1+grad2/(ff(rr)*rr*rr))**(-mp.mpf(1)/4)
                return 2*mp.pi*mp.quad(integrand,[-1,0,1])
            predicted = 4*mp.pi/(2*ell+1)*(u2-ell*(ell+1)/(2*R**3*mp.sqrt(ff(R))))
            errors = []
            for eps in [mp.mpf('1e-3'),mp.mpf('1e-4'),mp.mpf('1e-5')]:
                second = (action(eps)+action(-eps)-2*action(0))/(eps*eps)
                errors.append(float(abs((second-predicted)/predicted)))
                assert second < 0
            assert errors[-1] < 1e-8
            assert errors[-1] < errors[0]/9000
            rows.append({'a':1,'a2_over_L2':float(zeta),'ell':ell,
                         'predicted_action_second_derivative':float(predicted),
                         'relative_errors_h_1e_3_1e_4_1e_5':errors})
    return rows


def comparison():
    # Extremal asymptotically flat RN, a=1, R=2, f(R)=1/4.
    # The collisionless norm-growth exponent follows from Hamiltonian flow.
    out = []
    for ell in (0,2,5,10,50):
        material_growth2 = mp.mpf(1)/32+mp.mpf(1)/32*ell*(ell+1)
        out.append({'ell':ell,'material_coordinate_time_growth':float(mp.sqrt(material_growth2)),
                    'collisionless_phase_space_norm_growth':float(mp.sqrt(mp.mpf(1)/32))})
    # Only a matter-model comparison; no microscopic cutoff is assigned.
    return out


def main():
    result = {'date':'2026-10-12','status':'probe_completion_comparison_verified',
              'python':platform.python_version(),'sympy':s.__version__,'mpmath':mp.__version__,
              'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'symbolic':symbolic_checks(),'nonlinear_action_checks':finite_difference_checks(),
              'growth_comparison':comparison(),
              'not_computed':['full constrained gravitational Hessian','microscopic radiation wall state',
                              'collisionless Euclidean determinant','instanton negative-mode count',
                              'wormhole nucleation rate','wormhole stability']}
    path = Path('evidence/ads-wall-response.json')
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: fluid action, radial Routhian equality, geodesic Hessian, RN-AdS formula, thermodynamic constraint and 18 nonlinear-action checks.')
    print('No gravitational instanton mode count or formation rate accepted.')


if __name__=='__main__':
    main()
