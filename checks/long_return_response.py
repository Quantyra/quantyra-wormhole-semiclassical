"""Independent checks of the prescribed free-channel long-return response.

Run with Python and the repository's pinned numpy/scipy/sympy dependencies.
No Lean, cloud, or four-dimensional state construction is performed.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import roots_legendre
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
PI = np.pi
KAPPA = 1 + np.arcsinh(np.sqrt(3)) / (2 * np.sqrt(3))
T, H = PI, 0.04
WLO, WHI = 0.6, 1.4


def weight(n, length, throat=T, collar=H):
    omega = 2 * PI * n / length
    return (2 / (3 * length**2) * (n*n-1) / n
            * np.sin(omega*throat/2)**2 * np.exp(-(collar*omega)**2))


def window(omega):
    z = (2*omega-WLO-WHI)/(WHI-WLO)
    return np.maximum(1-z*z, 0)**4 * (np.abs(z) < 1)


def dipole_scaled_length(u):
    """Return u*f(u), evaluated without singular 1/u factors."""
    def integrand(eta):
        aa = 1-eta*eta
        qq = np.hypot(aa, 2*u*eta)
        umu = (aa+qq)/2
        return np.sqrt(umu + 4*umu*eta*eta*(umu-u*u)/(qq*qq))
    # Resolve the eta=1 layer, whose width is O(u). An unsplit adaptive
    # rule can miss its O(u) contribution while reporting tiny error.
    points=sorted({1-scale*u for scale in [1,10,100,1000] if 0<1-scale*u<1})
    return quad(integrand, 0, 1, points=points,
                epsabs=2e-11, epsrel=2e-11, limit=150)[0]


def dipole_length_from_coordinates(u):
    """Independent arclength: differentiate z=mu*eta/2, r^2=(mu^2-1)(1-eta^2)/4."""
    def integrand(eta):
        aa = 1-eta*eta
        qq = np.hypot(aa, 2*u*eta)
        mu = (aa+qq)/(2*u)
        mup = -2*eta*(mu-u)/qq
        zprime = (mu+eta*mup)/2
        rr = (mu*mu-1)*aa/4
        rrprime = (2*mu*mup*aa-2*eta*(mu*mu-1))/4
        return 2*np.sqrt(zprime*zprime+rrprime*rrprime/(4*rr))
    points=sorted({1-scale*u for scale in [1,10,100,1000] if 0<1-scale*u<1})
    return quad(integrand, 0, 1, points=points,
                epsabs=1e-9, epsrel=2e-11, limit=150)[0]


def coefficient_density(omega, nmax=4000):
    """rho/(q*d) in the inverse-tail model; no frequency delta broadening."""
    n = np.arange(2, nmax+1, dtype=float)
    value = np.sum((1-1/(n*n))/(2*PI*n/omega-T)**2)
    return KAPPA/3*np.sin(omega*T/2)**2*np.exp(-(H*omega)**2)*value


def frequency_integral(nmax=4000):
    value = quad(lambda w: window(w)*coefficient_density(w,nmax), WLO,WHI,
                 epsabs=3e-12, epsrel=3e-11)[0]
    # For x=omega*T/(2*pi)<1, sum_{n>N}(n-x)^-2 <= 1/(N-x).
    tail = quad(lambda w: window(w)*KAPPA*w*w/(12*PI*PI)
                *np.sin(w*T/2)**2*np.exp(-(H*w)**2)/(nmax-w*T/(2*PI)),
                WLO,WHI,epsabs=1e-14)[0]
    return value,tail


def loop_integral(d, nodes, exact_dipole):
    """Integrate pi*sum W_n*g(omega_n) over u directly after u=d*v.

    The rigorous bound u*f(u)<=sqrt(5) limits the contributing v interval.
    This path never uses a spectral Jacobian or a derivative of f.
    """
    vmax = np.sqrt(5)/(4*PI/WHI-T)
    x,w = roots_legendre(nodes)
    values=[]
    for v in (x+1)*vmax/2:
        delay=(dipole_scaled_length(d*v) if exact_dipole else KAPPA)/v
        length=T+delay
        n=np.arange(max(2,int(np.ceil(WLO*length/(2*PI)))),
                    int(np.floor(WHI*length/(2*PI)))+1,dtype=float)
        values.append(PI*np.sum(weight(n,length)*window(2*PI*n/length)))
    return float(np.dot(w,values)*vmax/2)


def discrete_integral(d,q):
    """Midpoint Riemann sum over independent inverse-tail channels."""
    vmax = KAPPA/(4*PI/WHI-T)
    count=min(q,int(np.floor(q*d*vmax+0.5)))
    total=0.0
    for j in range(count):
        u=(j+0.5)/q
        length=T+KAPPA*d/u
        n=np.arange(max(2,int(np.ceil(WLO*length/(2*PI)))),
                    int(np.floor(WHI*length/(2*PI)))+1,dtype=float)
        total+=PI*np.sum(weight(n,length)*window(2*PI*n/length))
    return total/(q*d)


def main():
    n,j=sy.symbols('n j',integer=True,positive=True)
    length,throat=sy.symbols('L T',positive=True)
    exact={}
    def check(name,expression):
        residual=sy.simplify(expression)
        assert residual==0,(name,residual)
        exact[name]=str(residual)
    check('one_chirality_particle_hole_sum',
          sy.summation((j-sy.Rational(1,2)-n/2)**2,(j,1,n))-n*(n*n-1)/12)
    rational=(n*n-1)/(n*n*(2*n-1)**2)
    check('resonance_sum_partial_fractions',rational-(8/(2*n-1)-3/(2*n-1)**2-4/n-1/n**2))
    check('resonance_sum_closed_form',
          4*(sy.digamma(2)-sy.digamma(sy.Rational(3,2)))
          -3*(sy.pi**2/8-1)-(sy.pi**2/6-1)-(8*sy.log(2)-13*sy.pi**2/24))
    check('low_frequency_coefficient',
          (sy.zeta(2)-sy.zeta(4))/(48*sy.pi**2)-(15-sy.pi**2)/4320)
    omega,dd,kk=sy.symbols('omega d kappa',positive=True)
    root=kk*dd/(2*sy.pi*n/omega-throat)
    check('inverse_tail_spectral_jacobian',sy.diff(root,omega)
          -kk*dd*(2*sy.pi*n/omega)**2/(2*sy.pi*n*(2*sy.pi*n/omega-throat)**2))
    u,mu,e=sy.symbols('u mu e',positive=True)
    aa=1-e
    qq=2*u*mu-aa
    coordinate_speed=(mu-2*e*(mu-u)/qq)**2+e*(-2*mu*(mu-u)*aa/qq-(mu*mu-1))**2/((mu*mu-1)*aa)
    compact_speed=mu/u+4*mu*e*(mu-u)/qq**2
    check('dipole_arclength_identity_on_field_line',
          sy.factor((coordinate_speed-compact_speed).subs(e,mu*(1-u*mu)/(mu-u))))
    x=sy.symbols('x',nonnegative=True)
    check('dipole_tail_integral',sy.integrate(sy.sqrt(1+3*x*x),(x,0,1))
          -(1+sy.asinh(sy.sqrt(3))/(2*sy.sqrt(3))))
    # Countercheck the normalization directly with complex Fourier integrals.
    modes=[]
    for ll,tt,hh in [(5.0,3.0,0.04),(9.0,2.1,0.2),(4.2,3.9,0.01)]:
        for nn in [1,2,3,9,17]:
            wave=2*PI*nn/ll
            f=(quad(lambda x:np.cos(wave*x),0,tt)[0]
               -1j*quad(lambda x:np.sin(wave*x),0,tt)[0])/ll
            f*=np.exp(-(hh*wave)**2/2)
            filled=np.arange(nn,dtype=float)-nn+0.5
            direct=2*np.sum(np.abs((2*PI/ll)*(filled+nn/2)*f)**2)
            closed=weight(nn,ll,tt,hh)
            error=abs(direct-closed)
            assert error<2e-13
            modes.append(dict(L=ll,T=tt,h=hh,n=nn,direct=float(direct),closed=float(closed),error=float(error)))
    # Uniform lapse gives only n=0, and the n=1 vacuum weight vanishes.
    assert weight(1,5.0)==0
    geometry=[]
    for u in [0.9,0.4,0.1,0.02,0.002,0.0002,0.000002]:
        scaled=dipole_scaled_length(u)
        other=u*dipole_length_from_coordinates(u)
        assert abs(scaled-other)<2e-9
        assert 0<scaled<=np.sqrt(5)
        geometry.append(dict(u=u,u_f=scaled,coordinate_check=other,error=abs(scaled-other),tail_error=abs(scaled-KAPPA)))
    assert geometry[-1]['tail_error']<1e-5
    freq,tail=frequency_integral()
    inverse=[dict(nodes=k,value=loop_integral(0.1,k,False)) for k in [240,480]]
    # The missing mode tail in the frequency route has an analytic upper bound.
    assert freq-3e-9<=inverse[-1]['value']<=freq+tail+3e-9
    assert abs(inverse[-1]['value']-inverse[0]['value'])<3e-9
    dipole=[]
    for d in [0.4,0.2,0.1,0.05,0.02,0.005]:
        coarse,fine=(loop_integral(d,k,True) for k in [240,480])
        assert abs(fine-coarse)<5e-9
        dipole.append(dict(d=d,value=fine,quadrature_change=abs(fine-coarse),difference_from_tail=fine-inverse[-1]['value']))
    assert abs(dipole[-1]['difference_from_tail'])<abs(dipole[0]['difference_from_tail'])/20
    discrete=[]
    for q in [1000,4000,16000]:
        value=discrete_integral(0.1,q)
        discrete.append(dict(q=q,d=0.1,value=value,error=abs(value-inverse[-1]['value'])))
    assert discrete[-1]['error']<discrete[0]['error']/20
    # Low-frequency and resonance checks use extrapolated tail-bounded sums.
    resonance_exact=8*np.log(2)-13*PI*PI/24
    modes_n=np.arange(2,100001,dtype=float)
    resonance_sum=float(np.sum((modes_n*modes_n-1)/(modes_n*modes_n*(2*modes_n-1)**2)))
    assert 0<resonance_exact-resonance_sum<1/(4*(100000-0.5))
    low=[]
    for w in [0.1,0.02,0.004]:
        ratio=coefficient_density(w,100000)/(KAPPA*T*T*w**4*(15-PI*PI)/4320)
        low.append(dict(omega=w,ratio=ratio))
    assert abs(low[-1]['ratio']-1)<0.01
    report=dict(status='passed_prescribed_long_return_response',
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sy.__version__),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        note_sha256=hashlib.sha256((ROOT/'notes/long-return-response.md').read_bytes()).hexdigest(),
        exact_checks=exact,fermion_mode_cases=modes,dipole_geometry_cases=geometry,
        weak_spectral_test=dict(T=T,h=H,omega_window=[WLO,WHI],
            frequency_integral_partial=freq,analytic_omitted_mode_bound=tail,
            inverse_tail_direct_loop_integrals=inverse,dipole_limits=dipole,finite_channel_limits=discrete),
        constants=dict(kappa=KAPPA,resonance_sum=resonance_exact,
            conditional_gamma_over_omega_coefficient=4*KAPPA*resonance_exact/(3*PI*PI)),
        low_frequency_checks=low,model='Independent NS Dirac circles, smooth prescribed relative lapse; continuum point-dipole length distribution',
        four_dimensional_quantum_state_constructed=False,physical_collective_mode_coupling_verified=False,
        physical_damping_proved=False,novelty_established=False,lean_invocations=0)
    (ROOT/'evidence/long-return-response.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=report['status'],exact_checks=len(exact),mode_cases=len(modes),
          geometry_cases=len(geometry),frequency_integral=freq,mode_tail_bound=tail,
          direct_loop_integral=inverse[-1]['value'],dipole_last=dipole[-1],
          finite_q_last=discrete[-1],constants=report['constants']),indent=2))


if __name__=='__main__':
    main()
