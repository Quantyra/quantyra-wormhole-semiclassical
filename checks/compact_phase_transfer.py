"""S072 compact electric/rotation phase and RN-AdS capacitance gate.

The angular result is conditional on weak dressing of the KN inertia.
The scale family is not a four-dimensional dynamical error bound.
"""
from pathlib import Path
import hashlib
import json
import warnings
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad, IntegrationWarning
import sympy as s


def symbolic_checks():
    r,b,L,G,g,q,nf,eta,z=s.symbols('r b L G g q nf eta z',positive=True)
    X=b*b/L**2; C=1+6*X
    Q2=b*b*(1+3*X)
    M=(r+r**3/L**2+Q2/r)/(2*G)
    f=1-2*G*M.subs(r,b)/r+Q2/r**2+r**2/L**2
    factor=(r-b)**2*(r*r+2*b*r+L*L+3*b*b)/(L*L*r*r)
    assert s.factor(f-factor)==0
    assert s.simplify((1-b/r)**2/f-L*L/(r*r+2*b*r+L*L+3*b*b))==0
    # Extremal electric-charge curvature at fixed magnetic flux: the
    # extremal-radius variation drops out of dM/d(Q_e^2) at this order.
    assert s.simplify(s.diff(M,r).subs(r,b))==0
    Ec=g*g/(8*s.pi*b)
    phir=b**3/(G*C); N=2*nf*q
    T=N*eta*z/(2*s.pi*phir)
    flux_G=g*g*b*b*(1+3*X)/(s.pi*q*q)
    x0=s.simplify((Ec/T).subs(G,flux_G))
    assert s.simplify(x0-s.pi*q/(8*nf*eta*z*C*(1+3*X)))==0
    # Kerr-Newman-AdS: J=a M_ext and Omega=a(1+X)/b^2 to first order.
    I=b**3*(1+2*X)/(G*(1+X))
    m0=b*(1+2*X)
    m2=(1+X-2*Q2/L**2)/(2*b)
    # At fixed physical magnetic charge, the KN charge parameter is Q*Xi.
    # Extremality makes the O(a^2) radius shift drop out of Delta_r=0.
    mass2=(m2+2*m0/L**2)/G
    assert s.simplify(mass2/(m0/G)**2-1/(2*I))==0
    j=(q-1)/2
    rot=s.simplify(j*(j+1)/(2*I*T))
    gap=(q+1)/C
    assert s.simplify(rot-s.pi*(q-1)*gap*(1+X)/(8*nf*q*eta*z*(1+2*X)))==0
    margin=s.Rational(3,32)*N*N*eta*eta/phir
    assert s.simplify((Ec/margin).subs(G,flux_G)-1/(3*nf*nf*eta*eta*C*(1+3*X)))==0
    return dict(electric_thermal_parameter=str(x0),rotation_thermal_exponent=str(rot),
                first_LL_mass_squared_in_AdS2=str(gap),fixed_charge_curvature='pass')


def D(X):
    return np.arctan(np.sqrt(1+2*X)/(2*np.sqrt(X)))/np.sqrt(X*(1+2*X))


def capacitance(X,kappa,cut):
    nu=(1+np.sqrt(1+4*kappa/(1+6*X)))/2
    def rhs(u,v):
        poly=X+2*X*u+(1+3*X)*u*u
        return [v[1],kappa*v[0]/((1-u)**2*poly)]
    sol=solve_ivp(rhs,(1-cut,0),[1,-nu/cut],method='DOP853',
                  rtol=3e-12,atol=1e-13,max_step=.01)
    assert sol.success
    return -sol.y[1,-1]/sol.y[0,-1]


def electric_kernel_checks():
    mp.mp.dps=70
    rows=[]
    for xx in ('.01','.1','1','5','50','500'):
        x=mp.mpf(xx)
        limit=int(mp.sqrt(170/x))+4
        # Direct normalized TFD amplitudes, versus the theta-function result.
        amps=[mp.exp(-x*n*n/2) for n in range(-limit,limit+1)]
        norm=mp.fsum(v*v for v in amps)
        direct=mp.fsum(a*b for a,b in zip(amps,amps[1:]))/norm
        theta=mp.exp(-x/4)*mp.jtheta(2,0,mp.exp(-x))/mp.jtheta(3,0,mp.exp(-x))
        assert abs(direct-theta)<mp.mpf('1e-55')*theta
        if x>=1:
            assert direct<=2*mp.exp(-x/2)/(1-mp.exp(-2*x))
        if x<=1:
            R=mp.exp(-mp.pi**2/x)
            assert direct>=1-x/4-4*R/(1-R)
        rows.append(dict(x=float(x),coherence=float(direct),log_coherence=float(mp.log(direct)),
                         theta_error=str(abs(direct-theta))))
    return rows


def spin_kernel(a,j,integrate=False):
    twol=np.arange(int(2*(j+np.sqrt(180/a)))+4)
    l=twol/2; d=twol+1; cas=l*(l+1)
    Z=np.sum(d*d*np.exp(-a*cas))
    weights=d*np.exp(-a*cas/2)
    allowed=((np.abs(l[:,None]-l[None,:])<=j)&(l[:,None]+l[None,:]>=j)
             &((twol[:,None]+twol[None,:]+int(2*j))%2==0))
    val=float(np.sum(weights[:,None]*weights[None,:]*allowed)/((2*j+1)*Z))
    assert 0<=val<=1+1e-12
    assert val<=12*np.exp(-a*j*(j+1)/8)*(1+1e-10)+1e-300
    ratio=float(np.sum(d*d*np.exp(-a*cas/2))/Z)
    assert ratio<12
    error=None
    if integrate:
        def integrand(angle):
            sn=np.sin(angle)
            if abs(sn)<1e-14:
                chars=d if angle<np.pi/2 else d*(-1.)**twol
                target=1. if angle<np.pi/2 else (-1.)**int(2*j)
            else:
                chars=np.sin(d*angle)/sn
                target=np.sin((2*j+1)*angle)/((2*j+1)*sn)
            K=np.dot(weights,chars)
            return (2/np.pi)*sn*sn*K*K*target/Z
        direct=quad(integrand,0,np.pi,epsabs=2e-11,epsrel=2e-11,
                    points=[min(np.pi/2,4*np.sqrt(a))],limit=160)[0]
        error=abs(val-direct)
        assert error<2e-9
    return dict(a=a,j=j,coherence=val,heat_trace_ratio=ratio,
                character_quadrature_error=error,
                small_a_approximation=float(np.exp(-a*j*(j+1)/4)))


def family(k):
    mp.mp.dps=65
    k=mp.mpf(k); G=mp.mpf(1); g=k**-4; L=k**5; q=k**3; nf=k**3
    eta=1/k; z=mp.mpf('.5'); N=2*nf*q
    X=(mp.sqrt(1+12*mp.pi*q*q/(g*g*L*L))-1)/6
    C=1+6*X; b=L*mp.sqrt(X); phir=b**3/C
    T=N*eta*z/(2*mp.pi*phir); Ec=g*g/(8*mp.pi*b)
    I=b**3*(1+2*X)/(1+X); j=(q-1)/2
    electric=Ec/T; angular=j*(j+1)/(2*I*T)
    DD=mp.atan(mp.sqrt(1+2*X)/(2*mp.sqrt(X)))/mp.sqrt(X*(1+2*X))
    screening=nf*q*g*g*DD/(4*mp.pi**2)
    llgap=(q+1)/C
    vals=dict(k=k,flavors=nf,magnetic_flux=q,g=g,AdS4_radius_over_lp=L,eta_f=eta,
              radius_squared_over_AdS4_squared=X,electric_thermal_parameter=electric,
              angular_thermal_exponent=angular,spherical_screening_upper=screening,
              first_LL_gap_squared=llgap,Schwinger_mass_squared_AdS2=nf*q*g*g/(4*mp.pi**2*C),
              angular_polarization_scale_proxy=nf*g*g*q,
              species_ratio_N_G_over_b_squared=N/(b*b),area_ratio_N_eta_over_Phi_ext=2*N*eta/(b*b),
              N_eta=N*eta,temperature_times_radius=T*b,
              entropy_curvature_over_electric_thermal_parameter=g*g/(4*C*electric))
    assert electric<mp.mpf('.01') and angular<mp.mpf('.01')
    assert screening<mp.mpf('.01') and llgap>1 and eta<=mp.mpf('0.0005')
    return {key:float(value) for key,value in vals.items()}


if __name__=='__main__':
    warnings.simplefilter('error',IntegrationWarning)
    symbolic=symbolic_checks()
    cap=[]
    for X in (.25,1,100,10000):
        integral=quad(lambda u:1/(X+2*X*u+(1+3*X)*u*u),0,1,epsabs=1e-13)[0]
        assert abs(integral-D(X))<1e-12
        for kappa in (.001,.1,1):
            v=capacitance(X,kappa,1e-7)
            fine=capacitance(X,kappa,1e-9)
            assert 1-1e-10<=fine<=1+kappa*D(X)+1e-10
            assert abs(v-fine)<2e-8
            cap.append(dict(X=X,kappa=kappa,capacitance_ratio=fine,
                            variational_upper=1+kappa*D(X),refinement_error=abs(v-fine)))
    electric=electric_kernel_checks()
    spin=[spin_kernel(a,j,True) for a in (.01,.1,1) for j in (.5,1,5,20)]
    warm=[spin_kernel(a,5) for a in (.01,.001,.0001)]
    assert abs(warm[-1]['coherence']-warm[-1]['small_a_approximation'])<1e-7
    cold_min=np.pi*.99/(16*.0005*.75)
    cold_upper=12*np.exp(-cold_min/(8*1.01))
    assert cold_upper<2e-27
    family_rows=[family(k) for k in (2000,4000,10000)]
    out=dict(date='2026-10-13',status='compact_phase_state_gate_and_spherical_capacitance_verified',
             code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             numpy=np.__version__,scipy=scipy.__version__,sympy=s.__version__,mpmath=mp.__version__,
             symbolic=symbolic,capacitance=cap,electric_kernel=electric,spin_kernel=spin,
             spin_warm_refinement=warm,single_flavor_spin_coherence_upper=cold_upper,
             candidate_many_flavor_family=family_rows,
             scope=['quadratic compact boundary modes at leading RN/KN susceptibility',
                    'static spherical massless-LLL Gaussian polarization with regular horizon',
                    'angular obstruction conditional on inertia dressing at most one percent'],
             not_proved=['full dynamical gauge-invariant switching kernel','matched angular inertia correction',
                         'higher angular sectors and finite-cutoff dynamics','initial-pair preparation',
                         'full four-dimensional stability or physical realization'])
    (Path(__file__).resolve().parents[1]/'evidence/compact-phase-transfer.json').write_text(
        json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: RN/KN compact modes, 12 capacitance cases, 6 electric and 12 spin kernels, scale family')
