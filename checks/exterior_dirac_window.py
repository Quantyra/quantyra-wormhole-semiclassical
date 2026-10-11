"""Exact exterior Weyl reduction and norm-stable finite-window comparison."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
import sympy as sp


sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex)
sz=np.diag([1.,-1.])


def potential(r,z,q,a=1.):
    rm=np.hypot(r,z+a); rp=np.hypot(r,z-a)
    return q*(1-((z+a)/rm-(z-a)/rp)/2)


def field(r,z,q,a=1.):
    rm=np.hypot(r,z+a); rp=np.hypot(r,z-a)
    return q/2*np.array([r/rm**3-r/rp**3,(z+a)/rm**3-(z-a)/rp**3])


def verify():
    r,z,ph,ell=sp.symbols('r z ph ell',real=True)
    # Use exact Pauli matrices, without floating symbolic coefficients.
    SX=sp.Matrix([[0,1],[1,0]]); SY=sp.Matrix([[0,-sp.I],[sp.I,0]]); SZ=sp.diag(1,-1)
    U=sp.diag(sp.exp(sp.I*(ell-sp.Rational(1,2))*ph),sp.exp(sp.I*(ell+sp.Rational(1,2))*ph))
    SR=sp.Matrix([[0,sp.exp(-sp.I*ph)],[sp.exp(sp.I*ph),0]])
    SP=sp.Matrix([[0,-sp.I*sp.exp(-sp.I*ph)],[sp.I*sp.exp(sp.I*ph),0]])
    M=sp.Function('M')(r,z); f=sp.Function('f')(r,z); g=sp.Function('g')(r,z)
    vector=sp.Matrix([f,g])
    def H(v):return -sp.I*SX*v.diff(r)-sp.I*SZ*v.diff(z)+SY*M*v
    square=H(H(vector))-(-vector.diff(r,2)-vector.diff(z,2)+M*M*vector+(SZ*sp.diff(M,r)-SX*sp.diff(M,z))*vector)
    identities={
        'radial_spin_frame':U.conjugate().T*SR*U-SX,
        'azimuthal_spin_frame':U.conjugate().T*SP*U-SY,
        'radial_measure_cancellation':sp.I*SX/(2*r)-sp.I*U.conjugate().T*SP*U.diff(ph)/r-SY*ell/r,
        'meridional_square':square,
    }
    theta,k,n=sp.symbols('theta k n',real=True)
    S=sp.cos(theta/2)*sp.eye(2)-sp.I*sp.sin(theta/2)*SY
    tangent=sp.sin(theta)*SX+sp.cos(theta)*SZ
    normal=-sp.cos(theta)*SX+sp.sin(theta)*SZ
    identities['tangent_spin_frame']=S.conjugate().T*tangent*S-SZ
    identities['normal_spin_frame']=S.conjugate().T*normal*S+SX
    identities['fermi_normal_connection_cancellation']=-sp.I*SZ*(-k)*S.conjugate().T*S.diff(theta)/(1-k*n)+sp.I*SX*k/(2*(1-k*n))
    beta=sp.symbols('beta',positive=True)
    gaussian=sp.Matrix([sp.exp(-beta*n*n/2),0])
    identities['local_oscillator_zero']=sp.I*SX*gaussian.diff(n)+SY*beta*n*gaussian
    mu,eta=sp.symbols('mu eta',real=True)
    radius2=(mu*mu-1)*(1-eta*eta); height=mu*eta
    numerator=((mu+eta)**4+(mu-eta)**4-2*(radius2+height**2-1)*(mu*mu-eta*eta))/4
    N=eta**4+3*eta**2*mu**2-eta**2+mu**2
    identities['dipole_field_strength']=sp.Matrix([sp.expand(numerator-N)])
    identities['field_lower_factor']=sp.Matrix([sp.expand(N-mu**2-eta**2*(eta**2+3*mu**2-1))])
    y=sp.symbols('y',positive=True)
    F=y/(sp.exp(y)-1)-sp.log(1-sp.exp(-y))
    identities['window_error_integral']=sp.Matrix([sp.diff(F,y)+y*sp.exp(-y)/(1-sp.exp(-y))**2])
    d,rp,u=sp.symbols('d rp u',positive=True)
    identities['finite_port_chord_margin']=sp.Matrix([(d/u-2*rp)-(d-2*rp)/u-2*rp*(1/u-1)])
    exact={name:all(sp.simplify(v)==0 for v in expr) for name,expr in identities.items()}
    assert all(exact.values()),exact

    # Independent Cartesian finite differences of the original Weyl equation.
    cartesian=[]
    def chi(r,z):
        e=np.exp(-r*r/3-z*z/5)
        return e*np.array([1+.2*r+1j*.3*z,.4*r-.1*z+1j*.2])
    def chi_r(r,z):return -2*r/3*chi(r,z)+np.exp(-r*r/3-z*z/5)*np.array([.2,.4])
    def chi_z(r,z):return -2*z/5*chi(r,z)+np.exp(-r*r/3-z*z/5)*np.array([.3j,-.1])
    for Q in [2,9,25]:
        for L in [.5,Q-.5]:
            for R,Z,P in [(.4,-.3,.2),(1.1,.4,.7),(2.,-.5,-.8)]:
                def psi(v):
                    x,y,z=v; r=np.hypot(x,y); phi=np.arctan2(y,x)
                    return np.exp(1j*np.array([L-.5,L+.5])*phi)*chi(r,z)/np.sqrt(r)
                point=np.array([R*np.cos(P),R*np.sin(P),Z]); h=1e-6
                deriv=[(psi(point+h*np.eye(3)[j])-psi(point-h*np.eye(3)[j]))/(2*h) for j in range(3)]
                A=potential(R,Z,Q); Avec=A/R*np.array([-np.sin(P),np.cos(P),0])
                physical=sum(-1j*mat@dv-mat@psi(point)*av for mat,dv,av in zip([sx,sy,sz],deriv,Avec))
                reduced=-1j*sx@chi_r(R,Z)-1j*sz@chi_z(R,Z)+sy@chi(R,Z)*(L-A)/R
                expected=np.exp(1j*np.array([L-.5,L+.5])*P)*reduced/np.sqrt(R)
                error=np.linalg.norm(physical-expected)/max(1,np.linalg.norm(expected))
                assert error<2e-8,(Q,L,R,error)
                cartesian.append(dict(q=Q,ell=L,r=R,z=Z,relative_residual=float(error)))

    geometry=[]
    for u in [.025,.1,.3,.6,.9,.975]:
        for eta in [-.85,-.4,0.,.3,.8]:
            mu=((1-eta*eta)+np.sqrt((1-eta*eta)**2+4*u*u*eta*eta))/(2*u)
            r=np.sqrt((mu*mu-1)*(1-eta*eta)); z=mu*eta
            q=13.; ell=q*(1-u); A=potential(r,z,q); B=field(r,z,q)
            mass=lambda r,z:(ell-potential(r,z,q))/r
            refinements=[]
            for step in [1e-5,3e-6,1e-6]:
                grad=np.array([(mass(r+step,z)-mass(r-step,z))/(2*step),(mass(r,z+step)-mass(r,z-step))/(2*step)])
                error=np.linalg.norm(grad-np.array([-B[1],B[0]]))/max(1,np.linalg.norm(B))
                refinements.append(dict(step=step,residual=float(error)))
            assert abs(A-ell)<1e-12 and error<2e-8,(u,eta,error)
            assert np.linalg.norm(B)>=q*u**3-1e-12
            geometry.append(dict(u=u,eta=eta,mass_gradient_residual=float(error),
                derivative_refinements=refinements,field_ratio_to_midplane=float(np.linalg.norm(B)/(q*u**3))))
    assert max(c['derivative_refinements'][-1]['residual'] for c in geometry)<max(c['derivative_refinements'][0]['residual'] for c in geometry)

    # Resolvent comparison on passive, noncommuting perturbed matrices.
    rng=np.random.default_rng(20261010)
    matrix=[]
    for Q in [3,6]:
        labels=(np.arange(Q)+.5)/Q
        raw=rng.normal(size=(Q,Q))+1j*rng.normal(size=(Q,Q)); K=(raw+raw.conj().T)/2
        G=expm(.7j*K)
        raw=rng.normal(size=(Q,Q))+1j*rng.normal(size=(Q,Q)); V=expm(1j*(raw+raw.conj().T)/2)
        for xi in [.05,.2,1.,4.]:
            half=np.diag(np.exp(-xi*(1/labels-1)/2))
            S0=half@G@half; epsilon=.001
            R0=np.exp(-xi)*S0
            R=np.exp(-xi)*((1-epsilon/2)*S0+epsilon/2*V)
            assert np.linalg.norm(R-R0,2)<=epsilon*np.exp(-xi)+1e-13
            for T in [.1,.5,1.]:
                X=np.exp(-xi*T)*R; X0=np.exp(-xi*T)*R0
                tr=np.trace(np.linalg.solve(np.eye(Q)+X,X))
                tr0=np.trace(np.linalg.solve(np.eye(Q)+X0,X0))
                bound=Q*epsilon*np.exp(-xi*(T+1))/(-np.expm1(-xi*(T+1)))**2
                assert abs(tr-tr0)<=bound+1e-12
                matrix.append(dict(q=Q,xi_d=xi,T_over_d=T,trace_difference=float(abs(tr-tr0)),bound=float(bound)))

    window=[]
    for y in [.02,.05,.1,.4,1.]:
        analytic=y/np.expm1(y)-np.log(-np.expm1(-y))
        integral=quad(lambda t:t*np.exp(-t)/(-np.expm1(-t))**2,y,np.inf,epsabs=1e-11)[0]
        assert abs(integral-analytic)<1e-10
        window.append(dict(y=y,integral=integral,residual=float(abs(integral-analytic))))
    # Explicit conservative exact certificate at delta=.05, y=.05,
    # x=4 and epsilon=.001. F(.05)<=1+log(21)<81/20.
    log2lo=2*sum(sp.Rational(1,3)**(2*n+1)/sp.Integer(2*n+1) for n in range(8))
    exp4lo=sum(sp.Integer(4)**n/sp.factorial(n) for n in range(24))
    exp305lo=sum(sp.Rational(61,20)**n/sp.factorial(n) for n in range(20))
    assert exp305lo>21
    certificate=sp.factor(sp.Rational(314159,100000)**2*(8*log2lo-5-sp.Rational(1,20))
        -240/exp4lo-sp.Rational(24,800)-sp.Rational(24*81,10000))
    assert certificate>0
    y=.05; epsilon=.001; delta=.05
    Fnum=y/np.expm1(y)-np.log(-np.expm1(-y))
    margin=np.pi*(8*np.log(2)-5-delta)/24-10/(np.pi*np.exp(4))-y*y/(2*np.pi)-2*epsilon/np.pi*Fnum
    return dict(status='passed_exact_reduction_and_supplementary_checks',exact_identities=exact,
        cartesian_reduction_cases=cartesian,field_line_geometry_cases=geometry,
        matrix_perturbation_cases=matrix,window_integral_cases=window,
        positive_rational_margin_certificate=str(certificate),log21_bound_certificate=str(exp305lo-21),
        accepted_window_example=dict(xi_d_min=y,xi_d_max=4.,return_operator_error=epsilon,
                                     delay_cdf_error=delta,margin_times_d_squared_over_q=float(margin)),
        exact_flat_exterior_differential_reduction=True,
        physical_port_return_derived=False,physical_return_error_bounded=False,
        complete_four_dimensional_energy_derived=False,novelty_established=False,
        full_numerical_error_certified=False,lean_certified=False,lean_invocations=0)


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output',required=True); args=p.parse_args()
    report=verify(); report['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if not k.endswith('_cases')},indent=2))
