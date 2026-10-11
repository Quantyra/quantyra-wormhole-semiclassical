"""Finite spherical ports of the flat point dipole; no scattering certification."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize_scalar
import sympy as sp


def verify():
    x, e, q = sp.symbols('x e q', positive=True)
    c = 1-2*x
    S = sp.sqrt(1+e*e-2*e*c)
    D = S*(S+1-e*c)
    w = sp.log((S+1-e)/((1-e)*(S+1+e)))/2
    delta = q*(1-(1-e*c)/S)/2
    A = q*x+delta
    th, ell = sp.symbols('theta ell', real=True)
    af = sp.Function('A')(th)
    s = (ell-af)/sp.sin(th)
    t = s+sp.cot(th)
    f = sp.Function('f')(th)
    SX = sp.Matrix([[0,1],[1,0]])
    SY = sp.Matrix([[0,-sp.I],[sp.I,0]])
    SZ = sp.diag(1,-1)
    R = sp.cos(th/2)*sp.eye(2)-sp.I*sp.sin(th/2)*SY
    n, k = sp.symbols('n k', positive=True)
    B = sp.Matrix([[0,sp.sqrt(2*n+2)],[sp.sqrt(2*n+2),2*k]])
    Binv = sp.Matrix([[-2*k,sp.sqrt(2*n+2)],[sp.sqrt(2*n+2),0]])/(2*n+2)
    m = sp.symbols('m', integer=True, nonnegative=True)
    p = x**m*(1-x)**(q-m-1)
    identities = {
        'gauge_rationalization':delta-2*q*e*e*x*(1-x)/D,
        'mode_log_derivative':sp.diff(w,x)-e*e/D,
        'zero_mode_equation':sp.diff(w,x)*2*x*(1-x)-delta/q,
        'normal_magnetic_field':sp.diff(A,x)-q*(1+e*e*(c-e)/S**3),
        'gap_factorization':s*s-sp.diff(s,th)-(t*t+sp.diff(t,th)+1+2*sp.diff(af,th)/sp.sin(th)),
        'gap_positive_factor':(sp.diff(-sp.diff(f,th)+t*f,th)+t*(-sp.diff(f,th)+t*f))-(-sp.diff(f,th,2)+(t*t+sp.diff(t,th))*f),
        'D_x':sp.diff(D,x)-e*(3*S*S+4*S+1-e*e)/S,
        'beta_stein':sp.diff(-x*(1-x)/(q+1),x)-x*(1-x)/(q+1)*(m/x-(q-m-1)/(1-x))-(x-(m+1)/(q+1)),
        # Positive square-root endpoints use 0<e<1 explicitly.
        'north_W_endpoint':sp.log(((1-e)+1-e)/((1-e)*((1-e)+1+e)))/2,
        'south_exponential_W_endpoint':((1+e)+1-e)/((1-e)*((1+e)+1+e))*(1-e*e)-1,
        'radial_spin_rotation':R.conjugate().T*(sp.sin(th)*SX+sp.cos(th)*SZ)*R-SZ,
        'polar_spin_rotation':R.conjugate().T*(sp.cos(th)*SX-sp.sin(th)*SZ)*R-SX,
        'radial_measure_connection':-sp.I*SX*R.conjugate().T*sp.diff(R,th)+sp.I*SZ/2,
        'curved_edge_block_inverse':B*Binv-sp.eye(2),
    }
    exact = {}
    for name, expr in identities.items():
        if 'endpoint' in name:
            expr = sp.expand_log(sp.factor(expr),force=True)
            expr = sp.refine(expr,sp.Q.positive(1-e))
            expr = sp.expand_log(expr,force=True)
        vals = list(expr) if isinstance(expr,sp.MatrixBase) else [expr]
        exact[name] = all(sp.simplify(z)==0 for z in vals)
    assert all(exact.values()), exact
    rational = {
        'second_derivative_coefficient':sp.Rational(225750,59049)<4,
        'mixed_derivative_coefficient':2/sp.Rational(81,50)+sp.Rational(1,10)*sp.Rational(77,15)/sp.Rational(81,50)**2<sp.Rational(3,2),
        'convexity_q_at_least_three':1-sp.Rational(1,3)-8*sp.Rational(1,10)**3>sp.Rational(1,2),
        'normal_field_lower':1-(sp.Rational(1,10)/sp.Rational(9,10))**2==sp.Rational(80,81),
    }
    assert all(rational.values()),rational
    wf = sp.lambdify((x,e),w,'numpy')
    we = sp.lambdify((x,e),sp.diff(w,e),'numpy')

    def distribution(Q,M,E,tau):
        def logp(z):
            return M*np.log(z)+(Q-M-1)*np.log1p(-z)-tau*Q*wf(z,E)
        peak = minimize_scalar(lambda z:-logp(z),bounds=(1e-12,1-1e-12),method='bounded',options={'xatol':1e-13}).x
        scale = max(logp(peak),logp(1e-15),logp(1-1e-15))
        points=sorted(set([peak]+[v for v in [peak-5/np.sqrt(Q),peak+5/np.sqrt(Q)] if 0<v<1]))
        integ=lambda fun:quad(lambda z:np.exp(logp(z)-scale)*fun(z),0,1,points=points,epsabs=2e-12,epsrel=2e-11,limit=250)[0]
        norm=integ(lambda z:1.)
        return scale+np.log(norm),lambda fun:integ(fun)/norm

    profiles=[]
    for Q in [3,8,31,101]:
        for M in [0,Q//2,Q-1]:
            for E in [.02,.05,.1]:
                z0,_=distribution(Q,M,E,0)
                z1,_=distribution(Q,M,E,1)
                z2,expect=distribution(Q,M,E,2)
                overlap=np.exp(z1-(z0+z2)/2)
                distance=np.sqrt(max(0.,2-2*overlap))
                bound=50*np.sqrt(2*Q)/81*E*E
                mean=expect(lambda z:Q*we(z,E))
                derivative=np.sqrt(expect(lambda z:(Q*we(z,E)-mean)**2))
                radial_bound=3/np.sqrt(2)*np.sqrt(Q)*E # d=1
                assert distance<=bound+1e-10,(Q,M,E,distance,bound)
                assert derivative<=radial_bound+1e-10
                # Independent finite difference of the normalized probability amplitude.
                de=1e-5
                zp,_=distribution(Q,M,E+de,2);zm,_=distribution(Q,M,E-de,2)
                fd=expect(lambda z:((np.exp(-Q*(wf(z,E+de)-wf(z,E))-(zp-z2)/2)-np.exp(-Q*(wf(z,E-de)-wf(z,E))-(zm-z2)/2))/(2*de))**2)
                assert abs(np.sqrt(fd)-derivative)<2e-7
                profiles.append(dict(q=Q,m=M,epsilon=E,overlap=float(overlap),mode_distance=float(distance),distance_bound=float(bound),radial_derivative_norm_d1=float(derivative),radial_derivative_bound_d1=float(radial_bound),finite_difference_derivative_norm=float(np.sqrt(fd))))

    distortion=[]
    for Q in [100,1000,10000,100000,1000000]:
        M=Q//2;E=.1
        zs=[distribution(Q,M,E,tau)[0] for tau in [0,1,2]]
        overlap=np.exp(zs[1]-(zs[0]+zs[2])/2)
        distortion.append(dict(q=Q,m=M,epsilon=E,overlap=float(overlap)))
    assert all(a['overlap']>b['overlap'] for a,b in zip(distortion,distortion[1:]))
    assert distortion[-1]['overlap']<.1

    afun=sp.lambdify((x,e,q),A,'numpy')
    axfun=sp.lambdify((x,e,q),sp.diff(A,x),'numpy')
    spectra=[]
    for Q in [3,7,15]:
        for M in [0,Q//2,Q-1]:
            for E in [0.,.1]:
                values=[]
                for N in [300,600,1200]:
                    h=np.pi/(N+1);theta=h*np.arange(1,N+1);X=np.sin(theta/2)**2
                    ss=(M+.5-afun(X,E,Q))/np.sin(theta)
                    sprime=-axfun(X,E,Q)/2-ss/np.tan(theta)
                    val=eigh_tridiagonal(2/h**2+ss*ss-sprime,np.full(N-1,-1/h**2),select='i',select_range=(0,0))[0][0]
                    values.append(float(val))
                lower=1+Q*(1-(E/(1-E))**2)
                assert values[-1]>=lower-.001*(Q+1)
                assert abs(values[2]-values[1])<abs(values[1]-values[0])
                if E==0:assert abs(values[-1]-(Q+1))<.001*(Q+1)
                spectra.append(dict(q=Q,m=M,epsilon=E,grid_sizes=[300,600,1200],lowest_lower_spin_squared_eigenvalues=values,proved_squared_gap_lower=float(lower)))

    # Test the whole 3D differential operator on independent smooth spinors.
    sx=np.array(SX).astype(complex);sy=np.array(SY).astype(complex);sz=np.array(SZ).astype(complex)
    def rotate(theta):return np.cos(theta/2)*np.eye(2)-1j*np.sin(theta/2)*sy
    def eta(rho,theta):return np.exp(-rho)*np.array([1+.2*np.cos(theta)+.3j*rho,.4*np.sin(theta)+.1j*rho])
    def eta_r(rho,theta):return -eta(rho,theta)+np.exp(-rho)*np.array([.3j,.1j])
    def eta_t(rho,theta):return np.exp(-rho)*np.array([-.2*np.sin(theta),.4*np.cos(theta)])
    cartesian=[]
    for Q in [3,11]:
        for E in [.05,.1]:
            for theta in [.65,1.2,2.2]:
                L=Q//2+.5;phi=.7
                def psi(v):
                    X,Y,Z=v;rho=np.linalg.norm(v);tt=np.arccos(Z/rho);pp=np.arctan2(Y,X)
                    return np.exp(1j*np.array([L-.5,L+.5])*pp)*(rotate(tt)@eta(rho,tt))/(rho*np.sqrt(np.sin(tt)))
                point=E*np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
                X,Y,Z=point;rr=np.hypot(X,Y)
                # Independent Cartesian point-dipole potential, centers z=0,1.
                aval=Q*(1-(Z/E-(Z-1)/np.linalg.norm(point-np.array([0,0,1])))/2)
                avec=aval/rr*np.array([-np.sin(phi),np.cos(phi),0.])
                errors=[]
                for h in [3e-6,1e-6,3e-7]:
                    deriv=[(psi(point+h*np.eye(3)[j])-psi(point-h*np.eye(3)[j]))/(2*h) for j in range(3)]
                    actual=sum(-1j*mat@dv-mat@psi(point)*av for mat,dv,av in zip([sx,sy,sz],deriv,avec))
                    ss=(L-afun(np.sin(theta/2)**2,E,Q))/np.sin(theta)
                    reduced=-1j*sz@eta_r(E,theta)+(-1j*sx@eta_t(E,theta)+sy@eta(E,theta)*ss)/E
                    expected=np.exp(1j*np.array([L-.5,L+.5])*phi)*(rotate(theta)@reduced)/(E*np.sqrt(np.sin(theta)))
                    errors.append(float(np.linalg.norm(actual-expected)/max(1,np.linalg.norm(expected))))
                assert errors[-1]<1e-7,(Q,E,theta,errors)
                cartesian.append(dict(q=Q,epsilon=E,theta=theta,steps=[3e-6,1e-6,3e-7],relative_residuals=errors))

    return dict(status='passed_finite_spherical_mouth_modes_and_source_inverse_checks',exact_identities=exact,rational_inequalities={k:bool(v) for k,v in rational.items()},profile_cases=profiles,large_flux_distortion_cases=distortion,angular_spectrum_cases=spectra,cartesian_operator_cases=cartesian,source_specific_corrector_bound='27/sqrt(160) * (rho/d)^2 for q>=3 and 0<rho/d<=0.1',numerical_checks_are_supplementary=True,physical_port_return_derived=False,physical_return_error_bounded=False,novelty_established=False,lean_certified=False,lean_invocations=0,versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);args=parser.parse_args()
    result=verify();result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:len(result[k]) for k in ['exact_identities','profile_cases','large_flux_distortion_cases','angular_spectrum_cases','cartesian_operator_cases']}))
