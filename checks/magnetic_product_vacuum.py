"""Massless Dirac determinant and homogeneous magnetic-product pressures.

Zeta subtraction at fixed mu, Euclidean L2 heat kernel, unit charge.
Analytic UV/mode-tail bounds supplement (not certify) adaptive quadrature.
No wormhole construction or local Lean execution.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import zeta
import sympy as sp

PI=math.pi
EM_ORDER=10
H_ORDER=9
MODE_EXPONENT=80.0


def symbolic_checks():
    v,q,eta=sp.symbols('v q eta',positive=True)
    sphere=2/v
    for k in range(1,5):
        sphere+=2*sp.bernoulli(2*k)*sum((-v)**(2*k-j)*q**(2*k-2*j)
            /sp.factorial(2*k-2*j)/sp.factorial(j) for j in range(k+1))/v
    sphere/=4*sp.pi
    hyper=1/(2*sp.pi*v)+sum(sp.bernoulli(2*k+2)*v**k
        /((2*k+2)*sp.factorial(k)*sp.pi*eta**(k+1)) for k in range(3))
    A0=1/(4*sp.pi**2)
    A1=(1/eta-1)/(24*sp.pi**2)
    A2=(15*q*q-3-5/eta-3/eta**2)/(720*sp.pi**2)
    checks={}
    def zero(name,expression):
        assert sp.simplify(expression)==0,name
        checks[name]='exact_zero'
    series=sp.series(sphere*hyper,v,0,1).removeO().expand()
    zero('spectral_a0',series.coeff(v,-2)-A0)
    zero('spectral_a2',series.coeff(v,-1)-A1)
    zero('spectral_a4_vs_four_dimensional_anomaly',series.coeff(v,0)-A2)
    zero('pressure_subtraction_constant',A2+eta*sp.diff(A2,eta)
         -(15*q*q-3+3/eta**2)/(720*sp.pi**2))
    a,b,mu=sp.symbols('a b mu',positive=True)
    F=sp.Function('F')
    w=b**-4*F(a*a/b**2,sp.log(mu*mu*b*b))
    pA=-w-a*sp.diff(w,a)/2
    pS=-w-b*sp.diff(w,b)/2
    zero('independent_radius_variations_trace',2*pA+2*pS
         +2*b**-4*sp.Subs(sp.Derivative(F(eta,v),v),(eta,v),(a*a/b**2,sp.log(mu*mu*b*b))))
    u,z,eps=sp.symbols('u z epsilon',positive=True)
    eqs=sp.Matrix([-1/z+1/z**2,1/u-1/z**2])
    jac=eqs.jacobian([u,z]).subs({u:1,z:1})
    assert jac.det()==-1
    checks['classical_two_radius_Jacobian_determinant']='-1'
    PA,PS,lam=sp.symbols('P_A P_S loop')
    ub=1-lam*(2*PA+PS);zb=1-lam*PA
    e=eqs.subs({u:ub,z:zb})-lam*sp.Matrix([PA,PS])
    zero('sphere_radius_first_order',sp.diff(e[0],lam).subs(lam,0))
    zero('AdS_radius_first_order',sp.diff(e[1],lam).subs(lam,0))
    return checks


def coefficients(q,eta):
    """Laurent arrays: index 0 is v^-1, for each 2D factor."""
    sphere=np.zeros(2*EM_ORDER+1,dtype=np.longdouble)
    sphere[0]=1/(2*PI)
    for k in range(1,EM_ORDER+1):
        B=float(sp.bernoulli(2*k))/(2*PI)
        for j in range(k+1):
            power=2*k-j-1
            sphere[power+1]+=np.longdouble(B)*(-1)**(2*k-j)*np.longdouble(q)**(2*k-2*j)/math.factorial(2*k-2*j)/math.factorial(j)
    hyper=np.zeros(H_ORDER+2,dtype=np.longdouble);hyper[0]=1/(2*PI)
    phyper=hyper.copy()
    for k in range(H_ORDER+1):
        value=float(sp.bernoulli(2*k+2))/((2*k+2)*math.factorial(k)*PI*eta**(k+1))
        hyper[k+1]=value;phyper[k+1]=-k*value
    product=np.convolve(sphere,hyper);pressure=np.convolve(sphere,phyper)
    A=np.array([1/(4*PI**2),(1/eta-1)/(24*PI**2),
                (15*q*q-3-5/eta-3/eta**2)/(720*PI**2)])
    P=np.array([A[0],-1/(24*PI**2),(15*q*q-3+3/eta**2)/(720*PI**2)])
    assert np.max(np.abs(product[:3]-A))<2e-12*max(1,q*q)
    assert np.max(np.abs(pressure[:3]-P))<2e-12*max(1,q*q)
    return product,pressure,A,P


def uv_error_bounds(q,eta,v0):
    N=EM_ORDER;M=H_ORDER
    hermite_integral_bound=math.sqrt(PI)/2*math.sqrt(2**(2*N+1)*math.factorial(2*N+1))
    sphere_bound=zeta(2*N,1)*hermite_integral_bound*math.exp(v0*q*q/4)/(PI*(2*PI)**(2*N))
    h_bound=2*math.factorial(2*M+3)*zeta(2*M+4,1)/(math.factorial(M+1)*(2*PI)**(2*M+4)*PI*eta**(M+2))
    C_s=(2+3*q*v0+2*math.sqrt(2*v0/math.e))/(4*PI)
    C_sp=C_s+sphere_bound*v0**N
    C_h=1/(2*PI)+v0/(12*PI*eta)
    C_ph=1/(2*PI)+v0*v0/(120*PI*eta*eta)
    def bound(C,err):return C*sphere_bound*v0**(N-2)/(N-2)+C_sp*err*v0**M/M
    return float(bound(C_h,h_bound)),float(bound(C_ph,(M+1)*h_bound))


def kernels(v,q,eta):
    """Positive spectral factors. Inner quadrature tails beyond z=12 negligible."""
    N=max(1,math.ceil((-q+math.sqrt(q*q+4*MODE_EXPONENT/v))/2))
    assert v*(2*N+q)**2>=2
    n=np.arange(1,N+1,dtype=float)
    sphere=(q+2*np.sum((q+2*n)*np.exp(-v*n*(n+q))))/(4*PI)
    if v<1:
        # Split the exact flat term so its leading 1/v behavior is explicit.
        def thermal(lam,power):
            if lam==0:return 1/(2*PI) if power==1 else 0.0
            return lam**power*math.exp(-v*lam*lam/eta)/math.expm1(2*PI*lam)
        h=1/(2*PI*v)+2*quad(lambda l:thermal(l,1),0,12,epsabs=2e-13,epsrel=2e-13)[0]/(PI*eta)
        hp=1/(2*PI*v)+2*v*quad(lambda l:thermal(l,3),0,12,epsabs=2e-13,epsrel=2e-13)[0]/(PI*eta*eta)
    else:
        c=PI*math.sqrt(eta/v)
        def g(z,power):
            if z==0:return 1/c if power==1 else 0.0
            return z**power/math.tanh(c*z)*math.exp(-z*z)
        h=quad(lambda z:g(z,1),0,12,epsabs=2e-13,epsrel=2e-13)[0]/(PI*v)
        hp=quad(lambda z:g(z,3),0,12,epsabs=2e-13,epsrel=2e-13)[0]/(PI*v)
    return np.array([h*sphere,hp*sphere])


def vacuum(q,eta=1.0,log_mu_b_squared=0.0,cutoff_factor=1.0):
    assert q>=1 and eta>0
    v0=cutoff_factor*0.03/(q+1)
    products,pressures,A,P=coefficients(q,eta)
    uv=np.array([sum(float(c)*v0**(j-2)/(j-2) for j,c in enumerate(poly) if j>=3)
                 for poly in [products,pressures]])
    values=uv.copy();errors=np.zeros(2)
    for channel,sub in enumerate([A,P]):
        def integrand(v):
            return (kernels(v,q,eta)[channel]-sub[0]/(v*v)-sub[1]/v-sub[2])/v
        for lo,hi in [(v0,0.1),(0.1,1.0)]:
            out,err=quad(integrand,lo,hi,epsabs=3e-10,epsrel=3e-11,limit=150)
            values[channel]+=out;errors[channel]+=err
        def ir(h):
            if h==0:return q/(4*PI**2.5*math.sqrt(eta))*(1 if channel==0 else 0.5)
            return 2*kernels(1/(h*h),q,eta)[channel]/h
        out,err=quad(ir,0,1,epsabs=3e-11,epsrel=3e-11,limit=120)
        values[channel]+=out;errors[channel]+=err
    finite=np.array([values[0]-A[0]/2-A[1]+A[2]*(np.euler_gamma+log_mu_b_squared),
                     values[1]-P[0]/2-P[1]+P[2]*(np.euler_gamma+log_mu_b_squared)])/2
    # F=b^4 w_E; pressure channel gives F+eta F_eta.
    pA=-finite[1];pS=finite[1]-A[2]/2
    uvbounds=uv_error_bounds(q,eta,v0)
    tail=math.exp(-MODE_EXPONENT)/(2*PI)*(1/(4*PI*v0*v0)+1/(3*PI**1.5*math.sqrt(eta)*v0**1.5))
    return {'q':q,'eta':eta,'log_mu_b_squared':log_mu_b_squared,'uv_cutoff':v0,
            'b4_Euclidean_action_density':float(finite[0]),'b4_pressure_AdS':float(pA),
            'b4_pressure_sphere':float(pS),'b4_trace':float(-A[2]),
            'uv_remainder_bounds_before_half':list(uvbounds),'angular_tail_integral_bound_before_half':tail,
            'outer_quadrature_error_estimates_before_half':errors.tolist(),
            'quadrature_error_estimates_are_not_certified_bounds':True}


def numerical_checks():
    samples=[]
    for q in [3,8,16,64]:
        item=vacuum(q);comparison=vacuum(q,cutoff_factor=0.65)
        residual=max(abs(item[k]-comparison[k]) for k in ['b4_Euclidean_action_density','b4_pressure_AdS','b4_pressure_sphere'])
        assert residual<3e-7
        item['cutoff_change_max_residual']=residual
        # Classical b^2=kappa*q^2/(8g_c^2). Coefficients of g_c^2.
        item['sphere_fractional_radius_squared_shift_per_g_c_squared']=-8*item['b4_pressure_AdS']/(q*q)
        item['AdS_fractional_radius_squared_shift_per_g_c_squared']=-8*(2*item['b4_pressure_AdS']+item['b4_pressure_sphere'])/(q*q)
        samples.append(item)
    def w(a,b):return vacuum(3,a*a/(b*b),math.log(b*b))['b4_Euclidean_action_density']/b**4
    base=vacuum(3);derivatives=[]
    for h in [0.004,0.002]:
        da=(w(1+h,1)-w(1-h,1))/(2*h)
        db=(w(1,1+h)-w(1,1-h))/(2*h)
        derivatives.append([-base['b4_Euclidean_action_density']-da/2,
                            -base['b4_Euclidean_action_density']-db/2])
    extrapolated=(4*np.array(derivatives[1])-np.array(derivatives[0]))/3
    expected=np.array([base['b4_pressure_AdS'],base['b4_pressure_sphere']])
    residual=float(np.max(np.abs(extrapolated-expected)))
    assert residual<2e-8
    return {'samples':samples,'independent_two_radius_finite_difference_pressures':extrapolated.tolist(),
            'finite_difference_vs_direct_pressure_max_residual':residual}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data={'status':'passed_magnetic_product_vacuum_feasibility',
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'sympy_version':sp.__version__,'scipy_version':scipy.__version__,
          'symbolic_checks':symbolic_checks(),'numerical_checks':numerical_checks(),
          'state':'Euclidean complete H2 L2 spectral prescription, continued to invariant AdS2 vacuum',
          'renormalization':'zeta at fixed mu; finite additional local terms zero in this stated scheme',
          'scope':'Homogeneous prescribed product vacuum and first-order stationary reference',
          'full_numerical_error_certified':False,'global_wormhole_constructed':False,
          'physical_stability_proved':False,'lean_certified':False,'novelty_established':False}
    args.output.write_bytes((json.dumps(data,indent=2)+'\n').encode())


if __name__=='__main__':main()
