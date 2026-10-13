"""S072: fixed-source axial Einstein-Maxwell response and projected LLL bound.

Ordinary variational proof: notes/angular-boundary-response.md.
Numerics are supplementary Ritz upper bounds, not a four-dimensional quench.
"""
from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import scipy.linalg as la
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.integrate import quad
from scipy.special import beta
import sympy as sy

warnings.filterwarnings("error")
ROOT = Path(__file__).resolve().parents[1]


def symbolic_checks():
    u, x = sy.symbols("u x", positive=True)
    h = x + 2*x*u + (1+3*x)*u*u
    w = (2*(1+2*x)*u**3-(1+3*x)*u**4)/(1+x)
    v = (u*u+x)/(1+x)
    mix = 2*(1+3*x)/(3*(1-u)**2*h)
    ew = sy.diff(sy.diff(w,u)/(3*u*u),u)+2*mix*(v-w)
    ev = sy.diff(2*(1+3*x)*sy.diff(v,u)/3,u)-2*mix*(v-w)
    assert sy.simplify(ew) == sy.simplify(ev) == 0
    energy = sy.diff(w,u)**2/(6*u*u)+(1+3*x)*sy.diff(v,u)**2/3+mix*(v-w)**2
    kn = sy.factor(sy.integrate(sy.cancel(energy),(u,0,1)))
    assert sy.simplify(kn-(1+2*x)/(1+x)) == 0
    assert sy.simplify((v-w)-(1-u)**2*h/(1+x)) == 0
    # The trial independently minimizes the two radial derivative terms.
    derivative = sy.integrate(9*u**2/6+(1+3*x)/3,(u,0,1))
    assert sy.simplify(derivative-(5+6*x)/6) == 0
    # LLL <cos(theta)> from its normalized beta density.
    q,m = sy.symbols("q m",positive=True)
    j=(q-1)/2
    mean=1-2*(j+m+1)/(2*j+2)
    assert sy.simplify(mean+2*m/(q+1)) == 0
    return {"KN_Euler_Lagrange_residuals":[0,0], "KN_functional":str(kn),
            "fixed_source_radial_lower_bound":str(derivative),
            "LLL_cosine_projection":"-2 m/(q+1)"}


def bounds(x,q,gamma,tau=0.0):
    y=1+3*x-tau
    assert y>0 and tau>=0
    h=lambda u:x+2*x*u+(1+3*x)*u*u
    d=np.arctan(np.sqrt(1+2*x)/(2*np.sqrt(x)))/np.sqrt(x*(1+2*x))
    b=quad(lambda u:u*u*(1+u)**2/h(u),0,1,epsabs=1e-12)[0]
    p=quad(lambda u:(1+(u+u*u)/(q+1))**2/h(u),0,1,epsabs=1e-12)[0]
    lo=.5+y/3
    hi=lo+2*y*b/3+gamma*p
    simple=lo+2*y/(9*x)+gamma*(1+2/(q+1))**2*d
    assert b <= 1/(3*x)*(1+1e-12)
    assert hi <= simple*(1+1e-12)
    return lo,hi,simple,d,b,p


def ritz(x,q,gamma,degree,nu=0.0,tau=0.0):
    """Minimize the full positive quadratic form over a nested polynomial space.

    u=b/r; W(0)=0, W(1)=V(1)=1, V(0)=nu. Gauss integrates the
    derivative polynomials exactly and resolves the smooth factored weights.
    """
    nodes,weights=leggauss(max(160,4*degree+24))
    u=(nodes+1)/2; weights=weights/2
    p=np.column_stack([Legendre.basis(i)(nodes) for i in range(degree)])
    dp=np.column_stack([2*Legendre.basis(i).deriv()(nodes) for i in range(degree)])
    z=np.zeros_like(p)
    bw=u[:,None]**3*(1-u[:,None])*p
    bv=u[:,None]*(1-u[:,None])*p
    dw=(3*u*u-4*u**3)[:,None]*p+(u**3*(1-u))[:,None]*dp
    dv=(1-2*u)[:,None]*p+(u*(1-u))[:,None]*dp
    ymag=1+3*x-tau
    h=x+2*x*u+(1+3*x)*u*u+tau*u**3/(1-u)
    # (V-W)/(1-u) and [1-(W+q V)/(q+1)]/(1-u), without cancellation.
    base=[3*u*u, np.full_like(u,1-nu), nu+(1+u)*u,
          (1+u+u*u+q*(1-nu))/(q+1)]
    mat=[np.column_stack([dw,z]),np.column_stack([z,dv]),
         np.column_stack([-u[:,None]**3*p,u[:,None]*p]),
         np.column_stack([-u[:,None]**3*p,-q*u[:,None]*p])/(q+1)]
    scale=[weights/(6*u*u),weights*ymag/3,
           weights*2*ymag/(3*h),weights*gamma/h]
    a=np.vstack([m*np.sqrt(s)[:,None] for m,s in zip(mat,scale)])
    y=np.concatenate([b*np.sqrt(s) for b,s in zip(base,scale)])
    # Least squares avoids squaring the condition number of the basis.
    coef,_,_,_=la.lstsq(a,-y,lapack_driver="gelsy")
    residual=a@coef+y
    val=float(residual@residual)
    station=float(np.max(np.abs(a.T@residual)))/max(1.0,val)
    return val,station


def projection_checks():
    rows=[]
    for q in (2,3,8,19):
        j=(q-1)/2
        for m in (-j,0.5 if q%2==0 else 0,j):
            aa=j+m+1; bb=j-m+1
            den=beta(aa,bb)
            mean=quad(lambda t:(1-2*t)*t**(aa-1)*(1-t)**(bb-1)/den,
                      0,1,epsabs=2e-12)[0]
            err=abs(mean+2*m/(q+1))
            assert err<2e-12
            rows.append({"q":q,"m":m,"cosine_error":err})
    return rows


def scale_family(k):
    q=k**4; nf=1; eta=1/k; g=k**-4; L=k**5; z0=.5
    # X(1+3X)=pi k^6; the rationalized root avoids small-X cancellation.
    x=2*np.pi*k**6/(1+np.sqrt(1+12*np.pi*k**6))
    c=1+6*x; b2=x*L**2; n=2*q
    gamma=nf*g*g*(q*q-1)*(1+3*x)/(12*np.pi**2*q)
    lo,hi,simple,d,_,_=bounds(x,q,gamma)
    xe=np.pi*q/(8*nf*eta*z0*c*(1+3*x))
    chi=np.pi*(q*q-1)/(8*nf*q*eta*z0*c*lo)
    rows={"k":k,"N_f":nf,"q":q,"X":x,"eta":eta,"g":g,
          "L_over_lp":L,"I_G_over_b3_interval":[lo,hi],
          "relative_interval_width":(hi-lo)/lo,
          "relative_width_safe_bound":(2*(1+3*x)/(9*x)+gamma*(1+2/(q+1))**2*d)/lo,
          "x_e_bare_upper":xe,"chi_J_static_upper":chi,
          "LLL_gap_squared":(q+1)/c,
          "electric_screening_upper":nf*q*g*g*d/(4*np.pi**2),
          "angular_LLL_trial_increment_over_lower":gamma*(1+2/(q+1))**2*d/lo,
          "species_NG_over_b2":n/b2,"Neta_over_Phi_ext":2*n*eta/b2,
          "Neta":n*eta,"Tb":n*eta*c/(2*np.pi*b2)}
    assert xe<1/k and chi<1/k and rows["LLL_gap_squared"]>k/7
    return rows


def main():
    symbolic=symbolic_checks()
    projection=projection_checks()
    numerical=[]; kn=[]
    for x in (.25,1,10,10000):
        value,station=ritz(x,10,0,4,x/(1+x))
        expected=(1+2*x)/(1+x)
        assert abs(value-expected)<2e-10
        kn.append({"X":x,"Ritz":value,"exact":expected,"stationarity":station})
        for q,gamma in ((3,0),(10,.1),(100,10)):
            for tau in (0,.1):
                lo,hi,simple,_,_,_=bounds(x,q,gamma,tau)
                vals=[ritz(x,q,gamma,dg,tau=tau) for dg in (12,24,48)]
                assert all(lo-1e-8<=v<=hi+1e-8 for v,_ in vals)
                assert vals[2][0]<=vals[1][0]+1e-8<=vals[0][0]+2e-8
                assert max(t for _,t in vals)<2e-9
                rel=abs(vals[-1][0]-vals[-2][0])/vals[-1][0]
                assert rel<2e-5
                numerical.append({"X":x,"q":q,"gamma":gamma,"tau_4pi_rT":tau,"lower":lo,
                    "trial_upper":hi,"simple_upper":simple,
                    "Ritz_12_24_48":[v for v,_ in vals],"relative_refinement":rel,
                    "stationarity":max(t for _,t in vals)})
    receipt={"date":"2026-10-13","scope":"static axial Einstein-Maxwell Dirichlet problem at zero and positive temperature, with projected Gaussian LLL chemical polarization",
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "symbolic":symbolic,"LLL_projection":projection,
        "KN_wrong_boundary_comparator":kn,"fixed_source_Ritz_cases":numerical,
        "single_species_family":[scale_family(k) for k in (2000,4000,8000)],
        "limitations":["Ritz numbers are upper approximations; the interval proof is analytic",
            "Thermal/background stress, full KK determinant and real-time response not established",
            "Using this static susceptibility in the compact rotor still requires dynamical matching",
            "No initial-pair preparation or four-dimensional stability result"],"passed":True}
    (ROOT/"evidence/angular-boundary-response.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8",newline="\n")
    print("PASS: magnetic KN boundary source, coupled static form, 12 LLL projections, 24 Ritz cases, single-species scales")


if __name__=="__main__":
    main()
