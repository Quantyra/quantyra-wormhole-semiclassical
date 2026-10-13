"""S072 causal Gaussian boundary-transfer response on a fixed RN exterior.

The analytic half-line argument is in notes/causal-compact-response.md.
Finite passive regulators below check the spectral and time-domain identities;
they are not claimed to reproduce an ingoing horizon at all times.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import warnings

import numpy as np
import scipy.linalg as la
from numpy.polynomial.legendre import leggauss
import sympy as sy

warnings.filterwarnings("error")
ROOT=Path(__file__).resolve().parents[1]


def algebra():
    r,G,g,q,n=sy.symbols("r G g q n",positive=True)
    P=q/2; mu=8*sy.pi/(3*g*g); lam=n*q*(q*q-1)/(12*sy.pi)
    kg=r**4/(6*G); ke=4*sy.pi*r*r/(3*g*g)
    c=sy.Matrix([-P,1]); d=sy.Matrix([1,2])/(q+1); nv=sy.Matrix([1,P])
    kinv=sy.diag(1/kg,1/ke)
    mass=mu*c*c.T+lam*d*d.T
    aa=sy.factor(mu*(c.T*kinv*c)[0])
    dd=sy.factor(lam*(d.T*kinv*d)[0])
    cross=(c.T*kinv*d)[0]
    det=sy.factor(aa*dd-mu*lam*cross**2)
    assert sy.simplify(mass.det()-mu*lam)==0
    assert sy.simplify(c.dot(nv))==0 and sy.simplify(d.dot(nv)-1)==0
    assert sy.simplify(aa-(2/r**2+4*sy.pi*G*q*q/(g*g*r**4)))==0
    assert sy.simplify(det-mu*lam/(kg*ke))==0
    assert sy.simplify(dd-mu*lam*cross**2/aa-lam/(kg+P*P*ke))==0
    # Kodama-Ishibashi (74), n_sphere=2, K=1: 2/r^2+4Q^2/r^4.
    return {"heavy_mass_squared":str(aa),"light_diagonal_mass_squared":str(dd),
            "mass_determinant":str(det),"boundary_source_projection":[0,1],
            "mixing_upper_squared":str(sy.factor(4*lam/(mu*(q+1)**2)))}


def finite_regulator(x,q,G,tau,kind,cells):
    """Linear finite elements in u=1/r, with reflecting passive endpoints.

    The horizon cut u=1-1e-5 is explicit. All matrix integrals use independent
    Gauss quadrature, and the two source/response calculations are compared.
    """
    y=1+3*x-tau; g2=np.pi*G*q*q/y
    lam=q/np.pi if kind=="electric" else q*(q*q-1)/(12*np.pi)
    mu=8*np.pi/(3*g2); P=q/2
    eps2=0.0 if kind=="electric" else 4*lam/(mu*(q+1)**2)
    assert 0<=eps2<1
    nf=1 if kind=="electric" else 2
    # Map resolves the near-horizon region while retaining u=0 exactly.
    grid=1-np.exp(np.linspace(0,np.log(1e-5),cells+1))
    size=nf*(cells+1)
    stiff=np.zeros((size,size)); mass=np.zeros_like(stiff)
    gx,gw=leggauss(10)
    for j in range(cells):
        left,right=grid[j:j+2]; h=right-left
        u=left+(gx+1)*h/2; weight=gw*h/2
        basis=np.stack([(1-gx)/2,(1+gx)/2],axis=1)
        derivative=np.array([-1,1])/h
        F=(1-u)**2*(x+2*x*u+(1+3*x)*u*u)+tau*u**3*(1-u)
        p=np.einsum("k,ki,kj->ij",weight/F,basis,basis)
        kgrad=np.outer(derivative,derivative)*np.dot(weight,F)
        if kind=="electric":
            pot=np.full((len(u),1,1),lam*g2/(4*np.pi))
        else:
            pot=np.empty((len(u),2,2))
            pot[:,0,0]=2+4*y*u*u
            pot[:,0,1]=pot[:,1,0]=np.sqrt(mu*lam)/(q+1)*(-6*G*P*u*u+3*g2/(2*np.pi))
            pot[:,1,1]=lam/(q+1)**2*(6*G*u*u+3*g2/np.pi)
        for a in range(2):
            for b in range(2):
                sl1=slice(nf*(j+a),nf*(j+a+1));sl2=slice(nf*(j+b),nf*(j+b+1))
                stiff[sl1,sl2]+=kgrad[a,b]*np.eye(nf)+np.einsum("k,k,k,kij->ij",weight,basis[:,a],basis[:,b],pot)
                mass[sl1,sl2]+=p[a,b]*np.eye(nf)
    source=np.zeros(size);source[nf-1]=1
    static_vec=la.solve(stiff,source,assume_a="pos")
    static=lam*np.dot(source,static_vec)
    ev,vec=la.eigh(stiff,mass,driver="gvd")
    assert min(ev)>0
    # The small static component can be sensitive to eigenvalue conditioning.
    spectral=lam*(source@vec)**2/ev
    identity_error=abs(sum(spectral)-static)/static
    assert identity_error<2e-7
    prob=spectral/static
    mixing=0.0 if nf==1 else float(static_vec[0]/static_vec[1])
    assert mixing*mixing<=eps2/(1-eps2)+1e-8
    response=[]
    a_time=static/(lam*(1-eps2))
    for ratio in (.01,.1,1,10):
        s=ratio/a_time
        direct=lam*np.dot(source,la.solve(stiff+s*s*mass,source,assume_a="pos"))/static
        by_modes=float(np.sum(prob*ev/(ev+s*s)))
        assert abs(direct-by_modes)<2e-7
        lower=1/(1+s*a_time)
        assert lower-2e-5<=direct<=1+1e-8
        response.append({"s_times_a":ratio,"response_ratio":direct,"lower":lower,
                         "spectral_error":abs(direct-by_modes)})
    # Explicit causal convolution, independently of Laplace solve.
    drive_duration=100*a_time
    scaled_freq=np.sqrt(ev)*drive_duration
    times=np.linspace(0,4,121)
    errors=np.array([abs(np.dot(prob,cubic_memory(scaled_freq,t))) for t in times])
    analytic=4*np.sqrt(6)*a_time/drive_duration
    assert max(errors)<=analytic+2e-7
    return {"X":x,"q":q,"G_over_R2":G,"tau_4pi_RT":tau,"kind":kind,
            "cells":cells,"horizon_cut_1_minus_u":1e-5,"static_coefficient":float(static),
            "epsilon_squared":eps2,"static_mixing_ratio":mixing,
            "spectral_identity_relative_error":identity_error,
            "minimum_squared_frequency":float(min(ev)),"Laplace_comparisons":response,
            "cubic_max_sampled_error":float(max(errors)),"cubic_analytic_bound":analytic}


def cubic_memory(w,t):
    b=min(t,1.0)
    result=np.empty_like(w,dtype=float)
    low=w<.2
    if np.any(low):
        nodes,weights=leggauss(24);u=(nodes+1)*b/2
        result[low]=np.cos(w[low,None]*(t-u))@(weights*b/2*(6*u-6*u*u))
    high=~low
    wh=w[high];iw=1j*wh; ex=np.exp(-iw*b)
    i0=(1-ex)/iw; i1=(-b*ex+i0)/iw;i2=(-b*b*ex+2*i1)/iw
    result[high]=np.real(np.exp(1j*wh*t)*(6*i1-6*i2))
    return result


def envelope_check():
    nodes,weights=leggauss(160)
    cases=[]
    for w in (.01,.1,.5,2,10,30,80):
        for t in (.1,.5,1,1.3,3):
            b=min(t,1);u=(nodes+1)*b/2
            exact=float(np.dot(weights*b/2*(6*u-6*u*u),np.cos(w*(t-u))))
            calculated=float(cubic_memory(np.array([w]),t)[0])
            assert abs(exact-calculated)<2e-12
            assert abs(calculated)<=min(1,24/w**2)+2e-12
            cases.append(abs(exact-calculated))
    return {"cases":len(cases),"max_independent_quadrature_error":max(cases)}


def candidates():
    spec=importlib.util.spec_from_file_location("angular_boundary_response",Path(__file__).with_name("angular_boundary_response.py"))
    dependency=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dependency)
    bounds=dependency.bounds
    rows=[]
    for k in (2000,4000,8000):
        q=k**4;g=k**-4;eta=1/k
        x=2*np.pi*k**6/(1+np.sqrt(1+12*np.pi*k**6))
        c=1+6*x;y=1+3*x
        gamma=g*g*(q*q-1)*y/(12*np.pi**2*q)
        _,dhi,_,dx,_,_=bounds(x,q,gamma)
        eps2=g*g*q*(q-1)/(8*np.pi**2*(q+1))
        epse=q*g*g*dx/(4*np.pi**2)
        for delta in (.25,.5):
            ratio_j=24*np.pi*eta*c*dhi/(delta*(q*q-1)*(1-eps2))
            ratio_e=8*np.pi*eta*c*y*(1+epse)/(delta*q*q)
            bound_j=4*np.sqrt(6)*ratio_j;bound_e=4*np.sqrt(6)*ratio_e
            assert bound_j<2.33e-6*(2000/k)**3*(.25/delta)
            assert bound_e<2.33e-6*(2000/k)**3*(.25/delta)
            rows.append({"k":k,"q":q,"delta":delta,"epsilon_angular_squared":eps2,
                "electric_screening_upper":epse,"a_over_drive_duration_angular_upper":ratio_j,
                "a_over_drive_duration_electric_upper":ratio_e,
                "cubic_uniform_relative_error_angular_upper":bound_j,
                "cubic_uniform_relative_error_electric_upper":bound_e})
    return rows


def main():
    alg=algebra();envelope=envelope_check();runs=[]
    for x,q,G,tau in ((.5,3,.01,.1),(3,7,.005,.2),(12,11,.002,.1)):
        for kind in ("electric","angular"):
            for cells in (100,180):
                runs.append(finite_regulator(x,q,G,tau,kind,cells))
    refinements=[]
    for coarse,fine in zip(runs[::2],runs[1::2]):
        rel=abs(coarse['static_coefficient']/fine['static_coefficient']-1)
        assert rel<.004
        refinements.append(rel)
    receipt={"date":"2026-10-13","scope":"retarded Gaussian charge and axial boundary-transfer response on fixed RN-AdS backgrounds",
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "dependency_sha256":hashlib.sha256((ROOT/'checks/angular_boundary_response.py').read_bytes()).hexdigest(),
        "symbolic":alg,"cubic_memory_envelope":envelope,
        "passive_finite_regulators":runs,"static_grid_refinement_relative":refinements,
        "single_species_candidates":candidates(),
        "limitations":["Finite regulators have reflecting horizon cuts; the half-line theorem is proved separately",
            "The prescribed phase/chemical drive is not the actual nonlinear fermion-transfer switch",
            "Quantum transfer covariance, noise, moving geometry and other sectors remain unmatched"],"passed":True}
    (ROOT/'evidence/causal-compact-response.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: constrained master operator, 35 cubic integrals, 12 passive regulators, 6 physical scale bounds')


if __name__=='__main__':
    main()
