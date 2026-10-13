"""S072: finite-q quadratic electric feedback with quantum noise.

Checks the static common subtraction, oscillator realization, positive-real
response, switching stiffness and finite-band symplectic quantum evolution.
No all-order junction or renormalized four-dimensional stress is computed.
"""
from pathlib import Path
import hashlib
import json
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp, simpson

ROOT=Path(__file__).resolve().parents[1]


def spectra(q,limit):
    def coeff(a):
        f=np.zeros(limit+1);f[0]=1.;f[1]=4*a
        for e in range(2,limit+1):f[e]=(4*a*f[e-1]+(e-2)*f[e-2])/e
        return f
    fa=coeff(1-1/q);fc=coeff(1/q);e=np.arange(limit+1)
    v=q/(4*np.pi**2)*(fa+(q-1)*fc);w=q/(4*np.pi**2)*(fa-(q-1)*fc)
    v[(e==0)|(e%2!=0)]=0;w[e%2==0]=0;w[1]=0
    return v,w


def static_delta(q):
    with mp.workdps(75):
        Q=mp.mpf(q)
        return Q/mp.pi**2*(Q/2*(mp.digamma(mp.mpf('.5'))-mp.digamma(mp.mpf('.5')-1/Q))-Q/(Q-2))


def subtraction_checks():
    rows=[];mp.mp.dps=70
    for q in [3,5,31,127,1000,2000**4]:
        Q=mp.mpf(q);r=Q/(Q-2)
        def integrand(y):
            if y==0:return r*(Q-1)
            if y==1:return mp.mpf(0)
            ly=mp.log(y)
            return r*(-1+Q*mp.expm1((r-1)*ly)/mp.expm1(2*r*ly))
        direct=Q/mp.pi**2*mp.quad(integrand,[0,.25,.5,.9,1])
        exact=static_delta(q);err=abs(direct/exact-1)
        assert exact>0 and err<mp.mpf('1e-40')
        v,w=spectra(q,1025)
        # Equal numbers of even and odd modes: this is deliberately not
        # accepted as the continuum static subtraction at large q.
        even=np.arange(2,1025,2);odd=even+1
        naive=np.sum(2*v[even]/even-2*w[odd]/odd)
        rows.append(dict(q=q,delta_chi=mp.nstr(exact,30),ratio_to_mu=mp.nstr(exact/(Q/(2*mp.pi)),30),
                         integral_relative_error=mp.nstr(err,8),paired_512_mode_ratio=float(naive/(q/(2*np.pi)))))
    assert abs(rows[-1]['paired_512_mode_ratio']-float(static_delta(2000**4)/(mp.mpf(2000**4)/(2*mp.pi))))>.6
    return rows


def hdata(t,hf,T):
    if t<=0:return 0.,0.
    if t>=T:return hf,0.
    x=t/T
    return hf*(3*x*x-2*x**3),hf*6*x*(1-x)/T


def cosine_drive(e,t,hf,T):
    if t<=T:
        z=hf*(6/T**2*(1-np.cos(e*t))/e**2-12/T**3*(t/e**2-np.sin(e*t)/e**3))
        zd=hf*(6/T**2*np.sin(e*t)/e-12/T**3*(1-np.cos(e*t))/e**2)
    else:
        y=e*T;f=12*(2*np.sin(y/2)-y*np.cos(y/2))/y**3
        z=hf*f*np.cos(e*(t-T/2));zd=-hf*f*e*np.sin(e*(t-T/2))
    return z,zd


def switch_checks():
    rows=[]
    for q,hf,T in [(31,.03,.5),(127,.08,.6),(1000,.08,.5)]:
        limit=16384;v,_=spectra(q,limit);e=np.arange(2,limit+1,2,dtype=float)
        a=2*v[e.astype(int)]/e;mu=q/(2*np.pi);delta=float(static_delta(q));A=24*hf/T**2
        error_bound=16*float(mp.zeta(3))/np.pi**2*A
        trunc_bound=q*A/(np.pi**2*limit)
        ratios=[];errors=[]
        times=np.r_[np.linspace(.002,2.9,70),np.linspace(3.2,10*np.pi+.8,90),
                    np.arange(1,11)*np.pi+T/2]
        for t in times:
            h,hd=hdata(t,hf,T);z,_=cosine_drive(e,t,hf,T)
            returned_derivative=hd if t<np.pi else 2*hdata(t%np.pi,hf,T)[1]
            R=float(a@z);Rinf=mu*returned_derivative-2*mu*h/np.pi
            assert abs(R-Rinf)<=error_bound+trunc_bound
            stiffness=mu*h+h*h*delta-h*R
            ratios.append(stiffness/(mu*h));errors.append(abs(R-Rinf))
        rows.append(dict(q=q,hf=hf,T=T,minimum_stiffness_ratio=float(min(ratios)),
                         largest_R_comparison_error=float(max(errors)),R_error_bound=error_bound,
                         R_truncation_bound=trunc_bound,
                         pre_return_stiffness_lower_ratio=1-1.5*hf/T-error_bound/mu,
                         all_time_stiffness_lower_ratio=1-3*hf/T-error_bound/mu))
    candidate=[]
    for k in [2000,4000,8000]:
        q=k**4;hf=4*np.pi/k;T=.25
        eps=768*float(mp.zeta(3))/np.pi*hf/(q*T*T)
        lower=1-1.5*hf/T-eps
        assert lower>.962
        candidate.append(dict(k=k,T=T,maximum_h_derivative=1.5*hf/T,
                              finite_q_error=eps,pre_return_stiffness_lower_ratio=lower,
                              all_time_stiffness_lower_ratio=1-3*hf/T-eps))
    return dict(reference=rows,candidate=candidate)


def oscillator_model(q,h,I,number=24):
    # Two retained radial modes illustrate the independent positive G_alpha.
    wr=np.array([.7,2.1]);sig=np.array([.04,.01])
    e=np.arange(3,2*number+2,2,dtype=float)
    _,w=spectra(q,int(e[-1]));c=np.sqrt(2*e*w[e.astype(int)])
    n=3+number;alpha=np.zeros(n);alpha[0]=1;alpha[1:3]=np.sqrt(sig)
    base=np.diag(np.r_[0.,wr*wr,e*e]);inv_mass=np.r_[1/I,np.ones(n-1)]
    delta=float(static_delta(q));mu=q/(2*np.pi);kappa=mu*h+h*h*delta
    K=base+(kappa+h*h*np.sum(c*c/(e*e)))*np.outer(alpha,alpha)
    for j,cj in enumerate(c):
        K[:,j+3]-=h*cj*alpha;K[j+3,:]-=h*cj*alpha
    return K,inv_mass,alpha,e,c,wr,sig,kappa


def resolvent_checks():
    rows=[]
    for q,h,I in [(3,.08,1.7),(31,.12,.8),(1000,.08,2.)]:
        K,D,a,e,c,wr,sig,kappa=oscillator_model(q,h,I,32)
        eig=np.linalg.eigvalsh(np.sqrt(D)[:,None]*K*np.sqrt(D)[None,:])
        assert min(eig)>0
        errors=[];real_parts=[]
        for re in [.03,.2,1.,5.]:
            for im in [-30.,-3.1,0.,.7,9.]:
                s=re+1j*im
                G=1/(I*s*s)+np.sum(sig/(s*s+wr*wr))
                impedance=1/(s*G)+kappa/s+h*h*np.sum(c*c/(e*e)*s/(s*s+e*e))
                direct=a@np.linalg.solve(s*s*np.diag(1/D)+K,a)
                errors.append(abs(direct*s*impedance-1));real_parts.append(impedance.real)
                assert impedance.real>0
        assert max(errors)<2e-10
        # The omitted infinite neutral tower has a convergent impedance tail.
        _,w=spectra(q,8193);ee=np.arange(3,8194,2,dtype=float)
        for s in [.2+3j,1.+8j]:
            L=255.;keep=ee>L
            tail=h*h*np.sum(2*w[ee[keep].astype(int)]/ee[keep]*s/(s*s+ee[keep]**2))
            bound=4*h*h*q*abs(s)/(3*np.pi**2*L)
            assert abs(tail)<=bound
        rows.append(dict(q=q,minimum_normal_frequency_squared=float(min(eig)),
                         maximum_resolvent_error=float(max(errors)),minimum_impedance_real_part=float(min(real_parts))))
    return rows


def gaussian_switch(q,hf,T,I,rtol=2e-10,with_floquet=False):
    number=12;nx=128
    _,D,a,e,c,wr,sig,_=oscillator_model(q,hf,I,number)
    n=len(a);mu=q/(2*np.pi);delta=float(static_delta(q))
    vx,_=spectra(q,2*nx);ex=np.arange(2,2*nx+1,2,dtype=float);ax=2*vx[ex.astype(int)]/ex
    base=np.diag(np.r_[0.,wr*wr,e*e]);aa=np.outer(a,a);sumct=np.sum(c*c/(e*e))
    cross=np.zeros((n,n))
    for j,cj in enumerate(c):cross[:,j+3]-=cj*a;cross[j+3,:]-=cj*a
    def matrices(t):
        h,hd=hdata(t,hf,T);z,zd=cosine_drive(ex,t,hf,T)
        R=ax@z;Rd=ax@zd;k=mu*h+h*h*delta-h*R
        kd=mu*hd+2*h*hd*delta-hd*R-h*Rd
        K=base+(k+h*h*sumct)*aa+h*cross
        Kd=(kd+2*h*hd*sumct)*aa+hd*cross
        return K,Kd,k
    frequencies=np.r_[wr,e]
    initial=np.r_[.5,1/(2*frequencies),.5,frequencies/2]
    def rhs(t,y):
        S=y.reshape(2*n,2*n);K,_,_=matrices(t)
        return np.vstack((D[:,None]*S[n:],-K@S[:n])).ravel()
    stop=1.4
    sol=solve_ivp(rhs,(0,stop),np.eye(2*n).ravel(),method='DOP853',rtol=rtol,atol=rtol/50,
                  max_step=.003,dense_output=True)
    assert sol.success
    J=np.block([[np.zeros((n,n)),np.eye(n)],[-np.eye(n),np.zeros((n,n))]])
    output=[];full=[];bare=[];symplectic=[]
    times=np.linspace(0,stop,101)
    for t,y in zip(times,sol.sol(times).T):
        S=y.reshape(2*n,2*n);A=a@S[:n];P=S[n]
        va=np.sum(A*A*initial);vp=np.sum(P*P*initial);cov=np.sum(A*P*initial)
        full.append(va*vp-cov*cov)
        indices=np.array([0,n]);v0=initial[indices]
        ba=A[indices];bp=P[indices]
        bare.append(np.sum(ba*ba*v0)*np.sum(bp*bp*v0)-np.sum(ba*bp*v0)**2)
        symplectic.append(np.max(abs(S@J@S.T-J)))
        output.append([float(va),float(vp),float(cov)])
    assert min(full)>=.25-2e-8
    assert min(bare)<.245
    assert max(symplectic)<2e-8
    def energy(t):
        S=sol.sol(t).reshape(2*n,2*n);K,_,_=matrices(t)
        qq=(S[:n]*initial)@S[:n].T
        pp=(S[n:]*initial)@S[n:].T
        return .5*(np.sum(D*np.diag(pp))+np.sum(K*qq.T))
    works=[]
    for points in [2001,4001]:
        ts=np.linspace(0,stop,points);integrand=[]
        for t,y in zip(ts,sol.sol(ts).T):
            Sq=y.reshape(2*n,2*n)[:n];_,Kd,_=matrices(t)
            integrand.append(.5*np.sum(Kd*((Sq*initial)@Sq.T).T))
        works.append(float(simpson(integrand,x=ts)))
    change=float(energy(stop)-energy(0));err=abs(change-works[-1])
    assert err<2e-5
    result=dict(q=q,hf=hf,T=T,I=I,neutral_bath_modes=number,radial_modes=2,X_response_modes=nx,
                minimum_full_uncertainty_determinant=float(min(full)),
                minimum_homogeneous_only_determinant=float(min(bare)),maximum_symplectic_error=float(max(symplectic)),
                final_covariance=output[-1],finite_band_energy_change=change,integrated_switch_work=works,
                energy_balance_error=err,rtol=rtol)
    if with_floquet:
        start=T+.07
        periodic=solve_ivp(rhs,(start,start+np.pi),np.eye(2*n).ravel(),method='DOP853',
                           rtol=rtol,atol=rtol/50,max_step=.003)
        assert periodic.success
        monodromy=periodic.y[:,-1].reshape(2*n,2*n)
        radius=float(max(abs(np.linalg.eigvals(monodromy))))
        result['floquet']=dict(period=float(np.pi),largest_multiplier_modulus=radius,
                                inferred_growth_rate=float(np.log(radius)/np.pi),
                                symplectic_error=float(np.max(abs(monodromy@J@monodromy.T-J))))
    return result


def main():
    switches=[gaussian_switch(127,.08,.6,.7),gaussian_switch(1000,.08,.5,2.,with_floquet=True)]
    fine=gaussian_switch(1000,.08,.5,2.,2e-12,with_floquet=True)
    refinement=max(abs(np.array(switches[-1]['final_covariance'])-fine['final_covariance']))
    assert refinement<2e-7
    result=dict(date='2026-10-13',status='PASS',python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                subtraction=subtraction_checks(),switch_stiffness=switch_checks(),resolvent=resolvent_checks(),
                quantum_switches=switches,quantum_refinement=fine,covariance_refinement_error=float(refinement),
                limitations=['Phase-linear and second-order neutral influence, resummed as a Gaussian model; angular sector held fixed',
                             'Stationary no-growing-pole result and all-time switched stiffness bound are different statements',
                             'Finite-band covariance and energy only; no renormalized continuum stress, nonlinear stability or physical capture'])
    (ROOT/'evidence/passive-electric-feedback.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: finite-q subtraction, positive-real feedback, switching stiffness, quantum noise and finite-band energy balance.')
    print('Candidate:',result['switch_stiffness']['candidate'][0])
    print('Uncertainty pairs:',[(x['minimum_full_uncertainty_determinant'],x['minimum_homogeneous_only_determinant']) for x in switches])
    print('Covariance refinement:',refinement)
    print('Floquet refinement:',switches[-1]['floquet'],fine['floquet'])


if __name__=='__main__':main()
