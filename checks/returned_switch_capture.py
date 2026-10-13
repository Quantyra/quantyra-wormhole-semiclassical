"""S072: exact returned free-fermion stress and a longer JT opening.

No gauged four-dimensional state or physical creation protocol is computed.
"""
from pathlib import Path
import hashlib
import json
import platform
import warnings

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp, IntegrationWarning
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import Polynomial
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
warnings.simplefilter('error', IntegrationWarning)
TAU=.25
STEP=Polynomial([0,0,0,10,-15,6])
ISQ=(STEP*STEP).integ()
GLX,GLW=leggauss(24)


def S(x):
    return STEP(np.clip(x,0.,1.))


def primitive(t,tau=TAU,squared=False):
    tt=np.asarray(t)
    x=np.clip(tt/tau,0.,1.)
    if squared:
        val=tau*ISQ(x)
        endpoint=tau*ISQ(1.)
    else:
        val=tau*(2.5*x**4-3*x**5+x**6)
        endpoint=tau/2
    return np.where(tt<=0,0.,np.where(tt<tau,val,endpoint+tt-tau))


def drive(t,eta):
    a=(primitive(t)-primitive(np.asarray(t)-np.pi))/np.pi
    r=4*eta/np.pi*(primitive(t,squared=True)-primitive(np.asarray(t)-np.pi,squared=True))
    ap=(S(np.asarray(t)/TAU)-S((np.asarray(t)-np.pi)/TAU))/np.pi
    rp=4*eta/np.pi*(S(np.asarray(t)/TAU)**2-S((np.asarray(t)-np.pi)/TAU)**2)
    return a,r,ap,rp


def hnorm(t,eta):
    v=S(np.asarray(t)/TAU)
    return (-v+4*eta*v*v)/(2*np.pi)


def H(t,eta):
    x=np.asarray(t)
    end=np.clip(x,0.,TAU)
    u=end[...,None]*(GLX+1)/2
    weights=end[...,None]*GLW/2
    src=-hnorm(u,eta)
    p=np.sum(weights*src,axis=-1)
    c=np.sum(weights*src*np.cos(u),axis=-1)
    s=np.sum(weights*src*np.sin(u),axis=-1)
    value=p-np.cos(x)*c-np.sin(x)*s
    deriv=np.sin(x)*c-np.cos(x)*s
    tail=np.maximum(x-TAU,0.)
    const=(1-4*eta)/(2*np.pi)
    return value+const*(tail-np.sin(tail)),deriv+const*(1-np.cos(tail))


def boundary(t,eta,z0,v0):
    return z0*np.cos(t)+v0*np.sin(t)+H(t,eta)[0]-H(np.asarray(t)-np.pi,eta)[0]


def bulk(t,x,eta,z0,v0):
    hm,hmd=H(t-x,eta)
    hp,hpd=H(t+x-np.pi,eta)
    return (z0*np.cos(t)+v0*np.sin(t))/np.sin(x)+hmd+hpd+(hm-hp)/np.tan(x)


def symbolic_proof_checks():
    R=sp.Rational
    # pi>25/8, sqrt(3)<7/4 imply atan(14pi/25)>pi/3.
    assert R(14,25)*R(25,8)==R(7,4) and R(7,4)**2>3
    assert R(1,3)-R(7,25)==R(4,75)
    ramp_lower=R(4,75)-R(65,1600)-R(65,32000)-R(1,2000)
    assert ramp_lower>R(1,100)
    assert (R(7,25)**2+2*R(7,25)*R(1,8)*R(8,25)+R(8,25)**2)==R(127,625)
    assert R(113,10)**2>127
    amplitude=R(113,250)+R(1,1000)+R(1,2000)
    center=R(499,1000)
    assert amplitude==R(907,2000)
    assert center-amplitude>R(1,25)
    assert amplitude**2-center**2 < -R(1,25)
    assert R(7,25)+R(1,2000)+R(65,64)<R(13,10)
    x=sp.symbols('x',real=True)
    step=10*x**3-15*x**4+6*x**5
    assert sp.integrate(step,(x,0,1))==R(1,2)
    return dict(ramp_lower=str(ramp_lower),late_amplitude_upper=str(amplitude),
                late_redshift_lower=str(center-amplitude),energy_upper=str(amplitude**2-center**2),
                uniform_redshift_lower='1/100',uniform_redshift_upper='13/10',area_lower='-49/100')


def point_split_checks():
    mp.mp.dps=65
    J=mp.matrix([[0,-1],[1,0]])
    def rot(x):
        return mp.eye(2)*mp.cos(x)+J*mp.sin(x)
    def step(x):
        x=max(mp.mpf(0),min(mp.mpf(1),x))
        return 10*x**3-15*x**4+6*x**5
    def prim(t):
        tau=mp.mpf('.25')
        if t<=0:return mp.mpf(0)
        if t>=tau:return t-tau/2
        x=t/tau
        return tau*(mp.mpf('2.5')*x**4-3*x**5+x**6)
    rows=[]
    for eta in [mp.mpf('0.0005'),mp.mpf('.02')]:
        theta=lambda t:4*eta*(prim(t)-prim(t-mp.pi))
        accum=lambda t:sum(theta(t-j*mp.pi) for j in range(max(0,int(mp.floor(t/mp.pi))+1)))
        for t in [mp.mpf('.12'),mp.mpf('.8'),mp.pi+mp.mpf('.12'),mp.pi+mp.mpf('.4'),5*mp.pi+mp.mpf('.1')]:
            D=4*eta*step(t/mp.mpf('.25'))
            Dprev=4*eta*step((t-mp.pi)/mp.mpf('.25'))
            errors=[]
            for eps in [mp.mpf('1e-3'),mp.mpf('1e-4'),mp.mpf('1e-5')]:
                def covariance_trace(d,offset):
                    G=mp.eye(2)/(4*mp.pi*mp.j*mp.sin(d/2))+mp.j*J/(4*mp.pi*mp.cos(d/2))
                    diff=rot(accum(t+d/2-offset))*G*rot(accum(t-d/2-offset)).T-G
                    return diff[0,0]+diff[1,1]
                out=-mp.j*(covariance_trace(eps,0)-covariance_trace(-eps,0))/(4*eps)
                inc=-mp.j*(covariance_trace(eps,mp.pi)-covariance_trace(-eps,mp.pi))/(4*eps)
                predicted=(D*D-D-Dprev*Dprev+Dprev)/(4*mp.pi)
                errors.append(abs(out-inc-predicted))
            assert errors[-1]<mp.mpf('1e-9')
            assert abs(accum(t)-4*eta*prim(t))<mp.mpf('1e-60')
            # Closed spectral sum at a finite imaginary time: both towers.
            d=mp.mpf('.21')-mp.j*mp.mpf('.3')
            G=mp.eye(2)/(4*mp.pi*mp.j*mp.sin(d/2))+mp.j*J/(4*mp.pi*mp.cos(d/2))
            exact=rot(4*eta*d)*G
            cp=mp.exp(-mp.j*(mp.mpf('1.5')-4*eta)*d)/(mp.pi*(1-mp.exp(-2*mp.j*d)))
            cm=mp.exp(-mp.j*(mp.mpf('.5')+4*eta)*d)/(mp.pi*(1-mp.exp(-2*mp.j*d)))
            spectral=(cp+cm)/2*mp.eye(2)+(cp-cm)/(2*mp.j)*J
            spectrum_error=max(abs(exact[i,j]-spectral[i,j]) for i in range(2) for j in range(2))
            assert spectrum_error<mp.mpf('1e-60')
            rows.append(dict(eta=float(eta),time=float(t),stress_errors=list(map(float,errors)),
                             spectral_covariance_error=float(spectrum_error)))
    return rows


def quadrature_piecewise(fun,a,b,cuts,**kwargs):
    split=sorted(set([a,b]+[x for x in cuts if a<x<b]))
    return sum(quad(fun,l,r,epsabs=2e-12,epsrel=2e-12,**kwargs)[0] for l,r in zip(split[:-1],split[1:]))


def work_checks():
    rows=[]
    for shape,T in [('cubic',.25),('cubic',4.),('shaped',np.pi+TAU)]:
        def f(t):
            if shape=='shaped':return float(drive(t,0.)[2])
            if 0<t<T:return 6*t*(T-t)/T**3
            return 0.
        cuts=[0.,TAU,np.pi,np.pi+TAU,T]
        def period(t):
            return sum(f(t+j*np.pi) for j in range(1+int(T/np.pi)))
        Q=quadrature_piecewise(lambda t:period(t)**2,0.,np.pi,[TAU,T,T-np.pi])
        direct=quadrature_piecewise(lambda t:f(t)*(f(t)+2*sum(f(t-j*np.pi) for j in range(1,1+int(T/np.pi)))),0.,T,cuts)
        assert abs(direct-Q)<2e-12
        excess=Q-1/np.pi
        assert excess>=-2e-12
        parseval=0.
        for n in range(1,257):
            re=quadrature_piecewise(f,0.,T,cuts,weight='cos',wvar=2*n)
            im=quadrature_piecewise(f,0.,T,cuts,weight='sin',wvar=2*n)
            parseval+=2/np.pi*(re*re+im*im)
        # Cubic transform bounded by 12/(T^2 E^2)+24/(T^3 E^3).
        tail=0. if shape=='shaped' else (2/np.pi)*(3/T**2+3/(256*T**3))**2/(3*256**3)
        assert abs(excess-parseval)<=tail+3e-11
        if shape=='shaped':assert abs(excess)<2e-12
        rows.append(dict(shape=shape,duration=T,full_work_coefficient=Q,
                         excess_coefficient=excess,parseval_256= parseval,tail_bound=tail,
                         naive_pre_return_excess=quadrature_piecewise(lambda t:f(t)**2,0,T,cuts)-1/np.pi))
    return rows


def geometry_case(eta,z0,v0,tight=False):
    finish=np.pi+TAU
    def rhs(t,Y):
        a,r,_,_=drive(t,eta)
        return [Y[1],-Y[0]+(a-r)/2]
    sol=solve_ivp(rhs,(0,finish),[z0,v0],method='DOP853',rtol=2e-13 if tight else 2e-11,
                  atol=2e-15 if tight else 2e-13,max_step=.005 if tight else .01,dense_output=True)
    assert sol.success
    ts=np.linspace(0,finish,501)
    minimum=float(np.min(sol.sol(ts)[0]))
    assert minimum> .01
    error=float(np.max(abs(sol.sol(ts)[0]-boundary(ts,eta,z0,v0))))
    assert error<1e-10
    zf,vf=sol.y[:,-1]
    A=1-4*eta
    ef=vf*vf+zf*zf-A*zf
    amp=np.hypot(zf-A/2,vf)
    assert ef<-.04 and A/2-amp>.04
    # Independent original nonlinear equations in s=eta*u.
    def physical(s,Y):
        z,w,t=Y
        a,r,_,_=drive(t,eta)
        return [z*w,-z*z+(a-r)*z/2,z]
    def stop(s,Y):return Y[2]-finish
    stop.terminal=True
    nonlinear=solve_ivp(physical,(0,400),[z0,v0,0.],events=stop,method='DOP853',
                        rtol=2e-12,atol=2e-14,max_step=.08)
    assert nonlinear.success and len(nonlinear.t_events[0])==1
    match=float(max(abs(nonlinear.y[:2,-1]-sol.y[:,-1])))
    assert match<2e-9
    times=np.linspace(0,finish+4*np.pi,89)
    xs=np.linspace(.03,np.pi-.03,65)
    tt,xx=np.meshgrid(times,xs,indexing='ij')
    values=bulk(tt,xx,eta,z0,v0)
    assert values.min()>-.49
    reflected=float(np.max(abs(values-bulk(tt,np.pi-xx,eta,z0,v0))))
    assert reflected<2e-11
    late_t=np.linspace(finish,finish+2*np.pi,31)[:,None]
    periodic=float(np.max(abs(bulk(late_t,xs,eta,z0,v0)-bulk(late_t+2*np.pi,xs,eta,z0,v0))))
    assert periodic<2e-11
    spatial=[]
    for t in [.5,np.pi/2+.02,np.pi+.12,finish+.1,finish+4.]:
        B=boundary(t,eta,z0,v0)
        for x in [.17,.63,np.pi/2]:
            def f(r):
                k=-hnorm(t-r,eta)-hnorm(t+r-np.pi,eta)
                return k*(1-np.cos(min(x,r)))
            cuts=[x,t,t-TAU,np.pi-t,np.pi+TAU-t]
            W=B-quadrature_piecewise(f,0,np.pi/2,cuts)
            spatial.append(abs(W/np.sin(x)-bulk(t,x,eta,z0,v0)))
    assert max(spatial)<2e-11
    return dict(eta=eta,z0=z0,v0=v0,tight=tight,max_step=.005 if tight else .01,minimum_ramp_redshift=minimum,
                final_energy=ef,late_redshift_lower=A/2-amp,
                final_state=[float(zf),float(vf)],boundary_reconstruction_error=error,
                nonlinear_endpoint_error=match,duration_eta_u=float(nonlinear.t[-1]),
                minimum_sampled_dilaton=float(values.min()),reflection_error=reflected,
                late_periodicity_error=periodic,spatial_reconstruction_error=float(max(spatial)))


def main():
    symbolic=symbolic_proof_checks()
    stress=point_split_checks()
    work=work_checks()
    cases=[geometry_case(eta,z0,v0) for eta in [1e-5,.0005]
           for z0 in [.25,.265,.28] for v0 in [-.0005,0.,.0005]]
    fine=geometry_case(.0005,.28,-.0005,True)
    ref=next(x for x in cases if x['eta']==.0005 and x['z0']==.28 and x['v0']==-.0005)
    refinement=float(max(abs(np.array(fine['final_state'])-ref['final_state'])))
    assert refinement<1e-10
    result=dict(date='2026-10-13',status='PASS',python=platform.python_version(),numpy=np.__version__,
                scipy=scipy.__version__,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                symbolic=symbolic,point_split=stress,work=work,geometry=cases,refined=fine,
                final_state_refinement=refinement,
                scope=['Exact free-Majorana nondispersive junction and symmetric semiclassical JT constraints',
                       'All-mode returned covariance and renormalized stress; final matter vacuum, oscillating geometry',
                       'Prepared pair and imposed nonlocal boundary interaction; no physical preparation or four-dimensional matching',
                       'Electric-neutral h and free-junction theta retain distinct physical normalization gates'])
    (ROOT/'evidence/returned-switch-capture.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: returned covariance, exact final vacuum, capture with preparation margin and conditional area bound.')
    print('Largest final energy:',max(x['final_energy'] for x in cases))
    print('Smallest ramp/late redshift:',min(x['minimum_ramp_redshift'] for x in cases),min(x['late_redshift_lower'] for x in cases))
    print('Refinement:',refinement)
    print('Work:',work)


if __name__=='__main__':main()
