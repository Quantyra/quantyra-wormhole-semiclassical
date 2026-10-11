"""Independent checks of transparent JT/channel mean response, not 4D closure."""
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import brentq
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
NOTE=ROOT/'notes/transparent-clock-flux.md'
OUT=ROOT/'evidence/transparent-clock-flux.json'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def exact_checks():
    x,A,c,D,q=sy.symbols('x A c D q',positive=True)
    T=sy.pi/x;L=T+D
    energy=A*x*x+q*c*x/24-q*c*sy.pi/(6*L)
    derivative=x*sy.diff(energy,x)
    stiffness=x*sy.diff(derivative,x)
    checks={
        'static_stiffness_at_stationarity':sy.simplify(stiffness-derivative-(2*A*x*x+q*c*sy.pi*T*D/(3*L**3))),
        'total_boundary_plus_exterior_energy':sy.simplify(energy-derivative-(-A*x*x-q*c*sy.pi*D/(6*L**2))),
    }
    h,r1,r2,r3=sy.symbols('h r1 r2 r3',positive=True)
    rp=r1+r2*h/2+r3*h*h/8
    rm=r1-r2*h/2+r3*h*h/8
    diff=r1*h+r3*h**3/24
    coeff=sy.expand(sy.series(sy.sqrt(rp*rm)/diff,h,0,2).removeO()).coeff(h,1)
    checks['point_split_schwarzian']=sy.simplify(coeff-(r3/r1-sy.Rational(3,2)*(r2/r1)**2)/12)
    k,t=sy.symbols('k t',positive=True)
    p=sy.exp(k*t)
    sch=sy.diff(p,t,3)/sy.diff(p,t)-sy.Rational(3,2)*(sy.diff(p,t,2)/sy.diff(p,t))**2
    checks['positive_hawking_sign']=sy.simplify(-sch/(24*sy.pi)-k*k/(48*sy.pi))
    for name,value in checks.items():assert value==0,(name,value)
    return [dict(identity=name,residual=str(value)) for name,value in checks.items()]


def point_split():
    mp.mp.dps=55
    rows=[]
    for L,t in [(mp.mpf('3.7'),mp.mpf('.4')),(mp.mpf('8.2'),mp.mpf('1.3'))]:
        R=lambda v:v+mp.mpf('.12')*mp.sin(mp.mpf('.9')*v)
        rp=mp.diff(R,t);rpp=mp.diff(R,t,2);rppp=mp.diff(R,t,3)
        expected=-mp.pi*rp**2/(12*L**2)-(rppp/rp-mp.mpf('1.5')*(rpp/rp)**2)/(24*mp.pi)
        values=[]
        for h in [mp.mpf('.004'),mp.mpf('.002'),mp.mpf('.001')]:
            gap=R(t+h/2)-R(t-h/2)
            kernel=mp.sqrt(mp.diff(R,t+h/2)*mp.diff(R,t-h/2))*mp.pi/(L*mp.sin(mp.pi*gap/L))
            value=-(kernel-1/h)/(2*mp.pi*h)
            values.append(dict(h=float(h),stress=float(value),error=float(abs(value-expected))))
        rich=(4*values[-1]['stress']-values[-2]['stress'])/3
        assert abs(rich-float(expected))<2e-12
        rows.append(dict(length=float(L),t=float(t),expected=float(expected),refinements=values,richardson_error=abs(rich-float(expected))))
    return rows


tt=sy.symbols('t',real=True)
bump=(1-sy.cos(2*sy.pi*tt/8))**4
bump_derivatives=[sy.lambdify(tt,sy.diff(bump,tt,n),'numpy') for n in range(4)]


def pulse(t,n):
    return float(bump_derivatives[n](t)) if 0<t<8 else 0.0


def clock(t,amp,n=0):
    return (t if n==0 else 1.0 if n==1 else 0.0)+amp*pulse(t,n)


def ray(t,amp):
    p=brentq(lambda s:clock(s,amp)-clock(t,amp)+math.pi,t-5,t,xtol=8e-14)
    pt,pp=clock(t,amp,1),clock(p,amp,1)
    d1=pt/pp
    d2=(clock(t,amp,2)-clock(p,amp,2)*d1*d1)/pp
    d3=(clock(t,amp,3)-clock(p,amp,3)*d1**3-3*clock(p,amp,2)*d1*d2)/pp
    return p,d1,d2,d3


def state(t,D,amp):
    L=math.pi+D;initial=-math.pi/(12*L*L)
    s=t;steps=[]
    while s>0:
        p,d1,d2,d3=ray(s,amp)
        steps.append((d1,d2,d3));s=p-D
        assert len(steps)<100
    R=s+len(steps)*L
    a1,a2,a3=1.0,0.0,0.0
    stress=initial
    for d1,d2,d3 in reversed(steps):
        stress=d1*d1*stress-(d3/d1-1.5*(d2/d1)**2)/(24*math.pi)
        a1,a2,a3=a1*d1,a2*d1*d1+a1*d2,a3*d1**3+3*a2*d1*d2+a1*d3
    direct=initial*a1*a1-(a3/a1-1.5*(a2/a1)**2)/(24*math.pi)
    assert abs(stress-direct)<3e-13
    assert a1>0
    return R,a1,a2,a3,stress


def linear_stress(t,D):
    L=math.pi+D;e0=-math.pi/(12*L*L);out=0.0
    while t>0:
        out+=2*e0*(pulse(t,1)-pulse(t-math.pi,1))-(pulse(t,3)-pulse(t-math.pi,3))/(24*math.pi)
        t-=L
    return out


def causal_cases():
    rows=[]
    times=np.linspace(.21,27.3,83)
    for D in [0,.3,2.7,8.0]:
        target=np.array([linear_stress(t,D) for t in times])
        initial=-math.pi/(12*(math.pi+D)**2)
        refinements=[]
        for amp in [.001,.0005,.00025]:
            data=[state(t,D,amp) for t in times]
            if D==0:
                for t,v in zip(times,data):
                    for n in range(4):assert abs(v[n]-clock(t,amp,n))<2e-11
            response=np.array([(v[-1]-initial)/amp for v in data])
            error=float(np.max(np.abs(response-target)))
            refinements.append(dict(amplitude=amp,error=error,minimum_phase_derivative=min(v[1] for v in data)))
        ratios=[old['error']/new['error'] for old,new in zip(refinements,refinements[1:])]
        assert all(1.8<r<2.2 for r in ratios),ratios
        # Check the composed phase derivatives against actual ray inversion.
        t=11.31;amp=.0005;R,d1,d2,d3,_=state(t,D,amp)
        fd=[]
        for h in [.02,.01,.005]:
            vals=[state(t+k*h,D,amp)[0] for k in [-2,-1,0,1,2]]
            v3=(vals[4]-2*vals[3]+2*vals[1]-vals[0])/(2*h**3)
            fd.append(dict(h=h,third_derivative=v3,error=abs(v3-d3)))
        assert fd[-1]['error']<2e-6
        rows.append(dict(return_length=D,refinements=refinements,error_ratios=ratios,
                         phase_third_derivative=d3,phase_finite_difference=fd,
                         scope='prescribed clock, causal free state; not a solution of backreaction'))
    return rows


def transfer(w,T,D,eps=0.0):
    if eps:
        x=math.pi/T;b=math.sqrt(1-eps*eps*x*x)
        Te=2*math.acos(eps*x)/x
        B=-np.expm1(1j*w*Te)-1j*w*eps/b*(1+np.exp(1j*w*Te))
    else:Te=T;B=-np.expm1(1j*w*T)
    L=Te+D
    return (w*w-(2*math.pi/L)**2)*B/(-np.expm1(1j*w*L))/(24*math.pi)


def matter(w,T,D,eps=0.0):
    return 2j*np.expm1(1j*w*D)*transfer(w,T,D,eps)/w


def spectral(w,T,D,N):
    L=T+D;n=np.arange(2,N+1,dtype=float);E=2*np.pi*n/L
    W=2*(n*n-1)*np.sin(np.pi*n*T/L)**2/(3*L*L*n)
    value=np.sum(2*w*w*W/(E*(w*w-E*E)))
    assert 2*np.pi*N/L>2*abs(w)
    bound=2*abs(w)**2*L/(9*np.pi**3*N)
    return complex(value),float(bound)


def positive_spectral_factor(w,T,D,N):
    L=T+D;n=np.arange(2,N+1,dtype=float);E=2*np.pi*n/L
    W=2*(n*n-1)*np.sin(np.pi*n*T/L)**2/(3*L*L*n)
    return float(np.sum(2*W/(E*np.abs(w*w-E*E)**2)))


def response_cases():
    rows=[];equilibria=[];axis=[]
    for Ds in [[0,0,0],[.1,.3,2.7],[3,10,100],[.01,1000,.5]]:
        A=len(Ds)/16
        f=lambda x:2*A*x+len(Ds)/24-sum(np.pi**2/(6*(np.pi+D*x)**2) for D in Ds)
        x=brentq(f,0,1.01);T=math.pi/x
        K0=2*A*x*x+sum(math.pi*T*D/(3*(T+D)**3) for D in Ds)
        mp.mp.dps=40
        def es(phi):
            y=mp.mpf(str(x))*mp.exp(phi)
            return mp.mpf(str(A))*y*y+len(Ds)*y/24-sum(mp.pi/(6*(mp.pi/y+D)) for D in Ds)
        numerical_stiffness=float(mp.diff(es,0,2))
        boundary=-A*x*x-sum(math.pi*D/(6*(T+D)**2) for D in Ds)
        assert abs(numerical_stiffness-K0)<2e-11 and abs(float(es(0))-boundary)<2e-11
        equilibria.append(dict(returns=Ds,clock_rate=x,stationarity=abs(f(x)),stiffness=K0,
                              independent_stiffness=numerical_stiffness,energy_identity_error=abs(float(es(0))-boundary)))
        for wrel in [.7+.2j,1+.1j,2.3+.3j,-1+.2j]:
            w=wrel*x
            K=-2*A*w*w+2*A*x*x+sum(matter(w,T,D) for D in Ds)
            trig=-2*A*w*w+2*A*x*x-sum((w-4*np.pi**2/(w*(T+D)**2))*np.sin(w*T/2)*np.sin(w*D/2)/(6*np.pi*np.sin(w*(T+D)/2)) for D in Ds)
            assert abs(K-trig)<2e-12
            spec,bounds=zip(*(spectral(w,T,D,16384) for D in Ds))
            comparator=-2*A*w*w+K0+sum(spec)
            error=abs(K-comparator);bound=sum(bounds)
            assert error<=bound*1.001+2e-12
            z=w*w
            imaginary_identity=-z.imag*(K0/abs(z)**2+sum(positive_spectral_factor(w,T,D,16384) for D in Ds))
            finite_identity_error=abs((comparator/z).imag-imaginary_identity)
            assert finite_identity_error<2e-11
            assert abs((K/z).imag-imaginary_identity)<=bound/abs(z)+2e-11
            assert (K/(w*w)).imag*(w*w).imag<0
            rows.append(dict(returns=Ds,omega=[w.real,w.imag],kernel=[K.real,K.imag],
                flux_trigonometric_error=abs(K-trig),spectral_error=error,spectral_tail_bound=bound,
                independent_pole_identity_error=finite_identity_error,
                pole_exclusion_sign=(K/(w*w)).imag*(w*w).imag))
        for gamma in [.2*x,x,4*x]:
            K=2*A*gamma*gamma+2*A*x*x+sum(matter(1j*gamma,T,D) for D in Ds)
            assert K.real>0 and abs(K.imag)<2e-12
            axis.append(dict(returns=Ds,gamma=gamma,kernel=K.real))
    return equilibria,rows,axis


def cutoff_cases():
    rows=[];T=math.pi
    for eps in [.04,.02,.01]:
        Te=2*math.acos(eps);b=math.sqrt(1-eps*eps)
        for D in [0,.2,2,100,1e7]:
            for w in [1+.03j,1+.1j,2+.2j]:
                L=T+D;Le=Te+D
                B=-np.expm1(1j*w*T)
                Be=-np.expm1(1j*w*Te)-1j*w*eps/b*(1+np.exp(1j*w*Te))
                d0=-math.expm1(-w.imag*T/2)
                adiff=abs((2*math.pi/Le)**2-(2*math.pi/L)**2)
                Aden=abs(w*w-(2*math.pi/L)**2)
                den_diff=abs(np.expm1(1j*w*Le)-np.expm1(1j*w*L))
                tri=4/abs(w)/(24*math.pi)*(adiff*abs(Be)/d0+Aden*(abs(Be-B)/d0+abs(B)*den_diff/d0**2))
                err=abs(matter(w,T,D,eps)-matter(w,T,D))
                assert err<=tri+3e-12
                assert abs(Be-B)<=6*eps*abs(w)
                rows.append(dict(epsilon=eps,return_length=D,omega=[w.real,w.imag],error=err,explicit_triangle_bound=tri))
    # Direct ray inversion independently checks the endpoint displacement.
    endpoints=[]
    for eps in [.03,.1]:
        Te=2*math.acos(eps)
        for t in [3.2,5.1,8.7]:
            target=pulse(t,0)-pulse(t-Te,0)+eps/math.sqrt(1-eps*eps)*(pulse(t,1)+pulse(t-Te,1))
            vals=[]
            for amp in [.0004,.0002,.0001]:
                # s=acos(eps*tau') has the exact first variation of the
                # induced curve; its velocity correction starts quadratically.
                ss=lambda u:math.acos(eps*clock(u,amp,1))
                p=brentq(lambda u:clock(t,amp)-clock(u,amp)-ss(t)-ss(u),t-5,t,xtol=7e-14)
                response=(p-(t-Te))/amp
                vals.append(dict(amplitude=amp,response=response,error=abs(response-target)))
            rich=2*vals[-1]['response']-vals[-2]['response']
            assert abs(rich-target)<2e-5
            endpoints.append(dict(epsilon=eps,t=t,linear_displacement=target,refinements=vals,richardson_error=abs(rich-target)))
    return rows,endpoints


def main():
    start=time.monotonic()
    eq,res,axis=response_cases();cut,ends=cutoff_cases()
    data=dict(status='passed_transparent_clock_flux_checks',date='2026-10-11',exact_checks=exact_checks(),
        point_split_cases=point_split(),causal_state_cases=causal_cases(),equilibrium_cases=eq,
        numerical_pulse_regularities='Compact C7 pulse; enough for all displayed derivative checks. The patch Hadamard argument assumes a smooth history and is analytic, not certified by these samples.',
        response_cases=res,imaginary_axis_cases=axis,cutoff_cases=cut,endpoint_cases=ends,
        leading_free_transparent_response_derived=True,causal_free_state_for_prescribed_clock_constructed=True,
        full_four_dimensional_influence_action_derived=False,interacting_stress_bound_closed=False,
        finite_channel_decay_window_proved=False,physical_damping_proved=False,novelty_established=False,
        lean_invocations=0,note_sha256=sha(NOTE),script_sha256=sha(Path(__file__)),
        dependency_hashes={f:sha(ROOT/f) for f in ['notes/collective-clock-response.md','notes/proper-clock-anomaly.md',
            'notes/collar-conservative-response.md','notes/weak-coupling-window.md']},
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sy.__version__,mpmath=mp.__version__))
    data['elapsed_seconds']=time.monotonic()-start
    OUT.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:data[k] for k in ['status','elapsed_seconds','physical_damping_proved']}))
    print(json.dumps({'maximum_flux_trigonometric_error':max(v['flux_trigonometric_error'] for v in res),
        'maximum_endpoint_error':max(v['richardson_error'] for v in ends),
        'causal_linearization_ratios':[v['error_ratios'] for v in data['causal_state_cases']]}))


if __name__=='__main__':main()
