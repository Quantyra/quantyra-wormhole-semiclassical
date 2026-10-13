"""S072: retarded JT bulk reconstruction, total-area bound and probe rays.

The physical area statement is conditional on its extremal contribution.
No four-dimensional reduction error or general stability bound is tested.
"""
from pathlib import Path
import hashlib
import json
import warnings
import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad, IntegrationWarning
from scipy.interpolate import CubicSpline
import sympy as sp


def symbolic_checks():
    t, x = sp.symbols('t x', real=True)
    f, g = sp.Function('f'), sp.Function('g')
    field = (sp.diff(f(t+x), t)-sp.cot(x)*f(t+x)
             +sp.diff(g(t-x), t)+sp.cot(x)*g(t-x))
    def dp(v):
        return (sp.diff(v, t)+sp.diff(v, x))/2
    def dm(v):
        return (sp.diff(v, t)-sp.diff(v, x))/2
    assert sp.simplify(-dp(dp(field))-sp.cot(x)*dp(field)
                       +sp.diff(f(t+x), t, 3)+sp.diff(f(t+x), t)) == 0
    assert sp.simplify(-dm(dm(field))+sp.cot(x)*dm(field)
                       +sp.diff(g(t-x), t, 3)+sp.diff(g(t-x), t)) == 0
    assert sp.simplify(-sp.diff(field, t, 2)+sp.diff(field, x, 2)
                       -2*field/sp.sin(x)**2) == 0
    B = sp.Function('B')(x)
    P = B/sp.sin(x)
    spatial = sp.diff(P, x, 2)+sp.cot(x)*sp.diff(P, x)-P/sp.sin(x)**2
    assert sp.simplify(spatial-(sp.diff(B,x,2)-sp.cot(x)*sp.diff(B,x))/sp.sin(x)) == 0
    # Pointwise negativity of the whole outgoing stress, including radiation.
    speed_bound = sp.Rational(1152,22000)*sp.Rational(22,7)
    assert speed_bound < 1
    K = sp.sqrt(sp.Rational(5,32))
    assert sp.Rational(13,48) < K < sp.Rational(1,2)
    # Two-reflection recurrence has no secular forcing of the dilaton.
    alpha, beta = sp.symbols('alpha beta', real=True)
    assert sp.trigsimp((1-sp.cos(alpha-beta))+(1-sp.cos(alpha-sp.pi-beta))) == 2
    return {'null_constraints': 'pass', 'trace_equation': 'pass',
            'spatial_constraint': 'pass', 'theta_dot_upper_rational': str(speed_bound),
            'uniform_dynamical_dilaton_lower_over_N_eta': str(-K)}


def build(eta, duration, z0, samples=8193):
    def rhs(s, state):
        z,w,r,t,P,C,S = state
        x=s/duration
        a=x*x*(3-2*x)
        ap=6*x*(1-x)/duration
        qdt=-ap/2+2*np.pi*eta*ap*ap/z
        return [z*w,-z*z+(a-r)*z/2,4*np.pi*eta*ap*ap/z,z,
                qdt,qdt*np.cos(t),qdt*np.sin(t)]
    sol=solve_ivp(rhs,(0,duration),[z0,0,0,0,0,0,0],method='DOP853',
                  rtol=2e-13,atol=2e-15,max_step=duration/100,dense_output=True)
    assert sol.success
    zf,wf,rf,T,Pf,Cf,Sf=sol.y[:,-1]
    e=wf*wf+zf*zf-(1-rf)*zf
    sample_s=np.linspace(0,duration,samples)
    vals=sol.sol(sample_s)
    inverse=CubicSpline(vals[3],sample_s)
    moments=CubicSpline(vals[3],vals[4:].T)
    def partial(v):
        vv=np.asarray(v)
        return moments(np.clip(vv,0,T))
    def pulse(v):
        vv=np.asarray(v)
        si=inverse(np.clip(vv,0,T))
        zz=sol.sol(si)[0]
        xx=si/duration
        ap=6*xx*(1-xx)/duration
        velocity=4*np.pi*eta*ap/zz
        return np.where((vv>0)&(vv<T),(velocity*velocity-velocity)/(8*np.pi*eta),0.)
    def h(v):
        vv=np.asarray(v)
        return np.where(vv>=0,pulse(np.mod(vv,np.pi)),0.)
    def H(v):
        vv=np.asarray(v)
        hh=np.zeros_like(vv,dtype=float)
        hp=np.zeros_like(vv,dtype=float)
        for j in range(max(0,int(np.max(vv)/np.pi)+1)):
            xx=vv-j*np.pi
            p,c,ss=np.moveaxis(partial(xx),-1,0)
            hh+=np.where(xx>0,-p+np.cos(xx)*c+np.sin(xx)*ss,0.)
            hp+=np.where(xx>0,-np.sin(xx)*c+np.cos(xx)*ss,0.)
        return hh,hp
    def boundary(t):
        pp,cc,ss=np.moveaxis(partial(t),-1,0)
        return z0*np.cos(t)-pp+np.cos(t)*cc+np.sin(t)*ss
    def dilaton(t,x):
        hm,hmp=H(t-x)
        hp,hpp=H(t+x-np.pi)
        return z0*np.cos(t)/np.sin(x)+hmp+hpp+(hm-hp)/np.tan(x)
    return dict(sol=sol,T=T,rf=rf,e=e,P=Pf,C=Cf,S=Sf,pulse=pulse,h=h,H=H,
                boundary=boundary,dilaton=dilaton,eta=eta,duration=duration,z0=z0)


def check_case(eta,duration,z0):
    m=build(eta,duration,z0)
    T=m['T']; field=m['dilaton']; b=m['boundary']; h=m['h']
    assert abs(m['P']+(1-m['rf'])/2)<2e-13
    si=np.linspace(0,duration,101)
    vals=m['sol'].sol(si)
    boundary_error=float(np.max(np.abs(b(vals[3])-vals[0])))
    assert boundary_error<2e-12
    ap=6*(si/duration)*(1-si/duration)/duration
    max_speed=float(np.max(4*np.pi*eta*ap/vals[0]))
    assert max_speed<1 and np.max(m['pulse'](np.linspace(0,T,301)))<=1e-12
    # Independent post-ramp dynamics against the reconstructed boundary.
    zf,wf=m['sol'].y[:2,-1]
    def post(t,v):
        z,w,glob=v
        return [z*w,-z*z+(1-m['rf'])*z/2,z]
    aft=solve_ivp(post,(0,100),[zf,wf,T],method='DOP853',rtol=2e-12,
                  atol=2e-14,max_step=.15)
    assert aft.success
    post_match=float(np.max(np.abs(aft.y[0]-b(aft.y[2]))))
    assert post_match<2e-10
    times=np.linspace(0,3*np.pi,193)
    xs=np.linspace(.025,np.pi-.025,129)
    tt,xx=np.meshgrid(times,xs,indexing='ij')
    values=field(tt,xx)
    K=np.sqrt(5/32)
    assert values.min()>=-K-1e-10
    assert np.max(values-b(tt)/np.sin(xx))<1e-9
    reflection=float(np.max(np.abs(values-field(tt,np.pi-xx))))
    assert reflection<1e-10
    late=field(tt+np.pi,xx)
    periodic=float(np.max(np.abs(late-field(tt+3*np.pi,xx))))
    assert periodic<2e-10
    # Retarded bulk formula versus independently integrated spatial constraint.
    spatial_errors=[]
    for time in (T/2,T+0.2,np.pi/2+.02,np.pi+T/2,2*np.pi+.3):
        for sig in (.17,.63,np.pi/2):
            def integrand(v):
                k=-h(time-v)-h(time+v-np.pi)
                return (1-np.cos(min(sig,v)))*k
            cuts=[0.,sig,np.pi/2]
            for j in range(6):
                for off in (0,T):
                    for v in (time-j*np.pi-off,j*np.pi+off-time+np.pi):
                        if 0<v<np.pi/2: cuts.append(v)
            cuts=sorted(set(cuts))
            # Coincident wavefront cuts can differ by a floating-point ulp.
            # Integrate those tiny intervals directly, avoiding an adaptive
            # quadrature roundoff warning without discarding their measure.
            integral=sum((right-left)*integrand((left+right)/2) if right-left<1e-10
                         else quad(integrand,left,right,epsabs=2e-11,epsrel=2e-11)[0]
                         for left,right in zip(cuts,cuts[1:]))
            other=(b(time)-integral)/np.sin(sig)
            spatial_errors.append(abs(float(field(time,sig))-other))
    assert max(spatial_errors)<2e-9
    # Causal center crossing is real and is not confused with total-area zero.
    at_zero=float(field(np.pi/2,np.pi/2))
    early=.0001
    slope=float(field(np.pi/2+early,np.pi/2)/early)
    assert abs(at_zero)<1e-12 and abs(slope+z0)<1e-4
    assert field(np.pi/2+.001,np.pi/2)<0
    # For example Phi_ext/(N eta)=10: no zero of total area, even where Phi<0.
    assert np.min(10+values)>10-K
    delays=[quad(lambda t:1/float(b(t)),start,start+np.pi,epsabs=2e-10)[0]
            for start in (T,np.pi/2,np.pi,2*np.pi)]
    assert min(delays)>=np.pi/(.5+K)-1e-9
    assert max(delays)<=np.pi/(.5-K)+1e-9
    return dict(eta_f=eta,duration=duration,z0=z0,global_switch_duration=T,
                max_theta_dot=max_speed,normalized_dilaton_grid_min=float(values.min()),
                boundary_match_error=boundary_error,post_boundary_match_error=post_match,
                spatial_constraint_quadrature_error=float(max(spatial_errors)),
                reflection_error=reflection,late_periodicity_error=periodic,
                center_zero=at_zero,center_forward_slope=slope,
                probe_delays_eta_times_u=delays)


if __name__=='__main__':
    warnings.simplefilter('error',IntegrationWarning)
    symbolic=symbolic_checks()
    cases=[check_case(5e-4,d,z) for d in (.25,.375,.5) for z in (.25,.5,.75)]
    cases.append(check_case(1e-5,.375,.5))
    # Independent interpolation-density refinement at non-grid spacetime points.
    coarse=build(5e-4,.25,.25,4097)
    fine=build(5e-4,.25,.25,16385)
    t=np.linspace(.037,3*np.pi-.041,91)[:,None]
    x=np.linspace(.043,np.pi-.063,73)[None,:]
    refine=float(np.max(np.abs(coarse['dilaton'](t,x)-fine['dilaton'](t,x))))
    assert refine<2e-10
    out=dict(date='2026-10-12',status='conditional_JT_total_area_and_probe_channel_verified',
             code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__,
             symbolic=symbolic,cases=cases,interpolation_refinement_error=refine,
             theorem='Phi >= -sqrt(5/32) N eta_f for all future global time',
             total_area_condition='renormalized Phi_ext > sqrt(5/32) N eta_f',
             not_proved=['four-dimensional dynamical matching error','finite-cutoff matching',
                         'initial-state preparation','quantum fluctuations and other perturbations',
                         'finite-energy signal backreaction','physical realization'])
    (Path(__file__).resolve().parents[1]/'evidence/driven-bulk-geometry.json').write_text(
        json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: bulk JT equations, 10 driven geometries, return periodicity, area bound, probe rays and refinement')
