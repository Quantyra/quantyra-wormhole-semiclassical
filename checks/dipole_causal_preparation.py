"""Sharp initial-time obstruction and smooth diagnostic polar response bounds."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import mpmath as mp
import numpy as np
import sympy as s
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'checks'))
from dipole_vacuum_response import rho, norm


def geometry(C):
    x=1/C;r=mp.sqrt(1+12*x+12*x*x)
    sn=mp.sqrt(6*x*(1+x)/(r*(r+1)));cs=mp.sqrt((r+1)/(2*r))
    lo=24*x*x/(1+6*x+r);hi=1+6*x+r
    d=(1+mp.sqrt(9+16*x))/2;h=(1+mp.sqrt(1+4*lo))/2
    c=1/mp.sqrt(3*(C+1))
    return x,d,h,lo,hi,sn,cs,c


def amplitudes(C,lam):
    x,d,h,lo,hi,sn,cs,c=geometry(C)
    den=lam*lam-(2+12*x)*lam+24*x*x
    fx=mp.mpf(9)/4+(3*lam*(1-x)+36*x+72*x*x)/den
    fa=-mp.mpf(3)/4+(9*lam*x-12*x-60*x*x)/den
    light=-sn*c*fx+cs*fa  # symmetrized master divided by G E0
    chi_light=sn*sn*fx-sn*cs*fa/c
    chi_heavy=cs*cs*fx+sn*cs*fa/c
    assert abs(chi_light+chi_heavy-fx)<mp.mpf('1e-50')
    assert max(abs(chi_light),abs(chi_heavy))<=mp.mpf(24)/5
    return light,chi_light,chi_heavy


def sharp_checks():
    rows=[]
    for C in map(mp.mpf,['1','4','100000000']):
        x,d,h,lo,hi,sn,cs,c=geometry(C)
        v=-sn*c*mp.mpf(9)/4-cs*mp.mpf(3)/4
        assert v<-3/(4*mp.sqrt(2))
        limit=16*mp.sqrt(2*mp.pi)/15*norm(h,0)**2*v*v
        terms=[]
        for j in [100,10000,1000000]:
            D=2*d+2*j;lam=D*(D-1)
            J=norm(h,0)*norm(D,0)*mp.sqrt(mp.pi)*mp.gamma((h+D+1)/2)/mp.gamma((h+D+2)/2)
            f=amplitudes(C,lam)[0]
            velocity=(16*mp.pi/15)*rho(d,j)*f*f*(D/2)*J*J
            B=mp.exp(-1j*D*mp.pi/(2*h))+1j*D/h
            retarded=(16*mp.pi/15)*rho(d,j)*f*f*J*J/(2*D)*abs(B)**2
            terms.append({'j':j,'velocity_scaled_to_limit':str(velocity*mp.sqrt(j)/limit),
                          'sharp_response_scaled_to_limit':str(retarded*mp.sqrt(j)/(limit/h**2))})
        assert abs(mp.mpf(terms[-1]['velocity_scaled_to_limit'])-1)<mp.mpf('1e-4')
        assert abs(mp.mpf(terms[-1]['sharp_response_scaled_to_limit'])-1)<mp.mpf('1e-4')
        rows.append({'C':str(C),'light_weight':str(h),'local_light_coefficient':str(v),
                     'positive_velocity_asymptotic_constant':str(limit),'primary_asymptotics':terms})
    return rows


def switching_checks():
    t,O,w=s.symbols('t O w',positive=True,real=True)
    g=s.Function('g')(t);u=s.exp(-s.I*O*t)
    assert s.simplify(s.diff(g*u,t,2)+w*w*g*u-g*(s.diff(u,t,2)+w*w*u)
                      -(s.diff(g,t,2)*u+2*s.diff(g,t)*s.diff(u,t)))==0
    z=s.symbols('z',real=True);S=10*z**3-15*z**4+6*z**5
    v=s.diff(S,z)
    assert s.factor(v)==30*z*z*(z-1)**2
    assert s.integrate(v,(z,0,1))==1
    assert 2*v.subs(z,s.Rational(1,2))==s.Rational(15,4)
    assert all(s.diff(S,z,k).subs(z,0)==s.diff(S,z,k).subs(z,1)==0 for k in [1,2])
    rows=[]
    for T,omega,Omega,stop in [(.25,1.1,4.7,.95),(.25,4.,4.,.95),(.25,9.,4.7,.95),
                              (1.,2.3,17.,1.7),(2.,1.,31.,2.7),(.25,1.1,4.7,.1),
                              (.25,4.,4.,.2),(1.,2.3,17.,.7)]:
        # Independent ODE integration includes the complete multiplied original source.
        def switch(time):
            x=min(max(time/T,0.),1.)
            return 10*x**3-15*x**4+6*x**5
        def rhs(time,y):return np.array([y[1],-omega**2*y[0]+(omega**2-Omega**2)*switch(time)*np.exp(-1j*Omega*time)])
        sol=solve_ivp(rhs,(0,stop),np.array([0j,0j]),method='DOP853',rtol=2e-11,atol=2e-12,max_step=min(T/35,.03))
        assert sol.success
        tt=mp.mpf(str(T));ww=mp.mpf(str(omega));oo=mp.mpf(str(Omega));ts=mp.mpf(str(stop))
        upper=min(ts/tt,mp.mpf(1));gt=10*upper**3-15*upper**4+6*upper**5
        def F(nu):return mp.quad(lambda x:30*x*x*(1-x)**2*mp.exp(-1j*nu*tt*x),[0,upper])
        expected=gt*mp.exp(-1j*oo*ts)-(oo+ww)/(2*ww)*mp.exp(-1j*ww*ts)*F(oo-ww)+(oo-ww)/(2*ww)*mp.exp(1j*ww*ts)*F(oo+ww)
        err=abs(complex(expected)-sol.y[0,-1])
        assert err<2e-8
        K=1+2*max(1,15/(4*T))
        assert abs(expected)<=K
        if Omega==omega:assert abs(expected)<mp.mpf('1e-50')
        rows.append({'T':T,'observation_time':stop,'master_frequency':omega,'pair_frequency':Omega,'ODE_error':err,'amplitude':str(abs(expected)),'bound':K})
    return {'product_commutator':'exact','quintic_total_variation':'15/(4T)','independent_ODE_cases':rows}


def constant_checks():
    x=s.symbols('x',nonnegative=True)
    assert s.expand((1+12*x+12*x*x)-(1+x)**2)==10*x+11*x*x
    # Positive light/heavy masses are bounded by 2 and 12, respectively.
    assert s.expand(4-2*(2+12*x)+24*x*x)==24*x*x-24*x
    z,mu=s.symbols('z mu',real=True)
    f=8/(3*s.pi)*s.cos(z)**4
    Hf=s.trigsimp(-s.diff(f,z,2)+mu/s.cos(z)**2*f)
    target=8/(3*s.pi)*(16*s.cos(z)**4+(mu-12)*s.cos(z)**2)
    assert s.trigsimp(Hf-target)==0
    I4=s.Rational(3,8)*s.pi;I6=s.Rational(5,16)*s.pi;I8=s.Rational(35,128)*s.pi
    Hnorm=64/(9*s.pi**2)*(256*I8+32*(mu-12)*I6+(mu-12)**2*I4)
    assert s.simplify(Hnorm-64/(9*s.pi)*(s.Rational(3,8)*mu*mu+mu+4))==0
    pi_lo=s.Rational(157,50);pi_hi=s.Rational(22,7)
    A0_hi=8*pi_hi/15*(pi_hi**2/6-1)
    assert A0_hi<s.Rational(11,10)
    assert s.Rational(11,10)*s.Rational(5,4)<s.Rational(6,5)**2
    assert 5/pi_lo<s.Rational(13,10)**2
    assert 160/(3*pi_lo)<s.Rational(21,5)**2
    assert 4480/(9*pi_lo)<s.Rational(63,5)**2
    assert s.Rational(21,5)+s.Rational(63,5)<17
    coefficient=s.Rational(13,10)*s.Rational(24,5)*s.Rational(6,5)*31*17
    assert coefficient<4000
    candidates=[]
    for k0 in [2000,4000,10000]:
        k=mp.mpf(k0);C=mp.sqrt(1+12*mp.pi*k**6);b2=k**10*(C-1)/6
        candidates.append({'k':k0,'all_mode_diagnostic_sphere_RMS_upper':str(4000/b2),'physical_preparation_bound':False})
    return {'all_mode_certificate':'RMS spatially averaged fractional sphere response < 4000 G/b^2, T=1/4, every t>=0; all retained radial modes',
            'spatial_test':'8 cos(sigma)^4/(3 pi)','intermediate_rational_constant':str(coefficient),
            'candidate':candidates}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    mp.mp.dps=65
    out={'status':'PASS','scope':'Sharp-zero-data obstruction and all-mode smooth diagnostic source envelope on the fixed global H2 vacuum product.',
         'sharp_initial_data':sharp_checks(),'switching':switching_checks(),'all_mode_bound':constant_checks(),
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'dependencies':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
             ['checks/dipole_vacuum_response.py','notes/dipole-vacuum-response.md','notes/dipole-quadrupole-source.md','notes/dipole-causal-preparation.md']},
         'limitations':['Envelope multiplying master sources is not a conserved four-dimensional preparation.',
                        'No construction of the dressed constrained quantum state or its intrinsic metric noise.',
                        'No H3/H4, fermion, contact, massive-mode or physical-boundary error bound.',
                        'No renormalized energy/stress, capture, creation, stability or transport conclusion.']}
    args.output.write_bytes((json.dumps(out,indent=2)+'\n').encode())
    print(json.dumps({'status':'PASS','ODE_cases':len(out['switching']['independent_ODE_cases']),
                      'candidate':out['all_mode_bound']['candidate'][0]},indent=2),flush=True)
