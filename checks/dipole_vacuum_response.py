"""Stationary vacuum-to-pair particular response, not the driven physical state."""
from pathlib import Path
import argparse
import hashlib
import json
import math

import mpmath as mp
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def primary(delta,N):
    out=[mp.mpf(1)]
    for p in range(N):
        out.append(-out[-1]*mp.sqrt((N-p)*(N-p-1+2*delta)/((p+1)*(p+2*delta))))
    norm=mp.sqrt(mp.fsum(v*v for v in out))
    return [v/norm for v in out]


def central_weight(delta,n):
    return mp.gamma(n+mp.mpf('.5'))*mp.gamma(n+delta)/(2*mp.pi*mp.gamma(n+1)*mp.gamma(n+delta+mp.mpf('.5')))


def rho(delta,j):
    return (4*delta+4*j-1)/(4*mp.pi)*mp.gamma(j+mp.mpf('.5'))*mp.gamma(j+delta)**2*mp.gamma(j+2*delta-mp.mpf('.5'))/(mp.gamma(j+1)*mp.gamma(j+delta+mp.mpf('.5'))**2*mp.gamma(j+2*delta))


def representation_checks():
    mp.mp.dps=80;weights=[];grams=[]
    for delta in map(mp.mpf,['2','2.3','3']):
        for j in range(9):
            N=2*j;c=primary(delta,N)
            assert abs(mp.fsum(x*x for x in c)-1)<mp.mpf('1e-70')
            for p in range(N):
                residual=c[p+1]*mp.sqrt((p+1)*(p+2*delta))+c[p]*mp.sqrt((N-p)*(N-p-1+2*delta))
                assert abs(residual)<mp.mpf('1e-65')
            B=mp.fsum(c[p]*(-1)**j*mp.sqrt(central_weight(delta,p//2)*central_weight(delta,(N-p)//2)) for p in range(0,N+1,2))
            actual=B*B/central_weight(2*delta+N,0)
            err=abs(actual/rho(delta,j)-1)
            assert err<mp.mpf('1e-65')
            weights.append({'Delta':str(delta),'j':j,'relative_weight_error':str(err)})
        dd=float(delta)
        for M in range(19):
            basis=[]
            for N in range(0,M+1,2):
                v=np.array([float(a) for a in primary(delta,N)])
                for level in range(N,M):
                    raised=np.zeros(level+2)
                    for p in range(level+1):
                        q=level-p
                        raised[p+1]+=math.sqrt((p+1)*(p+2*dd))*v[p]
                        raised[p]+=math.sqrt((q+1)*(q+2*dd))*v[p]
                    step=level-N
                    v=raised/math.sqrt((step+1)*(step+4*dd+2*N))
                assert np.max(abs(v-v[::-1]))<1e-12
                basis.append(v)
            mat=np.array(basis)
            error=float(np.max(abs(mat@mat.T-np.eye(len(basis)))))
            assert len(basis)==M//2+1 and error<1e-11
            grams.append(error)
    return {'independent_primary_weights':weights,'symmetric_energy_levels':len(grams),'max_descendant_Gram_error':max(grams)}


def Wcenter(delta,time):
    return central_weight(delta,0)*mp.exp(-delta*time)*mp.hyp2f1(mp.mpf('.5'),delta,delta+mp.mpf('.5'),mp.exp(-2*time))


def KL_checks():
    rows=[]
    for delta in map(mp.mpf,['2','2.3','3']):
        for time in map(mp.mpf,['.45','.8','1.2']):
            lhs=Wcenter(delta,time)**2
            rhs=mp.fsum(rho(delta,j)*Wcenter(2*delta+2*j,time) for j in range(100))
            error=abs(rhs/lhs-1)
            assert error<mp.mpf('1e-35')
            rows.append({'Delta':str(delta),'Euclidean_time':str(time),'terms':100,'relative_error':str(error)})
    return rows


def algebra_checks():
    u,e,lam=s.symbols('u e lam',nonnegative=True)
    M=s.Matrix([[2+6*u,6*(1+u)],[2*u,6*u]])
    v=s.Matrix([s.Rational(9,4),-s.Rational(3,4)])
    R=s.Matrix([3*(1-u),9*u])
    D=lam**2-(2+12*u)*lam+24*u*u
    cx=s.Rational(9,4)+(3*lam*(1-u)+36*u+72*u*u)/D
    ca=-s.Rational(3,4)+(9*lam*u-12*u-60*u*u)/D
    assert (v+(lam*s.eye(2)-M).inv()*R-s.Matrix([cx,ca])).applyfunc(s.factor)==s.zeros(2,1)
    assert s.factor((1+6*u)**2-(1+12*u+12*u*u))==24*u*u
    # lambda >= 12+16u and 0<=u<=1 suffice for all following signs.
    px=lam**2-(12+2*u)*lam-120*u-216*u*u
    assert s.expand(px.subs(lam,12+16*u+e)-(e*e+(12+30*u)*e+48*u+8*u*u))==0
    pa=lam**2-(2+24*u)*lam+16*u+104*u*u
    assert s.expand(pa.subs(lam,12+16*u+e)-(e*e+(22+8*u)*e+120+56*u+24*u*(1-u)))==0
    assert s.factor((9*lam*u-12*u-60*u*u).subs(lam,12+16*u+e)-3*u*(32+28*u+3*e))==0
    return {'particular_response':'exact rational inverse','Casimir_gap_lower':'10+4/C',
            'metric_coefficient_interval':['9/4','51/20'],'Maxwell_coefficient_interval':['-3/4','0'],
            'method':'exact polynomial nonnegativity certificates on 0<=1/C<=1'}


def norm(delta,n):
    return mp.sqrt(mp.factorial(n)*(n+delta)*mp.gamma(delta)**2/(mp.pi*2**(1-2*delta)*mp.gamma(n+2*delta)))


def overlap_checks():
    rows=[]
    for C in map(mp.mpf,['1','4','61','1000000','49119681981.431447']):
        dd=(1+mp.sqrt(9+16/C))/2
        muplus=(C+6+mp.sqrt(C*C+12*C+12))/C
        dh=(1+mp.sqrt(1+4*muplus))/2
        power=dh+2*dd-2
        I0=mp.sqrt(mp.pi)*mp.gamma((power+1)/2)/mp.gamma((power+2)/2)
        factor=norm(dh,2)*norm(dd,0)**2*dh*I0/(dh+2*dd)
        analytic=factor*(dh+2-2*dd)
        def integrand(x):
            return norm(dh,2)*mp.cos(x)**dh*(2*dh*(dh+1)*mp.sin(x)**2-dh)*norm(dd,0)**2*mp.cos(x)**(2*dd-2)
        numerical=mp.quad(integrand,[-mp.pi/2,0,mp.pi/2])
        assert abs(numerical-analytic)<mp.mpf('1e-65')
        detuning=(2*dd)**2-(dh+2)**2
        continuation=-factor/(2*dd+dh+2)
        if C!=1:assert abs(analytic/detuning-continuation)<mp.mpf('1e-60')
        else:assert analytic==0 and detuning==0
        rows.append({'C':str(C),'overlap':str(analytic),'frequency_squared_detuning':str(detuning),
                     'ratio_continuation':str(continuation),'quadrature_error':str(abs(numerical-analytic)),
                     'exact_resonance_forces_no_mode':C==1})
    return rows


def all_branch_overlap_checks():
    """Independent spatial integrals test the general Green identity, not just n=2."""
    rows=[]
    for C in map(mp.mpf,['1','4','1000000']):
        delta=(1+mp.sqrt(9+16/C))/2
        root=mp.sqrt(C*C+12*C+12)
        for branch,mu in [('light',24/(C*(C+6+root))),('heavy',(C+6+root)/C)]:
            h=(1+mp.sqrt(1+4*mu))/2
            for j,k,n in [(0,0,0),(0,0,2),(1,1,1),(2,2,4)]:
                D=2*delta+2*j;lam=D*(D-1)
                # Finite Gegenbauer recurrence avoids hypergeometric evaluation at x=0.
                def gegen(alpha,degree,x):
                    if degree==0:return mp.mpf(1)
                    a,b=mp.mpf(1),2*alpha*x
                    for ell in range(1,degree):
                        a,b=b,(2*(ell+alpha)*x*b-(ell+2*alpha-1)*a)/(ell+1)
                    return b
                def product(x,power):
                    return norm(h,n)*norm(D,k)*mp.cos(x)**(h+D-power)*gegen(h,n,mp.sin(x))*gegen(D,k,mp.sin(x))
                J=mp.quad(lambda x:product(x,0),[-mp.pi/2,0,mp.pi/2])
                I=mp.quad(lambda x:product(x,2),[-mp.pi/2,0,mp.pi/2])
                detuning=(D+k)**2-(h+n)**2
                residual=detuning*J-(lam-mu)*I
                assert abs(residual)<mp.mpf('1e-60')
                rows.append({'C':str(C),'branch':branch,'j':j,'descendant':k,'master_mode':n,
                             'Green_identity_error':str(abs(residual))})
    t,Omega,omega=s.symbols('t Omega omega',positive=True,real=True)
    B=s.exp(-s.I*Omega*t)-s.cos(omega*t)+s.I*Omega/omega*s.sin(omega*t)
    assert s.simplify(s.diff(B,t,2)+omega**2*B-(omega**2-Omega**2)*s.exp(-s.I*Omega*t))==0
    assert B.subs(t,0)==0 and s.diff(B,t).subs(t,0)==0
    assert s.simplify(s.expand_complex(B.subs(Omega,omega)))==0
    return {'spatial_identity':rows,'retarded_mode_identity':'exact ODE and zero initial data; exact resonance gives zero forcing/response',
            'continuum_initial_data_bound':False}


def variance_checks():
    # Rational upper bounds use e^-x <= 1/sum_(k=0)^32 x^k/k!.
    def exp_upper(x):return 1/sum(x**k/s.factorial(k) for k in range(33))
    partial=sum((m+1)*exp_upper(s.Rational((m+2)**2,4)) for m in range(5))
    ratio=s.Rational(7,6)*exp_upper(s.Rational(15,4))
    upper=partial+6*exp_upper(s.Rational(49,4))/(1-ratio)
    # Machin's identity and alternating atan series certify pi>157/50.
    atan_lower=sum((-1)**k*s.Rational(1,5)**(2*k+1)/(2*k+1) for k in range(4))
    pi_lower=16*atan_lower-s.Rational(4,239)
    assert pi_lower>s.Rational(157,50)
    var_bound=64*upper/(135*s.Rational(157,50))
    assert var_bound<s.Rational(1,10)
    assert s.Rational(51,20)**2/(2*s.Rational(157,50))<s.Rational(1018,1000)**2
    rows=[]
    for C in map(mp.mpf,['1','4','100000000']):
        delta=(1+mp.sqrt(9+16/C))/2
        dvar=mp.fsum(central_weight(delta,p)*central_weight(delta,q)*mp.exp(-mp.mpf('.25')**2*(2*p+2*q+2*delta)**2) for p in range(22) for q in range(22))*16*mp.pi/15
        kvar=mp.fsum(rho(delta,j)*mp.fsum(central_weight(2*delta+2*j,n)*mp.exp(-mp.mpf('.25')**2*(2*n+2*delta+2*j)**2) for n in range(22)) for j in range(22))*16*mp.pi/15
        assert abs(dvar-kvar)<mp.mpf('1e-55') and dvar<mp.mpf('.1')
        rows.append({'C':str(C),'Q_variance_truncated':str(dvar),'KL_double_mode_difference':str(abs(dvar-kvar))})
    candidates=[]
    for k0 in [2000,4000,10000]:
        k=mp.mpf(k0);C=mp.sqrt(1+12*mp.pi*k**6);b2=k**10*(C-1)/6
        rms=mp.mpf('1.018')/b2
        candidates.append({'k':k0,'b_squared_Planck_units':str(b2),'smeared_fractional_sphere_metric_RMS_upper':str(rms),
                           'physical_capture_bound':False})
    return {'uniform_Q_variance_upper':'1/10','rational_bound_value':str(s.N(var_bound,22)),
            'Gaussian_time_width':'1/4','point':'global throat center',
            'numerical_covariance_checks':rows,'candidate':candidates}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    out={'status':'PASS','scope':'stationary vacuum-to-two-dipole particular solution on global product with Friedrichs boundaries',
         'representations':representation_checks(),'KL_normalization':KL_checks(),
         'response_bounds':algebra_checks(),'near_resonant_overlap':overlap_checks(),
         'all_branch_overlap':all_branch_overlap_checks(),'smeared_variance':variance_checks(),
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'dependencies':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
                         ['notes/dipole-quadrupole-source.md','notes/quadrupole-product-response.md','notes/dipole-vacuum-response.md']},
         'versions':{'mpmath':mp.__version__,'numpy':np.__version__,'sympy':s.__version__},
         'limitations':['A regular particular solution does not select causal initial data or a dressed interacting state.',
                        'No bound on the intrinsic homogeneous quadrupole field, H3/H4, fermions or response contacts.',
                        'Global vacuum and Gaussian measurement are not physical horizon or finite switching conditions.',
                        'No full physical capture, creation, preparation, stability or transport result.']}
    args.output.write_bytes((json.dumps(out,indent=2)+'\n').encode())
    print(json.dumps({'status':'PASS','max_Gram_error':out['representations']['max_descendant_Gram_error'],
                      'Q_variance_upper':out['smeared_variance']['uniform_Q_variance_upper'],
                      'candidate':out['smeared_variance']['candidate'][0]},indent=2),flush=True)
