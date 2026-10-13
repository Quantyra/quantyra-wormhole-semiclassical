"""S072: massive one-loop electric response on R^(1,1) x S^2.

The sphere spectrum is exact. The longitudinal space is flat, the state
is vacuum, and the source in the pulse test is prescribed and classical.
No AdS2 remainder, interacting quantum switch or physical capture is tested.
"""
from pathlib import Path
import hashlib
import json
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.special import roots_legendre, jv

ROOT = Path(__file__).resolve().parents[1]


def harmonic2(q):
    return mp.zeta(2)-mp.polygamma(1,q+1)


def delta_coefficient(q,C,Z):
    """Subtracted Euclidean F^2 coefficient, in units with a=1."""
    q,C,Z = map(mp.mpf,(q,C,Z))
    def integrand(x):
        u=x*(1-x);v=C*Z*u
        rm=2*v/(q+mp.sqrt(q*q-4*v))
        return u*mp.re(mp.digamma(1)+mp.digamma(q+1)
                       -mp.digamma(1+rm)-mp.digamma(q+1-rm))
    return 2*C/mp.pi*mp.quad(integrand,[0,mp.mpf('.5'),1])


def sum_checks():
    rows=[]
    gx,gw=roots_legendre(96);u=(gx+1)/2;weights=gw/2;v=u*(1-u)
    for q,C,Z in [(1,2,.3),(3,7,20),(17,31,200),(100,61,300)]:
        n=np.arange(1,2001,dtype=float);mass2=n*(n+q)/C;d=2*(q+2*n)
        partial=np.sum(d[:,None]*(-Z*v[None,:]**2)
                       /(mass2[:,None]*(mass2[:,None]+Z*v[None,:]))*weights)/np.pi
        exact=float(delta_coefficient(q,C,Z))
        tail4=2*C*C/q*sum(1/j**2 for j in range(2001,2001+q))
        tail_bound=Z*tail4/(30*np.pi)
        assert 0 <= partial-exact <= tail_bound+2e-11
        full4=2*mp.mpf(C)**2/q*harmonic2(q)
        partial4=np.sum(d/mass2**2)
        assert abs(float(full4)-partial4-tail4)<2e-11
        bound=float(Z*full4/(30*mp.pi))
        assert 0 < -exact <= bound
        # Common physical angular-mass cutoff, b=mu=1.
        target=(2*mp.euler-(mp.digamma(q+1)+mp.euler))/(3*mp.pi)
        errors=[]
        for N in [100,1000,10000]:
            val=(mp.harmonic(N)+mp.harmonic(N+q)-mp.harmonic(q)
                 -mp.log(N*(N+q)))/(3*mp.pi)
            errors.append(float(abs(val-target)))
        assert errors[2]<errors[1]<errors[0]
        rows.append(dict(q=q,C=C,Z=Z,subtracted_coefficient=exact,
                         partial_sum_error=float(partial-exact),tail_bound=tail_bound,
                         euclidean_absolute_bound=bound,static_subtraction_errors=errors))
    return rows


def cut_checks():
    rows=[]
    sigma3=np.diag([1.,-1.])
    for m in [.3,1.,4.]:
        for ratio in [1.001,1.2,3.,20.]:
            omega=2*m*ratio;p=np.sqrt(omega*omega/4-m*m)
            energies,vectors=np.linalg.eigh(np.array([[p,m],[m,-p]]))
            element=abs(vectors[:,1]@sigma3@vectors[:,0])**2
            # Two roots in the pair cut; Jacobian d(2 E_p)/dp.
            im_kappa=element*(omega/2)/(2*p)/(omega*omega)
            formula=2*m*m/(omega**4*np.sqrt(1-4*m*m/omega**2))
            assert abs(im_kappa/formula-1)<2e-12
            rows.append(dict(m=m,omega_over_threshold=ratio,
                             relative_cut_error=float(abs(im_kappa/formula-1))))
    # The primitive is independently integrated after s=4m^2+t^2.
    primitive_errors=[]
    for top in [mp.mpf('4.01'),mp.mpf('7'),mp.mpf('100')]:
        integral=mp.quad(lambda t:4/(4+t*t)**mp.mpf('1.5'),[0,mp.sqrt(top-4)])
        primitive_errors.append(float(abs(integral-mp.sqrt(1-4/top))))
    assert max(primitive_errors)<1e-45
    return dict(pair_cuts=rows,primitive_errors=primitive_errors)


def cumulative_checks():
    rows=[]
    for q in [1,3,100,10000]:
        largest=0.
        for top in np.geomspace(4*(q+1)*1.000001,4*(q+1)*1e4,50):
            N=int(np.floor(top/(2*(q+np.sqrt(q*q+top)))))
            # N = floor[(sqrt(q^2+top)-q)/2], evaluated without cancellation.
            n=np.arange(1,N+1,dtype=float)
            value=np.sum(2*(q+2*n)*np.sqrt(np.maximum(0,1-4*n*(n+q)/top)))
            bound=top/2*(1+1/q)
            assert value<=bound*(1+1e-12)
            largest=max(largest,value/bound)
        rows.append(dict(q=q,largest_cumulative_fraction=largest))
    return rows


def pulse_fourier(y):
    y=np.asarray(y)
    return 6*np.sqrt(np.pi)*(2/y)**3.5*jv(3.5,y)


def pulse_sum(q,tau,order,N=300):
    x,w=roots_legendre(order);t=(x+1)/2;w=w/2
    n=np.arange(1,N+1,dtype=float);mass=np.sqrt(n*(n+q))
    vals=tau*pulse_fourier(2*mass[:,None]*tau/np.sqrt(1-t[None,:]**2))
    work=np.sum(2*(q+2*n)*np.sum(vals*vals*w[None,:],axis=1))/(2*np.pi)
    Y=2*np.sqrt((N+1)*(N+1+q))*tau
    D=96*(1+6/Y+15/Y**2+15/Y**3)
    tail=D*D/(1890*np.pi*tau**6*N**6)
    return work,tail


def pulse_checks():
    rows=[]
    a0=32/35;a2=384/(25*np.sqrt(5))
    # Independent time integral checks the Fourier normalization.
    fourier_errors=[]
    for y in [.2,2.,10.,30.]:
        direct=mp.quad(lambda x:(1-x*x)**3*mp.cos(mp.mpf(y)*x),[-1,0,1])
        error=float(abs(mp.mpf(float(pulse_fourier(y)))-direct))
        assert error<1e-14
        fourier_errors.append(error)
    for q,tau in [(1,.1),(3,.2),(100,.05)]:
        w1,tail=pulse_sum(q,tau,512)
        w2,_=pulse_sum(q,tau,1024)
        refined,_=pulse_sum(q,tau,2048)
        bound=(1+1/q)*a0*a2/(2*np.pi)
        assert abs(refined-w2)<2e-8
        assert refined+tail<bound
        rows.append(dict(q=q,tau=tau,computed_work=refined,
                         refinement_error=abs(refined-w2),previous_refinement=abs(w2-w1),
                         omitted_level_tail_bound=tail,all_frequency_work_bound=bound))
    return dict(fourier_normalization_errors=fourier_errors,pulses=rows)


def pole_checks():
    # Use a visible coupling so a root shift can be independently resolved.
    rows=[]
    for q,C,g in [(1,2,1.),(3,7,.7),(17,31,.5)]:
        K=4*mp.pi*C/g**2;s0=q/(mp.pi*K);m12=mp.mpf(q+1)/C
        # Z=-s is below the first pair cut, so the same digamma sum is real.
        root=mp.findroot(lambda s:s*(K+delta_coefficient(q,C,-s))-q/mp.pi,(s0*.99,s0))
        bound=g**4*harmonic2(q)/(240*mp.pi**4*(1-s0/(4*m12)))
        shift=1-root/s0
        assert 0 < shift <= bound
        rows.append(dict(q=q,C=C,g=g,relative_pole_shift=float(shift),bound=float(bound)))
    return rows


def candidate(ki):
    k=mp.mpf(ki);q=k**4;g=k**-4
    X=(mp.sqrt(1+12*mp.pi*k**6)-1)/6;C=1+6*X
    Hq=mp.digamma(q+1)+mp.euler
    # Specified flat-product spectral subtraction scheme, mu=1/b.
    static_fraction=g*g*(2*mp.euler-Hq)/(12*mp.pi**2)
    gstar2=g*g/(1+static_fraction)
    K=4*mp.pi*C/gstar2
    om2=q*q/(C*(1+3*X));gam=4*mp.pi/k*om2
    estimates=[]
    for Z in [om2,4*gam*gam]:
        fraction=-delta_coefficient(q,C,Z)/K
        bound=gstar2*C*Z*harmonic2(q)/(60*mp.pi**2*q)
        assert 0 < fraction < bound
        estimates.append(dict(euclidean_invariant=float(Z),fraction=float(fraction),bound=float(bound)))
    threshold_ratio=gstar2*q/(16*mp.pi**2*(q+1))
    pole_bound=gstar2**2*harmonic2(q)/(240*mp.pi**4*(1-threshold_ratio))
    pulse_bound=3072*gstar2*(1+1/q)/(875*mp.sqrt(5)*mp.pi**2)
    return dict(k=ki,q=int(q),renormalized_static_fraction=float(static_fraction),
                euclidean_comparisons=estimates,light_pole_to_pair_threshold_ratio=float(threshold_ratio),
                relative_light_pole_shift_bound=float(pole_bound),
                prescribed_pulse_work_over_peak_field_energy_bound=float(pulse_bound))


def main():
    mp.mp.dps=60
    result=dict(date='2026-10-13',status='PASS',python=platform.python_version(),
                numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                subtracted_sums=sum_checks(),absorptive_cut=cut_checks(),
                cumulative_measure=cumulative_checks(),prescribed_pulses=pulse_checks(),
                light_poles=pole_checks(),candidate=[candidate(k) for k in [2000,4000,8000]],
                limitations=['Exact sphere spectrum but flat longitudinal R^(1,1), zero-temperature vacuum',
                             'Quadratic one-loop response; prescribed classical pulse, not the interacting compact source',
                             'Static coefficient uses a stated common subtraction; its AdS/RN matching remains open',
                             'No curved-throat, thermal, quantum-noise, full-stress or physical capture error bound'])
    (ROOT/'evidence/massive-bulk-polarization.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: common subtraction, tower sums, pair cut, cumulative measure, all-frequency pulse and light pole.')
    print(json.dumps(result['candidate'][0],indent=2))


if __name__=='__main__':
    main()
