"""Finite-channel linear clock decay: independent checks, not a 4D theorem."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import eigh
from scipy.optimize import brentq
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
PI = np.pi
KAPPA = 1 + np.arcsinh(np.sqrt(3)) / (2*np.sqrt(3))
GAMMA = 4*KAPPA*(8*np.log(2)-13*PI**2/24)/(3*PI**2)


def lengths(q, delta):
    return KAPPA*delta*q/(np.arange(q)+0.5)


def force(x, delta):
    a = KAPPA*delta
    # Integral of the static energy, differentiated analytically.
    return x/8-1/8+a/(6*PI)*(2*x*np.log1p(PI/(a*x))-PI*x/(a*x+PI))


def equilibrium(q, delta):
    ds = lengths(q, delta)
    return brentq(lambda x:x/8+1/24-PI**2*np.mean((PI+ds*x)**-2)/6,
                  0.05, 1.01, xtol=5e-15)


def circle(ds, w, tt):
    aa = np.exp(1j*w*tt)
    bb = np.exp(1j*w*ds)
    return -1j/(12*PI)*(w-4*PI**2/(w*(tt+ds)**2))*(1-aa)*(1-bb)/(1-aa*bb)


def kernel(q, delta, w, x=None):
    if x is None:
        x = equilibrium(q, delta)
    return -w*w/8+x*x/8+np.mean(circle(lengths(q, delta), w, PI/x))


def variation_bound(w, tt):
    r = -np.expm1(-w.imag*tt)
    a = abs(1-np.exp(1j*w*tt))
    return (8*PI**2*a/(abs(w)*tt**2*r)
            +a*a*(abs(w)**2+4*PI**2/tt**2)/(w.imag*r))/(12*PI)


def continuum_complex(delta, w, x, order=80):
    """Gauss quadrature on return-phase cells, with exponentially small tail.

    Refinement is numerical evidence. The independent mathematical bound is
    the variation estimate in the note, not this quadrature error estimate.
    """
    a = KAPPA*delta
    period = 2*PI/w.real
    cutoff = 24/w.imag
    nodes, weights = np.polynomial.legendre.leggauss(order)
    # Half-period cells put poles at endpoints; refinement resolves them.
    cuts = np.unique(np.r_[a, np.arange(period/4, cutoff, period/2), cutoff])
    cuts = cuts[cuts >= a]
    mid = (cuts[1:]+cuts[:-1])/2
    half = np.diff(cuts)/2
    dd = mid[:, None]+half[:, None]*nodes
    integral = np.sum(half*np.sum(weights*circle(dd, w, PI/x)/dd**2, axis=1))
    tt = PI/x
    const = -1j*w*(1-np.exp(1j*w*tt))/(12*PI)
    # Leading tail integrated exactly; omitted inverse-length term bounded.
    integral += const/cutoff
    tail_bound = (abs(1-np.exp(1j*w*tt))*PI/(9*abs(w)*cutoff**3)
                  +4*abs(const)*(1+4*PI**2/(abs(w)**2*(tt+cutoff)**2))
                  *np.exp(-w.imag*cutoff)/(w.imag*cutoff**2))
    return -w*w/8+x*x/8+a*integral, a*tail_bound


def continuum_real(delta, w, x, cells=160):
    """Real-axis PV with explicit pole subtraction and separate smooth part."""
    tt = PI/x
    start = KAPPA*delta
    ss, cc = np.sin(w*tt/2), np.cos(w*tt/2)
    def g(dd):
        return (w-4*PI**2/(w*(tt+dd)**2))*ss**2/(6*PI*dd**2)
    def gp(dd):
        return 8*PI**2*ss**2/(6*PI*w*dd**2*(tt+dd)**3)-2*g(dd)/dd
    smooth = quad(lambda v:(w-4*PI**2/(w*(tt+start/v)**2)),
                  0, 1, epsabs=1e-12, epsrel=1e-12)[0]/start
    total = -ss*cc*smooth/(6*PI)
    # A possible initial partial cell around the nonpositive n=0 pole.
    first = PI/w-tt
    if first > start:
        total += quad(lambda d:g(d)/np.tan(w*(tt+d)/2),start,first,
                      epsabs=1e-11,epsrel=1e-11)[0]
    for nn in range(1, cells+1):
        center = 2*PI*nn/w-tt
        lo = max(start, (2*nn-1)*PI/w-tt)
        hi = (2*nn+1)*PI/w-tt
        if hi <= lo:
            continue
        gc = g(center)
        def regular(dd):
            diff = dd-center
            if abs(diff) < 1e-7:
                return 2*gp(center)/w
            return (g(dd)-gc)/np.tan(w*diff/2)
        total += quad(regular,lo,hi,points=[center] if lo<center<hi else None,
                      epsabs=3e-11,epsrel=3e-10)[0]
        total += 2*gc/w*np.log(abs(np.sin(w*(hi-center)/2)/np.sin(w*(lo-center)/2)))
    return -w*w/8+x*x/8+start*total


def continuum_center(delta, cells=160):
    x = brentq(lambda x:force(x,delta),0.05,1.01,xtol=5e-15)
    w = brentq(lambda w:continuum_real(delta,w,x,cells),0.6,1.01,xtol=1e-12)
    return x,w


def real_kernel_and_derivative(ds, x, w):
    tt = PI/x
    ll = tt+ds
    ss, cc = np.sin(w*tt/2), np.cos(w*tt/2)
    phase = w*ll/2
    cot = 1/np.tan(phase)
    rr = ss*cc-ss*ss*cot
    rp = tt/2*np.cos(w*tt)-tt*ss*cc*cot+ll/2*ss**2/np.sin(phase)**2
    kk = 4*PI**2/ll**2
    value = -w*w/8+x*x/8-np.mean((w-kk/w)*rr)/(6*PI)
    deriv = -w/4-np.mean((1+kk/w**2)*rr+(w-kk/w)*rp)/(6*PI)
    return value,deriv


def spectral_time_case(q, delta, cutoff=1.5, lower=0.4):
    x = equilibrium(q,delta)
    ds = lengths(q,delta)
    ll = PI/x+ds
    poles = np.sort(np.concatenate([2*PI*np.arange(2,int(cutoff*l/(2*PI))+1)/l
                                   for l in ll]))
    # Select a fixed frequency band; positivity bounds ALL omitted atoms,
    # including very soft roots too close to bath poles for double precision.
    poles=poles[poles>lower]
    # The chosen examples have no coincident positive-weight poles.
    assert len(poles)>0 and np.min(np.diff(poles))>1e-11
    roots=[];mass=[]
    edges=np.r_[lower,poles,cutoff]
    for lo,hi in zip(edges[:-1],edges[1:]):
        gap=hi-lo
        a=max(1e-7,lo+max(1e-13,gap*1e-9))
        b=hi-max(1e-13,gap*1e-9)
        fa=real_kernel_and_derivative(ds,x,a)[0]
        fb=real_kernel_and_derivative(ds,x,b)[0]
        if fa*fb>=0:
            # The first/last partial interval need not contain a zero.
            assert hi==cutoff or lo==lower, (q,lo,hi,fa,fb)
            continue
        root=brentq(lambda w:real_kernel_and_derivative(ds,x,w)[0],a,b,
                    xtol=5e-15,rtol=1e-14)
        deriv=real_kernel_and_derivative(ds,x,root)[1]
        weight=-root/(4*deriv) # -2 M omega / K'(omega), per-channel K used.
        assert weight>0
        roots.append(root);mass.append(weight)
    roots=np.array(roots);mass=np.array(mass)
    total=float(np.sum(mass)); assert 0<total<1+2e-8
    _,center=continuum_center(delta)
    times=np.linspace(0,3/GAMMA,151)
    envelope=np.exp(1j*center*times/delta)*np.sum(
        mass[:,None]*np.exp(-1j*roots[:,None]*times/delta),axis=0)
    error=np.abs(envelope-np.exp(-GAMMA*times))
    transform_checks=[]
    for probe in [.8+.07j,1.1+.13j]:
        partial=np.sum(mass/(roots**2-probe**2))
        exact=1/(8*kernel(q,delta,probe,x))
        bound=(1-total)/abs((probe**2).imag)
        assert abs(exact-partial)<=bound+1e-9
        transform_checks.append(dict(probe=[probe.real,probe.imag],
                                     error=float(abs(exact-partial)),positive_tail_bound=float(bound)))
    k0=x*x/8+np.mean(PI*(PI/x)*ds/(3*ll**3))
    inverse_moment=1/(8*k0)
    missing_inverse_moment=inverse_moment-np.sum(mass/roots**2)
    assert missing_inverse_moment>=-2e-8
    h_tail=np.sqrt(max(0,1-total)*max(0,missing_inverse_moment))
    hh=np.sum((mass/roots)[:,None]*np.sin(roots[:,None]*times/delta),axis=0)
    h_target=np.exp(-GAMMA*times)*np.sin(center*times/delta)/center
    return dict(q=q,delta=delta,frequency_band=[lower,cutoff],poles=len(poles),clock_atoms=len(roots),
                captured_probability=total,omitted_probability_bound=1-total,
                max_envelope_error_with_tail_bound=float(np.max(error)+1-total),
                velocity_response_tail_bound=float(h_tail),
                max_velocity_response_error_with_tail_bound=float(np.max(abs(hh-h_target))+h_tail),
                sampled_lifetime_checks=[dict(lifetimes=m,
                    envelope_error_with_tail_bound=float(np.max(error[times<=m/GAMMA+1e-10])+1-total),
                    velocity_error_with_tail_bound=float(np.max(abs(hh-h_target)[times<=m/GAMMA+1e-10])+h_tail))
                    for m in [1,2,3]],
                transform_checks=transform_checks,
                rescaled_times=times.tolist(),complex_envelope=[[float(z.real),float(z.imag)] for z in envelope])


def main():
    exact={}
    w,e,W,M,K0=sy.symbols('w e W M K0',positive=True)
    h2=2*W/(M*e**3)
    ex=K0-M*w*w-M*w*w*e*e*h2/(e*e-w*w)
    exact['Schur_subtracted_response']=str(sy.simplify(ex-(K0-M*w*w+2*w*w*W/(e*(w*w-e*e)))))
    a,b=sy.symbols('a b')
    exact['return_length_derivative']=str(sy.simplify(sy.diff((1-b)/(1-a*b),b)+(1-a)/(1-a*b)**2))
    exact['throat_length_derivative']=str(sy.simplify(sy.diff((1-a)*(1-b)/(1-a*b),a)+(1-b)**2/(1-a*b)**2))
    om,z,dd=sy.symbols('omega z delta')
    exact['scaled_Stieltjes_partial_fraction']=str(sy.simplify(dd/(om-w)-2*w*dd/(om*om-w*w)-dd/(om+w)))
    assert all(v=='0' for v in exact.values())

    oscillator=[]
    for q,delta in [(3,.08),(7,.2)]:
        x=equilibrium(q,delta);ds=lengths(q,delta);tt=PI/x;ll=tt+ds
        k0=x*x/8+np.mean(PI*tt*ds/(3*ll**3))
        for nmax in [16,48]:
            nn=np.arange(2,nmax+1)[None,:];ll2=ll[:,None]
            ee=(2*PI*nn/ll2).ravel()
            ww=(2*(nn*nn-1)*np.sin(PI*nn*tt/ll2)**2/(3*ll2**2*nn)).ravel()
            hh=np.sqrt(2*ww/((q/8)*ee**3))
            matrix=np.diag(np.r_[8*k0,ee**2]);matrix[0,0]+=np.sum(ee**2*hh**2)
            matrix[0,1:]=matrix[1:,0]=-ee**2*hh
            ev,vec=eigh(matrix);weights=vec[0,:]**2
            assert np.min(ev)>0 and abs(sum(weights)-1)<3e-13
            for omega in [.87+.13j,1.03+.07j]:
                spectral=np.sum(weights/(ev-omega**2))
                direct=k0-omega**2/8+np.sum(2*omega**2*ww/(ee*(omega**2-ee**2)))/q
                err=abs(spectral-1/(8*direct));assert err<2e-9
                closed=kernel(q,delta,omega,x)
                assert nmax>max(abs(omega)*ll/(PI))
                tail=abs(omega)**2*np.mean(ll)/(6*PI**3)*4/(3*nmax)
                assert abs(direct-closed)<=tail
                oscillator.append(dict(q=q,delta=delta,nmax=nmax,omega=[omega.real,omega.imag],
                     resolvent_error=float(err),closed_kernel_error=float(abs(direct-closed)),tail_bound=float(tail),
                     inverse_moment_error=float(abs(np.sum(weights/ev)-1/(8*k0)))))
            assert abs(np.sum(weights/ev)-1/(8*k0))<2e-9

    quadrature=[];centers=[];scaling=[]
    for delta in [.08,.04,.02]:
        x,center=continuum_center(delta,100)
        _,fine=continuum_center(delta,200)
        assert abs(center-fine)<2e-7
        centers.append(dict(delta=delta,x=x,center=fine,pv_refinement=abs(center-fine)))
        for z0 in [.15+.3j,-.2+.6j]:
            omega=fine+delta*z0
            cv0,b0=continuum_complex(delta,omega,x,160)
            cv,bound=continuum_complex(delta,omega,x,320)
            refine=abs(cv-cv0)
            assert refine<3e-9
            for q in [1024,8192,65536]:
                xq=equilibrium(q,delta)
                assert abs(xq-x)<=2/(3*q)+1e-13
                fixed=kernel(q,delta,omega,x)
                err=abs(fixed-cv)
                vb=variation_bound(omega,PI/x)/(2*q)
                assert err<=vb+bound+refine
                actual=kernel(q,delta,omega,xq)
                quadrature.append(dict(q=q,delta=delta,z=[z0.real,z0.imag],root_difference=abs(xq-x),
                     fixed_root_kernel_error=err,variation_bound=vb,continuum_refinement=refine,
                     continuum_tail_bound=bound,total_kernel_error=abs(actual-cv)))
            scaling.append(dict(delta=delta,z=[z0.real,z0.imag],
                continuum_scaled_kernel=[float((8*cv/delta).real),float((8*cv/delta).imag)],
                predicted=[float((-2*(z0+1j*GAMMA)).real),float((-2*(z0+1j*GAMMA)).imag)],
                error=float(abs(8*cv/delta+2*(z0+1j*GAMMA)))))
    # Analytic endpoint counterexample: rho(y)=(1+|y|)/3 on [-1,1].
    # It is a probability density and is Lipschitz near E0=0 with rho(0)>0.
    endpoint=[]
    for eta in [1e-2,1e-3,1e-4,1e-5]:
        exact_im=(2*np.arctan(1/eta)+eta*np.log((1+eta*eta)/(eta*eta)))/3
        numeric=2*quad(lambda y:eta*(1+y)/(3*(y*y+eta*eta)),0,1,
                       points=[eta,10*eta],epsabs=1e-12,epsrel=1e-12)[0]
        assert abs(numeric-exact_im)<1e-11
        endpoint.append(dict(eta=eta,poisson_error=abs(exact_im-PI/3),
                             error_over_eta=abs(exact_im-PI/3)/eta,quadrature_error=abs(numeric-exact_im)))
    assert all(a['error_over_eta']<b['error_over_eta'] for a,b in zip(endpoint,endpoint[1:]))
    print('Algebra, oscillator, quadrature and source-endpoint checks passed.',flush=True)
    times=[]
    for q in [128,512,2048,8192]:
        row=spectral_time_case(q,.08);times.append(row)
        print(json.dumps({k:row[k] for k in ['q','poles','captured_probability','max_envelope_error_with_tail_bound']}),flush=True)
    assert times[-1]['max_envelope_error_with_tail_bound']<times[0]['max_envelope_error_with_tail_bound']
    assert times[-1]['sampled_lifetime_checks'][0]['envelope_error_with_tail_bound']<.1
    receipt=dict(status='passed_finite_channel_decay_checks',date='2026-10-11',
        note_sha256=sha(ROOT/'notes/finite-channel-decay.md'),script_sha256=sha(Path(__file__)),
        dependency_hashes={f:sha(ROOT/f) for f in ['notes/transparent-clock-flux.md','notes/collective-clock-response.md']},
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sy.__version__),
        gamma=GAMMA,exact_checks=exact,oscillator_cases=oscillator,quadrature_cases=quadrature,
        continuum_centers=centers,scaling_cases=scaling,source_endpoint_cases=endpoint,spectral_time_cases=times,
        finite_free_inverse_tail_decay_interval_proved=True,condition='delta -> 0 and q delta^3 -> infinity',
        fixed_finite_system_recurrent=True,quantitative_time_error_threshold_proved=False,
        time_numerics_are_sampled_not_uniform_certificates=True,
        physical_damping_proved=False,interacting_stress_bound_closed=False,actual_dipole_full_response_transferred=False,
        nonlinear_time_window_controlled=False,novelty_established=False,lean_invocations=0)
    (ROOT/'evidence/finite-channel-decay.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'],gamma=GAMMA,counts={k:len(receipt[k]) for k in
        ['exact_checks','oscillator_cases','quadrature_cases','continuum_centers','scaling_cases','source_endpoint_cases','spectral_time_cases']})))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__=='__main__':
    main()
