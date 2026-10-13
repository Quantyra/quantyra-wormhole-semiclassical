"""S072: ordered Gaussian electric vacuum coupled to the massive tower.

Fixed global product, normalizable charged scalar, quadratic cross-sector
response. Tests do not identify the coupling ramp with the physical junction.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import platform

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.special import eval_gegenbauer, gammaln, roots_jacobi, jv

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('curved_massive_response', ROOT/'checks/curved_massive_response.py')
_curved = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_curved)
B, modes, grid, pulse_transform = _curved.B, _curved.modes, _curved.grid, _curved.pulse_transform


def scalar_modes(delta, count, r, derivatives=False):
    s = np.sin(r); c = np.cos(r)
    values = []; second = []
    for j in range(count):
        logh = (np.log(np.pi) + (1-2*delta)*np.log(2)
                + gammaln(j+2*delta)-gammaln(j+1)-np.log(j+delta)-2*gammaln(delta))
        factor = np.exp(-logh/2)*c**delta
        P = eval_gegenbauer(j, delta, s)
        values.append(factor*P)
        if derivatives:
            P1 = 0 if j == 0 else 2*delta*eval_gegenbauer(j-1, delta+1, s)
            P2 = 0 if j < 2 else 4*delta*(delta+1)*eval_gegenbauer(j-2, delta+2, s)
            second.append(factor*((delta**2*np.tan(r)**2-delta/c**2)*P
                                  -(2*delta+1)*s*P1+c*c*P2))
    return np.array(values), np.array(second)


def scalar_checks():
    rows = []
    for delta in [1., 1.2, 2.3, 4.]:
        count = 24
        r, w = grid(320)
        v, dd = scalar_modes(delta, count, r, True)
        omega = delta+np.arange(count)
        residual = np.sqrt(np.sum(w*(-dd+delta*(delta-1)/np.cos(r)**2*v-omega[:,None]**2*v)**2,axis=1))
        assert max(residual) < 2e-9
        # Gauss-Jacobi integrals independently remove the endpoint powers.
        x, wx = roots_jacobi(128, delta-.5, delta-.5)
        pol = scalar_modes(delta, count, np.arcsin(x))[0]/(1-x*x)**(delta/2)
        gram = (pol*wx)@pol.T
        x, wx = roots_jacobi(128, delta-1.5, delta-1.5)
        pol = scalar_modes(delta, count, np.arcsin(x))[0]/(1-x*x)**(delta/2)
        hf = np.sum(wx*pol*pol,axis=1)
        target = 2*omega/(2*delta-1)
        assert np.max(abs(gram-np.eye(count))) < 2e-10
        assert np.max(abs(hf/target-1)) < 2e-10
        s = delta*(delta-1); step = 1e-5
        def energy2(m2):return (np.arange(count)+(1+np.sqrt(1+4*m2))/2)**2
        fd = (energy2(s+step)-energy2(s-step))/(2*step)
        assert np.max(abs(fd/target-1)) < 2e-8
        rows.append(dict(delta=delta,scalar_residual=float(max(residual)),
                         gram_error=float(np.max(abs(gram-np.eye(count)))),
                         hellmann_feynman_relative_error=float(np.max(abs(hf/target-1)))))
    return rows


def ultraviolet_checks():
    delta = 1.2; rows = []
    for count in [10,100,1000,10000]:
        omega = delta+np.arange(count,dtype=float)
        # Divide out K=s^2/[lambda(2delta-1)]. These are positive vacuum norms.
        rows.append(dict(cutoff_modes=count,equal_time_norm_squared_over_K=count,
                         second_derivative_norm_squared_over_K=float(np.sum(omega**4)),
                         fourth_sum_over_N5=float(np.sum(omega**4)/count**5)))
    assert abs(rows[-1]['fourth_sum_over_N5']-.2)<1e-4
    return rows


def pulse_data(T):
    F0 = 32*T/35; F2 = 384/(25*np.sqrt(5)*T)
    return F0,F2


def ramp_transform(y):
    # Magnitude up to the harmless midpoint phase exp(iy/2).
    y=np.asarray(y);out=np.empty_like(y,dtype=float);small=abs(y)<1e-6
    out[small]=1-y[small]**2/56
    x=abs(y[~small])/2
    out[~small]=15/8*np.sqrt(np.pi)*(2/x)**2.5*jv(2.5,x)
    return out


def ramp_normalization_checks():
    roots=[0.,(3-np.sqrt(3))/6,(3+np.sqrt(3))/6,1.]
    norm=sum(quad(lambda x:abs(60-360*x+360*x*x),a,b,epsabs=1e-12)[0]
             for a,b in zip(roots[:-1],roots[1:]))
    assert abs(norm-40/np.sqrt(3))<2e-12
    errors=[]
    for y in [.2,3.,12.]:
        direct=quad(lambda x:30*x*x*(1-x)**2*np.cos(y*(x-.5)),0,1,epsabs=2e-14)[0]
        error=abs(direct-float(ramp_transform(np.array([y]))[0]));assert error<2e-14
        errors.append(error)
    return dict(third_derivative_norm=float(norm),fourier_errors=errors)


def pair_source_sum(q,C,delta,T,cutoff,order):
    r,w = grid(order); source_count = int(np.floor(cutoff-delta))+1
    v,_ = scalar_modes(delta,source_count,r)
    boson = delta+np.arange(source_count)
    # Each profile below has exactly Q=1 by Hellmann-Feynman.
    derivative = v/np.cos(r)**2*np.sqrt((2*delta-1)/(2*boson[:,None]))
    top = C*((cutoff-delta)/2-.5)**2
    N = int(np.floor(2*top/(q+np.sqrt(q*q+4*top))))
    total_work = 0.; massive_work = 0.; ramp_work = 0.; triples = 0; cumulative = [0.,0.,0.]
    bands = [cutoff*.6,cutoff*.8,cutoff]
    for n in range(1,N+1):
        M = np.sqrt(n*(n+q)/C)
        count = int(np.floor(cutoff-delta-2*M-1))+1
        if count <= 0:continue
        F,G,E = modes(M,count,r); pair_E = E[:,None]+E[None,:]
        for j,omega in enumerate(boson):
            U = pair_E+omega; keep = U<=cutoff
            if not np.any(keep):break
            comm = (F*(w*derivative[j]))@G.T+(G*(w*derivative[j]))@F.T
            probability = 2*(q+2*n)*(comm[keep]/pair_E[keep])**2
            transform = T*pulse_transform(T*U[keep])
            total_work += np.sum(U[keep]*probability*transform**2)
            massive_work += np.sum(pair_E[keep]*probability*transform**2)
            ramp_work += np.sum(probability/U[keep]*ramp_transform(T*U[keep])**2)
            triples += int(np.sum(keep))
            for k,band in enumerate(bands):
                selected = U[keep]<=band
                cumulative[k] += float(np.sum(U[keep][selected]*probability[selected]))
    D = B*C*(1+1/q); F0,F2 = pulse_data(T)
    total_bound = 8*D/3*np.sqrt(F0)*F2**1.5
    massive_bound = total_bound/2
    y = cutoff*T; envelope = 96*(1+6/y+15/y**2+15/y**3)
    tail_massive = 8*D*envelope**2/(15*T**6*cutoff**5)
    tail_total = 2*tail_massive
    G2=40/(np.sqrt(3)*T*T);gap=2*np.sqrt((q+1)/C)+1+delta
    ramp_bound=min(8*D/3*np.sqrt(G2),4*D/3*G2**2/gap**3)
    ramp_tail=4*D/3*G2**2/cutoff**3
    for value,band in zip(cumulative,bands):assert value <= 2*D*band**3/3
    assert total_work+tail_total < total_bound
    assert massive_work+tail_massive < massive_bound
    assert ramp_work+ramp_tail < ramp_bound
    return dict(q=q,C=C,delta=delta,T=T,cutoff=cutoff,allowed_triples=triples,
                work_over_K=float(total_work),massive_work_over_K=float(massive_work),
                tail_over_K=float(tail_total),massive_tail_over_K=float(tail_massive),
                bound_over_K=float(total_bound),massive_bound_over_K=float(massive_bound),
                ramp_excess_over_K=float(ramp_work),ramp_tail_over_K=float(ramp_tail),ramp_bound_over_K=float(ramp_bound),
                cumulative_total_over_K=cumulative,cumulative_cutoffs=bands)


def vacuum_checks():
    rows = []
    for q,C,delta,T in [(31,1.,1.,.3),(127,4.,1.2,.4),(3,.1,2.3,.5)]:
        coarse = pair_source_sum(q,C,delta,T,40.,192)
        fine = pair_source_sum(q,C,delta,T,40.,384)
        err = abs(coarse['work_over_K']-fine['work_over_K'])
        assert err < 2e-9
        fine['quadrature_refinement_error'] = err
        rows.append(fine)
    return rows


def thermal_checks():
    rows = []; delta = 1.2; T=.3; F0,F2=pulse_data(T)
    for beta in [.5,2.,2*np.pi]:
        omega = delta+np.arange(2000)
        # Beyond omega=1400 the exponential is below floating-point relevance.
        z = np.exp(-beta*omega); occ = z/(1-z)
        moments = [float(np.sum(occ*omega**r)) for r in range(3)]
        ell = np.arange(1,2001); z = np.exp(-beta*ell); pref = np.exp(-beta*ell*delta)
        geometric = [float(np.sum(pref/(1-z))),
                     float(np.sum(pref*(delta/(1-z)+z/(1-z)**2))),
                     float(np.sum(pref*(delta**2/(1-z)+2*delta*z/(1-z)**2+z*(1+z)/(1-z)**3)))]
        assert max(abs(np.array(moments)-geometric)) < 2e-12
        bound = 4*F0*F2*moments[0]+8/3*F0**1.5*np.sqrt(F2)*moments[1]+F0**2*moments[2]
        rows.append(dict(beta=beta,occupation_moments=moments,massive_thermal_bound_over_DK=float(bound)))
    # Test the shifted envelope against the saturating positive measure dmu=2E dE.
    shifted = []
    for omega in [.3,3.,30.]:
        L0 = np.sqrt(F2/F0)
        def integrand(E):
            z=abs(E-omega)
            return 2*E*(F0**2 if z<=L0 else F2**2/z**4)
        cuts=sorted(set([0.,max(0.,omega-L0),omega+L0,np.inf]))
        direct=sum(quad(integrand,a,b,epsabs=2e-10,epsrel=2e-11)[0] for a,b in zip(cuts[:-1],cuts[1:]))
        bound=2*F0*F2+8/3*omega*F0**1.5*np.sqrt(F2)+omega**2*F0**2
        assert direct<=bound*(1+2e-11)
        shifted.append(dict(omega=omega,direct=float(direct),bound=float(bound)))
    return dict(moments=rows,shifted_envelope=shifted)


def dressed_ramp_checks():
    rows = []
    for U,T,start,finish in [(7.,.7,0.,1.),(13.,.5,.2,.9),(19.,.8,1.,.4)]:
        v=.3
        def f(t):
            x=t/T;return start+(finish-start)*(10*x**3-15*x**4+6*x**5)
        def fd(t):
            x=t/T;return (finish-start)*30*x*x*(1-x)**2/T
        transform = complex(quad(lambda t:fd(t)*np.cos(U*t),0,T,epsabs=2e-14)[0],
                            quad(lambda t:fd(t)*np.sin(U*t),0,T,epsabs=2e-14)[0])
        target=v*v/U*abs(transform)**2
        errors=[]; excess=[]
        for eps in [.16,.08,.04]:
            def H(t):return np.array([[0,eps*v*f(t)],[eps*v*f(t),U]])
            e0,u0=np.linalg.eigh(H(0));ef,uf=np.linalg.eigh(H(T))
            sol=solve_ivp(lambda t,y:-1j*H(t)@y,(0,T),u0[:,0].astype(complex),
                          method='DOP853',rtol=1e-12,atol=1e-14)
            assert sol.success
            # Positive final excitation avoids subtracting nearly equal ground energies.
            W=(ef[1]-ef[0])*abs(np.vdot(uf[:,1],sol.y[:,-1]))**2
            excess.append(float(W));errors.append(float(abs(W/eps**2-target)))
        assert 3.8<errors[0]/errors[1]<4.2 and 3.8<errors[1]/errors[2]<4.2
        rows.append(dict(U=U,T=T,start=start,finish=finish,second_order_target=float(target),
                         epsilon=[.16,.08,.04],exact_excess=excess,coefficient_errors=errors))
    return rows


def candidate_checks():
    rows=[]
    for k in [2000,4000,8000]:
        q=float(k)**4; g=float(k)**-4
        X=(np.sqrt(1+12*np.pi*k**6)-1)/6; C=1+6*X
        lam=q/np.pi; s=q*g*g/(4*np.pi**2*C)
        root=np.sqrt(1+4*s); delta=1+2*s/(1+root)
        K=s*s/(lam*root); D=B*C*(1+1/q)
        gap=2*np.sqrt((q+1)/C)+1+delta
        T=.25; F0,F2=pulse_data(T); G2=40/(np.sqrt(3)*T*T)
        pulse=8*D*K/3*np.sqrt(F0)*F2**1.5
        ramp=min(8*D*K/3*np.sqrt(G2),4*D*K/3*G2**2/gap**3)
        rows.append(dict(k=k,C=float(C),schwinger_mass_squared=float(s),K=float(K),
                         vacuum_pulse_total_bound=float(pulse),dressed_ramp_excess_bound=float(ramp),
                         T=T,threshold=float(gap)))
    return rows


def main():
    result=dict(date='2026-10-13',status='PASS',python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                dependency_sha256=hashlib.sha256((ROOT/'checks/curved_massive_response.py').read_bytes()).hexdigest(),
                scalar=scalar_checks(),ultraviolet=ultraviolet_checks(),vacuum=vacuum_checks(),
                thermal=thermal_checks(),ramp_normalization=ramp_normalization_checks(),
                dressed_ramp=dressed_ramp_checks(),candidate=candidate_checks(),
                limitations=['Fixed global Gaussian electric sector and free massive vacuum; no boundary compact rotor or driven junction',
                             'Quadratic cross-sector response only, not all massive Coulomb interactions or gravitational stress',
                             'Dressed ramp at finite regulator with a uniform energy bound; continuum interacting vacuum not constructed'])
    (ROOT/'evidence/ordered-electric-vacuum.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: constrained scalar modes, divergent raw moments, ordered full-spectrum work, thermal envelopes and dressed ramps.')
    print('Candidate:',result['candidate'][0])
    print('Refinements:',[x['quadrature_refinement_error'] for x in result['vacuum']])


if __name__=='__main__':main()
