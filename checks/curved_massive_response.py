"""S072: fixed global-AdS2 massive current measure and quadratic work.

Checks exact Jacobi modes, the Ward matrix identity, the Bessel heat
comparison, complete finite-energy pair sums and quantum source ordering.
It does not evolve an RN geometry or prove a full interacting remainder.
"""
from pathlib import Path
import hashlib
import json
import platform

import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad
from scipy.special import eval_jacobi, gammaln, roots_legendre, ive, jv

ROOT=Path(__file__).resolve().parents[1]
A=16/(9*np.sqrt(np.pi))*(2*np.e/3)**1.5
B=3*np.pi**2*A/8
PULSE_PRODUCT=(32/35)*(384/(25*np.sqrt(5)))


def modes(M,count,r,derivatives=False):
    s=np.sin(r);c=np.cos(r);F=[];G=[];FD=[];GD=[]
    for j in range(count):
        norm=np.exp((gammaln(j+1)+gammaln(j+2*M+1))/2
                    -(M+.5)*np.log(2)-gammaln(j+M+.5))
        pf=norm*c**M*np.sqrt(1+s);pg=norm*c**M*np.sqrt(1-s)
        p=eval_jacobi(j,M-.5,M+.5,s);q=eval_jacobi(j,M+.5,M-.5,s)
        F.append(pf*p);G.append(pg*q)
        if derivatives:
            dp=0 if j==0 else .5*(j+2*M+1)*eval_jacobi(j-1,M+.5,M+1.5,s)
            dq=0 if j==0 else .5*(j+2*M+1)*eval_jacobi(j-1,M+1.5,M+.5,s)
            FD.append(pf*((-M*np.tan(r)+c/(2*(1+s)))*p+c*dp))
            GD.append(pg*((-M*np.tan(r)-c/(2*(1-s)))*q+c*dq))
    result=(np.asarray(F),np.asarray(G),M+.5+np.arange(count))
    return result+(np.asarray(FD),np.asarray(GD)) if derivatives else result


def grid(order):
    x,w=roots_legendre(order)
    return x*np.pi/2,w*np.pi/2


def profile(r,kind):
    if kind=='linear':return r,np.ones_like(r)
    if kind=='sine':return np.sin(r),np.cos(r)
    epsilon=float(kind)
    c=np.cos(r)
    return np.sin(r)/(epsilon+c),(1+epsilon*c)/(epsilon+c)**2


def mode_checks():
    r,w=grid(320);rows=[]
    for M in [4.,4.5,7.2,18.04813031688602]:
        F,G,E,FD,GD=modes(M,18,r,True)
        gram=(F*w)@F.T+(G*w)@G.T
        residual=np.sqrt(np.sum(w*((-GD+M/np.cos(r)*G-E[:,None]*F)**2
                                  +(FD+M/np.cos(r)*F-E[:,None]*G)**2),axis=1))
        error=float(np.max(abs(gram-np.eye(len(E)))))
        assert error<2e-11 and max(residual)<2e-10
        ward=[]
        for kind in ['linear','sine','.1','.025']:
            V,D=profile(r,kind)
            pair=(F*(w*V))@F.T-(G*(w*V))@G.T
            comm=(F*(w*D))@G.T+(G*(w*D))@F.T
            err=float(np.max(abs((E[:,None]+E[None,:])*pair-comm)))
            assert err<2e-9
            ward.append(dict(profile=kind,maximum_error=err))
        pure_gauge=(F*w)@F.T-(G*w)@G.T
        assert np.max(abs(pure_gauge))<2e-11
        rows.append(dict(M=M,gram_error=error,dirac_residual=float(max(residual)),ward=ward,
                         constant_potential_pair_error=float(np.max(abs(pure_gauge)))))
    return rows


def heat_checks():
    rows=[];weber=[]
    for M,delta,t in [(4.,1.,.05),(7.2,1.4,.2),(18.048,1.2,.02)]:
        direct,error=quad(lambda k:delta*k*np.exp(-t*k*k)*jv(M-.5,delta*k)**2,
                          0,np.inf,epsabs=2e-14,epsrel=2e-12)
        exact=delta/(2*t)*ive(M-.5,delta*delta/(2*t))
        assert abs(direct/exact-1)<2e-10
        weber.append(dict(M=M,delta=delta,t=t,relative_error=float(abs(direct/exact-1))))
    for M in [4.,7.2,18.04813031688602]:
        r=np.array([0.,.8,1.3,1.52]);delta=np.pi/2-abs(r)
        F,G,E=modes(M,160,r)
        for t in [.002,.02,.2]:
            partial=np.sum(np.exp(-t*E[:,None]**2)*(F*F+G*G),axis=0)
            z=delta*delta/(2*t)
            bessel=delta/(2*t)*ive(M-.5,z)
            rational=(4*np.pi*t)**-.5*(1+(M-1)*t/delta**2)**(-M)
            power=(4*np.pi*t)**-.5*(1+3*M*M*t/(16*delta**2))**(-4)
            assert np.all(partial<=bessel*(1+2e-10))
            assert np.all(bessel<=rational*(1+2e-10))
            assert np.all(rational<=power*(1+2e-10))
            rows.append(dict(M=M,t=t,max_partial_over_bessel=float(max(partial/bessel)),
                             max_bessel_over_rational=float(max(bessel/rational))))
    return dict(weber_normalization=weber,comparison=rows)


def tower_checks():
    rows=[]
    for q,C in [(1,.1),(31,1.),(127,4.)]:
        assert np.sqrt((q+1)/C)>=4
        n=np.arange(1,100001,dtype=float);M2=n*(n+q)/C;d=2*(q+2*n)
        fractions=[]
        for a in [.001,.1,10.]:
            partial=np.sum(d/(1+a*M2)**4)
            bound=2*C*(1+1/q)/(3*a)
            assert partial<=bound
            fractions.append(float(partial/bound))
        rows.append(dict(q=q,C=C,majorant_fractions=fractions))
    return rows


def pulse_transform(y):
    y=np.asarray(y);out=np.empty_like(y,dtype=float);small=abs(y)<1e-6
    out[small]=32/35*(1-y[small]**2/18)
    x=abs(y[~small]);out[~small]=6*np.sqrt(np.pi)*(2/x)**3.5*jv(3.5,x)
    return out


def pair_sum(q,C,kind,tau,cutoff,order):
    r,w=grid(order);V,D=profile(r,kind)
    Q=float(np.sum(w*np.cos(r)**2*D*D))
    # Normalize the radial Maxwell seminorm to one.
    V=V/np.sqrt(Q);D=D/np.sqrt(Q)
    top=C*(cutoff/2-.5)**2
    N=int(np.floor(2*top/(q+np.sqrt(q*q+4*top))))
    totals=[];energies=[];work=0.
    for n in range(1,N+1):
        M=np.sqrt(n*(n+q)/C);count=int(np.floor(cutoff-2*M-1))+1
        if count<=0:continue
        F,G,E=modes(M,count,r)
        pair=(F*(w*V))@F.T-(G*(w*V))@G.T
        en=E[:,None]+E[None,:];keep=en<=cutoff
        measure=2*(q+2*n)*en[keep]*pair[keep]**2
        omega=en[keep]
        energies.extend(omega.tolist());totals.extend(measure.tolist())
        work+=np.sum(measure*(tau*pulse_transform(omega*tau))**2)
    energies=np.array(energies);totals=np.array(totals)
    Dcoef=B*C*(1+1/q)
    cumulative=[]
    for omega in [20.,40.,cutoff]:
        value=float(np.sum(totals[energies<=omega]));bound=Dcoef*omega**2
        assert value<=bound
        cumulative.append(dict(omega=omega,value=value,bound=bound))
    y=cutoff*tau;envelope=96*(1+6/y+15/y**2+15/y**3)
    tail=4*Dcoef*envelope**2/(3*tau**6*cutoff**6)
    full_bound=2*Dcoef*PULSE_PRODUCT
    assert work+tail<full_bound
    return dict(q=q,C=C,profile=kind,tau=tau,cutoff=cutoff,radial_normalization=Q,
                finite_pair_count=len(totals),computed_work=float(work),tail_bound=float(tail),
                full_work_bound=float(full_bound),cumulative=cumulative)


def pulse_checks():
    rows=[]
    for q,C,kind,tau in [(31,1.,'linear',.2),(127,4.,'.05',.3),(3,.1,'sine',.4)]:
        base=pair_sum(q,C,kind,tau,60.,256)
        refined=pair_sum(q,C,kind,tau,60.,512)
        error=abs(base['computed_work']-refined['computed_work'])
        assert error<2e-9
        refined['quadrature_refinement_error']=error
        rows.append(refined)
    return rows


def quantum_source_checks():
    # One pair sector of energy E coupled to a two-level quantum source.
    # The full source+pair Hamiltonian splits into two invariant 2x2 blocks.
    E=10.;gap=7.;v=.2;tau=.3;rows=[]
    for initial in ['excited','ground']:
        detuning=E-gap if initial=='excited' else E+gap
        target=E*v*v*(tau*float(pulse_transform(np.array([detuning*tau]))[0]))**2
        errors=[];works=[]
        for eps in [.1,.05,.025]:
            def rhs(t,state):
                f=(1-(t/tau)**2)**3
                H=np.array([[0,eps*v*f],[eps*v*f,detuning]])
                return -1j*H@state
            sol=solve_ivp(rhs,(-tau,tau),np.array([1.+0j,0j]),method='DOP853',rtol=2e-12,atol=2e-14)
            assert sol.success
            work=E*abs(sol.y[1,-1])**2
            errors.append(abs(work/eps**2-target));works.append(float(work))
        assert 3.8<errors[0]/errors[1]<4.2 and 3.8<errors[1]/errors[2]<4.2
        rows.append(dict(initial_source=initial,ordered_second_order_coefficient=float(target),
                         epsilon=[.1,.05,.025],exact_bath_work=works,coefficient_errors=errors))
    assert rows[0]['ordered_second_order_coefficient']>rows[1]['ordered_second_order_coefficient']*10
    return rows


def main():
    result=dict(date='2026-10-13',status='PASS',python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                constants=dict(local_projector=A,pair_measure=B,separable_pulse=B/np.pi*PULSE_PRODUCT,
                               general_dyadic=8/np.sqrt(3)),
                modes=mode_checks(),heat=heat_checks(),tower=tower_checks(),pulses=pulse_checks(),
                quantum_source=quantum_source_checks(),
                candidate=[dict(k=k,separable_pulse_bound=B/np.pi*PULSE_PRODUCT*k**-8*(1+k**-4)) for k in [2000,4000,8000]],
                limitations=['Fixed global AdS2 x S2, normalizable massive vacuum, M1 >= 4',
                             'Quadratic bath work; quantum extension assumes initially factorized bath vacuum and finite ordered time-derivative moments',
                             'No evolved RN exterior, interacting remainder, actual source moment bound or complete stress'])
    (ROOT/'evidence/curved-massive-response.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: exact modes, Ward identity, heat comparison, massive tower, pulse work and ordered quantum source.')
    print('Constants:',result['constants'])
    print('Candidate:',result['candidate'][0])
    print('Pulse refinements:',[r['quadrature_refinement_error'] for r in result['pulses']])


if __name__=='__main__':main()
