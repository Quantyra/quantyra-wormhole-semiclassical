"""S072 projected phase action, positive bath and quantum recoil checks.

RN regulators verify the radial identities; they do not replace the
analytic half-line proof or compute the interacting fermion stress.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
import scipy
from scipy.linalg import solve, eigh, expm
from scipy.integrate import quad, solve_ivp
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]


def matrices(X,q,G,tau,cells,kind,s=0,primal=False):
    Y=1+3*X-tau;g2=np.pi*G*q*q/Y
    lam=q/np.pi if kind=='electric' else q*(q*q-1)/(12*np.pi)
    nf=1 if kind=='electric' else 2
    mu=8*np.pi/(3*g2);P=q/2
    grid=1-np.exp(np.linspace(0,np.log(1e-3),cells+1))
    size=nf*(cells+1)
    A=np.zeros((size,size));M=np.zeros_like(A)
    nodes,weights=leggauss(8)
    for i in range(cells):
        left,right=grid[i:i+2];h=right-left
        u=left+(nodes+1)*h/2;wt=weights*h/2
        basis=np.stack([(1-nodes)/2,(1+nodes)/2],axis=1)
        der=np.array([-1,1])/h
        F=(1-u)**2*(X+2*X*u+(1+3*X)*u*u)+tau*u**3*(1-u)
        pot=np.zeros((len(u),nf,nf))
        if nf==1:
            pot[:,0,0]=lam*g2/(4*np.pi)
        else:
            pot[:,0,0]=2+4*Y*u*u
            pot[:,0,1]=pot[:,1,0]=np.sqrt(mu*lam)/(q+1)*(-6*G*P*u*u+3*g2/(2*np.pi))
            pot[:,1,1]=lam/(q+1)**2*(6*G*u*u+3*g2/np.pi)
        for k in range(len(u)):
            if primal:
                grad=F[k]*np.linalg.inv(F[k]*pot[k]+s*s*np.eye(nf))
                local_pot=np.eye(nf)/F[k]
            else:
                grad=F[k]*np.eye(nf);local_pot=pot[k]
            for a in range(2):
                for b in range(2):
                    aa=slice(nf*(i+a),nf*(i+a+1));bb=slice(nf*(i+b),nf*(i+b+1))
                    A[aa,bb]+=wt[k]*(der[a]*der[b]*grad+basis[k,a]*basis[k,b]*local_pot)
                    M[aa,bb]+=wt[k]*basis[k,a]*basis[k,b]*np.eye(nf)/F[k]
    return A,M,lam,nf


def radial_case(X,q,G,tau,cells,kind):
    A,M,lam,nf=matrices(X,q,G,tau,cells,kind)
    source=np.zeros(len(A));port=nf-1;source[port]=1
    response=solve(A,source,assume_a='pos')
    z0=response/response[port];K0=1/response[port]
    keep=np.delete(np.arange(len(A)),port)
    Ai=A[np.ix_(keep,keep)];Mi=M[np.ix_(keep,keep)]
    ev,vec=eigh(Ai,Mi,driver='gvd')
    assert min(ev)>0
    overlap=vec.T@(M@z0)[keep]
    sigma=ev*overlap**2/lam
    # Finite elements give their boundary basis a nonzero L2 norm. This
    # produces a contact term absent in the continuum; retain, do not hide it.
    contact=(z0@M@z0-np.sum(overlap**2))/lam
    assert contact>=-2e-10
    rows=[]
    for ratio in (.15,.7,2):
        s=ratio*K0
        inverse=solve(A+s*s*M,source,assume_a='pos')[port]
        direct=1/(lam*s*s*inverse)
        spectral=K0/(lam*s*s)+contact+np.sum(sigma/(ev+s*s))
        spectral_error=abs(direct-spectral)/direct
        assert spectral_error<3e-8
        primal,_,_,_=matrices(X,q,G,tau,cells,kind,s=s,primal=True)
        boundary=np.r_[np.arange(nf),np.arange(len(A)-nf,len(A))]
        free=np.setdiff1d(np.arange(len(A)),boundary)
        b=np.zeros(len(A));b[port]=1
        b[free]=solve(primal[np.ix_(free,free)],-primal[np.ix_(free,boundary)]@b[boundary],assume_a='pos')
        primal_I=lam*(b@primal@b)
        dual_I=lam*inverse
        error=abs(primal_I/dual_I-1)
        rows.append(dict(s_over_K0=ratio,spectral_relative_error=spectral_error,
                         primal_dual_relative_error=error))
    return dict(X=X,q=q,G=G,tau=tau,kind=kind,cells=cells,
                horizon_cut_1_minus_u=.001,static_I=lam/K0,
                finite_element_contact=contact,minimum_bath_weight=float(min(sigma)),
                comparisons=rows)


def oscillator_checks():
    rows=[]
    n=100
    annih=np.diag(np.sqrt(np.arange(1,n)),1)
    for omega,beta,weight in ((.7,1.2,.04),(1.3,2,.5),(3,.8,1.2)):
        energies=omega*(np.arange(n)+.5)
        probs=np.exp(-beta*(energies-energies[0]));probs/=sum(probs)
        X=np.sqrt(weight)*(annih+annih.T)/np.sqrt(2*omega)
        V=expm(1j*X)
        kick=np.sum(probs*(np.sum(abs(V)**2*energies[:,None],axis=0)-energies))
        assert abs(kick-weight/2)<2e-11
        tests=[]
        for t in (0,.2,1):
            evolved=np.exp(1j*energies[:,None]*2*t)*V*np.exp(-1j*energies[None,:]*2*t)
            midpoint=np.sum(np.sqrt(probs[:,None]*probs[None,:])*evolved*V.conj())
            W=weight/(2*omega)*(1/np.tanh(beta*omega/2)-np.cos(2*omega*t)/np.sinh(beta*omega/2))
            err=abs(midpoint-np.exp(-W))
            assert err<3e-12
            # Ordinary Wightman function includes the commutator phase.
            value=np.sum(probs[None,:]*np.exp(-1j*(energies[:,None]-energies[None,:])*t)*abs(V)**2)
            J=weight/(2*omega)*((np.cos(omega*t)-1)/np.tanh(beta*omega/2)-1j*np.sin(omega*t))
            werr=abs(value-np.exp(J))
            assert werr<3e-12
            tests.append(dict(t=t,TFD_error=float(err),Wightman_error=float(werr)))
        rows.append(dict(omega=omega,beta=beta,weight=weight,
                         mean_recoil=kick,expected_recoil=weight/2,comparisons=tests))
    return rows


def measure_bounds():
    # d sigma = 2 B d omega saturates the cumulative envelope. It need
    # not saturate the stronger Stieltjes bound; it checks the integration.
    rows=[]
    for beta in (.3,2,10):
        for cutoff in (.2,4,100):
            B=.02
            W0=quad(lambda w:B*np.tanh(beta*w/4)/w,0,cutoff,epsabs=1e-12)[0]
            t=beta/(8*np.pi)
            extra=quad(lambda w:B*(1-np.cos(2*w*t))/(w*np.sinh(beta*w/2)) if beta*w<700 else 0,
                       0,cutoff,epsabs=1e-12)[0]
            mid_bound=B*(beta*cutoff/4 if beta*cutoff<=4 else 1+np.log(beta*cutoff/4))
            time_bound=2*np.pi**2*B*(t/beta)**2
            assert W0<=mid_bound+1e-12 and extra<=time_bound+1e-12
            rows.append(dict(beta=beta,cutoff=cutoff,W0=W0,W0_bound=mid_bound,
                             free_TFD_increment=extra,increment_bound=time_bound))
    s,w=sy.symbols('s w',positive=True)
    # An atom of weight 2 B w obeys the full Stieltjes envelope.
    assert sy.factor(1-2*w*s/(s*s+w*w)) == (s-w)**2/(s*s+w*w)
    assert abs(quad(lambda v:v/np.sinh(v) if v<700 else 0,0,np.inf)[0]-np.pi**2/4)<1e-12
    return rows


def interacting_current_check(levels,h_final):
    """A Gauss-constrained single-transfer block, not a 4D stress model."""
    omega=1.3;beta=2.;weight=.2;duration=1.2
    annih=np.diag(np.sqrt(np.arange(1,levels)),1)
    boson_H=omega*(annih.T@annih+.5*np.eye(levels))
    X=np.sqrt(weight)*(annih+annih.T)/np.sqrt(2*omega)
    vertex=expm(1j*X)
    raise_charge=np.array([[0.,0.],[1.,0.]])
    transfer=np.kron(raise_charge,vertex)
    interaction=transfer+transfer.conj().T
    Q=np.kron(np.diag([0.,1.]),np.eye(levels))
    # Q=n_R eliminates the independent charge rotor. The diagonal energy
    # is epsilon n_R+Q^2/(2 I_0), with I_0=4 and epsilon=.3.
    H0=np.kron(np.eye(2),boson_H)+.425*Q
    a=np.kron(np.eye(2),annih)
    dim=len(Q);total=dim*dim
    def rhs(t,y):
        U=y[:total].reshape(dim,dim)
        h=h_final*(3*(t/duration)**2-2*(t/duration)**3)
        H=H0+h*interaction
        j=1j*h*(-transfer+transfer.conj().T)
        return np.r_[(-1j*H@U).ravel(),(np.exp(1j*omega*t)*(U.conj().T@j@U)).ravel()]
    init=np.r_[np.eye(dim,dtype=complex).ravel(),np.zeros(total,dtype=complex)]
    result=solve_ivp(rhs,(0,duration),init,method='DOP853',rtol=2e-11,atol=2e-12)
    assert result.success
    U=result.y[:total,-1].reshape(dim,dim)
    current_integral=result.y[total:,-1].reshape(dim,dim)
    probs=np.exp(-beta*omega*np.arange(levels));probs/=sum(probs)
    rho=np.kron(np.ones((2,2))/2,np.diag(probs))
    rho_half=np.kron(np.ones((2,2))/2,np.diag(np.sqrt(probs)))
    free=np.exp(-1j*omega*duration)*a
    delta=1j*np.sqrt(weight/(2*omega))*np.exp(-1j*omega*duration)*current_integral
    evolved=U.conj().T@a@U
    state_error=np.linalg.norm((evolved-free-delta)@rho_half)
    assert state_error<2e-7
    actual=omega*np.trace(rho@(evolved.conj().T@evolved-a.conj().T@a)).real
    positive=omega*np.trace(rho@(delta.conj().T@delta)).real
    mixed=2*omega*np.trace(rho@(free.conj().T@delta)).real
    assert positive>=0 and abs(mixed)>1e-7
    energy_error=abs(actual-positive-mixed)
    assert energy_error<2e-8
    return dict(boson_levels=levels,h_final=h_final,duration=duration,
                exact_bath_energy_change=actual,current_square_term=positive,
                initial_noise_current_term=mixed,energy_identity_error=energy_error,
                state_weighted_operator_error=state_error)


def scale_checks():
    rows=[]
    for k in (2000,4000,8000):
        q=k**4
        X=2*np.pi*k**6/(1+np.sqrt(1+12*np.pi*k**6))
        b=np.sqrt(X)*k**5;C=1+6*X
        beta=2*np.pi*b**3/(C*(2*q)*(1/k)*.5)
        B=np.pi/q
        assert beta<2*k**13.5
        W0=B*(1+np.log(beta/4))
        uniform=B*(1-np.log(2)+13.5*np.log(k)+1/32)
        assert W0+B/32<uniform<2.03e-11
        rows.append(dict(k=k,beta_over_lp=beta,frequency_cutoff_times_lp=1,
                         electric_B=B,midpoint_bath_exponent_upper=W0,
                         free_pair_increment_upper=B/32,
                         uniform_family_exponent_upper=uniform))
    return rows


def main():
    radial=[]
    for kind in ('electric','angular'):
        for cells in (60,120,240):
            radial.append(radial_case(3,7,.005,.2,cells,kind))
    refinement=[]
    for start in (0,3):
        cases=radial[start:start+3]
        errors=[max(r['primal_dual_relative_error'] for r in c['comparisons']) for c in cases]
        assert errors[2]<errors[1]<errors[0]
        assert errors[2]<.004
        assert cases[2]['finite_element_contact']<cases[1]['finite_element_contact']<cases[0]['finite_element_contact']
        refinement.append(dict(kind=cases[0]['kind'],primal_dual_max_errors=errors))
    interacting=[interacting_current_check(n,h) for n,h in ((18,.2),(18,.7),(24,.7))]
    interacting_refinement=abs(interacting[1]['exact_bath_energy_change']-interacting[2]['exact_bath_energy_change'])
    assert interacting_refinement<1e-10
    receipt=dict(date='2026-10-13',scope='off-shell projected Gaussian phase action and exact bath current identity; physical many-channel switching stress not computed',
                 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 numpy=np.__version__,scipy=scipy.__version__,radial_regulators=radial,
                 refinement=refinement,quantum_oscillators=oscillator_checks(),
                 interacting_current_checks=interacting,interacting_energy_refinement=interacting_refinement,
                 cumulative_measure_checks=measure_bounds(),single_species_electric_band=scale_checks(),
                 limitations=['Finite RN regulators have reflecting horizon cuts and an explicit vanishing mesh contact term',
                              'The compact winding/sector weights and nonlinear angular completion are not determined',
                              'The Planck bound is an upper cutoff comparison for the retained Gaussian band, not UV matching',
                              'Small vertex overlap loss does not bound switching energy; a saturating one-mode counterexample is proved'],passed=True)
    (ROOT/'evidence/phase-influence-action.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: 6 RN regulators, 3 equilibrium oscillators, 3 interacting current identities, 9 measure bounds, 3 electric scale bounds')


if __name__=='__main__':
    main()
