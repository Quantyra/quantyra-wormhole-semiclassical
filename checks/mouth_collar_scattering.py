"""Prepared-input scattering in a finite flat dipole annulus, supplementary checks."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_bvp
from scipy.special import eval_jacobi, roots_jacobi, betaln
import sympy as sp


class Angular:
    def __init__(self,q,m,n,quadrature=64):
        self.q=q;self.m=m;self.n=n
        z,w=roots_jacobi(quadrature,q-m-1,m)
        self.x=(z+1)/2;self.weights=w/w.sum()
        P=np.array([eval_jacobi(k,m,q-m-1,1-2*self.x) for k in range(n+1)])
        norms=np.sqrt(np.sum(P*P*self.weights,axis=1))
        self.upper=P/norms[:,None]
        lam=np.sqrt(np.arange(1,n+1)*(np.arange(1,n+1)+q))
        self.lower=np.array([(k+q)/lam[k-1]/norms[k]*np.sqrt(self.x*(1-self.x))*eval_jacobi(k-1,m+1,q-m,1-2*self.x) for k in range(1,n+1)])
        self.K0=np.zeros((n+1,n));self.K0[np.arange(1,n+1),np.arange(n)]=lam
        assert np.max(abs((self.upper*self.weights)@self.upper.T-np.eye(n+1)))<2e-12
        assert np.max(abs((self.lower*self.weights)@self.lower.T-np.eye(n)))<2e-12
        self.gauss_z,self.gauss_w=np.polynomial.legendre.leggauss(64)

    def fields(self,x,e):
        c=1-2*x;S=np.sqrt(1+e*e-2*e*c);D=S*(S+1-e*c)
        W=self.q*np.log((S+1-e)/((1-e)*(S+1+e)))/2
        # Differentiate the closed formula independently of the angular inverse.
        Se=(e-c)/S
        We=self.q/2*((Se-1)/(S+1-e)+1/(1-e)-(Se+1)/(S+1+e))
        return W,We,D

    def profiles(self,e):
        x=self.x;W,We,_=self.fields(x,e)
        Z=np.sum(self.weights*np.exp(-2*W));g=np.exp(-W)/np.sqrt(Z) # g/g0
        mu=np.sum(self.weights*g*g*We)
        gp=-(We-mu)*g # d=1
        # Direct cumulative-integral inverse, integrating from the nearer pole.
        pivot=(self.m+1)/(self.q+1)
        left=x<pivot;lo=np.where(left,0.,x);hi=np.where(left,x,1.)
        xx=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]/2*self.gauss_z
        wi,wei,_=self.fields(xx,e)
        logp=self.m*np.log(xx)+(self.q-self.m-1)*np.log1p(-xx)-betaln(self.m+1,self.q-self.m)-2*wi-np.log(Z)
        I=(hi-lo)/2*np.sum(self.gauss_w*(wei-mu)*np.exp(logp),axis=1)*np.where(left,1.,-1.)
        p=np.exp(self.m*np.log(x)+(self.q-self.m-1)*np.log1p(-x)-betaln(self.m+1,self.q-self.m))*g*g
        ratio=e*I/(p*np.sqrt(x*(1-x)))
        h=ratio*g # h/g0
        return g,gp,h

    def K(self,e):
        E=np.atleast_1d(e);x=self.x[:,None]
        _,_,D=self.fields(x,E[None,:])
        ds=-self.q*E[None,:]**2*np.sqrt(x*(1-x))/D
        return self.K0[:,:,None]+np.einsum('ia,ja,ak,a->ijk',self.upper,self.lower,ds,self.weights)

    def project(self,g,h):
        return self.upper@(self.weights*g),self.lower@(self.weights*h)

    def solve(self,e0,e1,Xi,tol=2e-9):
        g0,_,h0=self.profiles(e0);g1,_,h1=self.profiles(e1)
        u0,v0=self.project(g0,h0);u1,v1=self.project(g1,h1)
        N=self.n+1
        def fun(r,y):
            K=self.K(r)
            return np.vstack([-Xi*y[:N]-np.einsum('ijk,jk->ik',K,y[N:])/r,Xi*y[N:]-np.einsum('ijk,ik->jk',K,y[:N])/r])
        def bc(ya,yb):return np.r_[ya[:N]-u0,yb[N:]]
        mesh=np.linspace(e0,e1,40)
        guess=np.zeros((2*self.n+1,len(mesh)));guess[:N]=u0[:,None]*np.exp(-Xi*(mesh-e0))
        sol=solve_bvp(fun,bc,mesh,guess,tol=tol,max_nodes=1600)
        assert sol.success,sol.message
        outgoing=np.r_[sol.y[N:,0],sol.y[:N,-1]]
        trial=np.r_[v0,np.exp(-Xi*(e1-e0))*u1]
        error=np.linalg.norm(outgoing-trial)
        # Independent bulk current balance and a posteriori trial residual bound.
        bulk=quad(lambda r:np.sum(sol.sol(r)**2),e0,e1,epsabs=1e-12,epsrel=1e-10)[0]
        current_defect=np.dot(outgoing,outgoing)+2*Xi*bulk-np.dot(u0,u0)
        assert abs(current_defect)<2e-8
        assert np.dot(outgoing,outgoing)<=np.dot(u0,u0)+2e-10
        return outgoing,dict(q=self.q,m=self.m,angular_order=self.n,epsilon_inner=e0,epsilon_outer=e1,xi_d=Xi,error_to_corrected_trial=float(error),current_balance_residual=float(current_defect),mesh_nodes=len(sol.x),max_collocation_residual=float(np.max(sol.rms_residuals)),input_projection_norm=float(np.linalg.norm(u0)))


def analytic_error(q,e0,e1,Xi):
    integrand=lambda e:np.exp(-2*Xi*(e-e0))*(27/5*e+12*np.sqrt(q)*e**3+22/5*Xi*e**2)**2
    return np.sqrt(729/160*e1**4*np.exp(-2*Xi*(e1-e0))+quad(integrand,e0,e1,epsabs=1e-18,epsrel=1e-11)[0]/(2*Xi))


def verify():
    x,e,q=sp.symbols('x e q',positive=True);c=1-2*x;S=sp.sqrt(1+e*e-2*e*c)
    delta=q*(1-(1-e*c)/S)/2
    exact={
        'remote_potential_derivative':sp.simplify(sp.diff(delta,e)-2*q*e*x*(1-x)/S**3)==0,
        'distance_second_derivative':sp.simplify(sp.diff(S,e,2)-(1-c*c)/S**3)==0,
    }
    a,v,b=sp.symbols('a v b',real=True)
    exact['centered_fourth_moment_cancellation']=sp.expand((a*a-2*v)**2-a**4+4*v*(a*a-v))==0
    # Current identity for arbitrary complex finite-dimensional fields/source.
    rng=np.random.default_rng(87231);current=[]
    for N in [1,2,7]:
        for Xi in [.05,1.,5.]:
            C=rng.normal(size=(N,N))+1j*rng.normal(size=(N,N));A=np.block([[np.zeros((N,N)),C],[C.conj().T,np.zeros((N,N))]])
            J=np.diag(np.r_[np.ones(N),-np.ones(N)])
            psi=rng.normal(size=2*N)+1j*rng.normal(size=2*N);f=rng.normal(size=2*N)+1j*rng.normal(size=2*N)
            derivative=1j*J@f-Xi*J@psi-1j*J@A@psi
            lhs=2*np.real(np.vdot(psi,J@derivative));rhs=-2*Xi*np.vdot(psi,psi).real-2*np.imag(np.vdot(psi,f))
            assert abs(lhs-rhs)<1e-11
            current.append(dict(dimension=2*N,xi_d=Xi,residual=float(abs(lhs-rhs))))
    D=F(81,50);D1=F(77,15);D2=F(2420,243);emax=F(1,10)
    mixed=2/D+4*emax*D1/D**2+emax**2*D2/D**2+2*emax**2*D1**2/D**3
    rational=dict(mixed_second_derivative_bound=mixed<F(9,4),first_corrector_derivative_constant=F(135,8)**2<F(27,5)**2*10,second_corrector_derivative_constant=F(25,12)**2<2*(12-F(81,8))**2,corrector_amplitude_constant=F(729,160)<F(11,5)**2)
    E=F(1,12500);aa=F(27,5);bb=F(22);cc=F(12000000)
    poly=aa*aa*E**3/3+2*aa*bb*E**4/4+(bb*bb+2*aa*cc)*E**5/5+2*bb*cc*E**6/6+cc*cc*E**7/7
    uniform_squared=F(729,160)*E**4+10*poly
    rational['uniform_collar_error_less_than_8e_minus6']=uniform_squared<F(1,125000)**2
    assert all(exact.values()) and all(rational.values())

    derivatives=[]
    for Q in [3,11,31]:
        for M in [0,Q//2,Q-1]:
            angular=Angular(Q,M,8)
            for E in [.01,.05,.1]:
                g,gp,h=angular.profiles(E)
                hnorm=np.sqrt(np.sum(angular.weights*h*h))
                estimates=[]
                for step in [2e-5,1e-5]:
                    hp=angular.profiles(E+step)[2];hm=angular.profiles(E-step)[2]
                    estimates.append(float(np.sqrt(np.sum(angular.weights*((hp-hm)/(2*step))**2))))
                bound=27/5*E+12*np.sqrt(Q)*E**3
                assert abs(estimates[1]-estimates[0])<3e-8
                assert estimates[-1]<bound and hnorm<27/np.sqrt(160)*E*E
                u,hh=angular.project(gp,h)
                inverse_residual=np.linalg.norm(angular.K(E)[:,:,0]@hh+E*u)
                assert inverse_residual<2e-9
                derivatives.append(dict(q=Q,m=M,epsilon=E,corrector_norm=float(hnorm),derivative_norms=estimates,derivative_steps=[2e-5,1e-5],proved_derivative_bound=float(bound),angular_inverse_projection_residual=float(inverse_residual)))

    scattering=[]
    for Q in [3,11]:
        for M in [0,Q//2,Q-1]:
            for Xi in [.05,4.]:
                rows=[];outputs=[]
                for N in [4,7,10]:
                    angular=Angular(Q,M,N)
                    out,row=angular.solve(.01,.05,Xi)
                    bound=analytic_error(Q,.01,.05,Xi)
                    assert row['error_to_corrected_trial']<bound
                    row['proved_error_bound']=float(bound)
                    rows.append(row);outputs.append(out)
                # Compare separately the reflected and transmitted coefficients.
                differences=[]
                for j in [0,1]:
                    n=rows[j]['angular_order'];nn=rows[j+1]['angular_order']
                    a=outputs[j];b=outputs[j+1]
                    diff=np.linalg.norm(np.r_[a[:n]-b[:n],a[n:]-b[nn:nn+n+1]])
                    tail=np.linalg.norm(np.r_[b[n:nn],b[nn+n+1:]])
                    differences.append(float(np.hypot(diff,tail)))
                assert differences[-1]<2e-8
                scattering.append(dict(q=Q,m=M,xi_d=Xi,angular_refinements=rows,successive_output_differences=differences))
    return dict(status='passed_finite_annulus_prepared_input_scattering_checks',exact_identities=exact,rational_inequalities=rational,mixed_second_derivative_rational=str(mixed),uniform_error_squared_rational=str(uniform_squared),uniform_error_parameters=dict(q_max=1000000000000,epsilon_outer_max=1/12500,xi_d_min=.05,xi_d_max=5.,error_bound=8e-6),current_identity_cases=current,corrector_derivative_cases=derivatives,scattering_cases=scattering,finite_annulus_prepared_input_bound=True,full_exterior_return_derived=False,full_exterior_return_error_bounded=False,curved_throat_matching_derived=False,novelty_established=False,lean_certified=False,lean_invocations=0,numerical_checks_are_supplementary=True,versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);args=parser.parse_args()
    result=verify();result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:len(result[k]) for k in ['exact_identities','rational_inequalities','current_identity_cases','corrector_derivative_cases','scattering_cases']}))
