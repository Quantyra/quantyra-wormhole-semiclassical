"""Joint small-port/large-flux checks for the matched exterior trial.

Tests overlap and bulk residuals, not a physical loaded PDE solution.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import cumulative_simpson, simpson, solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_hermite, roots_legendre
import sympy as sp


class Dipole:
    def __init__(self):
        r,z,u=sp.symbols('r z u',real=True)
        rm=sp.sqrt(r*r+(z+sp.Rational(1,2))**2)
        rp=sp.sqrt(r*r+(z-sp.Rational(1,2))**2)
        flux=((z+sp.Rational(1,2))/rm-(z-sp.Rational(1,2))/rp)/2
        mass=(flux-u)/r
        B=sp.Matrix([r/rm**3-r/rp**3,(z+sp.Rational(1,2))/rm**3-(z-sp.Rational(1,2))/rp**3])/2
        self.bf=sp.lambdify((r,z),B,'numpy',cse=True)
        self.jf=sp.lambdify((r,z),B.jacobian([r,z]),'numpy',cse=True)
        self.hf=sp.lambdify((r,z,u),sp.hessian(mass,[r,z]),'numpy',cse=True)
        self.mf=sp.lambdify((r,z,u),mass,'numpy',cse=True)
        self.uf=sp.lambdify((r,z),flux,'numpy',cse=True)

    def geometry(self,points,U):
        r,z=np.asarray(points).reshape(2,-1)
        B=np.asarray(self.bf(r,z),float).reshape(2,-1)
        J=np.asarray(self.jf(r,z),float).reshape(2,2,-1)
        H=np.asarray(self.hf(r,z,U),float).reshape(2,2,-1)
        b=np.linalg.norm(B,axis=0);t=B/b;n=np.array([-t[1],t[0]])
        bs=np.einsum('in,ijn,jn->n',t,J,t)
        k=np.einsum('in,ijn,jn->n',n,J,t)/b
        c=np.einsum('in,ijn,jn->n',n,H,n)
        return dict(b=b,t=t,n=n,omega=bs/(4*b),k=k,c=c)

    def on_sphere(self,rho,U):
        # Select the source-to-midplane branch; a large sphere can also
        # intersect the return branch and enclose the opposite monopole.
        lower=max(.0001,np.arccos(min(1.,.5/rho)))
        theta=brentq(lambda th:float(self.uf(rho*np.sin(th),-.5+rho*np.cos(th)))-U,
                     lower,np.pi-.0001,xtol=5e-15)
        return theta,np.array([rho*np.sin(theta),-.5+rho*np.cos(theta)])

    def section(self,rho,U,Q,Xi,step_fraction=1e-4,normal_order=256):
        theta,p=self.on_sphere(rho,U);g=self.geometry(p,U)
        t=g['t'][:,0];n=g['n'][:,0];step=step_fraction*rho
        shifted=[self.geometry(p+j*step*t,U) for j in [-2,-1,1,2]]
        for name in ['omega','k','c']:
            g[name+'s']=(shifted[0][name]-8*shifted[1][name]+8*shifted[2][name]-shifted[3][name])/(12*step)
        b=float(g['b'][0]);c=float(g['c'][0]);a=float(g['omega'][0]);k=float(g['k'][0])
        cs=float(g['cs'][0]);ks=float(g['ks'][0]);aas=float(g['omegas'][0])
        def width(point,curvature):
            R,Z=point
            return .15/np.sqrt(1/R**2+1/(R*R+(Z+.5)**2)+1/(R*R+(Z-.5)**2)+1+curvature**2)
        tube=width(p,k)
        widths=[width(p+j*step*t,float(gg['k'][0])) for j,gg in zip([-2,-1,1,2],shifted)]
        tube_s=(widths[0]-8*widths[1]+8*widths[2]-widths[3])/(12*step)
        y,w=roots_hermite(normal_order);w/=np.sqrt(np.pi);N=y/np.sqrt(Q*b)
        inside=abs(N)<tube;N=N[inside];w=w[inside]
        h=1-k*N;hs=-ks*N
        assert min(h)>.8
        mm=self.mf(*(p[:,None]+n[:,None]*N),U)
        P=1-Q*c*N**3/6;Z=a*N;G=a*(1-2*Q*b*N*N)
        delta=mm-b*N;rem=delta-c*N*N/2
        upper=G*(-k*N+Q*c*N**3/6)/h+Q*cs*N**3/(6*h)+Xi*(1/h-1)*P+hs*P/(2*h*h)-Q*a*N*delta
        lower=Q*rem-Q*Q*c*N**3*delta/6+aas*N/h+Z*G/h-Xi*Z*(1/h+1)-hs*Z/(2*h*h)
        # Direct Pauli-operator evaluation, independent of the expanded residual.
        sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.])
        vec=np.array([P,Z])
        dn=np.array([-Q*c*N*N/2-Q*b*N*P,a-Q*b*N*Z])
        ds=np.array([-Q*cs*N**3/6+P*(G-Xi),aas*N+Z*(G-Xi)])
        direct=1j*sx@dn-1j*sz@ds/h+1j*sz@vec*hs/(2*h*h)+sy@vec*(Q*mm)-1j*Xi*vec
        assert np.max(abs(direct-1j*np.array([upper,lower])))<1e-7
        chi=np.ones_like(N);chi_y=np.zeros_like(N);chi_s=np.zeros_like(N)
        transition=abs(N)>tube/2;angle=np.pi*(abs(N[transition])/tube-.5)
        chi[transition]=np.cos(angle)**2
        chi_y[transition]=-np.pi/tube*np.sin(2*angle)*np.sign(N[transition])
        chi_s[transition]=np.pi*abs(N[transition])*tube_s/tube**2*np.sin(2*angle)
        upper,lower=chi*upper+chi_y*Z-chi_s*P/h,chi*lower+chi_y*P+chi_s*Z/h
        norm=np.sqrt(np.sum(w*(upper*upper+lower*lower)))
        radial_cos=float(t@np.array([np.sin(theta),np.cos(theta)]))
        assert radial_cos>0
        return float(norm),radial_cos

    def overlap(self,rho,U,Q,Xi,order):
        theta0,p=self.on_sphere(rho,U);gg=self.geometry(p,U)
        t0=gg['t'][:,0]
        cosine=float(t0@np.array([np.sin(theta0),np.cos(theta0)]))
        A=rho*rho*float(gg['b'][0])*cosine
        T=np.linspace(-9,9,order)
        theta=theta0+T/np.sqrt(Q*A)
        assert min(theta)>0 and max(theta)<np.pi
        def rhs(_,point):
            B=np.asarray(self.bf(*point),float).ravel()
            return B/np.linalg.norm(B)
        opts=dict(rtol=3e-13,atol=1e-15*rho,dense_output=True,max_step=rho/40)
        fw=solve_ivp(rhs,[0,.6*rho],p,**opts)
        bw=solve_ivp(rhs,[0,-.6*rho],p,**opts)
        assert fw.success and bw.success
        def curve(s):
            answer=np.empty((2,len(s)));pos=s>=0
            if np.any(pos):answer[:,pos]=fw.sol(s[pos])
            if np.any(~pos):answer[:,~pos]=bw.sol(s[~pos])
            return answer
        points=np.array([rho*np.sin(theta),-.5+rho*np.cos(theta)])
        ss=np.einsum('i,in->n',t0,points-p[:,None])
        for _ in range(8):
            center=curve(ss);g=self.geometry(center,U)
            diff=points-center;yy=np.sum(diff*g['n'],axis=0)
            ss+=np.sum(diff*g['t'],axis=0)/(1-g['k']*yy)
        center=curve(ss);g=self.geometry(center,U)
        yy=np.sum((points-center)*g['n'],axis=0)
        assert np.max(np.linalg.norm(points-center-g['n']*yy,axis=0))/rho<1e-10
        metric=1-g['k']*yy
        gauss=(Q*g['b']/np.pi)**.25*np.exp(-Q*g['b']*yy*yy/2)
        P=1-Q*g['c']*yy**3/6;V=g['omega']*yy
        angle=np.arctan2(g['t'][0],g['t'][1])-theta
        scalar=np.sqrt(rho/metric)*gauss*np.exp(-Xi*ss)
        tube_in=scalar*(np.cos(angle/2)*P-np.sin(angle/2)*V)
        tube_out=scalar*(np.sin(angle/2)*P+np.cos(angle/2)*V)

        # Independent exact angular zero mode and its integral inverse.
        th=np.asarray(theta,np.longdouble);rr=np.longdouble(rho)
        S=np.sqrt(1+rr*rr-2*rr*np.cos(th));ell=(1-U)*Q
        W=Q/2*np.log((S+1-rr)/((1-rr)*(S+1+rr)))
        logg=ell*np.log(np.sin(th/2))+(Q-ell)*np.log(np.cos(th/2))-W
        ge=np.exp(np.asarray(logg-max(logg),float));ge/=np.sqrt(simpson(ge*ge,x=theta))
        Se=(rr-np.cos(th))/S
        We=np.asarray(Q/2*((Se-1)/(S+1-rr)+1/(1-rr)-(Se+1)/(S+1+rr)),float)
        centered=We-simpson(We*ge*ge,x=theta)
        integrand=centered*ge*ge
        left=cumulative_simpson(integrand,x=theta,initial=0)
        right=-cumulative_simpson(integrand[::-1],x=-theta[::-1],initial=0)[::-1]
        cumulative=np.where(T<=0,left,right)
        he=rho*cumulative/ge
        # The unitary radial operator gives lower component h with this sign.
        norm=lambda v:float(np.sqrt(simpson(v*v,x=theta)))
        mismatch=np.sqrt(norm(tube_in-ge)**2+norm(tube_out-he)**2)
        reflected_in=scalar*np.exp(2*Xi*ss)*(np.cos(angle/2)*P-np.sin(angle/2)*V)
        reflected_out=-scalar*np.exp(2*Xi*ss)*(np.sin(angle/2)*P+np.cos(angle/2)*V)
        receiver_mismatch=np.sqrt(norm(reflected_in-ge)**2+norm(reflected_out+he)**2)
        return dict(mismatch=float(mismatch),receiver_mismatch=float(receiver_mismatch),
                    scaled=float(max(mismatch,receiver_mismatch)/(Q**(-.5)+rho*rho)),
                    corrector_norm=norm(he),corrector_upper=27/np.sqrt(160)*rho*rho,
                    exact_mode_norm=norm(ge))


def verify():
    # Exponents are exact rational consequences of the complete error bound.
    outer=F(1,4);port=F(3,4)
    exponents={'interior_squared':-1+outer,'collar_squared':-3*outer,
               'collar_derivative_squared':1-7*outer,
               'overlap_squared':-1+outer,
               'overlap_tilt_squared':-3*outer,
               'basis_weighted_squared':F(3,2)-3*port,
               'corrector_weighted_squared':F(1,2)-3*port}
    assert all(v<=F(-3,4) for v in exponents.values())
    assert port>outer
    model= Dipole();overlaps=[];bulk=[]
    for U in [.25,.5,.75]:
        qs=[10001,160001,2560001] if U==.5 else [10002,160002,2560002]
        for Xi in [.05,4.]:
            scaled_rows=[]
            for Q in qs:
                assert abs((1-U)*Q-(np.floor((1-U)*Q)+.5))<1e-10
                a=.4*Q**(-.25);p=Q**(-.75)
                for factor in [1.,1.5,2.]:
                    rho=factor*a
                    co=model.overlap(rho,U,Q,Xi,1025)
                    fi=model.overlap(rho,U,Q,Xi,2049)
                    refinement=max(abs(co[k]-fi[k]) for k in co)
                    assert refinement<2e-5
                    assert fi['corrector_norm']<=fi['corrector_upper']+2e-8
                    assert fi['scaled']<5
                    overlaps.append(dict(q=Q,flux_label=U,xi_d=Xi,port=p,
                        overlap_radius=rho,refinement=refinement,**fi))
                integrals=[]
                # Integrate the exact normal-section residual over the half line,
                # then reflect. Omit attenuation to obtain a larger bulk norm.
                for order,step,normal_order in [(64,1e-4,256),(96,5e-5,512)]:
                    z,w=roots_legendre(order)
                    lo=np.log(a);hi=np.log(1/(2*U))
                    rr=np.exp((lo+hi)/2+(hi-lo)/2*z)
                    vals=[model.section(r,U,Q,Xi,step,normal_order) for r in rr]
                    integral=2*(hi-lo)/2*np.sum(w*rr*np.array([v[0]**2/v[1] for v in vals]))
                    integrals.append(float(integral))
                assert abs(integrals[0]-integrals[1])<.002*integrals[1]+1e-9,(U,Xi,Q,integrals)
                scaled=np.sqrt(integrals[1]*Q*a)
                scaled_rows.append(float(scaled))
                bulk.append(dict(q=Q,flux_label=U,xi_d=Xi,overlap_scale=a,
                    bulk_residual_squared=integrals[1],scaled_residual=float(scaled),
                    integration_refinement=abs(integrals[0]-integrals[1])))
            # The q-dependent cutoff term can dominate at the lowest flux;
            # retain that case and test the asymptotic upper order later.
            assert scaled_rows[-1]<=1.15*scaled_rows[-2],(U,Xi,scaled_rows)
    budgets=[]
    for Q in [10**4,10**6,10**8,10**10,10**12]:
        a=.4*Q**(-.25);p=Q**(-.75);xi=.05
        beta=1/p;ell=p/4;gap=np.sqrt(Q+1)
        d0=-np.expm1(-2*(xi+gap*beta)*ell)
        M=Q*p*p/(4*(1-p)**2)
        basis=M*M/d0*(1/(Q+1)+beta/(xi*gap))
        hbound=27/np.sqrt(160)*p*p
        Lhbound=3/np.sqrt(2)*np.sqrt(Q)*p*p+M*hbound
        corrector=(hbound*hbound+beta/xi*hbound*Lhbound)/d0
        F2=1/(Q*a)+a**3+Q*a**7
        budgets.append(dict(q=Q,port=p,overlap=a,basis_cost=float(basis),
                            corrector_cost=float(corrector),
                            bulk_order_expression=F2,
                            basis_scaled=float(basis*Q**.75),
                            bulk_scaled=float(F2*Q**.75)))
    assert max(b['basis_scaled'] for b in budgets)<2
    return dict(status='passed_joint_matching_checks',exact_exponents={k:str(v) for k,v in exponents.items()},
        overlap_cases=overlaps,bulk_cases=bulk,analytic_order_budget_cases=budgets,
        ordinary_matched_trial_constructed=True,conditional_trimmed_return_order='q^(-3/8)',
        finite_flux_threshold_certified=False,loaded_exterior_existence_proved=False,
        physical_metric_and_gauge_matching_proved=False,complete_energy_derived=False,
        novelty_established=False,lean_invocations=0,numerical_checks_are_supplementary=True,
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
    result=verify();result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:len(result[k]) for k in ['exact_exponents','overlap_cases','bulk_cases','analytic_order_budget_cases']}))
