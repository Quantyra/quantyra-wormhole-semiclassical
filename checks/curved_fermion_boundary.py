"""Curved paired-Dirac boundary heat coefficients, Weyl transfer and contacts."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import sympy as s
from scipy.integrate import quad
import scipy
from scipy.linalg import expm


def zero(expr):
    items=list(expr) if isinstance(expr,s.MatrixBase) else [expr]
    for item in items:
        assert s.simplify(s.expand(item))==0,item


def clifford_checks():
    I=s.eye(2); p1=s.Matrix([[0,1],[1,0]]);p2=s.Matrix([[0,-s.I],[s.I,0]]);p3=s.diag(1,-1)
    k=s.kronecker_product
    normal=k(p3,p1,I)
    tangent=[-k(I,p3,I),-k(I,p2,p1),-k(I,p2,p2)]
    gam=[normal]+tangent
    gen=s.I*k(p2,p1,I)
    eye=s.eye(8)
    for i in range(4):
        zero(gen*gam[i]+gam[i]*gen)
        for j in range(4):zero(gam[i]*gam[j]+gam[j]*gam[i]-2*(i==j)*eye)
    zero(gen.H+gen);zero(gen**2+eye)
    Cs=[normal*g*gen for g in tangent]
    for sign in [1,-1]:
        P=(eye+sign*normal)/2
        for a in range(3):
            zero(Cs[a].H-Cs[a]);zero(Cs[a]*P-P*Cs[a])
            for b in range(3):zero(Cs[a]*Cs[b]+Cs[b]*Cs[a]-2*(a==b)*eye)
        # All minimal-Dirac curvature/gauge endomorphism and connection bases.
        basis=[eye]+[gam[i]*gam[j] for i in range(4) for j in range(i+1,4)]
        for E in basis:
            zero(gen*E-E*gen);zero(s.trace(gen*E))
            for c in Cs:zero(s.trace(c*P*E))
    l0,l1,l2,l01,l02,l12=s.symbols('L00 L11 L22 L01 L02 L12',real=True)
    L=s.Matrix([[l0,l01,l02],[l01,l1,l12],[l02,l12,l2]])
    uv=s.Matrix(s.symbols('u0 u1 u2',real=True)); traceL=s.trace(L);u2=(uv.T*uv)[0]
    Lu=(uv.T*L*uv)[0]
    results=[]
    for sign in [1,-1]:
        # Adjoint domain has reversed base projector and angle.
        angle_sign=sign
        chi=sign*normal;P=(eye+chi)/2
        derivatives=[sign*(-sum((L[a,b]*tangent[b] for b in range(3)),s.zeros(8))+angle_sign*uv[a]*gen*normal) for a in range(3)]
        S0=-traceL*P/2
        S1=angle_sign*sum((uv[a]*Cs[a] for a in range(3)),s.zeros(8))*P/2
        S=S0+S1
        induced=sum((P*normal*tangent[a]*derivatives[a]*P/2 for a in range(3)),s.zeros(8))
        zero(induced-S)
        T=sum((d*d for d in derivatives),s.zeros(8))
        zero(s.trace(T)-8*(s.trace(L*L)+u2))
        zero(s.trace(S1));zero(s.trace(S1*S1)-u2)
        zero(s.trace(T*S)-s.trace(L*L)*s.trace(S0)+2*traceL*u2+4*Lu)
        coeff3=s.trace(192*(S*S-S0*S0)-12*(T-s.trace(L*L)*eye))/384
        zero(coeff3-u2/4)
        # Complete theta-derivative part of the unsmeared mixed a4 formula.
        term=480*(S*S-S0*S0)*traceL+480*(S**3-S0**3)
        term-=12*(T-s.trace(L*L)*eye)*traceL
        term-=24*sum((L[a,b]*(derivatives[a]*derivatives[b]-
                    sum((L[a,c]*tangent[c] for c in range(3)),s.zeros(8))*
                    sum((L[b,c]*tangent[c] for c in range(3)),s.zeros(8)))
                    for a in range(3) for b in range(3)),s.zeros(8))
        term-=120*(T*S-s.trace(L*L)*S0)
        coeff4=s.expand(s.trace(term)/360)
        expected=s.Rational(4,5)*(Lu-traceL*u2/3)
        zero(coeff4-expected)
        normal_smear=s.trace(240*(S*S-S0*S0)-18*(T-s.trace(L*L)*eye))/360
        zero(normal_smear-s.Rational(4,15)*u2)
        # The full cross term is necessary: dropping the L.u.C contribution
        # in chi_:a chi_:a S gives a nonzero spurious umbilic coefficient.
        missing_cross=coeff4-s.Rational(4,3)*Lu
        zero(missing_cross.subs({l0:1,l1:1,l2:1,l01:0,l02:0,l12:0})+s.Rational(4,3)*u2)
        results.append(dict(domain='Dirac' if sign==1 else 'adjoint',a3=str(s.factor(coeff3)),a4=str(s.factor(coeff4)),a4_normal_smearing=str(s.factor(normal_smear))))
    return dict(Clifford_and_minimal_curvature_gauge_traces='PASS',arbitrary_symmetric_second_fundamental_form='PASS',induced_Robin_both_domains='PASS',rank=8,coefficients=results,umbilic_omission_control='PASS')


def weyl_and_asymptotic_checks():
    sig,ns,L,dsq,hat=s.symbols('sigma normal_sigma L theta_squared Lhat_theta_theta',real=True)
    # Orthonormal tensor contractions have weights -3 (hat) and -2 (dsq).
    zero(s.exp(3*sig)*s.exp(-3*sig)*hat-hat)
    zero(s.exp(3*sig)*s.exp(-sig)*s.exp(-2*sig)*dsq-dsq)
    transformedJ=s.exp(3*sig)*s.exp(-sig)*(L-3*ns)*s.exp(-2*sig)*dsq
    zero(transformedJ-L*dsq+3*ns*dsq)
    sigma2,n2=s.symbols('sigma2 normal_sigma2',real=True)
    b=s.Rational(4,5);d=s.Rational(4,15)
    W=b*sig*hat+d*ns*dsq
    W2=s.exp(3*sig)*(b*sigma2*s.exp(-3*sig)*hat+d*s.exp(-sig)*n2*s.exp(-2*sig)*dsq)
    zero(W+W2-(b*(sig+sigma2)*hat+d*(ns+n2)*dsq))
    R,ell,m,Q2=s.symbols('R ell m Q2',positive=True)
    f=1-2*m/R+Q2/R**2+R**2/ell**2
    B2=R**2/f;w=s.diff(f,R)/2-f/R
    zero(w-(-1/R+3*m/R**2-2*Q2/R**3))
    typeB=s.Rational(2,3)*B2*w
    normalW=-B2*s.diff(f,R)/2
    # Relative reference uses the same asymptotic clock and round sphere.
    base={m:0,Q2:0}
    z=s.symbols('z',positive=True)
    expansionB=s.series((typeB-typeB.subs(base)).subs(R,1/z),z,0,4).removeO()
    expansionN=s.series((normalW-normalW.subs(base)).subs(R,1/z),z,0,4).removeO()
    zero(expansionB-(2*m*ell**2*z**2-s.Rational(4,3)*Q2*ell**2*z**3))
    zero(expansionN-(-3*m*ell**2*z**2+2*Q2*ell**2*z**3))
    zero(s.limit(normalW+R,R,s.oo))
    zero(s.limit(typeB,R,s.oo))
    return dict(conformal_weights_and_WZ_cocycle='PASS',finite_L_source_term_shift='delta(J)=-3 integral n(omega)*(Dtheta)^2',
                RN_typeB_relative_series=str(expansionB),RN_normal_relative_series=str(expansionN),
                RN_WZ_relative_limit='0; typeB is O(log(R)/R^2), normal term O(1/R^2), at fixed parameters')


def contact_variations():
    t=s.symbols('t',real=True);zeta=s.Rational(7,9)
    # Lorentz boundary -N(t)^2 dt^2+R(t)^2 dOmega^2, prescribed theta(t).
    N=s.exp(s.sin(t)/9);R=2+s.cos(t)/7;theta=s.sin(t)+s.sin(2*t)/11
    derivatives=[s.sin(2*t),s.cos(3*t)/5,s.cos(t)+s.sin(3*t)/7]
    # Periodic integration removes temporal endpoint terms without fixing
    # lapse/radius variations; angular directions are retained by 4*pi R^2.
    lag=lambda nn,rr,th:2*s.pi*zeta*rr**2*s.diff(th,t)**2/nn
    euler_theta=-4*s.pi*zeta*s.diff(R**2*s.diff(theta,t)/N,t)
    dN=-2*s.pi*zeta*R**2*s.diff(theta,t)**2/N**2
    dR=4*s.pi*zeta*R*s.diff(theta,t)**2/N
    arrays=[[1,0,0],[0,1,0],[0,0,1],[1,1,1]]
    cases=[]
    for signs in [1,-1]:
        for mask in arrays:
            var=[signs*a*b for a,b in zip(mask,derivatives)]
            expected=s.lambdify(t,dN*var[0]+dR*var[1]+euler_theta*var[2],'numpy',cse=True)
            e=s.symbols('eps',real=True)
            direct=s.lambdify((t,e),lag(N+e*var[0],R+e*var[1],theta+e*var[2]),'numpy',cse=True)
            target=quad(expected,0,2*np.pi,epsabs=1e-11,epsrel=1e-11)[0]
            estimates=[]
            for step in [1e-3,5e-4]:
                def action(eps):return quad(lambda v:direct(v,eps),0,2*np.pi,epsabs=1e-11,epsrel=1e-11)[0]
                estimates.append((-action(2*step)+8*action(step)-8*action(-step)+action(-2*step))/(12*step))
            err=abs(estimates[-1]-target);ref=abs(estimates[-1]-estimates[0])
            assert err<2e-7 and ref<2e-7
            cases.append(dict(mask=mask,sign=signs,discrepancy=err,step_refinement=ref))
    # Explicit homogeneous Ward identity with a changing boundary radius:
    # rho=p=zeta*(theta_dot/N)^2/2, O=-zeta/(N R^2)d_t(R^2 theta_dot/N).
    rho=zeta*s.diff(theta,t)**2/(2*N**2)
    source=-zeta*s.diff(R**2*s.diff(theta,t)/N,t)/(N*R**2)
    zero(s.diff(rho,t)/N+4*s.diff(R,t)*rho/(N*R)+source*s.diff(theta,t)/N)
    return cases


def direct_surface_cases():
    eye2=np.eye(2);p1=np.array([[0,1],[1,0]],complex);p2=np.array([[0,-1j],[1j,0]],complex);p3=np.diag([1.,-1.])
    k=lambda a,b,c:np.kron(np.kron(a,b),c)
    gn=k(p3,p1,eye2);ga=[-k(eye2,p3,eye2),-k(eye2,p2,p1),-k(eye2,p2,p2)]
    gen=1j*k(p2,p1,eye2);identity=np.eye(8);rng=np.random.default_rng(74261013)
    cases=[]
    for th in [-.6,0.,4*np.pi/2000,.2,.7,1.]:
        for kind in ['umbilic','anisotropic']:
            raw=rng.normal(size=(3,3))/10
            L=np.eye(3)/7 if kind=='umbilic' else (raw+raw.T)/2
            u=rng.normal(size=3)/5
            # Boundary graph r=y^T L y/2 in flat ambient space. Differentiate
            # its unit normal and rotating projector directly in coordinates.
            def chi(y,gradient):
                slope=L@y
                gamma_n=(gn-sum((slope[a]*ga[a] for a in range(3)),np.zeros((8,8),complex)))/np.sqrt(1+slope@slope)
                V=expm((th+gradient@y)*gen/2)
                return V@gamma_n@V.conj().T
            base_chi=chi(np.zeros(3),u);P=(identity+base_chi)/2
            Cs=[gn@g@gen for g in ga]
            targetS=(-np.trace(L)*identity+sum((u[a]*Cs[a] for a in range(3)),np.zeros((8,8),complex)))@P/2
            targets=np.array([u@u/4,4*(u@L@u-np.trace(L)*(u@u)/3)/5,4*(u@u)/15])
            values=[];s_errors=[]
            for step in [1e-3,5e-4]:
                def coefficient(gradient):
                    derivatives=[]
                    for a in range(3):
                        e=np.eye(3)[a]*step
                        derivatives.append((-chi(2*e,gradient)+8*chi(e,gradient)-8*chi(-e,gradient)+chi(-2*e,gradient))/(12*step))
                    S=sum((P@gn@ga[a]@derivatives[a]@P/2 for a in range(3)),np.zeros((8,8),complex))
                    T=sum((d@d for d in derivatives),np.zeros((8,8),complex))
                    tensor=sum((L[a,b]*derivatives[a]@derivatives[b] for a in range(3) for b in range(3)),np.zeros((8,8),complex))
                    a3=np.trace(192*S@S-12*T)/384
                    a4=np.trace(480*np.trace(L)*S@S+480*S@S@S-12*np.trace(L)*T-24*tensor-120*T@S)/360
                    smear=np.trace(240*S@S-18*T)/360
                    return np.array([a3,a4,smear]),S
                active,SS=coefficient(u);static,_=coefficient(np.zeros(3))
                values.append(active-static);s_errors.append(np.linalg.norm(SS-targetS))
            error=np.max(np.abs(values[-1]-targets));refinement=np.max(np.abs(values[-1]-values[0]))
            assert error<2e-8 and refinement<2e-8 and s_errors[-1]<2e-8,(th,kind,error,refinement,s_errors[-1],values[-1],targets)
            cases.append(dict(theta=th,shape=kind,heat_discrepancy=float(error),step_refinement=float(refinement),Robin_discrepancy=float(s_errors[-1]),unsmeared_derivative_a4=float(targets[1])))
    return cases


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=dict(status='PASS',clifford=clifford_checks(),weyl=weyl_and_asymptotic_checks(),contact_variations=contact_variations(),direct_surfaces=direct_surface_cases(),
                versions=dict(python=platform.python_version(),sympy=s.__version__,numpy=np.__version__,scipy=scipy.__version__),
                scope='Matched identical Euclidean collars for a free minimal Dirac field; local derivative heat/Weyl sector and prescribed-source contact distinction. No full state, determinant phase, controller or opening proof.',
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))
