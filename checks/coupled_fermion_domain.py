"""Coupled four-dimensional Dirac domain and driven product-collar heat data."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
import sympy as s
from scipy.linalg import expm
from scipy.optimize import brentq

I2=s.eye(2)
p1=s.Matrix([[0,1],[1,0]])
p2=s.Matrix([[0,-s.I],[s.I,0]])
p3=s.diag(1,-1)
kron=s.kronecker_product
I4=s.eye(4)
A=kron(p3,p2)
B=-A
G=s.I*kron(p2,p1)
C=s.I*A*G
chi0=kron(p3,p1)


def zero(expr):
    values=list(expr) if isinstance(expr,s.MatrixBase) else [expr]
    for v in values:
        assert s.simplify(s.trigsimp(s.expand_trig(v)))==0,v


def analytic():
    th,dth=s.symbols('theta theta_prime',real=True)
    chi=s.cos(th)*chi0-s.sin(th)*kron(p1,I2)
    Pi=(I4+chi)/2
    P0=(I4+chi0)/2
    V=s.cos(th/2)*I4+s.sin(th/2)*G
    dotPi=s.diff(Pi,th)*dth
    S=C*Pi*dth/2
    for expr in [G.H+G,G*G+I4,C.H-C,C*C-I4,chi.H-chi,
                 chi*chi-I4,Pi*Pi-Pi,chi*B+B*chi,Pi*B*Pi,
                 V.H*V-I4,V*P0*V.H-Pi,A*G-G*A,
                 G*kron(I2,p1)-kron(I2,p1)*G,
                 S-s.I*A*dotPi*Pi,S.H-S,Pi*S-S,S*Pi-S,
                 -s.I*P0*V.H*B*s.diff(V,th)*P0-P0*C*P0/2]:
        zero(expr)
    assert s.trace(Pi)==2
    # Include the two angular spinor components: full folded rank eight.
    angular2=lambda z: 2*s.trace(z)
    chip=s.diff(chi,th)*dth
    zero(angular2(chip*chip)-8*dth**2)
    zero(angular2(S*S)-dth**2)
    for z in [S,S**3,chip*chip*S,chi]: zero(angular2(z))
    # Product E is an arbitrary angular endomorphism; radial/flavor trace
    # of S and chi kills every such mixed term, including the magnetic part.
    e0,e1,e2,e3=s.symbols('e0 e1 e2 e3')
    Eang=e0*I2+e1*p1+e2*p2+e3*p3
    zero(s.trace(kron(S,Eang)))
    zero(s.trace(kron(chi,Eang)))
    a3=s.simplify(angular2(192*S*S-12*chip*chip)/384)
    a4normal=s.simplify(angular2(240*S*S-18*chip*chip)/360)
    zero(a3-dth**2/4)
    zero(a4normal-s.Rational(4,15)*dth**2)
    omitted=s.simplify(angular2(-12*chip*chip)/384)
    zero(omitted+dth**2/4)
    assert s.simplify(a3-omitted)==dth**2/2
    # Differentiate Pi Psi=Psi before squaring D_E. The minus projector
    # maps D_E Psi to dotPi Psi-i A Pi d_r Psi; the Robin term cancels it.
    zero(dotPi*Pi+s.I*A*S)
    # Direct interval propagation from the two original projector rows.
    h,z=s.symbols('h z',real=True)
    prop=s.cos(z)*I2+s.I*p2*s.sin(z)
    plus=s.Matrix([[1,1]])/s.sqrt(2)
    minus=s.Matrix([[1,-1]])/s.sqrt(2)
    condition=(plus*(prop+h*I2)).col_join(minus*(I2+h*prop))
    zero(condition.det()+((1+h*h)*s.cos(z)+2*h))
    # Original trace equations versus the orthogonal projector, without
    # assuming the free spectral formula.
    va=s.Matrix([1,1,-h,-h])
    vb=s.Matrix([-h,h,1,-1])
    chih=((1-h*h)*chi0-2*h*kron(p1,I2))/(1+h*h)
    zero(chih*va-va); zero(chih*vb-vb)
    zero(va.H*B*vb); zero(va.H*B*va); zero(vb.H*B*vb)
    return dict(rank=8,allowed_rank=4,domain_Green_and_work='PASS',
                static_collar_intertwiner='PASS for identical spherical collars',
                induced_mixed_Robin='PASS',product_Delta_a3_without_4pi='theta_prime**2/4',
                product_Delta_a4_normal_smearing_without_4pi='4*theta_prime**2/15',
                product_Delta_a4_unsmeared=0,omitting_Robin_a3='-theta_prime**2/4 (incorrect)',
                interval_propagation_determinant='-(1+h^2)*cos(pi*omega)-2*h')


def boundary_work_cases():
    rng=np.random.default_rng(742026)
    arr=lambda m:np.kron(np.array(m,dtype=complex),np.eye(2))
    an,bn,gn,cn,chi=map(arr,[A,B,G,C,chi0])
    identity=np.eye(8)
    p0=(identity+chi)/2
    cases=[]
    for theta in [-0.6,-0.1,0.,4*np.pi/2000,.2,.65]:
        for dtheta in [-.3,.4]:
            v=p0@(rng.normal(size=8)+1j*rng.normal(size=8))
            vd=p0@(rng.normal(size=8)+1j*rng.normal(size=8))
            V=expm(theta*gn/2); psi=V@v
            allowed=V@p0@V.conj().T
            forbidden=identity-allowed
            m=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
            angular=(m+m.conj().T)/2
            tangent=np.kron(np.array(kron(I2,p1),dtype=complex),angular)
            robin=dtheta*cn@allowed/2
            free_normal=forbidden@(rng.normal(size=8)+1j*rng.normal(size=8))
            analytic=dtheta*np.vdot(v,cn@v)/2
            def sample(e):return expm((theta+dtheta*e)*gn/2)@(v+e*vd)
            estimates=[]
            for eps in [1e-3,5e-4]:
                deriv=(-sample(2*eps)+8*sample(eps)-8*sample(-eps)+sample(-2*eps))/(12*eps)
                estimates.append(-1j*np.vdot(psi,bn@deriv))
            discrepancy=abs(estimates[-1]-analytic)
            flux=abs(np.vdot(psi,bn@psi))
            # Independently differentiate an allowed boundary trace and
            # evaluate the first-order Euclidean operator with its normal jet.
            normal=free_normal-robin@psi
            mixed_residual=np.linalg.norm(forbidden@(deriv-1j*an@normal+tangent@psi))
            omitted_residual=np.linalg.norm(forbidden@(deriv-1j*an@free_normal+tangent@psi))
            assert discrepancy<2e-9 and flux<2e-13
            assert mixed_residual<2e-9 and omitted_residual>1e-3
            assert abs(analytic)<=abs(dtheta)*np.vdot(psi,psi).real/2+1e-13
            cases.append(dict(theta=theta,theta_dot=dtheta,discrepancy=float(discrepancy),
                              refinement=float(abs(estimates[-1]-estimates[0])),charge_flux=float(flux),
                              boundary_power=float(analytic.real),mixed_domain_residual=float(mixed_residual),
                              domain_residual_if_Robin_omitted=float(omitted_residual)))
    return cases


def static_spectrum_cases():
    cases=[]
    al=np.array(p2,dtype=complex); identity=np.eye(2)
    ep=np.array([1,1])/np.sqrt(2); em=np.array([1,-1])/np.sqrt(2)
    for theta in [-.6,0.,4*np.pi/2000,.2,.65]:
        h=np.tan(theta/2)
        def determinant(w):
            U=expm(1j*np.pi*w*al)
            rows=np.vstack([ep@(U+h*identity),em@(identity+h*U)])
            return float(np.linalg.det(rows).real)
        err=0.; residual=0.
        for k in range(-8,9):
            root=brentq(determinant,k+.1,k+.9,xtol=1e-13)
            expected=k+.5+(-1.)**k*theta/np.pi
            err=max(err,abs(root-expected)); residual=max(residual,abs(determinant(expected)))
        assert err<2e-11 and residual<2e-12
        cases.append(dict(theta=theta,roots=17,max_frequency_error=err,max_direct_determinant_residual=residual))
    return cases


def heat_cases():
    mp.mp.dps=65
    cases=[]
    for th in ['-0.6','0.006283185307179586476925286766559','.2','.65']:
        for width in ['.12','.25','.5','1.0']:
            theta=mp.mpf(th); heat=mp.mpf(width); offset=mp.mpf('.5')+theta/mp.pi
            direct=2*mp.fsum(mp.exp(-heat*(2*n+offset)**2) for n in range(-80,81))
            baseline=2*mp.fsum(mp.exp(-heat*(2*n+mp.mpf('.5'))**2) for n in range(-80,81))
            factor=mp.pi**2/(4*heat)
            poisson=mp.sqrt(mp.pi/heat)*(1+2*mp.fsum(mp.exp(-factor*m*m)*mp.cos(mp.pi*m*offset) for m in range(1,81)))
            difference=abs(direct-baseline)
            bound=4*mp.sqrt(mp.pi/heat)*mp.exp(-factor)/(1-mp.exp(-3*factor))
            error=abs(direct-poisson)
            assert error<mp.mpf('1e-55') and difference<=bound
            cases.append(dict(theta=th,heat_time=width,poisson_discrepancy=str(error),
                              relative_heat_difference=str(difference),exponential_bound=str(bound)))
    assert mp.mpf(cases[-1]['relative_heat_difference'])>mp.mpf('.01')
    return cases


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=dict(status='PASS',analytic=analytic(),boundary_work=boundary_work_cases(),
                spectrum=static_spectrum_cases(),heat=heat_cases(),
                versions=dict(python=platform.python_version(),sympy=s.__version__,numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__),
                scope='Four-dimensional free domain and local ultraviolet theorem; driven heat coefficients only on the stated Euclidean product collar; no physical state or opening bound.',
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))
