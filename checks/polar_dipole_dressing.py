"""Reciprocal dipole-pair bubble: a specified quadratic comparison, not full stability."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'checks'))
from dipole_vacuum_response import rho


def algebra():
    C,b,G,e,z,Q,y1,y2=s.symbols('C b G e z Q y1 y2',positive=True)
    k2=8*s.pi*G;a2=b*b/C
    W=s.Matrix([[2*b*b/k2*(1+C/3),2*b*b*e],[2*b*b*e,6]])
    K=s.Matrix([[(6+2*C)/b**2,6*k2*e/b**2],[2*e,6/b**2]])
    v=s.Matrix([9*G/(4*b*b),-3*G*e/4])
    s0=s.Matrix([3*G*(C-4)/b**4,9*G*e/b**2])
    R=s.Matrix([3*G*(C-1)/b**4,9*G*e/b**2])
    back=lambda a:s.factor(s.expand(a).subs(e**2,(C+1)/(k2*b*b)))
    assert (K*v+s0-R).applyfunc(back)==s.zeros(2,1)
    assert (W*K-K.T*W).applyfunc(back)==s.zeros(2)
    A=27*G/(32*s.pi*b*b);B=9*G*(2*C-5)/(16*s.pi*b**4)
    assert back((v.T*W*v)[0]-A)==0
    assert back((v.T*W*R)[0]-(v.T*W*K*v)[0]/2-B)==0
    # Scalar z represents the self-adjoint Box in the integrated identity.
    y=s.Matrix([y1,y2]);u=y+v*Q
    original=(u.T*W*(z*s.eye(2)-K)*u)[0]/2-(u.T*W*(v*z+s0))[0]*Q
    transformed=(y.T*W*(z*s.eye(2)-K)*y)[0]/2-(y.T*W*R)[0]*Q-A*z*Q**2/2-B*Q**2
    assert back(original-transformed)==0
    weighted=3*G*(C*C+25*C+28)/(8*s.pi*C*b*b)
    norm=3*G*(C**3+19*C*C+76*C+66)/(4*s.pi*C*C*b*b)
    assert back(a2*(R.T*W*K.inv()*R)[0]-weighted)==0
    assert back(a2*a2*(R.T*W*R)[0]-norm)==0
    x,j,t=s.symbols('x j t',nonnegative=True)
    assert s.expand((1+6*x)**2-(1+12*x+12*x*x)-24*x*x)==0
    assert s.expand((2*t+2*j)*(2*t+2*j-1)-2*t*(2*t-1)-4*j*j-(8*t-2)*j)==0
    beta_bound=s.Rational(3,2)+s.Rational(176,675)
    assert beta_bound<2
    return {'field_shift_quartic_derivative':str(A/2),'field_shift_quartic_potential':str(-B),'g_Minv_g':str(weighted),'g_squared':str(norm),'beta_rational_upper_bound':str(beta_bound),'identities':'exact'}


def data(C,gamma,N):
    C=mp.mpf(C);gamma=mp.mpf(gamma)
    radical=mp.sqrt(C*C+12*C+12)
    masses=[24/(C*(C+6+radical)),(C+6+radical)/C]
    sn=mp.sqrt(6*(C+1)/(radical*(radical+C)))
    cs=mp.sqrt((radical+C)/(2*radical))
    d=mp.sqrt(gamma/(8*mp.pi))/C**mp.mpf('1.5')*mp.matrix([mp.sqrt(3)*(C-1),9*mp.sqrt(C+1)])
    U=mp.matrix([[-sn,cs],[cs,sn]])
    g=mp.matrix([mp.sqrt(C*masses[r])*(U[:,r].T*d)[0] for r in range(2)])
    delta=(1+mp.sqrt(9+16/C))/2
    lam=[(2*delta+2*j)*(2*delta+2*j-1) for j in range(N)]
    weight=[16*mp.pi/15*rho(delta,j) for j in range(N)]
    q=3*gamma/(8*mp.pi)*(1+25/C+28/C**2)
    assert abs(mp.fsum(g[r]**2/masses[r] for r in range(2))/q-1)<mp.mpf('1e-55')
    return masses,g,delta,lam,weight,q


def spectral_bounds():
    mp.mp.dps=75;rows=[];candidate=None
    k=mp.mpf(2000);candidateC=mp.sqrt(1+12*mp.pi*k**6)
    for C in [mp.mpf(1),mp.mpf(4),mp.mpf(100),candidateC]:
        gamma=6*C/(k**10*(C-1)) if C==candidateC else mp.mpf('.001')
        masses,g,delta,lam,w,q=data(C,gamma,256)
        rho0=15*w[0]/(16*mp.pi)
        assert rho0<=11/(6*mp.sqrt(3*mp.pi))
        for j in range(1,256):
            assert 15*w[j]/(16*mp.pi)<=15/(4*mp.pi*j)
            assert lam[j]-masses[1]>=10+4/C+14*j+4*j*j
        for z in [mp.mpf(0),masses[0],masses[1]]:
            val=mp.fsum(wj/(lj-z) for wj,lj in zip(w,lam))
            deriv=mp.fsum(wj/(lj-z)**2 for wj,lj in zip(w,lam))
            # rho_j<=15/(4pi j), gap>=4j^2; integral remainder from N-1.
            tail=1/(2*mp.mpf(255)**2)
            derivative_tail=tail/10
            assert val+tail<2 and deriv+derivative_tail<mp.mpf('.2')
            assert (val+tail)*q<2*q
            rows.append({'C':str(C),'z':str(z),'terms':256,'beta_lower':str(val),'beta_upper':str(val+tail),'beta_prime_upper':str(deriv+derivative_tail)})
        if C==candidateC:
            candidate={'C':str(C),'G_over_a2':str(gamma),'mu_minus':str(masses[0]),'mu_plus':str(masses[1]),'epsilon':str(2*q),'uniform_14G_over_a2':str(14*gamma),'light_absolute_shift_bound':str(2*q*masses[0]),'g_minus':str(g[0]),'g_plus':str(g[1])}
    return {'series':rows,'candidate_k2000':candidate,'tail_proof':'beta_N <= beta <= beta_N + 1/[2(N-1)^2], uniformly 0<=z<=mu_plus; derivative tail <= this /10.'}


def bisect(f,left,right):
    fl=f(left);fr=f(right);assert fl*fr<0
    for _ in range(200):
        mid=(left+right)/2;fm=f(mid)
        if fl*fm>0:left=mid;fl=fm
        else:right=mid
    return (left+right)/2


def finite_bath_checks():
    rows=[]
    for C in [1,4,25]:
        for gamma in ['.001','.01']:
            for N in [16,48]:
                masses,g,delta,lam,w,q=data(C,gamma,N);epsilon=2*q
                D=np.array([float(v) for v in masses+lam]);coupling=np.outer(np.array(list(g),dtype=float).reshape(2),np.sqrt(np.array(w,dtype=float)))
                A=np.diag(D);A[:2,2:]=coupling;A[2:,:2]=coupling.T
                eigen,V=np.linalg.eigh(A)
                assert eigen[0]>0
                scaled=A/np.sqrt(D[:,None]*D[None,:])
                b0=mp.fsum(wj/lj for wj,lj in zip(w,lam))
                expected=1-math.sqrt(float(b0*q))
                assert abs(np.linalg.eigvalsh(scaled)[0]-expected)<2e-13
                def det(z):
                    B=mp.fsum(wj/(lj-z) for wj,lj in zip(w,lam))
                    return (masses[0]-z)*(masses[1]-z)-B*(g[0]**2*(masses[1]-z)+g[1]**2*(masses[0]-z))
                roots=[bisect(det,(1-epsilon)*mu,mu) for mu in masses]
                error=max(abs(eigen[r]-float(roots[r])) for r in range(2))
                assert error<2e-9
                norm_errors=[]
                for r in range(2):
                    z=eigen[r];u=V[:2,r]
                    Bp=sum(float(wj)/(float(lj)-z)**2 for wj,lj in zip(w,lam))
                    predicted=float(u@u+Bp*(np.array(list(g),dtype=float).reshape(2)@u)**2)
                    norm_errors.append(abs(predicted-1))
                    assert abs(predicted-1)<2e-10
                    assert float((1-epsilon)*masses[r])-2e-9<=z<=float(masses[r])+2e-9
                rows.append({'C':C,'G_over_a2':gamma,'bath_terms':N,'epsilon':str(epsilon),'mass_roots':[str(x) for x in roots],'matrix_root_error':error,'residue_norm_error':max(norm_errors),'positive_form_lowest_eigenvalue':expected})
    return rows


def gaussian_state_checks():
    # A separate spatial Galerkin Hamiltonian tests the stationary correlated state.
    rows=[];nr=7;N=6
    Vrad=np.zeros((nr,nr))
    for i in range(nr):
        for j in range(nr):
            def integrand(x):
                return 2/math.pi*math.sin((i+1)*(x+math.pi/2))*math.sin((j+1)*(x+math.pi/2))/math.cos(x)**2
            Vrad[i,j]=quad(integrand,-math.pi/2,math.pi/2,epsabs=1e-10,epsrel=1e-10)[0]
    T=np.diag(np.arange(1,nr+1,dtype=float)**2)
    for C in [1,4]:
        masses,g,delta,lam,w,q=data(C,'.01',N)
        A=np.diag([float(x) for x in masses+lam]);c=np.outer(np.array(list(g),dtype=float).reshape(2),np.sqrt(np.array(w,dtype=float)))
        A[:2,2:]=c;A[2:,:2]=c.T
        H=np.kron(T,np.eye(N+2))+np.kron(Vrad,A)
        eig,E=np.linalg.eigh(H);assert eig[0]>0
        qcov=(E*(1/(2*np.sqrt(eig))))@E.T
        pcov=(E*(np.sqrt(eig)/2))@E.T
        purity=float(np.max(abs(qcov@pcov-np.eye(len(eig))/4)))
        stationary=float(np.max(abs(H@qcov-pcov)))
        assert purity<2e-12 and stationary<2e-10
        # Independent internal-mass rotation into uncoupled radial Hamiltonians.
        ma,U=np.linalg.eigh(A);reconstructed=np.zeros_like(qcov)
        for j in range(N+2):
            freq,R=np.linalg.eigh(T+ma[j]*Vrad)
            spatial=(R*(1/(2*np.sqrt(freq))))@R.T
            reconstructed+=np.kron(spatial,np.outer(U[:,j],U[:,j]))
        difference=float(np.max(abs(reconstructed-qcov)))
        assert difference<2e-12
        cross=float(np.max(abs(qcov.reshape(nr,N+2,nr,N+2)[:,:2,:,2:])))
        assert cross>1e-8
        rows.append({'C':C,'bath_terms':N,'radial_basis_size':nr,'stationarity_error':stationary,'purity_error':purity,'independent_covariance_error':difference,'nonzero_polar_bath_cross_covariance':cross})
    return rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    out={'story':'S072','scope':'Specified reciprocal nondifferentiated cubic representative, free dipole pair bubble and quadratic auxiliary-field comparison; omitted local action terms remain required.',
         'algebra':algebra(),'spectral':spectral_bounds(),'finite_bath':finite_bath_checks(),'Gaussian_state':gaussian_state_checks()}
    out['versions']={'sympy':s.__version__,'mpmath':mp.__version__,'numpy':np.__version__,'scipy':scipy.__version__}
    out['sources']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['checks/polar_dipole_dressing.py','notes/polar-dipole-dressing.md','checks/dipole_vacuum_response.py','notes/dipole-quadrupole-source.md','notes/polar-canonical-vacuum.md']}
    out['status']='PASS';out['goal_complete']=False
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_bytes((json.dumps(out,indent=2)+'\n').encode())
    print('PASS: reciprocal bubble, relative light-mass bound, retained quartic shift and correlated quadratic vacuum',flush=True)


if __name__=='__main__':main()
