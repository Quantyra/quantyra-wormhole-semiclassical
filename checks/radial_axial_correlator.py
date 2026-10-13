"""S072 radial electric plus one commuting angular generator.

Matrix DtN, independent ODE, actual vector-charge Fock traces and joint
charge lattice. The noncommuting angular theory is outside this model.
"""
from pathlib import Path
from itertools import product
import hashlib
import json
import math
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import eigh, solve
from scipy.integrate import solve_ivp, quad
from scipy.special import polygamma

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'checks'))
from radial_electric_correlator import fem


def matrix_fem(length, count, potential):
    nodes = np.linspace(-length, 0, count+1)
    h = length/count
    M = np.zeros((2*(count+1),)*2)
    H = M.copy()
    points, weights = np.polynomial.legendre.leggauss(6)
    shape = np.stack(((1-points)/2, (1+points)/2))
    for n in range(count):
        xx = (nodes[n]+nodes[n+1])/2+h*points/2
        local_m = np.kron(h*np.array([[2., 1.], [1., 2.]])/6, np.eye(2))
        local_h = np.kron(np.array([[1., -1.], [-1., 1.]])/h, np.eye(2))
        for s, x in enumerate(xx):
            local_h += h/2*weights[s]*np.kron(np.outer(shape[:, s], shape[:, s]), potential(x))
        sl = slice(2*n, 2*n+4)
        M[sl, sl] += local_m
        H[sl, sl] += local_h
    # Both fields vanish at the inner cut. The final light value is fixed;
    # the heavy outer value is varied, giving its natural Neumann condition.
    M, H = M[2:, 2:], H[2:, 2:]
    z = np.r_[-solve(H[:-1, :-1], H[:-1, -1], assume_a='pos'), 1.]
    K0 = float(z@H@z)
    vals, vecs = eigh(H[:-1, :-1], M[:-1, :-1])
    coeffs = vecs.T@(M@z)[:-1]
    contact = float(z@M@z-coeffs@coeffs)
    return M, H, K0, vals, coeffs, contact


def ode_dtn(length, potential, s):
    eps = 1e-6
    V = potential(-length)
    initial = np.eye(2)/eps+(V+s*s*np.eye(2))*eps/3
    def rhs(x, yy):
        y = yy.reshape(2, 2)
        return (potential(x)+s*s*np.eye(2)-y@y).reshape(-1)
    sol = solve_ivp(rhs, (-length+eps, 0), initial.reshape(-1), rtol=3e-11, atol=2e-12)
    assert sol.success
    Y = sol.y[:, -1].reshape(2, 2)
    assert abs(Y[0, 1]-Y[1, 0]) < 1e-11
    return float(Y[1, 1]-Y[1, 0]*Y[0, 1]/Y[0, 0])


def rn_profile(q=3, R=2., G=.2, g=1., X=.6, length=3.):
    Y = math.pi*G*q*q/(g*g*R*R)
    assert 0 < Y < 1+3*X
    F = lambda u: X+u*u-(1+X+Y)*u**3+Y*u**4
    radial = solve_ivp(lambda x, u: [-F(u[0])/R], (0, -length), [0.],
                       rtol=2e-12, atol=2e-13, dense_output=True)
    assert radial.success
    lam = q*(q*q-1)/(12*math.pi)
    mu, P = 8*math.pi/(3*g*g), q/2
    def potential(x):
        u = float(radial.sol(float(x))[0])
        pref = F(u)/(R*R)
        A = pref*(2+4*Y*u*u)
        B = pref*math.sqrt(mu*lam)/(q+1)*(-6*P*G*u*u/(R*R)+3*g*g/(2*math.pi))
        D = pref*lam/(q+1)**2*(6*G*u*u/(R*R)+3*g*g/math.pi)
        return np.array([[A, B], [B, D]])
    def electric(x):
        u = radial.sol(x)[0]
        return F(u)/(R*R)*(q/math.pi)*g*g/(4*math.pi)
    rcut = R/float(radial.sol(-length)[0])
    W = lam/(q+1)**2*(2*G/rcut**3+3*g*g/(math.pi*rcut))
    exact = quad(lambda x: potential(x)[1, 1], -length, 0, epsabs=1e-12)[0]
    assert abs(exact-W) < 2e-12
    return potential, electric, lam, W, rcut


def thermal_difference(current, reference, beta, lam):
    M, H, K, v, c, b = current
    M0, H0, K0, v0, c0, b0 = reference
    assert abs(b-b0) < 5e-13
    w, w0 = np.sqrt(v), np.sqrt(v0)
    direct = beta*(K-K0)/(8*lam)+(np.sum(w*c*c*np.tanh(beta*w/4))
                -np.sum(w0*c0*c0*np.tanh(beta*w0/4)))/(2*lam)
    j = np.arange(1, 20002, 2, dtype=float)
    s = 2*math.pi*j/beta
    dk = np.full_like(s, K-K0)
    for start in range(0, len(s), 500):
        ss = s[start:start+500, None]**2
        dk[start:start+500] += np.sum(ss*v*c*c/(ss+v), axis=1)
        dk[start:start+500] -= np.sum(ss*v0*c0*c0/(ss+v0), axis=1)
    mats = 4/(beta*lam)*np.sum(dk/s**2)
    dk_inf = K-K0+np.sum(v*c*c)-np.sum(v0*c0*c0)
    next_odd = j[-1]+2
    tail2 = float(polygamma(1, next_odd/2))/4
    tail4 = float(polygamma(3, next_odd/2))/96
    mats += beta*dk_inf/(math.pi**2*lam)*tail2
    rem = beta**3/(4*math.pi**4*lam)*tail4*(np.sum(v*v*c*c)+np.sum(v0*v0*c0*c0))
    assert abs(direct-mats) <= rem+8e-11
    return float(direct), float(mats), float(rem), float(abs(b-b0))


def joint_theta(q, c, ye, ya, tag, t=.5, cutoff=10):
    N = np.array(list(product(range(-cutoff, cutoff+1), repeat=q)), dtype=float)
    m = np.arange(q)-(q-1)/2
    A = c*np.eye(q)+ye*np.ones((q,q))+ya*np.outer(m,m)
    E = np.einsum('ni,ij,nj->n', N, A, N)
    shifted = N.copy(); shifted[:, tag] += t
    Ep = np.einsum('ni,ij,nj->n', shifted, A, shifted)
    return float(np.exp(-Ep).sum()/np.exp(-E).sum())


def radial_checks():
    rows, joint = [], []
    length, beta = 3., 1.7
    for q, G, X in [(2, .2, .6), (3, .2, .6)]:
        pot, epot, lam, W, rcut = rn_profile(q=q, G=G, X=X, length=length)
        ref = lambda x: np.diag([pot(x)[0, 0], 0.])
        samples = [0., .2, 1., 4.]
        exact = [ode_dtn(length, pot, s) for s in samples]
        delta = [v-(1/length if s == 0 else s/np.tanh(s*length)) for v,s in zip(exact,samples)]
        assert all(0 <= v <= W for v in delta)
        previous = None
        for count in [32,64,128]:
            current, reference = matrix_fem(length,count,pot), matrix_fem(length,count,ref)
            M,H,K,vals,coeffs,contact=current
            primal_errors, continuum_errors=[] ,[]
            for s, truth in zip(samples,exact):
                F=H+s*s*M
                primal=F[-1,-1]-F[-1,:-1]@solve(F[:-1,:-1],F[:-1,-1],assume_a='pos')
                spectral=K+s*s*contact+np.sum(s*s*vals*coeffs**2/(s*s+vals))
                primal_errors.append(float(abs(primal-spectral)))
                continuum_errors.append(float(abs(primal-truth)))
            assert max(primal_errors)<3e-11
            D, mats, rem, contact_err=thermal_difference(current,reference,beta,lam)
            assert 0<D<beta*W/(8*lam)
            row=dict(q=q,elements=count,K0=K,D_midpoint=D,matsubara=mats,
                     thermal_error=abs(D-mats),tail_bound=rem,
                     spectral_error=max(primal_errors),continuum_error=max(continuum_errors),
                     contact_difference=contact_err,radial_integral=W,
                     bare_midpoint_cap=beta*W/(8*lam),continuum_delta_K=delta,
                     refinement=None if previous is None else abs(D-previous))
            rows.append(row); previous=D
            print('radial',row,flush=True)
        recent=rows[-3:]
        assert recent[2]['refinement']<.4*recent[1]['refinement']
        assert recent[2]['continuum_error']<.4*recent[1]['continuum_error']
        le=q/math.pi
        electric=fem(length,128,epot)
        efree=fem(length,128,lambda x:np.zeros_like(x))
        De,_,_,_=thermal_difference(electric,efree,beta,le)
        c=math.pi*beta/(2*length)
        ye=beta*(electric[2]-1/length)/(2*le)
        ya=beta*(K-1/length)/(2*lam)
        be=beta/(8*math.pi*rcut)
        ba=beta*W/(2*lam)
        for tag in range(q):
            m=tag-(q-1)/2
            theta=joint_theta(q,c,ye,ya,tag)/joint_theta(q,c,0,0,tag)
            refined=joint_theta(q,c,ye,ya,tag,cutoff=12)/joint_theta(q,c,0,0,tag,cutoff=12)
            assert abs(theta-refined)<2e-12
            ratio=math.exp(-De-m*m*D)*theta
            lamtag=ye+m*m*ya; btag=be+m*m*ba
            assert lamtag<=btag+2e-6 # FEM upper error shrinks with refinement.
            lower=math.exp(-btag/2)
            upper=math.exp(btag/2)/math.cosh(btag/2)
            assert lower<=ratio<=upper
            joint.append(dict(q=q,m=m,ratio=ratio,lower=lower,upper=upper,
                              screened_increment=lamtag,bare_increment=btag))
    return rows,joint


def vector_charge_checks():
    rows=[]
    for q, energies, ye, ya in [(2,[.3,1.2],.2,.7),(3,[0.,.8],.4,.2),(3,[1.1],.05,1.2)]:
        m=np.arange(q)-(q-1)/2
        charges=np.array([[1.,mm] for u in energies for mm in m])
        u=np.repeat(energies,q)
        D=len(u)
        # Adjacent particle/antiparticle pairs; all actual occupations retained.
        ee=np.repeat(u,2)
        vv=np.stack([v for v in charges for v in (v,-v)])
        states=np.array(list(product([0.,1.],repeat=2*D)))
        Q=states@vv
        H=states@ee+ye*Q[:,0]**2+ya*Q[:,1]**2
        Z=np.exp(-H).sum()
        for tag in [0,D-1]:
            col=2*tag
            initial=states[:,col]==0
            target=states[initial].copy(); target[:,col]=1
            Qt=target@vv
            Ht=target@ee+ye*Qt[:,0]**2+ya*Qt[:,1]**2
            b=np.exp(-(H[initial]+Ht)/2).sum()/Z
            z=math.exp(-u[tag]); ratio=b/(math.sqrt(z)/(1+z))
            v=charges[tag]; lam=ye*v[0]**2+ya*v[1]**2
            assert math.exp(-lam/2)<=ratio<=math.exp(lam/4)/math.cosh(lam/2)
            # Pair removed, then direct vector-Gaussian averaging gives a
            # separate formula for the insertion including both tagged occupations.
            keep=[i for i in range(2*D) if i not in [col,col+1]]
            rest=np.array(list(product([0.,1.],repeat=len(keep))))
            qr=rest@vv[keep]; weight=np.exp(-rest@ee[keep]); weight/=weight.sum()
            def f(s):
                shifted=qr+s*v
                return np.sum(weight*np.exp(-ye*shifted[:,0]**2-ya*shifted[:,1]**2))
            formula=(1+z)**2*math.exp(-lam/4)*f(.5)/((1+z*z)*f(0)+2*z*f(1))
            assert abs(formula-ratio)<3e-13
            rows.append(dict(q=q,energies=energies,tag=tag,ye=ye,ya=ya,
                             exact_ratio=float(ratio),formula=float(formula),error=float(abs(formula-ratio))))
    return rows


def matrix_counterexample():
    V=np.array([[1.,1.],[1.,1.0001]])
    w,U=eigh(V); length=40.
    def K(s):
        roots=np.sqrt(w+s*s)
        Y=(U*(roots/np.tanh(length*roots)))@U.T
        return float(Y[1,1]-Y[0,1]**2/Y[0,0])
    d0=K(0)-1/length
    d1=K(.4)-.4/np.tanh(.4*length)
    assert w[0]>0 and d1>5*d0>0
    return dict(potential=V.tolist(),length=length,delta_K_zero=d0,
                delta_K_at_point4=d1,scalar_monotonicity_extension=False)


def candidate_checks():
    mp.mp.dps=70
    rows=[]
    for kk in [2000,10000,1000000]:
        for zz in ['0.25','0.28']:
            k,z=mp.mpf(kk),mp.mpf(zz)
            q=k**4; j=(q-1)/2; g=k**-4; L=k**5
            X=(mp.sqrt(1+12*mp.pi*k**6)-1)/6
            b=L*mp.sqrt(X); Y=1+3*X
            beta=mp.pi*b**3/((1+6*X)*k**3*z)
            be=beta*g*g/(8*mp.pi*b)
            ba=beta/(2*(q+1)**2)*(2/b**3+3*g*g/(mp.pi*b))
            full=be+j*j*ba
            alternative=be*(1+j*j/(q+1)**2*(12+8*Y/q**2))
            assert abs(full-alternative)<mp.mpf('1e-60')
            assert be<1/(48*k*z) and 2*Y/q**2<mp.mpf('1e-8')
            cap=(4+mp.mpf('1e-8'))/(24*k)
            loss=-mp.expm1(-full/2)
            assert full/2<cap and loss<mp.mpf('8.334e-5')
            rows.append(dict(k=kk,z0=zz,bare_exponent=str(full/2),loss_bound=str(loss),uniform_cap=str(cap)))
    return rows


def main():
    counter=matrix_counterexample()
    print('matrix static-cap counterexample',counter,flush=True)
    vector=vector_charge_checks()
    radial,joint=radial_checks()
    result=dict(date='2026-10-13',status='PASS',
                scope='Fixed radial electric plus one commuting angular generator, quadratic Einstein-Maxwell sector, neutral inner cut',
                matrix_counterexample=counter,vector_charge_Fock=vector,
                radial_matrix=radial,joint_midpoints=joint,candidate=candidate_checks(),
                versions=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                limitations=['No noncommuting SU(2) interaction bound','No physical horizon/gluing or 4D vertex match','No driven stress, evolving geometry, stability, preparation or transport'])
    (ROOT/'evidence/radial-axial-correlator.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('PASS',len(vector),'vector traces;',len(radial),'radial meshes;',len(joint),'joint midpoints; 6 candidate cases',flush=True)


if __name__=='__main__':
    main()
