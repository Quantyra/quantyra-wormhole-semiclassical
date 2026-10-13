"""S072 classical RN charge matching and uniform finite-cut radial midpoint.

No quantum gravitational density operator, driven stress or creation is
certified by these checks. See the ordinary proof for assumptions.
"""
from pathlib import Path
import hashlib
import json
import math
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.special import erfc
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'checks'))
from radial_electric_correlator import fem, lattice, rotor_theta


def identities():
    r, G, g, L, P, w, T = sy.symbols('r G g L P w T', positive=True)
    alpha = G*g**2/(4*sy.pi)
    M = r/(2*G)*(1+r*r/L**2+(P*P+alpha*w)/r**2)
    S = sy.pi*r*r/G
    temp = (1+3*r*r/L**2-(P*P+alpha*w)/r**2)/(4*sy.pi*r)
    assert sy.simplify(sy.diff(M, r)-temp*sy.diff(S, r)) == 0
    assert sy.simplify(sy.diff(M, w)-g*g/(8*sy.pi*r)) == 0
    poly = 3*r**4/L**2-4*sy.pi*T*r**3+r*r-P*P-alpha*w
    C = 1+6*r*r/L**2-6*sy.pi*T*r
    assert sy.simplify(sy.diff(poly, r)-2*r*C) == 0
    # Extremal entropy response: db/dw=alpha/(2b C_b).
    entropy_response = sy.diff(S, r)*alpha/(2*r*(1+6*r*r/L**2))
    assert sy.simplify(entropy_response-g*g/(4*(1+6*r*r/L**2))) == 0
    return dict(on_shell_radius_cancellation=True, electric_free_energy_derivative=True,
                horizon_equation_derivative=True, extremal_entropy_response=True)


def rn_numeric():
    mp.mp.dps = 70
    rows = []
    for radius, ads, temperature, coupling in ((3, 2, '.02', '.6'),
                                               (5, 1, '.03', '1.1'),
                                               (2, 2, '.01', '.3')):
        R, L, T, g = map(mp.mpf, (radius, ads, temperature, coupling))
        alpha = g*g/(4*mp.pi)
        P2 = R*R*(1+3*R*R/L**2-4*mp.pi*T*R)
        C = 1+6*R*R/L**2-6*mp.pi*T*R
        a, d = g*g/(8*mp.pi*R*T), alpha/(2*R*R*C)
        assert C > 0 and R >= mp.pi*T*L*L/2

        def r_of_w(w):
            fun = lambda r: 3*r**4/L**2-4*mp.pi*T*r**3+r*r-P2-alpha*w
            derivative = lambda r: 2*r*(1+6*r*r/L**2-6*mp.pi*T*r)
            return mp.findroot(fun, R, df=derivative, solver='newton') if w else R

        def free(r, w):
            return r/2*(1+r*r/L**2+(P2+alpha*w)/r**2)-mp.pi*T*r*r

        for Q in (1, 4, 15, 60):
            w = mp.mpf(Q)**2
            r = r_of_w(w)
            action = (free(r, w)-free(R, 0))/T
            integral = mp.quad(lambda u: g*g/(8*mp.pi*T*r_of_w(u)), [0, w])
            lower, upper = a/d*mp.log1p(d*w), a*w
            err, cap = upper-action, a*d*w*w/2
            assert abs(action-integral) < mp.mpf('1e-60')
            assert lower <= action <= upper and 0 <= err <= cap
            rows.append(dict(R=str(R), L=str(L), T=str(T), g=str(g), Q=Q,
                             action=str(action), integral_error=str(abs(action-integral)),
                             log_lower=str(lower), quadratic_upper=str(upper),
                             quadratic_error=str(err), error_bound=str(cap)))
    return rows


def classical_distribution_checks():
    rows = []
    for R, L, T, g in ((2., 2., .02, .5), (3., 2., .03, .8), (5., 2., .02, 1.)):
        alpha = g*g/(4*np.pi)
        P2 = R*R*(1+3*R*R/L**2-4*np.pi*T*R)
        C = 1+6*R*R/L**2-6*np.pi*T*R
        a, d = g*g/(8*np.pi*R*T), alpha/(2*R*R*C)
        p = a/d
        assert a <= 1 and p >= 100
        stop = math.ceil(math.sqrt(80/a))
        Q = np.arange(0, 3*stop+1, dtype=float)

        def horizon(w):
            if not w:
                return R
            return brentq(lambda r: 3*r**4/L**2-4*np.pi*T*r**3+r*r-P2-alpha*w,
                          R, R*(1+d*w)+1.e-10, xtol=1.e-13)

        # Independent thermodynamic charge integration avoids cancellation.
        action = np.array([quad(lambda u: a*R/horizon(u), 0, q*q,
                                epsabs=2.e-11, epsrel=2.e-12)[0] for q in Q])
        multiplicity = np.where(Q == 0, 1., 2.)
        actual = multiplicity*np.exp(-action)
        gaussian = multiplicity*np.exp(-a*Q*Q)
        actual /= sum(actual)
        gaussian /= sum(gaussian)
        sampled_tail = float(sum(actual[Q > stop]))
        bound = (2*np.sqrt(a/np.pi)*np.exp(-a*stop*stop/2)
                 +np.sqrt(2)*erfc(np.sqrt(a/2)*stop)
                 +2**(1-p)/np.sqrt(np.pi)*(np.sqrt(a)+np.sqrt(p)/(p-1)))
        # The analytic bound, not the finite charge sample, controls all Q.
        assert sampled_tail <= bound < np.exp(-39)
        central_error = max((a*Q*Q-action)[Q <= stop])
        assert central_error <= 5000/p
        tv = float(sum(abs(actual-gaussian))/2)
        assert tv <= np.expm1(5000/p)+2*np.exp(-39)
        rows.append(dict(R=R, L=L, T=T, g=g, a=a, p=p, Q_window=stop,
                         sampled_tail=sampled_tail, rigorous_tail_upper=bound,
                         sampled_tv=tv, central_log_error=central_error,
                         central_log_error_upper=5000/p))
    return rows


def radial_cut_checks():
    rows = []
    # Direct charge traces, not a large-n or infinite-volume replacement.
    for n in (1, 2, 3):
        for length, beta, mass in ((.7, 1.2, .8), (2., 1.1, 1.7), (5., 1., .3)):
            lam = n/np.pi
            c = np.pi*beta/(2*length)
            K0 = mass/np.tanh(mass*length)
            x = beta*K0/(2*lam)
            y = x-c/n
            j = np.arange(1, 100002, 2, dtype=float)
            s = 2*np.pi*j/beta
            sv = np.sqrt(s*s+mass*mass)
            D = 4/(beta*lam)*np.sum((sv/np.tanh(sv*length)-s/np.tanh(s*length))/s**2)
            ratio = np.exp(-D)*lattice(n,c,x,.5,16)/rotor_theta(c,.5)
            lo, hi = np.exp(-y/2), np.exp(y/2)/np.cosh(y/2)
            assert lo-2.e-10 <= ratio <= hi+2.e-10
            rows.append(dict(profile='constant', n=n, L_optical=length, beta=beta,
                             y=y, D=D, ratio=ratio, lower=lo, upper=hi))
    profile = lambda x: .8*np.exp(1.2*x)*(1+.2*np.sin(2*x)**2)
    length, beta = 3., 1.7
    M,H,K0,vals,co,contact = fem(length, 160, profile)
    _,_,Kfree,vfree,cfree,_ = fem(length, 160, lambda x: np.zeros_like(x))
    integral = quad(profile,-length,0)[0]
    assert K0-1/length <= integral
    for n in (1,2,3):
        lam=n/np.pi; c=np.pi*beta/(2*length); x=beta*K0/(2*lam); y=x-c/n
        w, w0 = np.sqrt(vals), np.sqrt(vfree)
        osc = np.sum(w*co**2*np.tanh(beta*w/4)-w0*cfree**2*np.tanh(beta*w0/4))/(2*lam)
        D = y/4+osc
        ratio=np.exp(-D)*lattice(n,c,x,.5,16)/rotor_theta(c,.5)
        lo,hi=np.exp(-y/2),np.exp(y/2)/np.cosh(y/2)
        assert lo <= ratio <= hi
        rows.append(dict(profile='nonconstant_finite_mesh',n=n,L_optical=length,beta=beta,
                         y=y,D=D,ratio=ratio,lower=lo,upper=hi,
                         K0_minus_free=K0-1/length,integrated_potential_upper=integral))
    return rows


def candidate_checks():
    mp.mp.dps = 400
    rows=[]
    for kk in (2000,10000,1000000):
        for zz in ('0.25','0.28'):
            k,z=mp.mpf(kk),mp.mpf(zz)
            g,L=k**-4,k**5
            X=(mp.sqrt(1+12*mp.pi*k**6)-1)/6
            b=L*mp.sqrt(X); Cb=1+6*X
            beta=2*mp.pi*b**3/(Cb*2*k**3*z)
            T=1/beta; alpha=g*g/(4*mp.pi); P2=mp.pi*k**16
            # Scale the root by b so all residuals remain well-conditioned.
            def equation(v,w=0):
                return 3*X*v**4-4*mp.pi*T*b*v**3+v*v-(1+3*X)-alpha*w/b**2
            v=mp.findroot(equation,(1,1+4*mp.pi*T*b/Cb))
            R=b*v; C=1+6*R*R/L**2-6*mp.pi*T*R
            a=beta*g*g/(8*mp.pi*R); d=alpha/(2*R*R*C); p=a/d
            p_lower=mp.pi**2*k**23/(12*z)
            assert b<=R<=2*b and C>=Cb/2 and R>=mp.pi*T*L**2/2
            assert 0<a<1/(12*k) and p>=p_lower
            Qstar=mp.ceil(mp.sqrt(80/a)); w=Qstar**2
            vw=mp.findroot(lambda v:equation(v,w),(v,v*(1+d*w)))
            Rw=b*vw
            def free(r, w):
                return r/2*(1+r*r/L**2+(P2+alpha*w)/r**2)-mp.pi*T*r*r
            action=beta*(free(Rw,w)-free(R,0))
            error=a*w-action
            cap=a*d*w*w/2
            assert 0<error<=cap and cap<5000/p
            assert abs(error/cap-1)<mp.mpf('1e-50')
            entropy=g*g/(4*Cb)
            extremal_a=beta*g*g/(8*mp.pi*b)
            # First thermal correction equals -dS_ext/d(Q^2).
            thermal_ratio=(extremal_a-a)/entropy
            assert abs(thermal_ratio-1)<mp.mpf('1e-28')
            rows.append(dict(k=kk,z0=zz,R_over_b=str(v),a=str(a),p=str(p),
                             p_lower=str(p_lower),Q_window=str(Qstar),
                             actual_central_action_error=str(error),
                             pointwise_error_bound=str(cap),window_error_upper=str(5000/p),
                             entropy_cancellation_ratio=str(thermal_ratio),
                             finite_cut_loss_upper=str(-mp.expm1(-1/(24*k)))))
    return rows


def main():
    receipt=dict(date='2026-10-13',status='PASS',
                 scope='Classical fixed-temperature RN charge thermodynamics and neutral-inner-cut radial LLL comparison',
                 identities=identities(),classical_RN=rn_numeric(),
                 classical_charge_probabilities=classical_distribution_checks(),
                 finite_cut_radial=radial_cut_checks(),candidate=candidate_checks(),
                 versions=dict(numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,sympy=sy.__version__),
                 limitations=['No full gravity/matter density operator or charged matrix-element matching',
                              'No quantum determinant or angular-sector acceptance',
                              'Finite-cut bound uses the same free reference and neutral inner flux',
                              'No driven stress, stability, preparation or physical creation'])
    files=['notes/electric-thermodynamic-matching.md','checks/electric_thermodynamic_matching.py',
           'checks/radial_electric_correlator.py','notes/radial-electric-correlator.md',
           'notes/electric-reservoir-midpoint.md']
    receipt['source_hashes']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}
    (ROOT/'evidence/electric-thermodynamic-matching.json').write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('PASS: 4 symbolic identities, 12 independent RN action comparisons,')
    print('3 classical charge distributions, 12 finite-cut radial comparisons, 6 candidate cases')
    print('largest classical numerical action error',max(float(r['integral_error']) for r in receipt['classical_RN']))
    print('k=2000 window log-weight bound',mp.nstr(mp.mpf(receipt['candidate'][1]['window_error_upper']),8))


if __name__=='__main__':
    main()
