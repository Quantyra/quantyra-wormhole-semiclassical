"""Exact and independent integral checks for the dipole midplane audit.

No four-dimensional reduction, novelty verdict, or Lean invocation.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import betaln, expit
import sympy as sp


def verify():
    r, a, q, j, t, x = sp.symbols('r a q j t x', positive=True)
    m = q-j-1
    R = sp.sqrt(r*r+a*a)
    u = a/R
    A = q*(1-u)
    Phi = q*sp.log((R+a)/(2*a))
    logf = m*sp.log(r)-Phi
    Eu_inv = (2*q*q-j)/(q*(2*j-1))
    variance_weighted = j*(q*q-j*j)/(q**3*(2*j-1))
    identities = {
        'magnetic_field': sp.diff(A,r)/r-q*a/R**3,
        'scalar_potential': sp.diff(Phi,r)-A/r,
        'zero_mode': sp.diff(logf,r)-(m-A)/r,
        'dilation_derivative': r*sp.diff(logf,r)+1-(q*u-j),
        'dirac_corrector': (1+r*sp.diff(logf,r)+m+1-A)/2-(q*u-j),
        'weighted_variance': (j/q)**2*Eu_inv-j/q-variance_weighted,
        'inverse_moment_excess': Eu_inv-q/j-(q*q-j*j)/(q*j*(2*j-1)),
        'variance_upper_bound': q*q*(2*j-1)-j*(q*q-j*j)-(q*q*(j-1)+j**3),
        'radial_coordinate': (2*a*sp.sqrt(t)/(1-t))**2+a*a-(a*(1+t)/(1-t))**2,
        'radial_jacobian': sp.diff((2*a*sp.sqrt(t)/(1-t))**2,t)/2-2*a*a*(1+t)/(1-t)**3,
    }
    assert all(sp.simplify(v)==0 for v in identities.values())

    # Beta recurrences written as rational products, valid for q>j>0.
    alpha, beta = q-j, 2*j
    zratio = 1+alpha/(alpha+beta)
    eu = (beta/(alpha+beta))/zratio
    einv = ((alpha+beta-1)/(beta-1)+2*alpha/(beta-1)
            +alpha*(alpha+1)/((beta-1)*(alpha+beta)))/zratio
    er2 = 4*(alpha*(alpha+beta-1)/((beta-1)*(beta-2))
             +alpha*(alpha+1)/((beta-1)*(beta-2)))/zratio
    beta_identities = {
        'normalization_ratio':zratio-2*q/(q+j),
        'mean_u':eu-j/q,
        'mean_inverse_u':einv-Eu_inv,
        'second_radius_moment':er2-4*(q*q-j*j)/((2*j-1)*(2*j-2)),
    }
    assert all(sp.factor(v)==0 for v in beta_identities.values())
    bump=(x-1)**2*(2-x)**2
    bump_norm=sp.integrate(x*bump**2,(x,1,2))
    bump_gradient=sp.integrate(x*sp.diff(bump,x)**2,(x,1,2))
    gap_constant=sp.factor(2*bump_gradient/bump_norm)

    cases=[]
    for Q in [2,3,8,25,100]:
        for J in sorted({1,min(2,Q-1),max(1,Q//2),Q-1}):
            M=Q-J-1
            logZ=np.log(2*Q/(Q+J))+betaln(Q-J,2*J)
            logden=np.log(2)+M*np.log(4)+logZ
            def density(y):
                if y==0: return 0.
                return np.exp((2*M+1)*np.log(y)-2*Q*np.log((np.hypot(1,y)+1)/2)-logden)
            funcs=[lambda y:1.,lambda y:1/np.hypot(1,y),lambda y:np.hypot(1,y)]
            expected=[1.,J/Q,(2*Q*Q-J)/(Q*(2*J-1))]
            if J>=2:
                funcs.append(lambda y:y*y)
                expected.append(4*(Q*Q-J*J)/((2*J-1)*(2*J-2)))
            values=[quad(lambda y:f(y)*density(y),0,np.inf,epsabs=2e-9,epsrel=2e-11)[0] for f in funcs]
            rel=max(abs(v-w)/max(1,abs(w)) for v,w in zip(values,expected))
            assert rel<2e-9,(Q,J,values,expected)
            cases.append({'q':Q,'j':J,'radial_integrals':values,'max_relative_residual':rel})

    # Critical-mode cutoff norm and error, using exact radial integration.
    tail_checks=[]
    for Q in [2,8,25]:
        logZ=np.log(2/(Q*Q-1))
        def pt(v):
            if v==0: return 0.
            return np.exp((Q-2)*np.log(v)+np.log1p(-v)+np.log1p(v)-logZ)
        for L in [10.,100.,1000.]:
            tL=(np.hypot(1,L)-1)/(np.hypot(1,L)+1)
            tail=quad(pt,tL,1,epsabs=2e-12)[0]
            def cutoff(y):
                if y<=L:return 1.
                if y>=2*L:return 0.
                return np.cos(np.pi*(y/L-1)/2)**2
            t2=(np.hypot(1,2*L)-1)/(np.hypot(1,2*L)+1)
            norm2=quad(lambda v:cutoff(2*np.sqrt(v)/(1-v))**2*v/(1-v)**2*pt(v),0,t2,
                       points=[tL],epsabs=2e-8,epsrel=2e-10)[0]
            assert tail<=2*(Q*Q-1)/L**2+1e-12
            assert norm2<=(Q*Q-1)*np.log1p(L)+1e-8
            tail_checks.append({'q':Q,'L_over_a':L,'tail_probability':tail,'cutoff_corrector_norm_squared_over_a_squared':norm2})

    # Finite imaginary-frequency band: sharp low-band delay bound plus
    # generic high-band front bound. Piecewise envelopes are not claimed
    # to be physical analytic scattering matrices.
    band_checks=[]
    kappa=8*np.log(2)-5
    for Q in [21,50,100]:
        labels=np.arange(1,Q+1,dtype=float)/Q
        labels[0]=1.  # q-1 modes, critical j=1 reserved, missing mode at u=1
        delays=1/labels
        for T in [.1,.5,1.]:
            low=quad(lambda xi:2/np.pi*xi*np.sum(expit(-xi*(T+delays))),0,4,epsabs=1e-10)[0]
            high=quad(lambda xi:2*Q/np.pi*xi*expit(-xi*(T+1)),4,np.inf,epsabs=1e-10)[0]
            full=np.pi/6*np.sum(1/(T+delays)**2)
            uv=2*Q/np.pi*np.exp(-4*(T+1))*(4/(T+1)+1/(T+1)**2)
            assert low+high<=full+uv+1e-10
            force=-Q*np.pi/(24*T*T)+low+high
            uniform=-Q*(np.pi*(kappa-2/Q)/24-10/(np.pi*np.exp(4)))
            assert force<=uniform+1e-10 and uniform<0
            band_checks.append({'q':Q,'T_over_d':T,'envelope_force_times_d_squared':force,'uniform_bound':uniform})
    # Exact rational certificate for q>=21 at Lambda*d=4, delta=2/q.
    log2_lower=2*sum(sp.Rational(1,3)**(2*n+1)/sp.Integer(2*n+1) for n in range(8))
    exp4_lower=sum(sp.Rational(4)**n/sp.factorial(n) for n in range(24))
    pi_lower=sp.Rational(314159,100000)
    certificate=sp.factor(pi_lower**2*(8*log2_lower-5-sp.Rational(2,21))-240/exp4_lower)
    assert certificate>0
    return {
        'status':'passed_exact_identities_and_supplementary_quadratures',
        'scope':'Exact dipole midplane and conditional finite-band Gaussian energy; no physical 4D reduction',
        'exact_identities':list(identities)+list(beta_identities),
        'gap_test_bump_norm':str(bump_norm),
        'gap_test_bump_gradient':str(bump_gradient),
        'gap_test_gradient_constant':str(gap_constant),
        'radial_quadrature_cases':cases,
        'critical_tail_and_cutoff_cases':tail_checks,
        'finite_band_envelope_cases':band_checks,
        'q21_margin_positive_rational_certificate':str(certificate),
        'band4_continuum_margin_times_d_squared_over_q':float(np.pi*kappa/24-10/(np.pi*np.exp(4))),
        'q21_reserved_two_channel_margin_times_d_squared_over_q':float(np.pi*(kappa-2/21)/24-10/(np.pi*np.exp(4))),
        'no_positive_transverse_gap':True,
        'exact_zero_mode_count':'max(q-1,0), integer q>0; full plane',
        'physical_four_dimensional_projection_derived':False,
        'physical_energy_corrections_bounded':False,
        'novelty_established':False,
        'full_numerical_error_certified':False,
        'lean_certified':False,
        'lean_invocations':0,
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    report=verify()
    report['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('radial_quadrature_cases','critical_tail_and_cutoff_cases','finite_band_envelope_cases')},indent=2))
