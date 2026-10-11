"""Checks for the conditional dipole-delay comparison; no Lean execution."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
import sympy as sp


def field_speed(u, eta):
    """Arc length derivative / d in the common prolate eta coordinate."""
    if u == 1:
        return .5
    disc = np.sqrt((1-eta*eta)**2+4*u*u*eta*eta)
    mu = ((1-eta*eta)+disc)/(2*u)
    mup = 2*eta*(u-mu)/disc
    zprime = (mu+eta*mup)/2
    root = np.sqrt(mu/u-1)
    rprime = -eta*root+(1-eta*eta)*mup/(4*u*root)
    return np.hypot(zprime, rprime)


def field_length(u):
    """Length/d from prolate coordinates; u=nu/2, point dipole, d=1."""
    return 2*quad(lambda eta:field_speed(u,eta),0,1,epsabs=2e-10,epsrel=2e-10)[0]


def energy_force(ret, T):
    q = len(ret(1.0))
    def terms(xi):
        X = np.exp(-xi*T)*ret(xi)
        _, logabs = np.linalg.slogdet(np.eye(q)+X)
        tr = np.trace(np.linalg.solve(np.eye(q)+X, X)).real
        return -2*logabs/np.pi, 2*xi*tr/np.pi
    return tuple(quad(lambda x: terms(x)[i], 0, np.inf,
                      epsabs=2e-10, epsrel=2e-10)[0] for i in (0, 1))


def verify():
    t,u,x,mu,eta = sp.symbols('t u x mu eta', positive=True)
    J = 1+1/(1+t)-2*sp.log(1+t)/t
    primitive = u-1/(t*(1+t*u))-2*sp.log(1+t*u)/t
    identities = {
        'integral_primitive': sp.diff(primitive,u)-(t*u/(1+t*u))**2,
        'integral_endpoints': primitive.subs(u,1)-primitive.subs(u,0)-J,
        'logistic_convexity': sp.diff(sp.exp(x)/(1+sp.exp(x)),x,2)-sp.exp(x)*(1-sp.exp(x))/(1+sp.exp(x))**3,
        'point_dipole_relation': u*(mu*mu-1)-(mu-u)*(1-eta*eta)-(u*(mu*mu-eta*eta)-mu*(1-eta*eta)),
        'gap_at_one': 1-4*J.subs(t,1)-(8*sp.log(2)-5),
        'gap_at_two': 1-4*J.subs(t,2)-(4*sp.log(3)-sp.Rational(13,3)),
        'J_integrand_derivative': sp.diff((t*u/(1+t*u))**2,t)-2*t*u*u/(1+t*u)**3,
        'local_speed_lower_bound':
            (u*mu+4*mu*mu*(1-u*mu)*(mu-u)**2/(mu*mu-2*u*mu+1)**2-1)
            -(1-u*mu)*(mu*mu-1)*(3*mu*mu-4*u*mu+1)/(mu*mu-2*u*mu+1)**2,
    }
    for k,v in identities.items():
        assert sp.simplify(v)==0,k
    log3_lower=2*sum(sp.Rational(1,2)**(2*n+1)/sp.Integer(2*n+1) for n in range(5))
    log2_lower=2*sum(sp.Rational(1,3)**(2*n+1)/sp.Integer(2*n+1) for n in range(2))
    assert 4*log3_lower-sp.Rational(13,3)-sp.Rational(16,270)>0
    assert 8*log2_lower-5>0
    jf = lambda a: 1+1/(1+a)-2*np.log1p(a)/a
    root = brentq(lambda a:4*jf(a)-1,2,3,xtol=1e-13)
    integral_checks=[]
    for a in [.2,1,2,root,5]:
        integral = quad(lambda z:(a*z/(1+a*z))**2,0,1,epsabs=1e-13)[0]
        assert abs(jf(a)-integral)<2e-14
        integral_checks.append({'t':a,'residual':abs(jf(a)-integral)})
    # Fixed passive return families with the same unequal half delays.
    q=5
    delays=1/((np.arange(q)+.5)/q)
    rng=np.random.default_rng(281109)
    unitary=lambda: np.linalg.qr(rng.normal(size=(q,q))+1j*rng.normal(size=(q,q)))[0]
    U,V=unitary(),unitary()
    poles=np.array([.2,.4,.9,1.4,2.2])
    cores={
        'ideal':lambda xi:np.eye(q),
        'unitary_mixing':lambda xi:U@V,
        'lossy_mixing':lambda xi:U@np.diag([1,.9,.6,.3,.1])@V,
        'resonant_mixing':lambda xi:U@np.diag((xi-poles)/(xi+poles))@V,
    }
    permutation=np.zeros((q,q))
    # One three-cycle and one two-cycle; both have antiperiodic total gluing.
    for i,j,v in [(0,1,1),(1,2,1),(2,0,1),(3,4,1),(4,3,-1)]:
        permutation[i,j]=v
    cores['signed_cycles']=lambda xi:permutation
    def ret(core):
        def result(xi):
            H=np.diag(np.exp(-xi*delays/2))
            return H@core(xi)@H
        return result
    pointwise=0; products=0; integrated=[]; variations=[]
    for name,core in cores.items():
        R=ret(core)
        for T in [.4,1,2.2]:
            for xi in [.001,.05,.3,1,4]:
                X=np.exp(-xi*T)*R(xi)
                eig=np.sort(np.abs(np.linalg.eigvals(X)))[::-1]
                envelope=np.sort(np.exp(-xi*(T+delays)))[::-1]
                assert np.all(np.cumprod(eig)<=np.cumprod(envelope)+3e-13)
                products+=q
                trace=np.trace(np.linalg.solve(np.eye(q)+X,X)).real
                bound=np.sum(envelope/(1+envelope))
                assert trace<=bound+2e-12
                pointwise+=1
            energy,force=energy_force(R,T)
            bound=np.pi*np.sum(1/(T+delays)**2)/6
            assert force<=bound+2e-10
            h=2e-4
            fp=(energy_force(R,T+h)[0]-energy_force(R,T-h)[0])/(2*h)
            residual=abs(force-fp)
            assert residual<2e-7
            if name=='ideal':
                assert abs(force-bound)<2e-11
                assert abs(energy+np.pi*np.sum(1/(T+delays))/6)<2e-11
            if name=='signed_cycles':
                lengths=[3*T+sum(delays[:3]),2*T+sum(delays[3:])]
                # Independent real massless-circle spectra after unfolding cycles.
                expected_e=-np.pi*sum(1/a for a in lengths)/6
                expected_f=np.pi*(3/lengths[0]**2+2/lengths[1]**2)/6
                assert abs(energy-expected_e)<2e-11
                assert abs(force-expected_f)<2e-11
            integrated.append({'map':name,'T':T,'energy':energy,'derivative':force,'upper_bound':bound})
            variations.append({'map':name,'T':T,'residual':residual})
    lengths=[{'u':float(a),'length_over_d':field_length(float(a)),'triangle_lower_bound':1/a}
             for a in [.025,.1,.25,.5,.75,.9,.99,1.]]
    for row in lengths: assert row['length_over_d']+1e-10>=row['triangle_lower_bound']
    speed_cases=0
    for a in [.025,.1,.25,.5,.75,.9,.99,1.]:
        for z in np.linspace(-1,1,21):
            assert field_speed(a,z)+2e-11>=1/(2*a)
            speed_cases+=1
    # Distributed Hermitian mixing. This does not assume a midpoint factorization.
    transport_cases=[]
    for n in [4,6]:
        labels=(np.arange(n)+.5)/n
        Z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        K0=(Z+Z.conj().T)/5
        Z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        K1=(Z+Z.conj().T)/7
        for xi in [.07,.3,.7]:
            def rhs(eta,flat):
                rates=np.array([field_speed(a,eta) for a in labels])
                B=np.diag(rates)
                K=K0*np.cos(2*eta)+K1*np.sin(3*eta)
                return ((-xi*B-1j*K)@flat.reshape(n,n)).ravel()
            sol=solve_ivp(rhs,(-1,1),np.eye(n,dtype=complex).ravel(),rtol=1e-10,atol=2e-12)
            assert sol.success
            R=sol.y[:,-1].reshape(n,n)
            sv=np.linalg.svd(R,compute_uv=False)
            lhs=np.cumprod(sv); bound=np.exp(-xi*np.cumsum(np.sort(1/labels)))
            assert np.all(lhs<=bound*(1+2e-8))
            for T in [.4,1,2]:
                X=np.exp(-xi*T)*R
                trace=np.trace(np.linalg.solve(np.eye(n)+X,X)).real
                envelope=1/(np.exp(xi*(T+1/labels))+1)
                assert trace<=sum(envelope)+2e-10
            transport_cases.append({'q':n,'xi':xi,'maximum_product_ratio':float(max(lhs/bound)),
                                    'rhs_evaluations':sol.nfev,'local_mixing_commutator':float(np.linalg.norm(K0@K1-K1@K0))})
    # Independent reproduction of the established MMP numerical limit.
    reproductions=[]
    for order in [32,64,128]:
        nodes,weights=leggauss(order); nodes=(nodes+1)/2; weights=weights/2
        fs=np.array([field_length(a) for a in nodes])
        threshold=brentq(lambda a:4*np.sum(weights*(a/(a+fs))**2)-1,2,3,xtol=1e-12)
        reproductions.append({'outer_order':order,'MMP_large_separation_ratio':threshold})
    assert abs(reproductions[-1]['MMP_large_separation_ratio']-reproductions[-2]['MMP_large_separation_ratio'])<2e-8
    assert 2.3<reproductions[-1]['MMP_large_separation_ratio']<2.4
    discretizations=[]
    for n in [1,2,5,15,32,128]:
        quantiles=(np.arange(n)+.5)/n
        delta=1/(2*n)
        for a in [1,2,root]:
            value=np.mean((a*quantiles/(1+a*quantiles))**2)
            error=abs(value-jf(a)); bound=delta*(a/(1+a))**2
            assert error<=bound+1e-14
            discretizations.append({'q':n,'t':a,'quadrature_error':error,'discrepancy_bound':bound})
    kappa=8*np.log(2)-5; kappa2=4*np.log(3)-13/3
    assert kappa2-16/(9*30)>0
    sample=ret(cores['resonant_mixing'])(.3)
    nonnormal=float(np.linalg.norm(sample@sample.conj().T-sample.conj().T@sample))
    assert nonnormal>.01
    # Passivity and a common front alone do not imply unequal-delay control.
    xi=.3
    front_only_product=float(np.exp(-xi*q))
    geometric_product=float(np.exp(-xi*sum(delays)))
    assert front_only_product>geometric_product
    return {'symbolic':{k:'exact_zero' for k in identities},
            'exact_positive_margin_certificates':{
                'kappa_lower_bound':str(8*log2_lower-5),
                'q15_ratio_two_bracket_lower_bound':str(4*log3_lower-sp.Rational(13,3)-sp.Rational(16,270))},
            'integral_checks':integral_checks,
            'pointwise_matrix_cases':pointwise,'exterior_power_product_cases':products,
            'integrated_cases':integrated,'energy_variations':variations,
            'nonnormal_commutator_norm':nonnormal,'field_line_lengths':lengths,
            'local_speed_cases':speed_cases,'distributed_transport_cases':transport_cases,
            'geometric_threshold_ratio':root,'MMP_reproductions':reproductions,
            'uniform_shortcut_gap_coefficient':kappa,'ratio_two_gap_coefficient':kappa2,
            'finite_quantile_checks':discretizations,
            'front_only_counterexample':{'product':front_only_product,'required_product_bound':geometric_product}}


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True); args=ap.parse_args()
    result={'status':'passed_conditional_model_checks','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'model':'Leading MMP Gaussian energy with exterior-power delay control from factored networks or specified local chiral propagation; point-dipole continuum for geometry',
            'checks':verify(),'ordinary_conditional_delay_bound':True,
            'actual_MMP_return_factorization_derived':False,'actual_four_dimensional_projection_derived':False,
            'complete_four_dimensional_energy_derived':False,
            'actual_omitted_derivatives_bounded':False,'physical_stability_proved':False,
            'novelty_established':False,'full_numerical_error_certified':False,'lean_certified':False,'lean_invocations':0}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'sha256':result['script_sha256'],
                      'geometric_ratio':result['checks']['geometric_threshold_ratio'],
                      'MMP_ratio':result['checks']['MMP_reproductions'][-1],
                      'max_variation_error':max(x['residual'] for x in result['checks']['energy_variations'])}))
