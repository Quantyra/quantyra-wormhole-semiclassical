"""Cross-check the conditional passive-return travel-time obstruction.

The note proves the all-frequency bound. Finite matrix samples do not prove
passivity of an actual exterior or certify every numerical error.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp


def loop_energy_force(return_map, T):
    q = len(return_map(1.0))
    def integrands(xi):
        X = np.exp(-xi*T)*return_map(xi)
        _, logabs = np.linalg.slogdet(np.eye(q)+X)
        force = np.trace(np.linalg.solve(np.eye(q)+X, X)).real
        return -2*logabs/np.pi, 2*xi*force/np.pi
    energy, error_e = quad(lambda x: integrands(x)[0], 0, np.inf, epsabs=2e-10, epsrel=2e-10)
    force, error_f = quad(lambda x: integrands(x)[1], 0, np.inf, epsabs=2e-10, epsrel=2e-10)
    return energy, force, error_e, error_f


def resonant_real_spectrum(length, pole):
    """Independent real roots for R(s)=exp(-D*s)*(s-a)/(s+a).

    Positive frequencies solve w*L+2*atan(w/a)=2*pi*n. The
    infinite-line plus isolated-resonance density uses the same regulator.
    Differentiate the real roots implicitly to obtain the vacuum force.
    """
    epsilons = [0.12, 0.06, 0.03, 0.015]
    cutoff = 45/min(epsilons)
    count = int(np.ceil((cutoff*length+np.pi)/(2*np.pi)))+2
    roots = np.array([brentq(lambda w: w*length+2*np.arctan(w/pole)-2*np.pi*n,
                            (2*np.pi*n-np.pi)/length, 2*np.pi*n/length,
                            xtol=2e-12) for n in range(1, count+1)])
    forces = []
    for eps in epsilons:
        slope = -roots/(length+2*pole/(pole*pole+roots*roots))
        value = -2*np.sum(slope*(1-eps*roots)*np.exp(-eps*roots))+1/(np.pi*eps*eps)
        forces.append(float(value))
    extrap = [(forces[i]-20*forces[i+1]+64*forces[i+2])/45 for i in range(2)]
    return {'positive_roots': count, 'regulator_cutoffs': epsilons,
            'regulated_forces': forces, 'extrapolated_forces': extrap}


def verify():
    r,t,c,T,D,A,q,x = sp.symbols('r t c T D A q x', real=True)
    gap = 2*A/T**3+q*sp.pi/(24*T**2)-q*sp.pi/(6*(T+D)**2)
    f = 4*T**3/(T+D)**2-T
    identities = {
        'disk_bound_numerator': r+(r-1)*t*c-t*t-((r-t)*(1+t)+(1-r)*t*(1-c)),
        'short_throat_gap': gap-(2*A/T**3+q*sp.pi*(D-T)*(D+3*T)/(24*T**2*(T+D)**2)),
        'stationary_threshold': (q*sp.pi/(6*(T+D)**2)-q*sp.pi/(24*T**2)-2*A/T**3)*24*T**3/(q*sp.pi)-(f-48*A/(q*sp.pi)),
        'threshold_monotonicity': (sp.diff(f,T)*(T+D)**3).subs(T,D+x)-(8*D**3+24*D**2*x+18*D*x*x+3*x**3),
        'MMP_shortest_return': (sp.pi**2*A/T**2+q*sp.pi/(24*T)-q*sp.pi/(6*(T+D))).subs(T,sp.pi*x) - (A/x**2+q/(24*x)-q*sp.pi/(6*(sp.pi*x+D))),
    }
    for name,value in identities.items(): assert sp.simplify(value)==0,name
    fermi, fermi_error = quad(lambda u: u*np.exp(-u)/(1+np.exp(-u)),0,np.inf)
    assert abs(fermi-np.pi**2/12)<2e-12

    rng=np.random.default_rng(20261011)
    def unitary(n):
        return np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]
    dimension=4; delay=0.7; area_energy=0.13
    U,V=unitary(dimension),unitary(dimension)
    extra=np.array([0.0,0.2,0.5,1.1]); poles=np.array([0.4,0.8,1.5,2.5])
    P=np.roll(np.eye(dimension),1,axis=0)
    maps={
        'ideal': lambda s: np.exp(-delay*s)*np.eye(dimension),
        'periodic': lambda s: -np.exp(-delay*s)*np.eye(dimension),
        'cyclic': lambda s: np.exp(-delay*s)*P,
        'unequal_delays_mixed': lambda s: np.exp(-delay*s)*U@np.diag(np.exp(-s*extra))@V,
        'resonant_mixed': lambda s: np.exp(-delay*s)*U@np.diag((s-poles)/(s+poles))@V,
        'lossy_mixed': lambda s: np.exp(-delay*s)*U@np.diag(poles/(s+poles))@V,
    }
    pointwise=[]; integrated=[]; variations=[]
    for name,fn in maps.items():
        for tt in [0.3,delay,1.4]:
            for xi in [0.001,0.08,0.5,3.0,10.0]:
                R=fn(xi); X=np.exp(-xi*tt)*R; radius=np.exp(-xi*(tt+delay))
                norm=np.linalg.norm(R,2)
                trace=float(np.trace(np.linalg.solve(np.eye(dimension)+X,X)).real)
                upper=dimension*radius/(1+radius)
                assert norm<=np.exp(-xi*delay)*(1+3e-14)
                assert trace<=upper+3e-12
                pointwise.append({'map':name,'T':tt,'xi':xi,'trace':trace,'upper_bound':upper})
            energy,force,err_e,err_f=loop_energy_force(fn,tt)
            bound=dimension*np.pi/(6*(tt+delay)**2)
            lower=-dimension*np.pi/(3*(tt+delay)**2)
            total=-2*area_energy/tt**3-dimension*np.pi/(24*tt**2)+force
            assert lower-2e-10<=force<=bound+2e-10
            if tt<=delay: assert total < 0
            if name=='ideal':
                assert abs(energy+dimension*np.pi/(6*(tt+delay)))<1e-10
                assert abs(force-bound)<1e-10
            if name=='periodic':
                assert abs(energy-dimension*np.pi/(3*(tt+delay)))<1e-10
                assert abs(force-lower)<1e-10
            integrated.append({'map':name,'T':tt,'energy':energy,'energy_derivative':force,
                               'derivative_bound':bound,'derivative_lower_bound':lower,'total_energy_derivative':total,
                               'quadrature_error_estimates':[err_e,err_f]})
            h=2e-4
            ep=loop_energy_force(fn,tt+h)[0]; em=loop_energy_force(fn,tt-h)[0]
            diff=(ep-em)/(2*h)
            assert abs(diff-force)<3e-7
            variations.append({'map':name,'T':tt,'finite_difference_residual':abs(diff-force)})
    fn=maps['unequal_delays_mixed']
    nonnormal=np.linalg.norm(fn(0.6)@fn(0.6).conj().T-fn(0.6).conj().T@fn(0.6))
    noncommuting=np.linalg.norm(fn(0.3)@fn(1.1)-fn(1.1)@fn(0.3))
    assert nonnormal>0.01 and noncommuting>0.01

    independent=[]
    for pole,length in [(0.7,1.3),(1.8,2.2),(5.0,1.1)]:
        # D=0.4 and T=L-D. The derivative is with respect to T at fixed D,a.
        outside=0.4
        fn=lambda s: np.array([[np.exp(-outside*s)*(s-pole)/(s+pole)]])
        _, force, _, error=loop_energy_force(fn,length-outside)
        spectrum=resonant_real_spectrum(length,pole)
        residuals=[abs(value-force) for value in spectrum['extrapolated_forces']]
        assert residuals[-1]<3e-7 and residuals[-1]<residuals[0]/8
        independent.append({'pole':pole,'total_length':length,'imaginary_frequency_force':force,
                            'quadrature_error_estimate':error,**spectrum,'residuals':residuals})

    thresholds=[]
    for coeff,charge,outside in [(0.13,4,0.7),(2.0,8,1.0),(0.02,16,3.0)]:
        d0=16*coeff/(charge*np.pi)
        root=brentq(lambda length: 4*length**3/(length+outside)**2-length-3*d0,
                    outside, outside+4*d0+4*outside)
        assert root>max(outside,d0)
        thresholds.append({'A':coeff,'q':charge,'D':outside,'d0':d0,'minimum_stationary_T':root,
                           'uniform_positive_error_derivative_budget':2*coeff/outside**3})
    # Audit a generic inference used after Eq. (51) of arXiv:2505.13682v2.
    # This does not construct that paper's physical operator/eigenvector.
    signed_laplace=sp.exp(-x)-sp.Rational(3,4)*sp.exp(-2*x)
    value=sp.simplify(signed_laplace.subs(x,sp.log(sp.Rational(4,3))))
    slope=sp.simplify(sp.diff(signed_laplace,x).subs(x,sp.log(sp.Rational(4,3))))
    assert value==sp.Rational(21,64) and slope==sp.Rational(3,32)
    return {'symbolic':{name:'exact_zero' for name in identities},
            'fermi_integral_residual':abs(fermi-np.pi**2/12),
            'pointwise_matrix_bound_cases':len(pointwise),
            'integrated_cases':integrated,'energy_variations':variations,
            'nonnormal_matrix_commutator_norm':float(nonnormal),
            'different_frequency_commutator_norm':float(noncommuting),
            'independent_real_spectra':independent,'sharp_thresholds':thresholds,
            'signed_Laplace_inference_audit':{'positive_difference':str(value),'positive_derivative':str(slope),
              'generic_monotonicity_inference_refuted':True,'physical_force_theorem_refuted':False}}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    result={'status':'passed_passive_return_bound_cross_checks',
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'checks':verify(),
            'ordinary_conditional_no_shortcut_theorem':True,
            'model':'MMP leading length energy with q free complex channels and T-independent causal passive Gaussian return of front time D',
            'actual_MMP_return_map_derived':False,'complete_four_dimensional_energy_derived':False,
            'universal_wormhole_no_go':False,'physical_stability_proved':False,'novelty_established':False,
            'full_numerical_error_certified':False,'lean_certified':False,'lean_invocations':0}
    args.output.write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print(json.dumps({'status':result['status'],'script_sha256':result['script_sha256'],
                      'real_spectrum_residuals':[s['residuals'] for s in result['checks']['independent_real_spectra']]}))


if __name__=='__main__': main()
