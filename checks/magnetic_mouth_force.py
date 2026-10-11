"""Direct linear Einstein-Maxwell check of a force-balanced magnetic RN mouth.

Reproduces the weak-field Ernst/C-metric response, not a binary wormhole.
All perturbations are coefficients of a common bookkeeping parameter.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def text_hash(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').encode()).hexdigest()


def tensor_checks():
    t, r, c, phi = sp.symbols('t r c phi', real=True)
    M, P, acc, B = sp.symbols('M P acc B', real=True, nonzero=True)
    coords = [t, r, c, phi]
    f = 1-2*M/r+P**2/r**2
    sphere = 1-c*c
    diagonal = [-f, 1/f, r*r/sphere, r*r*sphere]
    relative = [-2*c*(acc*r+B*P), -2*c*(acc*r+B*P),
                -2*c*(acc*r+acc*M+B*P), -2*c*(acc*r-acc*M-B*P)]
    metric, inverse = sp.diag(*diagonal), sp.diag(*[1/x for x in diagonal])
    h = sp.diag(*[g*v for g, v in zip(diagonal, relative)])
    dinverse = -inverse*h*inverse
    gamma, dgamma = {}, {}
    for a in range(4):
        for b in range(4):
            for d in range(4):
                dg = sp.diff(metric[a,d],coords[b])+sp.diff(metric[a,b],coords[d])-sp.diff(metric[b,d],coords[a])
                dh = sp.diff(h[a,d],coords[b])+sp.diff(h[a,b],coords[d])-sp.diff(h[b,d],coords[a])
                gamma[a,b,d] = sp.factor(inverse[a,a]*dg/2)
                dgamma[a,b,d] = sp.factor((inverse[a,a]*dh+dinverse[a,a]*dg)/2)
    ricci, dricci = sp.zeros(4), sp.zeros(4)
    for a in range(4):
        for b in range(a,4):
            base = linear = sp.S.Zero
            for k in range(4):
                base += sp.diff(gamma[k,a,b],coords[k])-sp.diff(gamma[k,a,k],coords[b])
                linear += sp.diff(dgamma[k,a,b],coords[k])-sp.diff(dgamma[k,a,k],coords[b])
                for j in range(4):
                    base += gamma[k,a,b]*gamma[j,k,j]-gamma[j,a,k]*gamma[k,b,j]
                    linear += (dgamma[k,a,b]*gamma[j,k,j]+gamma[k,a,b]*dgamma[j,k,j]
                               -dgamma[j,a,k]*gamma[k,b,j]-gamma[j,a,k]*dgamma[k,b,j])
            ricci[a,b] = ricci[b,a] = sp.factor(base)
            dricci[a,b] = dricci[b,a] = sp.factor(linear)
    scalar = sp.factor(sp.trace(inverse*ricci))
    dscalar = sp.factor(sp.trace(dinverse*ricci+inverse*dricci))
    einstein = ricci-metric*scalar/2
    deinstein = dricci-h*scalar/2-metric*dscalar/2
    field, dfield = sp.zeros(4), sp.zeros(4)
    field[2,3], field[3,2] = -P, P
    potential = B*(r*r+P*P)*sphere/2
    for i in [1,2]:
        dfield[i,3], dfield[3,i] = sp.diff(potential,coords[i]), -sp.diff(potential,coords[i])
    raised = inverse*field*inverse
    draised = dinverse*field*inverse+inverse*dfield*inverse+inverse*field*dinverse
    square = sum(field[i,j]*raised[i,j] for i in range(4) for j in range(4))
    dsquare = sum(dfield[i,j]*raised[i,j]+field[i,j]*draised[i,j] for i in range(4) for j in range(4))
    stress = 2*(field*inverse*field.T-metric*square/4)
    dstress = 2*(dfield*inverse*field.T+field*dinverse*field.T+field*inverse*dfield.T
                 -h*square/4-metric*dsquare/4)
    checks = {}
    def zero(name, value):
        result = sp.factor(value)
        assert result == 0, (name,str(result))
        checks[name] = 'exact_zero'
    for i in range(4):
        for j in range(i,4):
            zero(f'background_Einstein_Maxwell_{i}_{j}',einstein[i,j]-stress[i,j])
            zero(f'linear_Einstein_Maxwell_{i}_{j}',deinstein[i,j]-dstress[i,j])
    volume = r*r
    dvolume = volume*sp.trace(inverse*h)/2
    for j in range(4):
        zero(f'background_Maxwell_{j}',sum(sp.diff(volume*raised[i,j],coords[i]) for i in range(4)))
        zero(f'linear_Maxwell_{j}',sum(sp.diff(dvolume*raised[i,j]+volume*draised[i,j],coords[i]) for i in range(4)))
    zero('linear_trace_free',dscalar)
    print('direct_tensor_equations_passed',flush=True)
    # In a fixed RN metric the coupled magnetic field does not solve Maxwell.
    frozen_density = inverse*dfield*inverse
    frozen_div = sp.factor(sum(sp.diff(volume*frozen_density[i,3],coords[i]) for i in range(4)))
    zero('nonzero_frozen_metric_Maxwell_residual',frozen_div+2*B*P*P/r**2)
    eps, z = sp.symbols('eps z',real=True)
    # Circumference / (2 pi proper radius), before choosing a common period.
    north = (1+2*M*eps*acc+P*P*eps**2*acc**2)/(1-eps*B*P/2)**4
    south = (1-2*M*eps*acc+P*P*eps**2*acc**2)/(1+eps*B*P/2)**4
    zero('linear_north_axis_ratio',sp.diff(north,eps).subs(eps,0)-2*(M*acc+B*P))
    zero('linear_south_axis_ratio',sp.diff(south,eps).subs(eps,0)+2*(M*acc+B*P))
    # Exact equality of both cone factors in the extremal accelerating family.
    balanced_acc = -B/(1+(B*P/2)**2)
    north_ext = north.subs({eps:1,M:P,acc:balanced_acc})
    south_ext = south.subs({eps:1,M:P,acc:balanced_acc})
    zero('exact_extremal_axis_balance',north_ext-south_ext)
    zero('exact_extremal_period_factor',1/north_ext-(1+(B*P/2)**2)**2)
    # Static source paper's Lambda at each pole: do not import its printed
    # expanded C expressions, which do not equal the displayed metric.
    lam = 1-B*P*c+B*B*(r*r*sphere+P*P*c*c)/4
    zero('north_period_from_metric',lam.subs(c,1)**2-(1-B*P/2)**4)
    zero('south_period_from_metric',lam.subs(c,-1)**2-(1+B*P/2)**4)
    # Extremal first-order angular area and flux perturbations at fixed r.
    # Their l=1 difference is invariant under angular relabeling, but not
    # under arbitrary radial gauge transformations away from the product limit.
    area_coeff = 2*B*r
    flux_coeff = B*(r*r+P*P)/P
    zero('area_flux_mismatch',flux_coeff-area_coeff-B*(r-P)**2/P)
    zero('horizon_area_flux_agree',(flux_coeff-area_coeff).subs(r,P))
    zero('radial_field_in_balanced_frame',P*(flux_coeff-area_coeff)/r**2-B*(1-P/r)**2)
    angular_shift=-B*P*(1-c*c)
    for value in [P*P/(1-c*c),P*P*(1-c*c)]:
        lie=angular_shift*sp.diff(value,c)
        if value==P*P/(1-c*c):
            lie+=2*value*sp.diff(angular_shift,c)
        zero('horizon_angular_relabeling_'+str(value),lie-2*B*P*c*value)
    zero('horizon_flux_relabeling',-P*sp.diff(angular_shift,c)+2*B*P*P*c)
    # Coupled scalar mass-six radial equation on the Poincare AdS2 throat.
    zero('mass_six_regular_mode',sp.diff(z*z*sp.diff(z*z,z),z)-6*z*z)
    zero('mass_six_singular_mode',sp.diff(z*z*sp.diff(z**-3,z),z)-6*z**-3)
    return dict(checks=checks, count=len(checks), frozen_metric_Maxwell_residual=str(frozen_div),
                axis_balance='M acc + B P = 0 at first order',
                extremal_exact_balance='acc = -B/(1+(B P/2)^2)',
                area_flux_mismatch='B (r-P)^2/P in the displayed areal-r gauge')


def cone_checks():
    cases=[]
    for mass in [1.,1.4,2.]:
        for field in [-.02,-.001,.001,.02]:
            def pole(c,a):
                return (1+2*mass*a*c+a*a)/(1-field*c/2)**4
            def imbalance(a):
                return pole(1,a)-pole(-1,a)
            endpoints=sorted([0.,-2*field/mass])
            acc=brentq(imbalance,*endpoints,xtol=1e-15)
            period=1/pole(1,acc)
            assert abs(imbalance(acc))<1e-13
            radius=2*(mass+np.sqrt(mass*mass-1))
            assert abs(acc*radius)<.1
            def coefficients(theta):
                c=np.cos(theta)
                sn=np.sin(theta)
                om=1+acc*radius*c
                ang=1+2*mass*acc*c+acc*acc*c*c
                lam=(1-field*c/2)**2+field*field*radius*radius*ang*sn*sn/(4*om*om)
                return lam*radius/(om*np.sqrt(ang)),radius*np.sqrt(ang)*sn/(om*lam)
            errors=[]
            for cutoff in [1e-3,5e-4]:
                pair=[]
                for south in [False,True]:
                    def proper_integrand(x):
                        return coefficients(np.pi-x if south else x)[0]
                    length=quad(proper_integrand,0,cutoff,epsabs=1e-14,epsrel=1e-12)[0]
                    angle=np.pi-cutoff if south else cutoff
                    ratio=period*coefficients(angle)[1]/length
                    pair.append(abs(ratio-1))
                errors.append(max(pair))
            assert errors[-1]<1e-7 and errors[-1]<.3*errors[0]
            cases.append(dict(M=mass,P=1.,B=field,acc=acc,period_factor=period,
                              balance_residual=abs(imbalance(acc)),
                              linear_acceleration_error=abs(acc+field/mass),
                              cone_ratio_error=errors[-1],refinement_ratio=errors[-1]/errors[0]))
    return cases


def orbital_checks():
    q, d, ell, G, coupling, P=sp.symbols('q d ell G coupling P',positive=True)
    omega=2*sp.sqrt(P/d**3)
    pem=4*sp.pi/coupling**2*sp.Rational(2,3)*(q*d/2)**2*omega**4
    pgrav=sp.Rational(2,15)*P**2*d**4*omega**6/G
    einstein=sp.pi*G*q*q/coupling**2
    assert sp.simplify((pgrav/pem).subs(G,P*P*coupling**2/(sp.pi*q*q))-sp.Rational(4,5)*P/d)==0
    lifetime=sp.simplify((P*P/(G*d)/pem).subs(G,P*P*coupling**2/(sp.pi*q*q)))
    assert sp.simplify(lifetime-sp.Rational(3,32)*d**3/P**2)==0
    # Fix ell=1 as the unit, then P=1/(16 pi q^3), d=log(q)^-2.
    selected={coupling:1/q,P:1/(16*sp.pi*q**3),d:sp.log(q)**-2}
    assert sp.simplify(pem.subs(selected)-sp.log(q)**8/(24*sp.pi*q*q))==0
    assert sp.simplify(omega.subs(selected)**2-sp.log(q)**6/(4*sp.pi*q**3))==0
    assert sp.simplify((omega*d)**2*q-sp.Integer(4)*P*q/d)==0
    for power, logpower in [(sp.Rational(-3,2),3),(-1,1),(-2,8)]:
        expression=q**power*sp.log(q)**logpower
        assert sp.simplify(q*sp.diff(expression,q)/expression-power-logpower/sp.log(q))==0
    results=[]
    for flux in [1e3,1e4,1e6,1e12]:
        delta=np.log(flux)**-2
        pe=1/(16*np.pi*flux**3)
        orbital=2*np.sqrt(pe/delta**3)
        emission=np.log(flux)**8/(24*np.pi*flux**2)
        results.append(dict(q=flux,delta=delta,orbit_frequency_times_ell=orbital,
                            unscreened_drift_spacing_ratio=orbital*delta*np.sqrt(flux),
                            emitted_energy_over_one_quantum=emission,
                            emitted_energy_over_q_over_ell=emission/flux,
                            lifetime_over_ell=3*delta**3/(32*pe**2),
                            local_magnetic_metric_size=(pe/delta)**2,
                            photon_to_gravity_power_ratio=1/(.8*pe/delta)))
    assert results[-1]['orbit_frequency_times_ell']<results[0]['orbit_frequency_times_ell']
    assert results[-1]['emitted_energy_over_one_quantum']<results[0]['emitted_energy_over_one_quantum']
    return dict(samples=results, normalization='geometric P^2=pi G q^2/g^2; g=1/q; ell/P=16 pi q^3',
                orbital_frequency='Omega ell=(log q)^3/(2 sqrt(pi) q^(3/2))',
                unscreened_drift_spacing_ratio='Omega d sqrt(q)=(log q)/(2 sqrt(pi) q)',
                emission_per_crossing_over_one_quantum='P_EM ell^2=(log q)^8/(24 pi q^2)',
                lifetime_over_ell='24 pi^2 q^6/(log q)^6',
                state_adiabatic_theorem_proved=False)


def verify():
    result = dict(status='passed_magnetic_mouth_force_checks', tensor=tensor_checks(),
                  cone_cases=cone_checks(),orbital_scale_checks=orbital_checks(),
                  full_binary_metric_matched=False, quantum_open_throat_matched=False,
                  quantum_state_derived=False, physical_transport_proved=False,
                  novelty_established=False, lean_invocations=0,
                  input_sha256={p:text_hash(ROOT/p) for p in [
                      'notes/magnetic-mouth-force.md','notes/local-maxwell-stress.md',
                      'notes/logarithmic-long-return.md']},
                  versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__))
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    result=verify()
    result['script_sha256']=text_hash(Path(__file__))
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],identities=result['tensor']['count'])))
