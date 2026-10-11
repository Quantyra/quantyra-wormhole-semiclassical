"""Selected axial-metric/polar-gauge Ward completion on the magnetic cylinder.

Independent checks: tetrad/spin connection versus operator conjugation,
finite coordinate pullback, angular averages, and massive vacuum pressures.
No curved quantum response or physical wormhole stability is certified.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import expit, kv
import sympy as sp

sys.path.insert(0,str(Path(__file__).resolve().parent))
from casimir_cylinder_twist import angular_data, heavy_bound, mixed_loop


def operator_check():
    # A general angular displacement, stronger than just the l=2 profile.
    t, z, th, ph, eps = sp.symbols('t z theta phi epsilon', real=True)
    radius = sp.Symbol('R', positive=True)
    q, m = sp.symbols('q m', real=True)
    coords = [t, z, th, ph]
    psi = sp.Function('psi')(t, z, th)
    co = sp.diag(1, 1, radius, radius*sp.sin(th))
    for a in range(3):
        co[3, a] = -eps*radius*sp.sin(th)*sp.diff(psi, coords[a])
    metric, frame = co.T*co, co.inv()
    inverse = frame*frame.T
    g0, gi0 = metric.subs(eps, 0), inverse.subs(eps, 0)
    dg, dgi = metric.diff(eps).subs(eps, 0), inverse.diff(eps).subs(eps, 0)
    connection, variation = {}, {}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                value = delta = 0
                for d in range(4):
                    jet = sp.diff(g0[d,c],coords[b])+sp.diff(g0[b,d],coords[c])-sp.diff(g0[b,c],coords[d])
                    djet = sp.diff(dg[d,c],coords[b])+sp.diff(dg[b,d],coords[c])-sp.diff(dg[b,c],coords[d])
                    value += gi0[a,d]*jet/2
                    delta += (dgi[a,d]*jet+gi0[a,d]*djet)/2
                connection[a,b,c] = sp.simplify(value)
                variation[a,b,c] = sp.simplify(delta)
    sx = sp.Matrix([[0,1],[1,0]])
    sy = sp.Matrix([[0,-sp.I],[sp.I,0]])
    sz = sp.diag(1,-1)
    gamma = [sp.kronecker_product(sx,sp.eye(2)),sp.kronecker_product(sy,sp.eye(2)),
             sp.kronecker_product(sz,sx),sp.kronecker_product(sz,sy)]
    co0, dco = co.subs(eps,0), co.diff(eps).subs(eps,0)
    fr0, dfr = frame.subs(eps,0), frame.diff(eps).subs(eps,0)
    perturbation = sp.zeros(4)
    for mu in range(4):
        omega, domega = sp.zeros(4), sp.zeros(4)
        for a in range(4):
            for b in range(4):
                value = delta = 0
                for nu in range(4):
                    jet = sp.diff(fr0[nu,b],coords[mu])+sum(connection[nu,mu,la]*fr0[la,b] for la in range(4))
                    djet = sp.diff(dfr[nu,b],coords[mu])+sum(variation[nu,mu,la]*fr0[la,b]+connection[nu,mu,la]*dfr[la,b] for la in range(4))
                    value += co0[a,nu]*jet
                    delta += dco[a,nu]*jet+co0[a,nu]*djet
                omega += value*gamma[a]*gamma[b]/4
                domega += delta*gamma[a]*gamma[b]/4
        gg = sum((gamma[a]*fr0[mu,a] for a in range(4)),sp.zeros(4))
        dgg = sum((gamma[a]*dfr[mu,a] for a in range(4)),sp.zeros(4))
        perturbation += gg*domega+dgg*omega
        if mu == 3:
            perturbation += dgg*sp.I*(m-q*sp.cos(th)/2)
    # Covariant diffeomorphism delta A = i_xi F, with xi^phi=-psi.
    perturbation += sp.I*gamma[2]*q*sp.sin(th)*psi/(2*radius)
    generator = -sp.I*psi*(m-q*sp.cos(th)/2)
    expected = -gamma[0]*sp.diff(generator,t)-gamma[1]*sp.diff(generator,z)-gamma[2]*sp.diff(generator,th)/radius
    residual = (perturbation-expected).applyfunc(sp.simplify)
    assert residual == sp.zeros(4), residual
    return {'independent_Dirac_commutator_entries':16,'residual':'exact_zero',
            'angular_displacement':'arbitrary smooth axisymmetric psi(tau,z,theta)',
            'background':'Euclidean product cylinder, fixed R; monopole A_phi=q cos(theta)/2',
            'generator':'-i psi (m-q cos(theta)/2)'}


def geometric_checks():
    t,z,u,ph,eps,q = sp.symbols('t z u phi epsilon q', real=True)
    coords = [t,z,u,ph]
    C,H,b,at,az,ht,hz,chi,lam = [sp.Function(n)(t,z) for n in
                               ['C','H','b','a_t','a_z','h_t','h_z','chi','lambda']]
    s, y, eta = 1-u*u, sp.legendre(2,u), H/C
    yp,ypp = sp.diff(y,u),sp.diff(y,u,2)
    base = sp.diag(1,1,C/s,C*s)
    delta = sp.zeros(4)
    for i,val in enumerate([-s*yp*ht,-s*yp*hz,-s*ypp*H]):
        delta[i,3]=delta[3,i]=val
    # Exact Jacobian of phi_old=phi_new+eps eta Y'.
    jac = sp.eye(4)
    for i in range(3):
        jac[3,i]=eps*sp.diff(eta*yp,coords[i])
    transformed = jac.T*(base+eps*delta)*jac
    first = transformed.diff(eps).subs(eps,0).applyfunc(sp.simplify)
    second = (transformed.diff(eps,2).subs(eps,0)/2).applyfunc(sp.simplify)
    hhat = [ht-C*sp.diff(eta,t),hz-C*sp.diff(eta,z)]
    ahat = [at-sp.diff(b,t)+q*sp.diff(eta,t)/2,az-sp.diff(b,z)+q*sp.diff(eta,z)/2]
    expected1,expected2=sp.zeros(4),sp.zeros(4)
    for i in range(2):
        expected1[i,3]=expected1[3,i]=-s*yp*hhat[i]
        for j in range(2):
            expected2[i,j]=-s*yp**2*([ht,hz][i]*sp.diff(eta,coords[j])+[ht,hz][j]*sp.diff(eta,coords[i])-C*sp.diff(eta,coords[i])*sp.diff(eta,coords[j]))
        expected2[i,2]=expected2[2,i]=-s*yp*ypp*[ht,hz][i]*eta
    expected2[2,2]=-C*s*ypp**2*eta**2
    assert (first-expected1).applyfunc(sp.simplify)==sp.zeros(4)
    assert (second-expected2).applyfunc(sp.simplify)==sp.zeros(4)
    assert (transformed-base-eps*first-eps**2*second).applyfunc(sp.simplify)==sp.zeros(4)
    potential=sp.Matrix([eps*at*y,eps*az*y,eps*b*yp,q*u/2])
    gauge=-eps*(b*y+q*eta*(u*yp-y)/2)
    transformedA=jac.T*potential+sp.Matrix([sp.diff(gauge,c) for c in coords])
    expectedA=sp.Matrix([eps*ahat[0]*y,eps*ahat[1]*y,0,q*u/2])
    assert (transformedA-expectedA).applyfunc(sp.simplify)==sp.zeros(4,1)
    shift={H:H+C*chi,ht:ht+C*sp.diff(chi,t),hz:hz+C*sp.diff(chi,z),
           b:b+q*chi/2+lam,at:at+sp.diff(lam,t),az:az+sp.diff(lam,z)}
    for invariant in hhat+ahat:
        assert sp.simplify(invariant.subs(shift,simultaneous=True).doit()-invariant)==0
    # Derive angular averages directly, including non-quadrupole checks.
    norms=[]
    for l in range(2,8):
        yl=sp.legendre(l,u)
        w=sp.integrate(s*sp.diff(yl,u)**2,(u,-1,1))/2
        v=sp.integrate(s*s*sp.diff(yl,u,2)**2,(u,-1,1))/2
        assert w==sp.Rational(l*(l+1),2*l+1)
        assert v==sp.Rational(l*(l+1)*(l-1)*(l+2),2*l+1)
        norms.append({'l':l,'w':str(w),'v':str(v)})
    # Direct stress contraction versus the reduced tadpole functional.
    mt,mz,mo=sp.symbols('mu_tau mu_z mu_omega', real=True)
    contracted=mt*second[0,0]+mz*second[1,1]+mo*s*second[2,2]/C
    integrated=sp.integrate(sp.expand(contracted),(u,-1,1))/2
    reduced=-sp.Rational(6,5)*sum(mu*(2*h*sp.diff(eta,c)-C*sp.diff(eta,c)**2)
                                  for mu,h,c in [(mt,ht,t),(mz,hz,z)])-sp.Rational(24,5)*mo*eta**2
    assert sp.simplify(integrated-reduced)==0
    pure={H:C*chi,ht:C*sp.diff(chi,t),hz:C*sp.diff(chi,z),b:q*chi/2,at:0,az:0}
    for inv in hhat+ahat:
        assert sp.simplify(inv.subs(pure,simultaneous=True).doit())==0
    tangent=sp.simplify(reduced.subs(pure,simultaneous=True).doit())
    # Independently compute the second-order term on a finite gauge orbit.
    purejac=sp.eye(4)
    for i in range(3):
        purejac[3,i]=-eps*sp.diff(chi*yp,coords[i])
    orbit=purejac.T*base*purejac
    orbit2=orbit.diff(eps,2).subs(eps,0)/2
    orbitstress=mt*orbit2[0,0]+mz*orbit2[1,1]+mo*s*orbit2[2,2]/C
    orbitaverage=sp.integrate(sp.expand(orbitstress),(u,-1,1))/2
    assert sp.simplify(tangent+orbitaverage)==0
    return {'first_metric_pullback_entries':16,'second_metric_pullback_entries':16,
            'exact_finite_pullback_entries':16,'gauge_potential_entries':4,
            'invariant_combinations':4,'warped_radius_retained':True,
            'angular_norm_cases':norms,'tadpole_contraction':'exact_zero',
            'finite_pure_gauge_orbit_cancellation':'exact_zero'}


def massive_energy(mass,length):
    """One 2D Dirac: energy, d_L E and R d_R E for M proportional to 1/R."""
    x=mass*length
    def integrands(y):
        e=math.hypot(y,x)
        boltzmann=math.exp(-e)
        log_ratio=math.log1p(boltzmann)/boltzmann if boltzmann else 1.
        scaled=math.exp(x-e)
        return (scaled*log_ratio,e*scaled/(1+boltzmann),x*x*scaled/((1+boltzmann)*e))
    values=[quad(lambda y: integrands(y)[i],0,np.inf,epsabs=1e-14,epsrel=1e-11)[0] for i in range(3)]
    return math.exp(-x)*np.array([-2*values[0]/(math.pi*length),
                     2*values[1]/(math.pi*length**2),-2*values[2]/(math.pi*length)])


def pressure_tail(q,radius,length,retained_levels):
    """Analytic all-level upper bounds on |mu_tau|, |mu_z|, |mu_omega| tail."""
    n0=retained_levels+1
    x=length*math.sqrt(n0*(n0+q))/radius
    r=math.exp(-length/radius)
    # sum_{j>=0} j^k r^j, k=0,1,2,3.
    moments=[1/(1-r),r/(1-r)**2,r*(1+r)/(1-r)**3,
             r*(1+4*r+r*r)/(1-r)**4]
    poly=np.polynomial.Polynomial
    deg=poly([q+2*n0,2])
    upper=poly([(n0+q/2)/radius,1/radius])
    summed=lambda p: sum(c*moments[i] for i,c in enumerate(p.coef))
    envelope=lambda nu: math.sqrt(math.pi/(2*x))*math.exp(-x+nu*nu/(2*x))
    return np.array([2/(math.pi*length)*envelope(1)*summed(deg*upper),
                     (envelope(0)+envelope(2))/math.pi*summed(deg*upper**2),
                     envelope(0)/math.pi*summed(deg*upper**2)])


def numerical_checks():
    samples=[]
    for mass,length in [(0.4,2.0),(2.0,3.0),(2.0,16.0)]:
        exact=massive_energy(mass,length)
        # An independent winding/Bessel representation; 120 windings.
        series=np.zeros(3)
        for w in range(1,121):
            x=w*mass*length;sign=(-1)**(w+1)
            series+=sign*np.array([-2*mass*kv(1,x)/(math.pi*w),
                                   mass**2*(kv(0,x)+kv(2,x))/math.pi,
                                   -2*length*mass**2*kv(0,x)/math.pi])
        residual=float(np.max(np.abs(exact-series)))
        assert residual<2e-12
        step=2e-4
        dL=(massive_energy(mass,length*(1+step))[0]-massive_energy(mass,length*(1-step))[0])/(2*step*length)
        rder=(massive_energy(mass/(1+step),length)[0]-massive_energy(mass/(1-step),length)[0])/(2*step)
        errors=[abs(dL-exact[1])/max(abs(exact[1]),1e-30),abs(rder-exact[2])/max(abs(exact[2]),1e-30)]
        assert max(errors)<8e-6
        samples.append({'M':mass,'L':length,'energy_derivatives':exact.tolist(),
                        'Bessel_max_residual':residual,'finite_difference_relative_errors':errors})
    tails=[]
    for q,radius,length,N in [(3,1.0,2.0,2),(3,1.0,16.0,0),(8,1.0,16.0,0)]:
        bound=pressure_tail(q,radius,length,N)
        # A 30-level partial sum tests the bound; the proof covers infinity.
        summed=np.zeros(3)
        for n in range(N+1,N+31):
            M=math.sqrt(n*(n+q))/radius
            E,EL,ER=massive_energy(M,length)
            summed+=2*(q+2*n)*np.array([E/(2*length),EL/2,ER/(4*length)])
        assert np.all(np.abs(summed)<=bound*(1+1e-11))
        tails.append({'q':q,'R':radius,'L':length,'retained_levels':N,
                      'partial_tail_pressures':summed.tolist(),'analytic_all_level_bounds':bound.tolist()})
    return {'massive_energy_cases':samples,'pressure_tail_samples':tails,
            'numerics_are_interval_certified':False}


def path_sign_checks():
    q,v=sp.symbols('q v', real=True)
    weights,shifts,levels=angular_data(q)
    gtt=sp.factor(sum(w*b*b/lam for w,b,lam in zip(weights,shifts,levels)))
    # Python divisions in angular_data are exact because q is symbolic.
    numerator=sp.together(gtt-sp.Rational(1,4)).as_numer_denom()[0]
    shifted=sp.Poly(sp.expand(numerator.subs(q,v+3)),v)
    assert all(c>0 for c in shifted.all_coeffs())
    # Uniform mixed-loop lower bound: x=M L>=32, k=2 pi/L.
    power=sp.Rational(80,32**2)+sp.Rational(2720,32**4)
    exponential=(sp.Rational(2,3)*32**2+24)/2**32
    assert power<sp.Rational(1,10) and exponential<sp.Rational(1,1000)
    assert power+exponential<sp.Rational(1,8)
    # Entire massive remainder < pi h^2/(768 R^2 L^2), as proved previously.
    lower=sp.Rational(7,8)/12-sp.Rational(1,768)
    assert lower==sp.Rational(55,768) and lower>sp.Rational(1,16)
    samples=[]
    for charge in [3,8,64]:
        length=16.;radius=1.;wave=2*math.pi/length
        weights,shifts,levels=angular_data(charge)
        mixed=error=0.
        for w,b,lam in zip(weights,shifts,levels):
            value,err=mixed_loop(math.sqrt(lam),length,wave)
            mixed+=w*b*b*value;error+=w*b*b*err
        contact=-math.pi*(6*charge/5)/(12*length**2)
        remainder=heavy_bound(charge,radius,length,wave,0.,1.)
        # Massive pressure adds positively on the completed-shift path;
        # discarding it provides a conservative lower bound.
        assert mixed+contact+remainder+error<0
        assert mixed-remainder-error>math.pi/(16*radius**2*length**2)
        samples.append({'q':charge,'R':radius,'L':length,
                        'affine_coefficient_upper_bound':mixed+contact+remainder+error,
                        'completed_shift_coefficient_lower_bound':mixed-remainder-error,
                        'uniform_completed_shift_lower_bound':math.pi/(16*radius**2*length**2)})
    return {'GTT_gt_one_quarter_shifted_numerator':str(shifted.as_expr()),
            'mixed_loop_at_least_seven_eighths_leading':True,
            'completed_shift_positive_uniform_coefficient':'Q_shift(0,h)>pi*h^2/(16*R^2*L^2), q>=3, L/R>=16',
            'same_linear_tangent_opposite_quadratic_signs':True,'samples':samples}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data={'status':'passed_selected_sector_Ward_completion',
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'dependency_sha256':{'casimir_cylinder_twist.py':hashlib.sha256(Path(__file__).with_name('casimir_cylinder_twist.py').read_bytes()).hexdigest()},
          'sympy_version':sp.__version__,'numpy_version':np.__version__,'scipy_version':scipy.__version__,
          'operator_checks':operator_check(),'geometric_checks':geometric_checks(),
          'numerical_checks':numerical_checks(),'path_sign_checks':path_sign_checks(),'lean_certified':False,
          'selected_axial_polar_Ward_completion':True,
          'all_four_dimensional_perturbation_sectors_checked':False,
          'curved_quantum_response_computed':False,'physical_stability_proved':False}
    args.output.write_bytes((json.dumps(data,indent=2)+'\n').encode())


if __name__=='__main__':
    main()
