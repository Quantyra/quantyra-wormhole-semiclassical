"""Euclidean axial response with the full curl vertex and all-level error bound.

The mixed sector is exact. The massive/massive sector is bounded, not computed.
This is not the full four-dimensional Ward audit or a stability certificate.
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
from scipy.special import expit
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from casimir_cylinder_twist import angular_data, heavy_bound, mixed_loop


def coefficients(p, mass, omega, wave, c, curl):
    c0, c1 = c
    n2 = 2*(c1*c1-c0*c0)+mass*mass*curl*curl/8
    n1 = (4*c0*c1*(2*p+omega)+2*wave*(c1*c1-c0*c0)
          -mass*mass*curl*c0+mass*mass*curl*curl*wave/8)
    n0 = (2*p*(p+omega)*(c0*c0-c1*c1)+4*p*wave*c0*c1
          +mass*mass*curl*p*c1+mass*mass*curl*curl*p*(p+omega)/8)
    return n0, n1, n2


def fractions(p, mass, omega, wave, c, curl):
    n0, n1, n2 = coefficients(p, mass, omega, wave, c, curl)
    delta = mass*mass+omega*omega+wave*wave+2*p*omega
    den = delta*delta+4*wave*wave*p*p
    A = (delta*n1-2*wave*(n0-n2*p*p))/den
    B = (delta*(n0-n2*p*p)+2*wave*p*p*n1)/den
    return A, B, n2-wave*A-B


def mixed_raw(mass, length, omega, wave, c, curl):
    """Compact-minus-line bubble, including both mixed-sector orientations."""
    assert wave != 0
    def integrand(x):
        p = x/length
        _, B, D = fractions(p, mass, omega, wave, c, curl)
        energy = math.hypot(p+omega, mass)
        low = -B*expit(-abs(x))/abs(p) if p else 0.
        high = -D*expit(-length*energy)/energy
        return (low+high)/(2*math.pi*length)
    breaks = sorted(set([-omega*length, 0.]))
    endpoints = [-np.inf]+breaks+[np.inf]
    values = [quad(integrand, a, b, epsabs=2e-14, epsrel=2e-11)
              for a, b in zip(endpoints[:-1], endpoints[1:])]
    return sum(v for v, _ in values), sum(e for _, e in values)


def thermal_tadpole(mass, length):
    def integrand(x):
        energy = math.hypot(x/length, mass)
        return expit(-length*energy)/(math.pi*length*energy)
    return quad(integrand, 0, np.inf, epsabs=1e-16, epsrel=2e-11)


def mixed_reduced(mass, length, omega, wave, c, curl):
    # The exact Pauli-vertex identity leaves a transverse derivative shift
    # and one purely massive tadpole. Odd massless and massive residues cancel.
    effective = [c[0]-wave*curl/4, c[1]+omega*curl/4]
    value, error = mixed_raw(mass, length, omega, wave, effective, 0.)
    tadpole, tadpole_error = thermal_tadpole(mass, length)
    return value-mass*mass*curl*curl*tadpole/8, error+mass*mass*curl*curl*tadpole_error/8


def response(q, radius, length, omega, wave, source):
    aa, hh = np.asarray(source[:2]), np.asarray(source[2:])
    curl = omega*hh[1]-wave*hh[0]
    weights, shifts, lambdas = angular_data(q)
    value = math.pi*(6*q/5)/(12*radius**2*length**2)*(hh[0]**2-hh[1]**2)
    error = 0.
    for weight, shift, lam in zip(weights, shifts, lambdas):
        c = aa+shift*hh/radius**2
        loop, err = mixed_reduced(math.sqrt(lam)/radius, length, omega, wave, c, curl)
        value += weight*loop
        error += weight*err
    remainder = heavy_bound(q, radius, length, wave,
                           np.linalg.norm(aa)+3*abs(curl)/(8*radius), np.linalg.norm(hh))
    return value, error, remainder


def symbolic_checks():
    p, r, o, k, M, c0, c1, F = sp.symbols('p r o k M c0 c1 F', real=True)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.diag(1, -1)
    gp, gpp, gc = sx*p+sy*r, sx*(p+o)+sy*(r+k), sx*c0+sy*c1
    ident = sp.eye(2)
    # Order angular blocks as (+,-), with D_sphere=i*M*sigma_x.
    # S_L(p) V(-Q) S_H(p+Q) V(Q): derivative signs follow this routing.
    left = (M*F*ident/4).row_join(-sp.I*gc)
    right = (-M*F*ident/4).col_join(-sp.I*gc)
    numerator_H = (-sp.I*gpp).row_join(-sp.I*M*sz).col_join(
        (-sp.I*M*sz).row_join(-sp.I*gpp))
    raw = sp.expand(sp.trace(-sp.I*gp*left*numerator_H*right))
    vv = lambda x, y: 4*(p*x+r*y)*((p+o)*x+(r+k)*y)-2*(p*(p+o)+r*(r+k))*(x*x+y*y)
    expected = vv(c0,c1)+M*M*F*(p*c1-r*c0)+M*M*F*F*(p*(p+o)+r*(r+k))/8
    checks = {}
    def zero(name, expr):
        assert sp.expand(expr) == 0, name
        checks[name] = 'exact_zero'
    zero('full_spinor_trace_with_curl', raw-expected)
    zero('angular_commutator_curl_identity', (sy-sp.I*sx*sz)[0,1])

    # Verify the general partial fraction identity before substituting sources.
    n0,n1,n2 = sp.symbols('n0 n1 n2')
    delta = M*M+o*o+k*k+2*p*o
    den = delta*delta+4*k*k*p*p
    An = delta*n1-2*k*(n0-n2*p*p)
    Bn = delta*(n0-n2*p*p)+2*k*p*p*n1
    Dn = n2*den-k*An-Bn
    zero('general_partial_fractions',
         den*(n2*r*r+n1*r+n0)-(An*r+Bn)*((r+k)**2+(p+o)**2+M*M)
         -(-An*(r+k)+Dn)*(r*r+p*p))
    difference = sp.Poly(sp.expand(expected-vv(c0-k*F/4,c1+o*F/4)),r)
    sub = {n0:difference.nth(0), n1:difference.nth(1), n2:difference.nth(2)}
    odd = F*(F*o+8*c1)/8
    zero('curl_reduction_massless_residue', Bn.subs(sub)-den*odd*p)
    zero('curl_reduction_massive_residue', Dn.subs(sub)-den*(M*M*F*F/8-odd*(p+o)))
    zero('vector_kernel_traceless', vv(1,0)+vv(0,1))

    # The pure LLL occupation sum has finite support, exactly k*L/(2*pi)
    # for a right-moving channel; the left-moving channel gives its negative.
    for mode in range(1,10):
        occupied = lambda n: sp.Integer(n < 0)
        total = sum(occupied(n)-occupied(n+mode) for n in range(-mode,0))
        assert total == mode
    checks['pure_LLL_occupation_sum_cases'] = 9
    return checks


def numerical_checks():
    matsubara = []
    for length, mass, p, omega, mode, c, curl in [
        (2.,.7,.4,.3,1,[.2,1.1],.6), (3.,1.2,-.9,.4,2,[1.,-.4],-.7),
        (5.,2.,.3,-.2,1,[.7,.8],1.3), (4.,1.7,-.4,.8,1,[0.,1.],-.2),
    ]:
        wave = 2*math.pi*mode/length
        n0,n1,n2 = coefficients(p,mass,omega,wave,c,curl)
        n = np.arange(-40000,40000)
        rr = 2*math.pi*(n+.5)/length
        func = (n2*rr*rr+n1*rr+n0)/((rr*rr+p*p)*((rr+wave)**2+(p+omega)**2+mass*mass))
        discrete = np.sum(func-n2/(rr*rr))/length+n2*length/4
        _,B,D = fractions(p,mass,omega,wave,c,curl)
        energy = math.hypot(p+omega,mass)
        continuous = B/(2*abs(p))+D/(2*energy)
        expected = -B*expit(-length*abs(p))/abs(p)-D*expit(-length*energy)/energy
        error = abs(discrete-continuous-expected)
        assert error < 3e-11
        matsubara.append(error)

    reduction_errors = []
    for length,mass,omega,mode,c,curl in [
        (2.,.7,.3,1,[.2,1.1],.6), (3.,1.2,.4,2,[1.,-.4],-.7),
        (16.,2.,.5,1,[1.3,.8],1.), (16.,3.,1.,2,[-.4,2.],.7),
    ]:
        wave=2*math.pi*mode/length
        raw,er=mixed_raw(mass,length,omega,wave,c,curl)
        reduced,ee=mixed_reduced(mass,length,omega,wave,c,curl)
        error=abs(raw-reduced)
        assert error < 3e-12+er+ee
        reduction_errors.append(error)
    static_errors=[]
    for mass,length in [(2.,16.),(3.,32.),(.7,2.)]:
        wave=2*math.pi/length
        value,err=mixed_raw(mass,length,0.,wave,[0.,1.],0.)
        old,olderr=mixed_loop(mass,length,wave)
        error=abs(value-old)
        assert error < 2e-12+err+olderr
        static_errors.append(error)

    samples=[]
    for q,omega in [(3,0.),(3,.5),(3,1.),(8,.5),(8,1.5)]:
        radius,length=1.,16.
        wave=2*math.pi/length
        matrix=np.zeros((4,4))
        bounds=[]
        for i in range(4):
            v,e,b=response(q,radius,length,omega,wave,np.eye(4)[i])
            matrix[i,i]=v
            bounds.append({'quadrature_error_estimate':e,'all_massive_sector_bound':b})
        for i in range(4):
            for j in range(i):
                v,_,_=response(q,radius,length,omega,wave,np.eye(4)[i]+np.eye(4)[j])
                matrix[i,j]=matrix[j,i]=(v-matrix[i,i]-matrix[j,j])/2
        assert abs(matrix[0,0]+matrix[1,1]) < 1e-11
        assert matrix[2,2] > bounds[2]['all_massive_sector_bound']+bounds[2]['quadrature_error_estimate']
        assert matrix[3,3] < -bounds[3]['all_massive_sector_bound']-bounds[3]['quadrature_error_estimate']
        samples.append({'q':q,'length_over_radius':length,'omega_E':omega,'spatial_mode':1,
                        'source_order':['a_tau','a_z','h_tau','h_z'],
                        'mixed_plus_LLL_contact_matrix':matrix.tolist(),
                        'unit_source_errors':bounds,
                        'continuum_mixed_pair_threshold':math.sqrt(q+1+wave*wave),
                        'circle_mixed_pair_threshold':math.pi/length+math.sqrt(q+1+(wave-math.pi/length)**2)})
    return {'independent_Matsubara_cases':len(matsubara), 'max_Matsubara_residual':max(matsubara),
            'raw_vs_reduced_curl_cases':len(reduction_errors), 'max_curl_reduction_residual':max(reduction_errors),
            'static_limit_cases':len(static_errors), 'max_static_limit_residual':max(static_errors),
            'response_samples':samples}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    dependency=Path(__file__).with_name('casimir_cylinder_twist.py')
    result={'status':'passed_euclidean_cylinder_response_checks',
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'dependency_sha256':{dependency.name:hashlib.sha256(dependency.read_bytes()).hexdigest()},
            'sympy_version':sp.__version__,'scipy_version':scipy.__version__,
            'symbolic_checks':symbolic_checks(),'numerical_checks':numerical_checks(),
            'scope':'Free Dirac compact-minus-line axial Euclidean quadratic response, l=2, q>=3, nonzero circle momentum',
            'mixed_response_exact':True,'entire_massive_sector_bounded':True,
            'full_four_dimensional_Ward_audit_complete':False,
            'curved_background_response_complete':False,'physical_stability_proved':False,
            'lean_certified':False}
    args.output.write_bytes((json.dumps(result,indent=2)+'\n').encode())


if __name__=='__main__':
    main()
