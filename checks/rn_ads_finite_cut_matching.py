"""S134.1 static RN-AdS contacts, exterior delay and mass envelope."""
import hashlib
import json
from pathlib import Path
import platform
import sys

import mpmath as mp
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'evidence/rn-ads-finite-cut-matching.json'
DPS=60
QUAD_TOL='1e-30'
GEOMETRY_TOL='1e-40'
BOUND_ALLOWANCE='1e-30'


def zero(value):
    if s.simplify(value)!=0:
        raise AssertionError(str(value))


def symbolic():
    r,b,L4,X,delta=s.symbols('r b L4 X delta',positive=True)
    C=1+6*X
    f=1+X*r**2/b**2-2*b*(1+2*X)/r+b**2*(1+3*X)/r**2
    J=C+4*X*delta+X*delta**2
    fs=s.factor(f.subs(r,b*(1+delta)))
    zero(fs-delta**2*J/(1+delta)**2)
    zero(f.subs(r,b))
    zero(s.diff(f,r).subs(r,b))
    zero(s.diff(f,r,2).subs(r,b)/2-C/b**2)
    # Differentiate the metric to derive its static-surface Christoffels.
    U,angle,azimuth=s.symbols('U angle azimuth',real=True)
    ffun=s.Function('f')(r)
    metric=s.diag(-ffun,1/ffun,r**2,r**2*s.sin(angle)**2)
    inverse=metric.inv()
    coordinates=(U,r,angle,azimuth)
    def christoffel(a,i,j):
        return s.simplify(sum(inverse[a,k]*(s.diff(metric[k,j],coordinates[i])+s.diff(metric[k,i],coordinates[j])-s.diff(metric[i,j],coordinates[k]))/2 for k in range(4)))
    normal_cov=1/s.sqrt(ffun)
    Lu=-christoffel(1,0,0)*normal_cov
    Ltheta=-christoffel(1,2,2)*normal_cov
    kt=s.simplify(Lu/metric[0,0])
    ks=s.simplify(Ltheta/metric[2,2])
    zero(kt-s.diff(ffun,r)/(2*s.sqrt(ffun)))
    zero(ks-s.sqrt(ffun)/r)
    lhat_t=kt-(kt+2*ks)/3
    zero(lhat_t-s.Rational(2,3)*(kt-ks))
    area_volume=4*s.pi*r**2*s.sqrt(ffun)
    source= -area_volume*(-1/ffun)/2
    curvature=area_volume*lhat_t*(-1/ffun)
    A=r**2/s.sqrt(ffun)
    B=r**2*s.diff(ffun,r)/(2*ffun)-r
    zero(source-2*s.pi*A)
    zero(curvature+8*s.pi*B/3)
    zero(r**2*(kt+2*ks)/s.sqrt(ffun)-(B+3*r))
    # Normalize contact coefficients by their unequal-radius product factors.
    normalized_A=s.simplify((r**2/s.sqrt(f)).subs(r,b*(1+delta))/(b**2/(s.sqrt(C)*delta)))
    target_A=(1+delta)**3/s.sqrt(J/C)
    zero(normalized_A-target_A)
    normalized_B=s.simplify((r**2*s.diff(f,r)/(2*f)-r).subs(r,b*(1+delta))/(b/delta))
    target_B=1-delta**2+X*delta*(1+delta)**2*(delta+2)/J
    zero(normalized_B-target_B)
    # Nonnegative polynomial differences prove the denominator comparisons.
    zero(C*(1+delta)**2-J-(2+8*X)*delta-(1+5*X)*delta**2)
    zero(J-X*(1+delta)**2-(1+5*X)-2*X*delta)
    # r*f_ext is increasing: its derivative numerator has positive coefficients.
    monotone=s.factor(s.diff(delta**2*J/(1+delta),delta))
    numerator=s.Poly(s.together(monotone).as_numer_denom()[0],delta)
    if any(s.simplify(coefficient).is_positive is not True for coefficient in numerator.all_coeffs() if coefficient!=0):
        raise AssertionError('monotonic denominator polynomial not positive')
    return {'factorization':str(fs),'AdS2_radius':'b/sqrt(1+6*X)','sphere_radius':'b','principal_curvatures':[str(kt),str(ks)],'contact_A_normalized':str(target_A),'contact_B_normalized':str(target_B),'monotone_denominator_derivative':str(monotone)}


def candidate(k):
    k=mp.mpf(k)
    X=(mp.sqrt(1+12*mp.pi*k**6)-1)/6
    C=1+6*X
    b=k**5*mp.sqrt(X) # Planck units only for the synthetic numeric report.
    z0=mp.mpf(7)/25
    phi=b**3/C
    kap=2*k**4*(1/k)*z0/phi
    lam=kap*b/C
    return X,C,b,kap,lam,k**(-4)


def numerical():
    mp.mp.dps=DPS
    tol,geom_tol,allowance=map(mp.mpf,(QUAD_TOL,GEOMETRY_TOL,BOUND_ALLOWANCE))
    entries=[('reference-X1',mp.mpf(1),mp.mpf(1)/10,None),('reference-X100',mp.mpf(100),mp.mpf(1)/100,None)]
    for kval in (2000,4000):
        X,C,b,kap,lam,delta=candidate(kval)
        entries.append(('candidate-k'+str(kval),X,delta,(kval,b,kap,lam)))
    results=[]
    max_quad=mp.mpf(0)
    max_geometry=mp.mpf(0)
    for name,X,delta,parameters in entries:
        C=1+6*X
        J=lambda value:C+4*X*value+X*value**2
        # I=C*delta*T_ext/b is dimensionless and close to one at a deep cut.
        # Variable w=delta/s integrates from conformal infinity to the cut.
        def w_integrand(w,signed_mass=mp.mpf(0)):
            if w==0:
                return C/X
            depth=delta/w
            D=signed_mass*delta**2*J(delta)/(1+delta)
            relative=D*(1+depth)/(depth**2*J(depth))
            return C*(1+depth)**2/J(depth)/(1-relative)
        splits=[0,delta,mp.sqrt(delta),1]
        dimensionless=mp.quad(w_integrand,splits)
        # Independent direct emblackening function at s=delta*t.
        def t_integrand(t,signed_mass=mp.mpf(0)):
            depth=delta*t
            D=signed_mass*delta**2*J(delta)/(1+delta)
            f=1+X*(1+depth)**2-2*(1+2*X)/(1+depth)+(1+3*X)/(1+depth)**2-D/(1+depth)
            # Avoid cancellation of the double zero in the large-k candidate.
            # Its independently expanded polynomial is exactly r^2*f.
            poly=depth**2*(C+4*X*depth+X*depth**2)-D*(1+depth)
            f_stable=poly/(1+depth)**2
            if abs(f-f_stable)>mp.mpf('1e-40')*max(1,abs(f_stable)):
                raise AssertionError('direct RN and factorization disagree numerically')
            return C*delta**2/f_stable
        direct=mp.quad(t_integrand,[1,2,10,1/mp.sqrt(delta),1/delta,mp.inf])
        error=abs(dimensionless-direct)
        max_quad=max(max_quad,error)
        if error>tol:
            raise AssertionError('independent exterior quadratures disagree')
        upper=1+delta*(2*mp.log(1/delta)-delta+C/X)
        if dimensionless<1-allowance or dimensionless>upper+allowance:
            raise AssertionError('analytic travel bound failed')
        # Direct finite-cut curvature calculation in r units with b=1.
        r=1+delta
        f=delta**2*J(delta)/r**2
        fp=2*X*r+2*(1+2*X)/r**2-2*(1+3*X)/r**3
        kt=fp/(2*mp.sqrt(f))
        ks=mp.sqrt(f)/r
        Bdirect=r**2*(kt-ks)/mp.sqrt(f)
        Anorm=(r**2/mp.sqrt(f))/(1/(mp.sqrt(C)*delta))
        Bnorm=Bdirect*delta
        Aclosed=r**3/mp.sqrt(J(delta)/C)
        Bclosed=1-delta**2+X*delta*r**2*(delta+2)/J(delta)
        geometry_error=max(abs(Anorm-Aclosed),abs(Bnorm-Bclosed))
        max_geometry=max(max_geometry,geometry_error)
        if geometry_error>geom_tol:
            raise AssertionError('independent finite-cut curvature disagrees')
        if abs(Anorm-1)>7*delta+allowance or abs(Bnorm-1)>3*delta+allowance:
            raise AssertionError('contact matching error bound failed')
        mass_results=[]
        for sign in (-1,1):
            signed_mass=sign*mp.mpf(1)/100
            shifted=mp.quad(lambda w:w_integrand(w,signed_mass),splits)
            if not dimensionless/(1+abs(signed_mass))-allowance<=shifted<=dimensionless/(1-abs(signed_mass))+allowance:
                raise AssertionError('signed mass delay envelope failed')
            mass_results.append({'signed_relative_mass_at_cut':str(signed_mass),'normalized_delay':str(shifted)})
        record={'case':name,'X':str(X),'delta':str(delta),'normalized_delay_CdeltaT_over_b':str(dimensionless),'normalized_delay_upper':str(upper),'quadrature_error':str(error),'normalized_contact_A':str(Anorm),'normalized_contact_B':str(Bnorm),'geometry_error':str(geometry_error),'synthetic_signed_mass_cases':mass_results}
        if parameters:
            kval,b,kap,lam=parameters
            ep=lam/delta
            clock_delay=ep*dimensionless
            clock_upper=ep*upper
            seam_excess=2*(clock_delay-mp.atan(ep))
            seam_upper=2*(lam*(2*mp.log(1/delta)-delta+C/X)+ep**3/3)
            if seam_excess<-allowance or seam_excess>seam_upper+mp.mpf('1e-55'):
                raise AssertionError('leading-piece seam bound failed')
            # Exact thermal horizon equation, rescaled by lambda to avoid subtraction.
            equation=lambda w:(1+lam*w)**3-w-lam*(1+18*X)/(2*C)*w**2-lam**2*6*X/C*w**3-lam**3*3*X/(2*C)*w**4
            horizon_w=mp.findroot(equation,(mp.mpf('0.9'),mp.mpf('1.1')))
            sh=lam*horizon_w
            D=sh**2*J(sh)/(1+sh)
            thermal_m=D*(1+delta)/(delta**2*J(delta))
            thermal_bound=256*ep**2
            if sh>8*lam or thermal_m>thermal_bound or thermal_bound>=mp.mpf('2e-38'):
                raise AssertionError('initial thermal mass envelope failed')
            if clock_upper>=mp.mpf('8.76e-21') or seam_upper>=mp.mpf('7.5e-32'):
                raise AssertionError('candidate uniform upper bound failed')
            record['candidate']={'k':kval,'C':str(C),'sphere_to_AdS2_radius_ratio':str(mp.sqrt(C)),'b_Planck_units':str(b),'kappa0':str(kap),'lambda0':str(lam),'epsilon_cut':str(ep),'initial_clock_delay':str(clock_delay),'initial_clock_delay_upper':str(clock_upper),'seam_model_excess':str(seam_excess),'seam_model_excess_upper':str(seam_upper),'initial_thermal_horizon_s':str(sh),'initial_thermal_mass_relative_at_cut':str(thermal_m),'initial_thermal_mass_upper':str(thermal_bound)}
        results.append(record)
    # Uniform k>=2000 estimates reduce to elementary monotone functions.
    k0=mp.mpf(2000)
    correction=9*mp.log(k0)/k0**4
    delay_uniform=(mp.mpf(14)/25)/k0**6*(1+correction)
    seam_uniform=2*(mp.mpf(14)/25)/k0**10*9*mp.log(k0)+mp.mpf(2)/3*((mp.mpf(14)/25)/k0**6)**3
    if delay_uniform>=mp.mpf('8.76e-21') or seam_uniform>=mp.mpf('7.5e-32'):
        raise AssertionError('uniform all-k elementary bound failed')
    return {'cases':results,'max_quadrature_error':str(max_quad),'max_geometry_error':str(max_geometry),'uniform_k_ge_2000':{'clock_delay_upper':str(delay_uniform),'seam_model_excess_upper':str(seam_uniform),'A_relative_error_upper':'4.375e-13','B_relative_error_upper':'1.875e-13','thermal_mass_relative_upper':'2e-38','proof_scope':'fixed initial static reference; delta=k^-4 and z0<=7/25'}}


def main():
    identities=symbolic()
    checks=numerical()
    paths=['notes/rn-ads-finite-cut-matching.md','checks/rn_ads_finite_cut_matching.py','notes/cutoff-clock-contact-matching.md','notes/axisymmetric-axial-action.md','notes/boundary-vertex-matching.md','notes/causal-compact-response.md','notes/source-clock-joint.md']
    receipt={'status':'pass','milestone':'E003/S134.1','symbolic':identities,'numerical':checks,'predeclared':{'decimal_digits':DPS,'quadrature_absolute_tolerance':QUAD_TOL,'geometry_absolute_tolerance':GEOMETRY_TOL,'bound_rounding_allowance':BOUND_ALLOWANCE},'physical_scope':{'reference_geometry':'static RN-AdS4; unequal AdS2 and sphere radii','initial_clock_only':True,'mass_envelope_conditional':True,'seam_model_is_full_geometry':False,'moving_clock_map_supplied':False,'physical_source_coefficients_assigned':False,'S074_complete':False,'S076_margin_inferred':False,'physical_opening_claimed':False},'primary_source':{'url':'https://arxiv.org/html/2308.00871v2','equations':['2.18-2.27','2.34-2.43'],'use':'geometry and static coordinate dictionary; finite-cut bounds derived here'},'inputs_sha256':{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},'runtime':{'python':sys.version,'platform':platform.platform(),'sympy':s.__version__,'mpmath':mp.__version__}}
    OUTPUT.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('S134.1 RN-AdS finite-cut matching: PASS')
    print('Exact unequal-radius contacts, exterior delay, seam and conditional mass bounds verified.')
    print('Maximum independent quadrature error: '+checks['max_quadrature_error'])
    print('Maximum independent curvature error: '+checks['max_geometry_error'])
    print('Static initial-clock comparison only; evolving matching and full S074 remain open.')
    print('Receipt: '+str(OUTPUT))


if __name__=='__main__':
    main()
