"""Checks for constrained JT clock dynamics and the subtracted circle response.

The ordinary argument states the approximation boundaries. This script does
not certify a four-dimensional state, a physical decay law, or novelty.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
PI=np.pi
KAPPA=1+np.arcsinh(np.sqrt(3))/(2*np.sqrt(3))


def subtracted_circle(delay,omega,throat=PI):
    length=throat+delay
    a=omega*throat/2
    ratio=np.sin(a)*np.cos(a)-np.sin(a)**2/np.tan(omega*length/2)
    return -PI*throat*delay/(3*length**3)-(omega-4*PI**2/(omega*length**2))*ratio/(6*PI)


def direct_modes(delay,omega,throat=PI,nmax=30000):
    length=throat+delay
    n=np.arange(2,nmax+1,dtype=float)
    freq=2*PI*n/length
    weights=2*(n*n-1)/(3*length**2*n)*np.sin(PI*n*throat/length)**2
    value=np.sum(2*omega**2*weights/(freq*(omega**2-freq**2)))
    zz=omega*length/(2*PI)
    assert nmax>2*abs(zz)
    bound=abs(omega)**2*length/(6*PI**3)*4/(3*nmax)
    return value,bound


def dipole_g_and_log_derivative(u):
    def terms(eta):
        aa=1-eta*eta
        qq=np.hypot(aa,2*u*eta)
        rr=(aa+qq)/2
        z=aa/qq
        term0=4*eta*eta*rr*rr/(qq*qq)
        term1=4*eta*eta*rr*u*u/(qq*qq)
        ff=rr+term0-term1
        logdiff=rr*(1-z)-2*term0*z*(1-z)-term1*(1-z+2*z*z)
        return np.sqrt(ff),logdiff/(2*np.sqrt(ff))
    points=sorted({1-k*u for k in [1,10,100,1000] if 0<1-k*u<1})
    return tuple(quad(lambda eta:terms(eta)[j],0,1,points=points,
                      epsabs=1e-11,epsrel=1e-11,limit=150)[0] for j in [0,1])


def dipole_density(d,omega=1.0,nmax=60):
    # Finite partial sums suffice for comparison: each omitted summand has
    # the same uniform relative tail correction in the ordinary proof.
    terms=[]
    for n in range(2,nmax+1):
        delay=2*PI*n/omega-PI
        f=lambda u:d*dipole_g_and_log_derivative(u)[0]/u-delay
        root=brentq(f,0.5*d/delay,2.3*d/delay,xtol=1e-15,rtol=1e-13)
        gg,ugp=dipole_g_and_log_derivative(root)
        exact=(1-1/n**2)*d*gg**2/(delay**2*(gg-ugp))
        tail=(1-1/n**2)*KAPPA*d/delay**2
        terms.append((exact,tail))
    arr=np.array(terms)
    pref=np.sin(omega*PI/2)**2/3
    return dict(d=d,modes=nmax-1,actual_partial=float(pref*np.sum(arr[:,0])),
                inverse_partial=float(pref*np.sum(arr[:,1])),
                max_relative_term_difference=float(np.max(np.abs(arr[:,0]/arr[:,1]-1))))


def principal_value_tail(d,cells=100):
    """PV int_{kappa*d}^infinity Delta chi(D,1)/D^2 dD, T=pi.

    Pole subtraction is performed in complete cotangent periods. No finite
    imaginary regulator is used. The infinite remainder has an explicit bound.
    """
    start=KAPPA*d
    assert start<PI
    def gg(x):return (1-4*PI**2/(PI+x)**2)/(6*PI*x*x)
    def gp(x):return -2*gg(x)/x+8*PI**2/(6*PI*x*x*(PI+x)**3)
    total=0.0
    for n in range(1,cells+1):
        center=(2*n-1)*PI
        a=max(start,2*(n-1)*PI);b=2*n*PI
        at_center=gg(center)
        def regular(x):
            dx=x-center
            cc=(gg(x)-at_center)/np.tan(dx/2) if abs(dx)>1e-7 else 2*gp(center)
            return -PI**2/(3*x*(PI+x)**3)+cc
        value=quad(regular,a,b,points=[center],epsabs=2e-11,epsrel=2e-10,limit=180)[0]
        value+=2*at_center*np.log(abs(np.sin((b-center)/2)/np.sin((a-center)/2)))
        total+=value
    radius=2*PI*cells
    bound=(4/3)/(2*PI)**3*(cells**-3+0.5*cells**-2)
    bound+=(16*PI**2/3)/(2*PI)**5*(cells**-5+0.25*cells**-4)
    bound+=PI**2/(9*radius**3)
    return total,bound


def main():
    exact={}
    def check(name,value):
        value=s.simplify(value)
        assert value==0,(name,value)
        exact[name]=str(value)
    re,gn,ph,eps=s.symbols('r_e G Phi eps',positive=True)
    curvature=s.symbols('R2')
    radial_squared=re**2+2*gn*ph
    reduced=(radial_squared*curvature+2-2*re**2/radial_squared)/(4*gn)
    check('linear_JT_reduction',s.diff(reduced,ph).subs(ph,0)-(curvature+2/re**2)/2)
    rp=re*(1+eps);rm=re/(1+eps)
    mass=(rp+rm)/(2*gn)
    temp=(rp-rm)/(4*s.pi*rp**2)
    check('RN_thermodynamic_Schwarzian_coefficient',s.limit((mass-re/gn)/temp**2,eps,0)-2*s.pi**2*re**3/gn)
    p,b,c=s.symbols('p b c',positive=True)
    kk=(p**3-eps**2*p*p*c)/(p*p-eps**2*b*b)**s.Rational(3,2)
    check('Lorentzian_extrinsic_curvature',s.expand(s.series(kk,eps,0,3).removeO()).coeff(eps,2)+c/p-s.Rational(3,2)*(b/p)**2)
    z3=-b*b/(2*p)
    induced=(-p*p+eps**2*b*b)/(eps*p+eps**3*z3)**2
    check('fixed_induced_boundary_metric',s.limit(induced+eps**-2,eps,0))
    aa=s.symbols('A',positive=True)
    hh=s.Function('H')(p)
    lag=aa*(b*b/(p*p)-p*p)-hh
    momentum=s.diff(lag,p)-s.diff(s.diff(lag,b),p)*b-s.diff(s.diff(lag,b),b)*c
    check('cyclic_clock_charge',momentum+2*aa*(p+c/p**2-b*b/p**3)+s.diff(hh,p))
    energy=momentum*p+s.diff(lag,b)*b-lag
    check('constrained_energy',energy-momentum*p-(aa*(b/p)**2+aa*p*p+hh))
    n,zz=s.symbols('n z',positive=True)
    check('retarded_partial_fraction',(n*n-1)/(n*n*(n*n-zz*zz))-1/(zz*zz*n*n)-(1-1/zz**2)/(n*n-zz*zz))
    dd,ww,tt=s.symbols('D omega T',positive=True)
    ll=tt+dd
    ratio=s.sin(ww*tt/2)*s.cos(ww*tt/2)-s.sin(ww*tt/2)**2*s.cot(ww*ll/2)
    response=-s.pi*tt*dd/(3*ll**3)-(ww-4*s.pi**2/(ww*ll**2))*ratio/(6*s.pi)
    check('short_return_dispersive_coefficient',s.trigsimp(s.diff(response,dd).subs(dd,0)+ww**2/(12*s.pi)))
    # Geometric derivative lower bound, using the field-line relation.
    u,rr,ee=s.symbols('u r e',positive=True)
    qa=2*rr-1+ee
    ff=rr+4*ee*rr*(rr-u*u)/qa**2
    expression=ff-1-(1-rr)*(4*rr*rr-qa*qa)/qa**2
    check('dipole_integrand_lower_bound_factorization',s.factor(expression.subs(ee,rr*(1-rr)/(rr-u*u))))
    # Differentiate the exact averaged static energy; include the independent
    # dispersive inertia coefficient from (10), not only its curvature.
    delta,x,kappa=s.symbols('delta x kappa',positive=True)
    static=x*x/16-x/8+kappa*delta*x*x/(6*s.pi)*s.log(1+s.pi/(kappa*delta*x))
    force=s.diff(static,x)
    hessian=s.diff(force,x)
    equilibrium_coefficient=-8*s.limit(force.subs(x,1)/(delta*s.log(1/delta)),delta,0)
    curvature_coefficient=s.limit((hessian.subs(x,1)-s.Rational(1,8))/(delta*s.log(1/delta)),delta,0)
    inertia_coefficient=kappa/(12*s.pi)
    loaded_coefficient=4*(curvature_coefficient-inertia_coefficient)
    check('equilibrium_log_shift',equilibrium_coefficient+8*kappa/(3*s.pi))
    check('frequency_log_shift_at_loaded_clock',loaded_coefficient-kappa/s.pi)
    check('frequency_log_shift_at_unloaded_clock',equilibrium_coefficient+loaded_coefficient+5*kappa/(3*s.pi))
    mode_cases=[]
    for delay in [0.01,0.7,4.0,13.0]:
        for omega in [0.8+0.15j,1.1+0.08j]:
            direct,bound=direct_modes(delay,omega)
            closed=subtracted_circle(delay,omega)
            error=abs(direct-closed)
            assert error<=bound
            mode_cases.append(dict(D=delay,omega=[omega.real,omega.imag],error=float(error),analytic_tail_bound=float(bound)))
    # Independent isolated-pole normalization: compare the closed response
    # residue against the fermion spectral weight, including n=1 cancellation.
    pole_cases=[]
    for delay in [0.3,2.0,7.0]:
        length=PI+delay
        for nn in [1,2,3,7]:
            frequency=2*PI*nn/length
            expected=2*(nn*nn-1)*np.sin(PI*nn*PI/length)**2/(3*length**2*nn)
            estimates=[]
            for regulator in [1e-5,1e-6]:
                estimates.append(float(-regulator*subtracted_circle(delay,frequency+1j*regulator).imag))
            error=abs(estimates[-1]-expected)
            assert error<2e-10
            pole_cases.append(dict(D=delay,n=nn,weight=expected,residue=estimates[-1],error=error))
    derivative_cases=[]
    for u0 in [0.4,0.1,0.01,0.001,0.0001]:
        gg,ugp=dipole_g_and_log_derivative(u0)
        step=u0*0.001
        finite=u0*(dipole_g_and_log_derivative(u0+step)[0]-dipole_g_and_log_derivative(u0-step)[0])/(2*step)
        assert abs(ugp-finite)<5e-8
        derivative_cases.append(dict(u=u0,g=gg,u_gprime=ugp,finite_difference=finite,error=abs(ugp-finite),u2_minus_fprime=gg-ugp))
    density_cases=[dipole_density(d) for d in [0.2,0.05,0.01]]
    assert density_cases[-1]['max_relative_term_difference']<density_cases[0]['max_relative_term_difference']/10
    pv_cases=[]
    for d in [0.01,0.001,0.0001,0.00001]:
        coarse,bound0=principal_value_tail(d,60)
        fine,bound=principal_value_tail(d,120)
        assert abs(fine-coarse)<bound0+bound
        pv_cases.append(dict(d=d,integral=fine,analytic_tail_bound=bound,refinement=abs(fine-coarse),
                    renormalized_remainder=fine+np.log(1/(KAPPA*d))/(12*PI)))
    slope=(pv_cases[-1]['integral']-pv_cases[-2]['integral'])/np.log(pv_cases[-2]['d']/pv_cases[-1]['d'])
    assert abs(slope+1/(12*PI))<1e-5
    equilibrium_cases=[]
    static_numeric=s.lambdify((x,delta,kappa),static,'numpy')
    force_numeric=s.lambdify((x,delta,kappa),force,'numpy')
    for d in [0.01,0.001,0.0001,0.00001]:
        root=brentq(lambda xx:force_numeric(xx,d,KAPPA),0.5,1.1,xtol=2e-14)
        step=2e-5
        slope_check=(static_numeric(root+step,d,KAPPA)-static_numeric(root-step,d,KAPPA))/(2*step)
        assert abs(slope_check)<1e-10
        equilibrium_cases.append(dict(d=d,x=root,finite_difference_force=float(slope_check),
                    logarithmic_shift_ratio=float((1-root)/(d*np.log(1/d)))))
    # A finite quantum system checks the constrained action's energy transfer.
    energy_cases=[]
    for tolerance in [1e-8,2e-11]:
        A=1/16;B=1/8;gap=0.7;coupling=0.035
        def rhs(t,y):
            phi,vel=y[:2];psi=y[2:4]+1j*y[4:6]
            x=np.exp(phi)
            op=coupling*np.array([[0.,1.],[1.,0.]])
            ham=np.diag([0.,gap])+(x-1)*op
            force=float(np.vdot(psi,op@psi).real)
            deriv=-1j*ham@psi
            return np.r_[vel,-(2*A*x*x+x*(force-B))/(2*A),deriv.real,deriv.imag]
        ts=np.linspace(0,40,501)
        run=solve_ivp(rhs,(0,40),[0.03,0.,1.,0.,0.,0.],method='DOP853',rtol=tolerance,atol=tolerance/20,t_eval=ts)
        assert run.success
        energies=[];norms=[]
        for y in run.y.T:
            phi,vel=y[:2];psi=y[2:4]+1j*y[4:6];x=np.exp(phi)
            ham=np.diag([0.,gap])+(x-1)*coupling*np.array([[0.,1.],[1.,0.]])
            energies.append(A*vel*vel+A*x*x-B*x+np.vdot(psi,ham@psi).real)
            norms.append(np.vdot(psi,psi).real)
        drift=float(max(abs(np.array(energies)-energies[0])))
        norm_error=float(max(abs(np.array(norms)-1)))
        assert drift<2e-9 and norm_error<2e-9
        energy_cases.append(dict(rtol=tolerance,energy_drift=drift,norm_error=norm_error))
    assert energy_cases[-1]['energy_drift']<energy_cases[0]['energy_drift']/5
    out=dict(status='passed_collective_clock_response',versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=s.__version__),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        note_sha256=hashlib.sha256((ROOT/'notes/collective-clock-response.md').read_bytes()).hexdigest(),
        exact_checks=exact,complex_mode_comparisons=mode_cases,pole_residue_cases=pole_cases,dipole_derivative_cases=derivative_cases,
        pointwise_density_cases=density_cases,principal_value_cases=pv_cases,logarithmic_slope=slope,
        expected_logarithmic_slope=-1/(12*PI),equilibrium_cases=equilibrium_cases,quantum_energy_cases=energy_cases,
        constrained_JT_clock_action_verified=True,point_dipole_pointwise_absorption_limit_proved=True,
        inverse_tail_reduced_model_resonance_asymptotic_proved=True,
        four_dimensional_matching_and_state_errors_bounded=False,physical_damping_proved=False,
        novelty_established=False,lean_invocations=0)
    (ROOT/'evidence/collective-clock-response.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=out['status'],counts={k:len(out[k]) for k in ['exact_checks','complex_mode_comparisons','dipole_derivative_cases','pointwise_density_cases','principal_value_cases','quantum_energy_cases']},
          logarithmic_slope=slope,density_last=density_cases[-1],energy_cases=energy_cases),indent=2))


if __name__=='__main__':main()
