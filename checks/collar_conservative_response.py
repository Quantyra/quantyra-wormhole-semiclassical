"""Verify the conservative collar term and its separate dynamic limit.

Independent comparators: filled-fermion spectra, CFT ground-energy map,
Virasoro response sums, and spatial/heat-kernel quadratures. No Lean run.
"""
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.optimize import brentq
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
NOTE=ROOT/'notes/collar-conservative-response.md'
OUT=ROOT/'evidence/collar-conservative-response.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient(n,T,D,h):
    n=np.asarray(n,dtype=float)
    length=T+D
    k=2*np.pi*n/length
    out=np.full(n.shape,T/length,dtype=complex)
    keep=n!=0
    out[keep]=-np.expm1(-1j*k[keep]*T)/(1j*k[keep]*length)*np.exp(-h*h*k[keep]**2/2)
    return out


def profile(T,D,h,grid=8192):
    length=T+D
    indices=np.fft.fftfreq(grid)*grid
    k=2*np.pi*indices/length
    coeff=coefficient(indices,T,D,h)
    values=(grid*np.fft.ifft(coeff)).real
    derivative=(grid*np.fft.ifft(1j*k*coeff)).real
    assert values.min()>-2e-13 and values.max()<1+2e-13
    return values,derivative,length/grid


def ground_energy(phi,T,D,h,grid=8192):
    f,fp,ds=profile(T,D,h,grid)
    a=np.expm1(phi)
    velocity=1+a*f
    assert velocity.min()>0
    optical=ds*np.sum(1/velocity)
    gradient=ds*np.sum((a*fp)**2/velocity)
    return float(-np.pi/(6*optical)-gradient/(24*np.pi))


def fermion_energy(phi,T,D,h,cutoff):
    # NS modes +/-1/2,...,+/-(cutoff-1/2). Fill the negative energies.
    # Subtract the SAME original sea expectation of H[v], then restore
    # its fixed-circle Casimir constant. Both chiralities are included.
    r=np.arange(-cutoff,cutoff)+.5
    differences=r[:,None]-r[None,:]
    vcoeff=(differences==0).astype(float)+np.expm1(phi)*coefficient(differences,T,D,h)
    one=(2*np.pi/(T+D))*(r[:,None]+r[None,:])/2*vcoeff
    assert np.max(np.abs(one-one.conj().T))<1e-11
    eigenvalues=np.linalg.eigvalsh(one)
    assert np.sum(eigenvalues<0)==cutoff
    old_sea=np.trace(one[:cutoff,:cutoff]).real
    vmean=1+np.expm1(phi)*T/(T+D)
    energy=2*(eigenvalues[:cutoff].sum()-old_sea)-vmean*np.pi/(6*(T+D))
    return float(energy)


def spectral(T,D,h):
    length=T+D
    top=max(100,int(math.ceil(length*math.sqrt(80)/(2*np.pi*h))))
    n=np.arange(2,top+1,dtype=float)
    energies=2*np.pi*n/length
    weights=(2*np.pi/length)**2*n*(n*n-1)/6*np.abs(coefficient(n,T,D,h))**2
    return energies,weights


def static_parts(T,D,h,grid=8192):
    length=T+D
    f,fp,ds=profile(T,D,h,grid)
    first=ds*f.sum()
    second=ds*np.dot(f,f)
    gradient=ds*np.dot(fp,fp)
    chi=-gradient/(12*np.pi)+np.pi*(second-first*first/length)/(3*length*length)
    contact=-np.pi*first/(6*length*length)
    energies,weights=spectral(T,D,h)
    independent=-2*np.sum(weights/energies)
    assert abs(chi-independent)<2e-11
    return dict(first=float(first),second=float(second),gradient=float(gradient),
        chi_static=float(chi),chi_spectral=float(independent),contact=float(contact),
        matter_stiffness=float(chi+contact))


def exact_checks():
    phi,F,Q,I,L,c=sy.symbols('phi F Q I L c',real=True,nonzero=True)
    optical=L-F*phi+(Q-F/2)*phi**2
    energy=-sy.pi*c/(6*optical)-c*I*phi**2/(24*sy.pi)
    want=-sy.pi*c*F/(6*L**2)+sy.pi*c*(Q-F**2/L)/(3*L**2)-c*I/(12*sy.pi)
    assert sy.simplify(sy.diff(energy,phi,2).subs(phi,0)-want)==0
    # Map Schwarzian: y'=v_bar/v.
    x=sy.symbols('x');v=sy.Function('v')(x)
    yp=1/v
    schwarz=sy.diff(yp,x,2)/yp-sy.Rational(3,2)*(sy.diff(yp,x)/yp)**2
    assert sy.simplify(v*schwarz+sy.diff(v,x,2)-sy.diff(v,x)**2/(2*v))==0
    # Scalar curvature of diag(-v(t,s)^2,1), calculated from connections.
    t,s=sy.symbols('t s');vf=sy.Function('v')(t,s)
    coords=[t,s];metric=sy.diag(-vf**2,1);inverse=metric.inv()
    gamma=[[[sy.simplify(sum(inverse[a,d]*(sy.diff(metric[d,b],coords[k])+sy.diff(metric[d,k],coords[b])-sy.diff(metric[b,k],coords[d])) for d in range(2))/2) for k in range(2)] for b in range(2)] for a in range(2)]
    ricci=sy.Matrix(2,2,lambda a,b:sum(sy.diff(gamma[k][a][b],coords[k])-sy.diff(gamma[k][a][k],coords[b])+sum(gamma[k][a][b]*gamma[d][k][d]-gamma[d][a][k]*gamma[k][b][d] for d in range(2)) for k in range(2)))
    curvature=sy.simplify(sum(inverse[a,b]*ricci[a,b] for a in range(2) for b in range(2)))
    assert sy.simplify(curvature+2*sy.diff(vf,s,2)/vf)==0
    rho=-sy.log(vf)
    raw=sy.simplify(vf*(rho*curvature-sy.diff(rho,t)**2/vf**2+sy.diff(rho,s)**2))
    integrated=-sy.diff(vf,s)**2/vf-sy.diff(vf,t)**2/vf**3
    boundary=sy.diff(2*sy.log(vf)*sy.diff(vf,s),s)
    assert sy.simplify(raw-integrated-boundary)==0
    return ['clock_contact_plus_static_hessian','integrated_schwarzian_sign',
        'lapse_metric_scalar_curvature','lorentzian_weyl_action_difference']


def ground_cases():
    rows=[]
    for T,D,h,phi in [(np.pi,.7,.09,.25),(2.1,1.3,.07,-.3),(np.pi,.4,.05,.18)]:
        mapped=ground_energy(phi,T,D,h,16384)
        grid_error=abs(mapped-ground_energy(phi,T,D,h,8192))
        finite=[]
        for cutoff in [32,64,128]:
            value=fermion_energy(phi,T,D,h,cutoff)
            finite.append(dict(cutoff=cutoff,energy=value,error=abs(value-mapped)))
        assert finite[-1]['error']<2e-8 and grid_error<2e-12
        rows.append(dict(T=T,D=D,h=h,phi=phi,ground_energy=mapped,grid_error=grid_error,fermion_refinements=finite))
    return rows


def curvature_cases():
    rows=[]
    for h in [.08,.04,.02]:
        T,D=np.pi,.4
        parts=static_parts(T,D,h)
        at_zero=ground_energy(0,T,D,h)
        values=[]
        for step in [.02,.01,.005]:
            value=(ground_energy(step,T,D,h)-2*at_zero+ground_energy(-step,T,D,h))/step**2
            values.append(value)
        richardson=(4*values[-1]-values[-2])/3
        error=abs(richardson-parts['matter_stiffness'])
        assert error<2e-8
        rows.append(dict(h=h,parts=parts,steps=[.02,.01,.005],finite_difference_curvatures=values,
            richardson_curvature=richardson,richardson_error=error))
    return rows


def mobius_cases():
    rows=[]
    length=3.7
    theta=2*np.pi*np.arange(32768)/32768
    for a in [.2,.6,.85]:
        velocity=1+a*np.cos(theta)
        derivative=-2*np.pi*a*np.sin(theta)/length
        ds=length/len(theta)
        optical=ds*np.sum(1/velocity)
        gradient=ds*np.sum(derivative*derivative/velocity)
        actual=-np.pi/(6*optical)-gradient/(24*np.pi)
        expected=-np.pi/(6*length)
        assert abs(actual-expected)<2e-13
        rows.append(dict(a=a,energy=float(actual),unchanged_energy=expected,error=abs(actual-expected)))
    return rows


def heat_cases():
    rows=[]
    for T,D,h in [(np.pi,.2,.04),(np.pi,.2,.02),(np.pi,10.,.04),(2.1,.8,.09)]:
        gradient=static_parts(T,D,h,16384)['gradient']
        m=np.arange(-16,17,dtype=float)
        image=np.sum(np.exp(-(m*(T+D))**2/(4*h*h))-np.exp(-(T+m*(T+D))**2/(4*h*h)))/(np.sqrt(np.pi)*h)
        rel=abs(gradient*np.sqrt(np.pi)*h-1)
        bound=8*np.exp(-min(T,D)**2/(4*h*h))
        assert abs(gradient-image)<3e-11
        assert rel <= bound+3e-13
        rows.append(dict(T=T,D=D,h=h,fourier_gradient=gradient,image_gradient=float(image),
            comparison_error=abs(gradient-image),relative_leading_error=rel,analytic_bound=float(bound)))
    return rows


def sharp_delta(w,T,D):
    length=T+D
    return -np.pi*T*D/(3*length**3)-(w-4*np.pi**2/(w*length**2))/(6*np.pi)*np.sin(w*T/2)*np.sin(w*D/2)/np.sin(w*length/2)


def dynamic_cases():
    rows=[]
    for D,w in [(.2,1+.03j),(2.,1.15+.025j),(20.,.9+.02j)]:
        for h in [.04,.02,.01]:
            T=np.pi
            energies,weights=spectral(T,D,h)
            smooth=np.sum(2*w*w*weights/(energies*(w*w-energies*energies)))
            sharp=sharp_delta(w,T,D)
            difference=abs(smooth-sharp)
            bound=8*abs(w)**2*h/(9*np.pi**2)+8*abs(w)**5*h*h/(3*np.pi**2*w.imag*w.real)
            assert difference<bound
            rows.append(dict(D=D,h=h,w=[w.real,w.imag],smooth=[float(smooth.real),float(smooth.imag)],
                sharp=[float(sharp.real),float(sharp.imag)],difference=float(difference),bound=float(bound)))
    return rows


def instability_cases():
    rows=[]
    for h in [.04,.02,.01]:
        T,D=np.pi,.2
        parts=static_parts(T,D,h)
        mass=1/8
        static=7/24+parts['matter_stiffness']
        assert static<0
        energies,weights=spectral(T,D,h)

        def kernel(gamma):
            return mass*gamma*gamma+static+np.sum(2*gamma*gamma*weights/(energies*(energies*energies+gamma*gamma)))

        upper=2*math.sqrt(-static/mass)+1
        gamma=brentq(kernel,0,upper,xtol=1e-13,rtol=1e-14)
        derivative=2*mass*gamma+np.sum(4*gamma*weights*energies/(energies*energies+gamma*gamma)**2)
        assert abs(kernel(gamma))<2e-12 and derivative>0
        rows.append(dict(h=h,mass=mass,static_stiffness=static,growing_rate=gamma,
            root_residual=float(kernel(gamma)),positive_derivative=float(derivative),
            gradient_completion_stiffness=parts['gradient']/(12*np.pi),
            stiffness_after_gradient_completion=static+parts['gradient']/(12*np.pi)))
    return rows


def weyl_completion_cases():
    rows=[]
    for D,h in [(.2,.04),(.2,.02),(2.,.05)]:
        T=np.pi
        parts=static_parts(T,D,h)
        mass=1/8-parts['second']/(12*np.pi)
        stiffness=7/24+parts['matter_stiffness']+parts['gradient']/(12*np.pi)
        assert mass>=1/24-1e-13 and stiffness>=1/8-1e-13
        energies,weights=spectral(T,D,h)
        comparisons=[]
        for w in [.4+.07j,1+.02j,-1+.3j,2+.2j]:
            z=w*w
            delta=np.sum(2*z*weights/(energies*(z-energies*energies)))
            kernel=stiffness-mass*z+delta
            expected=-z.imag*(stiffness/abs(z)**2+np.sum(2*weights/(energies*abs(energies*energies-z)**2)))
            error=abs((kernel/z).imag-expected)
            assert error<2e-12 and expected*z.imag<0
            comparisons.append(dict(w=[w.real,w.imag],kernel=[float(kernel.real),float(kernel.imag)],imaginary_identity_error=float(error)))
        for gamma in [.01,.2,1.,5.]:
            kernel=mass*gamma*gamma+stiffness+np.sum(2*gamma*gamma*weights/(energies*(energies*energies+gamma*gamma)))
            assert kernel>0
        rows.append(dict(D=D,h=h,mass_after_weyl_completion=mass,
            stiffness_after_weyl_completion=stiffness,upper_half_plane_checks=comparisons))
    return rows


def weyl_action_quadrature():
    T,D,h=np.pi,.7,.08
    f,fp,ds=profile(T,D,h,2048)
    n=np.fft.fftfreq(len(f))*len(f)
    fpp=(len(f)*np.fft.ifft(-(2*np.pi*n/(T+D))**2*coefficient(n,T,D,h))).real
    rows=[]
    for nt in [128,256]:
        t=2*np.pi*np.arange(nt)/nt
        phi=.2*np.cos(t);phidot=-.2*np.sin(t)
        a=np.expm1(phi)[:,None]
        v=1+a*f[None,:]
        vs=a*fp[None,:];vss=a*fpp[None,:]
        vt=(np.exp(phi)*phidot)[:,None]*f[None,:]
        direct=ds*(2*np.pi/nt)*np.sum(2*np.log(v)*vss+vs*vs/v-vt*vt/v**3)/(24*np.pi)
        integrated=-ds*(2*np.pi/nt)*np.sum(vs*vs/v+vt*vt/v**3)/(24*np.pi)
        assert abs(direct-integrated)<2e-12
        rows.append(dict(time_points=nt,spatial_points=len(f),direct_action=float(direct),
            integrated_action=float(integrated),difference=abs(direct-integrated)))
    return rows


def scale_check():
    a,b,s=sy.Rational(1,16),sy.Integer(2),sy.Rational(1,4)
    powers={'h_over_d':a-s,'dynamic_error_first':a-s,'dynamic_error_second':2*(a-s),
        'static_error_over_clock_stiffness':s,'radius_over_h':-(1+b)+s,
        'collar_gradient_over_nominal_lll_gap':s-sy.Rational(1,2)}
    assert all(value<0 for key,value in powers.items() if key!='static_error_over_clock_stiffness')
    assert powers['static_error_over_clock_stiffness']>0
    return dict(a=str(a),b=str(b),h_power=str(-s),powers={k:str(v) for k,v in powers.items()})


def main():
    started=time.perf_counter()
    report=dict(status='passed_collar_conservative_response_checks',
        note_sha256=sha(NOTE),script_sha256=sha(Path(__file__)),
        dependency_hashes={name:sha(ROOT/name) for name in [
            'notes/collective-clock-response.md','notes/long-return-response.md','notes/weak-coupling-window.md']},
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sy.__version__),
        exact_checks=exact_checks(),ground_energy_cases=ground_cases(),curvature_cases=curvature_cases(),
        mobius_cases=mobius_cases(),heat_kernel_cases=heat_cases(),dynamic_cases=dynamic_cases(),
        instability_cases=instability_cases(),weyl_completion_cases=weyl_completion_cases(),
        weyl_action_quadrature=weyl_action_quadrature(),scale_check=scale_check(),
        canonical_static_collar_obstruction=True,subtracted_dynamic_collar_bound=True,
        same_absorption_different_conservative_metric_completions=True,
        prescribed_weyl_completed_clock_has_no_upper_half_plane_poles=True,
        defined_lapse_clock_instability=True,physical_wormhole_instability_proved=False,
        physical_conservative_completion_derived=False,physical_damping_proved=False,
        novelty_established=False,lean_invocations=0)
    report['elapsed_seconds']=time.perf_counter()-started
    OUT.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=report['status'],seconds=report['elapsed_seconds'],
        maximum_fermion_ground_error=max(row['fermion_refinements'][-1]['error'] for row in report['ground_energy_cases']),
        maximum_curvature_error=max(row['richardson_error'] for row in report['curvature_cases']),
        growing_rates=[row['growing_rate'] for row in report['instability_cases']])) )


if __name__=='__main__':
    main()
