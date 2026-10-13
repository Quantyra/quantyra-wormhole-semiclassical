"""S072: collective electric response from the projected neutral kernels.

Tests ordinary spectral/distribution arguments and a resummed deterministic
linear comparison. This is not a quantum-noise or higher-Landau simulation.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
_spec=importlib.util.spec_from_file_location('prior_neutral',ROOT/'checks/neutral_transfer_current.py')
_prior=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_prior)


def spectral_comparison():
    mp.mp.dps=75
    limit=512
    harmonic=sum(mp.mpf(1)/e for e in range(1,limit+1))
    rows=[]
    for qi in [2,3,5,10,100,2000**4]:
        q=mp.mpf(qi);c=1/q
        fa=_prior.power_coefficients(1-c,limit)
        fc=_prior.power_coefficients(c,limit)
        bound=mp.mpf('4.5')*harmonic/q+1/q**2
        errors=[]
        for e in range(2,limit+1):
            if e%2==0:
                normalized=(fa[e]+(q-1)*fc[e])/4
                target=mp.mpf(e)
            else:
                normalized=(fa[e]-(q-1)*fc[e])/4
                target=e-mp.mpf(1)/e
            loss=1-normalized/target
            assert -mp.mpf('1e-60') <= loss <= bound+mp.mpf('1e-60')
            errors.append(loss)
        rows.append(dict(q=qi,maximum_energy=limit,
                         largest_relative_spectral_loss=mp.nstr(max(errors),30),
                         analytic_relative_bound=mp.nstr(bound,30)))
    return rows


def laplace_comparison():
    # Sum convergent remainders after the common Abel divergence has
    # been removed analytically. Quantities are per complex channel.
    rows=[]
    for s in [0.3+0j,2+0j,1+2j]:
        target_x=-s/(2*np.pi)/np.tanh(np.pi*s/2)
        target_y=-(s+1/s)/(2*np.pi)*np.tanh(np.pi*s/2)
        errors=[]
        for eps in [0.02,0.01,0.005,0.0025]:
            e=np.arange(1,int(np.ceil(36/eps))+1,dtype=float)
            decay=np.exp(-eps*e)
            even=e%2==0
            sx=np.sum(decay[even]/(s*s+e[even]**2))
            sy=np.sum(decay[~even]/(s*s+e[~even]**2))
            xx=2/np.pi**2*(1/np.expm1(2*eps)-1/(2*eps)-s*s*sx)
            yy=2/np.pi**2*(1/(2*np.sinh(eps))-1/(2*eps)-(s*s+1)*sy)
            error=max(abs(xx-target_x),abs(yy-target_y))
            errors.append(error)
        assert all(b<a for a,b in zip(errors,errors[1:]))
        assert errors[-1] < 0.04
        rows.append(dict(s_real=s.real,s_imaginary=s.imag,
                         abel_eps=[0.02,0.01,0.005,0.0025],
                         largest_kernel_error_per_channel=[float(x) for x in errors],
                         limit_x=[float(target_x.real),float(target_x.imag)],
                         limit_y=[float(target_y.real),float(target_y.imag)]))
    return rows


def sin_monomials(e,t,degree):
    integrals=[(1-np.cos(e*t))/e]
    if degree>=1: integrals.append(t/e-np.sin(e*t)/e**2)
    for p in range(2,degree+1):
        integrals.append(t**p/e-p*(p-1)/e**2*integrals[-2])
    return integrals


def first_return_comparison():
    # h(s)=.08*[3(s/T)^2-2(s/T)^3]; alpha(s)=.3+.2s+.1s^2.
    # Integrate every oscillatory convolution analytically before summing.
    T=.5
    hc=np.array([0.,0.,.08*3/T**2,-.08*2/T**3])
    ac=np.array([.3,.2,.1])
    hac=np.polynomial.polynomial.polymul(hc,ac)
    rows=[]
    for t in [.17,.31,.46]:
        h=np.polynomial.polynomial.polyval(t,hc)
        alpha=np.polynomial.polynomial.polyval(t,ac)
        alphadot=ac[1]+2*ac[2]*t
        mem=sum(v*t**(p+1)/(p+1) for p,v in enumerate(hac))
        target=-h*h*alphadot-h*mem
        errors=[];bubble_only=[];free_errors=[]
        for eps in [1/256,1/512,1/1024,1/2048,1/4096]:
            e=np.arange(1,int(np.ceil(36/eps))+1,dtype=float)
            even=e%2==0
            integrals=sin_monomials(e,t,len(hac)-1)
            convh=sum(v*integrals[p] for p,v in enumerate(hc))
            convha=sum(v*integrals[p] for p,v in enumerate(hac))
            decay=np.exp(-eps*e)
            # chi/mu=4/pi times the displayed spectral sum.
            xx=4/np.pi*np.sum(decay[even]*e[even]*convh[even])
            yy=4/np.pi*np.sum(decay[~even]*(e[~even]-1/e[~even])*convha[~even])
            unprojected=4/np.pi*np.sum(decay[~even]*e[~even]*convha[~even])
            combined=-h*alpha*xx+h*yy
            errors.append(abs(combined-target))
            bubble_only.append(abs(h*yy-target))
            free_errors.append(abs(-h*alpha*xx+h*unprojected+h*h*alphadot))
        assert all(b<a for a,b in zip(errors,errors[1:]))
        assert all(b>a for a,b in zip(bubble_only,bubble_only[1:]))
        assert errors[-1] < 1.1e-5
        rows.append(dict(t=t,target_second_order_current_over_mu=target,
                         abel_eps=[1/256,1/512,1/1024,1/2048,1/4096],
                         combined_errors=errors,bubble_only_errors=bubble_only,
                         unprojected_no_memory_errors=free_errors))
    return rows


def candidate(k):
    mp.mp.dps=70
    k=mp.mpf(k);q=k**4;eta=1/k;z0=mp.mpf('.5')
    X=(-1+mp.sqrt(1+12*mp.pi*k**6))/6
    C=1+6*X
    x=mp.pi*q/(8*eta*z0*C*(1+3*X))
    ir=mp.pi/(2*x)
    h=4*mp.pi*eta
    mu=q/(2*mp.pi)
    omega2=mu*h/ir
    gamma=mu*h*h/ir
    landau2=(q+1)/C
    ratio=omega2/(landau2**2)
    assert abs(ratio-q*q*C/(2*z0*(q+1)**2*(1+3*X)))<mp.mpf('1e-55')
    discr=gamma*gamma-4*omega2
    assert discr>0
    fast=(gamma+mp.sqrt(discr))/2
    slow=2*omega2/(gamma+mp.sqrt(discr))
    # Single-mouth relative-response time in global t units. Dividing
    # both inertia and anomaly by two for the relative coordinate leaves it unchanged.
    arad=2*mp.pi*ir/q
    band=k*k
    spectral=mp.mpf('4.5')*(1+mp.log(band))/q+1/q**2
    return dict(k=int(k),q=int(q),X=float(X),h=float(h),x=float(x),
                relative_inertia=float(ir),omega=float(mp.sqrt(omega2)),
                gamma=float(gamma),first_higher_landau_mass=float(mp.sqrt(landau2)),
                omega_squared_over_landau_mass_fourth=float(ratio),
                local_fast_rate=float(fast),local_slow_rate=float(slow),
                slow_rate_over_landau_mass=float(slow/mp.sqrt(landau2)),
                fast_rate_over_landau_mass=float(fast/mp.sqrt(landau2)),
                radial_static_comparison_parameter_at_omega=float(arad*mp.sqrt(omega2)),
                illustrative_band_energy=int(band),spectral_loss_bound=float(spectral))


def collective_switch(k, damping=True, rtol=2e-10):
    p=candidate(k)
    om2=p['omega']**2;gamma=p['gamma'];ramp=.25;end=.6
    def a(t):
        if t>=ramp:return 1.,0.
        u=t/ramp
        return 3*u*u-2*u**3,(6*u-6*u*u)/ramp
    def rhs(t,state):
        A,Ad=a(t);out=[]
        for z in state.reshape(2,6):
            alpha,v,m,work,rad,mem=z
            acc=-om2*A*alpha-gamma*A*m
            if damping:acc-=gamma*A*A*v
            out.extend([v,acc,A*alpha,om2*Ad*alpha*alpha/2,
                        gamma*A*A*v*v,gamma*A*m*v])
        return out
    initial=np.array([[1.,0.,0.,0.,0.,0.],[0.,1.,0.,0.,0.,0.]]).reshape(-1)
    def half_event(t,state):return state[0]-.5
    half_event.direction=-1
    sol=solve_ivp(rhs,(0,end),initial,method='DOP853',rtol=rtol,
                  atol=rtol*1e-3,max_step=1/p['omega'],events=half_event)
    assert sol.success
    final=sol.y[:,-1].reshape(2,6)
    var=np.array([p['x']/2,2*p['x']/np.pi**2])
    weighted=lambda values:float(np.dot(var,values))
    e=weighted((final[:,1]**2+om2*final[:,0]**2)/2)
    work=weighted(final[:,3]);rad=weighted(final[:,4]);mem=weighted(final[:,5])
    initial_energy=var[1]/2
    defect=e+(rad if damping else 0)+mem-initial_energy-work
    assert abs(defect)<2e-6*max(1,abs(work))
    # The propagated initial component alone is not a full quantum state.
    # For these initial minimum-uncertainty moments, det Cov(alpha,p)
    # equals det(response map)^2/4. A missing noise term is observable.
    homogeneous_uncertainty=np.linalg.det(final[:,:2].T)**2/4
    if damping: assert homogeneous_uncertainty<0.01
    return dict(k=k,damping_retained=damping,ramp_time=ramp,end_time=end,
                phase_energy_over_I=e,restoring_work_over_I=work,
                integrated_friction_load_over_I=rad,memory_work_over_I=mem,
                initial_kinetic_energy_over_I=float(initial_energy),
                energy_identity_error_over_I=float(abs(defect)),
                friction_load_over_leading_negative_energy=rad/om2,
                friction_load_over_supplied_energy=rad/(initial_energy+work),
                initial_displacement_half_time=float(sol.t_events[0][0]),
                inverse_half_time_over_landau_mass=float(1/(sol.t_events[0][0]*p['first_higher_landau_mass'])),
                homogeneous_alpha_p_covariance_determinant=float(homogeneous_uncertainty),
                final_phase_variance=weighted(final[:,0]**2),
                scope='Deterministic quadratic electric comparison of the initial compact Gaussian component; no neutral quantum-noise contribution or higher-mode load included.')


def main():
    damped=[collective_switch(k) for k in [2000,4000,8000]]
    refined=collective_switch(2000,rtol=1e-11)
    refine=abs(refined['integrated_friction_load_over_I']-damped[0]['integrated_friction_load_over_I'])
    assert refine<1e-8
    undamped=collective_switch(2000,damping=False)
    assert undamped['friction_load_over_supplied_energy']>10
    result=dict(date='2026-10-13',status='PASS',python=platform.python_version(),
                numpy=np.__version__,scipy=scipy.__version__,
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                inherited_checker_sha256=hashlib.sha256((ROOT/'checks/neutral_transfer_current.py').read_bytes()).hexdigest(),
                spectra=spectral_comparison(),laplace_kernels=laplace_comparison(),
                first_return=first_return_comparison(),
                candidate_comparisons=[candidate(k) for k in [2000,4000,8000]],
                damped_switches=damped,refinement_error_over_I=refine,
                undamped_posthoc_load=undamped,
                limitations=['Only the electrically projected, angular-held-fixed reference sector',
                             'The kernel is second order in normalized h and linear in alpha; resummation is not an all-order remainder bound',
                             'Global-time and static-inertia calibration is conditional on physical matching',
                             'First-return comparison is not late-time stability or the full noisy switching stress'])
    (ROOT/'evidence/collective-electric-response.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: spectral bounds, Abel kernels, first-return cancellation and collective response.')
    print('k=2000 scales:',result['candidate_comparisons'][0])
    print('Damped switch:',damped[0])
    print('Undamped posthoc load / supplied energy:',undamped['friction_load_over_supplied_energy'])
    print('ODE refinement:',refine)


if __name__=='__main__':main()
