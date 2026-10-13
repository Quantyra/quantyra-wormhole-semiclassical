"""S072: full-spectrum return cancellation and its duration constraint.

Numerical checks support the ordinary Fourier proof; finite oscillator
evolution is not a continuum state, capture, or four-dimensional result.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import platform

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
DEP = ROOT / 'checks/passive_electric_feedback.py'
spec = importlib.util.spec_from_file_location('feedback', DEP)
feedback = importlib.util.module_from_spec(spec)
spec.loader.exec_module(feedback)


def smoothstep(x):
    x = np.clip(x, 0., 1.)
    return x**3 * (10. - 15.*x + 6.*x*x)


def integrated_step(t, tau):
    if t <= 0.:
        return 0.
    if t >= tau:
        return t - tau/2
    x = t/tau
    return tau*(2.5*x**4 - 3*x**5 + x**6)


def shaped(t, hf, tau):
    h = hf/np.pi*(integrated_step(t, tau)-integrated_step(t-np.pi, tau))
    hd = hf/np.pi*(smoothstep(t/tau)-smoothstep((t-np.pi)/tau))
    return h, hd


def transform(fun, omega, cuts):
    # Piecewise weighted quadrature, independent of Fourier factorization.
    c = sum(quad(fun, a, b, weight='cos', wvar=omega, epsabs=2e-13)[0]
            for a, b in zip(cuts[:-1], cuts[1:]))
    s = sum(quad(fun, a, b, weight='sin', wvar=omega, epsabs=2e-13)[0]
            for a, b in zip(cuts[:-1], cuts[1:]))
    return c-1j*s


def fourier_checks():
    rows = []
    for tau in [.15, .25, .7]:
        cuts = [0., tau, np.pi, np.pi+tau]
        f = lambda t: shaped(t, 1., tau)[1]
        area = transform(f, 0., cuts).real
        errors = [abs(transform(f, 2*n, cuts)) for n in range(1, 129)]
        assert abs(area-1) < 2e-12 and max(errors) < 2e-12
        detuned = []
        p = lambda t: 30*(t/tau)**2*(1-t/tau)**2/tau
        for n in [1, 3, 16, 64]:
            for delta in [-.15, .03, .2]:
                w = 2*n+delta
                direct = transform(f, w, cuts)
                factor = -np.expm1(-1j*np.pi*w)/(1j*np.pi*w)
                factored = factor*transform(p, w, [0., tau])
                bound = abs(delta)/(2*n-abs(delta))
                assert abs(direct-factored) < 2e-12 and abs(direct) <= bound+2e-12
                detuned.append(dict(n=n,delta=delta,amplitude=float(abs(direct)),bound=bound))
        rows.append(dict(tau=tau,duration=np.pi+tau,area=area,
                         maximum_even_transform= float(max(errors)),detuning=detuned))
    return rows


def fast_switch_checks():
    rows = []
    for T in [.1, .25, .8, 1.4]:
        for power in [1, 2, 4]:
            norm = quad(lambda x: (x*(1-x))**power, 0, 1)[0]
            f = lambda t: ((t/T)*(1-t/T))**power/(norm*T)
            value = transform(f, 2., [0., T])
            assert abs(value) >= np.cos(T)-2e-12
            rows.append(dict(T=T,power=power,first_amplitude=float(abs(value)),lower=float(np.cos(T))))
    k = 2000
    q = k**4
    hf = 4*np.pi/k
    T = .25
    cubic = 3*(np.sin(T)-T*np.cos(T))/T**3
    direct = abs(transform(lambda t: 6*t*(T-t)/T**3, 2., [0., T]))
    assert abs(cubic-direct) < 2e-13
    v, _ = feedback.spectra(31, 8)
    assert abs(v[2]-2*(31-1)/np.pi**2) < 1e-13
    eps = 300*float(feedback.mp.zeta(3))*hf/(np.pi**2*q*.25)
    return dict(monotone=rows,candidate=dict(k=k,T=T,
        first_stiffness_amplitude_over_mu_h=4*(1-1/q)*hf/np.pi*cubic,
        monotone_amplitude_lower=4*(1-1/q)*hf/np.pi*np.cos(T),
        shaped_stiffness_lower=1-2*hf/np.pi-eps,shaped_finite_q_loss=eps))


def oscillator_check(q, hf, tau, inertia, rtol=2e-10):
    number, nx = 8, 32
    Kstat, D, a, e, c, wr, sig, kstat = feedback.oscillator_model(q,hf,inertia,number)
    n = len(a)
    ex = np.arange(2,2*nx+1,2,dtype=float)
    vx, _ = feedback.spectra(q,2*nx)
    ax = 2*vx[ex.astype(int)]/ex
    mu = q/(2*np.pi)
    delta = float(feedback.static_delta(q))
    base = np.diag(np.r_[0.,wr*wr,e*e])
    aa = np.outer(a,a)
    counter = np.sum(c*c/(e*e))
    cross = np.zeros((n,n))
    for j,cj in enumerate(c):
        cross[:,j+3] -= cj*a
        cross[j+3,:] -= cj*a
    ns = (2*n)**2
    def rhs(t, state):
        S = state[:ns].reshape(2*n,2*n)
        z, s = state[ns:ns+nx], state[ns+nx:]
        h, hd = shaped(t,hf,tau)
        kappa = mu*h+h*h*delta-h*(ax@z)
        K = base+(kappa+h*h*counter)*aa+h*cross
        Sd = np.vstack((D[:,None]*S[n:],-K@S[:n])).ravel()
        return np.r_[Sd,hd-ex*s,ex*z]
    finish = np.pi+tau
    state0 = np.r_[np.eye(2*n).ravel(),np.zeros(2*nx)]
    sol = solve_ivp(rhs,(0,finish),state0,method='DOP853',rtol=rtol,atol=rtol/100,
                    max_step=.008,dense_output=True)
    assert sol.success
    state = sol.y[:,-1]
    S = state[:ns].reshape(2*n,2*n)
    residual = float(max(abs(state[ns:]))/hf)
    assert residual < 2e-8
    frequencies = np.r_[wr,e]
    variances = np.r_[.5,1/(2*frequencies),.5,frequencies/2]
    cov = (S*variances)@S.T
    B = np.block([[Kstat,np.zeros((n,n))],[np.zeros((n,n)),np.diag(D)]])
    J = np.block([[np.zeros((n,n)),np.eye(n)],[-np.eye(n),np.zeros((n,n))]])
    A = J@B
    energy = float(np.trace(B@cov)/2)
    errors, phase_ratios, uncertainty = [], [], []
    for t in [0.,.5,3.,17.,100.]:
        U = expm(A*t)
        evolved = U@cov@U.T
        errors.append(abs(float(np.trace(B@evolved)/2)/energy-1))
        va = float(a@evolved[:n,:n]@a)
        vp = float(evolved[n,n])
        ca = float(a@evolved[:n,n])
        uncertainty.append(va*vp-ca*ca)
        phase_ratios.append(kstat*va/(2*energy))
    assert max(errors) < 2e-10 and max(phase_ratios) <= 1+1e-10
    assert min(uncertainty) >= .25-1e-9
    # All-time stiffness proof with an explicit finite-band omission bound.
    Zbound = 75*hf/(8*np.pi*tau)
    model_error = 16*float(feedback.mp.zeta(3))/np.pi**2*Zbound
    tail = q*Zbound/(np.pi**2*(2*nx))
    lower = 1-2*hf/np.pi-(model_error+tail)/mu
    observed = []
    for t in np.linspace(.001,finish,250):
        h,_ = shaped(t,hf,tau)
        z = sol.sol(t)[ns:ns+nx]
        observed.append(1+h*delta/mu-(ax@z)/mu)
    assert min(observed) >= lower-1e-9
    return dict(q=q,hf=hf,tau=tau,I=inertia,neutral_modes=number,radial_modes=2,
                restoring_modes=nx,rtol=rtol,endpoint_return_over_hf=residual,
                final_phase_variance=float(a@cov[:n,:n]@a),
                maximum_energy_relative_error=float(max(errors)),
                maximum_phase_energy_ratio=float(max(phase_ratios)),
                minimum_uncertainty_determinant=float(min(uncertainty)),
                observed_stiffness_lower=float(min(observed)),proved_band_stiffness_lower=lower,
                symplectic_error=float(np.max(abs(S@J@S.T-J))))


def main():
    fourier = fourier_checks()
    fast = fast_switch_checks()
    dynamics = [oscillator_check(31,.03,.4,.8),oscillator_check(127,.08,.25,.7)]
    fine = oscillator_check(127,.08,.25,.7,2e-12)
    refinement = abs(fine['final_phase_variance']-dynamics[-1]['final_phase_variance'])
    assert refinement < 2e-7
    result = dict(date='2026-10-13',status='PASS',python=platform.python_version(),
                  numpy=np.__version__,scipy=scipy.__version__,
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  dependency_sha256=hashlib.sha256(DEP.read_bytes()).hexdigest(),
                  fourier=fourier,fast_switch=fast,dynamics=dynamics,refined=fine,
                  phase_variance_refinement=refinement,
                  limitations=['Exact even frequencies of fixed strip; second-order neutral influence and linear phase',
                    'Longer switch is outside the previous pre-return capture theorem',
                    'Conserved unsubtracted finite-band energy does not prove continuum quantum stress control',
                    'No creation, physical full stability, or novelty acceptance'])
    (ROOT/'evidence/electric-return-cancellation.json').write_text(
        json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: Fourier cancellation, short-switch lower bound, detuning, switched covariance and conserved final energy.')
    print('Candidate:',fast['candidate'])
    print('Endpoint residuals:',[x['endpoint_return_over_hf'] for x in dynamics])
    print('Phase variance refinement:',refinement)


if __name__ == '__main__':
    main()
