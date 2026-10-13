"""S072: correlated fermion switching flux and symmetric JT capture.

Free Majoranas, scale-invariant junction, cutoff-boundary limit. The
first-return gate is checked. No four-dimensional error bound is claimed.
"""
from pathlib import Path
import hashlib
import json
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import sympy as s


def symbolic_checks():
    d, v = s.symbols('d v', real=True)
    trace_difference = ((s.cos(v*d)-1)/(2*s.pi*s.I*s.sin(d/2))
                        -s.I*s.sin(v*d)/(2*s.pi*s.cos(d/2)))
    flux = s.simplify(-s.I*s.limit(trace_difference/d, d, 0)/2)
    assert s.simplify(flux-(v*v-v)/(4*s.pi)) == 0
    z, w, a, r, ap, eta = s.symbols('z w a r ap eta', positive=True)
    e = w*w+z*z-(a-r)*z
    de = (s.diff(e,z)*z*w+s.diff(e,w)*(-z*z+(a-r)*z/2)
          +s.diff(e,a)*ap+s.diff(e,r)*4*s.pi*eta*ap*ap/z)
    assert s.simplify(de-(-ap*z+4*s.pi*eta*ap*ap)) == 0
    x = s.symbols('x', real=True)
    ramp_cost = s.integrate(s.diff(3*x*x-2*x**3,x)**2,(x,0,1))
    assert ramp_cost == s.Rational(6,5)
    Z = s.Rational(96,127)
    assert Z*Z+Z/8 < s.Rational(2,3)
    assert max(xx*xx-s.Rational(11,12)*xx for xx in (s.Rational(1,4),s.Rational(3,4))) == -s.Rational(1,8)
    # pi<22/7 proves both strict rational inequalities in the bootstrap.
    assert s.Rational(192,11)*s.Rational(22,7)*s.Rational(24,5)/2000 < s.Rational(1,4)
    assert 4*s.Rational(22,7)*s.Rational(24,5)/2000 < s.Rational(1,32)
    t, te = s.symbols('t te', real=True)
    xl, xr = t-te, s.pi-(t-te)
    Q1 = s.cos(t)*s.cos(xl)+s.sin(t)*s.sin(xl)+s.cos(t)*s.cos(xr)-s.sin(t)*s.sin(xr)
    Q2 = s.sin(t)*s.cos(xl)-s.cos(t)*s.sin(xl)+s.sin(t)*s.cos(xr)+s.cos(t)*s.sin(xr)
    assert s.trigsimp(Q1) == 0 and s.trigsimp(Q2) == 0
    # RN-AdS thermodynamics fixes the physical Schwarzian coefficient.
    aa, L, G, rr = s.symbols('a L G rr', positive=True)
    charge2 = aa*aa*(1+3*aa*aa/(L*L))
    M = (rr+rr**3/(L*L)+charge2/rr)/(2*G)
    T = (1+3*rr*rr/(L*L)-charge2/(rr*rr))/(4*s.pi*rr)
    C = 1+6*aa*aa/(L*L)
    phi_r = aa**3/(G*C)
    coeff = s.simplify(s.diff(M,rr,2).subs(rr,aa)/2 / s.diff(T,rr).subs(rr,aa)**2)
    assert s.simplify(coeff-2*s.pi*s.pi*phi_r) == 0
    return {'flux_per_real_field_two_ends': str(flux), 'cubic_ramp_cost': str(ramp_cost),
            'RNAdS_Schwarzian_coefficient': str(phi_r), 'paired_nonenergy_SL2_charges': [0,0]}


def point_split_checks():
    mp.mp.dps = 65
    J = mp.matrix([[0,-1],[1,0]])
    def rotation(th):
        return mp.matrix([[mp.cos(th),-mp.sin(th)],[mp.sin(th),mp.cos(th)]])
    rows = []
    for theta_f in (mp.mpf('.01'),mp.mpf('.2')):
        for T in (mp.mpf('.15'),mp.mpf('.35')):
            for fraction in (mp.mpf('.2'),mp.mpf('.5'),mp.mpf('.8')):
                time = fraction*T
                def theta(tt):
                    xx=tt/T
                    return theta_f*(3*xx*xx-2*xx**3)
                velocity=theta_f*6*fraction*(1-fraction)/T
                expected=(velocity*velocity-velocity)/(4*mp.pi)
                def trace_at(eps):
                    Gin = mp.eye(2)/(4*mp.pi*mp.j*mp.sin(eps/2))+mp.j*J/(4*mp.pi*mp.cos(eps/2))
                    diff=rotation(theta(time+eps/2))*Gin*rotation(theta(time-eps/2)).T-Gin
                    return diff[0,0]+diff[1,1]
                errors=[]
                for eps in (mp.mpf('1e-3'),mp.mpf('1e-4'),mp.mpf('1e-5')):
                    observed=-mp.j*(trace_at(eps)-trace_at(-eps))/(4*eps)
                    errors.append(abs(observed-expected))
                assert errors[-1] < mp.mpf('1e-8')
                assert errors[-1] < errors[0]/9000
                rows.append({'theta_f':float(theta_f),'global_duration':float(T),'fraction':float(fraction),
                             'predicted_flux':float(expected),'absolute_errors':list(map(float,errors))})
    return rows


def evolve(eta, duration, z0, tight=False):
    def rhs(sigma, state):
        z,w,r,tt,work,radiation_work=state
        x=sigma/duration
        a=x*x*(3-2*x)
        ap=6*x*(1-x)/duration
        return [z*w,-z*z+(a-r)*z/2,4*np.pi*eta*ap*ap/z,z,-ap*z,4*np.pi*eta*ap*ap]
    sol=solve_ivp(rhs,(0,duration),[z0,0,0,0,0,0],method='DOP853',
                  rtol=2e-13 if tight else 2e-10,atol=2e-15 if tight else 2e-12,
                  max_step=duration/(120 if tight else 35))
    assert sol.success
    z,w,r,tt,work,radiation_work=sol.y[:,-1]
    ef=w*w+z*z-(1-r)*z
    assert r < .25 and ef <= -3/32 and tt < np.pi
    assert np.min(sol.y[0]) >= 11*z0/12-1e-11
    balance=abs(ef-z0*z0-work-radiation_work)
    assert balance<1e-10
    # Subtract the exact final coupled-strip vacuum, not the linear term.
    excitation_over_N_eta=r-4*eta
    assert excitation_over_N_eta > 0
    A=1-r
    lower,upper=(A-np.sqrt(A*A+4*ef))/2,(A+np.sqrt(A*A+4*ef))/2
    def post(t,v):
        zz,ww=v
        return [zz*ww,-zz*zz+A*zz/2]
    after=solve_ivp(post,(0,120),[z,w],method='DOP853',rtol=2e-11,atol=2e-13,max_step=.12)
    assert after.success
    zz,ww=after.y
    drift=float(np.max(np.abs(ww*ww+zz*zz-A*zz-ef)))
    assert drift<1e-9 and min(zz)>=lower-1e-9 and max(zz)<=upper+1e-9
    return {'eta_f':eta,'duration_eta_times_u':duration,'initial_y_over_eta':z0,
            'radiation_charge_R_over_N_eta':float(r),'global_time_elapsed':float(tt),
            'final_energy_over_N_eta_squared':float(ef),'work_balance_error':float(balance),
            'post_energy_drift':drift,'excitation_above_final_vacuum_over_N_eta':float(excitation_over_N_eta),
            'turning_interval_y_over_eta':[float(lower),float(upper)],'tight':tight}


if __name__ == '__main__':
    symbolic=symbolic_checks()
    stress=point_split_checks()
    cases=[evolve(eta,duration,z0) for eta in (1e-5,1e-4,5e-4)
           for duration in (.25,.375,.5) for z0 in (.25,.5,.75)]
    repeat=evolve(5e-4,.25,.75,True)
    ref=next(row for row in cases if row['eta_f']==5e-4 and row['duration_eta_times_u']==.25 and row['initial_y_over_eta']==.75)
    convergence=abs(repeat['final_energy_over_N_eta_squared']-ref['final_energy_over_N_eta_squared'])
    assert convergence<1e-10
    out={'date':'2026-10-12','status':'free_JT_switching_capture_verified',
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'sympy':s.__version__,'mpmath':mp.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
         'symbolic':symbolic,'point_split_checks':stress,'dynamics_checks':cases,
         'refinement':repeat,'refinement_energy_difference':convergence,
         'theorem_energy_margin':'3/32','uniform_redshift_lower_bound':float(.5-np.sqrt(5/32)),
         'scope':['free massless Majoranas and symmetric scale-invariant junction',
                  'initial AdS2 strip vacuum / prepared pair; no initial state creation',
                  'large-N semiclassical JT boundary constraint and zero relative charges',
                  'first-return condition proved during switching; constant coupling thereafter'],
         'not_proved':['four-dimensional matching error bounds','gauge and exterior response',
                       'full bulk dilaton/area positivity throughout the formation',
                       'initial-pair preparation','other perturbation sectors','physical realization']}
    (Path(__file__).resolve().parents[1]/'evidence/switching-radiation-capture.json').write_text(
        json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: correlated point-split stress, matter charges, capture bootstrap, 27 coupled runs and refinement')
