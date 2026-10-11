"""Independent reduced-model checks; no claim of physical evolution or novelty."""
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize_scalar
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
rho = s.symbols('rho', real=True)
a, r, G, m, alpha, beta, c = s.symbols('a r G m alpha beta c', positive=True)
k = s.symbols('k', real=True)
f, fa = 1 + rho**2, 1 + a**2
exact = {}


def zero(name, value):
    v = s.simplify(s.expand(value))
    assert v == 0, (name, v)
    exact[name] = str(v)


# Integrals over a coordinate delta, evaluated on its support.
volume = 4*s.pi*r**3/s.sqrt(fa)
lapse = r*s.sqrt(fa)
density = m*s.sqrt(fa)/(4*s.pi*r**3)
Ttt = density*lapse**2
zero('proper_mass_integral', density*volume-m)
zero('worldline_action_energy', Ttt/lapse*volume-m*r*s.sqrt(fa))
printed = m*fa/(4*s.pi*r)
zero('printed_source_actual_mass', printed/lapse**2*volume-m/s.sqrt(fa))
zero('printed_source_killing_energy', printed/lapse*volume-m*r)
# Varying the static action with the lapse gives the same density.
N = s.symbols('N', positive=True)
zero('lapse_action_density', -s.diff(-m*N,N)/volume-density)
# The proper acceleration uses Gamma^rho_tt = f f'/2 and u^t=1/(r sqrt f).
acc_rho = s.diff(f,rho)/(2*r**2)
acc_norm = s.sqrt(r**2/f*acc_rho**2).subs(rho,a)
mu = m*a/(r*s.sqrt(fa))
zero('worldline_acceleration_force', m*acc_norm-mu)
zero('distributional_radial_conservation', -mu/(4*s.pi*r**2)+a/fa*density)
zero('dimensionless_string_tension', (8*s.pi*G*mu).subs(m,beta*r/G)-8*s.pi*beta*a/s.sqrt(fa))

base_phi = alpha*(1+rho*s.atan(rho))
kappa = 2*beta/s.sqrt(fa)
phis = [base_phi+kappa*k*rho,
        base_phi-kappa*((rho-a)-k*rho)]
base_gamma = -alpha*((rho**2+3)*rho*s.atan(rho)+rho**2-s.log(f))
gammas = [base_gamma+beta/s.sqrt(fa)*((a-rho)*f+2*a*(rho-a)**2-3*a*rho**2+(1-2*k)*rho**3),
          base_gamma+beta/s.sqrt(fa)*((rho-a)*f-3*a*rho**2+(1-2*k)*rho**3)]
for i,(phi,gamma) in enumerate(zip(phis,gammas)):
    string = i*2*beta*a/s.sqrt(fa)
    zero(f'energy_equation_side_{i}', f*s.diff(phi,rho,2)+rho*s.diff(phi,rho)-phi-alpha/f+string)
    zero(f'radial_constraint_side_{i}', rho*s.diff(phi,rho)-phi+alpha/f+string)
    zero(f'angular_equation_side_{i}', s.diff(gamma,rho,2)+f*s.diff(phi,rho,2)+2*rho*s.diff(phi,rho)+4*phi)
zero('phi_continuity', (phis[1]-phis[0]).subs(rho,a))
zero('body_delta_jump', (f*(s.diff(phis[1],rho)-s.diff(phis[0],rho))).subs(rho,a)+2*beta*s.sqrt(fa))
zero('gamma_continuity', (gammas[1]-gammas[0]).subs(rho,a))
zero('angular_delta_jump', (s.diff(gammas[1]-gammas[0],rho)+f*s.diff(phis[1]-phis[0],rho)).subs(rho,a))
AR, AL = s.pi*alpha/2-kappa*(1-k), s.pi*alpha/2-kappa*k
zero('right_overlap_slope', s.limit(phis[1]/rho,rho,s.oo)-AR)
zero('left_overlap_slope', s.limit(phis[0]/(-rho),rho,-s.oo)-AL)
zero('centered_mass_limit', (s.pi*alpha/2-beta).subs({alpha:G*c/(4*s.pi*r**2),beta:G*c/(8*r**2)}))

# Independent matching versus energy minimization with fixed magnetic radius.
A, B, x, lam = s.symbols('A B x lam', positive=True)
energy = A*x**2+(r*m-B)*x
xm = (B-r*m)/(2*A)
K = B**2/(4*A)
zero('stationary_energy', s.diff(energy,x).subs(x,xm))
zero('loaded_minimum_energy', energy.subs(x,xm).subs(m,lam*B/r)+K*(1-lam)**2)
zero('total_added_energy', (energy.subs(x,xm)+K).subs(m,lam*B/r)-K*(2*lam-lam**2))
zero('body_energy', (m*r*xm).subs(m,lam*B/r)-2*K*lam*(1-lam))
zero('reversible_first_law', s.diff(energy.subs(x,xm),m)-r*xm)
zero('overlap_matches_minimum', ((s.pi*alpha/2-beta)/(2*r)).subs({alpha:G*c/(4*s.pi*r**2),beta:G*m/r})-xm.subs({A:r**3/G,B:c/8}))
zero('body_energy_maximum', (2*K*lam*(1-lam)).subs(lam,s.Rational(1,2))-K/2)

# The MQ charge constraint in the symmetric sector yields its known potential.
v, vp, vpp, eta, Delta, Et, nf = s.symbols('v vp vpp eta Delta Et nf', real=True)
V = s.exp(2*v)-eta*s.exp(2*Delta*v)+Et*s.exp(v)/nf
constraint = vpp+s.exp(2*v)-eta*Delta*s.exp(2*Delta*v)+Et*s.exp(v)/(2*nf)
zero('MQ_charge_to_potential', vpp+s.diff(V,v)/2-constraint)
zero('MQ_half_dimension_quadratic', V.subs(Delta,s.Rational(1,2))-(s.exp(2*v)+(Et/nf-eta)*s.exp(v)))

# Independent numerical initial-value integration, crossing the delta jump.
ode_cases=[]
for height, kval in [(0.,.5),(.4,.2),(1.7,.35),(3.,0.)]:
    av,bv=.006,.001
    kap=2*bv/np.sqrt(1+height**2)
    initial=[av*(1+(-6)*np.arctan(-6))+kap*kval*(-6),
             av*(np.arctan(-6)-6/37)+kap*kval]
    def rhs(side):
        def fun(z,y):
            fz=1+z*z
            return [y[1],(av/fz-side*2*bv*height/np.sqrt(1+height**2)-z*y[1]+y[0])/fz]
        return fun
    sols=[]
    errors=[]
    for tol in [2e-9,2e-11]:
        left=solve_ivp(rhs(0),(-6,height),initial,rtol=tol,atol=tol/100,dense_output=True,method='DOP853')
        assert left.success
        jump=left.y[:,-1].copy()
        jump[1]-=2*bv/np.sqrt(1+height**2)
        right=solve_ivp(rhs(1),(height,6),jump,rtol=tol,atol=tol/100,dense_output=True,method='DOP853')
        assert right.success
        z=np.linspace(-6,6,161)
        expected=av*(1+z*np.arctan(z))-kap*(np.maximum(z-height,0)-kval*z)
        actual=np.array([(left if t<=height else right).sol(t)[0] for t in z])
        errors.append(float(np.max(np.abs(expected-actual))))
    assert errors[-1]<2e-10
    ode_cases.append(dict(height=height,k=kval,max_errors=errors))

# Smooth normalized rest-density approximations, independently integrated.
mass_cases=[]
for height in [0.,.5,2.,4.]:
    eps=.025
    norm=quad(lambda z: np.exp(-((z-height)/eps)**2),height-10*eps,height+10*eps,epsabs=1e-13)[0]
    # epsilon(z) = sqrt(f(z)) times the normalized coordinate mollifier / (4 pi).
    mass=quad(lambda z: np.sqrt(1+z*z)*np.exp(-((z-height)/eps)**2)/norm/np.sqrt(1+z*z),height-10*eps,height+10*eps,epsabs=1e-13)[0]
    E=quad(lambda z: np.sqrt(1+z*z)*np.exp(-((z-height)/eps)**2)/norm,height-10*eps,height+10*eps,epsabs=1e-13)[0]
    assert abs(mass-1)<1e-12
    assert abs(E-np.sqrt(1+height**2))<eps**2
    mass_cases.append(dict(height=height,proper_mass=mass,killing_energy=E,delta_limit=float(np.sqrt(1+height**2))))

energy_cases=[]
for load in [0.,.1,.4,.5,.8,.95]:
    # Set x0=1,K=1: E=x^2-2(1-load)x; no analytic minimizer passed to scipy.
    fn=lambda y:y*y-2*(1-load)*y
    result=minimize_scalar(fn,bounds=(0,2),method='bounded',options={'xatol':1e-14})
    assert result.success
    err=abs(result.x-(1-load))
    assert err<2e-8
    energy_cases.append(dict(lambda_value=load,minimizer=float(result.x),minimum=float(result.fun),minimizer_error=err,
        total_added_energy=float(result.fun+1),body_energy=float(2*load*result.x)))

note=ROOT/'notes/static-load-normalization.md'
result=dict(status='passed_static_load_normalization',
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    note_sha256=hashlib.sha256(note.read_bytes()).hexdigest(),
    versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=s.__version__),
    exact_checks=exact,ode_cases=ode_cases,proper_mass_cases=mass_cases,energy_cases=energy_cases,
    source_comparison=dict(publisher_pdf_sha256='53cac622417261aac1fc7101096c5c8a2af0dd2176acd4fe4126068bf1c92d4d',
        publisher_url='https://inspirehep.net/files/bcd01c1c0d090e073c0ad34a5f460a9b',
        equations=['EGMM (3.4)-(3.10)','EGMM (3.15)-(3.23), Appendix B','MQ v3 (4.37)-(4.43), Appendix B']),
    fixed_rest_mass_normalization_verified=True,centered_static_limit_changed=False,
    conservative_lowering_envelope_changed=False,static_load_potential_is_known=True,
    quantum_state_constructed=False,finite_pulse_result=False,full_physical_stability_proved=False,
    final_primary_fallback_selected=False,novelty_established=False,lean_invocations=0)
out=ROOT/'evidence/static-load-normalization.json'
out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(status=result['status'],exact_checks=len(exact),ode_cases=len(ode_cases),
    mass_cases=len(mass_cases),energy_cases=len(energy_cases),max_final_ode_error=max(t['max_errors'][-1] for t in ode_cases),
    max_minimizer_error=max(t['minimizer_error'] for t in energy_cases))))
