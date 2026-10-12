"""Circular companion boundary data. Flat retarded checks are not a GR junction.

Run: python checks/binary_mouth_drive.py
Requires numpy, scipy, sympy; no cloud, Lean, or external service is invoked.
"""
from pathlib import Path
import hashlib
import json
import platform

import numpy as np
import scipy
from scipy.optimize import brentq
import sympy as S

ROOT = Path(__file__).resolve().parents[1]


def symbolic_checks():
    checks = {}

    def zero(name, expr):
        result = S.simplify(expr)
        assert result == 0, (name, result)
        checks[name] = str(result)

    m1, m2, Q1, Q2, d = S.symbols('m1 m2 Q1 Q2 d', positive=True)
    static = -m1*m2*(m1+m2)/(2*d**2) + Q1*Q2*(m1+m2)/d**2 \
        - (m1*Q2**2+m2*Q1**2)/(2*d**2)
    printed = static - 2*Q1*Q2*(m1+m2)/d**2
    zero('same_sign_extremal_static_cancellation', static.subs({Q1:m1,Q2:m2}))
    zero('printed_source_static_discrepancy',
         printed.subs({Q1:m1,Q2:m2}) + 2*m1*m2*(m1+m2)/d**2)
    # Independent exact RN test-particle action in harmonic radius rho=R-M.
    z, M, Q, m, e, eps = S.symbols('z M Q m e eps', real=True)
    R = 1/z+M
    exact = -m*S.sqrt(1-2*M/R+Q**2/R**2)-e*Q/R
    series = S.series(exact,z,0,3).removeO()
    test = S.diff(static.subs({m1:M,Q1:Q,m2:eps*m,Q2:eps*e,d:1/z}),eps).subs(eps,0)
    zero('RN_test_action_1PN_static',series+m-(m*M-e*Q)*z-test)

    P, O, chi = S.symbols('P O chi', positive=True)
    L0 = P*d**2*O**2/4+2*P**2/d
    L1 = P*d**4*O**4/64+S.Rational(7,4)*P**2*d*O**2-4*P**3/d**2
    Omega2 = 4*P/d**3*(1-8*P/d)
    circular = S.diff(L0+L1,d).subs(O**2,Omega2).subs(P,chi*d)
    zero('circular_1PN_radial_equation',S.series(circular,chi,0,4).removeO())
    energy = (O*S.diff(L0+L1,O)-L0-L1).subs(O**2,Omega2)
    expected = -P**2/d+S.Rational(15,4)*P**3/d**2
    zero('circular_1PN_binding_energy',S.series((energy-expected).subs(P,chi*d),chi,0,4).removeO())
    # Full two-body Cartesian variation, before imposing a circular orbit.
    x = S.Matrix(S.symbols('x0:3',real=True))
    va = S.Matrix(S.symbols('va0:3',real=True))
    vb = S.Matrix(S.symbols('vb0:3',real=True))
    r = S.sqrt(x.dot(x)); n = x/r
    L = P*(va.dot(va)+vb.dot(vb))/2+2*P**2/r \
        + P*(va.dot(va)**2+vb.dot(vb)**2)/8 \
        + P**2/(2*r)*(3*(va.dot(va)+vb.dot(vb))-8*va.dot(vb)-2*n.dot(va)*n.dot(vb)) \
        - 4*P**3/r**2
    v = S.symbols('v',real=True)
    sub = dict(zip(list(x)+list(va)+list(vb),[d,0,0,0,v,0,0,-v,0]))
    a1 = S.Matrix([-2*P/d**2,0,0]); a2 = -a1
    mom = S.diff(L,va[0])
    inertial_1pn = sum(S.diff(mom,va[j])*a1[j]+S.diff(mom,vb[j])*a2[j]
                      +S.diff(mom,x[j])*(va[j]-vb[j]) for j in range(3))
    # Leading acceleration corrected by +16 P^2/d^3 at body 1.
    residual = (inertial_1pn-S.diff(L,x[0])).subs(sub)+16*P**3/d**3
    zero('independent_Cartesian_circular_force',residual.subs(v**2,P/d))

    def tr(a): return S.series(a,v,0,6).removeO().expand()
    def vec(a): return S.Matrix([tr(b) for b in a])
    theta=2*v-v**3+S.Rational(13,12)*v**5
    zero('lightcone_series',tr(theta-2*v*S.cos(theta/2)))
    nr=vec([S.cos(theta/2),-S.sin(theta/2),0])
    beta=vec([-v*S.sin(theta),-v*S.cos(theta),0])
    acc=vec([2*v*v*S.cos(theta),-2*v*v*S.sin(theta),0])
    k=tr(1-nr.dot(beta)); distance=tr(S.cos(theta/2))
    B=-vec(vec((1-v*v)*(nr-beta))*tr(1/(k**3*distance**2))
           +vec(nr.cross(vec((nr-beta).cross(acc))))*tr(1/(k**3*distance)))
    E=-vec(nr.cross(B)); u=S.Matrix([0,v,0]); gamma=tr(1/S.sqrt(1-v*v))
    Bp=vec(gamma*vec(B-u.cross(E))-tr(gamma**2/(gamma+1))*u*tr(u.dot(B)))
    Ep=vec(gamma*vec(E+u.cross(B))-tr(gamma**2/(gamma+1))*u*tr(u.dot(E)))
    desired_B=S.Matrix([-1-5*v**4,-S.Rational(8,3)*v**3+S.Rational(112,15)*v**5,0])
    desired_E=S.Matrix([0,0,2*v+2*v**5])
    for i in range(3):
        zero('flat_rest_B_series_'+str(i),Bp[i]-desired_B[i])
        zero('flat_rest_E_series_'+str(i),Ep[i]-desired_E[i])
    return dict(count=len(checks),identities=checks,flat_B_series=[str(t) for t in Bp],
                flat_E_series=[str(t) for t in Ep],
                printed_source_same_sign_static=str(printed.subs({Q1:m1,Q2:m2})))


def source(t, v):
    angle=2*v*t
    return (np.array([-.5*np.cos(angle),-.5*np.sin(angle),0.]),
            np.array([v*np.sin(angle),-v*np.cos(angle),0.]),
            np.array([2*v*v*np.cos(angle),2*v*v*np.sin(angle),0.]))


def retarded(t, x, v):
    def equation(tau): return tau+np.linalg.norm(x-source(tau,v)[0])-t
    tau=brentq(equation,t-4,t,xtol=5e-15)
    position,beta,acc=source(tau,v)
    delta=x-position; radius=np.linalg.norm(delta); normal=delta/radius
    return tau,radius,normal,beta,acc,1-normal@beta


def potential(t, x, v):
    _,radius,_,beta,_,k=retarded(t,x,v)
    phi=-1/(k*radius)
    return np.concatenate(([phi],phi*beta))


def analytic_fields(t, x, v):
    _,radius,n,beta,acc,k=retarded(t,x,v)
    B=-((1-beta@beta)*(n-beta)/(k**3*radius**2)
        +np.cross(n,np.cross(n-beta,acc))/(k**3*radius))
    E=-np.cross(n,B)
    return E,B


def potential_fields(t, x, v, h):
    point=np.concatenate(([t],x))
    derivatives=[]
    for j in range(4):
        step=np.eye(4)[j]*h
        def evaluate(p): return potential(p[0],p[1:],v)
        derivatives.append((-evaluate(point+2*step)+8*evaluate(point+step)
                            -8*evaluate(point-step)+evaluate(point-2*step))/(12*h))
    deriv=np.array(derivatives)
    # Electric-source E, B followed by duality: magnetic-source E=-B_e, B=E_e.
    B=-deriv[1:,0]-deriv[0,1:]
    E=-np.array([deriv[2,3]-deriv[3,2],deriv[3,1]-deriv[1,3],deriv[1,2]-deriv[2,1]])
    return E,B


def boost(E,B,u):
    gamma=1/np.sqrt(1-u@u); longitudinal=gamma*gamma/(gamma+1)
    return (gamma*(E+np.cross(u,B))-longitudinal*u*(u@E),
            gamma*(B-np.cross(u,E))-longitudinal*u*(u@B))


def retarded_checks():
    rows=[]
    for v in [.003,.01,.03,.07,.1,.2]:
        for phase in [0.,.37,1.2]:
            n=np.array([np.cos(phase),np.sin(phase),0.])
            tangent=np.array([-np.sin(phase),np.cos(phase),0.])
            t=phase/(2*v); x=.5*n; u=v*tangent
            E,B=analytic_fields(t,x,v); Ep,Bp=boost(E,B,u)
            errors=[]
            for h in [.008,.004,.002]:
                Ei,Bi=potential_fields(t,x,v,h)
                Epi,Bpi=boost(Ei,Bi,u)
                errors.append(float(max(np.linalg.norm(Epi-Ep),np.linalg.norm(Bpi-Bp))))
            assert errors[-1]<2e-9, errors
            assert errors[-1]<errors[0]/80, errors
            invariants=max(abs((B@B-E@E)-(Bp@Bp-Ep@Ep)),abs(E@B-Ep@Bp))
            assert invariants<2e-14
            theta=brentq(lambda a:a-2*v*np.cos(a/2),0,2*v,xtol=5e-15)
            tau,radius,*_=retarded(t,x,v)
            assert abs(2*v*(t-tau)-theta)<1e-13
            series_B=(-1-5*v**4)*n+(-8*v**3/3+112*v**5/15)*tangent
            series_E=np.array([0,0,2*v+2*v**5])
            error_series=max(np.linalg.norm(Bp-series_B),np.linalg.norm(Ep-series_E))
            assert error_series<30*v**6
            rows.append(dict(v=v,phase=phase,derivative_errors=errors,
                             invariant_error=float(invariants),series_error=float(error_series),
                             rest_B_radial=float(Bp@n),rest_B_tangential=float(Bp@tangent),
                             rest_E_normal=float(Ep[2]),theta=float(theta)))
    return rows


def port_checks():
    nodes,weights=np.polynomial.legendre.leggauss(80)
    rows=[]
    for p in [.001,.01,.1,.3]:
        # Potential +1/|n+y| has gradient -n at the origin.
        potential=1/np.sqrt(1+2*p*nodes+p*p)
        coefficient=3*np.sum(weights*potential*nodes)/(2*p)
        assert abs(coefficient+1)<1e-10
        differences=[]
        bound=2*p/(1-p)**3
        for mu in np.linspace(-1,1,101):
            y=p*np.array([mu,np.sqrt(max(0,1-mu*mu)),0])
            b=(np.array([1.,0,0])+y)/np.linalg.norm(np.array([1.,0,0])+y)**3
            differences.append(np.linalg.norm(b-np.array([1.,0,0])))
        assert max(differences)<=bound
        rows.append(dict(port_fraction=p,l1_coefficient=float(coefficient),
                         max_sampled_field_change=float(max(differences)),proved_field_bound=bound))
    return rows


def scale_checks():
    rows=[]
    for q in [1e3,1e6,1e9,1e12]:
        logq=np.log(q); lam=1/(16*np.pi*q**3); delta=logq**-2; p=q**-1.5
        chi=lam/delta; v=np.sqrt(chi); sigma=2*v/delta
        row=dict(q=q,chi=float(chi),v=float(v),magnetic_drive=float(chi**2),
                 electric_to_magnetic=float(2*v),port_variation_bound=float(2*p/(1-p)**3),
                 regular_clock_shift=float(1.5*chi),acceleration_port_scale=float(chi*p),
                 gravity_tidal_port_scale=float(chi*p*p),Omega_ell=float(sigma),
                 candidate_rigid_LLL_rotation_norm=float((q-1)*sigma/2),
                 conditional_antipodal_master_D=float(-1.5*lam*lam*chi*chi))
        assert np.isclose(sigma,logq**3/(2*np.sqrt(np.pi)*q**1.5),rtol=1e-14)
        rows.append(row)
    return rows


def digest(path): return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def main():
    result=dict(status='passed_binary_drive_components_not_full_junction',
                scope='Known 1PN orbit and leading binary mouth drive; independent flat-retarded verification',
                python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=S.__version__,
                symbolic=symbolic_checks(),retarded_cases=retarded_checks(),port_cases=port_checks(),
                scale_cases=scale_checks(),script_sha256=digest(Path(__file__)),
                input_sha256={n:digest(ROOT/n) for n in ['notes/binary-mouth-drive.md',
                  'notes/dynamic-rn-dipole.md','notes/magnetic-mouth-force.md','notes/lll-return-holonomy.md']},
                full_four_dimensional_binary_junction=False,quantum_state_or_response_derived=False,
                electric_parity_response_derived=False,novelty_established=False,lean_invocations=0,goal_complete=False)
    path=ROOT/'evidence/binary-mouth-drive.json'
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],symbolic_checks=result['symbolic']['count'],
                         retarded_cases=len(result['retarded_cases']),
                         worst_derivative_error=max(x['derivative_errors'][-1] for x in result['retarded_cases']),
                         port_cases=len(result['port_cases']),scale_cases=len(result['scale_cases']))))


if __name__=='__main__': main()
