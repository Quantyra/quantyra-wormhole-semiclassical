"""S138 frozen-contract checks; supplementary to the reviewed ordinary proof.

Run with the designated Python from the satellite root. All created text is
UTF-8/LF. Raw input bytes are SHA-256 hashed before any calculation and checked
again at completion. Every attempt is retained in the JSON run list and log;
failed attempts append exact base64 input snapshots, without replacing history.
"""
from __future__ import annotations

import base64
import contextlib
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import traceback

import mpmath as mp
import numpy as np
import scipy
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "evidence/finite-pulse-higher-mode-work.json"
LOG = ROOT / "evidence/finite-pulse-higher-mode-work.log"
INPUT_NAMES = [
    "checks/finite_pulse_higher_mode_work.py",
    "notes/finite-pulse-higher-mode-work.md",
    "notes/rn-shifted-boundary-moments.md",
    "notes/rn-higher-mode-boundary-spectral-control.md",
    "notes/coupled-fermion-domain.md",
    "notes/laboratory-drive-clock-feedback.md",
    "notes/common-drive-validated-capture.md",
]
CONTRACT = {
    "scalar_cases": [[18, ".05", 1], [32, ".2", 3]],
    "scalar_orders": [2, 3, 4], "decimal_precisions": [60, 90],
    "scalar_tolerances": ["2e-8", "2e-9"],
    "scalar_tail": "8/[3*pi*(2*r-1)]*P**(1-2*r), P=max(100,4*mu)",
    "ode_refinements": [["2e-10", "2e-12"], ["2e-12", "2e-14"]],
    "finite_tolerances": ["2e-7", "2e-8"],
    "target": "abs(W_higher,U)/(nu*N*eta)<=1e-8; proved bound<1e-34",
    "scope": "ordinary continuum theorem with supplementary scalar/finite tests; no RN discretization",
}
CHECKS: list[dict] = []


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (F, mp.mpf, sp.Basic)):
        return str(value)
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def check(condition, label, **details):
    record = {"label": label, "passed": bool(condition), **clean(details)}
    CHECKS.append(record)
    if not record["passed"]:
        print("FAIL", label, json.dumps(clean(details), sort_keys=True), flush=True)
    return record["passed"]


def exact_algebra():
    x, y, s = sp.symbols("x y s", real=True)
    den = (x*x+s*s)**2*(y*y+s*s)**2
    g = lambda z: z/(z*z+s*s)**2
    a = lambda z: 1/(z*z+s*s)**2
    ng = s**4-x*y*(x*x+x*y+y*y+2*s*s)
    na = -(x+y)*(x*x+y*y+2*s*s)
    for label, fun, numerator in [("f", g, ng), ("a", a, na)]:
        residual = sp.factor((fun(x)-fun(y))*den-(x-y)*numerator)
        check(residual == 0, "rational_DD_"+label, residual=residual)
        diagonal = sp.simplify(numerator.subs(y, x)/den.subs(y, x)-sp.diff(fun(x), x))
        check(diagonal == 0, "rational_DD_diagonal_"+label, residual=diagonal)
    u, v = sp.symbols("u v", real=True)
    check(sp.simplify(sp.cosh(u)**2+v**2*sp.sinh(u)**2
                      -(1+(1+v**2)*sp.sinh(u)**2)) == 0,
          "below_barrier_density_denominator_identity")
    constants = {}
    for r, expected in [(2, sp.Rational(8, 3)), (3, sp.Rational(2, 3)), (4, sp.Rational(1, 3))]:
        cr = sp.simplify(8*sp.gamma(r-sp.Rational(3,2))/(3*sp.sqrt(sp.pi)*sp.gamma(r)))
        check(cr == expected, "moment_gamma_constant_r"+str(r), exact=cr)
        constants[str(r)] = str(cr)
    return {"moment_c_r": constants, "algebraic_residuals": "exact SymPy identities"}


def exact_inequalities():
    b = F(24, 25)
    assertions = {
        "f_local_250": 128/b**4+64/b**2+24 < 250,
        "f_tail_5": F(1,3)/b**5+F(4,3)/b**3+F(8,3)/b < 5,
        "f_mixed_3_squared": F(16,3)/b**3 < 9,
        "a_local_350": 64+256/b**2 < 350,
        "a_mixed_2_squared": F(8,3)/b**5 < 4,
        "a_mixed_5_squared": F(64,3)/b**3 < 25,
        "a_mixed_20_squared": F(1024,3)/b**3 < 400,
        "a_mixed_16_squared": F(512,3)/b**7 < 256,
        "a_mixed_43_squared": F(4096,3)/b**7 < 1849,
        "a_pure_tail_first_3": F(8,3)/b**2 < 3,
        "a_pure_tail_second_3_squared": F(32,9)/b**8 < 9,
        "a_mixed_sum_100": 2+5+20+16+43 < 100,
        "safe_sinh_slack": F(7,8)**2 > F(3,4),
        "exp4_greater8": sum(F(4)**j/math.factorial(j) for j in range(4)) > 8,
        "sqrt6_less_5over2_squared": 6 < F(5,2)**2,
        "sqrt18_less5_squared": 18 < 25,
        "exp18over5_greater32": sum(F(18,5)**j/math.factorial(j) for j in range(8)) > 32,
        "beta_local3_bound10": 2*(4+F(1,1)/F(5,4)) < 10,
        "beta_local4_bound10": 2*(4+F(1,1)/F(7,4)) < 10,
        "beta_tail1_bound8": 4+F(1,1)/F(1,4) <= 8,
        "beta_tail1half_bound6": 4+F(1,1)/F(1,2) <= 6,
        "beta_tail2_bound16over3": 4+F(1,1)/F(3,4) <= F(16,3),
        "a_integrated_constants": 3500 < 10000 and 6*F(16,3) < 100 and 100*6 < 1000,
        "f_integrated_constants": 2500 < 10000 and 5*8 < 100 and 3*6 < 1000,
        "pi_upper_allows_X": 3*F(41,40)**2 > F(22,7),
        "C_upper_implies_M1": F(123,20)+F(1,2000**3) < F(500,81),
        "delta_upper": F(3,1000*2000**9) < F(1,10**32),
        "angular_prefactor": 14*(1+F(1,2000**4)) < 15,
        "angular_pure": F(15,2000)*F(5,4) < F(1,100),
        "angular_tail2": F(15,2000)*2 < F(1,50),
        "angular_tail1half": F(15,2000)*F(10,32) < F(1,400),
        "dyadic_5power": 5**7 < 18**4,
        "dyadic_geometric": F(1,2) < F(3,4)**2,
        "dyadic_700": 8*18*4 < 700,
        "low_energy_4": 2**7 < 4**4,
    }
    delta = F(1,10**32)
    tower = 12500*F(1,10**80)+50*delta**2+30*F(1,10**48)
    assertions["full_angular_common_C_1e-45"] = tower < F(1,10**45)
    coefficient = (2000*(F(8,1000)+4000)+2*F(8,1000)+F(101,100))*1000
    assertions["normalized_work_1e-34"] = coefficient*F(1,10**45) < F(1,10**34)
    assertions["normalized_target_1e-8"] = F(1,10**34) < F(1,10**8)
    for label, passed in assertions.items():
        check(passed, "exact_"+label)
    return {"comparisons": len(assertions), "common_tower_upper": str(tower),
            "normalized_coefficient_upper": str(coefficient),
            "proof": "Fraction comparisons; exponential lower bounds are finite positive Taylor sums"}


def scalar_green(z, mu):
    a = mp.sqrt(z)
    k = mp.sqrt(z+mu*mu)
    return (mp.cosh(k)+(a/k)*mp.sinh(k))/(k*mp.sinh(k)+a*mp.cosh(k))


def scalar_spectral_moment(mu, a, r, P):
    def low(p):
        if not p:
            return mp.mpf(0)
        kap = mp.sqrt(max(mp.mpf(0), mu*mu-p*p))
        density = (2/mp.pi)/(mp.cosh(kap)**2+(kap/p)**2*mp.sinh(kap)**2)
        return density/(p*p+a*a)**r
    cuts = sorted(set([mp.mpf(0), min(a/2,mu), min(a,mu), min(2*a,mu), mu/2, mu]))
    below = mp.quadgl(low, cuts, maxdegree=9)
    qmax = mp.sqrt(P*P-mu*mu)
    def high(q):
        if not q:
            return mp.mpf(0)
        p = mp.sqrt(mu*mu+q*q)
        density = (2/mp.pi)/(mp.cos(q)**2+(q/p)**2*mp.sin(q)**2)
        return density*(q/p)/(p*p+a*a)**r
    cutsq = [mp.mpf(0)]
    while cutsq[-1]+mp.pi/2 < qmax:
        cutsq.append(cutsq[-1]+mp.pi/2)
    cutsq.append(qmax)
    above = mp.quadgl(high, cutsq, maxdegree=9)
    return below+above, below, above


def scalar_checks():
    results = []
    for precision, tolerance in [(60,"2e-8"),(90,"2e-9")]:
        with mp.workdps(precision):
            for M, delta_text, shift in CONTRACT["scalar_cases"]:
                nu = mp.mpf(delta_text)
                mu = mp.mpf(M)/(3*mp.sqrt(2))
                kap = mp.mpf(7)/25*nu
                a = mp.sqrt((nu*shift)**2-kap*kap)
                P = max(mp.mpf(100),4*mu)
                check(P >= 4*mu, f"scalar_tail_density_condition_{M}_{precision}")
                check(scalar_green(a*a,mu) <=2/mp.sqrt(mu*mu+a*a),
                      f"scalar_resolvent_bound_{M}_{precision}")
                for r in CONTRACT["scalar_orders"]:
                    exact = (-1)**(r-1)*mp.diff(lambda z: scalar_green(z,mu), a*a, r-1)/mp.factorial(r-1)
                    spectral, below, above = scalar_spectral_moment(mu,a,r,P)
                    tail = 8/(3*mp.pi*(2*r-1))*P**(1-2*r)
                    discrepancy = abs(exact-spectral)
                    allowed = mp.mpf(tolerance)*max(1,abs(exact))+tail
                    check(discrepancy <= allowed,
                          f"scalar_green_vs_dp_M{M}_r{r}_dps{precision}",
                          raw_difference=mp.nstr(exact-spectral,25),
                          analytic_tail=mp.nstr(tail,25), allowance=mp.nstr(allowed,25))
                    cr = 8*mp.gamma(r-mp.mpf("1.5"))/(3*mp.sqrt(mp.pi)*mp.gamma(r))
                    upper = 2*4**(r-1)*(mu*mu+a*a)**(mp.mpf(".5")-r) + cr*mp.exp(-mp.sqrt(3)*mu)*mu**-2*a**(3-2*r)
                    check(exact <= upper, f"scalar_moment_upper_M{M}_r{r}_dps{precision}")
                    normalized = 2*nu**(2*r-1)*exact
                    A = mu/nu
                    ks = mp.sqrt(A*A+(mp.mpf(24)/25*shift)**2)
                    normalized_upper = 4**r*ks**(1-2*r)+2*cr*mp.exp(-mp.mpf(M)/mp.sqrt(6))*A**-2*(mp.mpf(24)/25*shift)**(3-2*r)
                    check(normalized <= normalized_upper,
                          f"normalized_folded_moment_M{M}_r{r}_dps{precision}")
                    results.append({"M":M,"delta":delta_text,"shift":shift,"r":r,
                                    "dps":precision,"P":mp.nstr(P,25),
                                    "green_derivative":mp.nstr(exact,45),
                                    "spectral_dp":mp.nstr(spectral,45),
                                    "below":mp.nstr(below,30),"above":mp.nstr(above,30),
                                    "raw_difference":mp.nstr(exact-spectral,30),
                                    "tail_bound":mp.nstr(tail,30),"upper":mp.nstr(upper,30)})
                    print(f"scalar M={M} r={r} dps={precision}: diff={mp.nstr(exact-spectral,8)} tail={mp.nstr(tail,8)}",flush=True)
    for coarse in [v for v in results if v["dps"]==60]:
        refined = next(v for v in results if v["dps"]==90 and
                       (v["M"],v["r"])==(coarse["M"],coarse["r"]))
        with mp.workdps(90):
            drift = abs(mp.mpf(coarse["spectral_dp"])-mp.mpf(refined["spectral_dp"]))
            check(drift <=mp.mpf("2e-9")*max(1,abs(mp.mpf(refined["spectral_dp"]))),
                  f"scalar_precision_refinement_M{coarse['M']}_r{coarse['r']}", drift=mp.nstr(drift,20))
    return results


# SOURCE_DIAGNOSTIC_BEGIN: response_audit owns this section during authoring.
def trace_norm(matrix):
    return float(np.linalg.svd(matrix, compute_uv=False).sum())


def source_symbolic_checks():
    s = sp.symbols("s", real=True)
    p, w = sp.Function("p")(s), sp.Function("w")(s)
    P = sp.symbols("P0:5")
    W = sp.symbols("W0:5")
    replacements = {sp.diff(p,s,j):P[j] for j in range(5)}
    replacements.update({sp.diff(w,s,j):W[j] for j in range(5)})
    q = p*w
    expressions = []
    for j in range(1,5):
        q = sp.expand(w*sp.diff(q,s))
        expressions.append(sp.expand(q/w).xreplace(replacements))
    expected_lower=[P[1]*W[0]+P[0]*W[1],
        P[2]*W[0]**2+3*P[1]*W[0]*W[1]+P[0]*(W[0]*W[2]+W[1]**2),
        P[3]*W[0]**3+6*P[2]*W[0]**2*W[1]
        +P[1]*(7*W[0]*W[1]**2+4*W[0]**2*W[2])
        +P[0]*(W[1]**3+4*W[0]*W[1]*W[2]+W[0]**2*W[3])]
    for order,expected in enumerate(expected_lower,1):
        residual_lower=sp.expand(expressions[order-1]-expected)
        check(residual_lower==0,f"source_chain_identity_j{order}",residual=residual_lower)
    expected4 = (P[4]*W[0]**4+10*P[3]*W[0]**3*W[1]
        +P[2]*(25*W[0]**2*W[1]**2+10*W[0]**3*W[2])
        +P[1]*(15*W[0]*W[1]**3+30*W[0]**2*W[1]*W[2]+5*W[0]**3*W[3])
        +P[0]*(W[1]**4+11*W[0]*W[1]**2*W[2]+4*W[0]**2*W[2]**2
                 +7*W[0]**2*W[1]*W[3]+W[0]**3*W[4]))
    residual = sp.expand(expressions[3]-expected4)
    check(residual == 0,"source_fourth_chain_identity", residual=residual)
    Fj = [F(7,1000),F(1,250),F(3,50),F(4),F(123)]
    Wj = [F(2),F(2),F(3),F(3),F(5)]
    chain_bounds = []
    for j,expr in enumerate(expressions,1):
        bound = F(0)
        for powers,coeff in sp.Poly(expr,*P,*W).terms():
            p_powers, w_powers = powers[:5],powers[5:]
            check(sum(p_powers)==1,f"source_chain_linear_profile_j{j}_{powers}")
            term = abs(F(int(sp.numer(coeff)),int(sp.denom(coeff))))
            for idx,power in enumerate(p_powers):
                term *= Fj[idx]**power
            for idx,power in enumerate(w_powers):
                term *= Wj[idx]**power
            bound +=term
        chain_bounds.append(bound)
        check(bound < [F(1),F(1),F(100),F(3000)][j-1],
              f"source_variation_bound_a{j}", exact_upper=str(bound))
    S = 10*s**3-15*s**4+6*s**5
    check(sp.factor(sp.diff(S,s)-30*s**2*(1-s)**2)==0,
          "quintic_monotone_derivative")
    check(sp.integrate(sp.diff(S,s),(s,0,1))==1,"quintic_total_increment")
    primitive_S=sp.Rational(5,2)*s**4-3*s**5+s**6
    primitive_S2=sp.Rational(100,7)*s**7-sp.Rational(300,8)*s**8+sp.Rational(345,9)*s**9-18*s**10+sp.Rational(36,11)*s**11
    check(sp.expand(sp.diff(primitive_S,s)-S)==0,"source_theta_primitive_exact")
    check(sp.expand(sp.diff(primitive_S2,s)-S**2)==0,"source_Q_primitive_exact")
    source_assertions = {
        "theta_f_upper":F(22,7*500)<F(7,1000),
        "F1_exact":2*F(1,500)==F(1,250),
        "F2_exact":2*F(1,500)*4*F(15,4)==F(3,50),
        "F3_upper":2*F(1,500)*16*60<4,
        "F4_measure_upper":2*F(1,500)*64*(360+120)<123,
        "gprime_upper":(F(1,500)+F(1,500)**2)/(F(6,1000))<1,
        "gsecond_upper":(F(15,1000)+4*F(1,500)*F(15,1000))/F(6,1000)<3,
        "w_wprime_upper":F(53,200)+F(17,5)/2<2,
        "wsecond_upper":2+F(1,2)<3,
        "wthird_upper":2+F(1,2)<3,
        "wfourth_upper":3+F(3,2)<5,
    }
    for label,condition in source_assertions.items():
        check(condition,"source_exact_"+label)
    t = sp.symbols("t")
    aj = sp.symbols("a0:5")
    At = sum(aj[j]*t**j/sp.factorial(j) for j in range(5))
    b4 = sp.diff(At*sp.exp(At/2),t,4).subs(t,0)
    expected_b4 = sp.exp(aj[0]/2)*((1+aj[0]/2)*aj[4]
        +(1+aj[0]/4)*(4*aj[1]*aj[3]+3*aj[2]**2)
        +(sp.Rational(3,4)+aj[0]/8)*6*aj[1]**2*aj[2]
        +(sp.Rational(1,2)+aj[0]/16)*aj[1]**4)
    check(sp.simplify(b4-expected_b4)==0,"tagged_Dyson_B4_factorials_exact")
    a0,a1,a2,a3,a4 = F(7,1000),F(1),F(1),F(100),F(3000)
    exp_upper=1/(1-a0/2)
    B0=exp_upper*a0
    B1=exp_upper*(1+a0/2)*a1
    B4=exp_upper*((1+a0/2)*a4+(1+a0/4)*(4*a1*a3+3*a2*a2)
        +(F(3,4)+a0/8)*6*a1*a1*a2+(F(1,2)+a0/16)*a1**4)
    check(B0<F(8,1000) and B1<F(101,100) and B4<4000,
          "source_B0_B1_B4_uniform", B0=str(B0),B1=str(B1),B4=str(B4))
    theta=sp.Function("theta")(s)
    coeff=sp.Function("F")(theta)
    rho=sp.diff(theta,s)
    c=sp.symbols("c")
    lag=c*coeff*rho**2
    euler=sp.diff(lag,theta)-sp.diff(sp.diff(lag,rho),s)
    contact_residual=sp.simplify(-rho*euler-sp.diff(c*coeff*rho**2,s))
    check(contact_residual==0,"intrinsic_contact_Euler_work_exact", residual=contact_residual)
    return {"phase_variation_bounds":[str(v) for v in Fj],
            "laboratory_chain_upper":[str(v) for v in chain_bounds],
            "B_upper":[str(B0),str(B1),str(B4)],
            "source_fourth_derivative":"finite measure; endpoint atoms included in F4"}


def source_functions():
    eta=1/2000
    tau=.25
    def smooth(v):
        if v<=0: return 0.
        if v>=1: return 1.
        return v**3*(10+v*(-15+6*v))
    def smooth_prime(v):
        return 30*v*v*(1-v)**2 if 0<v<1 else 0.
    def primitive(v, squared=False):
        if v<=0: return 0.
        if squared:
            inside=lambda x:(100*x**7/7-300*x**8/8+345*x**9/9
                              -180*x**10/10+36*x**11/11)
        else:
            inside=lambda x:2.5*x**4-3*x**5+x**6
        return inside(v) if v<1 else inside(1)+(v-1)
    def p(s):
        return 4*eta*(smooth(s/tau)-smooth((s-math.pi)/tau))
    def pp(s):
        return (4*eta/tau)*(smooth_prime(s/tau)-smooth_prime((s-math.pi)/tau))
    def forcing(s):
        theta=4*eta*tau*(primitive(s/tau)-primitive((s-math.pi)/tau))
        Q=16*eta*eta*tau*(primitive(s/tau,True)-primitive((s-math.pi)/tau,True))
        return (theta-Q)/(4*math.pi*eta)
    return p,pp,forcing


def finite_evolution(rtol,atol,level):
    from scipy.integrate import solve_ivp,quad
    from scipy.linalg import expm
    p,pp,forcing=source_functions()
    cuts=[0.,.25,math.pi,math.pi+.25]
    reference_segments=[]
    ref_initial=np.array([.265,0.,0.])
    for left,right in zip(cuts[:-1],cuts[1:]):
        sol=solve_ivp(lambda s,y:[y[1],-y[0]+forcing(s)/2,1/y[0]],
                      (left,right),ref_initial,method="DOP853",rtol=rtol,
                      atol=atol,dense_output=True)
        check(sol.success,f"reference_ODE_success_level{level}_{left}",message=sol.message)
        reference_segments.append(sol)
        ref_initial=sol.y[:,-1]
    def ref(s):
        idx=min(2,max(0,int(np.searchsorted(cuts,s,side="right")-1)))
        return reference_segments[idx].sol(s)
    reference_samples=[]
    for phase in sorted(set(list(np.linspace(0,cuts[-1],13))+cuts)):
        integration_cuts=[v for v in cuts if 0<v<phase]
        integ=quad(lambda v:math.sin(phase-v)*forcing(v),0,phase,
                   points=integration_cuts,epsabs=2e-13,epsrel=2e-13)[0] if phase else 0.
        integp=quad(lambda v:math.cos(phase-v)*forcing(v),0,phase,
                    points=integration_cuts,epsabs=2e-13,epsrel=2e-13)[0] if phase else 0.
        exact=np.array([.265*math.cos(phase)+integ/2,
                        -.265*math.sin(phase)+integp/2])
        raw=ref(phase)[:2]-exact
        reference_samples.append({"phase":float(phase),"ODE":ref(phase)[:2].tolist(),
                                  "forced_Green":exact.tolist(),"raw_difference":raw.tolist()})
        check(float(np.max(np.abs(raw)))<[2e-7,2e-8][level],
              f"reference_vs_forced_Green_level{level}_{phase}",raw_difference=raw.tolist())
    angle=quad(p,0,cuts[-1],points=cuts[1:-1],epsabs=2e-14,epsrel=2e-14)[0]
    check(abs(angle-math.pi/500)<2e-13,f"actual_source_total_angle_level{level}",raw_difference=angle-math.pi/500)
    check(p(0)==0 and p(cuts[-1])==0,
          f"actual_source_contact_endpoint_zero_level{level}",endpoint_profile=[p(0),p(cuts[-1])])
    h=np.diag([-3.7,-1.3,-.15,.15,1.3,3.7]).astype(complex)
    j=np.arange(6)
    fourier=np.exp(2j*math.pi*np.outer(j,j)/6)/math.sqrt(6)
    K=fourier@np.diag([.5,.5,.5,-.5,-.5,-.5])@fourier.conj().T
    check(np.linalg.norm(h@K-K@h)>1,f"synthetic_noncommutation_level{level}")
    check(np.max(np.abs(np.linalg.eigvalsh(K)-np.array([-.5]*3+[.5]*3)))<2e-14,
          f"synthetic_K_eigenvalues_level{level}")
    initial=np.concatenate([np.eye(6,dtype=complex).ravel(),np.zeros(36,dtype=complex)])
    nfev=0
    def rhs(s,y):
        w,wp,_=ref(s)
        state=y[:36].reshape(6,6)
        rho_s=pp(s)*w+p(s)*wp
        ds=-1j*(h/w+p(s)*K)@state
        work=rho_s*(state.conj().T@K@state)
        return np.concatenate([ds.ravel(),work.ravel()])
    for left,right in zip(cuts[:-1],cuts[1:]):
        sol=solve_ivp(rhs,(left,right),initial,method="DOP853",rtol=rtol,atol=atol)
        check(sol.success,f"finite_ODE_success_level{level}_{left}",message=sol.message)
        initial=sol.y[:,-1]
        nfev+=sol.nfev
    state=initial[:36].reshape(6,6)
    work_integral=initial[36:].reshape(6,6)
    lab_end=float(ref(cuts[-1])[2])
    U=expm(1j*h*lab_end)@state
    delta=U.conj().T@h@U-h
    comm=h@U-U@h
    energies=np.diag(h).real
    avec=(1+energies**2)**(-1.75)
    a=np.diag(avec)
    f=np.diag(energies*avec)
    ca=a@U-U@a
    cf=f@U-U@f
    ad_ca=h@ca-ca@h
    identity=comm@a-(cf-ad_ca-ca@h)
    cb_f=trace_norm(f@K-K@f)
    cb_a=trace_norm(a@K-K@a)
    B0,B1,B4=.008,1.01,4000.
    ad4cf=cf*(energies[:,None]-energies[None,:])**4
    diagnostics={
        "unitarity":float(np.linalg.norm(U.conj().T@U-np.eye(6),ord=2)),
        "work_integral_matrix":float(np.linalg.norm(delta-work_integral,ord=2)),
        "thermal_algebra_matrix":float(np.linalg.norm(identity,ord=2)),
        "work_hermiticity":float(np.linalg.norm(delta-delta.conj().T,ord=2)),
    }
    tolerance=[2e-7,2e-8][level]
    for label,error in diagnostics.items():
        check(error<tolerance,f"finite_{label}_level{level}",raw_discrepancy=error,tolerance=tolerance)
    endpoint_bounds={"f_comm":{"observed":trace_norm(cf),"bound":B0*cb_f},
                     "a_comm":{"observed":trace_norm(ca),"bound":B0*cb_a},
                     "a_first_comm":{"observed":trace_norm(ad_ca),"bound":B1*cb_a},
                     "f_fourth_comm":{"observed":trace_norm(ad4cf),"bound":B4*cb_f}}
    for label,values in endpoint_bounds.items():
        check(values["observed"]<=values["bound"]+tolerance,
              f"finite_tagged_bound_{label}_level{level}",**values)
    minus=np.where(energies<0)[0]
    plus=np.where(energies>0)[0]
    B=U[np.ix_(plus,minus)]
    C=U[np.ix_(minus,plus)]
    vacuum=float(np.sum(energies[plus,None]*np.abs(B)**2)
                 +np.sum(-energies[minus,None]*np.abs(C)**2))
    Pminus=np.diag((energies<0).astype(float))
    vacuum_direct=float(np.trace(h@(U@Pminus@U.conj().T-Pminus)).real)
    covariance_cases=[]
    for zi in [.25,.28]:
        sigma=math.pi/zi
        for phase in [0.,.37]:
            X=np.zeros((6,6),complex)
            for neg,pos in [(0,5),(1,4),(2,3)]:
                occupation=1/(1+math.exp(2*sigma*energies[pos]))
                cross=math.sqrt(occupation*(1-occupation))*np.exp(1j*phase)
                X[pos,pos]=occupation
                X[neg,neg]=-occupation
                X[pos,neg]=cross
                X[neg,pos]=cross.conjugate()
            covariance=Pminus+X
            thermal_complex=np.trace(U.conj().T@comm@X)
            thermal=float(thermal_complex.real)
            total=vacuum+thermal
            direct_complex=np.trace(h@(U@covariance@U.conj().T-covariance))
            integrated_complex=np.trace(work_integral@covariance)
            weighted=np.linalg.norm(np.diag(1/avec)@X,ord=2)
            weighted_h=np.linalg.norm(np.diag(energies/avec)@X,ord=2)
            errs={"vacuum_crossing_identity":abs(vacuum-vacuum_direct),
                  "total_relative_trace":abs(total-direct_complex),
                  "integrated_work":abs(total-integrated_complex),
                  "thermal_trace_reality":abs(thermal_complex.imag),
                  "covariance_projector":np.linalg.norm(covariance@covariance-covariance,ord=2)}
            for label,error in errs.items():
                check(float(error)<tolerance,f"finite_{label}_level{level}_zi{zi}_phase{phase}",
                      raw_discrepancy=float(error),tolerance=tolerance)
            check(weighted<=1+1e-13 and weighted_h<=1+1e-13,
                  f"finite_TFD_weighted_norms_level{level}_zi{zi}_phase{phase}",
                  a_inverse_X=float(weighted),h_a_inverse_X=float(weighted_h))
            covariance_cases.append({"z_i":zi,"phase":phase,"vacuum":vacuum,
                "thermal_TFD":thermal,"total_relative":total,
                "direct_total":clean(complex(direct_complex)),
                "integrated_total":clean(complex(integrated_complex)),
                "weighted_norms":[float(weighted),float(weighted_h)],
                "raw_discrepancies":{key:float(value) for key,value in errs.items()}})
    return {"level":level,"rtol":rtol,"atol":atol,"nfev":nfev,
            "laboratory_completion":lab_end,"angle":angle,"reference_samples":reference_samples,
            "raw_discrepancies":diagnostics,"finite_static_commutator_norms":{"f":cb_f,"a":cb_a},
            "endpoint_bounds":endpoint_bounds,"TFD_cases":covariance_cases,
            "U_raw":U,"work_matrix_raw":work_integral}


def source_and_finite_checks():
    symbolic=source_symbolic_checks()
    results=[finite_evolution(2e-10,2e-12,0),finite_evolution(2e-12,2e-14,1)]
    refinement={"U_operator":float(np.linalg.norm(results[0]["U_raw"]-results[1]["U_raw"],ord=2)),
                "work_matrix_operator":float(np.linalg.norm(results[0]["work_matrix_raw"]-results[1]["work_matrix_raw"],ord=2)),
                "completion":abs(results[0]["laboratory_completion"]-results[1]["laboratory_completion"])}
    for label,error in refinement.items():
        check(error<2e-8,"finite_refinement_"+label,raw_discrepancy=error)
    for coarse,refined in zip(results[0]["TFD_cases"],results[1]["TFD_cases"]):
        error=abs(coarse["total_relative"]-refined["total_relative"])
        check(error<2e-8,f"finite_refinement_total_zi{coarse['z_i']}_phase{coarse['phase']}",raw_discrepancy=error)
    for result in results:
        for key in ["U_raw","work_matrix_raw"]:
            matrix=result[key]
            result[key]={"real":matrix.real.tolist(),"imag":matrix.imag.tolist()}
    return {"symbolic_source":symbolic,"finite_refinements":results,"refinement_raw":refinement,
            "scope":"finite noncommuting identity diagnostic; no RN C_n imposed on synthetic matrix"}
# SOURCE_DIAGNOSTIC_END


class Tee:
    def __init__(self, stream, logfile):
        self.stream, self.logfile = stream, logfile
    def write(self, value):
        self.stream.write(value)
        self.logfile.write(value.replace("\r\n","\n").replace("\r","\n"))
        self.logfile.flush()
        return len(value)
    def flush(self):
        self.stream.flush()
        self.logfile.flush()


def json_write(payload):
    RECEIPT.write_text(json.dumps(clean(payload),indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")


def main():
    snapshots = {name:(ROOT/name).read_bytes() for name in INPUT_NAMES}
    hashes = {name:hashlib.sha256(data).hexdigest() for name,data in snapshots.items()}
    existing = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs":[]}
    run = {"run_id":datetime.now(timezone.utc).isoformat(),"status":"running",
           "contract":CONTRACT,"execution_inputs_sha256":hashes,
           "command":list(sys.orig_argv),"cwd":str(Path.cwd()),
           "text_policy":"raw-byte SHA-256; current-turn source/notes and all outputs UTF-8/LF without BOM; historical inputs preserved raw",
           "versions":{"python":platform.python_version(),"numpy":np.__version__,
                       "scipy":scipy.__version__,"sympy":sp.__version__,"mpmath":mp.__version__}}
    existing["runs"].append(run)
    json_write(existing)
    exit_status = 1
    with LOG.open("a",encoding="utf-8",newline="\n") as logfile:
        with contextlib.redirect_stdout(Tee(sys.stdout,logfile)), contextlib.redirect_stderr(Tee(sys.stderr,logfile)):
            print("\nRUN",run["run_id"],flush=True)
            print("COMMAND",json.dumps(run["command"]),"CWD",run["cwd"],flush=True)
            print("INPUT_SHA256",json.dumps(hashes,sort_keys=True),flush=True)
            print("PREDECLARED_CONTRACT",json.dumps(CONTRACT,sort_keys=True),flush=True)
            try:
                for name in INPUT_NAMES[:3]:
                    data = snapshots[name]
                    check(not data.startswith(b"\xef\xbb\xbf") and b"\r" not in data,
                          "UTF8_LF_noBOM_"+name)
                    data.decode("utf-8")
                run["exact_algebra"] = exact_algebra()
                run["exact_inequalities"] = exact_inequalities()
                run["source_and_finite"] = source_and_finite_checks()
                run["scalar_comparisons"] = scalar_checks()
            except BaseException as exc:
                run["exception"] = {"type":type(exc).__name__,"message":str(exc)}
                traceback.print_exc()
                check(False,"execution_exception",error=str(exc))
            final_hashes = {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in INPUT_NAMES}
            check(final_hashes==hashes,"input_hashes_unchanged", final_sha256=final_hashes)
            run["final_inputs_sha256"] = final_hashes
            run["checks"] = CHECKS
            run["failed_checks"] = [v for v in CHECKS if not v["passed"]]
            exit_status = 1 if run["failed_checks"] else 0
            run["terminal_exit_status"] = exit_status
            run["status"] = "passed" if exit_status==0 else "failed"
            run["completed_utc"] = datetime.now(timezone.utc).isoformat()
            if exit_status:
                logfile.write("FAILED_RUN_EXACT_INPUT_SNAPSHOTS_BASE64_BEGIN\n")
                for name,data in snapshots.items():
                    logfile.write(json.dumps({"path":name,"sha256":hashes[name],
                                      "base64":base64.b64encode(data).decode("ascii")})+"\n")
                logfile.write("FAILED_RUN_EXACT_INPUT_SNAPSHOTS_BASE64_END\n")
                logfile.flush()
                print("Failed attempt's exact input snapshots retained in the raw log.",flush=True)
            print("TERMINAL_EXIT_STATUS",exit_status,"CHECKS",len(CHECKS),
                  "FAILURES",len(run["failed_checks"]),flush=True)
            json_write(existing)
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
