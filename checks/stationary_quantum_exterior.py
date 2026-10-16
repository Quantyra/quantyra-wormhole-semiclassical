"""S139 focused checks, supplementary to the ordinary global component proof.

Do not execute until the frozen contract and this source have been reviewed.
Raw input bytes are hashed before/after execution. Each attempt is retained;
failed attempts retain exact input snapshots. Outputs are strict UTF-8/LF.
"""
from __future__ import annotations

import base64
import contextlib
from datetime import datetime, timezone
from fractions import Fraction as Q
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
from scipy.integrate import solve_ivp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "evidence/stationary-quantum-exterior.json"
LOG = ROOT / "evidence/stationary-quantum-exterior.log"
INPUTS = [
    "checks/stationary_quantum_exterior.py",
    "notes/stationary-lll-quantum-exterior.md",
    "notes/stationary-exterior-conformal-domain.md",
    "notes/exact-spherical-em-reduction.md",
    "notes/spherical-anomaly-action.md",
    "notes/spherical-anomaly-sources.md",
    "notes/dirac-spherical-source-reduction.md",
]
CONTRACT = {
    "synthetic_parameters": {"epsilon": "1e-10", "xi": "3", "h": ["1e-3", "1e-4"]},
    "integration": "t=log((r-r_h)/h), r_h=1+h, r_max=20",
    "refinements": [
        {"rtol": "2e-9", "atol": "2e-11", "offset_over_h": "1e-6", "tolerance": "2e-4"},
        {"rtol": "2e-11", "atol": "2e-13", "offset_over_h": "1e-8", "tolerance": "2e-5"},
    ],
    "horizon_seed": "p=p_h+p'_h delta, w=w'_h delta+w''_h delta^2/2; integral p through delta^2",
    "sample_gaps": ["0.01h", "0.1h", "h", "10h", "0.01", "0.1", "1", "4", "19-h"],
    "constraint_scale": "abs(C)/(epsilon*max(1,F^2,kappa^2)); C is exact HH first integral",
    "refinement_scale": "abs(coarse-refined)/max(1,abs(refined)) separately for p,w,integral_p",
    "refinement_tolerance": "2e-5",
    "precision": [80, 110],
    "benchmark_scaled_root_tolerance": "1e-60",
    "scope": "actual isolated Polyakov-LLL ODE at synthetic parameters; supplementary to ordinary uniform proof",
}
CHECKS = []


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, (Q, mp.mpf, sp.Basic)):
        return str(value)
    return value


def check(condition, label, **details):
    row = {"label": label, "passed": bool(condition), **clean(details)}
    CHECKS.append(row)
    if not row["passed"]:
        print("FAIL", label, json.dumps(clean(details), sort_keys=True), flush=True)


def zero(expr, label):
    residual = sp.factor(sp.together(expr))
    check(residual == 0, label, residual=residual)


def exact_bulk():
    x = sp.symbols("x")
    r, psi = sp.Function("r")(x), sp.Function("psi")(x)
    g, kap = sp.symbols("gamma kappa")
    U = sp.Function("U")
    rx, vx = sp.diff(r, x), sp.diff(psi, x)
    lag = 4*r*rx*vx+6*rx**2+sp.exp(2*psi)*U(r)+2*g*vx**2
    er = sp.diff(lag, r)-sp.diff(sp.diff(lag, rx), x)
    ep = sp.diff(lag, psi)-sp.diff(sp.diff(lag, vx), x)
    zero(er+4*(r*sp.diff(psi, x, 2)+3*sp.diff(r, x, 2)-sp.exp(2*psi)*sp.diff(U(r), r)/4),
         "action_radius_EL")
    zero(ep+4*(g*sp.diff(psi, x, 2)+r*sp.diff(r, x, 2)+rx**2-sp.exp(2*psi)*U(r)/2),
         "action_metric_EL")
    energy = 4*r*rx*vx+6*rx**2+2*g*(vx**2-kap**2)-sp.exp(2*psi)*U(r)
    zero(sp.diff(energy, x)+er*rx+ep*vx, "HH_first_integral_Noether")
    omega = psi+sp.log(r)
    tt = kap**2+2*sp.diff(psi, x, 2)-vx**2
    xx = kap**2-vx**2
    zero(sp.diff(xx, x)+(sp.diff(omega, x)-rx/r)*(tt-xx), "stress_radius_Ward")
    trace = sp.symbols("trace")
    jx = -2*g*trace/r**2
    zero(trace+r**2*jx/(2*g), "stress_radius_trace_dictionary")
    R, F, gam, lam, charge = sp.symbols("R F gamma lambda Q", positive=True)
    D = R**2-3*gam
    H = R**2+3*lam*R**4-charge-gam*(1+6*lam*R**2)
    a = H/(R*D)
    b = R/D
    A = a-b*F
    B = (-2-3*lam*R**2+3*charge/R**2+3*F)/D
    S = 3*gam*(-F*(R**2+6*gam)-3*charge+3*gam+R**2
               +6*gam*lam*R**2+lam*R**4)/(R**2*D**2)
    formed = 2*B+2*A/R-2*F/R**2-sp.diff(a, R)+sp.diff(b, R)*F+b*A
    zero(formed-S, "exact_radial_S")
    ell = sp.symbols("ell")
    fp = A-F*ell
    lp = S/F-(b+2*A/F)*ell
    zero((fp*ell+F*lp)-(S-(b+A/F)*F*ell-F*ell**2), "quadratic_lapse_cancellation")
    zero(ell*F+fp+b*F-H/(R*D), "first_integrating_factor")
    zero(lp+2*ell**2+2*ell*fp/F+b*ell-S/F, "second_integrating_factor")
    xi, eps, t = sp.symbols("xi epsilon t", positive=True)
    P = xi*(t*t+2*t+3)+1
    Text = 2/(t-1)+sp.diff(P, t)/P-2/t
    zero((sp.Rational(3, 2)*Text-1/t-1/(t-1))*t*(t-1)*P
         -(xi*(t**3-t+12)+4-2*t), "barrier_positive_polynomial")
    zero((t**3-3*t+16)-14-(t-1)**2*(t+2), "barrier_polynomial_minimum")
    Qstar = 1+3*xi
    Hstar = t*t+3*xi*t**4-Qstar
    Hq = Hstar-6*eps*xi*(t*t-1)
    aq = Hq/(t*(t*t-3*eps))
    astar = Hstar/t**3
    theta = 6*eps/(1/xi+3*(t*t+1))
    zero(aq/astar-(1-theta)/(1-3*eps/t**2), "exact_reference_forcing_ratio")
    Fext = (t-1)**2*P/t**2
    zero(sp.diff(Fext, t)-2*xi*(t-1)
         -(2*xi*(t**3+2*t-3)+2*(t-1))/t**3, "reference_derivative_lower_bound")
    h, r0 = sp.symbols("h r0", positive=True)
    U0 = (2*R**2+6*lam*R**4-2*charge)
    zero(R*B+3*A-sp.diff(U0,R)/(4*R**2), "optical_radius_EL_elimination")
    zero(gam*B+R*A+F-U0/(2*R**2), "optical_metric_EL_elimination")
    zero((2*R*U0-gam*sp.diff(U0, R))/4-R*H, "exact_horizon_clock_numerator")
    zero(H.subs({R: r0, charge: r0**2+3*lam*r0**4-gam*(1+6*lam*r0**2)}),
         "exact_quantum_center")
    z = sp.symbols("z", positive=True)
    primitive = sp.log((z-1)/z)+1/z
    zero(sp.diff(primitive, z)-1/(z*z*(z-1)), "barrier_integral_primitive")
    u = sp.symbols("u", positive=True)
    primitive_F = lam*u**3+(1+3*gam*lam)*u+(charge+gam)/sp.sqrt(3*gam)*sp.atan(sp.sqrt(3*gam)/u)
    zero(sp.diff(primitive_F,u)-(3*lam*u*u+1+3*gam*lam-(charge+gam)/(u*u+3*gam)),
         "zero_lapse_Fbar_primitive")
    zero(lam*(R*R-3*gam)+1+3*gam*lam-(lam*R*R+1), "Fbar_constant_one")
    m1,m2 = sp.symbols("mass1 mass2")
    asymptotic_F = lam*R*R+1+m1/R+m2/R**2
    numerator = sp.cancel((S*(R*R*D*D)/(3*gam)).subs(F,asymptotic_F))
    check(sp.Poly(sp.expand(numerator*R**2),R).degree()<=3, "asymptotic_source_numerator_order_R")
    check(sp.limit(R*(b+2*A/F).subs(F,asymptotic_F),R,sp.oo)==5, "asymptotic_lapse_damping_five")
    d,a0,a1,b0,s0,s1,p0,p1,w1,w2,c0,c1 = sp.symbols("delta a0 a1 b0 s0 s1 p0 p1 w1 w2 c0 c1", nonzero=True)
    seedp0=s0/(2*a0)
    seedp1=(s1-2*a1*seedp0+a0*b0*seedp0)/(3*a0)
    residual_p = p1-((s0+s1*d)-2*(a0+a1*d)*(p0+p1*d))/(a0*d)-b0*p0
    zero(sp.limit(residual_p.subs({p0:seedp0,p1:seedp1}),d,0), "analytic_horizon_p_seed")
    zero((w2-(c1-b0*w1-seedp0*a0)).subs({w1:c0,w2:c1-b0*c0-seedp0*a0}),
         "analytic_horizon_w_seed")
    return {"normalization": "gamma=Gq/(12pi); stress equations divided by q/(24pi) in Ward check"}


def exact_domain_and_counterterm():
    y, GN, inverse_L2 = sp.symbols("y G inverse_L2")
    ads_trace = -sp.Rational(11,360)*24*y*y/(16*sp.pi**2)
    zero(ads_trace+11*y*y/(240*sp.pi**2), "conditional_outer_AdS_trace")
    zero(3*y-3*inverse_L2-8*sp.pi*GN*ads_trace/4
         -3*(y+11*GN*y*y/(360*sp.pi)-inverse_L2), "conditional_outer_AdS_Einstein")
    r,u,b0,m,Qmag,lam = sp.symbols("r u b0 mass Qmag lambda", positive=True)
    physical_F = 1+lam*r*r-2*m/r+Qmag/r**2
    quartic = lam*b0*b0+u*u-2*m*u**3/b0+Qmag*u**4/b0**2
    zero(u*u*physical_F.subs(r,b0/u)-quartic, "conformal_product_quartic")
    curvature = r*r/b0**2*(-sp.diff(physical_F,r,2)+2*sp.diff(physical_F,r)/r-2*physical_F/r**2)
    zero(curvature.subs(r,b0/u)+sp.diff(quartic,u,2)/b0**2, "conformal_orbit_curvature")
    zero(sp.diff(quartic,u,3).subs(u,0)+12*m/b0, "outer_third_geometry_jet")
    zero(-sp.diff(quartic,u,2).subs(u,0)/b0**2+2/b0**2, "outer_product_scalar_zero")
    s1 = sp.Matrix([[0, 1], [1, 0]])
    s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    s3 = sp.diag(1, -1)
    eye = sp.eye(2)
    gn = sp.kronecker_product(eye, s1, eye)
    tangent = [sp.kronecker_product(eye, s2, eye),
               sp.kronecker_product(eye, s3, s1),
               sp.kronecker_product(eye, s3, s2)]
    tau = sp.kronecker_product(s3, sp.eye(4))
    chi = tau*gn
    pp, pm = (sp.eye(8)+chi)/2, (sp.eye(8)-chi)/2
    check(chi*chi == sp.eye(8) and pp.rank() == 4, "physical_projector_rank_involution")
    check(chi*gn == gn*chi and pm*gn*pp == sp.zeros(8), "complementary_adjoint_Green_form")
    check(pp*gn*pp != sp.zeros(8), "first_order_domain_not_selfadjoint")
    ks = sp.symbols("K0:3")
    robin = sp.zeros(8)
    for j, ga in enumerate(tangent):
        check(chi*ga+ga*chi == sp.zeros(8), "tangential_projector_flip_"+str(j))
        robin += pp*ga*(tau*ga*ks[j]/2)*pp
    check(sp.simplify(robin-sum(ks)*gn*pp/2) == sp.zeros(8), "induced_Robin_K_over_2")
    sigma, sn, field, dn, khat = sp.symbols("sigma sigma_n field field_n Khat")
    transformed = sp.exp(-sigma)*(sp.exp(-sp.Rational(3,2)*sigma)*(dn-sp.Rational(3,2)*sn*field)
                  +sp.exp(-sp.Rational(3,2)*sigma)*(khat+3*sn)*field/2)
    zero(transformed-sp.exp(-sp.Rational(5,2)*sigma)*(dn+khat*field/2),
         "finite_cut_conformal_Robin_covariance")
    rf = sp.symbols("fiber_rank")
    zero(((240-120)*rf/2*(-sp.Rational(1,4))+(42-18)*rf/2)/360+rf/120,
         "boundary_a4_normal_flux")
    zero((60*(-rf/4)+12*rf)/360+rf/120, "bulk_a4_derivative_coefficient")
    r = sp.symbols("r", positive=True)
    f = sp.Function("F")(r)
    fp = sp.diff(f, r)
    A = -sp.diff(f, r, 2)/2-fp/r
    B = (1-f-r*fp)/r**2
    density = r**2*((2*(A+B))**2-2*(2*A*A+2*B*B))
    euler = sp.diff(density, f)-sp.diff(sp.diff(density, fp), r)+sp.diff(sp.diff(density, sp.diff(f,r,2)),r,2)
    zero(euler+4*(sp.diff(f,r,2)+2*(1-f)/r**2), "blind_counterterm_radial_variation")
    mass, charge, lam, c, G = sp.symbols("mass Q lambda c G")
    rn = 1+lam*r*r-2*mass/r+charge/r**2
    zero(euler.subs(f, rn).doit()+16*charge/r**4, "blind_counterterm_RN_source")
    alpha = -32*sp.pi*G*c*charge/r**4
    zero(sp.diff(alpha,r)-128*sp.pi*G*c*charge/r**5, "blind_counterterm_lapse")
    check(4*c+2*(-2*c) == 0, "blind_product_pressure_combination")
    return {"domain_scope": "DdaggerD modulus; no phase/eta cancellation asserted"}


def exact_inequalities():
    e = Q(1, 10**10)
    comparisons = {
        "forcing_upper": (1-3*e)*(1+4*e) > 1,
        "bootstrap_A_lower": Q(1,1)/Q(101,100)-(1-Q(1,1)/Q(101,100)+4*e) > Q(98,100),
        "bootstrap_A_upper": (1+4*e)/Q(99,100)+((1+4*e)/Q(99,100)-1) < Q(1022,1000),
        "barrier_damping": (Q(1,2)-Q(4,100))/Q(1022,1000) > Q(4,10),
        "source_numerator": Q(1602,100)+Q(2112,100)*e < Q(161,10),
        "source_50": 3*Q(161,10)/(1-3*e)**2 < 50,
        "source_ratio_5000": 50/(Q(97,100)*2) < 5000,
        "barrier_majorization": Q(1,10)*10**6 > 5000,
        "integrating_factor_tail": 3/(2*(1-3*e)) < 2,
        "geometry_first_jet": Q(2,1000)+6*Q(201,100)+4*Q(301,100) < 5*6,
        "geometry_third_jet": 12*Q(201,100)+24*Q(301,100) < 17*6,
        "geometry_fourth_jet": 24*Q(301,100) < 13*6,
        "curvature_vertex_majorant": Q(42,10**6)+Q(49,10**12) < 2,
    }
    k = 2000
    eb = Q(1,36*k**9)
    vb = 81*10**6*eb
    comparisons.update({
        "benchmark_epsilon": eb < Q(6,10**32),
        "benchmark_lapse_1e20": vb < Q(1,10**20),
        "component_target": Q(1,10**20) < Q(1,10**8),
        "bootstrap_F_lower": (1-vb)*(1-2*eb) > Q(99,100),
        "bootstrap_F_upper": (1+6*eb)/(1-vb) < Q(101,100),
        "bootstrap_variation": vb < Q(1,10),
        "benchmark_log80": sum(Q(80)**j/math.factorial(j) for j in range(201)) > 16*k**9*2000,
    })
    hb = Q(1,10**30)
    clock_lower = (1-eb)/(1+hb)**3*(1-vb)
    clock_upper = ((2+hb)*(1+hb+hb*hb/2)/(2*(1-3*eb)))/(1-vb)
    comparisons["temperature_factor_lower"] = clock_lower > Q(9,10)
    comparisons["temperature_factor_upper"] = clock_upper < Q(11,10)
    for label, value in comparisons.items():
        check(value, "rational_"+label)
    return {"epsilon_upper": str(eb), "lapse_upper": str(vb), "comparisons": len(comparisons)}


def benchmark():
    rows = []
    for dps in CONTRACT["precision"]:
        with mp.workdps(dps):
            k = mp.mpf(2000)
            xi = (mp.sqrt(1+12*mp.pi*k**6)-1)/6
            L2 = k**10
            b2 = xi*L2
            gamma = k**4/(12*mp.pi)
            charge = b2*(1+3*xi)
            coeff = 1-6*gamma/L2
            R2 = 2*(charge+gamma)/(mp.sqrt(coeff**2+12*(charge+gamma)/L2)+coeff)
            residual = abs(R2+3*R2**2/L2-charge-gamma*(1+6*R2/L2))/charge
            check(residual < mp.mpf("1e-60"), "benchmark_center_residual_"+str(dps), scaled=mp.nstr(residual,25))
            eps = gamma/R2
            nu = 2*k**4/mp.mpf(2000)*(1+6*xi)/mp.sqrt(b2)**3
            cases = []
            for zi in [mp.mpf("0.25"),mp.mpf("0.28")]:
                target = zi*nu
                hlo = target*mp.sqrt(R2)/(2*(1+6*R2/L2))
                hhi = 4*hlo
                bound = 10**6*eps*(1+mp.log(1/hlo))
                check(bound < mp.mpf("1e-20"), f"benchmark_lapse_dps{dps}_zi{zi}", upper=mp.nstr(bound,30))
                check(hlo >= 1/(mp.mpf(2000)*16*k**9) and mp.log(1/hlo)<80,
                      f"benchmark_horizon_bracket_dps{dps}_zi{zi}")
                cases.append({"zi":str(zi),"h_low":mp.nstr(hlo,40),"h_high":mp.nstr(hhi,40),"lapse_bound":mp.nstr(bound,40)})
            rows.append({"dps":dps,"epsilon":mp.nstr(eps,45),"R_over_b_minus1":mp.nstr(mp.sqrt(R2/b2)-1,45),
                         "root_scaled_residual":mp.nstr(residual,25),"cases":cases})
    return rows


def synthetic_run(htext, level):
    spec = CONTRACT["refinements"][level]
    eps, xi, h = 1e-10, 3., float(htext)
    rh, qstar, cc = 1+h, 1+3*xi, 1+6*xi
    charge = qstar-eps*cc
    with mp.workdps(80):
        em, xm, hm = mp.mpf("1e-10"),mp.mpf(3),mp.mpf(htext)
        rm, qm = 1+hm,1+3*xm-em*(1+6*xm)
        def reference(r):
            d = r-rm
            hh = hm*(2+hm)*(1+3*xm*(rm*rm+1))
            poly = hh/rm+d*(1+6*xm*rm*rm)+4*xm*rm*d*d+xm*d**3
            return d*poly/(r*r)
        def am(r):
            return ((r*r-1)*(1+3*xm*(r*r+1)-6*em*xm))/(r*(r*r-3*em))
        def bm(r):
            return r/(r*r-3*em)
        def sm(r, f):
            return 3*(-f*(r*r+6*em)-3*qm+3*em+r*r+6*em*xm*r*r+xm*r**4)/(r*r*(r*r-3*em)**2)
        def cm(r):
            return (3*mp.diff(reference,r)-6*xm*(r*r-1)/r)/(r*r-3*em)
        ah, bh = am(rm),bm(rm)
        p0 = sm(rm,mp.mpf(0))/(2*ah)
        s1 = mp.diff(lambda r:sm(r,mp.mpf(0)),rm)+mp.diff(lambda f:sm(rm,f),mp.mpf(0))*ah
        p1 = (s1-2*mp.diff(am,rm)*p0+ah*bh*p0)/(3*ah)
        w1 = cm(rm)
        w2 = mp.diff(cm,rm)-bh*w1-p0*ah
        offset = hm*mp.mpf(spec["offset_over_h"])
        initial = np.array([float(p0+p1*offset),float(w1*offset+w2*offset**2/2),
                            float(p0*offset+p1*offset**2/2)])
        kappa = ah/2
        seeds = {"p_h":mp.nstr(p0,30),"p_prime_h":mp.nstr(p1,30),
                 "w_prime_h":mp.nstr(w1,30),"w_second_h":mp.nstr(w2,30)}
    def geometry(t, state):
        delta = h*math.exp(t)
        r = rh+delta
        hs = h*(2+h)*(1+3*xi*(rh*rh+1))
        poly = hs/rh+delta*(1+6*xi*rh*rh)+4*xi*rh*delta**2+xi*delta**3
        pd = (1+6*xi*rh*rh)+8*xi*rh*delta+3*xi*delta**2
        fs = delta*poly/r**2
        fsp = (poly+delta*pd)/r**2-2*delta*poly/r**3
        p, w, integ = state
        f = fs+eps*w
        D = r*r-3*eps
        a = ((h+delta)*(2+h+delta)*(1+3*xi*(r*r+1)-6*eps*xi))/(r*D)
        b = r/D
        Sscaled = 3*(-f*(r*r+6*eps)-3*charge+3*eps+r*r+6*eps*xi*r*r+xi*r**4)/(r*r*D*D)
        wp = (3*fsp-6*xi*(h+delta)*(2+h+delta)/r)/D-b*w-p*f
        pp = (Sscaled-2*a*p)/f+b*p
        return delta,r,fs,f,np.array([pp,wp,p])
    def rhs(t, state):
        delta,_,_,_,derivative = geometry(t,state)
        return delta*derivative
    start = math.log(float(spec["offset_over_h"]))
    finish = math.log((20-rh)/h)
    sol = solve_ivp(rhs,(start,finish),initial,method="DOP853",
                    rtol=float(spec["rtol"]),atol=float(spec["atol"]),dense_output=True)
    check(sol.success,f"synthetic_ODE_success_h{htext}_level{level}",message=sol.message,nfev=sol.nfev)
    gaps = sorted(set([.01*h,.1*h,h,10*h,.01,.1,1.,4.,20-rh]))
    rows = []
    for gap in gaps:
        t = math.log(gap/h)
        state = sol.sol(t)
        delta,r,fs,f,derivative = geometry(t,state)
        with mp.workdps(80):
            rr = rm+hm*mp.exp(mp.mpf(str(t)))
            p,w,integ = [mp.mpf(str(float(v))) for v in state]
            fm = reference(rr)+em*w
            wp = cm(rr)-bm(rr)*w-p*fm
            fp = mp.diff(reference,rr)+em*wp
            Z = fp/2+fm*em*p-fm/rr
            V = 2+6*xm*rr*rr-2*qm/(rr*rr)
            constraint = 2*em*Z*Z+4*rr*fm*Z+6*fm*fm-V*fm-2*em*kappa*kappa*mp.exp(-2*em*integ)
            norm = abs(constraint)/(em*max(1,fm*fm,kappa*kappa))
            check(norm < mp.mpf(spec["tolerance"]),f"synthetic_HH_constraint_h{htext}_level{level}_gap{gap}",
                  raw=mp.nstr(constraint,25),normalized=mp.nstr(norm,25))
        check(f>0 and abs(eps*state[0])<=1e6*eps/(r*r*(r-1)),
              f"synthetic_positive_F_lapse_barrier_h{htext}_level{level}_gap{gap}")
        rows.append({"gap":gap,"r":r,"p":float(state[0]),"w":float(state[1]),
                     "integral_p":float(state[2]),"F":f,"F_reference":fs,
                     "constraint_normalized":mp.nstr(norm,25)})
    print("SYNTHETIC",htext,"LEVEL",level,"NFEV",sol.nfev,
          "MAX_NORMALIZED_CONSTRAINT",max(float(v["constraint_normalized"]) for v in rows),flush=True)
    return {"h":htext,"level":level,"spec":spec,"seeds":seeds,"samples":rows,"nfev":sol.nfev}


def finite_diagnostics():
    results = []
    for h in CONTRACT["synthetic_parameters"]["h"]:
        coarse,refined = synthetic_run(h,0),synthetic_run(h,1)
        for c,r in zip(coarse["samples"],refined["samples"]):
            for key in ["p","w","integral_p"]:
                error = abs(c[key]-r[key])/max(1,abs(r[key]))
                check(error<float(CONTRACT["refinement_tolerance"]),
                      f"synthetic_refinement_h{h}_{key}_gap{r['gap']}",normalized=error)
        results.extend([coarse,refined])
    return results


class Tee:
    def __init__(self, stream, logfile):
        self.stream,self.logfile = stream,logfile
    def write(self,value):
        self.stream.write(value)
        self.logfile.write(value.replace("\r\n","\n").replace("\r","\n"))
        self.logfile.flush()
        return len(value)
    def flush(self):
        self.stream.flush()
        self.logfile.flush()


def save(payload):
    RECEIPT.write_text(json.dumps(clean(payload),indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")


def main():
    snapshots = {name:(ROOT/name).read_bytes() for name in INPUTS}
    hashes = {name:hashlib.sha256(data).hexdigest() for name,data in snapshots.items()}
    history = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs":[]}
    run = {"run_id":datetime.now(timezone.utc).isoformat(),"status":"running","contract":CONTRACT,
           "execution_inputs_sha256":hashes,"command":list(sys.orig_argv),"cwd":str(Path.cwd()),
           "text_policy":"raw SHA256; new source/notes/outputs UTF-8 LF no BOM; historical bytes preserved",
           "versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
                       "sympy":sp.__version__,"mpmath":mp.__version__}}
    history["runs"].append(run)
    save(history)
    with LOG.open("a",encoding="utf-8",newline="\n") as logfile:
        with contextlib.redirect_stdout(Tee(sys.stdout,logfile)),contextlib.redirect_stderr(Tee(sys.stderr,logfile)):
            print("\nRUN",run["run_id"],"\nCOMMAND",json.dumps(run["command"]),"CWD",run["cwd"],flush=True)
            print("INPUT_SHA256",json.dumps(hashes,sort_keys=True),flush=True)
            print("FROZEN_CONTRACT",json.dumps(CONTRACT,sort_keys=True),flush=True)
            try:
                for name in INPUTS[:3]:
                    raw = snapshots[name]
                    raw.decode("utf-8")
                    check(not raw.startswith(b"\xef\xbb\xbf") and b"\r" not in raw,"UTF8_LF_"+name)
                run["bulk"] = exact_bulk()
                run["domain"] = exact_domain_and_counterterm()
                run["inequalities"] = exact_inequalities()
                run["benchmark"] = benchmark()
                run["finite_diagnostics"] = finite_diagnostics()
            except BaseException as error:
                run["exception"] = {"type":type(error).__name__,"message":str(error)}
                traceback.print_exc()
                check(False,"execution_exception",error=str(error))
            final_hashes = {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in INPUTS}
            check(hashes==final_hashes,"execution_inputs_unchanged",final_sha256=final_hashes)
            run["final_inputs_sha256"] = final_hashes
            run["checks"] = CHECKS
            run["failed_checks"] = [item for item in CHECKS if not item["passed"]]
            status = 1 if run["failed_checks"] else 0
            run["terminal_exit_status"] = status
            run["status"] = "failed" if status else "passed"
            run["completed_utc"] = datetime.now(timezone.utc).isoformat()
            if status:
                print("FAILED_INPUT_SNAPSHOTS_BASE64_BEGIN",flush=True)
                for name,raw in snapshots.items():
                    print(json.dumps({"path":name,"sha256":hashes[name],
                                      "base64":base64.b64encode(raw).decode("ascii")}),flush=True)
                print("FAILED_INPUT_SNAPSHOTS_BASE64_END",flush=True)
            print("TERMINAL_EXIT_STATUS",status,"CHECKS",len(CHECKS),"FAILURES",len(run["failed_checks"]),flush=True)
            save(history)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
