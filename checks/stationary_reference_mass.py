"""Exact identities and frozen margins for the S139 mass/horizon milestone.

The uniform heat estimates, convergence hypotheses and localized norm
counterexample remain ordinary arguments. No metric or charge is computed.
"""
from __future__ import annotations

import base64
import contextlib
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback

import sympy as S

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "evidence/stationary-reference-mass.json"
LOG = ROOT / "evidence/stationary-reference-mass.log"
INPUTS = [
    "checks/stationary_reference_mass.py",
    "notes/stationary-reference-finite-mass.md",
    "notes/stationary-horizon-product-matching.md",
    "notes/stationary-reference-higher-mode-clock.md",
    "notes/stationary-reference-nonbolt-stress.md",
    "notes/stationary-weyl-radial-sources.md",
    "notes/stationary-massive-boundary-layer.md",
    "notes/stationary-lll-quantum-exterior.md",
    "notes/uniform-dirac-reference.md",
    "notes/quantum-reference-clock-centering.md",
    "notes/stationary-exterior-solver-obstruction.md",
]
CONTRACT = {
    "acceptance": "all symbolic residuals exactly zero and all rational inequalities true",
    "numerical_tolerances": "none",
    "scope": "fixed classical RN/HH common Dirac source; exact necessary horizon equations",
    "mass_bound": "b integral abs(rho2_ren) dr < 10^16 q^2",
    "ordinary_arguments": ["uniform angular/heat remainder", "conditional product limit",
                           "auxiliary algebraic contraction", "localized horizon-jet norm witness"],
    "excluded": ["corrected exterior existence", "physical boundary Hamiltonian",
                 "metric-response bound", "physical gate completion"],
    "failure_policy": "append every run and log; preserve all failed input bytes",
}
ROWS = []
R = S.Rational


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)


def zero(label, expr):
    expr = S.factor(S.trigsimp(S.together(expr)))
    check(label, expr == 0, expr)


def angular_and_margins():
    x, q, C, tau, E = S.symbols("x q C tau E", positive=True)
    energy = x*(x+q)/C
    zero("angular_degeneracy_identity", 2*C*S.diff(energy,x)-2*(q+2*x))
    zero("angular_variable_jacobian", C**2*S.diff(energy,x)**2-(q*q+4*C*energy))
    g = E*S.sqrt(q*q+4*C*E)*S.exp(-tau*E)/C
    logprime = 1/E+2*C/(q*q+4*C*E)-tau
    zero("angular_log_derivative", S.diff(g,E)/g-logprime)
    zero("angular_unimodality_derivative", S.diff(logprime,E)
         +1/E**2+8*C*C/(q*q+4*C*E)**2)
    zero("continuous_angular_integral", S.integrate(E*S.exp(-tau*E),(E,0,S.oo))-1/tau**2)
    k = S.Integer(2000)
    q0 = k**4
    r0 = R(22,6250)
    Cs, Cb, K, K0, T = map(S.Integer, [10**15,10**16,10**15,10**11,4*10**15])
    check("xi_lower_bound_from_pi", 1+12*R(31,10)*k**6 > (1+6*k**3)**2)
    check("angular_ratio_bound", 1+12*R(22,7)*k**6 < r0**2*q0**2)
    check("proper_gap_C_bound", 1+12*R(22,7)*k**6 < (R(123,20)*k**3)**2)
    check("proper_gap_exceeds_18", k/R(123,20)>18**2)
    smooth = 8*Cs*r0**3
    long = 1400*10**10*r0
    local = r0*(R(3,8)+R(5,600))+500*r0**2
    check("short_smooth_gradient", smooth < 349000000)
    check("long_even_gradient", long <= 49280000000)
    check("local_restoration_gradient", local < 1)
    check("actual_cap_gradient", r0/5 < 1)
    check("universal_cubic_long_tail", 200*r0**2 < 1)
    check("smooth_long_restoration_ledger", smooth+long+3 < 5*10**10)
    check("angular_discrepancy_gradient", 12*r0+24/S.sqrt(q0)<1)
    check("log_trace_gradient_ledger", 96*Cb*r0**2+5*10**10+1 < 2*10**13)
    check("outer_geometry_quartic", R(1,1000)+R(402,100)*R(9,20)+R(301,100)/5<R(5,2))
    check("outer_Y4_less_than_one", R(8,9)<1)
    check("trace_anchor_to_log_tail_constant", K0+2*K+2<3*10**15)
    check("declared_trace_tail_coefficients", 3*10**15+1<T and 2*10**13<T)
    check("log_four_less_than_three_halves", sum(R(3,2)**j/S.factorial(j) for j in range(4))>4)
    check("outer_F_lower_bound", (1-R(1001,1000)/4)*(1-R(1,4))>R(1,2))
    check("outer_Bprime_bound", R(2,1000)+R(3,2)*R(201,100)+R(301,100)/4<4)
    check("pressure_integral_budget", K0+2*K+R(1,15)+R(3,8)*T<4*10**15)
    check("outer_trace_integral_budget", T*(2+R(3,2))/4<4*10**15)
    check("Weyl_log_coefficient", 432+288+864+432+60==2076)
    check("Weyl_constant_coefficient", 593+338+1408+504+588+288+11==3730)
    check("Weyl_remainder_majorant", 2076<2*2160 and 3730<2*2160)
    check("LLL_integral_budget", 3+R(3,4)<4 and 4*r0<1)
    mass_total = 3*(K+K0/2)+1+10**15+4*10**15+4+1
    check("complete_absolute_mass_budget", mass_total<9*10**15<10**16)
    check("normalized_geometric_mass", R(10**16,6*k**8)<R(6511,10**15))
    zero("unscaled_geometric_mass", R(10**16,k**5)-R(5,16))


def mass_identities():
    t, am, xi, Q, q = S.symbols("t a_m xi Q q", positive=True)
    F = 1+xi*t*t-2*am/t+Q/t**2
    B = F/t**2
    Bp = -2/t**3+6*am/t**4-4*Q/t**5
    zero("exact_Bprime",S.diff(B,t)-Bp)
    zero("Bprime_concavity",S.diff(t**5*Bp/2,t,2)+2)
    zero("B_four_extremal_upper",B.subs({t:4,am:1+2*xi,Q:1+3*xi})-xi-(9-13*xi)/256)
    zero("Bprime_four_extremal_lower",(t**5*Bp/2).subs({t:4,am:1+2*xi,Q:1+3*xi})-(-6+18*xi))
    zero("linear_trace_pressure_primitive",S.diff(2/t-3*am/t**2+4*Q/(3*t**3),t)-t*Bp)
    zero("log_pressure_integral",4*S.integrate((1+S.log(t))/t**3,(t,4,S.oo))-(3+2*S.log(4))/16)
    zero("log_density_integral",S.integrate((1+S.log(t))/t**2,(t,4,S.oo))-(2+S.log(4))/4)
    zero("global_Weyl_remainder_integral",S.integrate(2*(1+S.log(t))/t**2,(t,1,S.oo))-4)
    zero("LLL_tail_integral",S.integrate(4/t**3,(t,4,S.oo))-R(1,8))
    Alog = 432*Q**2+288*Q*xi*t**4-864*Q*am*t+432*Q*t**2-60*q*q*t**4
    Aconst = -593*Q**2-338*Q*xi*t**4+1408*Q*am*t-504*Q*t**2-588*am**2*t**2+288*am*t**3-11*t**4
    # Original frozen Weyl Euler numerator, before removing its two outer terms.
    Ealpha = (Alog*S.log(t)+Aconst-33*xi**2*t**8+144*xi*am*t**5)/(180*t**6)
    V = 11*xi**2*t**2/(240*S.pi)
    L = -am*xi/(5*S.pi)
    RW = -(Alog*S.log(t)+Aconst)/(720*S.pi*t**6)
    zero("exact_Weyl_vacuum_and_linear_split",-Ealpha/(4*S.pi)-V-L/t-RW)
    I, tail = S.symbols("I trace_remainder")
    zero("massive_Weyl_cancellation_inside_norm",I/(2*F)-(L*t+tail)/t**2
         +(-Ealpha/(4*S.pi))-V-(I/(2*F)-tail/t**2+RW))
    C, u = S.symbols("C u",positive=True)
    c = -am*S.sqrt(xi)/C**R(3,2)
    zero("universal_cubic_to_linear_trace",c*C*C/(5*S.pi*S.sqrt(C/xi)*u)-L/u)
    r, rh, G, Av, lam0, mu_h = S.symbols("r R G A_v lambda0 mu_h",positive=True)
    rho = S.Function("rho_ren")(r)
    m = S.Function("m")(r)
    mu = m-G*Av*r**3/3
    zero("vacuum_mass_derivative",S.diff(mu,r).subs(S.diff(m,r),G*(rho+Av*r*r))-G*rho)
    lam_eff = lam0-2*G*Av/3
    zero("vacuum_geometry_allocation",1+lam0*r*r-2*m/r+Q/r**2
         -(1+lam_eff*r*r-2*mu/r+Q/r**2))
    mh = rh*(1+lam0*rh**2+Q/rh**2)/2
    zero("vacuum_horizon_datum",mh-G*Av*rh**3/3-rh*(1+lam_eff*rh**2+Q/rh**2)/2)


def horizon_identities():
    # Compute Einstein components directly from the diagonal metric.
    U,r,th,ph = S.symbols("U r theta phi",real=True)
    coords = (U,r,th,ph)
    F, a = S.Function("F")(r), S.Function("alpha")(r)
    diag = (-S.exp(2*a)*F,1/F,r*r,r*r*S.sin(th)**2)
    metric, inv = S.diag(*diag), S.diag(*(1/x for x in diag))
    gamma = {}
    for k in range(4):
        for i in range(4):
            for j in range(4):
                gamma[k,i,j] = S.simplify(sum(inv[k,l]*(S.diff(metric[l,j],coords[i])
                    +S.diff(metric[l,i],coords[j])-S.diff(metric[i,j],coords[l])) for l in range(4))/2)
    ricci = []
    for i in range(4):
        ricci.append(S.simplify(sum(S.diff(gamma[k,i,i],coords[k])-S.diff(gamma[k,i,k],coords[i])
            +sum(gamma[k,i,i]*gamma[l,k,l]-gamma[l,i,k]*gamma[k,i,l] for l in range(4)) for k in range(4))))
    scalar = S.simplify(sum(inv[i,i]*ricci[i] for i in range(4)))
    mixed = [S.simplify(inv[i,i]*ricci[i]-scalar/2) for i in range(4)]
    zero("direct_orbit_Einstein_component",mixed[0]-(S.diff(F,r)/r+(F-1)/r**2))
    zero("direct_orbit_Einstein_difference",mixed[1]-mixed[0]-2*F*S.diff(a,r)/r)
    angular = S.diff(F,r,2)/2+R(3,2)*S.diff(F,r)*S.diff(a,r)+S.diff(F,r)/r
    zero("direct_angular_horizon_component",(mixed[2]-angular).subs(F,0))
    R0,f1,f2,A,eps = S.symbols("R0 f1 f2 A epsilon",positive=True)
    x = r-R0
    poly = f1*x+f2*x*x/2+eps*A*x**3
    sc = scalar.subs({F:poly,a:0}).doit()
    box = poly*S.diff(sc,r,2)+(S.diff(poly,r)+2*poly/r)*S.diff(sc,r)
    zero("cubic_metric_keeps_horizon_scalar",S.diff(sc,eps).subs({r:R0,eps:0}))
    zero("cubic_metric_radial_curvature_derivative",S.diff(sc,eps,r).subs({r:R0,eps:0})+6*A)
    zero("cubic_metric_horizon_box_variation",S.diff(box,eps).subs({r:R0,eps:0})+6*f1*A)
    G = S.symbols("G",positive=True)
    zero("retained_horizon_anomaly_witness",8*S.pi*G/(30*16*S.pi**2)*(-6*f1*A)+G*f1*A/(10*S.pi))
    tol = S.symbols("epsilon_defect",positive=True)
    zero("zero_base_jet_size_for_pressure_defects",(G*f1*A/(10*S.pi)-4*tol*f1/R0)
         *(10*S.pi*R0/f1)-(G*R0*A-40*S.pi*tol))
    zero("two_allowed_metrics_jet_difference",(G*f1*A/(10*S.pi)-8*tol*f1/R0)
         *(10*S.pi*R0/f1)-(G*R0*A-80*S.pi*tol))
    lam,m,Q,kap,rhoh = S.symbols("lambda m Q kappa rho_h",positive=True)
    rn = 1+lam*r*r-2*m/r+Q/r**2
    mh = R0*(1+lam*R0**2+Q/R0**2)/2
    Qtemp = R0**2+3*lam*R0**4-2*kap*R0**3
    D = (S.diff(rn,r,2)/2+kap/r).subs({r:R0,m:mh}).subs(Q,Qtemp)
    zero("finite_temperature_radius_denominator",D-(1+6*lam*R0**2)/R0**2+3*kap/R0)
    zero("horizon_mass_first_variation",S.diff(mh,R0).subs(Q,Qtemp)-kap*R0)
    Fprime_h = S.diff(1+lam*r*r-2*S.Function("m")(r)/r+Q/r**2,r)
    Fprime_h = Fprime_h.subs(S.diff(S.Function("m")(r),r),G*rhoh).subs(S.Function("m")(r),mh).subs(r,R0)
    zero("source_dependent_horizon_clock",Fprime_h-(R0**2+3*lam*R0**4-Q-2*G*R0**2*rhoh)/R0**3)
    z,u,alpha,e,Aprod,Sprod,a2,fT,eA,eS = S.symbols("z u alpha epsilon Aprod Sprod a2 fT eA eS")
    EA = alpha*z+(1-alpha)*z*z/2-(1+alpha)/2+e*Aprod
    ES = z*z/u-(1-alpha)*z*z/2-(1+alpha)/2-e*Sprod
    Tz = 1-e*Aprod-(1-alpha)*(z-1)**2/2
    den = 1-alpha+alpha/z-e*a2/(2*z*z)
    zero("forced_product_z_map",EA-(z-Tz))
    zero("forced_product_u_map",(ES-EA).subs(Sprod,-a2/2-Aprod)
         -(z*z/u-z*z*den))
    zero("forced_product_u_rhs",(-fT+eS)-(fT-eA)-(-2*fT+eA+eS))
    check("homogeneous_inverse_bound",1/(1-R(1,50))==R(50,49))
    check("auxiliary_thermal_root_displacement",R(21,10)/(1-R(1,50))<R(22,10))
    check("auxiliary_positive_z_gap",R(1,50)*R(22,10)<R(45,1000))
    d = R(1,10**20)
    check("auxiliary_gap_lower_margin",R(955,1000)*(1-d)**3/(1+d)**2>R(9,10))
    check("auxiliary_gap_upper_margin",R(1045,1000)*(1+d)**3/(1-d)**2<R(11,10))
    thermal_max = 4*R(7,25)*R(1,2000)*3*R(102,100)/S.Integer(2000)**9
    check("auxiliary_thermal_size_at_bounded_alpha",thermal_max<R(1,10**30))


def hashes(data):
    return {name:hashlib.sha256(raw).hexdigest() for name,raw in data.items()}


def main():
    before = {name:(ROOT/name).read_bytes() for name in INPUTS}
    existing = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs":[]}
    run = {"started_utc":datetime.now(timezone.utc).isoformat(),"contract":CONTRACT,
           "command":sys.orig_argv,"dont_write_bytecode":sys.dont_write_bytecode,"cwd":str(Path.cwd()),
           "python":platform.python_version(),"sympy":S.__version__,"execution_inputs_sha256":hashes(before)}
    with LOG.open("a",encoding="utf-8",newline="\n") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print("RUN "+json.dumps(run,sort_keys=True),flush=True)
        try:
            check("actual_bytecode_disabled",sys.dont_write_bytecode and "-B" in sys.orig_argv)
            angular_and_margins()
            mass_identities()
            horizon_identities()
        except Exception:
            traceback.print_exc()
            check("uncaught_exception",False)
        run["final_inputs_sha256"] = hashes({name:(ROOT/name).read_bytes() for name in INPUTS})
        print("FINAL_INPUTS "+json.dumps(run["final_inputs_sha256"],sort_keys=True),flush=True)
        check("frozen_inputs_unchanged",run["execution_inputs_sha256"]==run["final_inputs_sha256"])
        run["checks"] = ROWS
        run["failed_checks"] = [row for row in ROWS if not row["passed"]]
        code = int(bool(run["failed_checks"]))
        run["status"] = "failed" if code else "passed"
        run["terminal_exit_status"] = code
        run["finished_utc"] = datetime.now(timezone.utc).isoformat()
        if code:
            for name,raw in before.items():
                print(json.dumps({"path":name,"sha256":hashlib.sha256(raw).hexdigest(),
                                  "base64":base64.b64encode(raw).decode("ascii")},sort_keys=True))
        print("TERMINAL "+json.dumps({"status":run["status"],"checks":len(ROWS),
              "failures":len(run["failed_checks"]),"terminal_exit_status":code,
              "finished_utc":run["finished_utc"]},sort_keys=True),flush=True)
    existing["runs"].append(run)
    RECEIPT.write_text(json.dumps(existing,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":run["status"],"checks":len(ROWS),"failures":len(run["failed_checks"]),
                      "receipt":str(RECEIPT)},sort_keys=True))
    return code


if __name__ == "__main__":
    sys.exit(main())
