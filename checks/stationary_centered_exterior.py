"""Focused exact checks for the complete centered source and linear geometry.

Analytic heat/domain estimates are ordinary arguments, not numerical tests.
This contract does not certify metric response or a full-equation residual.
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
RECEIPT = ROOT / "evidence/stationary-centered-exterior.json"
LOG = ROOT / "evidence/stationary-centered-exterior.log"
INPUTS = [
    "checks/stationary_centered_exterior.py",
    "notes/stationary-centered-exterior-source.md",
    "notes/stationary-centered-first-order-geometry.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-horizon-source-defect.md",
    "notes/stationary-reference-cold-bolt.md",
    "notes/stationary-reference-nonbolt-stress.md",
    "notes/stationary-reference-finite-mass.md",
    "notes/stationary-weyl-radial-sources.md",
    "notes/stationary-horizon-product-matching.md",
    "notes/spherical-anomaly-action.md",
    "notes/uniform-dirac-reference.md",
]
CONTRACT = {
    "acceptance": "all symbolic residuals exactly zero and all rational inequalities true",
    "numerical_tolerances": "none",
    "scope": "complete fixed-seed source, prepared-temperature linear geometry and exact regular comparison completion",
    "ordinary_arguments": ["heat and common-subtraction transfer", "actual paired-domain localization",
                           "global integral and curvature bounds", "cold maximum-principle barrier"],
    "excluded": ["source metric/state response", "full retained-action residual bound",
                 "controlled spatial reduction", "corrected exterior existence",
                 "full boundary Hamiltonian", "physical gate closure"],
    "failure_policy": "append every run/log and preserve failed-run input bytes",
}
ROWS = []
R = S.Rational


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)


def zero(label, expr):
    value = S.factor(S.together(expr))
    check(label, value == 0, value)


def seed_and_source():
    t, u, x, d, delta = S.symbols("t u x d delta", positive=True)
    C = 1+6*x+d
    Q0 = 1+3*x
    a = 1+2*x+delta/2
    f = 1+x*t*t-2*a/t+(Q0+d)/t**2-2*d/t**3+d/t**4
    b = x+u*u-2*a*u**3+Q0*u**4+d*u**4*(1-u)**2
    be = b.subs(delta, 0)
    p = u**4*(1-u)**2
    zero("compactified_seed", u*u*f.subs(t,1/u)-b)
    zero("extremal_compact_factorization", be-(1-u)**2*(x*(1+2*u+3*u*u)+u*u+d*u**4))
    zero("compact_derivative_single_maximum", S.diff(be,u)-2*u*(1-u)*(1-2*Q0*u+d*u*u*(2-3*u)))
    zero("compact_second_derivative_gap", 2*C-S.diff(be,u,2)-12*(1-u)*(x+Q0*u)-d*(2-S.diff(p,u,2)))
    zero("compact_fifth_derivative", S.diff(b,u,5)-d*(-240+720*u))
    zero("compact_sixth_derivative", S.diff(b,u,6)-720*d)
    zero("monotone_seed_numerator", t**3*S.diff(f,t)/2-(x*t**4+a*t-Q0+d*(-1+3/t-2/t**2)))
    zero("strong_convexity_x_part", S.diff(f,t,2).diff(x)-2*(1-4/t**3+9/t**4))
    bracket = 1-4/t**3+9/t**4
    zero("convexity_unique_minimum", S.diff(bracket,t)-12*(t-3)/t**5)
    check("convexity_minimum_value", bracket.subs(t,3)==R(26,27))
    zero("pressure_integral_factor", S.diff(f/t**2,t)-(-2/t**3+6*a/t**4-4*(Q0+d)/t**5+10*d/t**6-6*d/t**7))
    zero("pressure_integral_primitive", S.diff(2/t-3*a/t**2+4*(Q0+d)/(3*t**3)-5*d/(2*t**4)+6*d/(5*t**5),t)-t*S.diff(f/t**2,t))
    F = S.Function("F")(t)
    th = S.Function("Theta")(t)
    j = S.Function("J")(t)
    pressure = th/(2*t*t)-j/(2*F)
    density = -th/(2*t*t)-j/(2*F)
    rules = {S.diff(j,t):F*S.diff(th,t)/(t*t)}
    zero("canonical_density_derivative", S.diff(density,t).subs(rules)
         -(th/t**3-S.diff(th,t)/t**2+S.diff(F,t)*j/(2*F**2)))
    zero("canonical_pressure_derivative", S.diff(pressure,t).subs(rules)
         -(-th/t**3+S.diff(F,t)*j/(2*F**2)))
    zero("canonical_trace_identity", pressure-density-th/t**2)
    zero("canonical_HH_sum", pressure+density+j/F)
    # Exact local heat identity using the compact orbit Laplacian.
    Bn = S.Function("Bn")(u)
    curv = -S.diff(Bn,u,2)
    boxcurv = S.diff(Bn*S.diff(curv,u),u)
    zero("common_local_polynomial_identity", S.diff(curv**2+2*boxcurv,u)
         +4*S.diff(Bn,u)*S.diff(Bn,u,4)+2*Bn*S.diff(Bn,u,5))
    r4 = -S.diff(f,t,2)-4*S.diff(f,t)/t+2*(1-f)/t**2
    zero("centered_four_scalar", r4-(-12*x+4*d/t**5-6*d/t**6))
    pa, pv, dp = S.symbols("pa pv dp")
    ref = -4*S.pi*((pa-pv)/t**2+dp*(-t**-2+4*t**-3-3*t**-4))
    zero("reference_finite_mass_integral", S.integrate(ref,(t,1,S.oo))+4*S.pi*(pa-pv))


def local_weyl():
    t, a, qc, x, d, q = S.symbols("t a qc x d q", positive=True)
    F = S.Function("F")(t)
    rb = -t*t*S.diff(F,t,2)+2*t*S.diff(F,t)-2*F
    tt = t*S.diff(F,t)
    rr, tr, fv, sig = S.symbols("Rbar T Fvalue sigma")
    H = sig*(rr*rr/60-rr/18+R(1,15)-q*q/3)+rr*tr/30-tr*tr/10+11*tr*fv/90-11*fv*fv/180-fv/18
    substitutions = {rr:rb,tr:tt,fv:F,sig:S.log(t)}
    L = H.subs(substitutions)/t**2
    EF = sum((-1)**n*S.diff(S.diff(L,S.diff(F,t,n)),t,n) for n in range(3))
    expected = (S.log(t)*t*t*S.diff(F,t,4)+(4*S.log(t)+2)*t*S.diff(F,t,3)
                +2*S.diff(F,t,2)-2*S.diff(F,t)/t)/30
    zero("Weyl_F_Euler_direct_variation", EF-expected)
    f = 1+x*t*t-2*a/t+qc/t**2-2*d/t**3+d/t**4
    ef = expected.subs(F,f).doit()
    alpha = (a/t**3+qc*(R(3,2)*S.log(t)-R(13,8))/t**4)/(15*S.pi)
    alpha += d*((R(22,25)-R(8,5)*S.log(t))/t**5+(2*S.log(t)-R(11,15))/t**6)/(2*S.pi)
    zero("integrated_complete_Weyl_lapse", S.diff(alpha,t)+ef/(2*S.pi*t))
    zero("Weyl_EF_explicit_RN_and_seed_terms", ef-2*(6*qc*S.log(t)-8*qc+3*a*t)/(15*t**4)
         -d*((6-8*S.log(t))/t**5+(12*S.log(t)-R(32,5))/t**6))
    # The lapse variation is independent of the preceding F variation.
    al = S.Function("al")(t)
    Tfull = t*(S.diff(F,t)+F*S.diff(al,t))
    Rfull = -t*t*S.diff(F,t,2)-3*t*t*S.diff(al,t)*S.diff(F,t) \
        -2*t*t*F*(S.diff(al,t,2)+S.diff(al,t)**2)+2*Tfull-2*F
    Lfull = S.exp(al)*H.subs({rr:Rfull,tr:Tfull,fv:F,sig:S.log(t)})/t**2
    off = {al:0,S.diff(al,t):0,S.diff(al,t,2):0}
    EA = sum((-1)**n*S.diff(S.diff(Lfull,S.diff(al,t,n)).subs(off,simultaneous=True),t,n) for n in range(3))
    HR = S.diff(H,rr).subs(substitutions)
    HT = S.diff(H,tr).subs(substitutions)
    identity = H.subs(substitutions)-rb*HR-(t*t*S.diff(F,t)+2*t*F)*S.diff(HR,t) \
        -2*t*t*F*S.diff(HR,t,2)-t*F*S.diff(HT,t)-(t*S.diff(F,t)-F)*HT
    zero("Weyl_full_lapse_Euler_identity", t*t*EA-identity)


def geometry():
    r, kap, sr, G, rh, ah, fpp = S.symbols("r kappa sR G rhoh ah fpp")
    horizon = S.symbols("R", positive=True)
    D = fpp/2+kap/horizon
    svalue = (G*rh/horizon-kap*ah)/D
    fh = -2*kap*sr
    fph = -2*G*rh/horizon-fh/horizon
    zero("temperature_selects_unique_displacement", (fph+fpp*sr+2*kap*ah).subs(sr,svalue))
    F = S.Function("F")(r)
    f = S.Function("f")(r)
    a = S.Function("a")(r)
    J = S.Function("J")(r)
    H = a+(f+sr*S.diff(F,r))/F
    xi = sr-F*J
    uprimeform = a+f/(2*F)+(S.diff(F,r)/(2*F)-1/r)*xi
    uform = a+(f+sr*S.diff(F,r))/(2*F)-sr/r-(S.diff(F,r)/2-F/r)*J
    zero("common_disk_conformal_variation", uprimeform-uform)
    zero("common_disk_displacement_ODE", S.diff(xi,r).subs(S.diff(J,r),-H/F)
         -S.diff(F,r)*xi/F-a-f/F)
    zero("regular_H_controls_displacement", S.diff(xi,r).subs(S.diff(J,r),-H/F)-H+S.diff(F,r)*J)
    rho = S.Function("rho")(r)
    p = S.Function("p")(r)
    AA = S.Function("Arest")(r)
    frule = -2*G*rho/r-f/r
    arule = G*(rho+p)/(r*F)
    f2rule = S.diff(frule,r).subs(S.diff(f,r),frule)
    a2rule = S.diff(arule,r)
    pS = (AA-(p-rho)/(4*S.pi*r*r))/2
    # The joint Ward identity supplies the angular equation.
    prule = 8*S.pi*r*pS-(rho+p)*S.diff(F,r)/(2*F)
    dR2 = -S.diff(f,r,2)-3*S.diff(F,r)*S.diff(a,r)-2*F*S.diff(a,r,2)
    subs = {S.diff(f,r,2):f2rule,S.diff(a,r,2):a2rule,S.diff(a,r):arule}
    dR2 = S.expand(dR2.subs(subs,simultaneous=True)).subs(S.diff(p,r),prule)
    zero("invariant_orbit_curvature_from_constraints", dR2-(-16*S.pi*G*pS+2*(-f+G*(p-rho))/r**2))
    dR4 = dR2-4*(S.diff(f,r)+F*S.diff(a,r))/r-2*f/r**2
    zero("full_linear_trace_from_constraints", dR4.subs({S.diff(f,r):frule,S.diff(a,r):arule})+8*S.pi*G*AA)
    dRbar = r*r*dR2+2*r*(S.diff(f,r)+F*S.diff(a,r))-2*f
    zero("product_frame_curvature_from_constraints", dRbar.subs({S.diff(f,r):frule,S.diff(a,r):arule})
         -(-8*S.pi*G*r*r*AA-6*f+6*G*(p-rho)))
    chi, Dc, h = S.symbols("chi Dc h", positive=True)
    W = h*((Dc+1)*S.cosh(chi)-chi*S.sinh(chi))
    zero("cold_barrier_exact_operator", S.trigsimp(-S.diff(W,chi,2)-S.coth(chi)*S.diff(W,chi)+2*W
         -h*(3*S.cosh(chi)+chi/S.sinh(chi))))
    lam, B = S.symbols("lambda B", positive=True)
    V = S.Function("V")(r)
    rlam = B+(r-B)*S.exp(lam*V)
    zero("completion_radius_tangent", S.diff(rlam,lam).subs(lam,0)-(r-B)*V)
    zero("completion_radius_monotonicity", S.diff(rlam,r)-S.exp(lam*V)*(1+lam*(r-B)*S.diff(V,r)))
    uu, R0, dR = S.symbols("u R0 dR")
    completed = S.exp(-2*lam*uu)*((1+2*lam*uu)*R0+lam*dR)
    zero("exact_conformal_curvature_tangent", S.diff(completed,lam).subs(lam,0)-dR)
    zero("exact_conformal_curvature_second_coefficient", S.diff(completed,lam,2).subs(lam,0)/2+2*uu*(uu*R0+dR))
    z, u3, v3 = S.symbols("z u3 v3")
    # z=1/r. Keeping B in the exact radius completion checks that its
    # constant center does not corrupt the cubic outer clock cancellation.
    ratio = B*z+(1-B*z)*S.exp(lam*v3*z**3)
    lapse = 2*S.log(ratio)+2*lam*u3*z**3-lam*v3*z**3 \
        -S.log(1-3*lam*(1-B*z)*v3*z**3)
    cubic = S.series(lapse,z,0,4).removeO().coeff(z,3)
    zero("exact_completion_outer_cubic", cubic-lam*(2*u3+4*v3))
    zero("outer_mass_clock_correlation", cubic.subs(u3,-2*v3))


def ledgers():
    k = S.Integer(2000)
    rat = R(22,6250)
    K = S.Integer(10)**15
    K0 = S.Integer(10)**11
    sig = R(51,10**13)
    L = R(33,10**11)
    Qr = Qp = R(13,10**11)
    hh = R(1,10**30)
    check("benchmark_convexity_guard", 54<R(1,100)*R(99,100)*k**3)
    check("nonbolt_angular_ratio_guard", R(123,20*k)*(1+R(1,10**21))<R(35,10000))
    check("massive_lapse_allocation", 80*R(1001,1000)/R(99,100)*K/k**8<R(32,10**11))
    check("complete_lapse_allocation", R(32,10**11)+R(3,10**28)+2/k**10<L)
    check("mass_source_absolute_ledger", 3*(K+K0/2)+4*K+K+10**4+7+4*R(22,7)*(3301+12)<10**16)
    check("density_seminorm", (31*K+10**6)*R(1001,1000)/k**8<Qr)
    check("pressure_seminorm", (29*K+10**6)*R(1001,1000)/k**8<Qp)
    check("reference_mass_profile", 4*R(22,7)*(3301+12)<50000)
    check("Weyl_polynomial_quadratic_ledger", 16+60+30+64+8+10<200)
    check("Weyl_extra_energy_ledger", (200*(2*9*3+9)+5*3)/(4*3)<10000)
    check("geometric_mass_source", 10**16*R(1001,1000)/(6*k**8)<R(7,10**12))
    check("prepared_temperature_gap", (sig+L)/(1-7*hh)*(1+12*hh)<R(34,10**11))
    sRoverB = 2*hh*(sig+L)
    check("moving_bolt_F_norm", L+3*Qr+80*sRoverB<R(75,10**11))
    H = L+R(75,10**11)
    check("common_disk_H_norm", H==R(108,10**11))
    check("common_disk_radius_norm", R(34,10**11)+10*H<R(12,10**9))
    check("common_disk_radius_monotone", 1-R(3,10**7)>0)
    check("inner_horizon_lapse_difference", (Qr+Qp)/10<R(3,10**11))
    check("cold_seam_radius_norm", R(34,10**11)*100*hh+R(175,100)*H<R(19,10**10))
    check("anomaly_derivative_coefficient", R(14884,20)+R(11*8832,360)+R(123708,30)+R(4,3)<5166)
    check("anomaly_derivative_full_budget", R(5166,144)+56000<100000)
    check("cold_curvature_Lie_budget", R(45,10**11)+R(12,10**9)*R(81,10)<R(1,10**7))
    check("cold_potential_barrier_slack", 3-20*R(15,1000)==R(27,10))
    check("cold_forcing_barrier", R(56,10**9)<R(27,10)*R(1,10**6))
    check("cold_seam_barrier", R(3,10**9)<R(1,10**6)*R(1,100))
    check("cold_horizon_conformal_ideal", R(1,10**6)*83==R(83,10**6))
    # The nonlinear horizon curvature coefficient follows from the exact
    # identity, with |Rhat0| A^2<3 and |deltaRhat| A^2<1e-7 h.
    uh = R(83,10**6)
    check("completion_horizon_quadratic", 3*uh**2*R(21,10)+3*uh*R(1,10**7)<R(1,10**7))
    check("completion_exponential_gap", R(34,10**11)/(1-R(34,10**11))<R(35,10**11))
    check("outer_log_tail_mass_integrable", S.integrate((1+S.log(S.Symbol('t',positive=True)))/S.Symbol('t',positive=True)**2,
          (S.Symbol('t',positive=True),1,S.oo))==2)


def hashes(data):
    return {name: hashlib.sha256(raw).hexdigest() for name, raw in data.items()}


def main():
    before = {name: (ROOT/name).read_bytes() for name in INPUTS}
    existing = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs": []}
    run = {"started_utc": datetime.now(timezone.utc).isoformat(), "contract": CONTRACT,
           "command": sys.orig_argv, "dont_write_bytecode": sys.dont_write_bytecode, "cwd": str(Path.cwd()),
           "python": platform.python_version(), "sympy": S.__version__, "execution_inputs_sha256": hashes(before)}
    with LOG.open("a", encoding="utf-8", newline="\n") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print("RUN "+json.dumps(run, sort_keys=True), flush=True)
        try:
            check("actual_bytecode_disabled", sys.dont_write_bytecode and "-B" in sys.orig_argv)
            check("inputs_strict_UTF8_LF", all(raw.decode("utf-8") is not None and b"\r" not in raw
                  and not raw.startswith(b"\xef\xbb\xbf") for raw in before.values()))
            seed_and_source()
            local_weyl()
            geometry()
            ledgers()
        except Exception:
            traceback.print_exc()
            check("uncaught_exception", False)
        run["final_inputs_sha256"] = hashes({name: (ROOT/name).read_bytes() for name in INPUTS})
        print("FINAL_INPUTS "+json.dumps(run["final_inputs_sha256"], sort_keys=True), flush=True)
        check("frozen_inputs_unchanged", run["execution_inputs_sha256"] == run["final_inputs_sha256"])
        run["checks"] = ROWS
        run["failed_checks"] = [row for row in ROWS if not row["passed"]]
        code = int(bool(run["failed_checks"]))
        run["status"] = "failed" if code else "passed"
        run["terminal_exit_status"] = code
        run["finished_utc"] = datetime.now(timezone.utc).isoformat()
        if code:
            run["failed_inputs_base64"] = {name: base64.b64encode(raw).decode("ascii") for name, raw in before.items()}
            for name, raw in before.items():
                print(json.dumps({"path": name, "sha256": hashlib.sha256(raw).hexdigest(),
                                  "base64": base64.b64encode(raw).decode("ascii")}, sort_keys=True))
        print("TERMINAL "+json.dumps({"status": run["status"], "checks": len(ROWS),
              "failures": len(run["failed_checks"]), "terminal_exit_status": code,
              "finished_utc": run["finished_utc"]}, sort_keys=True), flush=True)
    existing["runs"].append(run)
    RECEIPT.write_text(json.dumps(existing, indent=2, sort_keys=True)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS), "failures": len(run["failed_checks"]),
                      "receipt": str(RECEIPT)}, sort_keys=True))
    return code


if __name__ == "__main__":
    sys.exit(main())
