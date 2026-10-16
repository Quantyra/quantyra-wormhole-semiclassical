"""Exact algebra and finite budgets for the complete candidate residual.

Heat, domain and regularity estimates are reviewed ordinary arguments.
This checker does not certify them formally or prove an exact solution.
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
RECEIPT = ROOT / "evidence/stationary-centered-complete-residual.json"
LOG = ROOT / "evidence/stationary-centered-complete-residual.log"
INPUTS = [
    "checks/stationary_centered_complete_residual.py",
    "notes/stationary-centered-constraint-residual.md",
    "notes/stationary-centered-conormal-source.md",
    "notes/stationary-centered-local-response.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-complete-residual.md",
    "notes/stationary-centered-exterior-source.md",
    "notes/stationary-centered-first-order-geometry.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-reference-cold-bolt.md",
    "notes/stationary-reference-nonbolt-stress.md",
    "notes/stationary-reference-finite-mass.md",
    "notes/stationary-weyl-radial-sources.md",
    "notes/spherical-anomaly-sources.md",
    "notes/uniform-dirac-reference.md",
]
CONTRACT = {
    "acceptance": "all displayed symbolic residuals exactly zero and rational inequalities true",
    "numerical_tolerances": "none",
    "scope": "algebra and finite ledgers of the specified complete retained-equation defect",
    "ordinary_arguments": ["heat and actual-domain localization", "real C5 and weighted outer transfer",
                           "complete source response", "radial and horizon integration estimates"],
    "excluded": ["formal proof of analytic estimates", "full higher-order inverse",
                 "proximity to an exact stationary solution", "physical state matching",
                 "boundary Hamiltonian", "physical gate closure"],
    "failure_policy": "append every run/log and retain failed-run input bytes",
}
ROWS = []
R = S.Rational


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)


def zero(label, expression):
    value = S.simplify(S.factor(S.together(expression)))
    check(label, value == 0, value)


def exact_identities():
    ell, tau, V, b, u = S.symbols("ell tau V b u", real=True)
    T = 1+tau*(S.exp(ell*V)-1)
    Fratio = S.exp(2*ell*(V-u))*(1+ell*b)**2/T**2
    radial = S.exp(ell*V)*(1+ell*b)
    jm = S.exp(ell*V)*(1+ell*b)/T**2
    ja = S.exp(ell*(2*u-V))/(T*(1+ell*b))
    zero("hat_mass_Jacobian", radial/T**2-jm)
    zero("hat_lapse_Jacobian", radial/(T**3*Fratio)-ja)
    alpha = 2*S.log(T)+2*ell*u-ell*V-S.log(1+ell*b)
    a = (2*tau-1)*V+2*u-b
    remainder = 2*(S.log(T)-ell*tau*V)+ell*b-S.log(1+ell*b)
    zero("exact_lapse_remainder", alpha-ell*a-remainder)
    zero("lapse_linear_tangent", S.diff(alpha,ell).subs(ell,0)-a)
    zero("lapse_quadratic_coefficient", S.diff(remainder,ell,2).subs(ell,0)/2
         -(tau*(1-tau)*V**2+b**2/2))
    zero("source_mass_Jacobian_tangent", S.diff(jm,ell).subs(ell,0)-((1-2*tau)*V+b))
    zero("source_lapse_Jacobian_tangent", S.diff(ja,ell).subs(ell,0)-(2*u-(1+tau)*V-b))
    H = S.exp(2*ell*(V-u))*(1+ell*b)**2/T
    zero("mass_completion_H_tangent", S.diff(H,ell).subs(ell,0)-(2*(V-u+b)-tau*V))

    z, aseed, x, qc, d = S.symbols("z aseed x qc d")
    calB=x+z**2-2*aseed*z**3+qc*z**4-2*d*z**5+d*z**6
    polynomial=1-3*aseed*z+2*qc*z**2-5*d*z**3+3*d*z**4
    zero("Ward_weight_unimodality_polynomial", S.diff(calB,z)-2*z*polynomial)
    zero("Ward_weight_convexity_polynomial", S.diff(polynomial,z,2)-(4*qc-30*d*z+36*d*z**2))
    alpha0, ux, ut, wx, wt, thp = S.symbols("alpha0 ux ut wx wt thp")
    kernel=-S.exp(2*alpha0*(ut-ux))*wt*thp/wx
    zero("finite_Ward_metric_kernel_derivative", S.diff(kernel,alpha0)
         -2*S.exp(2*alpha0*(ut-ux))*wt*(ux-ut)*thp/wx)
    zero("constant_trace_has_zero_kernel_response", S.diff(kernel,alpha0).subs(thp,0))

    A0, B4, HH = S.symbols("A0 B4 HH")
    born1=(B4+A0**2+2*HH)/(120*S.pi)
    born2=-(2*A0**2+HH)/(60*S.pi)
    zero("general_Hessian_Born_cancellation", born1+born2-(B4-3*A0**2)/(120*S.pi))
    curv=-2*A0
    boxcurv=-2*B4+4*A0**2
    zero("common_second_heat_coefficient", born1+born2+(curv**2+2*boxcurv)/(480*S.pi))

    nu, lapnu, gradmix, lapphi, gradphi2, rhat, r4 = S.symbols(
        "nu lapnu gradmix lapphi gradphi2 rhat r4")
    varied_r4 = rhat-6*lapnu-12*gradmix+12*u*(lapphi+gradphi2)-2*nu*r4
    rhs=(rhat-varied_r4)/6+2*u*(lapphi+gradphi2)-nu*r4/3
    zero("invariant_radius_equation", lapnu+2*gradmix-rhs)

    Y=S.symbols("Y", positive=True)
    f=S.Function("f")(Y)
    D=lambda expr: Y*S.diff(expr,Y)
    zero("conormal_to_ordinary_second", Y**2*S.diff(f,Y,2)-(D(D(f))-D(f)))
    zero("conormal_to_ordinary_third", Y**3*S.diff(f,Y,3)-(D(D(D(f)))-3*D(D(f))+2*D(f)))
    j3=S.symbols("j3")
    jp=j3*Y**4/4
    w=-jp/(2*Y)
    v=-j3*Y**3/2
    zero("proper_mass_cubic_correlation", v-4*w)
    zero("proper_mass_radius_coefficient", w+j3*Y**3/8)

    t=S.symbols("t", real=True)
    eta7=1-35*t**4+84*t**5-70*t**6+20*t**7
    zero("C3_cap_derivative", S.diff(eta7,t)+140*t**3*(1-t)**3)
    check("C3_cap_endpoint_jets", eta7.subs(t,0)==1 and eta7.subs(t,1)==0 and
          all(S.diff(eta7,t,j).subs(t,e)==0 for j in (1,2,3) for e in (0,1)))
    beta=12012*t**6*(1-t)**6
    check("global_endpoint_map_cutoff_normalization", S.integrate(beta,(t,0,1))==1)
    check("global_endpoint_map_cutoff_C6", all(S.diff(beta,t,j).subs(t,e)==0
          for j in range(6) for e in (0,1)))
    N,q,C=S.symbols("N q C", positive=True)
    count=2*q*N+2*N*(N+1)
    zero("angular_count_margin", 4*N*(N+q)-count-2*N*(N+q-1))
    s,s0=S.symbols("s s0", positive=True)
    zero("middle_primitive_maximum_derivative", S.diff((s-s0)/s**2,s)-(2*s0-s)/s**3)
    zero("middle_primitive_maximum_value", ((s-s0)/s**2).subs(s,2*s0)-1/(4*s0))
    L=S.symbols("L")
    outer_primitive=-((S.log(t))**2+4*S.log(t)+5)/t
    zero("outer_absolute_integral_primitive", S.diff(outer_primitive,t)-(1+S.log(t))**2/t**2)


def finite_ledgers():
    ten=S.Integer(10)
    V=R(12,10**9); b=R(3,10**7)
    TVV=R(26,10**6); TVu=R(18,10**5); TVb=R(4,10**4)
    check("E0_TV_b", R(33,10**11)+TVV+2*V+2*TVu<TVb)
    check("E0_mass_q_primitive", (12*V**2+7*V*b)/R(14,100)<R(2,10**13))
    check("E0_signed_mass_total", (R(2,10**13)+R(65,10**14))/2<R(5,10**13))
    check("E0_lapse_integral", 3*V**2+3*V*TVV+b*TVb/(1-b)<R(13,10**11))
    check("E0_lapse_amplitude", V**2/4+b**2/(2*(1-b))<R(46,10**15))
    check("local_four_jet_coefficient", ten**7*S.factorial(4)*100**4<ten**20)
    check("local_five_jet_coefficient", ten**8*S.factorial(5)*100**5<ten**22)
    check("local_cold_oscillation_conversion", R(162,100)*R(11,10)<2)
    check("cold_normalizer_amplitude", R(1062,1000)*3*12*R(232,1000)<9)
    check("cold_normalizer_third_center", 6*ten**6*9<3*ten**9)
    check("cold_normalizer_third_extension", 3*ten**12+60*2*ten**10+1500*3*ten**9<ten**13)
    check("cold_inverse_matrix_third", 8*ten**11+6*8*7*ten**6*92
          +6*16*9200*92**2<ten**12)
    check("cold_ordered_inverse_third", 6*2*ten**6*ten**9+6*ten**11*2*ten**4
          +6*2*ten**6*(2*ten**4)**2+5*ten**13<ten**18)
    check("cold_third_glue", 60*2*400**2*ten**18<ten**26)
    check("cold_seed_angular_third", ten**26*4*R(352,100000)/(6*100**3)<R(24,10)*ten**17)
    check("cold_seed_total", R(24,10)*ten**17+ten**8<ten**20)
    check("middle_seed_row_coefficient", 4*(3*ten**4)**2*6*2000**3*2*256**2<R(23,10)*ten**25)
    check("middle_seed_row_normalization", R(23,10)*ten**25/(6*100**3)<R(39,10)*ten**18)
    check("middle_angular_moment", 768*65**5<R(892,100)*ten**11)
    check("middle_seed_complete", R(39,10)*ten**18*R(892,100)*ten**11
          *R(124,10**7)+2*ten**13+ten**12<ten**26)
    check("middle_curvature_force", 10*R(68,10**29)*ten**26==R(68,100))
    check("middle_fourth_metric", R(68,100)*2*100**2/(2*24*100**4)<R(284,10**8))
    check("middle_fifth_metric", R(68,100)*6*100**3/(2*120*100**5)<=R(17,10**7))
    check("finite_endpoint_composition_count", S.factorial(6)**2*2**6<34*ten**6)
    check("finite_endpoint_composition_smallness", 34*ten**6*R(1,10**10)<R(1,100))
    check("cold_Born_fifth_loop", S.factorial(5)*4**5*100==12288000)
    check("cold_Born_fourth_loop", S.factorial(4)*4**4*100==614400)
    check("cold_Born_value", 2*12288000*ten**12<ten**22)
    check("cold_Born_mixed", 2*614400*ten**14<ten**22)
    e3=sum(R(3)**j/S.factorial(j) for j in range(10))
    e2=sum(R(2)**j/S.factorial(j) for j in range(10))
    check("cold_exponential_second_moment", 32**2/e2<150)
    check("cold_exponential_third_moment", 48**3/e3<5600)
    check("middle_exponential_third_moment", 390**3/e3<3*ten**6)
    check("cold_angular_oscillation", (R(8,3)+32)*ten**22<ten**24)
    check("middle_two_norm_response", ten**22*R(1,10**4)+ten**24*R(1,10**6)==2*ten**18)
    check("outer_short_safe", 48*ten**21+96*ten**13+2*ten**17+4*ten**23+ten**16<ten**24)
    check("outer_mass_mode_retained", ten**24*(R(1,10**7)+R(4,10**11))<R(1001,1000)*ten**17)
    pref=R(85,10**34); P=R(5,10**32)
    Tc=R(11,10)*ten**17; Tosc=2*ten**18; Tm=3*ten**18; To=2*ten**17
    check("cold_full_value_budget", ten**23*R(1,10**6)+ten**20*R(1,10**5)+ten**12<Tc)
    check("cold_full_oscillation_budget", ten**24*R(1,10**6)+2*ten**22*R(1,10**5)+ten**14<Tosc)
    check("middle_full_value_budget", 2*ten**18+ten**16+ten**15<Tm)
    check("outer_full_value_budget", R(1001,1000)*ten**17+ten**15<To)
    check("physical_cold_seam_middle_integral", 175+(1000-100)/R(8,10)==1300)
    check("middle_primitive_coefficient", 2/(R(9,10)*4*R(1,1000)*R(1001,1000))<600)
    coldM=2*pref*Tc*81
    middleM=600*pref*Tm
    check("source_cold_mass_primitive", coldM<R(16,10**14))
    check("source_middle_mass_primitive", middleM<=R(153,10**13))
    check("source_total_mass_primitive", coldM+middleM+R(1,10**16)+R(3,10**18)<R(2,10**11))
    near=12*pref*Tosc*81
    neck=26300*pref*Tc
    middle=2600*pref*Tm
    check("source_near_bolt_lapse", near<R(166,10**13))
    check("source_long_cold_lapse", neck<R(246,10**13))
    check("source_middle_lapse", middle<=R(663,10**13))
    check("source_total_lapse", near+neck+middle+R(1,10**13)+R(1,10**17)+R(27,10**17)<R(12,10**11))
    e4=sum(R(4)**j/S.factorial(j) for j in range(10))
    check("outer_log40_less_than_four", e4>40)
    check("outer_log_integral_bound", R(4**2+4*4+5,40)<1)
    check("absolute_source_mass_coefficient", 2*Tm+2*Tm/40+4*To<7*ten**18)
    check("absolute_source_mass_image", 7*ten**18*pref<R(1,10**13))
    hor=P*Tc*81/10
    check("horizon_orbit_force", hor<R(45,10**15))
    angular=R(2,5)*R(22,7)*P*ten**15*81
    check("horizon_angular_half_trace", angular<R(51,10**16))
    check("horizon_complete_force", hor+R(1,10**14)+R(1,10**20)<R(1,10**13))
    check("full_retained_mass_primitive", R(5,10**13)+R(2,10**11)<R(21,10**12))
    check("full_retained_lapse_integral", R(13,10**11)+R(12,10**11)<R(26,10**11))
    check("full_retained_horizon_force", R(1,10**17)+R(1,10**13)<R(2,10**13))
    check("component_Einstein_zero_datum_F", 2*R(21,10**12)==R(42,10**12))
    check("component_Einstein_lapse_amplitude", R(46,10**15)+R(12,10**11)<R(13,10**11))
    check("component_Einstein_prepared_gap", R(3,2)*(R(12,10**11)+R(1,10**13))+R(1,10**17)<R(2,10**10))


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
            exact_identities()
            finite_ledgers()
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
