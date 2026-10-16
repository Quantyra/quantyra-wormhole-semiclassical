"""Focused exact verification of the complete fixed-RN source/clock proof.

The uniform heat, Born and localization estimates are ordinary proofs in
the three frozen notes. This checks their algebra and conservative margins;
it is not a stress sampler, a metric-response estimate or a nonlinear solver.
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
RECEIPT = ROOT / "evidence/stationary-reference-clock.json"
LOG = ROOT / "evidence/stationary-reference-clock.log"
INPUTS = [
    "checks/stationary_reference_clock.py",
    "notes/stationary-reference-higher-mode-clock.md",
    "notes/stationary-reference-cold-bolt.md",
    "notes/stationary-reference-nonbolt-stress.md",
    "notes/spherical-anomaly-action.md",
    "notes/stationary-exterior-conformal-domain.md",
    "notes/stationary-weyl-radial-sources.md",
    "notes/stationary-paired-dirac-phase.md",
    "notes/stationary-massive-boundary-layer.md",
    "checks/stationary_reference_clock.failed-first.py",
]
CONTRACT = {
    "acceptance": "exact symbolic zero; strict rational conservative margins; all provenance and scope guards true",
    "numerical_tolerances": "none",
    "fixed_parameters": {"k": 2000, "eta": "1/2000", "z_i": "[1/4,7/25]", "mu0_b": 1},
    "frozen_bounds": {"K_gradient": 10**15, "K0_anchor": 10**11,
                      "cold_gluing": "2.2e14", "total_reference_lapse": "3.126e-10"},
    "ordinary_arguments": ["uniform heat and reflected bridge bounds", "actual-domain localization",
                           "all-order Born remainder", "weighted cold gluing and row estimates",
                           "common four-dimensional angular subtraction"],
    "excluded": ["corrected-metric Frechet response", "nonlinear exterior solution",
                 "finite mass norm", "full boundary Hamiltonian charge", "physical completion"],
    "failure_policy": "append every execution and raw log; retain failed input bytes",
    "hash_policy": "SHA-256 of raw file bytes; no newline or encoding normalization",
    "retained_correction": "First run used pi>3 in the free-row rational upper bound; use the predeclared pi>31/10. Proof and 9/4 target unchanged. Failed source and run retained.",
}
ROWS = []


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    print(("PASS " if row["passed"] else "FAIL ") + json.dumps(row, sort_keys=True), flush=True)


def zero(label, expression):
    expression = S.factor(S.together(expression))
    check(label, expression == 0, expression)


def exp_lower(x, n=24):
    return sum(x**j / S.factorial(j) for j in range(n+1))


def exp_upper(x, n=24):
    # Positive x, with the successive tail ratios bounded geometrically.
    assert 0 <= x < n+2
    return exp_lower(x, n) + x**(n+1)/S.factorial(n+1)/(1-x/S.Integer(n+2))


def geometry_and_sources():
    t, v, xi, G, b = S.symbols("t v xi G b", positive=True)
    Q = 1+3*xi
    am = (v+xi*v**3+Q/v)/2
    F = 1+xi*t*t-2*am/t+Q/t**2
    factor = (t-v)/t**2 * (xi*t*(t*t+v*t+v*v)+t-Q/v)
    zero("horizon_factorization", F-factor)
    remainder = xi*(v+1)*t*t+xi*v*v*t+t-Q/v
    zero("positive_comparison_polynomial", F-xi*(t-v)*(t-1)-(t-v)*remainder/t**2)
    zero("comparison_corner", remainder.subs({t:1,v:1}))
    zero("monotonicity_numerator", S.diff(F,t)*t**3/2-(xi*t**4+am*t-Q))
    thermal = (1/v**2+1/v**3+3*xi*(1+1/v+1/v**2+1/v**3))/(2*(1+6*xi))
    zero("exact_temperature_factor", S.diff(F,t).subs(t,v)/(2*(1+6*xi))-(v-1)*thermal)
    zero("temperature_extremal_limit", thermal.subs(v,1)-1)
    primitive = S.log((t-1)/t)+1/t
    zero("clock_integral_primitive", S.diff(primitive,t)-1/(t*t*(t-1)))
    # The independent I variable carries its defining differential equation.
    theta = S.Function("Theta")(t)
    I = S.Function("I")(t)
    B = F/t**2
    pressure = I/(2*F)
    rho = pressure-theta/t**2
    ward = S.diff(pressure,t)+S.diff(F,t)*(rho+pressure)/(2*F)
    ward = ward.subs(S.diff(I,t),S.diff(B,t)*theta)
    zero("orbit_stress_Ward_identity", ward+theta/t**3)
    # Integrating I by parts is equivalent to I=B Theta-J, J'=B Theta'.
    J = S.Function("J")(t)
    zero("integration_by_parts_source", S.diff(B*theta-J,t).subs(S.diff(J,t),B*S.diff(theta,t))-S.diff(B,t)*theta)
    zero("pressure_reconstruction", pressure.subs(I,B*theta-J)-(theta/(2*t*t)-J/(2*F)))
    zero("density_reconstruction", rho.subs(I,B*theta-J)-(-theta/(2*t*t)-J/(2*F)))
    zero("lapse_reconstruction", (G*(rho+pressure)/(b*b*t*F)).subs(I,B*theta-J)+G*J/(b*b*t*F*F))
    c = S.symbols("constant_trace")
    zero("constant_trace_pressure_shift", (pressure.subs(I,I+B*c)-pressure)-c/(2*t*t))
    zero("constant_trace_density_shift", (rho.subs({I:I+B*c,theta:theta+c})-rho)+c/(2*t*t))
    u = S.symbols("u")
    coeff = S.symbols("c0:5")
    polynomial = sum(coeff[j]*u**j for j in range(5))
    R = -S.diff(polynomial,u,2)
    boxR = S.diff(polynomial*S.diff(R,u),u)
    zero("quartic_local_coefficient_derivative", S.diff(R*R+2*boxR,u)+4*S.diff(polynomial,u)*S.diff(polynomial,u,4))
    print("geometry, conservation and source identities complete", flush=True)


def angular_and_born():
    x, q = S.symbols("x q", positive=True)
    P = 2*(2*x+q)*x*(x+q)
    em = -(q+1)**2+P.subs(x,1)/2
    for j in range(1,4):
        em -= S.bernoulli(2*j)/S.factorial(2*j)*S.diff(P,x,2*j-1).subs(x,1)
    zero("angular_zeta_minus_one_EM", em-(S.Rational(1,30)-q*q/6))
    zero("EM_polynomial_sixth_derivative", S.diff(P,x,6))
    # Independent Leibniz envelopes for the differentiated EM remainder.
    def PL_bound(n):
        total = S.Integer(0)
        for a in range(min(3,n)+1):
            pb = 12*S.factorial(3)/S.factorial(3-a)
            lb = 32 if n-a == 0 else 2*S.factorial(n-a-1)
            total += S.binomial(n,a)*pb*lb
        return total
    for n,target in [(1,1176),(3,3000),(5,7056),(6,31968)]:
        zero("EM_Leibniz_order_"+str(n), PL_bound(n)-target)
    em_bound = 130+192+S.Rational(1176,12)+S.Rational(3000,720)+S.Rational(7056,30240)+S.Rational(8*15984,720)
    check("angular_zeta_derivative_bound", em_bound<604 and em_bound<1000, em_bound)
    n = S.symbols("n",positive=True,integer=True)
    zero("angular_inverse_square_telescope",
         2*(q+2*n)/(n*n*(n+q)**2)-2/q*(1/n**2-1/(n+q)**2))
    zero("angular_finite_part_partial_fraction",
         2*(q+2*n)/(n*(n+q))-2*(1/n+1/(n+q)))
    # Check the common heat subtraction's finite restoration algebra.
    z, Zm, Zmp, Z0, Fang, b, k0, k1, k2 = S.symbols("z Zm Zmp Z0 Fang b k0 k1 k2")
    local_finite = k0*Zmp/b**2+k1*Z0+b*b*(Fang+2)*k2
    anomaly = -k0*Zm/b**2+k1*Z0+2*b*b*k2
    zero("common_angular_finite_restoration", local_finite-anomaly-(k0*(Zmp+Zm)/b**2+b*b*Fang*k2))
    # Gamma evaluations of the radial K_nu squared moments.
    def bessel_square_moment(mu,nu):
        return 2**(mu-3)*S.gamma(mu/2)**2*S.gamma(mu/2+nu)*S.gamma(mu/2-nu)/S.gamma(mu)
    zero("Born_first_quadratic_moment", bessel_square_moment(S.Integer(4),0)-bessel_square_moment(S.Integer(4),1)+S.Rational(1,3))
    zero("Born_first_quartic_moment", bessel_square_moment(S.Integer(6),0)-bessel_square_moment(S.Integer(6),1)+S.Rational(8,15))
    zz, T, H2, B4 = S.symbols("z T H2 B4")
    terms = [(16*zz**2*(T*T+2*H2)*(3-zz),7),
             (-64*zz*T*T*(1-zz),6), (8*(1-3*zz)*T*T-16*zz*H2,5)]
    beta_integral = 0
    for numerator,power in terms:
        for (degree,),value in S.Poly(numerator,zz).terms():
            beta_integral += value*S.factorial(degree)*S.factorial(power-degree-2)/S.factorial(power-1)
    zero("Born_second_angular_beta_integral", beta_integral/(4*S.pi)+(2*T*T+H2)/(15*S.pi))
    born = (B4+T*T+2*H2)/(120*S.pi)-(2*T*T+H2)/(60*S.pi)
    R, boxR = -2*T,-2*B4+4*T*T
    zero("Born_common_second_heat_coefficient", born+(R*R+2*boxR)/(480*S.pi))
    print("angular finite parts and Born coefficients complete", flush=True)


def rational_margins():
    r = S.Rational
    k, q, eta = S.Integer(2000), S.Integer(2000)**4, r(1,2000)
    ratio = r(352,100000)
    C2upper = 1+12*r(22,7)*k**6
    check("benchmark_C_over_q", C2upper<ratio**2*q*q)
    check("benchmark_sharper_nonbolt_ratio", C2upper<(7*q/k)**2)
    check("benchmark_M1_squared_at_least_324", C2upper*324**2<(q+1)**2)
    check("benchmark_xi_at_least_k_cubed", k**3+3*k**6<r(31,10)*k**6)
    check("benchmark_xi_less_than_2k_cubed", 2*k**3+12*k**6>r(22,7)*k**6)
    check("horizon_log_bound_from_lower_h", exp_lower(S.Integer(80),160)>8*k**9/eta)
    check("horizon_upper_bound_from_temperature", 4*r(7,25)*eta/k**9<r(1,10**30))
    check("quintic_extension_Jacobian", r(58,1000)*(1+30*r(7,32))<r(1,2))
    check("normalizing_map_bootstrap", exp_upper(r(48,1000))-1<r(5,100))
    radius = r(5,32)
    check("H2_coordinate_distance_conversion", r(105,100)/(1-(r(105,100)*radius/2)**2)<r(106,100))
    check("reserved_weight_conversion", exp_upper(2*r(6,100)*radius)<r(102,100))
    delta = r(26,1000)
    check("flat_weighted_Neumann_ratio", 7*delta<r(1821,10000))
    check("flat_weighted_Born_row", 56*delta/(1-7*delta)<r(18,10))
    free_row_sq_upper = (8/r(31,10))/(r(5,8)*exp_lower(r(5,8)))
    check("flat_weighted_free_row", free_row_sq_upper<r(9,4))
    check("flat_weighted_L1_margin", (1+r(11,7))*18**2/r(23,2)**2<7)
    check("reserved_weight_row_budget", 4*33*r(102,100)<400)
    la = r(18,16)
    ordinary_sum = (1+2/la+2/la**2)/exp_lower(la)
    weighted_sum = (4/la**2)/exp_lower(S.Integer(2))+12/la**4
    check("angular_exponential_sum", 2*q*ordinary_sum+4*weighted_sum<4*q)
    cold_glue = r(12,10)*10**9*400**2*80*4*ratio
    check("cold_gluing_budget", cold_glue<r(22,10)*10**14,cold_glue)
    c2, c6, cp = r(1664,1000), S.Integer(10011), S.Integer(200)
    remainder = 2949120000*c2**3/(1-5*delta)+612*c6+589824000*cp
    check("local_Born_remainder_budget", remainder<r(14,10)*10**11,remainder)
    check("cold_local_radial_remainder", 7*10**13*ratio**3<r(31,10)*10**6)
    cold_total = r(22,10)*10**14+r(23,10)*10**8+r(31,10)*10**6+100
    check("cold_patch_total", cold_total<3*10**14,cold_total)
    nonbolt = r(196,100)*10**14+r(24,10)*10**13+r(14,10)*10**12+106
    check("nonbolt_patch_total", nonbolt<3*10**14,nonbolt)
    anchor_short = 8*10**15*ratio**3
    anchor_long = 1400*10**10*ratio
    check("anchor_long_even_cap_ratio", anchor_long<r(493,100)*10**10,anchor_long)
    check("trace_value_anchor_budget", anchor_short+anchor_long+171<10**11)
    check("anchor_crossing_value", 24*10**11*ratio**2*101/exp_lower(S.Integer(100),220)<1)
    K, K0 = S.Integer(10**15), S.Integer(10**11)
    check("individual_stress_budget", K+K0/2<r(11,10)*10**15)
    massive = 80*K/k**8
    zero("massive_lapse_decimal", massive-r(3125,10**13))
    lll, weyl = 140/k**9, 1/k**10
    check("LLL_reference_lapse", lll<r(2735,10**31))
    check("Weyl_reference_lapse", weyl<r(1,10**33))
    check("complete_reference_clock", massive+lll+weyl<r(3126,10**13))
    check("component_comparison_1e_minus8", massive+lll+weyl<r(1,10**8))
    print("frozen rational margins complete", flush=True)


def hashes(data):
    return {name:hashlib.sha256(value).hexdigest() for name,value in data.items()}


def main():
    before = {name:(ROOT/name).read_bytes() for name in INPUTS}
    existing = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs":[]}
    run = {"started_utc":datetime.now(timezone.utc).isoformat(),"contract":CONTRACT,
           "command":sys.orig_argv,"dont_write_bytecode":sys.dont_write_bytecode,"cwd":str(Path.cwd()),
           "python":platform.python_version(),"sympy":S.__version__,"execution_inputs_sha256":hashes(before)}
    with LOG.open("a",encoding="utf-8",newline="\n") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print("RUN "+run["started_utc"],flush=True)
        print("METADATA "+json.dumps(run,sort_keys=True),flush=True)
        try:
            check("actual_bytecode_disabled",sys.dont_write_bytecode and "-B" in sys.orig_argv)
            check("retained_first_failed_source",
                  bool(existing["runs"]) and
                  hashlib.sha256(before["checks/stationary_reference_clock.failed-first.py"]).hexdigest()
                  == existing["runs"][0]["execution_inputs_sha256"]["checks/stationary_reference_clock.py"] and
                  existing["runs"][0]["status"] == "failed" and
                  [r["label"] for r in existing["runs"][0]["failed_checks"]] == ["flat_weighted_free_row"])
            geometry_and_sources()
            angular_and_born()
            rational_margins()
            check("fixed_benchmark_scope",CONTRACT["fixed_parameters"]["k"]==2000)
            check("no_nonlinear_completion","nonlinear exterior solution" in CONTRACT["excluded"])
            check("no_metric_response_claim","corrected-metric Frechet response" in CONTRACT["excluded"])
            check("no_charge_claim","full boundary Hamiltonian charge" in CONTRACT["excluded"])
            check("new_notes_UTF8_LF",all(b"\r" not in before[name] and before[name].decode("utf-8")
                  for name in INPUTS[1:4]))
        except Exception:
            traceback.print_exc()
            check("uncaught_exception",False)
        run["final_inputs_sha256"] = hashes({name:(ROOT/name).read_bytes() for name in INPUTS})
        check("frozen_inputs_unchanged",run["execution_inputs_sha256"]==run["final_inputs_sha256"])
        run["checks"] = ROWS
        run["failed_checks"] = [row for row in ROWS if not row["passed"]]
        exit_code = int(bool(run["failed_checks"]))
        run["status"] = "failed" if exit_code else "passed"
        run["terminal_exit_status"] = exit_code
        run["finished_utc"] = datetime.now(timezone.utc).isoformat()
        print("FINAL_INPUTS "+json.dumps(run["final_inputs_sha256"],sort_keys=True),flush=True)
        if exit_code:
            for name,data in before.items():
                print(json.dumps({"path":name,"sha256":hashlib.sha256(data).hexdigest(),"base64":base64.b64encode(data).decode("ascii")},sort_keys=True))
        print(json.dumps({"status":run["status"],"checks":len(ROWS),"failures":len(run["failed_checks"]),"terminal_exit_status":exit_code},sort_keys=True),flush=True)
    existing["runs"].append(run)
    RECEIPT.write_text(json.dumps(existing,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":run["status"],"checks":len(ROWS),"failures":len(run["failed_checks"]),"receipt":str(RECEIPT)},sort_keys=True))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
