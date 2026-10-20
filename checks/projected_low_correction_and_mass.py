"""S167: ordinary-proof algebra and quantitative low correction/mass budgets."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "projected-low-correction-and-mass"
OWN = ["checks/projected_low_correction_and_mass.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("matched-residual-inverse-composition", 288),
    ("global-compatible-scalar-response", 344),
    ("strong-effective-tt-source", 350),
    ("tt-mass-numerical-bound", 398),
    ("weighted-york-source-conversion", 293),
    ("stationary-centered-complete-residual", 92),
]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual compatible seed inverse and projected comparison low response",
    "forcing": "D(f_ell)<2.2e-13 ell^2/g",
    "scalar_weighted_value": "||W tau||/sqrt(gamma0)<2.2e-5 ell^2",
    "metric_unweighted_value": "||(S,T)||/sqrt(gamma0)<1.2e-6 ell^2, including mass",
    "mass": "|m_c|<1.2e-4 ell^2, actual sphere TT coefficient",
    "fixture_scope": "exact noncommuting matrix algebra only; not a physical state",
    "candidate_Hessian_inverse": False, "nonlinear_family_preserved": False,
    "signed_compatibility_evaluated": False, "useful_Newton_certificate": False,
    "formal_verification": False,
}


def check(name, ok, detail=""):
    assert name not in {row["name"] for row in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    expression = S.cancel(S.expand(expression))
    check(name, expression == 0, expression)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    found = set(OWN)
    for stem, _ in PREDECESSORS:
        found.add("evidence/"+stem+".json")
        found.update(inputs(latest(stem)))
    return sorted(found)


def predecessors():
    expected = {}
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(row["passed"] for row in run["checks"]) and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for p, h in sorted(expected.items()):
        check("input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)
    for p in ["stationary-centered-first-order-geometry", "stationary-completed-metric-response-norm",
              "stationary-centered-covariant-seed", "comparison-path-complete-source-response",
              "comparison-path-complete-equation-residual", "stationary-scalar-outgoing-inverse"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def elliptic_fixture():
    hh = S.Matrix([[5, 1, 0], [1, 8, 1], [0, 1, 12]])
    vt = S.Matrix([[2, 1], [3, -1], [1, 4]])
    vs = S.Matrix([[1, 0], [0, 2], [-1, 1]])*R("1e-12")
    js = S.diag(1, -1)
    ee = R(23)
    ww = vt*vt.T-vs*js*vs.T
    kk = hh+ww/ee
    pp = S.Matrix([7, -5, 11])
    tau = (hh*(hh-ee*S.eye(3))-ww).inv()*pp
    zz = hh*tau
    qq = pp+ww*tau
    check("fixture_noncommuting", hh*ww != ww*hh)
    check("fixture_H0_at_least_four", (hh-4*S.eye(3)).is_positive_definite)
    check("fixture_small_local_form", (R("1e-21")*hh-vs*vs.T/ee).is_positive_definite)
    check("fixture_K_coercivity", (kk-(1-R("1e-21"))*hh-vt*vt.T/ee).is_positive_definite)
    check("fixture_reference_solution", (hh-ee*S.eye(3))*zz == qq)
    check("fixture_outgoing_partial_fraction", tau == (zz-hh.inv()*qq)/ee)
    check("fixture_form_regrouping", kk*tau == (hh*zz-pp)/ee)
    regrouped = zz/ee-kk.inv()*pp/ee-kk.inv()*ww*zz/ee**2
    check("fixture_exact_elliptic_extraction", tau == regrouped)
    check("fixture_vertex_term_cannot_be_dropped", tau != zz/ee-kk.inv()*pp/ee)


def inverse_budgets():
    sg, qg, ig = R("4e-17"), R("6e-9"), R("2e8")
    gm = sg**2
    check("benchmark_sqrt_guard", R(2000)**-10 < sg**2)
    check("benchmark_fourth_guard", R(2000)**-10 < qg**4)
    check("benchmark_inverse_fourth_guard", 1/(R(".999")*R(2000)**-10) < ig**4)
    check("frequency_sqrt_E", 60*R(22, 7)+R(11, 3)*gm < 15**2)
    check("total_source_old_A", 3001*30000+6100000*R(".2")*sg < R("9.004e7"))
    check("total_source_old_B", 9*30000+R("1.001")*R(".2") < 270001)
    check("z_value", R("3.2e9")+4*R("4.4e7")*sg < R("3.201e9"))
    # Bound sqrt(Aq*z) by 6e8*g^(1/4), without floating square roots.
    check("z_energy_cross", R("9.004e7")*R("3.201e9") < R("6e8")**2)
    check("z_gradient", 15*R("3.201e9")+R("6e8")*qg+270001*sg < R("4.9e10"))
    check("z_TT_vertex", 201*R("3.201e9")+R(".001")*R("4.9e10") < R("6.5e11"))
    check("z_local_vertex", 1000*R("3.201e9")+3*R("4.9e10") < R("3.4e12"))
    check("coercivity_reciprocal", 1/(1-R("1e-21")) < R("1.001"))
    check("coercivity_sqrt_reciprocal", 1/(1-R("1e-21")) < R("1.001")**2)
    check("local_sqrt_reserve", R("1e-21") < R("1e-10")**2)
    check("e0_three_half_reserve", 180**3 > 2400**2)
    value = R("3.201e9")*sg/180+R("1.001")*(R(30000, 4)+R(".1")*sg)/180
    value += R("1.001")*(R("6.5e11")+R("1e-10")*R("3.4e12"))*sg/(2*2400)
    check("tau_unweighted_43", value < 43, value)
    check("tau_unweighted_gradient", 43*R("3.2e9") < R("3.8e5")**2)
    check("tau_unweighted_Hessian", R("3.2e9")+4*R("3.8e5")*qg < R("3.3e9"))
    w = S.Symbol("w", positive=True)
    zero("decay_value_weight", w**-2*w**4-w**2)
    zero("decay_derivative_weight", w**-2*w**3-w)
    mixed = (5*((7+20*R("4e-27")/180+R("5e5")*gm/180)*43
                +R(200, 180)*R("3.3e9")*sg+R(8000, 180)*R("3.8e5")*sg*qg)
             +27*R("6e-22")*R("3.8e5")*ig)/24
    check("mixed_strong_source_63", mixed < 63, mixed)
    check("effective_TT_64", 1+63 == 64)
    check("TT_value_1664", 26*64 == 1664)
    check("metric_value_42000", 43+25*1664 < 42000)
    check("mass_4224000", R("6.6e4")*64 == R("4.224e6"))
    check("local_value_vertex", 40*R(22, 7)*R("1.3e-22") < R("1.3e-10")**2)
    check("local_divergence_vertex", 5 < R("2.24")**2)
    check("TT_auxiliary_8644", 201*43+R(".001")*R("3.8e5")*qg < 8644)
    check("local_auxiliary_006", R("1.3e-10")*43+R("2.24")*R("3.8e5")*qg < R(".006"))
    check("refined_total_A", 30000+201*8644+R("1.3e-10")*R(".006") < R("1.768e6"))
    check("refined_total_B", R(".2")+R(".001")*8644+R("2.24")*R(".006") < R("8.86"))
    check("refined_weighted_tau", (81*R("1.768e6")+1100*R("8.86")*sg)/180 < R("8e5"))
    check("refined_weighted_Ttau", 348*R("1.768e6")/13+701*R("8.86") < R("4.8e7"))
    check("refined_gradient_sqrt", R("8e5")*R("4.8e7") < R("6.2e6")**2)
    check("refined_weighted_gradient", R("6.2e6")+9*R("8e5")*qg < R("6.3e6"))


def geometry_low_jets():
    t, x, d, delta, y = S.symbols("t x d delta y")
    ss = t-1
    ff = x+2*x/t+(3*x+1)/t**2+d/t**4
    f0 = ss**2*ff-delta/t
    poly = t**3-2*t**2-5*t+12
    bracket = x*poly+4-3*t+d*(8-7*t)/t**2
    zero("seed_slope_identity", 2*ss*S.diff(f0, t)-3*f0-ss**2*bracket/t**3-delta*(5*t-2)/t**2)
    zero("polynomial_at_two", poly.subs(t, 2)-2)
    zero("polynomial_derivative_at_two", S.diff(poly, t).subs(t, 2)+1)
    check("polynomial_second_positive_interval", S.diff(poly, t, 2).subs(t, 1) > 0)
    square = x*y**3+4*x*(y-(x+3)/(8*x))**2+(31*x*x-38*x-9)/(16*x)
    zero("outer_completed_square", (x*poly+4-3*t).subs(t, y+2)-square)
    rem = S.Poly(S.expand(15*(1000+y)**2-38*(1000+y)-9), y)
    check("outer_remainder_above_x", all(c > 0 for c in rem.coeffs()))
    check("inner_bracket_positive", 1998-15*R("1e-12") > 0)
    check("outer_bracket_positive", 1000-15*R("1e-12") > 0)
    # |(8-7t)/t^2| <=8/t^2+7/t <=15 for t>=1.
    check("d_bracket_reserve", 8+7 == 15)
    zero("f_decreasing_numerator", -t**5*S.diff(ff, t)-(2*x*t**3+2*(3*x+1)*t**2+4*d))
    check("f_decreasing_margin", 2*1000+2*(3000+1)-4*R("1e-12") > 0)
    zero("C_is_f_at_one", ff.subs(t, 1)-(1+6*x+d))
    check("Fext_lower_C_over_7", 1000-1-R("1e-12") > 0)
    zero("terminal_integral_derivative", -S.diff(2*ss/f0, t)-1/f0-(2*ss*S.diff(f0, t)-3*f0)/f0**2)
    zero("terminal_boundary_zero", S.limit(2*ss/f0, t, S.oo))
    zero("thermal_slope_excess", ss*S.diff(f0, t)-2*f0-ss**3*S.diff(ff, t)-delta*(3*t-1)/t**2)
    thermal = R("1.001e-60")/R(".001")**2
    check("F0_lower_C_over_8", R(1, 7)-thermal > R(1, 8))
    check("actual_slope_upper_20001", 24*thermal < R(".0001"))
    shift = R("3.4e-40")
    check("radius_V_22e9", 2*R("1.08e-9")+shift/R(".001") < R("2.2e-9"))
    check("conformal_u_3e9", R("3.3e-10")+R("7.5e-10")/2+shift+R("2.0001")*R("1.08e-9") < R("3e-9"))
    check("Eulerian_f_ratio", R("7.5e-10")+shift*R("2.0001")/R(".001") < R("7.6e-10"))
    check("radial_b_33e9", R("3.3e-10")+R("7.6e-10")+R("1.0001")*R("2.2e-9") < R("3.3e-9"))
    check("laplacian_u_56e8", 2*R("3e-9")+R("5e-8") == R("5.6e-8"))


def residual_budgets():
    p, theta, dp = S.symbols("P Theta DeltaP")
    zero("Ward_orbit_pressure_cancellation", (2*p-theta)-2*(p-theta)-theta)
    zero("Ward_metric_response_cancellation", 2*dp-2*dp)
    eta = R("8.5e-33")
    for label, const, factor, limit in [("cold", "1.1e17", 12, "1.2e-14"),
                                       ("middle", "3e18", 12, "3.1e-13"),
                                       ("outer_log", "2e17", 24, "4.1e-14")]:
        check("endpoint_orbit_"+label, factor*eta*R(const) < R(limit))
    t, ss = R(4), R(3)
    u, n, w, b, lapu = map(R, ["3e-9", "2.2e-9", "2.2e-9", "3.3e-9", "5.6e-8"])
    fac = 2*u+w+R("1.001")*n
    cp = 10*(u+n)**2+5*(u+n)*w+2*w*w
    d1, d2 = R("2.1")*t**4, R("6.3")*t**4
    orbit = R("1.01")*fac*R("5.201e-9")+2*R("1.1")*t**4*cp
    orbit += R("2.02")*d1*n*fac+R("1.002")*d1*w*w+R("1.004")*d2*n*n
    orbit += R("1.01")*w*w*6*t**3+R("2.02")*t**3*ss*b*b+R("2.02")*ss**2*(w+t*b)**2
    angular = R("1.01")*fac*R("2.6e-9")+R(".51")*(t**4+1)*cp
    angular += R("1.01")*d1*n*fac+R(".501")*d1*w*w+R(".502")*d2*n*n
    angular += R("1.01")*t*w*lapu+R("1.01")*w*w*3*t**3
    angular += R("2.02")*t**3*ss*b*b+R("1.01")*ss**2*(w+t*b)**2
    check("middle_geometric_orbit", orbit < R("2.5e-13"), S.N(orbit, 18))
    check("middle_geometric_angular", angular < R("7.1e-14"), S.N(angular, 18))
    no, ns, ell, dt, da = S.symbols("NO NS ell DeltaTrace DeltaAngular")
    zero("full_trace_source_cancellation", no-2*ell*(dt-da)+2*(ns-ell*da)-(no+2*ns-2*ell*dt))
    for label, no, ns, tr, os, ph, vo in [
        ("cold", "3e-12", "2e-12", "3.2e-16", "1.2e-14", "7.1e-12", "3.1e-12"),
        ("middle", "2.5e-13", "7.1e-14", "3.2e-17", "3.1e-13", "4e-13", "5.7e-13"),
    ]:
        check("complete_Phi_"+label, R(no)+2*R(ns)+2*R(tr) < R(ph))
        check("complete_V_"+label, R(no)+R(os) < R(vo))
    check("complete_Phi_4_40", R("5e-15")+2*R("2e-15")+2*R("3.2e-17") < R("9.1e-15"))
    check("complete_Phi_outer", R("5e-15")+2*R("2e-15")+2*R("3.2e-14") < R("7.4e-14"))
    check("complete_V_outer_log", R("5e-15")+R("4.1e-14") == R("4.6e-14"))


def forcing_and_composition():
    # Positive series partial sums certify the needed upper bounds on logs.
    for t, upper in [(4, R("1.387")), (40, R("3.689"))]:
        check("log_upper_"+str(t), sum(upper**j/S.factorial(j) for j in range(30)) > t)
    for t, ll, n, bound in [(4, R("2.387"), 2, R("3.12")),
                            (4, R("2.387"), 4, R(60)),
                            (40, R("4.689"), 2, R(".835")),
                            (40, R("4.689"), 4, R(33))]:
        moment = sum(S.binomial(n, j)*S.factorial(j)*ll**(n-j) for j in range(n+1))/t
        check(f"log_moment_{t}_{n}", moment < bound, moment)
    t = S.Symbol("t", positive=True)
    for n in [0, 2, 4]:
        moment = sum(S.binomial(n, j)*S.factorial(j)*(1+S.log(t))**(n-j) for j in range(n+1))/t
        zero("moment_derivative_"+str(n), S.diff(moment, t)+(1+S.log(t))**n/t**2)
        zero("moment_terminal_"+str(n), S.limit(moment, t, S.oo))
    cold = (R("7.1e-12")**2+R("3.1e-12")**2)*R(".00331")**2*R(".001")
    middle = (R("4e-13")**2+R("5.7e-13")**2)*R(3, 4)
    def outer(a, b, p, i0, i2, i4):
        return a*a*i0+2*a*b*i2+(b*b+p*p)*i4
    midout = outer(R("3.1e-13"), R("5e-15"), R("9.1e-15"), R(1, 4), R("3.12"), 60)
    tail = outer(R("3.1e-13"), R("4.6e-14"), R("7.4e-14"), R(1, 40), R(".835"), 33)
    check("volume_exponential", 1/(1-2*R("1.2e-8")) < R("1.001"))
    integral = R("1.001")**2*(cold+middle+midout+tail)
    check("actual_forcing_integral", integral < R("8.5e-13")**2, S.N(integral, 20))
    tt, ell, u = S.symbols("T ell u", positive=True)
    zero("fixed_density_radius_cancellation", tt**4*S.exp(2*ell*u)*(t*tt)**-4-t**-4*S.exp(2*ell*u))
    zero("physical_weighted_measure", (t**2*t**-4)**2*t**2-t**-2)
    check("physical_forcing_22e13", R("8.5e-13")*R("6.1")/24 < R("2.2e-13"))
    for label, coefficient, target in [("metric", 42000, "1.2e-6"),
                                       ("scalar_fast", 800000, "2.2e-5"),
                                       ("mass", 4224000, "1.2e-4")]:
        value = coefficient*124*R("2.2e-13")
        check("actual_projected_"+label, value < R(target), value)
    derivative = R("6.3e6")*124*R("2.2e-13")*R("2e8")
    check("weighted_gradient_only_35e4", derivative < R("3.5e4"), derivative)
    check("gradient_upper_not_small_certificate", derivative > 1)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 low projected correction and mass",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        elliptic_fixture()
        inverse_budgets()
        geometry_low_jets()
        residual_budgets()
        forcing_and_composition()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 low projected correction and mass",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "sharper_actual_inverse_values": good, "actual_residual_physical_forcing_bound": good,
        "projected_low_correction_and_mass": good,
        "candidate_Hessian_inverse": False, "signed_compatibility_evaluated": False,
        "nonlinear_family_preserved": False, "useful_Newton_certificate": False,
        "S167_complete": False, "goal_complete": False,
        "unresolved": ["local derivative norms", "nonlinear family and remainder",
                       "projected nonlinear solve", "signed compatibility", "physical requirements"],
    }
    if not good:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_bytes()) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
