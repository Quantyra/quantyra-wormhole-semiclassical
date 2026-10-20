"""S167: regional pointwise values of the actual projected correction.

Exact radial norm arithmetic; fixtures are not physical solutions.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "projected-regional-pointwise-bound"
OWN = ["checks/projected_regional_pointwise_bound.py", "notes/"+STEM+".md"]
PREDECESSORS = [("fractional-interior-dirac-response", 318),
                ("projected-first-derivative-bound", 246),
                ("projected-low-correction-and-mass", 294)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual seed-projected full correction, including mass and fixed return",
    "pointwise_pair": "5.5e-6 ell^2/[t sqrt(t-1)] for t>=1+100h",
    "middle_value": "1.75e-4 ell^2 for t>=1.001",
    "forcing_projection": "only the established weighted physical L2 norm-one Pi0",
    "method": "separate scalar and TT physical value/gradient norms",
    "fixture_scope": "two-component power-law norm test; not a physical solution",
    "uniform_bolt_value": False, "pointwise_first_jets": False,
    "nonlinear_family_preserved": False, "candidate_Hessian_inverse": False,
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


def zero(name, expr):
    expr = S.cancel(expr)
    check(name, expr == 0, expr)


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
    for p in ["weighted-york-source-conversion", "stationary-weighted-tt-inverse",
              "projected-low-correction-and-mass", "projected-first-derivative-bound",
              "fractional-interior-dirac-response", "tt-form-and-scalar-fourth-graph"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def radial_and_seed_identities():
    t, x, d, delta = S.symbols("t x d delta", positive=True)
    f = x+2*x/t+(3*x+1)/t**2+d/t**4
    expanded = 1+x*t**2-(2+4*x+delta)/t+(1+3*x+d)/t**2-2*d/t**3+d/t**4
    zero("seed_factorization", (t-1)**2*f-delta/t-expanded)
    # Worst sign d=-1: its term is covered by the positive 1/t^2.
    zero("seed_lower_positive_decomposition", (f.subs(d, -1)-x)/x
         -(2/t+3/t**2+(t**2-1)/(x*t**4)))
    check("seed_endpoint_upper", 1+2+3+R(2, 1000) == R("6.002"))
    check("thermal_delta_bound", R("1.001")*R("6.002") < 7)
    check("thermal_factor_guard", R("1e-30") < R(".001"))
    check("regional_seed_lower_bound", 1-R(7, 100**2) == R(".9993") > R(".99"))
    ff = S.Function("F")(t)
    weight = t**2*S.sqrt(ff/x)
    zero("evaluation_weight_log_derivative", S.diff(weight, t)/weight-2/t-S.diff(ff, t)/(2*ff))
    value, derivative = S.symbols("v dv", real=True)
    zero("two_physical_weights", (t*value)*(t*S.sqrt(ff/x)*derivative)/weight-value*derivative)
    b, lam = S.symbols("B lambda", positive=True)
    zero("physical_gradient_conversion", ff/lam*(derivative/b)**2
         -(ff/x*derivative**2).subs(x, lam*b**2))
    aa, bb, cc, dd = S.symbols("a b c d", real=True)
    zero("vector_derivative_Cauchy", (aa*aa+bb*bb)*(cc*cc+dd*dd)
         -(aa*cc+bb*dd)**2-(aa*dd-bb*cc)**2)
    zero("regional_denominator_power", (t*S.sqrt(t-1))**2-t**2*(t-1))
    check("threshold_denominator_increasing", S.factor(S.diff(t**2*(t-1), t).subs(t, 1+S.Symbol("s", positive=True))).is_positive)


def actual_split_and_thresholds():
    d0 = 124*R("2.2e-13")
    check("physical_source_composition", d0 == R("2.728e-11"))
    check("TT_value_constant", 25*1664 == 41600)
    check("scalar_value_normalized", 43*d0 == R("1.17304e-9"))
    check("scalar_gradient_normalized", R("4.1e6")*d0 == R(".000111848"))
    check("TT_value_normalized", 41600*d0 == R("1.134848e-6"))
    check("TT_gradient_normalized", 391680*d0 == R("1.06850304e-5"))
    # Fourth powers avoid floating-point square roots, including the .99^(1/4).
    a = d0**2*2*43*R("4.1e6")
    b = d0**2*2*41600*391680
    check("scalar_pointwise_numerator", a*a < R(".99")*R("5.2e-7")**4, a)
    check("TT_pointwise_numerator", b*b < R(".99")*R("4.95e-6")**4, b)
    check("full_pair_pointwise_numerator", R("5.2e-7")+R("4.95e-6") < R("5.5e-6"))
    middle_t = R("1.001")
    check("middle_value_175e4", R("5.5e-6")**2 < R("1.75e-4")**2*middle_t**2*(middle_t-1))
    check("point001_value_threshold", R("5.5e-6")**2 < R(".001")**2*R("3.1e-5"))
    check("threshold_above_thermal_region", 100*R("1e-30") < R("3.1e-5") < R(".001"))
    check("outer_pointwise_numerator", R("5.5e-6")**2 < R("5.6e-6")**2*R(39, 40))
    t = S.Symbol("t", positive=True)
    zero("outer_t_minus_one_identity", t-1-R(39, 40)*t-(t-40)/40)
    combined_square = 2*R("1.2e-6")*R("1.3e-4")
    check("combined_numerator_exceeds_17e5", combined_square**2 > R(".99")*R("1.7e-5")**4)
    check("strictly_sharper_actual_split", R("5.5e-6") < R("1.7e-5"))
    ell = S.Symbol("ell", nonnegative=True)
    zero("homogeneous_ell_squared", S.sqrt(ell**4)-ell**2)
    check("no_cold_or_jet_promotion", not PROTOCOL["uniform_bolt_value"]
          and not PROTOCOL["pointwise_first_jets"] and not PROTOCOL["nonlinear_family_preserved"])


def independent_power_fixture():
    # F/x=t^2 on t>=1. This is only a norm fixture, not the retained seed.
    t, cutoff = S.symbols("t T", positive=True)
    v = S.Matrix([t**-3, 2*t**-4])
    vp = v.diff(t)
    aa = S.integrate(t**2*v.dot(v), (t, cutoff, S.oo))
    bb = S.integrate(t**4*vp.dot(vp), (t, cutoff, S.oo))
    zero("fixture_value_integral", aa-1/(3*cutoff**3)-4/(5*cutoff**5))
    zero("fixture_gradient_integral", bb-3/cutoff**3-64/(5*cutoff**5))
    direct = v.dot(v).subs(t, cutoff)
    signed_integral = -2*S.integrate(v.dot(vp), (t, cutoff, S.oo))
    zero("fixture_fundamental_theorem", direct-signed_integral)
    margin = S.Poly(S.cancel(4*aa*bb-cutoff**6*direct**2)*cutoff**10, cutoff)
    check("fixture_pointwise_Cauchy_margin", all(c > 0 for c in margin.all_coeffs() if c != 0), margin.as_expr())
    check("fixture_noncollinear_components", S.simplify(v[0]*vp[1]-v[1]*vp[0]) != 0)
    for endpoint in [R(1), R("1.001"), R(2), R(40), R(100)]:
        check("fixture_endpoint_"+str(endpoint), direct.subs(cutoff, endpoint)**2
              < (4*aa*bb/cutoff**6).subs(cutoff, endpoint))
    check("fixture_vanishing_outer_limit", all(S.limit(entry, t, S.oo) == 0 for entry in v))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 regional pointwise correction",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        radial_and_seed_identities()
        actual_split_and_thresholds()
        independent_power_fixture()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 regional pointwise correction",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "actual_regional_pointwise_bound": good, "separate_inverse_gain": good,
        "middle_value_below_point001": good,
        "uniform_bolt_value": False, "pointwise_first_jets": False,
        "candidate_Hessian_inverse": False, "signed_compatibility_evaluated": False,
        "nonlinear_family_preserved": False, "useful_Newton_certificate": False,
        "S167_complete": False, "goal_complete": False,
        "unresolved": ["cold and higher spatial correction norms", "nonlinear family",
                       "projected solve", "signed compatibility", "physical requirements"],
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
