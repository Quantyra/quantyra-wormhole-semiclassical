"""S167 actual projected first derivative and bolt norm-embedding audit."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "projected-first-derivative-bound"
OWN = ["checks/projected_first_derivative_bound.py", "notes/"+STEM+".md"]
PREDECESSORS = [("projected-low-correction-and-mass", 294),
                ("tt-form-and-scalar-fourth-graph", 238),
                ("global-compatible-scalar-response", 344)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "first derivatives of the actual compatible seed response and projected correction",
    "weighted_scalar_gradient": "||W gradbar tau||/sqrt(gamma0)<2e-4 ell^2",
    "unweighted_full_gradient": "||gradbar(S,T)||/sqrt(gamma0)<1.3e-4 ell^2",
    "unweighted_full_H1": "H1(S,T)/sqrt(gamma0)<1.31e-4 ell^2, including actual mass",
    "method": "differentiate separated resolvent terms; retain positive TT form and exact return",
    "no_source_derivative_required": True,
    "bolt_fixture_scope": "norm non-embedding only; not the actual correction",
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
    for p in ["stationary-weighted-tt-inverse", "weighted-york-source-conversion",
              "stationary-scalar-outgoing-inverse", "matched-residual-inverse-composition"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def source_and_resolvent_budgets():
    sg, qg = R("4e-17"), R("6e-9")
    check("strong_total_q", 36000+250*8644+R("1e-6") < R("2.2e6"))
    check("strong_z_weighted_value", R("4.8e7")+4*R("8e5")*sg < R("4.801e7"))
    e, z, q, c = S.symbols("E Z Q c", positive=True)
    # Complete the weighted energy square before the upper estimate.
    gg = S.Symbol("G", positive=True)
    zero("weighted_energy_completed_square", gg**2-e*z**2-q*z-2*c*gg*z
         -((gg-c*z)**2-(e+c*c)*z*z-q*z))
    check("energy_cross_square_reserve", R("2.2e6")*R("4.801e7") < R("1.03e7")**2)
    check("twice_log_weight_gradient", 4*20 < 9**2)
    zgrad = 15*R("4.801e7")+R("1.03e7")*qg+9*R("4.801e7")*sg
    check("weighted_z_gradient", zgrad < R("7.3e8"), zgrad)
    aa, bb = S.symbols("A B", real=True)
    zero("elliptic_pair_Cauchy", 10*(aa*aa+bb*bb)-(aa+3*bb)**2-(3*aa-bb)**2)
    check("elliptic_weighted_gradient_254", 80**2*10 < 254**2)
    weighted = (R("7.3e8")+254*R("2.2e6"))/180
    check("scalar_weighted_gradient_72e6", weighted < R("7.2e6"), weighted)
    check("TT_vertex_z", 201*R("4.801e7")+R(".001")*R("7.3e8") < R("9.7e9"))
    check("local_vertex_z", R("1.3e-10")*R("4.801e7")+R("2.24")*R("7.3e8") < R("1.7e9"))
    check("coercivity_reciprocal", 1/(1-R("1e-21")) < R("1.001"))
    check("coercivity_sqrt_reciprocal", 1/(1-R("1e-21")) < R("1.001")**2)
    check("local_sqrt_guard", R("1e-21") < R("1e-10")**2)
    check("e0_three_half_guard", 180**3 > 2400**2)
    unweighted = R("7.3e8")/180+R("1.001")*(15000+R(".2")*sg)/180
    unweighted += R("1.001")*(R("9.7e9")+R("1e-10")*R("1.7e9"))*sg/2400
    check("scalar_unweighted_gradient_41e6", unweighted < R("4.1e6"), unweighted)
    check("weighted_inverse_gain_over_interpolation", R("7.2e6")*qg < R("6.3e6"))


def noncommuting_form_fixture():
    hh = S.Matrix([[5, 1, 0], [1, 8, 1], [0, 1, 12]])
    vt = S.Matrix([[2, 1], [3, -1], [1, 4]])
    vs = S.Matrix([[1, 0], [0, 2], [-1, 1]])*R("1e-12")
    js = S.diag(1, -1)
    ee, b0 = R(25), R("1e-21")
    bt, bs, a = vt/5, vs/5, 1-b0
    ww = vt*vt.T-vs*js*vs.T
    kk = hh+ww/ee
    ki = kk.inv()
    check("fixture_H0_positive", (hh-4*S.eye(3)).is_positive_definite)
    check("fixture_noncommuting", hh*ww != ww*hh)
    check("fixture_local_form_small", (b0*hh-bs*bs.T).is_positive_definite)
    check("fixture_K_lower_form", (kk-a*hh-bt*bt.T).is_positive_definite)
    check("fixture_source_energy_factorization", (hh.inv()/a**2-ki*hh*ki).is_positive_definite)
    check("fixture_TT_energy_factorization", (S.eye(2)/a-bt.T*ki*hh*ki*bt).is_positive_definite)
    check("fixture_local_energy_factorization", (b0*S.eye(2)/a**2-bs.T*ki*hh*ki*bs).is_positive_definite)
    pp = S.Matrix([7, -5, 11])
    tau = (hh*(hh-ee*S.eye(3))-ww).inv()*pp
    zz = hh*tau
    check("fixture_exact_regrouping", tau == zz/ee-ki*pp/ee-ki*ww*zz/ee**2)


def return_derivatives():
    a, r, frr, b, c, cp, eta = S.symbols("A R Frr b c cp eta")
    ap, rp = frr/2-a*a, a*r-r*r
    bp, etap = -2*a*(b+c)-2*r*(b-c), a*eta+b+c
    kg = frr/2-2*a*r+r*r
    scp = cp/2-rp*eta-r*etap
    tcp = -(bp+3*cp)/2-(ap-rp)*eta-(a-r)*etap
    zero("exact_Sc_derivative", scp-(cp/2-(2*a*r-r*r)*eta-r*(b+c)))
    zero("exact_Tc_derivative", tcp-(-3*cp/2+2*r*b-kg*eta))
    check("sqrt5_reserve", 5 < R("2.24")**2)
    check("smooth_Sc_coefficient", 2*11+5 == 27)
    check("Sc_gradient_650", R(34, 2)+27*22+3*R("2.24")*R("5.1") < 650)
    check("Tc_gradient_980", R(3, 2)*34+4*R("2.24")*R("5.1")+40*22 < 980)
    check("disk_return_pair_1200", 650**2+980**2 < 1200**2)
    check("actual_TT_auxiliary_3264", R("5.1")*64 == R("326.4"))
    check("actual_return_gradient_391680", 1200*R("326.4") == 391680)
    check("full_disk_gradient_45e6", R("4.1e6")+391680 < R("4.5e6"))


def bolt_nonembedding():
    v = S.Symbol("v", real=True)
    derivative = 2772*v**5*(1-v)**5
    step = S.integrate(derivative, (v, 0, v))
    zero("bolt_step_left_value", step.subs(v, 0))
    zero("bolt_step_right_value", step.subs(v, 1)-1)
    for j in range(1, 6):
        zero("bolt_step_left_jet_"+str(j), S.diff(step, v, j).subs(v, 0))
        zero("bolt_step_right_jet_"+str(j), S.diff(step, v, j).subs(v, 1))
    zero("bolt_step_gradient", S.diff(step, v)-derivative)
    # v(1-v)<=1/4 on [0,1], so the derivative is nonnegative and <3.
    zero("bolt_step_midpoint_identity", R(1, 4)-v*(1-v)-(v-R(1, 2))**2)
    check("bolt_step_gradient_bound", R(2772, 4**5) < 3)
    energy = S.integrate(S.expand(derivative**2), (v, 0, 1))
    check("bolt_step_exact_energy", 0 < energy < 9, energy)
    rho, rho0, length = S.symbols("rho rho0 L", positive=True)
    yy = S.log(rho0/rho)
    # du/drho and rho drho=-rho^2 dy give the exact derivative cancellation.
    zero("bolt_log_chain_derivative", S.diff((yy-length)/length, rho)+1/(length*rho))
    zero("bolt_energy_measure_cancellation", (1/(length*rho))**2*rho*rho-1/length**2)
    zero("bolt_support_L2_integral", S.integrate(rho, (rho, 0, rho0*S.exp(-length)))
         -rho0**2*S.exp(-2*length)/2)
    zero("bolt_gradient_energy_limit", S.limit(energy/length, length, S.oo))
    zero("bolt_value_energy_limit", S.limit(rho0**2*S.exp(-2*length)/2, length, S.oo))
    check("bolt_value_remains_one", step.subs(v, 1) == 1)


def physical_compositions():
    for label, coefficient, target in [("scalar_weighted_gradient", "7.2e6", "2e-4"),
                                       ("full_unweighted_gradient", "4.5e6", "1.3e-4")]:
        value = R(coefficient)*124*R("2.2e-13")
        check("actual_projected_"+label, value < R(target), value)
    check("actual_projected_full_H1", R("1.2e-6")**2+R("1.3e-4")**2 < R("1.31e-4")**2)
    check("no_S167_acceptance_promoted", not PROTOCOL["nonlinear_family_preserved"]
          and not PROTOCOL["useful_Newton_certificate"])


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 projected first derivative bound",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        source_and_resolvent_budgets()
        noncommuting_form_fixture()
        return_derivatives()
        bolt_nonembedding()
        physical_compositions()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 projected first derivative bound",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "actual_compatible_first_derivative_bound": good,
        "actual_projected_full_H1_bound": good, "bolt_H1_pointwise_nonembedding": good,
        "candidate_Hessian_inverse": False, "signed_compatibility_evaluated": False,
        "nonlinear_family_preserved": False, "useful_Newton_certificate": False,
        "S167_complete": False, "goal_complete": False,
        "unresolved": ["local higher jets and cold weights", "nonlinear family and remainder",
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
