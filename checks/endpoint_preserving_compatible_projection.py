"""S167: compatible projection preserving source endpoint neighborhoods.

Exact norm and range algebra; no physical adjoint or quantile computation.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "endpoint-preserving-compatible-projection"
OWN = ["checks/endpoint_preserving_compatible_projection.py", "notes/"+STEM+".md"]
PREDECESSORS = [("projected-regional-pointwise-bound", 243),
                ("matched-residual-inverse-composition", 288),
                ("stationary-scalar-phase-moment", 80)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "fixed actual seed-compatible projection in physical weighted L2",
    "quantile_levels": ["1/200", "199/200"],
    "forcing_norm_bound": "1.006", "difference_from_orthogonal": ".101",
    "source_endpoint_germs_preserved": True,
    "source_cold_region_entirely_unchanged": False,
    "quantile_radii_evaluated": False, "higher_projection_norm_quantified": False,
    "fixture_scope": "exact finite projection and probability examples only",
    "cold_metric_bound": False, "nonlinear_family_preserved": False,
    "candidate_Hessian_inverse": False, "signed_compatibility_evaluated": False,
    "useful_Newton_certificate": False, "formal_verification": False,
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
    for p in ["stationary-scalar-phase-moment", "stationary-adjoint-phase",
              "stationary-shape-logarithmic-response", "matched-residual-inverse-composition",
              "projected-first-derivative-bound", "projected-regional-pointwise-bound"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def zero_matrix(name, matrix):
    check(name, all(S.cancel(entry) == 0 for entry in matrix), matrix)


def abstract_projection_algebra():
    a, b, alpha, lam = S.symbols("a b alpha lambda", real=True)
    n = S.Matrix([1, 0, 0])
    w = S.Matrix([0, a, b])
    p = n+w
    pp = S.eye(3)-p*n.T
    p0 = S.eye(3)-n*n.T
    zero_matrix("abstract_idempotence", pp*pp-pp)
    zero_matrix("abstract_range", n.T*pp)
    zero_matrix("abstract_compatible_identity", pp*S.Matrix([0, a, b])-S.Matrix([0, a, b]))
    zero("abstract_projection_norm_spectrum", (lam*S.eye(3)-pp*pp.T).det()
         -lam*(lam-1)*(lam-(1+a*a+b*b)))
    zero_matrix("abstract_difference_rank_one", pp-p0+w*n.T)
    zero("abstract_difference_norm", S.trace((pp-p0)*(pp-p0).T)-a*a-b*b)
    vv = alpha*p
    zero("cutoff_beta_formula", vv.dot(vv)-alpha**2*(1+w.dot(w)))
    zero_matrix("adjoint_sign_invariance", S.eye(3)-(-p)*(-n).T-pp)
    check("physical_forcing_norm_budget", R(100, 99) < R("1.006")**2)
    check("projection_difference_budget", R(1, 99) < R(".101")**2)
    check("core_quantile_mass", R(199, 200)-R(1, 200) == R(".99"))
    check("strict_margin_over_one", R("1.006") > 1)


def continuous_quantile_fixture():
    # Actual quantile radii depend on the physical adjoint and are not sampled.
    t = S.Symbol("t", nonnegative=True)
    density = 1/(1+t)**2
    cdf = t/(1+t)
    zero("fixture_cdf_derivative", S.diff(cdf, t)-density)
    check("fixture_probability_total", S.integrate(density, (t, 0, S.oo)) == 1)
    aa, bb = R(1, 199), R(199)
    check("fixture_inner_quantile", cdf.subs(t, aa) == R(1, 200))
    check("fixture_outer_quantile", cdf.subs(t, bb) == R(199, 200))
    check("fixture_core_mass", S.integrate(density, (t, aa, bb)) == R(".99"))
    check("fixture_inner_tail", S.integrate(density, (t, 0, aa)) == R(1, 200))
    check("fixture_outer_tail", S.integrate(density, (t, bb, S.oo)) == R(1, 200))
    left, right = aa/2, 2*bb
    check("fixture_compact_cutoff_order", 0 < left < aa < bb < right)
    check("fixture_inner_endpoint_neighborhood", left > 0)
    check("fixture_transition_intervals_finite", aa-left > 0 and right-bb > 0)
    # Flat cutoff endpoint proof: each derivative is exp(-1/u) times a
    # polynomial in 1/u. The finite checks support, but do not replace, it.
    u = S.Symbol("u", positive=True)
    flat = S.exp(-1/u)
    for j in range(5):
        check("fixture_flat_cutoff_endpoint_"+str(j), S.limit(S.diff(flat, u, j), u, 0, dir="+") == 0)


def finite_projection_and_inverse_fixture():
    n = S.Matrix([R(399, 401), R(40, 401), 0])
    mask = S.diag(1, 0, 0)
    v = mask*n
    alpha, b2 = n.dot(v), v.dot(v)
    p = v/alpha
    pc = S.eye(3)-p*n.T
    p0 = S.eye(3)-n*n.T
    check("fixture_unit_adjoint", n.dot(n) == 1)
    check("fixture_cutoff_mass_budget", R(".99") < alpha <= 1)
    check("fixture_cutoff_Cauchy_bounds", alpha**2 <= b2 <= alpha)
    check("fixture_direction_normalized", n.dot(p) == 1)
    zero_matrix("fixture_projection_idempotent", pc*pc-pc)
    zero_matrix("fixture_projection_range", n.T*pc)
    check("fixture_cut_projection_not_orthogonal", pc != pc.T)
    check("fixture_orthogonal_projection_symmetric", p0 == p0.T)
    lam = S.Symbol("lambda")
    zero("fixture_exact_norm_spectrum", (lam*S.eye(3)-pc*pc.T).det()
         -lam*(lam-1)*(lam-R(401, 399)**2))
    check("fixture_projection_norm", R(401, 399) < R("1.006"))
    zero("fixture_difference_norm_squared", S.trace((pc-p0)*(pc-p0).T)-R(40, 399)**2)
    check("fixture_difference_norm", R(40, 399) < R(".101"))
    f = S.Matrix(S.symbols("f1 f2 f3"))
    endpoint = S.diag(0, 1, 1)
    zero_matrix("fixture_endpoint_preservation", endpoint*(pc*f-f))
    core_source = S.Matrix([1, 0, 0])
    check("fixture_orthogonal_endpoint_injection", endpoint*p0*core_source != S.zeros(3, 1))
    zero_matrix("fixture_cutoff_endpoint_no_injection", endpoint*pc*core_source)
    # A has the desired range; it is only a finite algebra fixture.
    aa = S.Matrix([[-R(40, 401), 0], [R(399, 401), 0], [0, 1]])
    gg = aa.T
    zero_matrix("fixture_inverse_on_domain", gg*aa-S.eye(2))
    zero_matrix("fixture_range_inverse", aa*gg-p0)
    hc, h0 = -gg*pc*f, -gg*p0*f
    zero_matrix("fixture_actual_linear_residual", f+aa*hc-p*(n.dot(f)))
    zero_matrix("fixture_exact_correction_difference", hc-h0-gg*(p-n)*n.dot(f))
    zero_matrix("fixture_residual_projection_zero", pc*(f+aa*hc))
    zero("fixture_signed_coefficient_retained", n.dot(f+aa*hc)-n.dot(f))
    zero_matrix("fixture_full_equation_decomposition", pc*f+p*n.dot(f)-f)


def transferred_budgets():
    for label, old, new in [
        ("full_value", "1.2e-6", "1.21e-6"),
        ("full_H1", "1.31e-4", "1.32e-4"),
        ("regional_pointwise", "5.5e-6", "5.6e-6"),
        ("middle_pointwise", "1.75e-4", "1.8e-4"),
    ]:
        check("cut_projection_"+label, R("1.006")*R(old) < R(new))
    for label, old, new in [
        ("regional", "5.5e-6", "5.6e-7"),
        ("middle", "1.75e-4", "1.8e-5"),
    ]:
        check("correction_difference_"+label, R(".101")*R(old) < R(new))
    check("middle_value_still_below_point001", R("1.8e-4") < R(".001"))
    check("no_cold_or_nonlinear_promotion", not PROTOCOL["cold_metric_bound"]
          and not PROTOCOL["nonlinear_family_preserved"]
          and not PROTOCOL["higher_projection_norm_quantified"])


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 source endpoint compatible projection",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        abstract_projection_algebra()
        continuous_quantile_fixture()
        finite_projection_and_inverse_fixture()
        transferred_budgets()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 source endpoint compatible projection",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "actual_endpoint_preserving_projection": good,
        "physical_H_norm_below_1006": good, "correction_bounds_transferred": good,
        "source_endpoint_germs_preserved": good,
        "quantile_radii_evaluated": False, "higher_projection_norm_quantified": False,
        "cold_metric_bound": False, "candidate_Hessian_inverse": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["cold and higher spatial inverse/source norms", "fixed cutoff direction norms",
                       "nonlinear family", "projected solve", "signed compatibility", "physical requirements"],
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
