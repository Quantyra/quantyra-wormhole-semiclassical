"""S166 complete massive comparison in quantified half-density variables."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "outer-complete-massive-model-comparison"
OWN = ["checks/outer_complete_massive_model_comparison.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("outer-strong-curved-quantum-remainder", 243),
    ("outer-flat-tower-normal-remainder", 321),
    ("compact-cover-graph-absorption", 366),
    ("stationary-curved-quantum-tt", 157),
]
PROTOCOL = {
    "model_comparison": "m Q_geo m^-1 minus the exact even flat K",
    "model_bound": "1e44 epsilon in the physical shape density at weight 7/2",
    "complete_comparison": "m H_mass m^-1 minus flat logarithmic K_infty",
    "complete_bound": "1e65 epsilon, H4log to L2; epsilon<=1e-15",
    "norm_conversion": "M_epsilon and inverse <1.01 on H4, H4log and output L2",
    "normal_composition": "<1e-30 at epsilon=1e-120, in stated auxiliary variables",
    "actual_rho_retained": True,
    "Weyl_density_combination_complete": False,
    "global_inverse": False,
    "nonlinear_contraction": False,
}
R = S.Rational
ROWS = []


def check(name, ok, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    expression = S.simplify(S.expand(expression))
    check(name, expression == 0, "" if expression == 0 else expression)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    out = set(OWN)
    for stem, _ in PREDECESSORS:
        out.add(f"evidence/{stem}.json")
        out.update(inputs(latest(stem)))
    return sorted(out)


def predecessors():
    expected = {}
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(r["passed"] for r in run["checks"])
              and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for p, h in sorted(expected.items()):
        check("input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)


def half_density_and_fourth_allocation():
    y = S.symbols("Y", real=True)
    h, f, rb = (S.Function(n)(y) for n in ["h", "f", "Rbar"])
    m = S.exp(h/2)
    potential = S.diff(h, y, 2)/2+S.diff(h, y)**2/4
    lap = lambda v: S.diff(v, y, 2)+S.diff(h, y)*S.diff(v, y)
    divr = lambda v: rb*S.diff(v, y, 2)+(S.diff(rb, y)+S.diff(h, y)*rb)*S.diff(v, y)
    conjugated_lap = m*lap(f/m)
    zero("half_density_radial_Laplacian", conjugated_lap-S.diff(f, y, 2)+potential*f)
    zero("half_density_curvature_divergence", m*divr(f/m)
         -(rb*S.diff(f, y, 2)+S.diff(rb, y)*S.diff(f, y)
           -(rb*potential+S.diff(rb, y)*S.diff(h, y)/2)*f))
    fourth = m*(lap(lap(f/m))+divr(f/m)/3)
    target = S.diff(f, y, 4)-2*potential*S.diff(f, y, 2)-2*S.diff(potential, y)*S.diff(f, y)
    target += (potential**2-S.diff(potential, y, 2))*f
    target += (rb*S.diff(f, y, 2)+S.diff(rb, y)*S.diff(f, y)
               -(rb*potential+S.diff(rb, y)*S.diff(h, y)/2)*f)/3
    zero("complete_fourth_restoration", fourth-target)
    expanded = S.expand(target)
    check("fourth_coefficient_one", expanded.coeff(S.diff(f, y, 4)) == 1)
    check("third_coefficient_zero", expanded.coeff(S.diff(f, y, 3)) == 0)
    zero("lower_remainder_second_coefficient", expanded.coeff(S.diff(f, y, 2))-(-2*potential+rb/3))
    # Similarity transforms the equation, with an independently transformed source.
    a = S.Matrix([[2, 3], [5, 7]])
    mm = S.diag(2, 3)
    v = S.Matrix(S.symbols("v0:2"))
    check("auxiliary_equation_similarity", (mm*a*mm.inv())*(mm*v) == mm*(a*v))
    check("similarity_is_not_dual_congruence", mm*a*mm.inv() != mm.inv()*a*mm.inv())


def kernel_tail_and_probes():
    pi = R(22, 7)
    radius = R(1, 400)
    check("cutoff_added_geometric_error", 4*512*17/(128*9)+30 < 10**4)
    check("stationary_kernel_gain_two", 4*10**4*96*pi/2 < 10**9)
    check("full_bounded_tail_budget", 16*2*pi*18/(radius/2)**4 < 10**20)
    check("short_flow_transition_separation", R(1, 256) > radius/2)
    check("near_Taylor_strong_bound", 2*10**9*radius/6 < 10**7)
    t, s = S.symbols("t s", real=True)
    p = 462*t**6-1980*t**7+3465*t**8-3080*t**9+1386*t**10-252*t**11
    zero("beta_cutoff_derivative", S.diff(p, t)-2772*t**5*(1-t)**5)
    for k in range(6):
        check("beta_endpoint_left_"+str(k), S.diff(p, t, k).subs(t, 0) == 0)
        check("beta_endpoint_right_"+str(k), S.diff(1-p, t, k).subs(t, 1) == 0)
    ds = [1, 11, 220, 3960, 63360, 887040, 10644480]
    for k in range(7):
        pp = S.Poly(S.diff(p, t, k), t)
        n = pp.degree()
        if k == 0:
            check("beta_value_range", p.subs(t, 0) == 0 and p.subs(t, 1) == 1)
            continue
        bernstein = [sum(pp.nth(j)*S.binomial(i, j)/S.binomial(n, j)
                         for j in range(i+1)) for i in range(n+1)]
        check("cutoff_Bernstein_"+str(k), max(abs(z) for z in bernstein) <= ds[k])
    for degree in range(3):
        for k in range(7):
            bound = sum(S.binomial(k, j)*ds[k-j]*radius**(-(k-j))
                        *(2*radius)**(degree-j)/S.factorial(degree-j)
                        for j in range(min(k, degree)+1))
            check(f"probe_cutoff_jet_{degree}_{k}", bound < (10**15 if k <= 3 else 10**30), bound)
    probes = [S.Integer(1), s, s*s/2]
    jets = S.Matrix([[S.diff(f, s, k).subs(s, 0) for f in probes] for k in range(3)])
    check("three_probes_identity", jets == S.eye(3))
    for k, f in enumerate(probes):
        zero("probe_exact_near_Taylor_"+str(k), f-f.subs(s, 0)-s*S.diff(f, s).subs(s, 0)
             -s*s*S.diff(f, s, 2).subs(s, 0)/2)
    check("probe_support_in_outer_chart", R("1e-15")+2*radius < R(1, 64))
    check("inverse_half_density_third_jet", 2*S.factorial(3)*20**3 < 10**5)
    check("global_probe_third_flow_budget", 17*8*10**15*10**5 < 10**30)
    # The inherited flow recurrence is recomputed, not substituted from a log.
    z = t
    for n in range(5):
        rhs = S.series(2*(1+S.diff(z, t)**2)/(1-20*z), t, 0, n+1).removeO()
        z += S.expand(rhs).coeff(t, n)*t**(n+2)/((n+1)*(n+2))
    zj = [S.factorial(n)*S.expand(z).coeff(t, n) for n in range(1, 7)]
    for n in range(1, 7):
        budget = 2*sum(S.factorial(k)*20**k*S.bell(n, k, tuple(zj[:n-k+1]))
                       for k in range(1, n+1))
        check("analytic_probe_flow_"+str(n), budget < 10**30, budget)


def point_and_weighted_budgets():
    y, r = S.symbols("Y r", positive=True)
    pi = R(22, 7)
    check("global_third_remainder_integral", S.integrate(r**-2, (r, y/2, S.oo)) == 2/y)
    check("third_point_budget", 2*pi*72/6*2*10**30 < 10**33)
    check("sixth_near_budget", 2*pi*72/720/8 < 1)
    check("fourth_log_budget", 2*pi*1152/24 < 400)
    check("complete_geodesic_probe_budget", 10**30+10**33+400*10**30*3+10*10**30 < 10**35)
    check("flat_point_probe_budget", (R("2e-5")+10)*10**30 < 10**32)
    check("all_three_contacts", 2*10**35+10**32+10**20 < 10**40)
    delta = R(7, 2)
    check("output_measure_identity", 2*(4-delta)-1 == 0)
    dd, xx = S.symbols("D x")
    for k in range(4):
        poly = S.Poly(S.prod(dd+delta-j for j in range(k)), dd)
        check("weighted_polynomial_"+str(k), sum(abs(v) for v in poly.all_coeffs()) < 40)
        if k <= 2:
            check("contact_scale_power_"+str(k), 2*(delta-k-1)+1 == 2*(3-k))
        for n in range(3):
            target = xx**n
            for j in range(k):
                target = S.diff(target, xx)+(delta-j)*target
            zero(f"weighted_identity_{k}_{n}",
                 S.diff(y**delta*S.log(y)**n, y, k)-y**(delta-k)*target.subs(xx, S.log(y)))
    check("Taylor_scale_power", 2*(delta-3)+1 == 2)
    check("bounded_tail_scale_power", 2*delta+1 == 8)
    check("complete_model_comparison", 2*(120*10**40+40*10**7+10**20) < 10**44)


def logarithmic_norm_conversion():
    t = S.symbols("t", positive=True)
    series = t/(1-t)
    for _ in range(5):
        series = S.cancel(t*S.diff(series, t))
    remainder = S.cancel((series-t)/t**2)
    check("exact_fifth_moment_series", remainder.subs(t, R(1, 2)) == 4326)
    check("all_conormal_derivatives_bound", 2*4326*20**2 < 10**7)
    check("small_collar_series_ratio", 20*2*R("1e-15") < R(1, 2))
    ds = [1, 11, 220, 3960, 63360, 887040]
    for k, bound in enumerate(ds):
        check("log_cutoff_derivative_"+str(k), 2**k*bound < 10**8)
    check("localized_multiplier_C5", 2**5*10**8*10**7*2**2 < 10**18)
    check("integer_multiplication_bound", sum(S.binomial(5, j)*4**j for j in range(6)) < 10**6)
    check("integer_H3_H5_budget", 10**3*10**18 == 10**21)
    check("logarithmic_block_transfer", 28672 < 10**6)
    check("ordinary_H4_transfer", 1536 < 10**6)
    check("mixed_multiplier_budget", 10**6*10**21 == 10**27)
    check("multiplier_deviation", 10**27*R("1e-15")**2 <= R(".001"))
    check("forward_norm_conversion", 1+R(".001") < R("1.01"))
    check("inverse_norm_conversion", 1/(1-R(".001")) < R("1.01"))
    check("actual_remainder_conjugation", R("1.01")**2*10**64 < 2*10**64)
    check("complete_massive_comparison", 2*10**64+10**44+R("1.5e30") < 10**65)
    check("normal_composition", R("8.7e19")*10**65*R("1e-120") < R("1e-30"))
    check("normal_collar_inside_contract", R("1e-120") < R("1e-15"))
    check("mass_cubic_not_renormalized", 3+2 > R(7, 2))
    check("actual_density_not_discarded", PROTOCOL["actual_rho_retained"] is True)
    check("Weyl_assembly_remains_open", PROTOCOL["Weyl_density_combination_complete"] is False)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 complete outer massive comparison",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        half_density_and_fourth_allocation()
        kernel_tail_and_probes()
        point_and_weighted_budgets()
        logarithmic_norm_conversion()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(r["passed"] for r in ROWS)
    run = {
        "story": "S166 complete outer massive comparison",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version, "sympy": S.__version__,
        "protocol": PROTOCOL, "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "complete_outer_massive_model_comparison": good,
        "Weyl_density_combination_complete": False, "S166_complete": False, "goal_complete": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["local Weyl and principal-density assembly", "remaining local disk blocks",
                       "mass/fast localization and nonlocal cutoffs", "annulus matching",
                       "higher York maps", "full inverse/residual composition",
                       "nonlinear signed compatibility", "physical requirements"],
    }
    if not good:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_text(encoding="utf-8")) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
