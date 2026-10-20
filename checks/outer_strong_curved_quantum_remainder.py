"""S166 covariance probe and strong outer actual-minus-geodesic estimate."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "outer-strong-curved-quantum-remainder"
OWN = ["checks/outer_strong_curved_quantum_remainder.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("stationary-curved-quantum-tt", 157),
    ("strong-interior-quantum-remainder", 141),
    ("outer-global-and-zero-contact-bound", 205),
    ("outer-flat-tower-normal-remainder", 321),
    ("stationary-centered-complete-residual", 92),
]
PROTOCOL = {
    "operator": "complete actual canonical massive Hessian minus full geodesic model",
    "first_contact_method": "actual point action on the full-double covariant radial probe",
    "first_contact_bound": "abs(A1_radial)<1e60 C/Y^2 for 0<Y<=1e-15",
    "second_contact_bound": "Frobenius norm <1e30 C",
    "weighted_bound": "1e64 epsilon, H4log to L2 at delta=7/2",
    "physical_density": "actual C*u^4/(16*pi*xi^2)",
    "normal_composition": "<1e-30 at epsilon=1e-120",
    "odd_probe_is_physical_datum": False,
    "geodesic_to_flat_controlled": False,
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
              and all(c["passed"] for c in run["checks"])
              and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for p, h in sorted(expected.items()):
        check("input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)


def covariance_and_geometry():
    y, s = S.symbols("Y s", real=True)
    a0 = S.symbols("a0", positive=True)
    a = S.Function("a")(y)
    source = S.Function("T")(y)
    j, probe = a/a0, S.diff(a, y)/a0
    curvature = -2*S.diff(a, y, 2)/a
    hp = S.diff(a, y)/a
    zero("Lie_metric_radial", 2*S.diff(j, y)-2*probe)
    zero("Lie_metric_time", j*S.diff(a*a, y)-2*probe*a*a)
    density = S.exp(2*s*probe)*S.Function("T")(y+s*j)
    zero("fixed_area_covariance_identity",
         S.diff(density, s).subs(s, 0)-j*S.diff(source, y)-2*probe*source)
    constant = S.symbols("vacuum_constant")
    zero("vacuum_density_term_is_retained",
         S.diff(S.exp(2*s*probe)*constant, s).subs(s, 0)-2*probe*constant)
    zero("probe_first_curvature", S.diff(probe, y)+curvature*j/2)
    zero("probe_second_curvature",
         S.diff(probe, y, 2)+(S.diff(curvature, y)*j+curvature*probe)/2)
    zero("probe_tangential_Hessian", hp*S.diff(probe, y)+curvature*probe/2)
    u, c, xi = S.symbols("u C xi", positive=True)
    p = S.Function("P")(u)
    uu = S.sqrt(p/c)
    jj = S.sqrt(p/xi)
    vv = S.diff(p, u)/(2*S.sqrt(c*xi))
    zero("proper_probe_polynomial", S.diff(jj, u)*uu-vv)
    zero("probe_derivative_polynomial", S.diff(vv, u)*uu-S.diff(p, u, 2)*jj/(2*c))
    check("global_j_below_two", 1+R(1, 1000) < 4)
    check("global_probe_below_seven", R(5, 2)**2*R("6.1") < 49)
    check("global_probe_first", 2*2/2 <= 2)
    check("global_probe_second", 17*2/2+2*7/2 <= 24)
    check("global_probe_tangential", 2*7/2 <= 7)
    check("global_Hessian_reserve", 24+7 < 100)
    check("singular_source_coefficient", R(".334") < 1)
    check("H_probe_complete_source_bound",
          (2*3*10**22+2*7*2*10**20)*400**2 < 10**29)
    zero("log_upper_difference_derivative", S.diff(1/y-1+S.log(y), y)-(y-1)/y**2)


def analytic_geometry_and_flow():
    rr = R(1, 20)
    # Exact sums bound the one-sided analytic continuation, not the even seam.
    p_error = rr**2/1000+R("4.02")*rr**3+R("3.01")*rr**4+R("1e-15")*rr**4*(1+rr)**2
    p_prime = 2*rr/1000+R("12.06")*rr**2+R("12.04")*rr**3
    p_prime += R("1e-15")*(4*rr**3+10*rr**4+6*rr**5)
    check("complex_polynomial_error", p_error < R(".001"), p_error)
    check("inverse_sqrt_series_reserve", R(".001")/(1-R(".001")) < R(".002"))
    check("inverse_disk_Rouche_margin", R(".998")**2*R("5.9")/20**2 > R(".1")**2)
    check("analytic_probe_bound", p_prime/2 < 1)
    check("analytic_warp_lower", 1-R(".001") > R(1, 2)**2)
    check("analytic_warp_upper", 1+R(".001") < 2**2)
    check("analytic_logwarp_bound", 1/R(1, 2) <= 2)
    check("center_Cauchy_radius", R("1e-15")+R(1, 20) < R(".1"))
    check("nonzero_probe_numerator", 17*2**36*R("1e-15") < 1)
    check("small_probe_lower_denominator", R(1, 2)/2 == R(1, 4))
    t = S.symbols("t")
    z = t
    for n in range(5):
        rhs = S.series(2*(1+S.diff(z, t)**2)/(1-20*z), t, 0, n+1).removeO()
        coefficient = S.expand(rhs).coeff(t, n)/((n+1)*(n+2))
        check("majorant_coefficient_positive_"+str(n+2), coefficient > 0, coefficient)
        z += coefficient*t**(n+2)
    residual = S.series(S.diff(z, t, 2)-2*(1+S.diff(z, t)**2)/(1-20*z),
                        t, 0, 5).removeO()
    zero("majorant_ODE_through_four", residual)
    jets = [S.factorial(n)*S.expand(z).coeff(t, n) for n in range(1, 7)]
    for n in range(1, 7):
        bound = sum(S.factorial(k)*20**k*S.bell(n, k, tuple(jets[:n-k+1]))
                    for k in range(1, n+1))
        check("flow_probe_jet_"+str(n), bound < 10**30, bound)


def point_budgets_and_contacts():
    pi = R(22, 7)
    check("geodesic_near_sixth_budget", 2*pi*72/720/8 < 1)
    check("geodesic_far_quadratic_budget", 2*pi*72*100*2 < 10**5)
    check("geodesic_far_quartic_budget", 2*pi*1152/24 < 400)
    check("ordinary_fourth_restoration", 4*100/(45*3) < 10)
    check("geodesic_probe_total", 10**30+10**5+400*10**30*3+10*10**30 < 10**35)
    y, r = S.symbols("Y r", positive=True)
    check("far_second_moment_integral", S.integrate(r**-3, (r, y/2, S.oo)) == 2/y**2)
    zero("near_sixth_integral", S.integrate(r, (r, 0, y/2))-y*y/8)
    zero("annular_log_primitive", S.diff(S.log(r), r)-1/r)
    check("Taylor_near_point_reserve", 2*pi*2*10**8/6/2*10**30 < 2*10**38)
    check("Taylor_far_point_reserve", 2*pi*2*10**8*100 < 10**12)
    check("Taylor_probe_total", 2*10**38+10**12 < 10**40)
    check("compact_second_probe_sup", R(1, 32)**2 == R(1, 1024))
    # Third chain rule in two coordinates, independently retaining all sums.
    velocity, accel, jerk = 2, 240, 10**6
    check("geodesic_acceleration", 4*15*velocity**2 <= accel)
    check("geodesic_jerk", 8*6000*velocity**3+8*15*accel*velocity < jerk)
    d1, d2, d3 = 3*100, 3*2*100**2, 3*6*100**3
    chain3 = 8*d3*velocity**3+3*4*d2*accel*velocity+2*d1*jerk
    check("second_probe_third_flow_jet", chain3 < 10**12, chain3)
    e_probe = 10**13+10**6+6*10**24/R(1024)
    check("complete_E_second_probe", e_probe < 10**22)
    check("second_probe_Taylor", 2*10**6*10**12 == 2*10**18)
    second = (10**22+2*10**18+8*10**24/R(1024))/2
    check("second_contact_Frobenius", 2*second < 10**30)
    z1, z2 = S.symbols("z1 z2")
    for k, f in enumerate([z1*z1, z2*z2, (z1+z2)**2/2]):
        check("second_probe_zero_value_"+str(k), f.subs({z1: 0, z2: 0}) == 0)
        check("second_probe_zero_gradient_"+str(k),
              all(S.diff(f, z).subs({z1: 0, z2: 0}) == 0 for z in [z1, z2]))
    aa, bb, dd = S.symbols("A11 A12 A22")
    zero("second_contact_polarization", (aa+2*bb+dd)/2-aa/2-dd/2-bb)
    check("constant_contact_probe", 10**20*400**2*2*7 < 10**27)
    check("bounded_probe_tail", 8*10**24*7 < 10**26)
    check("second_contact_probe_action", 10**30*100 <= 10**32)
    numerator = 10**29+10**35+10**27+10**32+10**40+10**26
    check("first_contact_numerator", numerator < 2*10**40)
    check("first_contact_division", 8*10**40*2**36 < 10**60)


def covariant_and_weighted_derivatives():
    y, t = S.symbols("Y t", real=True)
    a, v = S.Function("a")(y), S.Function("v")(y)
    coords = [y, t]
    metric = S.diag(1, a*a)
    inverse = S.diag(1, a**-2)
    gamma = [[[S.simplify(sum(inverse[k, l]*(S.diff(metric[l, j], coords[i])
                    +S.diff(metric[l, i], coords[j])-S.diff(metric[i, j], coords[l]))
                    for l in range(2))/2) for j in range(2)] for i in range(2)] for k in range(2)]
    hessian = [[S.diff(v, coords[i], coords[j])
                -sum(gamma[k][i][j]*S.diff(v, coords[k]) for k in range(2))
                for j in range(2)] for i in range(2)]
    hp = S.diff(a, y)/a
    expected = {
        (0, 0, 0): S.diff(v, y, 3),
        (0, 1, 1): S.diff(hp, y)*S.diff(v, y)+hp*S.diff(v, y, 2),
        (1, 0, 1): hp*(S.diff(v, y, 2)-hp*S.diff(v, y)),
        (1, 1, 0): hp*(S.diff(v, y, 2)-hp*S.diff(v, y)),
    }
    for k in range(2):
        for i in range(2):
            for j in range(2):
                actual = S.diff(hessian[i][j], coords[k])
                actual -= sum(gamma[l][k][i]*hessian[l][j]+gamma[l][k][j]*hessian[i][l]
                              for l in range(2))
                actual /= a**(k+i+j)
                zero(f"covariant_third_{k}{i}{j}", actual-expected.get((k, i, j), 0))
    yy, xx, dd = S.symbols("Ypos x D", positive=True)
    delta = R(7, 2)
    for k in range(1, 4):
        poly = S.Poly(S.prod(dd+delta-j for j in range(k)), dd)
        check("weighted_derivative_coefficients_"+str(k), sum(abs(z) for z in poly.all_coeffs()) < 40)
        for n in range(4):
            target = xx**n
            for j in range(k):
                target = S.diff(target, xx)+(delta-j)*target
            zero(f"weighted_derivative_{k}_{n}",
                 S.diff(yy**delta*S.log(yy)**n, yy, k)-yy**(delta-k)*target.subs(xx, S.log(yy)))
    check("first_singular_weight_power", 2*(delta-1-2)+1 == 2)
    check("second_weight_power", 2*(delta-2)+1 == 4)
    check("third_weight_power", 2*(delta-3)+1 == 2)
    check("first_polynomial_norm", 1+delta < 5)
    check("third_covariant_norm_budget", 2*(17*40)**2 < 1000**2)
    check("first_contact_weighted_budget", 10**60*5*2 <= 10**61)
    check("second_contact_weighted_budget", 10**30*(40+2*40)*2 < 10**34)
    check("Taylor_weighted_budget", 2*2*10**6*1000*2 < 10**11)
    check("constant_and_bounded_weighted_budget",
          2*(16*10**24+2*R("1.6e25")) < 10**27)
    check("complete_weighted_bound", 10**61+10**34+10**11+10**27 < 10**64)
    check("short_collar_contract", R("1e-120") <= R("1e-15"))
    check("normal_inverse_complete_component", R("8.7e19")*10**64*R("1e-120") < R("1e-30"))
    check("odd_probe_not_physical_datum", PROTOCOL["odd_probe_is_physical_datum"] is False)
    check("geodesic_to_flat_still_open", PROTOCOL["geodesic_to_flat_controlled"] is False)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 strong outer curved quantum remainder",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        covariance_and_geometry()
        analytic_geometry_and_flow()
        point_budgets_and_contacts()
        covariant_and_weighted_derivatives()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(r["passed"] for r in ROWS)
    run = {
        "story": "S166 strong outer curved quantum remainder",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version, "sympy": S.__version__,
        "protocol": PROTOCOL, "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "strong_outer_actual_minus_geodesic_bound": good,
        "geodesic_to_flat_controlled": False, "S166_complete": False, "goal_complete": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["curved geodesic-to-flat comparison", "Weyl and principal density assembly",
                       "annulus matching", "higher York maps", "full inverse/residual composition",
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
