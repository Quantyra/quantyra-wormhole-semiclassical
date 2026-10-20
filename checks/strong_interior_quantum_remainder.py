"""Frozen S166 strong interior quantum remainder verification."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "strong-interior-quantum-remainder"
OWN = ["checks/strong_interior_quantum_remainder.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("stationary-curved-quantum-tt", 157),
    ("general-interior-dirac-response", 61),
    ("coincident-normal-dirac-response", 46),
]
PROTOCOL = {
    "scope": "actual canonical massive Hessian minus retained geodesic model, Y>=1/8",
    "norm": "sum of covariant H3 seminorms in gbar on one physical half",
    "constant": "1e48 C; intentionally unoptimized",
    "method": "first Cauchy, actual zero-value probes, finite-part Taylor subtraction",
    "quantum_domain": "full paired actual domain, angular tower and common contacts",
    "arithmetic": "exact symbolic and rational, continuum derivative bounds",
    "global_fast_TT_inverse": False,
    "Newton_contraction": False,
}
R = S.Rational
ROWS = []


def check(name, ok, detail=""):
    assert name not in {v["name"] for v in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    expression = S.factor(S.cancel(S.expand(expression)))
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
        out.update([f"notes/{stem}.md", f"evidence/{stem}.json"])
        out.update(inputs(latest(stem)))
    return sorted(out)


def predecessors():
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed"
              and len(run["checks"]) == count
              and all(c["passed"] for c in run["checks"])
              and inputs(run) == inputs(run, True))
        for path, digest in inputs(run).items():
            check(stem+"_input_"+path, sha256((ROOT/path).read_bytes()).hexdigest() == digest)


def cauchy_and_density():
    total = R(2, 3)*2*10**19 + R("1.4e14")*24*256**3 + 2*10**19
    check("first_Cauchy_complete_tower", total < 6*10**22, total)
    s, v, t0, t1 = S.symbols("s v t0 t1")
    first = S.diff(S.exp(2*s*v)*(t0+s*t1), s).subs(s, 0)
    zero("fixed_area_density_contact", first-(t1+2*v*t0))
    zero("zero_value_probe_density_contact", first.subs(v, 0)-t1)
    # Independent quadratic form polarization of the two-jet contact.
    x, y = S.symbols("x y")
    probes = [x, y, x*x/2, y*y/2, x*y]
    jet = lambda p: S.Matrix([S.diff(p, x), S.diff(p, y), S.diff(p, x, 2),
                              S.diff(p, y, 2), S.diff(p, x, y)]).subs({x: 0, y: 0})
    matrix = S.Matrix.hstack(*(jet(p) for p in probes))
    check("five_probes_extract_identity", matrix == S.eye(5))
    check("all_probes_have_zero_value", all(p.subs({x: 0, y: 0}) == 0 for p in probes))


def cutoff_and_flow():
    t = S.symbols("t", real=True)
    poly = 462*t**6-1980*t**7+3465*t**8-3080*t**9+1386*t**10-252*t**11
    zero("beta_integral_derivative", S.diff(poly, t)-2772*t**5*(1-t)**5)
    zero("cutoff_left_value", poly.subs(t, 0))
    zero("cutoff_right_value", poly.subs(t, 1)-1)
    d = [11, 220, 3960, 63360, 887040, 10644480]
    for j, bound in enumerate(d, 1):
        pp = S.Poly(S.diff(poly, t, j), t)
        degree = pp.degree()
        coeff = {i[0]: c for i, c in pp.terms()}
        bernstein = [sum(coeff.get(i, 0)*S.binomial(k, i)/S.binomial(degree, i)
                         for i in range(k+1)) for k in range(degree+1)]
        check("Bernstein_derivative_"+str(j), max(abs(b) for b in bernstein) <= bound)
        if j < 6:
            zero("left_join_"+str(j), pp.as_expr().subs(t, 0))
            zero("right_join_"+str(j), pp.as_expr().subs(t, 1))
    # Independent Cartesian differentiation with r as an auxiliary variable.
    x, y, radius = S.symbols("x y radius")
    deriv = lambda p, z: S.diff(p, z)+S.diff(p, radius)*z/radius
    radial = [1, 2, 6, 36, 300, 4000]
    for j in range(1, 7):
        for nx in range(j+1):
            value = radius
            for _ in range(nx):
                value = deriv(value, x)
            for _ in range(j-nx):
                value = deriv(value, y)
            unit = S.Poly(S.expand(value.subs(radius, 1)), x, y)
            check(f"radial_Cartesian_{j}_{nx}",
                  sum(abs(c) for c in unit.coeffs()) <= radial[j-1])
    cut = [S.Integer(1)]
    for j in range(1, 7):
        cut.append(2**j*sum(d[k-1]*S.bell(j, k, tuple(radial[:j-k+1]))
                            for k in range(1, j+1)))
    for j in range(7):
        linear = R(64)**(j-1)*(cut[j]+(j*cut[j-1] if j else 0))
        quadratic = R(64)**(j-2)*(cut[j]+(2*j*cut[j-1] if j else 0)
                                  +(j*(j-1)*cut[j-2] if j >= 2 else 0))
        bound = max(2, linear, quadratic)
        check("probe_six_jet_"+str(j), bound < 10**22)
        if j <= 5:
            # Metric conversion of the global gradient uses exp(.027)<2.
            check("probe_U_"+str(j), 2*bound < 10**18)
    # Christoffel analytic majorant: remaining bidisk and real-coordinate conversion.
    check("log_seed_margin", R(26, 1000)/(1-R(26, 1000)) < R(27, 1000))
    for j in range(5):
        upper = 3*R(27, 1000)*128*(j+1)**2*128**j
        check("Christoffel_majorant_"+str(j), upper < 100*(j+1)*128**j)
    # Solve only the finite positive formal jet needed by the geodesic equation.
    majorant = 2*t
    coefficients = []
    for degree in range(5):
        rhs = 400*(1-128*majorant)**-2*S.diff(majorant, t)**2
        coefficient = S.expand(S.series(rhs, t, 0, degree+1).removeO()).coeff(t, degree)
        majorant += coefficient*t**(degree+2)/((degree+1)*(degree+2))
        coefficients.append(coefficient)
        check("positive_geodesic_coefficient_"+str(degree), coefficient > 0)
    residual = S.series(S.diff(majorant, t, 2)
                        -400*(1-128*majorant)**-2*S.diff(majorant, t)**2, t, 0, 5).removeO()
    zero("geodesic_equation_through_degree_four", residual)
    jets = [S.diff(majorant, t, j).subs(t, 0) for j in range(1, 7)]
    for n, ceiling in [(3, 10**32), (4, 10**46), (6, 10**46)]:
        bound = 10**22*sum(2**k*S.bell(n, k, tuple(jets[:n-k+1])) for k in range(1, n+1))
        check("full_flow_probe_X"+str(n), bound < ceiling, bound)


def taylor_and_budgets():
    x, r, t = S.symbols("x r t")
    for degree in range(13):
        p = x**degree
        residual = p.subs(x, x+r)-p-r*S.diff(p, x)-r*r*S.diff(p, x, 2)/2
        integral = r**3*S.integrate((1-t)**2*S.diff(p, x, 3).subs(x, x+t*r), (t, 0, 1))/2
        zero("integral_Taylor_degree_"+str(degree), residual-integral)
    cutoff = R(1, 400)
    pi_upper = R(22, 7)
    check("strong_near_kernel_half_double", 2*4*pi_upper*10**8*cutoff/6 < 2*10**6)
    check("point_near_probe", 2*10**6*10**32 == 2*10**38)
    annulus = 4*pi_upper*10**8*S.integrate(r**-3, (r, cutoff, R(1, 128)))
    check("comparison_annulus_row", annulus < 10**15, annulus)
    tail = 2*pi_upper*18/cutoff**4
    check("complete_geodesic_tail_row", tail < 10**13)
    check("complete_bounded_kernel_rows", 3*10**24+10**15+10**13 < 4*10**24)
    check("row_subtraction_strong_Schur", 2*4*10**24 == 8*10**24)
    check("actual_constant_contact", 10**20*400**2*(1+8) < 2*10**26)
    check("actual_probe_Hessian", 6*10**22*10**18 == 6*10**40)
    check("geodesic_probe_complete_action", (1+R("2e-5"))*10**46 < 2*10**46)
    check("actual_minus_geodesic_probe", 6*10**40+2*10**46 < 3*10**46)
    check("all_five_contact_coefficients", 3*10**46+2*10**38+8*10**24 < 4*10**46)
    check("final_strong_H3_reserve",
          2*10**26+8*10**24+5*4*10**46+2*10**6 < 10**48)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 strong interior quantum remainder",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        cauchy_and_density()
        cutoff_and_flow()
        taylor_and_budgets()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 strong interior quantum remainder",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed",
        "python": sys.version, "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "S166_complete": False, "goal_complete": False,
        "strong_interior_remainder_bound": good,
        "global_fast_TT_inverse": False,
        "unresolved": ["principal graph comparison", "outer fast fourth/log estimate",
                       "higher York reconstruction and residual composition",
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
